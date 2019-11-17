#include <memory>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <Eigen/Core>
#include "MultivariateT.h"
#include "ExpectedInformationEstimator.h"
#include "SimpleExperiment.h"
#include "MossbauerExperiment.h"
#include "Utilities.h"
#include "RandomGenerator.h"
#include "MCMC.h"

using namespace std;
using namespace Eigen;
using namespace Utilities;

double sigma2marg(const shared_ptr<Experiment> experiment, const shared_ptr<MultivariateNormal> q_ML, const VectorXd& y, const VectorXd& design, const int n) {
  VectorXd theta;
  VectorXd logSamples = VectorXd::Zero(n);
  VectorXd logWeights = VectorXd::Zero(n);
  for (int i = 0; i < n; i++) {
    theta = q_ML->GetSample();
    logWeights(i) = experiment->PriorLogDensity(theta) - q_ML->LogDensity(theta);
    logSamples(i) = experiment->LogLikelihoodCached(y, theta, experiment->Evaluate(theta, design), design);
  }
  return WeightedVariance(logWeights, logSamples);
}

double sigma2cond(const shared_ptr<Experiment> experiment, const shared_ptr<MultivariateNormal> q_CL, const VectorXd& y, const VectorXd& thetaOuter, const VectorXd& design, const int n) {
  VectorXd theta = thetaOuter;
  VectorXd logSamples = VectorXd::Zero(n);
  VectorXd logWeights = VectorXd::Zero(n);
  for (int i = 0; i < n; i++) {
    theta.tail(experiment->nuisanceParameters) = q_CL->GetSample();
    logWeights(i) = experiment->CondPriorLogDensity(theta) - q_CL->LogDensity(theta.tail(experiment->nuisanceParameters));
    logSamples(i) = experiment->LogLikelihoodCached(y, theta, experiment->Evaluate(theta, design), design);
  }
  return WeightedVariance(logWeights, logSamples);
}

int main(int argc, char **argv)
{
    RandomGenerator::Initialize();

    int dim         = GetOption<int>(argc, argv, "-dim", 4);
    int poi         = GetOption<int>(argc, argv, "-poi", 1);
    double sigeps   = GetOption<double>(argc, argv, "-sigeps", 0.4);
    double k        = GetOption<double>(argc, argv, "-k", 5.0);
    auto experiment = make_shared<SimpleExperiment>(dim, poi, sigeps, k);

    string   designStr = GetOption<string>(argc, argv, "-design", "0.8");
    VectorXd design    = StringToVectorXd(designStr);


    int outer = GetOption<int>(argc, argv, "-outer", 100);
    int inner = GetOption<int>(argc, argv, "-inner", 100);

    VectorXd D1 = VectorXd::Zero(outer);
    VectorXd D2 = VectorXd::Zero(outer);
    VectorXd D3 = VectorXd::Zero(outer);

    VectorXd C1 = VectorXd::Zero(outer);
    VectorXd C2 = VectorXd::Zero(outer);

    VectorXd CESS1 = VectorXd::Zero(outer);
    VectorXd CESS2 = VectorXd::Zero(outer);

#pragma omp parallel for
    for (int i = 0; i < outer; i++) {

      VectorXd theta = experiment->GetPriorSample();
      VectorXd y = experiment->Evaluate(theta, design) + experiment->GetNoiseSample();

      VectorXd postMean = experiment->PosteriorMean(y, design);
      MatrixXd postCov = experiment->PosteriorCovariance(y, design);
      VectorXd condMean;
      MatrixXd condCov;
      tie(condMean, condCov) = MultivariateNormal::GetConditional(postMean, postCov, theta.head(experiment->parametersOfInterest));
      auto q_ML = make_shared<MultivariateNormal>(postMean, postCov);
      auto q_CL = make_shared<MultivariateNormal>(condMean, condCov);

      int n = 10;
      VectorXd x;
      VectorXd lml = VectorXd::Zero(n);
      VectorXd lcl = VectorXd::Zero(n);
      for (int j = 0; j < n; j++) {
        x = q_ML->GetSample();
        lml(j) = experiment->LogLikelihoodCached(y, x, experiment->Evaluate(x, design), design) + experiment->PriorLogDensity(x) - q_ML->LogDensity(x);
        x = theta;
        x.tail(experiment->nuisanceParameters) = q_CL->GetSample();
        lcl(j) = experiment->LogLikelihoodCached(y, x, experiment->Evaluate(x, design), design) + experiment->CondPriorLogDensity(x) - q_CL->LogDensity(x.tail(experiment->nuisanceParameters));
      }
      double logmarglike = LogSumExp(lml) - log(n);
      double logcondlike = LogSumExp(lcl) - log(n);

      VectorXd logLikelihoods = VectorXd::Zero(inner);
      VectorXd logWeights = VectorXd::Zero(inner);
      for (int j = 0; j < inner; j++) {
        x = q_ML->GetSample();
        logWeights(j) = experiment->PriorLogDensity(x) - q_ML->LogDensity(x);
        logLikelihoods(j) = experiment->LogLikelihoodCached(y, x, experiment->Evaluate(x, design), design);
      }
      logWeights = logWeights.array() - LogSumExp(logWeights);
      double s2m = exp(LogSumExp(2 * logLikelihoods + logWeights)) - exp(2 * logmarglike);
      CESS1[i] = CustomizedEffectiveSampleSize((2 * logLikelihoods).array().exp(), logWeights.array().exp());

      x = theta;
      for (int j = 0; j < inner; j++) {
        x.tail(experiment->nuisanceParameters) = q_CL->GetSample();
        logWeights(j) = experiment->CondPriorLogDensity(x) - q_CL->LogDensity(x.tail(experiment->nuisanceParameters));
        logLikelihoods(j) = experiment->LogLikelihoodCached(y, x, experiment->Evaluate(x, design), design);
      }
      logWeights = logWeights.array() - LogSumExp(logWeights);
      double s2c = exp(LogSumExp(2 * logLikelihoods + logWeights)) - exp(2 * logcondlike);
      CESS2[i] = CustomizedEffectiveSampleSize((2 * logLikelihoods).array().exp(), logWeights.array().exp());
      
      // double s2m = sigma2marg(experiment, q_ML, y, design, inner);
      // double s2c = sigma2cond(experiment, q_CL, y, theta, design, inner);

      D1[i] = s2m / exp(2 * logmarglike);
      D2[i] = s2c / exp(2 * logcondlike);

      D3[i] = logcondlike - logmarglike;

      C1[i] = s2m / exp(logmarglike);
      C2[i] = s2c / exp(logcondlike);
    }

    cout << std::setprecision(16) << "D1 = " << D1.sum() / outer << endl;
    cout << std::setprecision(16) << "D2 = " << D2.sum() / outer << endl;
    cout << std::setprecision(16) << "D3 (mean) = " << (D3.sum() / outer) << endl;
    cout << std::setprecision(16) << "D3 (var) = " << D3.dot(D3) / outer - (D3.sum() / outer) * (D3.sum() / outer) << endl;
    cout << std::setprecision(16) << "C1 = " << (C1.sum() / outer) << endl;
    cout << std::setprecision(16) << "C2 = " << (C2.sum() / outer) << endl;
    cout << std::setprecision(16) << "CESS1 = " << (CESS1.sum() / outer) << std::setprecision(2) << " (" << (CESS1.sum() / outer / inner) << ")" << endl;
    cout << std::setprecision(16) << "CESS2 = " << (CESS2.sum() / outer) << std::setprecision(2) << " (" << (CESS2.sum() / outer / inner) << ")" << endl;


    string outfile = GetOption<string>(argc, argv, "-o", "none");
    if (outfile != "none") {
      MatrixXd m = MatrixXd::Zero(outer, 7);
      m.col(0) = C1;
      m.col(1) = C2;
      m.col(2) = D1;
      m.col(3) = D2;
      m.col(4) = D3;
      m.col(5) = CESS1;
      m.col(6) = CESS2;
      WriteEigenAsciiFile(outfile, m);
    }

    return 0;
}

