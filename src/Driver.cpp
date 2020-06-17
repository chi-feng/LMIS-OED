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

int main(int argc, char **argv)
{
  shared_ptr<Experiment> experiment;
  string experimentType = GetOption<string>(argc, argv, "-experiment", "linear");
  if (experimentType == "linear") {
    int dim       = GetOption<int>(argc, argv, "-dim", 4);
    int poi       = GetOption<int>(argc, argv, "-poi", 1);
    double sigeps = GetOption<double>(argc, argv, "-sigeps", 0.4);
    double k      = GetOption<double>(argc, argv, "-k", 5.0);
    experiment = make_shared<SimpleExperiment>(dim, poi, sigeps, k);
  } else if (experimentType == "mossbauer")   {
    int dim       = GetOption<int>(argc, argv, "-dim", 1);
    int index     = GetOption<int>(argc, argv, "-index", 0);
    int poi       = GetOption<int>(argc, argv, "-poi", 1);
    double sigeps = GetOption<double>(argc, argv, "-sigeps", 0.4);
    experiment = make_shared<MossbauerExperiment>(dim, poi, index, sigeps);
  } else {
    cerr << "Unknown experiment type: " << experimentType << endl;
    return 1;
  }

  long seed = GetOption<long>(argc, argv, "-seed", 0L);
  if (seed > 0L) {
    RandomGenerator::SetSeed(seed);
  }
  RandomGenerator::Initialize();

  string   designStr = GetOption<string>(argc, argv, "-design", "0.0");
  VectorXd design    = StringToVectorXd(designStr);

  // do MCMC
  string   mcmcfile = GetOption<string>(argc, argv, "-mcmcfile", "none");
  if (mcmcfile != "none") {
    int      nmcmc = GetOption<int>(argc, argv, "-nmcmc", 1000);
    string   inputStr = GetOption<string>(argc, argv, "-input", "none");
    VectorXd input = (inputStr != "none") ? StringToVectorXd(inputStr) : experiment->GetPriorMean();
    VectorXd data = experiment->Evaluate(input, design);
    auto     logPosterior = [&](const VectorXd& input) {
                              return experiment->PriorLogDensity(input) + experiment->LogLikelihood(data, input, design);
                            };
    auto     mcmc = make_shared<MCMC>(logPosterior, experiment->GetInputDim());
    mcmc->SetAdaptBeta(0.1);
    mcmc->Run(nmcmc);
    auto     chain = mcmc->GetChain();
    WriteEigenAsciiFile(mcmcfile + ".chain", mcmc->GetChain());
    WriteEigenAsciiFile(mcmcfile + ".acf", mcmc->GetAutocorrelation(min(200, nmcmc)));
    return 0;
  }

  auto estimator = make_shared<ExpectedInformationEstimator>(experiment);
  estimator->N                       = GetOption<int>(argc, argv, "-N", 1000);
  estimator->M1                      = GetOption<int>(argc, argv, "-M1", 100);
  estimator->M2                      = GetOption<int>(argc, argv, "-M2", 100);
  estimator->maxComponents           = GetOption<int>(argc, argv, "-maxComponents", estimator->N);
  estimator->useIS                   = GetOption<bool>(argc, argv, "-useIS", false);
  estimator->useMIS                  = GetOption<bool>(argc, argv, "-useMIS", true);
  estimator->useExactPosterior       = GetOption<bool>(argc, argv, "-useExactPosterior", false);
  estimator->useSortHeuristic        = GetOption<bool>(argc, argv, "-useSortHeuristic", true);
  estimator->useReverseLikelihood    = GetOption<bool>(argc, argv, "-useReverseLikelihood", true);
  estimator->useMinSampleDistance    = GetOption<bool>(argc, argv, "-useMinSampleDistance", false);
  estimator->useBiasCorrection       = GetOption<bool>(argc, argv, "-useBiasCorrection", false);
  estimator->verbose                 = GetOption<bool>(argc, argv, "-verbose", false);
  estimator->biasingDistributionType = GetOption<string>(argc, argv, "-biasingDistributionType", "MVT");
  estimator->dof                     = GetOption<double>(argc, argv, "-dof", 2.5);
  estimator->nugget                  = GetOption<double>(argc, argv, "-nugget", 1e-3);
  estimator->debugCondLikelihood     = GetOption<bool>(argc, argv, "-debugCondLikelihood", false);
  estimator->useMarginal             = GetOption<bool>(argc, argv, "-useMarginal", false);

  double EIG = estimator->Evaluate(design);

  string dumpfile = GetOption<string>(argc, argv, "-dumpfile", "none");
  if (dumpfile != "none") {
    estimator->WriteToFile(dumpfile, design);
  }

  string outfile = GetOption<string>(argc, argv, "-outfile", "none");
  if (outfile != "none") {
    ofstream out;
    out.open(outfile);
    out << "design = " << design.transpose() << endl;
    out << std::setprecision(16) << "EIG = " << EIG << endl;
    out << "time = " << estimator->GetExecutionTime() << endl;
    out.close();
  }

  cout << std::setprecision(16) << "EIG = " << EIG << endl;
  cout << "time = " << estimator->GetExecutionTime() << endl;

  return 0;
}
