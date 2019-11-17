#include <cassert>
#include <iostream>
#include <Eigen/Core>
#include <cmath>
#include <limits>
#include "merr/Model.h"
#include "RandomGenerator.h"
#include "MultivariateNormal.h"
#include "Utilities.h"

using namespace std;
using namespace Eigen;
using namespace Utilities;

namespace merr {
void Model::SetData(const MatrixXd& x, const MatrixXd& y)
{
  this->x = x;
  this->y = y;
  this->N = x.rows();

  // initialize the PCE expansions at each x_i
  if (useNISP) {
    pce.resize(N);
    vector<PCExpansion::GermType> vars(dim, PCExpansion::GermType::Normal);
    for (int i = 0; i < N; ++i) {
      pce[i] = make_shared<PCExpansion>();
      pce[i]->SetVariables(vars);
      pce[i]->SetMaxOrder(pcdim);
    }
  }
}

VectorXd Model::GetLambda(const VectorXd& alpha, const VectorXd& xi)
{
  VectorXd lambda(dim);
  for (int i = 0, d = 0; d < dim; ++d) {
    lambda(d) = alpha(i++);
    for (int j = 0; j < d + 1; ++j) {
      lambda(d) += alpha(i++) * xi(j);
    }
  }
  return lambda;
}

pair<double, double> Model::GetMoments(const int i, const VectorXd& alpha)
{
  if (useNISP) {
    auto f_i = [&](const VectorXd& xi) { return this->Evaluate(this->x.row(i), alpha, xi); };
    pce[i]->SetFunction(f_i);
    return pce[i]->GetMoments();
  } else {
    VectorXd samples(nsamps);
    for (int i = 0; i < nsamps; ++i)
      samples(i) = Evaluate(x.row(i), GetLambda(alpha, RandomGenerator::GetNormalRandomVector(dim)));
    return make_pair(samples.mean(), (samples.array() - samples.mean()).square().sum() / nsamps);
  }
}

double Model::LogPrior(const VectorXd& alpha)
{
  if ((alpha.array() > 80).any()) {
    return -numeric_limits<double>::infinity();
  }
  if (inferNoise) {
    return (alpha.tail(dim + 1).head(dim).array() < 0).any() ? -numeric_limits<double>::infinity() : 0;
  } else {
    return (alpha.tail(dim).array() < 0).any() ? -numeric_limits<double>::infinity() : 0;
  }
}

double Model::LogLikelihood(const VectorXd& alpha)
{
  double logLikelihood = 0;
  double sigma         = inferNoise ? exp(alpha.tail(1)[0]) : this->datanoise;

  switch (likelihood) {
  case LikelihoodType::ABC: {
    VectorXd logLikelihoods(N);
    #pragma omp parallel for schedule (dynamic)
    for (int i = 0; i < N; ++i) {
      auto moments = GetMoments(i, alpha);
      logLikelihoods(i) = -pow(moments.first - y(i), 2) - pow(sqrt(moments.second) - abcgamma * fabs(moments.first - y(i)), 2);
    }
    logLikelihood = logLikelihoods.sum() - N * log(abceps * sqrt(2.0 * M_PI));
  } break;

  case LikelihoodType::GaussMarginal: {
    VectorXd logLikelihoods(N);
    #pragma omp parallel for schedule (dynamic)
    for (int i = 0; i < N; ++i) {
      auto moments = GetMoments(i, alpha);
      logLikelihoods(i) = NormalLogDensity(moments.first, y(i), moments.second + pow(sigma, 2));
    }
    logLikelihood = logLikelihoods.sum();
  } break;

  case LikelihoodType::Gauss: {
    VectorXd mean(N);
    MatrixXd covariance(N, N);
    if (useNISP) {
      #pragma omp parallel for schedule (dynamic)
      for (int i = 0; i < N; ++i) {
        auto moments = GetMoments(i, alpha);
        mean(i)          = moments.first;
        covariance(i, i) = moments.second;
      }
      // separate loop so that pce[i] has coefficients for cross coavariance
      for (int i = 0; i < N; ++i) {
        for (int j = 0; j < i; ++j) {
          covariance(i, j) = pce[i]->GetCrossCovariance(pce[j]);
          covariance(j, i) = covariance(i, j);
        }
      }
    } else { // use Monte Carlo
      MatrixXd samples(N, nsamps);
      #pragma omp parallel for schedule (dynamic)
      for (int j = 0; j < nsamps; ++j) {
        auto xi     = RandomGenerator::GetNormalRandomVector(dim);
        auto lambda = GetLambda(alpha, xi);
        for (int i = 0; i < N; ++i) {
          samples(i, j) = Evaluate(x.row(i), lambda);
        }
      }
      mean       = samples.rowwise().mean();
      covariance = (samples * samples.transpose() - mean * mean.transpose() / nsamps) / nsamps;
    }
    covariance   += pow(sigma, 2) * MatrixXd::Identity(N, N);
    logLikelihood = MultivariateNormal::LogDensity(y, mean, covariance);
  } break;

  case LikelihoodType::Marginal: {
    VectorXd lambda, xi;
    for (int i = 0; i < N; ++i) {
      VectorXd logLikelihoods = VectorXd::Zero(nsamps);
      for (int j = 0; j < nsamps; ++j) {
        xi                 = RandomGenerator::GetNormalRandomVector(dim);
        lambda             = GetLambda(alpha, xi);
        logLikelihoods(j) += NormalLogDensity(Evaluate(x.row(i), lambda), y(i), pow(sigma, 2));
      }
      logLikelihood += LogSumExp(logLikelihoods) - log(nsamps);
    }
  } break;

  case LikelihoodType::Full: {
    VectorXd logLikelihoods = VectorXd::Zero(nsamps);
    VectorXd lambda, xi;
    for (int j = 0; j < nsamps; ++j) {
      xi     = RandomGenerator::GetNormalRandomVector(dim);
      lambda = GetLambda(alpha, xi);
      for (int i = 0; i < N; ++i) {
        logLikelihoods(j) += NormalLogDensity(Evaluate(x.row(i), lambda), y(i), pow(sigma, 2));
      }
    }
    logLikelihood = LogSumExp(logLikelihoods) - log(nsamps);
  } break;

  default: break;
  }

  return logLikelihood;
}

ostream& operator<<(ostream& os, const shared_ptr<Model>& model)
{
  os << "name         = \"" << model->name << "\"" << endl;
  os << "dim          = " << model->dim << endl;
  os << "mcmc dim     = " << model->GetMCMCDim() << endl;
  os << "N            = " << model->N << endl;
  os << "x            = [" << model->x.rows() << " by " << model->x.cols() << "]" << endl;
  os << "y            = [" << model->y.rows() << " by " << model->y.cols() << "]" << endl << endl;
  os << "inferNoise   = " << (model->inferNoise ? "true" : "false") << endl;
  os << " - datanoise = " << model->datanoise << endl << endl;
  os << "likelihood   = Model::LikelihoodType::" << model->likelihood << endl;
  os << " - nsamps    = " << model->nsamps << endl;
  os << " - useNISP   = " << (model->useNISP ? "true" : "false") << endl;
  os << " - pcdim     = " << model->pcdim << endl;
  os << " - abceps    = " << model->abceps << endl;
  os << " - abcgamma  = " << model->abcgamma << endl;
  return os;
}
}
