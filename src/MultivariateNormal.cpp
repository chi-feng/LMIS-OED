#include <cmath>
#include <Eigen/Dense>
#include "MultivariateNormal.h"
#include "RandomGenerator.h"

#define SQRT2PI 2.5066282746310005

using namespace Eigen;
using namespace std;

bool MultivariateNormal::IsDiagonal()
{
  return isDiagonal;
}

void MultivariateNormal::SetMean(const VectorXd &mean)
{
  this->mean = mean;
  dim = mean.size();
  constant = -0.5 * static_cast<double>(dim) * log(2.0 * M_PI);
}

void MultivariateNormal::SetCovariance(const MatrixXd &covariance, const bool isDiagonal)
{
  this->covariance = covariance;
  this->isDiagonal = isDiagonal;
  cholesky = covariance.llt().matrixL();
  cholSum = log(cholesky.diagonal().array()).sum();
  if (this->isDiagonal)
  {
    diagLogNormalization = VectorXd::Zero(dim);
    for (int i = 0; i < dim; ++i)
      diagLogNormalization(i) = log(sqrt(covariance.coeff(i, i)) * SQRT2PI);
  }
}

MatrixXd MultivariateNormal::GetCovariance() { return covariance; }
VectorXd MultivariateNormal::GetMean() { return mean; }

MultivariateNormal::MultivariateNormal(const VectorXd &mean, const MatrixXd &covariance) : isDiagonal(false)
{
  SetMean(mean);
  SetCovariance(covariance, isDiagonal);
}

MultivariateNormal::MultivariateNormal(const VectorXd &mean, const MatrixXd &covariance, const bool isDiagonal)
{
  SetMean(mean);
  SetCovariance(covariance, isDiagonal);
}

MultivariateNormal::MultivariateNormal(const int dim) : isDiagonal(true)
{
  SetMean(VectorXd::Zero(dim));
  SetCovariance(MatrixXd::Identity(dim, dim), true);
}

double MultivariateNormal::LogDensity(const VectorXd &x)
{
  double logDensity = 0;
  if (isDiagonal)
  {
    for (int i = 0; i < x.size(); ++i)
      logDensity += -pow(x.coeff(i) - mean.coeff(i), 2) / (2.0 * covariance.coeff(i, i)) - diagLogNormalization(i);
  }
  else
  {
    logDensity = constant - cholSum - 0.5 * (cholesky.triangularView<Lower>().solve(x - mean)).squaredNorm();
  }
  return logDensity;
}

double MultivariateNormal::LogDensity(const VectorXd &x, const VectorXd &mean, const MatrixXd &covariance)
{
  MatrixXd cholesky = covariance.llt().matrixL();
  return -0.5 * static_cast<double>(mean.size()) * log(2.0 * M_PI) - log(cholesky.diagonal().array()).sum() - 0.5 * (cholesky.triangularView<Lower>().solve(x - mean)).squaredNorm();
}

VectorXd MultivariateNormal::GetSample()
{
  return mean + cholesky * RandomGenerator::GetNormalRandomVector(dim);
}

VectorXd MultivariateNormal::GetSamples(const int n)
{
  MatrixXd iidNormalVectors(dim, n);
  for (int i = 0; i < n; ++i)
    iidNormalVectors.col(i) = RandomGenerator::GetNormalRandomVector(dim);
  return mean + cholesky * iidNormalVectors;
}

pair<VectorXd, MatrixXd> MultivariateNormal::GetConditional(const VectorXd &mu, const MatrixXd &sigma, const VectorXd &upper)
{
  int n_up = upper.size(), n_low = mu.size() - upper.size();
  auto ldlt = sigma.topLeftCorner(n_up, n_up).selfadjointView<Lower>().ldlt();

  return make_pair(mu.tail(n_low) + sigma.topRightCorner(n_up, n_low).transpose() * ldlt.solve(upper - mu.head(n_up)),
                   sigma.bottomRightCorner(n_low, n_low) - sigma.topRightCorner(n_up, n_low).transpose() * ldlt.solve(sigma.topRightCorner(n_up, n_low)));
}

shared_ptr<MultivariateNormal> MultivariateNormal::GetConditionalDist(const VectorXd &upper)
{
  auto conditional = MultivariateNormal::GetConditional(this->mean, this->covariance, upper);
  return make_shared<MultivariateNormal>(conditional.first, conditional.second);
}