#include <cmath>
#include <Eigen/Dense>
#include "MultivariateT.h"
#include "RandomGenerator.h"

#define SQRT2PI 2.5066282746310005

using namespace Eigen;
using namespace std;

MultivariateT::MultivariateT(const VectorXd &mean, const MatrixXd &scale, const double dof) : mean(mean), scale(scale), dof(dof)
{
  dim = mean.size();
  constant = lgamma((dof + static_cast<double>(dim)) / 2.0) - lgamma(dof / 2.0) - (static_cast<double>(dim) / 2.0) * log(M_PI * dof);
  cholesky = scale.llt().matrixL();
  cholSum = cholesky.diagonal().array().log().sum();
}

VectorXd MultivariateT::GetSample()
{
  return mean + sqrt(dof / RandomGenerator::GetChiSq(dof)) * (cholesky * RandomGenerator::GetNormalRandomVector(dim));
}

VectorXd MultivariateT::GetSamples(const int n)
{
  MatrixXd iidVectors(dim, n);
  for (int i = 0; i < n; ++i)
    iidVectors.col(i) = RandomGenerator::GetNormalRandomVector(dim) * sqrt(dof / RandomGenerator::GetChiSq(dof));
  return mean + cholesky * iidVectors;
}

double MultivariateT::LogDensity(const VectorXd &x)
{
  return constant - cholSum - 0.5 * (dof + dim) * log1p(cholesky.triangularView<Lower>().solve(x - mean).squaredNorm() / dof);
}