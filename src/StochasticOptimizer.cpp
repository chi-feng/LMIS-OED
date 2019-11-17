#include "StochasticOptimizer.h"

#include "Utilities.h"
#include "RandomGenerator.h"

using namespace Eigen;
using namespace std;

void StochasticOptimizer::SetObjective(shared_ptr<StochasticObjective> objective)
{
  this->dim       = objective->dim;
  this->objective = objective;
}

void StochasticOptimizer::SetBoxConstraints(const MatrixXd& constraints)
{
  this->boxConstraints = constraints;
}

vector<VectorXd> StochasticOptimizer::GetTrajectory()
{
  return trajectory;
}

vector<double> StochasticOptimizer::GetUtilityTrajectory()
{
  return utility;
}

StochasticOptimizer::StochasticOptimizer()
{
  alpha = 0.602;
  gamma = 0.101;
  A     = 10.0;
  a0    = 0.16;
  c0    = 0.05; // stdev of noise in loss function
}

/// Reset the optimizer
void StochasticOptimizer::Reset()
{
  trajectory.resize(0);
  utility.resize(0);
}

/// Initialize the optimizer
void StochasticOptimizer::Initialize(const VectorXd& d)
{
  Reset();
  trajectory.push_back(d);
  if (objective->needInitialize) {
    objective->Initialize();
  }
  utility.push_back(objective->Evaluate(d));
}

/// Project box constraints in place
void StochasticOptimizer::ProjectBoxConstraints(VectorXd& d, const double c_k, const VectorXd& delta_k)
{
  for (int i = 0; i < boxConstraints.rows(); ++i) {
    if (d(i) + c_k * fabs(delta_k(i)) > boxConstraints(i, 1)) {
      d(i) = boxConstraints(i, 1) - c_k * fabs(delta_k(i));
    }
    if (d(i) - c_k * fabs(delta_k(i)) < boxConstraints(i, 0)) {
      d(i) = boxConstraints(i, 0) + c_k * fabs(delta_k(i));
    }
    assert(d(i) - c_k * delta_k(i) <= boxConstraints(i, 1));
    assert(d(i) - c_k * delta_k(i) >= boxConstraints(i, 0));
    assert(d(i) + c_k * delta_k(i) <= boxConstraints(i, 1));
    assert(d(i) + c_k * delta_k(i) >= boxConstraints(i, 0));
  }
}

/// Enforce box constraints in place
void StochasticOptimizer::EnforceBoxConstraints(VectorXd& d)
{
  for (int i = 0; i < boxConstraints.rows(); ++i) {
    d(i) = fmin(d(i), boxConstraints(i, 1));
    d(i) = fmax(d(i), boxConstraints(i, 0));
  }
}

/// Take a single SPSA step
void StochasticOptimizer::TakeStep()
{
  int k = (int)trajectory.size();
  double   a_k, c_k, g_plus, g_minus;
  VectorXd delta_k, d, d_plus, d_minus, g_k;

  a_k     = a0 / pow(A + k + 1.0, alpha);
  c_k     = c0 / pow(k + 1.0, gamma);
  delta_k = RandomGenerator::GetBernoulliRandomVector(dim, 1, -1);

  d = trajectory.back();

  ProjectBoxConstraints(d, c_k, delta_k);

  d_plus  = d + c_k * delta_k;
  d_minus = d - c_k * delta_k;

  if (objective->needInitialize) {
    objective->Initialize();
  }

  g_plus  = objective->Evaluate(d_plus);
  g_minus = objective->Evaluate(d_minus);
  g_k     = (g_plus - g_minus) / (2.0 * c_k) * delta_k.cwiseInverse();

  d = d + a_k * g_k;

  EnforceBoxConstraints(d);

  trajectory.push_back(d);
  utility.push_back(objective->Evaluate(d));
}