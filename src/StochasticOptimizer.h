#ifndef StochasticOptimizer_h
#define StochasticOptimizer_h

#include <memory>
#include <vector>
#include <Eigen/Core>

#include "Experiment.h"
#include "StochasticObjective.h"

using std::shared_ptr;
using std::vector;
using Eigen::VectorXd;
using Eigen::MatrixXd;

class StochasticOptimizer {
  public:

    shared_ptr<StochasticObjective> objective;

    /// SPSA parameters
    double   alpha, gamma, A, a0, c0;

    /// A dim x 2 matrix where the first and second column are the minimum and maximum box constraint
    MatrixXd boxConstraints;

    /// A row-wise collection of design vectors at each iteration
    vector<VectorXd> trajectory;

    /// A row-wise collection of expected utility at each iteration
    vector<double>   utility;

    /// Dimension of the design vector
    int dim;

  public:

    StochasticOptimizer();

    void             SetObjective(shared_ptr<StochasticObjective> objective);

    void             SetBoxConstraints(const MatrixXd& boxConstraints);

    vector<VectorXd> GetTrajectory();

    vector<double>   GetUtilityTrajectory();

    /// Reset the optimizer
    void Reset();

    /// Initialize the optimizer
    void Initialize(const VectorXd& d);

    /// Take a single SPSA step
    void TakeStep();

  protected:

    /// Project box constraints in place
    void ProjectBoxConstraints(VectorXd& d, const double c_k, const VectorXd& delta_k);

    /// Enforce box constraints in place
    void EnforceBoxConstraints(VectorXd& d);
};

#endif // ifndef StochasticOptimizer_h