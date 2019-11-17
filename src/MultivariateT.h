#ifndef MultivariateT_h
#define MultivariateT_h

#include <Eigen/Dense>
#include "Distribution.h"

class MultivariateT : public Distribution
{
private:
  Eigen::VectorXd mean;
  Eigen::MatrixXd scale;
  double dof;
  double constant;
  Eigen::MatrixXd cholesky;
  double cholSum;
  Eigen::VectorXd delta;

public:
  MultivariateT(const Eigen::VectorXd &mean, const Eigen::MatrixXd &scale, const double dof);
  Eigen::VectorXd GetSample() override;
  Eigen::VectorXd GetSamples(const int n) override;
  double LogDensity(const Eigen::VectorXd &x) override;
};

#endif // ifndef MultivariateT_h
