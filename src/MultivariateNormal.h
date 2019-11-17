#ifndef MultivariateNormal_h
#define MultivariateNormal_h

#include <memory>
#include <tuple>
#include <Eigen/Core>
#include "Distribution.h"

class MultivariateNormal : public Distribution
{
private:
  Eigen::VectorXd mean;
  Eigen::MatrixXd covariance;
  Eigen::MatrixXd cholesky;
  double cholSum;
  double constant;
  bool isDiagonal;
  Eigen::VectorXd diagLogNormalization;

public:
  bool IsDiagonal();
  void SetMean(const Eigen::VectorXd &mean);
  void SetCovariance(const Eigen::MatrixXd &covariance, const bool isDiagonal = false);
  Eigen::MatrixXd GetCovariance();
  Eigen::VectorXd GetMean();
  MultivariateNormal(const Eigen::VectorXd &mean, const Eigen::MatrixXd &covariance);
  MultivariateNormal(const Eigen::VectorXd &mean, const Eigen::MatrixXd &covariance, const bool isDiagonal);
  MultivariateNormal(const int dim);
  double LogDensity(const Eigen::VectorXd &x) override;
  static double LogDensity(const Eigen::VectorXd &x, const Eigen::VectorXd &mean, const Eigen::MatrixXd &covariance);
  Eigen::VectorXd GetSample() override;
  Eigen::VectorXd GetSamples(const int n) override;
  static std::pair<Eigen::VectorXd, Eigen::MatrixXd> GetConditional(const Eigen::VectorXd &mu, const Eigen::MatrixXd &sigma, const Eigen::VectorXd &upper);
  std::shared_ptr<MultivariateNormal> GetConditionalDist(const Eigen::VectorXd &upper);
};

#endif // ifndef MultivariateNormal_h
