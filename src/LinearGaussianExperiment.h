#ifndef LinearGaussianExperiment_h
#define LinearGaussianExperiment_h

#include <Eigen/Core>
#include <memory>

#include "MultivariateNormal.h"
#include "Experiment.h"

class LinearGaussianExperiment : public Experiment
{
public:
  std::shared_ptr<MultivariateNormal> priorDist;
  std::shared_ptr<MultivariateNormal> noiseDist;

public:
  LinearGaussianExperiment();

  virtual Eigen::MatrixXd GetMatrix(const Eigen::VectorXd &design);

  virtual Eigen::VectorXd Evaluate(const Eigen::VectorXd &input, const Eigen::VectorXd &design) override;

  Eigen::VectorXd GetNoiseSample() override;

  double LogLikelihood(const Eigen::VectorXd &observation, const Eigen::VectorXd &input, const Eigen::VectorXd &design) override;

  double LogLikelihoodCached(const Eigen::VectorXd &observation,
                             const Eigen::VectorXd &input,
                             const Eigen::VectorXd &cachedOutput,
                             const Eigen::VectorXd &design) override;

  Eigen::VectorXd GetPriorSample() override;

  Eigen::VectorXd GetCondPriorSample(const Eigen::VectorXd &input) override;

  double PriorLogDensity(const Eigen::VectorXd &input) override;

  double CondPriorLogDensity(const Eigen::VectorXd &input) override;

  Eigen::VectorXd PosteriorMean(const Eigen::VectorXd &output, const Eigen::VectorXd &design) override;

  Eigen::MatrixXd PosteriorCovariance(const Eigen::VectorXd &output, const Eigen::VectorXd &design) override;

  double Evidence(const Eigen::VectorXd &observation, const Eigen::VectorXd &design);
};

#endif // ifndef LinearGaussianExperiment_h