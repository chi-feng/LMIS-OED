#ifndef _Experiment_h
#define _Experiment_h

#include <string>
#include <memory>
#include <Eigen/Core>
#include "Distribution.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class Experiment
{
  friend class ExpectedInformationEstimator;

  public:

    std::string name;

    int inputDim, outputDim;

    int designDim;

    int parametersOfInterest, nuisanceParameters;

    bool hasPosterior;

  public:

    virtual VectorXd Evaluate(const VectorXd& input, const VectorXd& design) = 0;

    virtual VectorXd GetNoiseSample() = 0;

    virtual double   LogLikelihood(const VectorXd& output, const VectorXd& input, const VectorXd& design) = 0;

    virtual double   LogLikelihoodCached(const VectorXd& output, const VectorXd& input, const VectorXd& cachedOutput,
                                         const VectorXd& design) = 0;

    virtual VectorXd GetPriorSample() = 0;

    virtual double   PriorLogDensity(const VectorXd& input) = 0;

    virtual VectorXd GetCondPriorSample(const VectorXd& input) = 0;

    virtual double   CondPriorLogDensity(const VectorXd& nuisance) = 0;

    virtual VectorXd PosteriorMean(const VectorXd& output, const VectorXd& design) = 0;

    virtual MatrixXd PosteriorCovariance(const VectorXd& output, const VectorXd& design) = 0;

    virtual VectorXd GetPriorMean() = 0;

    inline int       GetInputDim() { return inputDim; }

    Experiment() : hasPosterior(false) { }

    virtual ~Experiment() = default;
};

#endif // ifndef _Experiment_h