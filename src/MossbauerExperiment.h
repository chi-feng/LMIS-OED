#ifndef MossbauerExperiment_h
#define MossbauerExperiment_h

#include <iostream>
#include <memory>
#include <vector>
#include <limits>
#include "LinearGaussianExperiment.h"
#include "Experiment.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

class MossbauerExperiment : public LinearGaussianExperiment
{
private:
  int index;

public:
  double sigeps;

  MossbauerExperiment(const int measurements, const int poi, const int index = 0, const double sigeps = 0.4) : index(index)
  {
    name = "Mossbauer";
    outputDim = measurements;
    designDim = measurements;
    inputDim = 4;
    parametersOfInterest = poi;
    nuisanceParameters = inputDim - parametersOfInterest;
    hasPosterior = false;

    VectorXd priorMean(4);
    VectorXd priorSigma(4);

    priorMean(0) = 0.0;
    priorSigma(0) = 1.0;
    priorMean(1) = 0.0;
    priorSigma(1) = 0.3; // log normal moments (1.05, 0.10)
    priorMean(2) = 0.0;
    priorSigma(2) = 0.3; // log normal moments (1.05, 0.10)
    priorMean(3) = 1.0;
    priorSigma(3) = 0.2; // log normal moments (2.77, 0.31)

    // swap poi so that it is in front
    if (index > 0)
    {
      double temp = priorMean(0);
      priorMean(0) = priorMean(index);
      priorMean(index) = temp;
      temp = priorSigma(0);
      priorSigma(0) = priorSigma(index);
      priorSigma(index) = temp;
    }

    VectorXd priorVariance = priorSigma.array().square();

    priorDist = std::make_shared<MultivariateNormal>(priorMean, priorVariance.asDiagonal(), true);
    noiseDist = std::make_shared<MultivariateNormal>(VectorXd::Zero(outputDim), pow(sigeps, 2) * MatrixXd::Identity(outputDim, outputDim), true);
  }

  inline VectorXd GetPriorMean() override
  {
    return priorDist->GetMean();
  }

  inline VectorXd Evaluate(const VectorXd &input, const VectorXd &d) override
  {
    VectorXd input_copy = input;

    if (index > 0)
    {
      input_copy(0) = input(index);
      input_copy(index) = input(0);
    }

    VectorXd y(outputDim);
    double center = input_copy(0);
    double width = exp(input_copy(1));
    double height = exp(input_copy(2));
    double offset = exp(input_copy(3));

    auto dev = (d.array() - center).array().square();
    double wsq = width * width;

    // if width is large, then use second-order taylor expansion at infinity
    // to avoid inf / inf = nan situations, truncation error is O(width^-4)
    if (width < 1e3)
    {
      for (int i = 0; i < outputDim; ++i)
      {
        y(i) = offset - height * wsq / (dev(i) + wsq);
      }
    }
    else
    {
      for (int i = 0; i < outputDim; ++i)
      {
        y(i) = offset - height * (1.0 - dev(i) / wsq);
      }
    }
    return y;
  }
};

#endif // ifndef MossbauerExperiment_h