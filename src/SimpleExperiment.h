#ifndef SimpleExperiment_h
#define SimpleExperiment_h

#include "LinearGaussianExperiment.h"

class SimpleExperiment : public LinearGaussianExperiment
{
private:
  double sigeps, k;

public:
  SimpleExperiment(const int inputDim, const int parametersOfInterest, const double sigeps, const double k = 5.0)
  {
    name = "LinearExperiment";
    designDim = 1;
    this->inputDim = inputDim;
    outputDim = inputDim;
    this->k = k;
    this->sigeps = sigeps;
    this->parametersOfInterest = parametersOfInterest;
    nuisanceParameters = inputDim - parametersOfInterest;
    priorDist = std::make_shared<MultivariateNormal>(VectorXd::Zero(inputDim), MatrixXd::Identity(inputDim, inputDim), true);
    noiseDist = std::make_shared<MultivariateNormal>(VectorXd::Zero(outputDim), MatrixXd::Identity(outputDim, outputDim) * sigeps * sigeps, true);
  }

  inline VectorXd GetPriorMean() override
  {
    return priorDist->GetMean();
  }

  inline MatrixXd GetMatrix(const VectorXd &d) override
  {
    MatrixXd G = MatrixXd::Ones(inputDim, inputDim);

    G(0, 0) = k * d(0);
    for (int i = 1; i < inputDim; ++i)
    {
      G(i, i) = k * (1.0 - d(0));
    }

    if (parametersOfInterest > 1)
    {
      // G.block(0, 0, parametersOfInterest, parametersOfInterest) = k * d(0) * MatrixXd::Identity(parametersOfInterest, parametersOfInterest);
    }

    return G;
  }
};

#endif // ifndef SimpleExperiment_h
