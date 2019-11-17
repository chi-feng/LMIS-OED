#include <Eigen/Dense>

#include "Utilities.h"
#include "LinearGaussianExperiment.h"

LinearGaussianExperiment::LinearGaussianExperiment()
{
  name = "LinearGaussianExperiment";
  hasPosterior = true;
}

MatrixXd LinearGaussianExperiment::GetMatrix(const VectorXd &design)
{
  return MatrixXd::Identity(inputDim, inputDim);
}

VectorXd LinearGaussianExperiment::Evaluate(const VectorXd &input, const VectorXd &design)
{
  return GetMatrix(design) * input;
}

VectorXd LinearGaussianExperiment::GetNoiseSample()
{
  return noiseDist->GetSample();
}

double LinearGaussianExperiment::LogLikelihood(const VectorXd &observation,
                                               const VectorXd &input,
                                               const VectorXd &design)
{
  return noiseDist->LogDensity(observation - Evaluate(input, design));
}

double LinearGaussianExperiment::LogLikelihoodCached(const VectorXd &observation,
                                                     const VectorXd &input,
                                                     const VectorXd &cachedOutput,
                                                     const VectorXd &design)
{
  return noiseDist->LogDensity(observation - cachedOutput);
}

VectorXd LinearGaussianExperiment::GetPriorSample()
{
  return priorDist->GetSample();
}

VectorXd LinearGaussianExperiment::GetCondPriorSample(const VectorXd &input)
{
  if (priorDist->IsDiagonal())
  {
    return priorDist->GetSample().tail(nuisanceParameters);
  }
  else
  {
    auto conditionalDist = priorDist->GetConditionalDist(input.head(parametersOfInterest));
    return conditionalDist->GetSample();
  }
}

double LinearGaussianExperiment::PriorLogDensity(const VectorXd &input)
{
  return priorDist->LogDensity(input);
}

double LinearGaussianExperiment::CondPriorLogDensity(const VectorXd &input)
{
  if (priorDist->IsDiagonal())
  {
    return Utilities::NormalLogDensity(priorDist->GetMean().tail(nuisanceParameters), input.tail(nuisanceParameters),
                                       priorDist->GetCovariance().diagonal().tail(nuisanceParameters));
  }
  else
  {
    auto conditionalDist = priorDist->GetConditionalDist(input.head(parametersOfInterest));
    return conditionalDist->LogDensity(input.tail(nuisanceParameters));
  }
}

VectorXd LinearGaussianExperiment::PosteriorMean(const VectorXd &observation, const VectorXd &design)
{
  MatrixXd G = GetMatrix(design);
  MatrixXd Gamma_pr = priorDist->GetCovariance();
  MatrixXd Gamma_obs = noiseDist->GetCovariance();
  MatrixXd Gamma_post = Gamma_pr - Gamma_pr * G.transpose() * (G * Gamma_pr * G.transpose() + Gamma_obs).ldlt().solve(G) * Gamma_pr;
  VectorXd temp = Gamma_obs.ldlt().solve(observation);
  return Gamma_post * G.transpose() * temp;
}

MatrixXd LinearGaussianExperiment::PosteriorCovariance(const VectorXd &observation, const VectorXd &design)
{
  MatrixXd G = GetMatrix(design);
  MatrixXd Gamma_pr = priorDist->GetCovariance();
  MatrixXd Gamma_obs = noiseDist->GetCovariance();
  return Gamma_pr - Gamma_pr * G.transpose() * (G * Gamma_pr * G.transpose() + Gamma_obs).ldlt().solve(G) * Gamma_pr;
}
