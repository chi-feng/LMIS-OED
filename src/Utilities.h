#ifndef Utilities_h
#define Utilities_h

#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include <Eigen/Core>

#include "Timer.h"

namespace Utilities
{
Eigen::MatrixXd SampleCovariance(const Eigen::MatrixXd& mat);

Eigen::VectorXd WeightedSampleMean(const Eigen::MatrixXd& mat, const Eigen::VectorXd& w);

Eigen::MatrixXd WeightedSampleCovariance(const Eigen::MatrixXd& mat, const Eigen::VectorXd& w);

double          Variance(const Eigen::VectorXd& x);

double          WeightedVariance(const Eigen::VectorXd& logWeights, const Eigen::VectorXd& logSamples);

inline double   LogSumExp(const Eigen::VectorXd& v)
{
  return v.maxCoeff() + log(exp(v.array() - v.maxCoeff()).array().sum());
}

double          NormalLogDensity(const double x, const double y, const double sigma2);

double          NormalLogDensity(const Eigen::VectorXd& x, const Eigen::VectorXd& y, const Eigen::VectorXd& sigma2);

double          EffectiveSampleSize(const Eigen::VectorXd& weights);

double          CustomizedEffectiveSampleSize(const Eigen::VectorXd& f, const Eigen::VectorXd& w);

void            WriteEigenBinaryFile(const std::string& path, const Eigen::MatrixXd& m);

void            WriteEigenAsciiFile(const std::string& path, const Eigen::MatrixXd& m);

Eigen::MatrixXd ReadEigenBinaryFile(const std::string& path);

Eigen::MatrixXd ReadEigenAsciiFile(const std::string& filename);

Eigen::VectorXd StringToVectorXd(const std::string& str);

template<typename T>
inline T GetOption(const int argc, char **argv, const std::string& option, const T defaultOption)
{
  char **begin = argv;
  char **end   = argv + argc;
  char **itr   = std::find(begin, end, option);
  if ((itr != end) && (++itr != end)) {
    std::stringstream ss(*itr);
    T result;
    if (ss >> result) {
      return result;
    }
  }
  return defaultOption;
}

inline bool HasOption(const int argc, char **argv, const std::string& option)
{
  char **begin = argv;
  char **end   = argv + argc;

  return std::find(begin, end, option) != end;
}
}
#endif // ifndef Utilities_h