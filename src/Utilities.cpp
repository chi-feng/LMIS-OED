#include "Utilities.h"

#include <iostream>
#include <fstream>

using namespace std;
using namespace Eigen;

namespace Utilities
{

MatrixXd SampleCovariance(const MatrixXd &mat)
{
  MatrixXd centered = mat.rowwise() - mat.colwise().mean();

  return (centered.adjoint() * centered) / static_cast<double>(mat.rows());
}

VectorXd WeightedSampleMean(const MatrixXd &mat, const VectorXd &w)
{
  // want row vector * matrix
  if (w.rows() > 1)
  {
    return w.transpose() * mat;
  }
  else
  {
    return w * mat;
  }
}

MatrixXd WeightedSampleCovariance(const MatrixXd &mat, const VectorXd &w)
{
  MatrixXd centered = mat.rowwise() - WeightedSampleMean(mat, w).transpose();

  return centered.transpose() * w.asDiagonal() * centered;
}

double WeightedVariance(const Eigen::VectorXd &logWeights, const Eigen::VectorXd &logSamples)
{
  VectorXd normalizedLogWeights = logWeights.array() - LogSumExp(logWeights);
  double EX = exp(LogSumExp(normalizedLogWeights + logSamples));
  double EX2 = exp(LogSumExp(normalizedLogWeights + 2 * logSamples));
  return EX2 - EX * EX;
}

/// Compute the log density of a 1D normal distribution
double NormalLogDensity(const double x, const double y, const double sigma2)
{
#define SQRT2PI 2.5066282746310005
  return -(x - y) * (x - y) / (2.0 * sigma2) - log(sqrt(sigma2) * SQRT2PI);
#undef SQRT2PI
}

double NormalLogDensity(const VectorXd &x, const VectorXd &y, const VectorXd &sigma2)
{
  double logDensity = 0;
  for (int i = 0; i < x.size(); ++i)
  {
    logDensity += NormalLogDensity(x.coeff(i), y.coeff(i), sigma2.coeff(i));
  }
  return logDensity;
}

double EffectiveSampleSize(const VectorXd &weights)
{
  return pow(weights.array().sum(), 2) / (weights.array().square().sum());
}

double CustomizedEffectiveSampleSize(const VectorXd &f, const VectorXd &w)
{
  VectorXd weights = f.cwiseAbs().array() * w.array();

  weights /= weights.array().sum();
  return 1.0 / weights.array().square().sum();
}

void WriteEigenBinaryFile(const string &path, const MatrixXd &m)
{
  ofstream file;
  file.open(path, ios::out | ios::binary);
  int rows = m.rows();
  int cols = m.cols();
  file.write(reinterpret_cast<char *>(&rows), sizeof(rows));
  file.write(reinterpret_cast<char *>(&cols), sizeof(cols));
  for (int i = 0; i < rows; ++i)
  {
    for (int j = 0; j < cols; ++j)
    {
      double entry = m(i, j);
      file.write(reinterpret_cast<char *>(&entry), sizeof(double));
    }
  }
  file.close();
}

MatrixXd ReadEigenBinaryFile(const string &path)
{
  ifstream file;

  file.open(path, ios::in | ios::binary);
  int rows, cols;
  file.read(reinterpret_cast<char *>(&rows), sizeof(rows));
  file.read(reinterpret_cast<char *>(&cols), sizeof(cols));
  MatrixXd m(rows, cols);
  for (int i = 0; i < rows; ++i)
  {
    for (int j = 0; j < cols; ++j)
    {
      double entry;
      file.read(reinterpret_cast<char *>(&entry), sizeof(double));
      m(i, j) = entry;
    }
  }
  file.close();
  return m;
}

void WriteEigenAsciiFile(const string &path, const MatrixXd &m)
{
  ofstream file;
  file.open(path);
  IOFormat full(FullPrecision, DontAlignCols, " ", "\n", "", "", "", "");
  file << m.format(full) << endl;
  file.close();
}

Eigen::MatrixXd ReadEigenAsciiFile(const std::string &filename)
{
  int cols = 0, rows = 0;
  std::vector<double> values;
  std::ifstream input(filename);
  if (input.is_open())
  {
    std::string line;
    while (!input.eof())
    {
      std::getline(input, line);
      if (line.empty())
        continue;
      std::stringstream stream(line);
      double value;
      int ncols = 0;
      while (!stream.eof())
      {
        stream >> value;
        values.push_back(value);
        ncols++;
      }
      if (ncols > 0)
        rows++;
      if (cols == 0)
        cols = ncols;
    }
    return Eigen::Map<Eigen::MatrixXd>(values.data(), cols, rows).transpose();
  }
  else
  {
    std::cerr << "Unable to open file " << filename << std::endl;
    return Eigen::MatrixXd::Zero(1, 1);
  }
}

VectorXd StringToVectorXd(const string &str)
{
  std::istringstream stream(str);
  std::string token;
  std::vector<double> values;
  double value;
  while (std::getline(stream, token, ','))
  {
    values.push_back(stod(token));
  }
  return Eigen::Map<VectorXd>(values.data(), values.size());
}

void ProgressBar(const int completed, const int total)
{
  double progress = static_cast<double>(completed) / total;
  int barWidth = 70;

  cout << "[";
  int pos = barWidth * progress;

  for (int i = 0; i < barWidth; ++i)
  {
    if (i < pos)
    {
      cout << "=";
    }
    else if (i == pos)
    {
      cout << ">";
    }
    else
    {
      cout << " ";
    }
  }
  cout << "] " << int(progress * 100.0) << " %%\r";
  cout.flush();
}
}
