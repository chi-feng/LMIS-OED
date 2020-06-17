#ifndef _Model_h
#define _Model_h
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <Eigen/Core>
#include "PCExpansion.h"

namespace merr {
using namespace Eigen;
using namespace std;

class Model {
  public:

    enum class LikelihoodType { ABC, Full, Gauss, Marginal, GaussMarginal };

    string name;
    MatrixXd x;                          /// matrix where each row is x_i
    VectorXd y;                          /// vector where each entry is y_i
    int N;                               /// number of data points
    int dim;                             /// dimension of lambda
    double datanoise;                    /// stdev of data noise when not inferring sigma
    bool inferNoise;                     /// infer data noise assume alpha.tail(1) is log(sigma)
    LikelihoodType likelihood;           /// type of likelihood, e.g. "full" or "gaussmarg"
    int nsamps;                          /// # of MC samples for estimating moments and likelihood
    bool initialized;                    /// flag to initialize shared PCE, etc.
    bool useNISP;                        /// use NISP to find model moments, otherwise use MC
    int pcdim;                           /// maximum order of PC expansion when using NISP
    vector<shared_ptr<PCExpansion>> pce; /// PCE expansions at each x_i
    double abceps;                       /// ABC epsilon
    double abcgamma;                     /// ABC gamma

    virtual void                 SetData(const MatrixXd& x, const MatrixXd& y);
    virtual VectorXd             GetLambda(const VectorXd& alpha, const VectorXd& xi);
    virtual double               Evaluate(const VectorXd& x_i, const VectorXd& lambda) = 0;
    inline virtual double        Evaluate(const VectorXd& x_i, const VectorXd& alpha, const VectorXd& xi) { return Evaluate(x_i, GetLambda(alpha, xi)); }
    virtual pair<double, double> GetMoments(const int x_i, const VectorXd& alpha);
    virtual double               LogPrior(const VectorXd& alpha);
    virtual double               LogLikelihood(const VectorXd& alpha);
    inline virtual double        LogPosterior(const VectorXd& alpha) { return LogPrior(alpha) + LogLikelihood(alpha); }
    inline virtual int           GetMCMCDim()                        { return dim + (dim * (dim + 1) / 2) + (inferNoise ? 1 : 0); }
    friend ostream             & operator<<(ostream& os, const shared_ptr<Model>& model);
};

inline std::ostream& operator<<(std::ostream& os, Model::LikelihoodType likelihood)
{
  switch (likelihood) {
  case Model::LikelihoodType::ABC: os << "ABC"; break;
  case Model::LikelihoodType::GaussMarginal: os << "GaussMarginal"; break;
  case Model::LikelihoodType::Marginal: os << "Marginal"; break;
  case Model::LikelihoodType::Gauss: os << "Gauss"; break;
  case Model::LikelihoodType::Full: os << "Full"; break;
  default: os << "Unknown"; break;
  }
  return os;
}
}

#endif // ifndef _Model_h