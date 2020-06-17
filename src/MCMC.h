#ifndef MCMC_h
#define MCMC_h

#include <iostream>
#include <iomanip>
#include <memory>
#include <functional>
#include <Eigen/Core>
#include "MultivariateNormal.h"
#include "Timer.h"

class MCMC
{
  protected:

    /// (unnormalized) log density of target distribution \f$\mathbb{R}^{\rm dim} -> \mathbb{R}\f$
    std::function<double(const Eigen::VectorXd&)> LogDensityFuncHandle;
    int dim;                      /// dimension of parameter
    Eigen::MatrixXd chain;        /// [dim x N] column-wise storage of MCMC chain
    Eigen::VectorXd logDensities; /// [N x 1] storage of log densities of each sample
    int accepted;                 /// number of accepted proposals
    std::shared_ptr<MultivariateNormal> mhProposalDist;
    std::shared_ptr<MultivariateNormal> amProposalDist;
    int adaptStride;              /// iterations between covariance updates
    double adaptBeta;             /// probability [0...1] of using nominal MH proposal
    Eigen::VectorXd chainSum;     /// running sum of chain samples for rank-1 update of sample mean
    Eigen::MatrixXd chainScatter; /// scatter matrix for rank-1 update of sample covariance
    int existing;

  public:

    inline void Reset()
    {
      chainScatter = Eigen::MatrixXd::Zero(dim, dim);
      chainSum     = Eigen::VectorXd::Zero(dim);
      accepted     = 0;
      existing     = 0;
    }

    MCMC(const std::function<double(const Eigen::VectorXd&)>& f, const int dim) : LogDensityFuncHandle(f), dim(dim)
    {
      adaptStride    = 10 * dim;
      adaptBeta      = 0.1;
      mhProposalDist = std::make_shared<MultivariateNormal>(Eigen::VectorXd::Zero(dim), 0.1 * Eigen::MatrixXd::Identity(dim, dim));
      amProposalDist = std::make_shared<MultivariateNormal>(Eigen::VectorXd::Zero(dim), 0.1 * Eigen::MatrixXd::Identity(dim, dim));
      Reset();
    }

    inline void Run(const int N)
    {
      // re-entrant?
      if (existing == 0) {
        chain           = Eigen::MatrixXd::Zero(dim, N);
        logDensities    = Eigen::VectorXd::Zero(N);
        logDensities(0) = LogDensityFuncHandle(chain.col(0));
      } else {
        chain.conservativeResize(dim, existing + N);
        logDensities.conservativeResize(existing + N);
      }
      Timer timer;
      timer.Start();
      for (int i = existing ? existing : 1; i < existing + N; ++i) {
        if (i % (N / 100) == 0) {
          std::cout << timer.ETAString(i - existing, N) << std::endl;
        }
        if (i % adaptStride == 0) {
          // adapt covariance
          auto newSamps = chain.block(0, i - adaptStride, dim, adaptStride);
          chainScatter += newSamps * newSamps.transpose();
          chainSum     += newSamps.rowwise().sum();
          amProposalDist->SetCovariance(2.4 / dim * (chainScatter - chainSum * chainSum.transpose() / i) / i);
        }
        // store proposal in-place
        chain.col(i) = (RandomGenerator::GetUniform() > adaptBeta)
                         ? chain.col(i - 1) + amProposalDist->GetSample()
                         : chain.col(i - 1) + mhProposalDist->GetSample();
        logDensities(i) = LogDensityFuncHandle(chain.col(i));
        // MH accept/reject
        if (log(RandomGenerator::GetUniform()) < logDensities(i) - logDensities(i - 1)) {
          accepted++;
        } else {
          chain.col(i)    = chain.col(i - 1);
          logDensities(i) = logDensities(i - 1);
        }
      }
      existing = existing + N;
    }

    inline Eigen::VectorXd GetAutocorrelation(const int lag)
    {
      Eigen::VectorXd mean           = chain.rowwise().mean();
      Eigen::VectorXd autocovariance = Eigen::VectorXd::Zero(lag + 1);
      for (int k = 0; k < lag + 1; ++k)
        for (int i = k; i < chain.cols(); ++i)
          autocovariance(k) += (chain.col(i) - mean).dot(chain.col(i - k) - mean);
      return autocovariance / autocovariance(0);
    }

    inline void            SetAdaptBeta(const double beta) { adaptBeta = beta; }
    inline Eigen::MatrixXd GetChain()                      { return chain.transpose(); }
    inline Eigen::VectorXd GetLogDensities()               { return logDensities; }
    inline double          GetAcceptanceRatio()            { return static_cast<double>(accepted) / chain.cols(); }
    inline Eigen::MatrixXd GetAdaptiveCov()                { return amProposalDist->GetCovariance(); }
};

#endif // ifndef MCMC_h