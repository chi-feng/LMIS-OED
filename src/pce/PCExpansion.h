#ifndef PCExpansion_H
#define PCExpansion_H

#include <memory>     // std::shared_ptr
#include <vector>     // std::vector
#include <tuple>      // std::pair
#include <numeric>    // std::accumulate
#include <functional> // std::function

#include <Eigen/Core>

#include "Quadrature1D.h"
#include "Polynomial1D.h"
#include "MultiindexSet.h"

class PCExpansion {
  public:

    enum class GermType { Normal, Uniform, Exponential };

  private:

    /// Function handle
    std::function<double(const Eigen::VectorXd&)> f;

    int maxOrder;                            /// Maximum per-dimension PCE order
    size_t dim;                              /// Dimension of inputs
    std::vector<GermType> germs;             /// Germ types
    std::shared_ptr<MultiindexSet> indexSet; /// Multiindex set
    Eigen::VectorXd coeffs;                  /// PCE coefficients
    bool hasCoefficients;                    /// Flag indicating need to recompute coefficients

    /// Precomputed univariate quadrature rules indexed by [dim][order]
    std::vector<std::vector<std::shared_ptr<Quadrature1D>>> univarQuads;

    /// Precomputed multivariate polynomial square norms \f$\langle\Psi_k^2\rangle\f$
    Eigen::VectorXd multivarNorms;

    /// Evaluate the multivariate polynomial parameterized by a multiindex
    inline double EvalTensorPoly(const std::vector<int>& mindex, const Eigen::VectorXd& germ)
    {
      double product = 1.0;
      for (size_t i = 0; i < dim; ++i) {
        switch (germs[i]) {
        case GermType::Normal: product      *= HermitePoly1D::Evaluate(mindex[i], germ(i)); break;
        case GermType::Uniform: product     *= LegendrePoly1D::Evaluate(mindex[i], germ(i)); break;
        case GermType::Exponential: product *= LaguerrePoly1D::Evaluate(mindex[i], germ(i)); break;
        default: break;
        }
      }
      return product;
    }

    /// Get the quadrature order required to integrate the polyomial for a given multiindex
    inline std::vector<int> GetQuadratureOrders(const std::vector<int>& mindex)
    {
      std::vector<int> orders(dim, maxOrder);
      return orders;
    }

    /// Computes coeficients \f$c_k=\langle f_k\Psi_k\rangle/\langle\Psi_k^2\rangle\f$
    inline void ComputeCoefficients()
    {
      // for every Multiindex k, flatten the tensorized univariate quadratures by using an "odometer"
      // where each digit indexes a univariate quadrature rule. This avoids storing quadrature points
      // in memory.
      coeffs = Eigen::VectorXd::Zero(indexSet->size());
      for (size_t k = 0; k < indexSet->size(); ++k) {
        std::vector<int> quad_idx(dim, 0);                                           // odometer digits
        std::vector<int> quad_ord = GetQuadratureOrders(indexSet->GetMultiindex(k)); // max digit value
        int npts                  = std::accumulate(quad_ord.begin(), quad_ord.end(), 1, std::multiplies<int>());
        Eigen::VectorXd  point(dim);
        double weight;
        for (int i = 0; i < npts; ++i) {
          weight = 1.0;
          for (size_t j = 0; j < dim; ++j) {
            point(j) = univarQuads[j][quad_ord[j]]->abscissae(quad_idx[j]);
            weight  *= univarQuads[j][quad_ord[j]]->weights(quad_idx[j]);
          }
          coeffs(k) += f(point) * EvalTensorPoly(indexSet->GetMultiindex(k), point) * weight;
          // increment last digit and roll-over odometer if needed
          if (++quad_idx.back() >= quad_ord.back()) {
            for (int j = dim - 1; quad_idx[j] >= quad_ord[j] && j > 0; --j) {
              quad_idx[j] = 0;   // reset current digit
              quad_idx[j - 1]++; // increment next digit
            }
          }
        }
        coeffs(k) /= multivarNorms(k);
      }
      hasCoefficients = true;
    }

  public:

    /// Update the function handle
    inline void SetFunction(std::function<double(const Eigen::VectorXd&)> const& f)
    {
      this->f         = f;
      hasCoefficients = false;
    }

    /// Set the type of each germ
    inline void SetVariables(const std::vector<GermType>& germs)
    {
      this->germs     = germs;
      dim             = germs.size();
      hasCoefficients = false;
    }

    /// Constructs total order multiindex set and precomputes univariate quadratures
    inline void SetMaxOrder(const int maxOrder)
    {
      this->maxOrder  = maxOrder;
      indexSet        = std::make_shared<TotalOrderMultiindexSet>(dim, maxOrder);
      hasCoefficients = false;
      // precompute univariate quadrature rules
      univarQuads.resize(dim);
      for (size_t i = 0; i < dim; ++i) {
        univarQuads[i].resize(maxOrder + 1);
        for (int order = 1; order <= maxOrder; ++order) {
          std::shared_ptr<Quadrature1D> quad;
          switch (germs[i]) {
          case GermType::Normal: quad      = std::make_shared<HermiteQuad1D>(order); break;
          case GermType::Uniform: quad     = std::make_shared<LegendreQuad1D>(order); break;
          case GermType::Exponential: quad = std::make_shared<LaguerreQuad1D>(order); break;
          default: break;
          }
          univarQuads[i][order] = quad;
        }
      }
      // precompute multivariate polynomial norms
      multivarNorms.resize(indexSet->size());
      for (size_t k = 0; k < indexSet->size(); ++k) {
        auto   mindex  = indexSet->GetMultiindex(k);
        double product = 1.0;
        for (size_t i = 0; i < dim; ++i) {
          switch (germs[i]) {
          case GermType::Normal: product      *= HermitePoly1D::Norm(mindex[i]); break;
          case GermType::Uniform: product     *= LegendrePoly1D::Norm(mindex[i]); break;
          case GermType::Exponential: product *= LaguerrePoly1D::Norm(mindex[i]); break;
          default: break;
          }
        }
        multivarNorms(k) = product;
      }
    }

    PCExpansion() : dim(0), hasCoefficients(false) { }

    ~PCExpansion() = default;

    /// Construct the PCE given the function, variables, and maximum order
    PCExpansion(std::function<double(const Eigen::VectorXd&)> const& f,
                const std::vector<GermType>& variables, const int maxOrder)
    {
      SetFunction(f);
      SetVariables(variables);
      SetMaxOrder(maxOrder);
    }

    /// Accessor for PCE coefficients
    inline Eigen::VectorXd GetCoefficients()
    {
      if (!hasCoefficients) ComputeCoefficients();
      return coeffs;
    }

    /// Evaluate PCE for a given xi
    inline double Evaluate(const Eigen::VectorXd& xi)
    {
      if (!hasCoefficients) ComputeCoefficients();
      double sum = 0;
      for (int k = 0; k < coeffs.size(); ++k)
        sum += coeffs(k) * EvalTensorPoly(indexSet->GetMultiindex(k), xi);
      return sum;
    }

    /// Computes first and second moments of the PCE
    inline std::pair<double, double> GetMoments()
    {
      if (!hasCoefficients) ComputeCoefficients();
      double mean     = coeffs(0);
      double variance = 0;
      for (int k = 1; k < coeffs.size(); ++k) {
        variance += pow(coeffs(k), 2) * multivarNorms(k);
      }
      return std::make_pair(mean, variance);
    }

    /// Compute \f$\Sigma_{ij}=\sum_{k=1}^{K-1}c_{ik}c_{jk}\langle\Psi_k^2\rangle\f$
    inline double GetCrossCovariance(std::shared_ptr<PCExpansion> other)
    {
      if (!hasCoefficients) ComputeCoefficients();
      if (!other->hasCoefficients) other->ComputeCoefficients();
      double covariance = 0;
      for (int k = 1; k < coeffs.size(); ++k)
        covariance += this->coeffs(k) * other->coeffs(k) * multivarNorms(k);
      return covariance;
    }
};

#endif // ifndef PCExpansion_H