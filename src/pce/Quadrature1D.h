#ifndef Quadrature1D_h
#define Quadrature1D_h

#include <cmath>

#include <Eigen/Core>
#include <Eigen/Eigenvalues> // for Golub-Welsch algorithm

/**
 * Univariate quadrature rules
 */
class Quadrature1D {
  public:

    Eigen::VectorXd abscissae;
    Eigen::VectorXd weights;
};

/**
 * Probabilists' Gauss-Hermite quadrature (standard normal density)
 */
class HermiteQuad1D : public Quadrature1D {
  public:

    HermiteQuad1D(const int order)
    {
      Eigen::MatrixXd companion = Eigen::MatrixXd::Zero(order, order);

      for (int i = 0; i < order; ++i) {
        if (i != 0) {
          companion(i, i - 1) = sqrt(0.5 * i);
          companion(i - 1, i) = sqrt(0.5 * i);
        }
      }
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
      solver.compute(companion);
      abscissae = sqrt(2.0) * solver.eigenvalues();
      weights   = solver.eigenvectors().row(0).array().square();
    }
};

/**
 * Gauss-Legendre quadrature (uniform density on [-1,1])
 */
class LegendreQuad1D : public Quadrature1D {
  public:

    LegendreQuad1D(const int order)
    {
      Eigen::MatrixXd companion = Eigen::MatrixXd::Zero(order, order);

      for (int i = 0; i < order; ++i) {
        if (i > 0) {
          companion(i, i - 1) = static_cast<double>(i) / sqrt(4.0 * i * i - 1.0);
          companion(i - 1, i) = static_cast<double>(i) / sqrt(4.0 * i * i - 1.0);
        }
      }
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
      solver.compute(companion);
      abscissae = solver.eigenvalues();
      weights   = 2.0 * solver.eigenvectors().row(0).array().square();
    }
};

/**
 * Gauss-Laguerre quadrature \$f(x)=x^\alpha e^{-x}\f$
 */
class LaguerreQuad1D : public Quadrature1D {
  public:

    LaguerreQuad1D(const int order, const double alpha = 0.0)
    {
      Eigen::MatrixXd companion = Eigen::MatrixXd::Zero(order, order);

      for (int i = 0; i < order; ++i) {
        companion(i, i) = (2.0 * i - 1.0) + alpha;
        if (i > 0) {
          companion(i, i - 1) = sqrt(static_cast<double>(i) * (static_cast<double>(i) + alpha));
          companion(i - 1, i) = sqrt(static_cast<double>(i) * (static_cast<double>(i) + alpha));
        }
      }
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
      solver.compute(companion);
      abscissae = solver.eigenvalues();
      weights   = tgamma(alpha + 1.0) * solver.eigenvectors().row(0).array().square();
    }
};

#endif // ifndef _Quadrature1D_h