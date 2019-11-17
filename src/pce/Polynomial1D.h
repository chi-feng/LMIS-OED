#ifndef Polynomial1D_H
#define Polynomial1D_H

#include <cmath>

#define SQRT2PI 2.5066282746310005

/**
 * Univariate orthogonal polynomials for spectral methods
 *
 * Note: it's important to evaluate the polyomials using recurrence relations instead of directly
 *       summing their terms to minimize numerical error.
 */
class Polynomial1D {
  public:

    static double Evaluate(const int order, const double x);
    static double Norm(const int order);
    static double Var(const int order);

  protected:

    static inline long Factorial(const long n)
    {
      return (n == 0 || n == 1) ? 1 : Factorial(n - 1) * n;
    }
};

/** Univariate probabilists' Hermite polynomial */
class HermitePoly1D : Polynomial1D {
  public:

    static inline double Evaluate(const int order, const double x)
    {
      if (order == 0) {
        return 1;
      }
      if (order == 1) {
        return x;
      }
      double yn2 = 1;
      double yn1 = x;
      double y   = 0;
      for (int k = 2; k <= order; ++k) {
        y   = x * yn1 - (k - 1) * yn2;
        yn2 = yn1;
        yn1 = y;
      }
      return y;
    }

    static inline double Norm(const int order)
    {
      return Factorial(order);
    }

    static inline double Var(const int order)
    {
      return static_cast<double>(Factorial(order));
    }
};

/** Univariate Legendre polynomials */
class LegendrePoly1D : Polynomial1D {
  public:

    static inline double Evaluate(const int order, const double x)
    {
      if (order == 0) {
        return 1;
      }
      if (order == 1) {
        return x;
      }
      double yn2 = 1.0;
      double yn1 = x;
      double y   = 0;
      for (int k = 2; k <= order; ++k) {
        y   = ((2.0 * k - 1.0) * x * yn1 - (k - 1) * yn2) / k;
        yn2 = yn1;
        yn1 = y;
      }
      return y;
    }

    static inline double Norm(const int order)
    {
      return 2.0 / (2.0 * order + 1);
    }

    static inline double Var(const int order)
    {
      return 1.0 / (2.0 * order + 1.0);
    }
};

/** Univariate Laguerre polynomials */
class LaguerrePoly1D : Polynomial1D {
  public:

    static inline double Evaluate(const int order, const double x)
    {
      if (order == 0) {
        return 1;
      }
      if (order == 1) {
        return 1.0 - x;
      }
      double yn2 = 1.0;
      double yn1 = 1.0 - x;
      double y   = 0;
      for (int k = 2; k <= order; ++k) {
        y   = ((2.0 * k - 1.0 - x) * yn1 - (k - 1) * yn2) / k;
        yn2 = yn1;
        yn1 = y;
      }
      return y;
    }

    static inline double Norm(const int order)
    {
      return 1;
    }

    static inline double Var(const int order)
    {
      return 1;
    }
};

#endif // ifndef Polynomial1D_H