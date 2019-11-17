#ifndef InverseWishart_h
#define InverseWishart_h

#include <memory>
#include <tuple>
#include <Eigen/Dense>

class InverseWishart
{
  private:

    int dim;
    double dof;
    Eigen::MatrixXd S;

  public:

    inline void SetSinv(const Eigen::MatrixXd& Sinv) {
      this->dim = Sinv.rows();
      this->S = Sinv.inverse();
      double lgamma_d = 0.25 * dim * (dim - 1) * log(M_PI);
      for (int i = 1; i <= d; ++i)
        lgamma_d += tgamma((dof + 1 - i) / 2);
      this->logZ = dof * log(S.llt().matrixL().diagonal().array()).sum();
      logZ -= (dof * dim / 2) * log(2.0) - lgamma_d;
    }

    InverseWishart(const Eigen::MatrixXd& Sinv, const double dof) {
      this->dof = dof;
      SetSinv(Sinv);
    }

    double LogDensity(const Eigen::MatrixXd& X) {
      return -logZ - (dof + dim + 1) * log(X.llt().matrixL().diagonal().array()).sum() - 0.5 * (S * X.inverse()).trace();
    }
};

#endif // ifndef InverseWishart_h
