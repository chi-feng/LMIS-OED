#ifndef PolynomialModel_h
#define PolynomialModel_h

#include <Eigen/Core>
#include "Model.h"

namespace merr {
template<unsigned int degree>
class PolynomialModel : public Model {
  public:

    PolynomialModel() : Model()
    {
      dim = degree + 1;
    }

    using Model::Evaluate;
    inline virtual double Evaluate(const Eigen::VectorXd& x_i, const Eigen::VectorXd& lambda) override
    {
      double polynomial = lambda(0);
      for (unsigned int i = 1; i < degree + 1; ++i) {
        polynomial += lambda(i) * pow(x_i(0), i);
      }
      return polynomial;
    }
};
}

#endif // ifndef PolynomialModel_h