#ifndef StochasticObjective_h
#define StochasticObjective_h

#include <memory>
#include <Eigen/Core>

class StochasticObjective {
  public:

    int dim;

    bool needInitialize;

    StochasticObjective() : needInitialize(false) { }

    inline virtual void Initialize() { }

    virtual double      Evaluate(const VectorXd& input) = 0;
};

#endif // ifndef StochasticObjective