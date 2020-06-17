#ifndef MossbauerObjective_h
#define MossbauerObjective_h

#include "StochasticObjective.h"
#include "MossbauerExperiment.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::shared_ptr;
using std::vector;

class MossbauerObjective : public StochasticObjective
{
  public:

    shared_ptr<MossbauerExperiment> experiment;
    int N, M1, M2;
    int measurements;

    MossbauerObjective(const int N, const int M1, const int M2, const int measurements) : N(N), M1(M1), M2(M2), measurements(measurements)
    {
      this->needInitialize = true;
      this->dim            = measurements;
      experiment           = make_shared<MossbauerExperiment>(measurements);
    }

    inline void Initialize() override
    {
      experiment->PrecomputeInputs(N);
    }

    inline double Evaluate(const VectorXd& input) override
    {
      experiment->ExpectedUtility(input, N, M1, M2)
    }
};

#endif // ifndef MossbauerObjective_h