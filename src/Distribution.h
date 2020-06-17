#ifndef Distribution_h
#define Distribution_h

#include <cmath>
#include <Eigen/Core>

class Distribution
{
  protected:

    int dim;

  public:

    virtual Eigen::VectorXd GetSample()                               = 0;
    virtual Eigen::VectorXd GetSamples(const int n)                   = 0;
    virtual double          LogDensity(const Eigen::VectorXd& sample) = 0;
    inline double           Density(const Eigen::VectorXd& sample)
    {
      return std::exp(LogDensity(sample));
    }

    inline int GetDimension() { return dim; }
};

#endif // ifndef Distribution_h
