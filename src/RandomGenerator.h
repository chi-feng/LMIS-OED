#ifndef RandomGenerator_h
#define RandomGenerator_h

#include <ctime>
#include <vector>
#include <random>

#include <Eigen/Core>

class RandomGenerator
{
  public:

    static std::mt19937 engine;
    static bool initialized;

    inline static void Initialize()
    {
      std::random_device random;

      RandomGenerator::engine      = std::mt19937(random());
      RandomGenerator::initialized = true;
    }

    inline static void SetSeed(size_t seed = std::mt19937::default_seed)
    {
      if (!RandomGenerator::initialized) {
        RandomGenerator::Initialize();
      }
      engine.seed(seed);
    }

    inline static double GetNormal()
    {
      std::normal_distribution<double> dist(0, 1);
      return dist(RandomGenerator::engine);
    }

    inline static double GetUniform()
    {
      std::uniform_real_distribution<double> dist(0, 1);
      return dist(RandomGenerator::engine);
    }

    inline static double GetExponential()
    {
      std::exponential_distribution<double> dist(1);
      return dist(RandomGenerator::engine);
    }

    inline static double GetGamma(const double alpha, const double beta)
    {
      std::gamma_distribution<double> dist(alpha, beta);
      return dist(RandomGenerator::engine);
    }

    inline static double GetChiSq(const double n)
    {
      std::chi_squared_distribution<double> dist(n);
      return dist(RandomGenerator::engine);
    }

    // helper functions below

    inline static double GetNormal(const double mu, const double sigma)
    {
      return RandomGenerator::GetNormal() * sigma + mu;
    }

    inline static double GetUniform(const double a, const double b)
    {
      return RandomGenerator::GetUniform() * (b - a) + a;
    }

    inline static Eigen::VectorXd GetBernoulliRandomVector(int const n, const double a, const double b)
    {
      Eigen::VectorXd result(n);

      for (int i = 0; i < n; ++i) {
        result(i) = RandomGenerator::GetUniform() < 0.5 ? a : b;
      }

      return result;
    }

    inline static Eigen::VectorXd GetNormalRandomVector(int const n)
    {
      Eigen::VectorXd result(n);

      for (int i = 0; i < n; ++i) {
        result(i) = RandomGenerator::GetNormal();
      }
      return result;
    }

    inline static Eigen::VectorXd GetUniformRandomVector(int const n, const double a, const double b)
    {
      Eigen::VectorXd result(n);

      for (int i = 0; i < n; ++i) {
        result(i) = RandomGenerator::GetUniform(a, b);
      }
      return result;
    }
};

#endif // ifndef RandomGenerator_h