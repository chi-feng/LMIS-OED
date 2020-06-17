#ifndef Timer_h
#define Timer_h

#include <tuple>
#include <chrono>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>

class Timer
{
  private:

    std::chrono::time_point<std::chrono::steady_clock> start;

  public:

    inline void Start() { start = std::chrono::steady_clock::now(); }

    Timer() { }

    inline double Elapsed()
    {
      return std::chrono::duration<double, std::ratio<1>>(std::chrono::steady_clock::now() - start).count();
    }

    inline double ETA(const int completed, const int total)
    {
      double rate = static_cast<double>(completed) / Elapsed();
      return (total - completed) / rate;
    }

    inline std::string ETAString(const int completed, const int total)
    {
      auto   eta = ETA(completed, total);
      int    h   = eta / 60 / 60;
      eta -= h * 60 * 60;
      int    m = eta / 60;
      eta -= m * 60;
      double s = eta;
      std::ostringstream stream;
      stream << std::setfill('0') << std::setw(2) << h << ':' << std::setfill('0') << std::setw(2) << m << ':' << std::setfill('0') << std::setw(2) << s;
      return stream.str();
    }
};

#endif // ifndef Timer_h
