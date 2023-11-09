#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <functional>
#include <iomanip>
#include <limits>
#include <ostream>
#include <random>
#include <thread>

namespace utility
{

using ul = unsigned long;
using ull = unsigned long long;

/// @brief Sample a random double from [0,1), drawn from uniform distribution.
/// @return double-point value
double random_uniform()
{
    // seed the generator by thread id
    thread_local std::mt19937_64 generator(
        std::hash<std::thread::id>()(std::this_thread::get_id()));
    thread_local std::uniform_real_distribution<double> uniform_dist(0, 1);
    return uniform_dist(generator);
}

void write_double(std::ostream& out, double value)
{
    out << std::setprecision(std::numeric_limits<double>::digits10) << value;
}

}  // namespace utility

#endif