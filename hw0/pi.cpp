#include <iostream>
#include <cmath>
#include <iomanip>
#include <random>
#include <limits>

using ull = unsigned long long;

const ull MAX_TOSSES = (1ULL << 53) - 1;
const ull LEAST_TOSSES = 1e8;

int main()
{
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> uniformDist(-1, 1);

    double estimatedPi = 0;
    ull numberOfTosses = 0;
    ull numberInCircle = 0;
    while (numberOfTosses < MAX_TOSSES)
    {
        // toss dart randomly
        auto x = uniformDist(generator);
        auto y = uniformDist(generator);

        auto distanceToCenter = std::pow(x, 2) + std::pow(y, 2);
        if (distanceToCenter <= 1)
        {
            numberInCircle++;
        }

        numberOfTosses++;

        if (numberOfTosses >= LEAST_TOSSES)
        {
            // check if estimatedPi is accurate to two decimal places (3.14)
            estimatedPi = 4 * (numberInCircle / static_cast<double>(numberOfTosses));
            if (std::trunc(estimatedPi * 100) == 314)
            {
                break;
            }
        }
    }

    std::cout << std::setprecision(std::numeric_limits<double>::digits10) << estimatedPi << std::endl;
    return 0;
}
