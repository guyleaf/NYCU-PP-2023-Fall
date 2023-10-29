#include <cmath>
#include <iostream>
#include <string>
#include <thread>
#include <utility.hpp>
#include <vector>

using utility::ull;

void toss_darts_in_circle(ull number_of_tosses, double radius = 1.0)
{
    ull number_of_hit = 0;
    while ((number_of_tosses--) > 0)
    {
        // toss dart randomly
        auto x = utility::random_uniform();
        auto y = utility::random_uniform();

        auto &&distance_to_center = std::pow(x, 2) + std::pow(y, 2);
        if (distance_to_center <= radius)
        {
            number_of_hit++;
        }
    }
}

void print_usage(const std::string program_name)
{
    std::cout << "Usage: " << program_name
              << "<Number of threads> <Number of tosses>" << std::endl;
}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        print_usage(*argv);
        return 1;
    }

    long number_of_threads = std::stol(*(++argv));
    ull number_of_tosses = std::stoull(*(++argv));
    if (number_of_threads < 1 || number_of_tosses < 1)
    {
        std::cerr << "The number of threads should be greater than 0."
                  << std::endl;
        return 1;
    }

    std::vector<std::thread> threads;

    double estimatedPi = 0;
    utility::write_double(std::cout, estimatedPi);
    std::cout << std::endl;
    return 0;
}