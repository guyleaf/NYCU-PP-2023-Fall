#include <atomic>
#include <cmath>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#include "utility.hpp"

using utility::ull;

const double CIRCLE_RADIUS = 1.0;

void toss_darts_in_circle(ull number_of_tosses,
                          std::atomic_ullong &number_of_hit)
{
    ull _number_of_hit = 0;
    while ((number_of_tosses--) > 0)
    {
        // toss dart randomly
        auto x = utility::random_uniform();
        auto y = utility::random_uniform();

        auto &&distance_to_center = std::pow(x, 2) + std::pow(y, 2);
        if (distance_to_center <= CIRCLE_RADIUS)
        {
            _number_of_hit++;
        }
    }

    number_of_hit += _number_of_hit;
}

void print_usage(const std::string program_name)
{
    std::cout << "Usage: " << program_name
              << " <Number of threads> <Number of tosses>" << std::endl;
}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        print_usage(*argv);
        return 1;
    }

    int number_of_threads = std::stoi(*(++argv));
    ull number_of_tosses = std::stoull(*(++argv));
    if (number_of_threads < 1 || number_of_tosses < 1)
    {
        std::cerr << "The number of threads should be greater than 0."
                  << std::endl;
        return 1;
    }

    std::atomic_ullong number_of_hit(0);
    std::vector<std::thread> threads(number_of_threads);

    // create threads
    ull number_of_tosses_per_thread = number_of_tosses / number_of_threads;
    for (int i = 0; i < number_of_threads - 1; i++)
    {
        threads[i] =
            std::thread(toss_darts_in_circle, number_of_tosses_per_thread,
                        std::ref(number_of_hit));
    }
    number_of_tosses_per_thread += number_of_tosses % number_of_threads;
    threads[number_of_threads - 1] =
        std::thread(toss_darts_in_circle, number_of_tosses_per_thread,
                    std::ref(number_of_hit));

    // join threads
    for (auto &thread : threads)
    {
        thread.join();
    }

    // print estimated of pi
    utility::write_double(
        std::cout, 4 * (number_of_hit / static_cast<double>(number_of_tosses)));
    std::cout << std::endl;
    return 0;
}