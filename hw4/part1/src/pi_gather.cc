#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#define FAST_RAND_MAX 0x7FFF

static unsigned int g_seed;

// Used to seed the generator.
inline void fast_srand(int seed) { g_seed = seed; }

// Compute a pseudorandom integer.
// Output value in range [0, 32767]
inline int fast_rand(void)
{
    g_seed = (214013 * g_seed + 2531011);
    return (g_seed >> 16) & FAST_RAND_MAX;
}

long long int estimate_tosses_in_circle(long long int tosses)
{
    long long int number_in_circle = 0;
    double x, y, distance_squared;

    while (tosses > 0)
    {
        x = (double)fast_rand() / FAST_RAND_MAX;
        y = (double)fast_rand() / FAST_RAND_MAX;
        distance_squared = x * x + y * y;
        if (distance_squared <= 1)
        {
            number_in_circle++;
        }
        tosses--;
    }

    return number_in_circle;
}

int main(int argc, char **argv)
{
    // --- DON'T TOUCH ---
    MPI_Init(&argc, &argv);
    double start_time = MPI_Wtime();
    double pi_result;
    long long int tosses = atoi(argv[1]);
    int world_rank, world_size;
    // ---

    // MPI init
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    fast_srand(world_rank);

    long long int tosses_per_process = tosses / world_size;
    if (world_rank == 0)
    {
        tosses_per_process += tosses % world_size;
    }

    long long int *partial_tosses = new long long int[world_size];
    partial_tosses[world_rank] = estimate_tosses_in_circle(tosses_per_process);

    // use MPI_Gather
    MPI_Gather(partial_tosses + world_rank, 1, MPI_LONG_LONG_INT,
               partial_tosses, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

    if (world_rank == 0)
    {
        // PI result
        tosses_per_process = 0;
        for (int i = 0; i < world_size; i++)
        {
            tosses_per_process += partial_tosses[i];
        }

        pi_result = 4 * (double)tosses_per_process / tosses;

        // --- DON'T TOUCH ---
        double end_time = MPI_Wtime();
        printf("%lf\n", pi_result);
        printf("MPI running time: %lf Seconds\n", end_time - start_time);
        // ---
    }

    delete[] partial_tosses;

    MPI_Finalize();
    return 0;
}
