#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

long long int estimate_tosses_in_circle(long long int tosses)
{
    long long int number_in_circle = 0;
    double x, y, distance_squared;

    while (tosses > 0)
    {
        x = (double)rand() / RAND_MAX;
        y = (double)rand() / RAND_MAX;
        distance_squared = x * x + y * y;
        if (distance_squared <= 1.0)
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

    MPI_Win win;

    // MPI init
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    srand(world_rank);

    long long int tosses_per_process = tosses / world_size;
    if (world_rank == 0)
    {
        tosses_per_process += tosses % world_size;
    }

    tosses_per_process = estimate_tosses_in_circle(tosses_per_process);

    long long int total_tosses = 0;
    if (world_rank == 0)
    {
        // Master
        MPI_Win_create(&total_tosses, sizeof(long long int),
                       sizeof(long long int), MPI_INFO_NULL, MPI_COMM_WORLD,
                       &win);
    }
    else
    {
        // Workers
        MPI_Win_create(NULL, 0, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    }

    MPI_Win_fence(0, win);
    MPI_Accumulate(&tosses_per_process, 1, MPI_LONG_LONG_INT, 0, 0, 1,
                   MPI_LONG_LONG_INT, MPI_SUM, win);
    MPI_Win_fence(0, win);

    MPI_Win_free(&win);

    if (world_rank == 0)
    {
        // handle PI result
        pi_result = 4 * total_tosses / (double)tosses;

        // --- DON'T TOUCH ---
        double end_time = MPI_Wtime();
        printf("%lf\n", pi_result);
        printf("MPI running time: %lf Seconds\n", end_time - start_time);
        // ---
    }

    MPI_Finalize();
    return 0;
}