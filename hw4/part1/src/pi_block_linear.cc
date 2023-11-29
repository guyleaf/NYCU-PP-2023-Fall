#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

double estimate_pi(long long int tosses, long long int total_tosses)
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

    return 4 * (double)number_in_circle / total_tosses;
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

    // init MPI
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    srand(world_rank);

    long long int tosses_per_process = tosses / world_size;
    long long int remainder = tosses % world_size;
    double partial_pi;

    if (world_rank > 0)
    {
        // handle workers
        if (world_rank <= remainder)
        {
            tosses_per_process++;
        }

        partial_pi = estimate_pi(tosses_per_process, tosses);
        MPI_Send(&partial_pi, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    else if (world_rank == 0)
    {
        // master
        partial_pi = estimate_pi(tosses_per_process, tosses);
    }

    if (world_rank == 0)
    {
        // process PI result
        pi_result = partial_pi;
        for (int rank = 1; rank < world_size; rank++)
        {
            MPI_Recv(&partial_pi, 1, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);

            pi_result += partial_pi;
        }

        // --- DON'T TOUCH ---
        double end_time = MPI_Wtime();
        printf("%lf\n", pi_result);
        printf("MPI running time: %lf Seconds\n", end_time - start_time);
        // ---
    }

    MPI_Finalize();
    return 0;
}