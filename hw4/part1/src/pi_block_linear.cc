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

    for (long long int toss = 0; toss < tosses; toss++)
    {
        x = (double)rand() / RAND_MAX;
        y = (double)rand() / RAND_MAX;
        distance_squared = x * x + y * y;
        if (distance_squared <= 1.0)
        {
            number_in_circle++;
        }
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

    // init MPI
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    srand(world_rank);

    long long int tosses_per_process = tosses / world_size;
    if (world_rank > 0)
    {
        // handle workers
        tosses_per_process = estimate_tosses_in_circle(tosses_per_process);
        MPI_Send(&tosses_per_process, 1, MPI_LONG_LONG_INT, 0, 0,
                 MPI_COMM_WORLD);
    }
    else if (world_rank == 0)
    {
        // master
        tosses_per_process =
            estimate_tosses_in_circle(tosses_per_process + tosses % world_size);
    }

    if (world_rank == 0)
    {
        MPI_Status status;

        // process PI result
        long long int partial_tosses;
        for (int rank = 1; rank < world_size; rank++)
        {
            MPI_Recv(&partial_tosses, 1, MPI_LONG_LONG_INT, rank, 0,
                     MPI_COMM_WORLD, &status);
            if (status.MPI_ERROR != MPI_SUCCESS)
            {
                printf("MPI_Recv error\n");
                MPI_Abort(MPI_COMM_WORLD, status.MPI_ERROR);
            }

            tosses_per_process += partial_tosses;
        }

        pi_result = 4 * tosses_per_process / (double)tosses;

        // --- DON'T TOUCH ---
        double end_time = MPI_Wtime();
        printf("%lf\n", pi_result);
        printf("MPI running time: %lf Seconds\n", end_time - start_time);
        // ---
    }

    MPI_Finalize();
    return 0;
}