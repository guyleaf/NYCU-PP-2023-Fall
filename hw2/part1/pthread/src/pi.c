#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#include "utility.h"

#define CIRCLE_RADIUS 1.0

typedef struct thread_data
{
    ll number_of_tosses;
    double pi_factor;
} thread_data_t;

void *toss_darts_in_circle(void *args)
{
    thread_data_t *data = (thread_data_t *)args;
    unsigned int seed = pthread_self();
    ll number_of_tosses = data->number_of_tosses;

    double x, y;
    double *partial_pi = (double *)calloc(1, sizeof(double));
    while (number_of_tosses > 0)
    {
        // toss dart randomly
        x = (double)rand_r(&seed) / RAND_MAX;
        y = (double)rand_r(&seed) / RAND_MAX;

        x = x * x + y * y;
        if (x <= CIRCLE_RADIUS)
        {
            (*partial_pi)++;
        }

        number_of_tosses--;
    }

    (*partial_pi) *= data->pi_factor;
    return (void *)partial_pi;
}

void allocate_threads(pthread_t **threads, int number_of_threads)
{
    // dynamic allocate threads and data
    *threads = (pthread_t *)malloc(sizeof(pthread_t) * number_of_threads);
    if (threads == NULL)
    {
        perror("Failed to allocate memory for threads.");
        exit(1);
    }
}

double estimate_pi(int number_of_threads, ll total_of_tosses)
{
    int i;
    double pi;
    double *result = NULL;

    thread_data_t data = {
        .number_of_tosses = total_of_tosses / number_of_threads,
        .pi_factor = 4.0 / total_of_tosses};
    thread_data_t main_thread_data = {
        .number_of_tosses =
            data.number_of_tosses + total_of_tosses % number_of_threads,
        .pi_factor = data.pi_factor};

    // allocate threads
    pthread_t *threads =
        (pthread_t *)malloc(sizeof(pthread_t) * number_of_threads - 1);
    if (threads == NULL)
    {
        perror("Failed to allocate memory for threads.");
        exit(1);
    }

    // create threads
    for (i = 0; i < number_of_threads - 1; i++)
    {
        if (pthread_create(&threads[i], NULL, toss_darts_in_circle,
                           (void *)&data) != 0)
        {
            perror("Failed to create thread.");
            exit(1);
        }
    }

    // main thread is also a worker
    result = (double *)toss_darts_in_circle((void *)&main_thread_data);
    pi = *result;
    free(result);

    // join threads
    for (i = 0; i < number_of_threads - 1; i++)
    {
        pthread_join(threads[i], (void **)&result);
        pi += *result;
        free(result);
    }

    free(threads);
    return pi;
}

void print_usage(const char *program_name)
{
    printf("Usage: %s <Number of threads> <Number of tosses>\n", program_name);
}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        print_usage(*argv);
        return 1;
    }

    int number_of_threads = atoi(*(++argv));
    ll total_of_tosses = atoll(*(++argv));

    if (number_of_threads < 1)
    {
        perror(
            "Invalid Arguments: the number of threads should be greater than "
            "0.");
        return 1;
    }
    if (total_of_tosses < 1)
    {
        perror(
            "Invalid Arguments: the number of tosses should be greater than "
            "0.");
        return 1;
    }

    double pi = estimate_pi(number_of_threads, total_of_tosses);
    printf("%f\n", pi);
    return 0;
}