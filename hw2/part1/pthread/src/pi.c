#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#include "utility.h"

#define CIRCLE_RADIUS 1.0

typedef struct thread_data
{
    ll number_of_tosses;
    ll total_of_tosses;
    double *pi;
    pthread_mutex_t *mutex;
} thread_data_t;

void *toss_darts_in_circle(void *args)
{
    unsigned int seed = time(NULL);
    thread_data_t *data = (thread_data_t *)args;
    ll number_of_tosses = data->number_of_tosses;

    ll *partial_pi = (ll *)calloc(1, sizeof(ll));
    double x, y;
    while (number_of_tosses > 0)
    {
        // toss dart randomly
        x = random_uniform_r(&seed);
        y = random_uniform_r(&seed);

        x = x * x + y * y;
        if (x <= CIRCLE_RADIUS)
        {
            (*partial_pi)++;
        }

        number_of_tosses--;
    }

    // partial_pi = 4 * partial_pi / data->total_of_tosses;

    // pthread_mutex_lock(data->mutex);
    // (*data->pi) += partial_pi;
    // pthread_mutex_unlock(data->mutex);

    return (void *)partial_pi;
}

void allocate_threads(pthread_t **threads, thread_data_t **data,
                      int number_of_threads)
{
    // dynamic allocate threads and data
    *threads = (pthread_t *)malloc(sizeof(pthread_t) * number_of_threads);
    if (threads == NULL)
    {
        perror("Failed to allocate memory for threads.");
        exit(1);
    }

    *data = (thread_data_t *)malloc(sizeof(thread_data_t) * number_of_threads);
    if (data == NULL)
    {
        perror("Failed to allocate memory for thread data.");
        exit(1);
    }
}

void free_threads(pthread_t **threads, thread_data_t **data)
{
    free(*threads);
    free(*data);
    *threads = NULL;
    *data = NULL;
}

void create_thread(pthread_t *thread, thread_data_t *data,
                   ll number_of_tosses_per_thread, ll total_of_tosses,
                   double *pi, pthread_mutex_t *mutex)
{
    data->number_of_tosses = number_of_tosses_per_thread;
    data->total_of_tosses = total_of_tosses;
    data->pi = pi;
    data->mutex = mutex;
    if (pthread_create(thread, NULL, toss_darts_in_circle, (void *)data) != 0)
    {
        perror("Failed to create thread.");
        exit(1);
    }
}

double estimate_pi(int number_of_threads, ll total_of_tosses)
{
    int i;
    double pi = 0;

    pthread_mutex_t mutex;
    pthread_mutex_init(&mutex, NULL);

    pthread_t *threads = NULL;
    thread_data_t *data = NULL;
    allocate_threads(&threads, &data, number_of_threads);

    // create threads
    ll number_of_tosses_per_thread = total_of_tosses / number_of_threads;
    for (i = 0; i < number_of_threads; i++)
    {
        if (i == number_of_threads - 1)
        {
            number_of_tosses_per_thread += total_of_tosses % number_of_threads;
        }
        create_thread(&threads[i], &data[i], number_of_tosses_per_thread,
                      total_of_tosses, &pi, &mutex);
    }

    // join threads
    ll *returnValue;
    ll count = 0;
    for (i = 0; i < number_of_threads; i++)
    {
        pthread_join(threads[i], (void **)&returnValue);
        count += *returnValue;
        free(returnValue);
    }

    free_threads(&threads, &data);
    pthread_mutex_destroy(&mutex);

    return 4 * count / (double)total_of_tosses;
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