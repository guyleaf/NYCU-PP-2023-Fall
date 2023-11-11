#include "page_rank.h"

#include <omp.h>

#include <cmath>
#include <utility>

#include "../common/CycleTimer.h"
#include "../common/graph.h"

// pageRank --
//
// g:           graph to process (see common/graph.h)
// solution:    array of per-vertex vertex scores (length of array is
// num_nodes(g)) damping:     page-rank algorithm's damping parameter
// convergence: page-rank algorithm's convergence threshold
//
void pageRank(Graph g, double *solution, double damping, double convergence)
{
    const int number_of_nodes = num_nodes(g);
    const double inverse_number_of_nodes = 1.0 / number_of_nodes;
    const double random_jump_probability =
        (1 - damping) * inverse_number_of_nodes;

    double *score_old = new double[number_of_nodes];
    int *outgoing_sizes = new int[number_of_nodes];

    // initialize vertex weights to uniform probability. Double
    // precision scores are used to avoid underflow for large graphs

#pragma omp parallel for
    for (int i = 0; i < number_of_nodes; ++i)
    {
        solution[i] = inverse_number_of_nodes;
        outgoing_sizes[i] = outgoing_size(g, i);
    }

    /*
       For PP students: Implement the page rank algorithm here.  You
       are expected to parallelize the algorithm using openMP.  Your
       solution may need to allocate (and free) temporary arrays.

       Basic page rank pseudocode is provided below to get you started:

       // initialization: see example code above
       score_old[vi] = 1/numNodes;

       while (!converged) {

         // compute score_new[vi] for all nodes vi:
         score_new[vi] = sum over all nodes vj reachable from incoming edges
                            { score_old[vj] / number of edges leaving vj  }
         score_new[vi] = (damping * score_new[vi]) + (1.0-damping) / numNodes;

         score_new[vi] += sum over all nodes v in graph with no outgoing edges
                            { damping * score_old[v] / numNodes }

         // compute how much per-node scores have changed
         // quit once algorithm has converged

         global_diff = sum over all nodes vi { abs(score_new[vi] -
       score_old[vi]) }; converged = (global_diff < convergence)
       }

     */

    double sum, global_diff;
    bool converged = false;
    while (!converged)
    {
        global_diff = 0;
        sum = 0;

#pragma omp parallel for reduction(+ : sum)
        for (int i = 0; i < number_of_nodes; ++i)
        {
            // copy new scores to old scores for next iteration
            score_old[i] = solution[i];

            // sum over all nodes v in graph with no outgoing edges
            if (outgoing_sizes[i] == 0)
            {
                sum += solution[i];
            }
        }

        sum *= (damping * inverse_number_of_nodes);
        sum += random_jump_probability;

        // compute score_new[vi] for all nodes vi
#pragma omp parallel for schedule(static, 1) reduction(+ : global_diff)
        for (int i = 0; i < number_of_nodes; ++i)
        {
            // sum over all nodes vj reachable from incoming edges
            double incoming_sum = 0;
            for (const Vertex *j = incoming_begin(g, i);
                 j != incoming_end(g, i); ++j)
            {
                incoming_sum += score_old[*j] / outgoing_sizes[*j];
            }

            solution[i] = damping * incoming_sum + sum;
            global_diff += std::abs(solution[i] - score_old[i]);
        }

        converged = (global_diff < convergence);
    }

    delete[] outgoing_sizes;
    delete[] score_old;
}