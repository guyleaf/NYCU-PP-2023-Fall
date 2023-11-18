#include "bfs.h"

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include <cstddef>

#include "common/CycleTimer.h"
#include "common/graph.h"

#define ROOT_NODE_ID 0
#define NOT_VISITED_MARKER -1

#define ALPHA 14
#define BETA 24

struct _vertex_set
{
    // # of vertices in the set
    int count;

    // array of vertex ids in set
    int *vertices;
};

// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void top_down_step(Graph g, _vertex_set *frontier, int *distances,
                   int current_iter)
{
    int new_frontier_count = 0;

#pragma omp parallel for reduction(+ : new_frontier_count)
    for (int node = 0; node < g->num_nodes; node++)
    {
        if (frontier->vertices[node] == current_iter)
        {
            int start_edge = g->outgoing_starts[node];
            int end_edge = (node == g->num_nodes - 1)
                               ? g->num_edges
                               : g->outgoing_starts[node + 1];

            // attempt to add all neighbors to the new frontier
            for (int neighbor = start_edge; neighbor < end_edge; neighbor++)
            {
                int outgoing = g->outgoing_edges[neighbor];

                if (frontier->vertices[outgoing] == NOT_VISITED_MARKER)
                {
                    distances[outgoing] = distances[node] + 1;
                    frontier->vertices[outgoing] = current_iter + 1;
                    new_frontier_count++;
                }
            }
        }
    }

    frontier->count = new_frontier_count;
}

// Implements top-down BFS.
//
// Result of execution is that, for each node in the graph, the
// distance to the root is stored in sol.distances.
void bfs_top_down(Graph graph, solution *sol)
{
    _vertex_set frontier;
    frontier.vertices = (int *)malloc(graph->num_nodes * sizeof(int));
    frontier.count = 0;

    // initialize all nodes to NOT_VISITED
#pragma omp parallel for
    for (int i = 0; i < graph->num_nodes; i++)
    {
        frontier.vertices[i] = NOT_VISITED_MARKER;
        sol->distances[i] = NOT_VISITED_MARKER;
    }

    // setup frontier edges with the root node
    int current_iter = 0;
    frontier.vertices[ROOT_NODE_ID] = current_iter;
    frontier.count++;

    sol->distances[ROOT_NODE_ID] = 0;

    do
    {
#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        top_down_step(graph, &frontier, sol->distances, current_iter);

#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();
        printf("iter=%d %.4f sec\n", current_iter, end_time - start_time);
#endif

        current_iter++;
    } while (frontier.count != 0);

    free(frontier.vertices);
}

void bottom_up_step(Graph g, _vertex_set *frontier, int *distances,
                    int current_iter)
{
    int new_frontier_count = 0;

#pragma omp parallel for reduction(+ : new_frontier_count)
    for (int node = 0; node < g->num_nodes; node++)
    {
        if (frontier->vertices[node] == NOT_VISITED_MARKER)
        {
            int start_edge = g->incoming_starts[node];
            int end_edge = (node == g->num_nodes - 1)
                               ? g->num_edges
                               : g->incoming_starts[node + 1];

            // find neighbor which is in frontier & add to frontier
            for (int neighbor = start_edge; neighbor < end_edge; neighbor++)
            {
                int incoming = g->incoming_edges[neighbor];
                if (frontier->vertices[incoming] == current_iter)
                {
                    distances[node] = distances[incoming] + 1;
                    frontier->vertices[node] = current_iter + 1;
                    new_frontier_count++;
                    break;
                }
            }
        }
    }

    frontier->count = new_frontier_count;
}

void bfs_bottom_up(Graph graph, solution *sol)
{
    // For PP students:
    //
    // You will need to implement the "bottom up" BFS here as
    // described in the handout.
    //
    // As a result of your code's execution, sol.distances should be
    // correctly populated for all nodes in the graph.
    //
    // As was done in the top-down case, you may wish to organize your
    // code by creating subroutine bottom_up_step() that is called in
    // each step of the BFS process.

    _vertex_set frontier;
    frontier.vertices = (int *)malloc(graph->num_nodes * sizeof(int));
    frontier.count = 0;

    // initialize all nodes to NOT_VISITED
#pragma omp parallel for
    for (int i = 0; i < graph->num_nodes; i++)
    {
        frontier.vertices[i] = NOT_VISITED_MARKER;
        sol->distances[i] = NOT_VISITED_MARKER;
    }

    // setup frontier edges with the root node
    int current_iter = 0;
    frontier.vertices[ROOT_NODE_ID] = current_iter;
    frontier.count++;

    sol->distances[ROOT_NODE_ID] = 0;

    do
    {
#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        bottom_up_step(graph, &frontier, sol->distances, current_iter);

#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();
        printf("iter=%d %.4f sec\n", current_iter, end_time - start_time);
#endif

        current_iter++;
    } while (frontier.count != 0);

    free(frontier.vertices);
}

void bfs_hybrid(Graph graph, solution *sol)
{
    // For PP students:
    //
    // You will need to implement the "hybrid" BFS here as
    // described in the handout.

    _vertex_set frontier;
    frontier.vertices = (int *)malloc(graph->num_nodes * sizeof(int));
    frontier.count = 0;

    // initialize all nodes to NOT_VISITED
#pragma omp parallel for
    for (int i = 0; i < graph->num_nodes; i++)
    {
        frontier.vertices[i] = NOT_VISITED_MARKER;
        sol->distances[i] = NOT_VISITED_MARKER;
    }

    // setup frontier edges with the root node
    int current_iter = 0;
    frontier.vertices[ROOT_NODE_ID] = current_iter;
    frontier.count++;

    sol->distances[ROOT_NODE_ID] = 0;

    int size_of_graph = graph->num_nodes;
    bool choose_top_down = true;
    do
    {
#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif
        int frontier_size = frontier.count;

        if (choose_top_down)
        {
            top_down_step(graph, &frontier, sol->distances, current_iter);
        }
        else
        {
            bottom_up_step(graph, &frontier, sol->distances, current_iter);
        }

        // reference:
        // https://github.com/gbossi/Pagerank-DOBFS/blob/master/parallel_dobfs.cpp
        int new_frontier_size = frontier.count;
        float frontier_growing_ratio =
            (float)(new_frontier_size - frontier_size) / frontier_size;

        size_of_graph -= new_frontier_size;
        if (choose_top_down && frontier_growing_ratio > 0)
        {
            choose_top_down =
                (new_frontier_size * frontier_growing_ratio <=
                 (size_of_graph * frontier_growing_ratio) / ALPHA);
        }
        else if (!choose_top_down && frontier_growing_ratio < 0)
        {
            choose_top_down = (new_frontier_size < (float)size_of_graph / BETA);
        }

#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();
        printf("iter=%d %.4f sec\n", current_iter, end_time - start_time);
#endif

        current_iter++;
    } while (frontier.count != 0);

    free(frontier.vertices);
}