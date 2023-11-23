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

    // array of vertices in set
    bool *vertices;
};

void vertex_set_init(_vertex_set *list, int max_count)
{
    list->count = 0;
    list->vertices = (bool *)calloc(max_count, sizeof(int));
}

// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void top_down_step(Graph g, _vertex_set *frontier, _vertex_set *new_frontier,
                   int *distances)
{
    // reset all vertices in the new frontier
    memset(new_frontier->vertices, 0, g->num_nodes * sizeof(bool));

    int new_frontier_count = 0;

#pragma omp parallel for schedule(static, 1) reduction(+ : new_frontier_count)
    for (int node = 0; node < g->num_nodes; node++)
    {
        if (frontier->vertices[node])
        {
            int start_edge = g->outgoing_starts[node];
            int end_edge = (node == g->num_nodes - 1)
                               ? g->num_edges
                               : g->outgoing_starts[node + 1];

            // attempt to add all neighbors to the new frontier
            for (int neighbor = start_edge; neighbor < end_edge; neighbor++)
            {
                int outgoing = g->outgoing_edges[neighbor];
                if (distances[outgoing] == NOT_VISITED_MARKER)
                {
                    distances[outgoing] = distances[node] + 1;
                    new_frontier->vertices[outgoing] = true;
                    // the outgoing node may be checked by multiple threads at
                    // the same time the approximation is not precise, but it is
                    // ok
                    new_frontier_count++;
                }
            }
        }
    }

    new_frontier->count = new_frontier_count;
}

// Implements top-down BFS.
//
// Result of execution is that, for each node in the graph, the
// distance to the root is stored in sol.distances.
void bfs_top_down(Graph graph, solution *sol)
{
    _vertex_set list1;
    _vertex_set list2;
    _vertex_set *frontier = &list1;
    _vertex_set *new_frontier = &list2;

    vertex_set_init(frontier, graph->num_nodes);
    vertex_set_init(new_frontier, graph->num_nodes);

    // initialize all nodes to NOT_VISITED
#pragma omp parallel for
    for (int i = 0; i < graph->num_nodes; i++)
    {
        sol->distances[i] = NOT_VISITED_MARKER;
    }

    // setup frontier edges with the root node
    frontier->vertices[ROOT_NODE_ID] = true;
    frontier->count++;

    sol->distances[ROOT_NODE_ID] = 0;

    do
    {
#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        top_down_step(graph, frontier, new_frontier, sol->distances);

#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();
        printf("iter=%d %.4f sec\n", current_level, end_time - start_time);
#endif

        // swap frontier and new frontier
        _vertex_set *tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    } while (frontier->count != 0);

    free(frontier->vertices);
    free(new_frontier->vertices);
}

void bottom_up_step(Graph g, _vertex_set *frontier, _vertex_set *new_frontier,
                    int *distances)
{
    // reset all vertices in the new frontier
    memset(new_frontier->vertices, 0, g->num_nodes * sizeof(bool));

    int new_frontier_count = 0;

#pragma omp parallel for reduction(+ : new_frontier_count)
    for (int node = 0; node < g->num_nodes; node++)
    {
        if (distances[node] == NOT_VISITED_MARKER)
        {
            int start_edge = g->incoming_starts[node];
            int end_edge = (node == g->num_nodes - 1)
                               ? g->num_edges
                               : g->incoming_starts[node + 1];

            // find neighbor which is in frontier & add node to frontier
            for (int neighbor = start_edge; neighbor < end_edge; neighbor++)
            {
                int incoming = g->incoming_edges[neighbor];
                if (frontier->vertices[incoming])
                {
                    distances[node] = distances[incoming] + 1;
                    new_frontier->vertices[node] = true;
                    new_frontier_count++;
                    break;
                }
            }
        }
    }

    new_frontier->count = new_frontier_count;
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

    _vertex_set list1;
    _vertex_set list2;
    _vertex_set *frontier = &list1;
    _vertex_set *new_frontier = &list2;

    vertex_set_init(frontier, graph->num_nodes);
    vertex_set_init(new_frontier, graph->num_nodes);

    // initialize all nodes to NOT_VISITED
#pragma omp parallel for
    for (int i = 0; i < graph->num_nodes; i++)
    {
        sol->distances[i] = NOT_VISITED_MARKER;
    }

    // setup frontier edges with the root node
    frontier->vertices[ROOT_NODE_ID] = true;
    frontier->count++;

    sol->distances[ROOT_NODE_ID] = 0;

    do
    {
#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        bottom_up_step(graph, frontier, new_frontier, sol->distances);

#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();
        printf("iter=%d %.4f sec\n", current_level, end_time - start_time);
#endif

        // swap frontier and new frontier
        _vertex_set *tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    } while (frontier->count != 0);

    free(frontier->vertices);
    free(new_frontier->vertices);
}

void bfs_hybrid(Graph graph, solution *sol)
{
    // For PP students:
    //
    // You will need to implement the "hybrid" BFS here as
    // described in the handout.

    _vertex_set list1;
    _vertex_set list2;
    _vertex_set *frontier = &list1;
    _vertex_set *new_frontier = &list2;

    vertex_set_init(frontier, graph->num_nodes);
    vertex_set_init(new_frontier, graph->num_nodes);

    // initialize all nodes to NOT_VISITED
#pragma omp parallel for
    for (int i = 0; i < graph->num_nodes; i++)
    {
        sol->distances[i] = NOT_VISITED_MARKER;
    }

    // setup frontier edges with the root node
    frontier->vertices[ROOT_NODE_ID] = true;
    frontier->count++;

    sol->distances[ROOT_NODE_ID] = 0;

    int size_of_graph = graph->num_nodes;
    bool choose_top_down = true;
    do
    {
#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        if (choose_top_down)
        {
            top_down_step(graph, frontier, new_frontier, sol->distances);
        }
        else
        {
            bottom_up_step(graph, frontier, new_frontier, sol->distances);
        }

        // reference:
        // https://github.com/gbossi/Pagerank-DOBFS/blob/master/parallel_dobfs.cpp
        // it has a little bit different implementation from the paper

        int frontier_size = frontier->count;
        int new_frontier_size = new_frontier->count;
        int frontier_growing_rate = new_frontier_size - frontier_size;

        // remove number of visited vertices from size of graph
        size_of_graph -= new_frontier_size;

        // if choose top down previsously, and the frontier growing rate is
        // increasing
        if (choose_top_down && frontier_growing_rate > 0)
        {
            choose_top_down =
                (new_frontier_size <= (float)size_of_graph / ALPHA);
        }
        // if choose bottom up previsously, and the frontier growing rate is
        // decreasing
        else if (!choose_top_down && frontier_growing_rate < 0)
        {
            // I think the repo is wrong. It should take new frontier size
            // instead of the old one
            choose_top_down = (new_frontier_size < (float)size_of_graph / BETA);
        }

#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();
        printf("iter=%d %.4f sec\n", current_level, end_time - start_time);
#endif

        // swap frontier and new frontier
        _vertex_set *tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    } while (frontier->count != 0);

    free(frontier->vertices);
    free(new_frontier->vertices);
}