#include "bfs.h"

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <climits>
#include <cmath>
#include <cstddef>

#include "common/CycleTimer.h"
#include "common/graph.h"

#define ROOT_NODE_ID 0
#define NOT_VISITED_MARKER -1
#define NO_COMING_EDGE_MARKER -2

void vertex_set_clear(vertex_set *list) { list->count = 0; }

void vertex_set_init(vertex_set *list, int count)
{
    list->max_vertices = count;
    list->vertices = (int *)calloc(list->max_vertices, sizeof(int));
    vertex_set_clear(list);
}

void vertex_set_deinit(vertex_set *list)
{
    free(list->vertices);
    list->vertices = NULL;
}

// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void top_down_step(Graph g, vertex_set *frontier_edges,
                   vertex_set *new_frontier_edges, int *distances)
{
    bool has_new_frontier = false;

#pragma omp parallel for
    for (int node = 0; node < g->num_nodes; node++)
    {
        int index = 2 * node;
        int *start_edge = &frontier_edges->vertices[index];
        int *end_edge = &frontier_edges->vertices[index + 1];

        // attempt to add all neighbors to the new frontier
        for (int neighbor = *start_edge; neighbor < *end_edge; neighbor++)
        {
            int outgoing = g->outgoing_edges[neighbor];

            if (distances[outgoing] == NOT_VISITED_MARKER)
            {
                distances[outgoing] = distances[node] + 1;

                index = 2 * outgoing;
                new_frontier_edges->vertices[index] =
                    g->outgoing_starts[outgoing];
                new_frontier_edges->vertices[index + 1] =
                    (outgoing == g->num_nodes - 1)
                        ? g->num_edges
                        : g->outgoing_starts[outgoing + 1];
                has_new_frontier = true;
            }
        }

        // reset frontier edges
        *start_edge = *end_edge = 0;
    }

    new_frontier_edges->count = has_new_frontier ? 1 : 0;
}

// Implements top-down BFS.
//
// Result of execution is that, for each node in the graph, the
// distance to the root is stored in sol.distances.
void bfs_top_down(Graph graph, solution *sol)
{
    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, 2 * graph->num_nodes);
    vertex_set_init(&list2, 2 * graph->num_nodes);

    vertex_set *frontier_edges = &list1;
    vertex_set *new_frontier_edges = &list2;

    // initialize all nodes to NOT_VISITED
#pragma omp parallel for
    for (int i = 0; i < graph->num_nodes; i++)
    {
        sol->distances[i] = NOT_VISITED_MARKER;
    }

    // setup frontier edges with the root node
    frontier_edges->vertices[2 * ROOT_NODE_ID] =
        graph->outgoing_starts[ROOT_NODE_ID];
    frontier_edges->vertices[2 * ROOT_NODE_ID + 1] =
        (ROOT_NODE_ID == graph->num_nodes - 1)
            ? graph->num_edges
            : graph->outgoing_starts[ROOT_NODE_ID + 1];
    frontier_edges->count++;

    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier_edges->count != 0)
    {
#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(frontier_edges);

        top_down_step(graph, frontier_edges, new_frontier_edges,
                      sol->distances);

#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();
        printf("iter=%d %.4f sec\n", current_iter, end_time - start_time);
#endif

        // swap pointers
        vertex_set *tmp = frontier_edges;
        frontier_edges = new_frontier_edges;
        new_frontier_edges = tmp;
    }

    vertex_set_deinit(frontier_edges);
    vertex_set_deinit(new_frontier_edges);
}

void bottom_up_step(Graph g, vertex_set *frontier, int *distances,
                    int current_iter)
{
    bool has_new_frontier = false;

#pragma omp parallel for
    for (int node = 0; node < g->num_nodes; node++)
    {
        if (distances[node] == NOT_VISITED_MARKER)
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
                    has_new_frontier = true;
                    break;
                }
            }

            // mark as no coming edge
            if (start_edge == end_edge)
            {
                distances[node] = NO_COMING_EDGE_MARKER;
            }
        }
    }

    frontier->count = has_new_frontier ? 1 : 0;
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

    vertex_set frontier;
    vertex_set_init(&frontier, graph->num_nodes);

    // initialize all nodes to NOT_VISITED
#pragma omp parallel for
    for (int i = 0; i < graph->num_nodes; i++)
    {
        sol->distances[i] = NOT_VISITED_MARKER;
        frontier.vertices[i] = NOT_VISITED_MARKER;
    }

    // setup frontier with the root node
    int current_iter = ROOT_NODE_ID;
    frontier.vertices[ROOT_NODE_ID] = current_iter;
    frontier.count++;

    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier.count != 0)
    {
#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(&frontier);

        bottom_up_step(graph, &frontier, sol->distances, current_iter);

#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();
        printf("iter=%d %.4f sec\n", current_iter, end_time - start_time);
#endif

        current_iter++;
    }

    // restore macro values other than NOT_VISITED_MARKER
#pragma omp parallel for
    for (int i = 0; i < graph->num_nodes; i++)
    {
        sol->distances[i] = std::max(sol->distances[i], NOT_VISITED_MARKER);
    }

    vertex_set_deinit(&frontier);
}

void bfs_hybrid(Graph graph, solution *sol)
{
    // For PP students:
    //
    // You will need to implement the "hybrid" BFS here as
    // described in the handout.
}