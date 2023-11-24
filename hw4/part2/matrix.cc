#include <mpi.h>

#include <cstdio>
#include <string>

// Read size of matrix_a and matrix_b (n, m, l) and whole data of matrixes from
// stdin
//
// n_ptr:     pointer to n
// m_ptr:     pointer to m
// l_ptr:     pointer to l
// a_mat_ptr: pointer to matrix a (a should be a continuous memory space for
// placing n * m elements of int) b_mat_ptr: pointer to matrix b (b should be a
// continuous memory space for placing m * l elements of int)
void construct_matrices(int *n_ptr, int *m_ptr, int *l_ptr, int **a_mat_ptr,
                        int **b_mat_ptr)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank != 0)
    {
        *a_mat_ptr = nullptr;
        *b_mat_ptr = nullptr;
        return;
    }

    int n, m, l;
    int *a_mat, *b_mat;
    if (std::scanf("%d %d %d", &n, &m, &l) == EOF)
    {
        printf("Scanf error: cannot parse stdin\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    *n_ptr = n;
    *m_ptr = m;
    *l_ptr = l;

    a_mat = new int[n * m];
    b_mat = new int[m * l];

    *a_mat_ptr = a_mat;
    *b_mat_ptr = b_mat;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            if (std::scanf(" %d", &a_mat[i * m + j]) == EOF)
            {
                printf("Scanf error: cannot parse stdin\n");
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
    }

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < l; j++)
        {
            if (std::scanf(" %d", &b_mat[i * l + j]) == EOF)
            {
                printf("Scanf error: cannot parse stdin\n");
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
    }
}

// Just matrix multiplication (your should output the result in this function)
//
// n:     row number of matrix a
// m:     col number of matrix a / row number of matrix b
// l:     col number of matrix b
// a_mat: a continuous memory placing n * m elements of int
// b_mat: a continuous memory placing m * l elements of int
void matrix_multiply(const int n, const int m, const int l, const int *a_mat,
                     const int *b_mat)
{
    printf("Test\n");
}

// Remember to release your allocated memory
void destruct_matrices(int *a_mat, int *b_mat)
{
    delete[] a_mat;
    delete[] b_mat;
}