#include <mpi.h>

#include <cmath>
#include <cstdio>
#include <string>

#define MPI_MASTER 0
#define A_TAG 0
#define B_TAG 1

struct mpi_comm_cart_t
{
    MPI_Comm world;
    MPI_Comm row_comm;
    MPI_Comm col_comm;
    int rows;
    int cols;
    int size;

    int row;
    int col;
    int global_rank;  // MPI_COMM_WORLD's rank
};
using mpi_comm_cart = mpi_comm_cart_t *;

mpi_comm_cart comm_cart = nullptr;

mpi_comm_cart initialize_cart_group()
{
    mpi_comm_cart _comm_cart = new mpi_comm_cart_t;

    MPI_Comm_rank(MPI_COMM_WORLD, &_comm_cart->global_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &_comm_cart->size);

    // initialize cart world
    int periods[2] = {0};
    int dims[2] = {0};
    MPI_Dims_create(_comm_cart->size, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &_comm_cart->world);
    _comm_cart->rows = dims[0];
    _comm_cart->cols = dims[1];

    int remain_dims[2];

    // create row comm for A
    remain_dims[0] = 0;
    remain_dims[1] = 1;
    MPI_Cart_sub(_comm_cart->world, remain_dims, &_comm_cart->row_comm);

    // create col comm for B
    remain_dims[0] = 1;
    remain_dims[1] = 0;
    MPI_Cart_sub(_comm_cart->world, remain_dims, &_comm_cart->col_comm);

    // determine my cart coords
    int coords[2];
    MPI_Cart_coords(_comm_cart->world, _comm_cart->global_rank, 2, coords);
    _comm_cart->row = coords[0];
    _comm_cart->col = coords[1];

    return _comm_cart;
}

void get_submatrix_size(int &m, int &n)
{
    int m_ = m / comm_cart->rows;
    int n_ = n / comm_cart->cols;
    if (comm_cart->row < (m % comm_cart->rows))
    {
        m_++;
    }
    if (comm_cart->col < (n % comm_cart->cols))
    {
        n_++;
    }
    m = m_;
    n = n_;
}

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
    int n, m, l;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == MPI_MASTER)
    {
        // master
        // read data from stdin
        if (std::scanf("%d %d %d", &n, &m, &l) == EOF)
        {
            printf("Scanf error: cannot parse stdin\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // broadcast n, m, l
    MPI_Request requests[3];
    MPI_Ibcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD, requests);
    MPI_Ibcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD, requests + 1);
    MPI_Ibcast(&l, 1, MPI_INT, 0, MPI_COMM_WORLD, requests + 2);

    if (rank == MPI_MASTER)
    {
        int *a_mat = new int[n * m];
        int *b_mat = new int[m * l];

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

        *a_mat_ptr = a_mat;
        *b_mat_ptr = b_mat;
    }

    MPI_Waitall(3, requests, MPI_STATUSES_IGNORE);

    *n_ptr = n;
    *m_ptr = m;
    *l_ptr = l;
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
    comm_cart = initialize_cart_group();

    int *m_a = new int[comm_cart->rows];
    int *n_a = new int[comm_cart->cols];
    int *m_b = new int[comm_cart->rows];
    int *n_b = new int[comm_cart->cols];

    // matrix a / b format: m x n
    // 1. every process should know sizes of submatrix of a and b
    // 2. scatter submatrix a and b

    m_a[comm_cart->row] = n;
    n_a[comm_cart->col] = m;
    get_submatrix_size(m_a[comm_cart->row], n_a[comm_cart->col]);
    m_b[comm_cart->row] = m;
    n_b[comm_cart->col] = l;
    get_submatrix_size(m_b[comm_cart->row], n_b[comm_cart->col]);

    // 1. every process should know sizes of submatrix of a and b
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, n_a, 1, MPI_INT,
                  comm_cart->row_comm);
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, n_b, 1, MPI_INT,
                  comm_cart->row_comm);
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, m_a, 1, MPI_INT,
                  comm_cart->col_comm);
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, m_b, 1, MPI_INT,
                  comm_cart->col_comm);

    // 2. scatter submatrix a and b
    int size_a_local = m_a[comm_cart->row] * n_a[comm_cart->col];
    int size_b_local = m_b[comm_cart->row] * n_b[comm_cart->col];
    int *a_local = new int[size_a_local];
    int *b_local = new int[size_b_local];

    MPI_Request receive_requests[2];
    MPI_Irecv(a_local, size_a_local, MPI_INT, 0, A_TAG, MPI_COMM_WORLD,
              receive_requests);
    MPI_Irecv(b_local, size_b_local, MPI_INT, 0, B_TAG, MPI_COMM_WORLD,
              receive_requests + 1);

    if (comm_cart->global_rank == MPI_MASTER)
    {
        MPI_Datatype a_block;
        MPI_Datatype b_block;
        const int sizes_a[2] = {n, m}, sizes_b[2] = {m, l};
        int sizes_a_local[2], sizes_b_local[2];
        int starts[2] = {0, 0};
        int coords[2];
        int row, col;
        int offset_a = 0;
        int offset_b = 0;

        // TODO: change to use scatter?
        MPI_Request *send_requests = new MPI_Request[2 * comm_cart->size];
        for (int rank = 0; rank < comm_cart->size; rank++)
        {
            MPI_Cart_coords(comm_cart->world, rank, 2, coords);
            row = coords[0];
            col = coords[1];
            sizes_a_local[0] = m_a[row];
            sizes_b_local[0] = m_b[row];
            sizes_a_local[1] = n_a[col];
            sizes_b_local[1] = n_b[col];

            MPI_Type_create_subarray(2, sizes_a, sizes_a_local, starts,
                                     MPI_ORDER_C, MPI_INT, &a_block);
            MPI_Type_create_subarray(2, sizes_b, sizes_b_local, starts,
                                     MPI_ORDER_C, MPI_INT, &b_block);
            MPI_Type_commit(&a_block);
            MPI_Type_commit(&b_block);

            MPI_Isend(&a_mat[offset_a], 1, a_block, rank, A_TAG, MPI_COMM_WORLD,
                      &send_requests[2 * rank]);
            MPI_Isend(&b_mat[offset_b], 1, b_block, rank, B_TAG, MPI_COMM_WORLD,
                      &send_requests[2 * rank + 1]);

            MPI_Type_free(&a_block);
            MPI_Type_free(&b_block);

            offset_a += n_a[col];
            offset_b += n_b[col];
            if ((rank + 1) % comm_cart->cols == 0)
            {
                offset_a += std::max(0, m_a[row] - 1) * m;
                offset_b += std::max(0, m_b[row] - 1) * l;
            }
        }
        MPI_Waitall(2 * comm_cart->size, send_requests, MPI_STATUSES_IGNORE);
        delete[] send_requests;
    }

    MPI_Waitall(2, receive_requests, MPI_STATUSES_IGNORE);

    if (comm_cart->global_rank == MPI_MASTER)
    {
        for (int i = 0; i < size_a_local; i++)
        {
            if (i % n_a[comm_cart->col] == 0)
            {
                printf("\n");
            }
            printf("%d ", a_local[i]);
        }
    }

    delete[] m_a;
    delete[] n_a;
    delete[] m_b;
    delete[] n_b;
    delete[] a_local;
    delete[] b_local;
}

// Remember to release your allocated memory
void destruct_matrices(int *a_mat, int *b_mat)
{
    delete[] a_mat;
    delete[] b_mat;
    delete comm_cart;
}