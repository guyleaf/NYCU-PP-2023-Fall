#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#define MPI_MASTER 0
#define A_TAG 0
#define B_TAG 1
#define PANEL_WIDTH 32

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

void naive_matmul(const int m, const int n, const int k, const int *a_mat,
                  const int lda, const int *b_mat, const int ldb, int *c_mat,
                  const int ldc)
{
    // outer product
    for (int kk = 0; kk < k; kk++)
    {
        for (int ii = 0; ii < m; ii++)
        {
            for (int jj = 0; jj < n; jj++)
            {
                c_mat[jj * ldc + ii] +=
                    a_mat[kk * lda + ii] * b_mat[jj * ldb + kk];
            }
        }
    }
}

void print_matrix(const int m, const int n, const int *mat)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%d ", mat[j * m + i]);
        }
        printf("\n");
    }
}

void scatter_data(const int n, const int m, const int l, const int *a_mat,
                  const int *b_mat, int **m_a, int **n_a, int **m_b, int **n_b,
                  int **a_local, int **b_local)
{
    int *m_a_ = new int[comm_cart->rows];
    int *n_a_ = new int[comm_cart->cols];
    int *m_b_ = new int[comm_cart->rows];
    int *n_b_ = new int[comm_cart->cols];

    // matrix a / b format: m x n
    // 1. every process should know sizes of submatrix of a and b
    // 2. scatter submatrix a and b

    m_a_[comm_cart->row] = n;
    n_a_[comm_cart->col] = m;
    get_submatrix_size(m_a_[comm_cart->row], n_a_[comm_cart->col]);
    m_b_[comm_cart->row] = m;
    n_b_[comm_cart->col] = l;
    get_submatrix_size(m_b_[comm_cart->row], n_b_[comm_cart->col]);

    // 1. every process should know sizes of submatrix of a and b
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, n_a_, 1, MPI_INT,
                  comm_cart->row_comm);
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, n_b_, 1, MPI_INT,
                  comm_cart->row_comm);
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, m_a_, 1, MPI_INT,
                  comm_cart->col_comm);
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, m_b_, 1, MPI_INT,
                  comm_cart->col_comm);

    // 2. scatter submatrix a and b
    int size_a_local = m_a_[comm_cart->row] * n_a_[comm_cart->col];
    int size_b_local = m_b_[comm_cart->row] * n_b_[comm_cart->col];
    *a_local = new int[size_a_local];
    *b_local = new int[size_b_local];

    MPI_Request receive_requests[2];
    MPI_Irecv(*a_local, size_a_local, MPI_INT, 0, A_TAG, MPI_COMM_WORLD,
              receive_requests);
    MPI_Irecv(*b_local, size_b_local, MPI_INT, 0, B_TAG, MPI_COMM_WORLD,
              receive_requests + 1);

    if (comm_cart->global_rank == MPI_MASTER)
    {
        MPI_Datatype a_block;
        MPI_Datatype b_block;
        const int sizes_a[2] = {n, m}, sizes_b[2] = {m, l};
        int sizes_a_local[2], sizes_b_local[2];
        int starts_a[2] = {0, 0};
        int starts_b[2] = {0, 0};
        int coords[2];
        int row, col;

        // TODO: change to use scatter?
        MPI_Request *send_requests = new MPI_Request[2 * comm_cart->size];
        for (int rank = 0; rank < comm_cart->size; rank++)
        {
            MPI_Cart_coords(comm_cart->world, rank, 2, coords);
            row = coords[0];
            col = coords[1];
            sizes_a_local[0] = m_a_[row];
            sizes_b_local[0] = m_b_[row];
            sizes_a_local[1] = n_a_[col];
            sizes_b_local[1] = n_b_[col];

            MPI_Type_create_subarray(2, sizes_a, sizes_a_local, starts_a,
                                     MPI_ORDER_FORTRAN, MPI_INT, &a_block);
            MPI_Type_create_subarray(2, sizes_b, sizes_b_local, starts_b,
                                     MPI_ORDER_FORTRAN, MPI_INT, &b_block);
            MPI_Type_commit(&a_block);
            MPI_Type_commit(&b_block);

            MPI_Isend(a_mat, 1, a_block, rank, A_TAG, MPI_COMM_WORLD,
                      &send_requests[2 * rank]);
            MPI_Isend(b_mat, 1, b_block, rank, B_TAG, MPI_COMM_WORLD,
                      &send_requests[2 * rank + 1]);

            MPI_Type_free(&a_block);
            MPI_Type_free(&b_block);

            starts_a[1] = (starts_a[1] + n_a_[col]) % m;
            starts_b[1] = (starts_b[1] + n_b_[col]) % l;
            if ((rank + 1) % comm_cart->cols == 0)
            {
                starts_a[0] += m_a_[row];
                starts_b[0] += m_b_[row];
            }
        }
        MPI_Waitall(2 * comm_cart->size, send_requests, MPI_STATUSES_IGNORE);
        delete[] send_requests;
    }

    MPI_Waitall(2, receive_requests, MPI_STATUSES_IGNORE);

#ifdef DEBUG

    if (comm_cart->global_rank == 0)
    {
        printf("Scatter data\n");

        int coords[2];
        for (int rank = 0; rank < comm_cart->size; rank++)
        {
            MPI_Cart_coords(comm_cart->world, rank, 2, coords);
            printf("Rank %d, A: %d x %d, B: %d x %d\n", rank, m_a_[coords[0]],
                   n_a_[coords[1]], m_b_[coords[0]], n_b_[coords[1]]);
        }

        printf("Rank %d, A: \n", comm_cart->global_rank);
        print_matrix(m_a_[comm_cart->row], n_a_[comm_cart->col], *a_local);
        printf("Rank %d, B: \n", comm_cart->global_rank);
        print_matrix(m_b_[comm_cart->row], n_b_[comm_cart->col], *b_local);

        printf("=====================\n");
    }

    // printf("Rank %d, A: \n", comm_cart->global_rank);
    // print_matrix(m_a_[comm_cart->row], n_a_[comm_cart->col], *a_local);
    // printf("Rank %d, B: \n", comm_cart->global_rank);
    // print_matrix(m_b_[comm_cart->row], n_b_[comm_cart->col], *b_local);

#endif

    *m_a = m_a_;
    *n_a = n_a_;
    *m_b = m_b_;
    *n_b = n_b_;
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

        // stored in column-major order
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                if (std::scanf(" %d", &a_mat[j * n + i]) == EOF)
                {
                    printf("Scanf error: cannot parse stdin\n");
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }
        }

        // stored in column-major order
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < l; j++)
            {
                if (std::scanf(" %d", &b_mat[j * m + i]) == EOF)
                {
                    printf("Scanf error: cannot parse stdin\n");
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }
        }

        *a_mat_ptr = a_mat;
        *b_mat_ptr = b_mat;

#ifdef DEBUG
        printf("Read data\n");
        printf("A: %d x %d\n", n, m);
        print_matrix(n, m, *a_mat_ptr);
        printf("B: %d x %d\n", m, l);
        print_matrix(m, l, *b_mat_ptr);
        printf("=====================\n");
#endif
    }
    else
    {
        *a_mat_ptr = *b_mat_ptr = nullptr;
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

    int *m_a = nullptr;
    int *n_a = nullptr;
    int *m_b = nullptr;
    int *n_b = nullptr;
    int *a_local = nullptr;
    int *b_local = nullptr;
    int *c_local = nullptr;
    int *a_tmp = nullptr;
    int *b_tmp = nullptr;

    // matrix a / b format: m x n

    // 1. every process should know sizes of submatrix of a and b
    // 2. scatter submatrix a and b
    scatter_data(n, m, l, a_mat, b_mat, &m_a, &n_a, &m_b, &n_b, &a_local,
                 &b_local);

    // 3. initialize two tmp blocks for broadcasting and c submatrix
    int size_a_tmp = m_a[comm_cart->row] * PANEL_WIDTH;
    int size_b_tmp = PANEL_WIDTH * n_b[comm_cart->col];
    int size_c_local = m_a[comm_cart->row] * n_b[comm_cart->col];
    a_tmp = new int[size_a_tmp];
    b_tmp = new int[size_b_tmp];
    c_local = new int[size_c_local]();

    // 4. execute main loop
    int nb;
    int current_row = 0, current_col = 0;
    int ii = 0, jj = 0;
    int start;
    MPI_Datatype b_tmp_block;
    int sizes_b_tmp[] = {PANEL_WIDTH, n_b[comm_cart->col]};
    int sizes_b_tmp_local[2] = {0, n_b[comm_cart->col]};
    int starts[] = {0, 0};
    for (int kk = 0; kk < m; kk += nb)
    {
        nb = std::min(PANEL_WIDTH, m_b[current_row] - ii);
        nb = std::min(nb, n_a[current_col] - jj);

        if (comm_cart->col == current_col)
        {
            // thanks to the contiguous of array in column-major order
            // elements of column are stored first
            // so, we can copy the submatrix directly without stride
            start = jj * m_a[comm_cart->row];
            std::copy_n(a_local + start, m_a[comm_cart->row] * nb, a_tmp);

            // #ifdef DEBUG

            //             if (comm_cart->global_rank == 0)
            //             {
            //                 printf("A tmp\n");
            //                 printf("ii: %d, jj: %d\n", ii, jj);
            //                 printf("nb: %d, rest_n_a: %d, rest_m_b: %d\n",
            //                 nb,
            //                        n_a[current_col] - jj, m_b[current_row] -
            //                        ii);
            //                 printf("current_row: %d, current_col: %d\n",
            //                 current_row,
            //                        current_col);

            //                 printf("Rank %d, A tmp: \n",
            //                 comm_cart->global_rank);
            //                 print_matrix(m_a[comm_cart->row], PANEL_WIDTH,
            //                 a_tmp);

            //                 printf("=====================\n");
            //             }

            // #endif
        }
        if (comm_cart->row == current_row)
        {
            for (int i = ii; i < ii + nb; i++)
            {
                for (int j = 0; j < n_b[comm_cart->col]; j++)
                {
                    b_tmp[j * PANEL_WIDTH + i - ii] =
                        b_local[j * m_b[comm_cart->row] + i];
                }
            }

            // #ifdef DEBUG

            //             if (comm_cart->global_rank == 0)
            //             {
            //                 printf("B tmp\n");
            //                 printf("ii: %d, jj: %d\n", ii, jj);
            //                 printf("nb: %d, rest_n_a: %d, rest_m_b: %d\n",
            //                 nb,
            //                        n_a[current_col] - jj, m_b[current_row] -
            //                        ii);
            //                 printf("current_row: %d, current_col: %d\n",
            //                 current_row,
            //                        current_col);

            //                 printf("Rank %d, B tmp: \n",
            //                 comm_cart->global_rank);
            //                 print_matrix(PANEL_WIDTH, n_b[comm_cart->col],
            //                 b_tmp);

            //                 printf("=====================\n");
            //             }

            // #endif
        }

        MPI_Bcast(a_tmp, m_a[comm_cart->row] * nb, MPI_INT, current_col,
                  comm_cart->row_comm);

        sizes_b_tmp_local[0] = nb;
        MPI_Type_create_subarray(2, sizes_b_tmp, sizes_b_tmp_local, starts,
                                 MPI_ORDER_FORTRAN, MPI_INT, &b_tmp_block);
        MPI_Type_commit(&b_tmp_block);
        MPI_Bcast(b_tmp, 1, b_tmp_block, current_row, comm_cart->col_comm);
        MPI_Type_free(&b_tmp_block);

#ifdef DEBUG

        if (comm_cart->global_rank == 1)
        {
            printf("A tmp\n");
            printf("ii: %d, jj: %d\n", ii, jj);
            printf("nb: %d, rest_n_a: %d, rest_m_b: %d\n", nb,
                   n_a[current_col] - jj, m_b[current_row] - ii);
            printf("current_row: %d, current_col: %d\n", current_row,
                   current_col);

            printf("Rank %d, A tmp: \n", comm_cart->global_rank);
            print_matrix(m_a[comm_cart->row], PANEL_WIDTH, a_tmp);

            printf("=====================\n");

            printf("B tmp\n");
            printf("ii: %d, jj: %d\n", ii, jj);
            printf("nb: %d, rest_n_a: %d, rest_m_b: %d\n", nb,
                   n_a[current_col] - jj, m_b[current_row] - ii);
            printf("current_row: %d, current_col: %d\n", current_row,
                   current_col);

            printf("Rank %d, B tmp: \n", comm_cart->global_rank);
            print_matrix(PANEL_WIDTH, n_b[comm_cart->col], b_tmp);

            printf("=====================\n");
        }

#endif

        // perform matrix multiplication
        naive_matmul(m_a[comm_cart->row], n_b[comm_cart->col], nb, a_tmp,
                     m_a[comm_cart->row], b_tmp, PANEL_WIDTH, c_local,
                     m_a[comm_cart->row]);

        // matmul_naive(m_a[comm_cart->row], nb, n_b[comm_cart->col], a_tmp,
        // b_tmp,
        //              c_local);

#ifdef DEBUG

        if (comm_cart->global_rank == 0)
        {
            printf("C submatrix\n");

            printf("Rank %d, C: \n", comm_cart->global_rank);
            print_matrix(m_a[comm_cart->row], n_b[comm_cart->col], c_local);

            printf("=====================\n");
        }

#endif

        ii += nb;
        jj += nb;
        if (ii >= m_b[current_row])
        {
            ii = 0;
            current_row++;
        }
        if (jj >= n_a[current_col])
        {
            jj = 0;
            current_col++;
        }
    }

#ifdef DEBUG

    if (comm_cart->global_rank == 0)
    {
        printf("C submatrix\n");

        printf("Rank %d, C: \n", comm_cart->global_rank);
        print_matrix(m_a[comm_cart->row], n_b[comm_cart->col], c_local);

        printf("=====================\n");
    }

#endif

    // 5. gather and print matrix
    if (comm_cart->global_rank == MPI_MASTER)
    {
        MPI_Datatype c_block;
        const int sizes_c[] = {n, l};
        int sizes_c_local[2];
        int starts[] = {0, 0};
        int coords[2];

        // gather submatrices to master process
        MPI_Request *receive_requests = new MPI_Request[comm_cart->size];
        int *c_mat = new int[n * l];
        int row, col;
        for (int rank = 0; rank < comm_cart->size; rank++)
        {
            MPI_Cart_coords(comm_cart->world, rank, 2, coords);
            row = coords[0];
            col = coords[1];
            sizes_c_local[0] = m_a[row];
            sizes_c_local[1] = n_b[col];

            MPI_Type_create_subarray(2, sizes_c, sizes_c_local, starts,
                                     MPI_ORDER_FORTRAN, MPI_INT, &c_block);
            MPI_Type_commit(&c_block);
            MPI_Irecv(c_mat, 1, c_block, rank, 0, MPI_COMM_WORLD,
                      receive_requests + rank);
            MPI_Type_free(&c_block);

            starts[1] = (starts[1] + n_b[col]) % l;
            if ((rank + 1) % comm_cart->cols == 0)
            {
                starts[0] += m_a[row];
            }
        }

        MPI_Send(c_local, size_c_local, MPI_INT, 0, 0, MPI_COMM_WORLD);

        MPI_Waitall(comm_cart->size, receive_requests, MPI_STATUSES_IGNORE);

        // print matrix
        print_matrix(n, l, c_mat);

        delete[] c_mat;
        delete[] receive_requests;
    }
    else
    {
        MPI_Send(c_local, size_c_local, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    delete[] m_a;
    delete[] n_a;
    delete[] m_b;
    delete[] n_b;
    delete[] a_local;
    delete[] b_local;
    delete[] c_local;
    delete[] a_tmp;
    delete[] b_tmp;
}

// Remember to release your allocated memory
void destruct_matrices(int *a_mat, int *b_mat)
{
    delete[] a_mat;
    delete[] b_mat;
    delete comm_cart;
}