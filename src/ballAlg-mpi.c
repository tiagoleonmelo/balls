#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include "gen_points.c"
#include <stddef.h>

typedef struct node
{
    int id;
    long subset_len;
    long *subset;
    double *center_coord;
    double radius;
    struct node *L;
    struct node *R;
} node_t;

typedef struct dbl_index_struct
{
    double val;
    long index;
} dbl_index;

int n_dims;
long n_points;
int seed;
double **pts;
node_t **nodes;
// node_t *_nodes;
int pid;
int p;
MPI_Datatype mpi_dbl_index;
MPI_Op mpi_maxloc_dbl_index;

void destroy_memory()
{
    free(nodes);
    // free(_nodes);
    free(pts[0]);
    free(pts);
}

void fill_node(int id, double *center_coord, double radius)
{
    // Retrieve node pointer
    node_t *node = nodes[id];

    node->id = id;
    node->center_coord = center_coord;
    node->radius = radius;

    node->L = NULL;
    node->R = NULL;
}

void maxloc_dbl_index(void *in, void *inout, int *len, MPI_Datatype *type)
{
    /* ignore type, just trust that it's our dbl_index type */
    dbl_index *invals = in;
    dbl_index *inoutvals = inout;

    for (int i = 0; i < *len; i++)
    {
        if (invals[i].val > inoutvals[i].val)
        {
            inoutvals[i].val = invals[i].val;
            inoutvals[i].index = invals[i].index;
        }
    }

    return;
}

double distance(long a, long b)
{

    double total = 0;
    double diff;
    double *pt_a = pts[a];
    double *pt_b = pts[b];

    for (int i = 0; i < n_dims; i++)
    {
        diff = pt_a[i] - pt_b[i];
        total += (diff * diff);
    }

    return sqrt(total);
}

long *furthest_apart(long *subset, long subset_len)
{

    // Don't forget: A and B are indexes relative to the subset!
    long first = subset[0];

    MPI_Bcast(pts[first], n_dims, MPI_DOUBLE, first % p, MPI_COMM_WORLD);


    double max_dist = -1;
    double curr_dist;

    long *ret = (long *)malloc(sizeof(long) * 2);
    if (ret == NULL)
        return NULL;

    for (long i = 1; i < subset_len; i++)
    {
        if (subset[i] % p == pid)
        {
            curr_dist = distance(first, subset[i]);

            if (curr_dist > max_dist)
            {
                max_dist = curr_dist;
                ret[0] = i;
            }
        }
    }


    dbl_index local, global;
    local.index = ret[0];
    local.val = max_dist;

    MPI_Allreduce(&local, &global, 1, mpi_dbl_index, mpi_maxloc_dbl_index, MPI_COMM_WORLD);

    printf("%f %f\n", max_dist, global.val);

    // root = index % num_p

    //MPI_Bcast(buf, len, type, root, comm);


    max_dist = -1;

    for (long i = 0; i < subset_len; i++)
    {
        if (ret[0] == i)
        {
            continue;
        }

        curr_dist = distance(subset[ret[0]], subset[i]);

        if (curr_dist > max_dist)
        {
            max_dist = curr_dist;
            ret[1] = i;
        }
    }

    return ret;
}

double inner_product(double *a, double *b)
{
    double total = 0;

    for (int dim = 0; dim < n_dims; dim++)
    {
        total += a[dim] * b[dim];
    }

    return total;
}

void difference(double *a, double *b, double *res)
{
    for (int dim = 0; dim < n_dims; dim++)
    {
        res[dim] = a[dim] - b[dim];
    }
}

double *orth_projection(long *subset, long subset_len, long a, long b, double *pt_a, double *b_minus_a_vec)
{

    double *proj = (double *)malloc(sizeof(double) * subset_len);
    if (proj == NULL)
        return NULL;

    double b_minus_a = b_minus_a_vec[0];

    double inner_prod_b_minus_a = inner_product(b_minus_a_vec, b_minus_a_vec);

    for (long p = 0; p < subset_len; p++)
    {

        if (p == a || p == b)
        {
            proj[p] = pts[subset[p]][0];
            continue;
        }

        double *p_minus_a_vec = (double *)malloc(sizeof(double) * n_dims);
        if (p_minus_a_vec == NULL)
        {
            free(proj);
            return NULL;
        }

        difference(pts[subset[p]], pt_a, p_minus_a_vec);

        double scalar = inner_product(p_minus_a_vec, b_minus_a_vec) / inner_prod_b_minus_a;

        proj[p] = scalar * b_minus_a + pt_a[0];

        free(p_minus_a_vec);
    }

    return proj;
}

double *single_projection(double orth_zero, double *pt_a, long *subset, double *b_minus_a_vec)
{
    double *new_point = (double *)malloc(sizeof(double) * n_dims);
    if (new_point == NULL)
        return NULL;

    // retrieve scalar knowing that orth_zero = scalar * (b-a)[0] + a[0]
    double scalar = (orth_zero - pt_a[0]) / b_minus_a_vec[0];

    new_point[0] = orth_zero;
    for (int i = 1; i < n_dims; i++)
    {
        new_point[i] = scalar * b_minus_a_vec[i] + pt_a[i];
    }

    return new_point;
}

int cmpfunc(const void *pt_1, const void *pt_2)
{
    double const p1 = *(double const *)pt_1;
    double const p2 = *(double const *)pt_2;

    if (p1 < p2)
        return -1;
    if (p1 > p2)
        return 1;
    return 0;
}

void swap(double *a, double *b)
{
    double temp = *a;
    *a = *b;
    *b = temp;
}

// Standard partition process of QuickSort.
long partition(double arr[], long l, long r)
{
    double x = arr[r];
    long i = l;
    for (long j = l; j <= r - 1; j++)
    {
        if (arr[j] <= x)
        {
            swap(&arr[i], &arr[j]);
            i++;
        }
    }
    swap(&arr[i], &arr[r]);
    return i;
}

// Picks a random pivot element between l and r and partitions
// arr[l..r] around the randomly picked element using partition()
long randomPartition(double arr[], long l, long r)
{
    long n = r - l + 1;
    long pivot = rand() % n;
    swap(&arr[l + pivot], &arr[r]);
    return partition(arr, l, r);
}

// Quickselect implementation that returns the smallest k'th element
double quickselect(double arr[], long l, long r, long k)
{
    // If k is smaller than number of elements in array
    if (k > 0 && k <= r - l + 1)
    {
        // Partition the array around a random element and
        // get position of pivot element in sorted array
        long pos = randomPartition(arr, l, r);

        // If position is same as k
        if (pos - l == k - 1)
            return arr[pos];
        if (pos - l > k - 1) // If position is more, recur for left subarray
            return quickselect(arr, l, pos - 1, k);

        // Else recur for right subarray
        return quickselect(arr, pos + 1, r, k - pos + l - 1);
    }

    // If k is more than the number of elements in the array
    return 0.0;
}

double *find_median(double *orth, long *subset, long subset_len, double *pt_a, double *b_minus_a_vec)
{
    double *new_pt;

    double *ordered_orth = (double *)malloc(sizeof(double) * subset_len);
    if (ordered_orth == NULL)
        return NULL;

    memcpy(ordered_orth, orth, sizeof(double) * subset_len);

    long mid_index = subset_len / 2;

    double median = quickselect(ordered_orth, 0, subset_len - 1, mid_index + 1);

    new_pt = single_projection(median, pt_a, subset, b_minus_a_vec); // Pass along point
    if (new_pt == NULL)
    {
        free(ordered_orth);
        return NULL;
    }

    if (subset_len % 2 == 0)
    {
        double *pt2;
        long mid_index2 = (subset_len / 2) - 1;
        double median2 = quickselect(ordered_orth, 0, subset_len - 1, mid_index2 + 1);
        pt2 = single_projection(median2, pt_a, subset, b_minus_a_vec);
        if (pt2 == NULL)
        {
            free(ordered_orth);
            free(new_pt);
            return NULL;
        }

        for (int i = 0; i < n_dims; i++)
        {
            new_pt[i] = (new_pt[i] + pt2[i]) / 2.0;
        }

        free(pt2);
    }

    free(ordered_orth);

    return new_pt;
}

void split(double *orth, double *median, long *subset, long subset_len, long *left, long *right)
{
    long ind_left = 0;
    long ind_right = 0;

    for (int i = 0; i < subset_len; i++)
    {
        if (median[0] > orth[i])
        {
            left[ind_left] = subset[i];
            ind_left++;
        }
        else
        {
            right[ind_right] = subset[i];
            ind_right++;
        }
    }
}

double find_radius(double *center, long *subset, long subset_len)
{

    double current_dist;
    double max_dist = -1;
    double total;
    double diff;
    double *pt;

    for (long i = 0; i < subset_len; i++)
    {
        current_dist = 0;
        total = 0;
        pt = pts[subset[i]];
        for (int d = 0; d < n_dims; d++)
        {
            diff = center[d] - pt[d];
            total += (diff * diff);
        }
        current_dist = sqrt(total);

        if (current_dist > max_dist)
        {
            max_dist = current_dist;
        }
    }

    return max_dist;
}

void build_tree(long id, int recursive)
{
    node_t *root = nodes[id];
    long subset_len = root->subset_len;
    long *subset = root->subset;

    if (subset_len > 1)
    {
        long *a_b;
        double *pt_a;
        double *pt_b;
        double *orth;
        double *median;
        double radius;

        long *subset_L = (long *)malloc(sizeof(long) * subset_len / 2);
        if (subset_L == NULL)
        {
            perror("Error allocating memory for subset_L");
            destroy_memory();
            exit(1);
        }
        long *subset_R = (long *)malloc(sizeof(long) * (subset_len / 2 + (subset_len % 2)));
        if (subset_R == NULL)
        {
            perror("Error allocating memory for subset_R");
            free(subset_L);
            destroy_memory();
            exit(1);
        }

        if (subset_len == 2)
        {
            pt_a = pts[subset[0]];
            pt_b = pts[subset[1]];

            orth = (double *)malloc(sizeof(double) * subset_len);
            if (orth == NULL)
            {
                perror("Error allocating memory for orth");
                free(subset_L);
                free(subset_R);
                destroy_memory();
                exit(1);
            }

            orth[0] = pt_a[0];
            orth[1] = pt_b[0];
            double total;
            total = 0;

            median = (double *)malloc(sizeof(double) * n_dims);
            if (median == NULL)
            {
                perror("Error allocating memory for median");
                free(subset_L);
                free(subset_R);
                free(orth);
                destroy_memory();
                exit(1);
            }

            double diff;
            for (int i = 0; i < n_dims; i++)
            {
                median[i] = (pt_a[i] + pt_b[i]) / 2.0;
                diff = median[i] - pt_a[i];
                total += (diff * diff);
            }
            radius = sqrt(total);

            // Create node with id, center coords and radius
            fill_node(id, median, radius);

            split(orth, median, subset, subset_len, subset_L, subset_R);
        }
        else
        {
            // Find A and B
            a_b = furthest_apart(subset, subset_len);
            if (a_b == NULL)
            {
                perror("Error allocating memory in furthest_apart");
                free(subset_L);
                free(subset_R);
                destroy_memory();
                exit(1);
            }

            pt_a = pts[subset[a_b[0]]];
            pt_b = pts[subset[a_b[1]]];

            double *b_minus_a_vec = (double *)malloc(sizeof(double) * n_dims);
            if (b_minus_a_vec == NULL)
            {
                perror("Error allocating memory for b_minus_a_vec");
                free(subset_L);
                free(subset_R);
                free(a_b);
                destroy_memory();
                exit(1);
            }

            difference(pt_a, pt_b, b_minus_a_vec);

            // Orthogonal projection
            orth = orth_projection(subset, subset_len, a_b[0], a_b[1], pt_a, b_minus_a_vec);
            if (orth == NULL)
            {
                perror("Error allocating memory in orth_projection");
                free(subset_L);
                free(subset_R);
                free(a_b);
                free(b_minus_a_vec);
                destroy_memory();
                exit(1);
            }

            // Find median point
            median = find_median(orth, subset, subset_len, pt_a, b_minus_a_vec);
            if (median == NULL)
            {
                perror("Error allocating memory in find_median");
                free(subset_L);
                free(subset_R);
                free(a_b);
                free(b_minus_a_vec);
                free(orth);
                destroy_memory();
                exit(1);
            }

            // Find radius
            radius = find_radius(median, subset, subset_len);

            // Split between Left and Right branches
            split(orth, median, subset, subset_len, subset_L, subset_R);

            // Create node with id, center coords and radius
            fill_node(id, median, radius);

            free(b_minus_a_vec);
            free(a_b);
        }

        long l_id = id + 1;
        long r_id = id + subset_len - (subset_len % 2);

        node_t *left = nodes[l_id];
        root->L = left;
        left->id = l_id;
        left->subset_len = subset_len / 2;
        left->subset = subset_L;

        node_t *right = nodes[r_id];
        root->R = right;
        right->id = r_id;
        right->subset_len = subset_len / 2 + (subset_len % 2);
        right->subset = subset_R;

        if (recursive)
        {
            build_tree(l_id, 1);
            build_tree(r_id, 1);
        }

        free(orth);
    }
    else
    {
        // Create leaf
        fill_node(id, pts[subset[0]], 0);
    }

    free(subset);
}

void dump_tree(node_t *root)
{
    if (root == NULL)
        return;

    node_t *L = root->L;
    node_t *R = root->R;

    if (!(L == NULL && R == NULL))
        free(root->center_coord);

    // Recursive call
    dump_tree(L);
    dump_tree(R);
}

void print_tree(node_t *root)
{
    if (root == NULL)
        return;

    node_t *L = root->L;
    node_t *R = root->R;

    if (!(L == NULL && R == NULL))
    {
        printf("%d %d %d %f ", root->id, L->id, R->id, root->radius);
    }
    else
    {
        printf("%d -1 -1 %f ", root->id, root->radius);
    }

    for (size_t i = 0; i < n_dims; i++)
    {
        printf("%.6f ", root->center_coord[i]);
    }

    printf("\n");

    // Recursive call
    print_tree(L);
    print_tree(R);
}

void print_node(int id)
{
    node_t *root = nodes[id];

    node_t *L = root->L;
    node_t *R = root->R;

    printf("%d %d %d %f ", root->id, L->id, R->id, root->radius);

    for (size_t i = 0; i < n_dims; i++)
    {
        printf("%.6f ", root->center_coord[i]);
    }

    printf("\n");
}

int get_id(int i, int level)
{
    if (level == 0)
        return 0;

    int parent_id = get_id(i / 2, level - 1);
    if (i % 2 == 0)
    {
        return parent_id + 1;
    }
    else
    {
        long len = nodes[parent_id]->subset_len;
        return parent_id + len - (len % 2);
    }
}

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        printf("Usage: %s <n_dims> <n_points> <seed>\n", argv[0]);
        exit(1);
    }

    // MPI_Status status;

    MPI_Init(&argc, &argv);

    n_dims = atoi(argv[1]);
    n_points = atol(argv[2]);
    seed = atoi(argv[3]);

    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    /* create our new data type */
    MPI_Datatype types[2] = {MPI_DOUBLE, MPI_LONG};
    MPI_Aint disps[2] = {offsetof(dbl_index, val),
                         offsetof(dbl_index, index)};

    int lens[2] = {1, 1};
    MPI_Type_create_struct(2, lens, disps, types, &mpi_dbl_index);
    MPI_Type_commit(&mpi_dbl_index);

    /* create our operator */
    MPI_Op_create(maxloc_dbl_index, 1, &mpi_maxloc_dbl_index);

    double exec_time;
    exec_time = -MPI_Wtime();

    pts = get_points(argc, argv, &n_dims, &n_points, p, pid);
    long *full_set = (long *)malloc(sizeof(long) * n_points);
    if (full_set == NULL)
    {
        perror("Error allocating memory for full_set");
        exit(1);
    }

    for (long i = 0; i < n_points; i++)
    {
        full_set[i] = i;
    }

    long num_nodes = 2 * n_points - 1;

    /* _nodes = (node_t *)malloc(num_nodes * sizeof(*_nodes));
    if (_nodes == NULL)
    {
        perror("Error allocating memory for _nodes");
        free(full_set);
        exit(1);
    } */

    nodes = (node_t **)malloc(num_nodes * sizeof(*nodes));
    if (nodes == NULL)
    {
        perror("Error allocating memory for nodes");
        free(full_set);
        // free(_nodes);
        exit(1);
    }

    // alloc space for root
    nodes[0] = (node_t *)malloc(sizeof(nodes[0]));
    node_t *root = nodes[0];
    root->subset = full_set;
    root->subset_len = n_points;
    build_tree()

    for (long i = 1; i < num_nodes; i++)
        if (i % p == pid)
            nodes[i] = (node_t *)malloc(sizeof(nodes[i]));
            build_tree(i, 0);



    /* int best_thr_level = ceil(log(p) / log(2));
    long num_complete_levels = (long)(log(num_nodes + 1) / log(2));
    int level_thr = best_thr_level < num_complete_levels - 1 ? best_thr_level : num_complete_levels - 1;

    long level_size = 1;

    for (int l = 0; l < level_thr; l++) // niveis
    {
        for (int i = 0; i < level_size; i++)
        {
            int id = get_id(i, l);
            build_tree(id, 0);
        }
        level_size *= 2;
    }

    for (int i = pid; i < level_size; i += p)
    {
        int id = get_id(i, level_thr);
        build_tree(id, 1);
    } */

    MPI_Barrier(MPI_COMM_WORLD);
    exec_time += MPI_Wtime();

    if (!pid)
    {
        fprintf(stderr, "%.1lf\n", exec_time);

        printf("%d %ld\n", n_dims, num_nodes);

        level_size = 1;
        for (int l = 0; l < level_thr; l++) // niveis
        {
            for (int i = 0; i < level_size; i++)
            {
                int id = get_id(i, l);
                print_node(id);
            }
            level_size *= 2;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < level_size; i++)
    {
        if (i % p == pid)
        {
            int id = get_id(i, level_thr);
            node_t *temp_root = nodes[id];
            print_tree(temp_root);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    dump_tree(root); // to the stdout!

    destroy_memory();

    MPI_Op_free(&mpi_maxloc_dbl_index);
    MPI_Type_free(&mpi_dbl_index);
    MPI_Finalize();
    return 0;
}