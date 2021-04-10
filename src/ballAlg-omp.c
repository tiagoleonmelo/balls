#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include "gen_points.c"

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

int n_dims;
long n_points;
int seed;
double **pts;
node_t **nodes;

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
    double max_dist = -1;
    double curr_dist;

    long *ret = (long *)malloc(sizeof(long) * 2);

    for (long i = 1; i < subset_len; i++)
    {
        curr_dist = distance(first, subset[i]);

        if (curr_dist > max_dist)
        {
            max_dist = curr_dist;
            ret[0] = i;
        }
    }

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

    // #pragma omp parallel for
    for (int dim = 0; dim < n_dims; dim++)
    {
        res[dim] = a[dim] - b[dim];
    }
}

double *orth_projection(long *subset, long subset_len, long a, long b, double *pt_a, double *b_minus_a_vec)
{

    double *proj = (double *)malloc(sizeof(double) * subset_len);

    double b_minus_a = b_minus_a_vec[0];

    double inner_prod_b_minus_a = inner_product(b_minus_a_vec, b_minus_a_vec);
    // #pragma omp parallel for
    for (long p = 0; p < subset_len; p++)
    {

        if (p == a || p == b)
        {
            proj[p] = pts[subset[p]][0];
            continue;
        }

        double *p_minus_a_vec = (double *)malloc(sizeof(double) * n_dims);
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
    memcpy(ordered_orth, orth, sizeof(double) * subset_len);

    long mid_index = subset_len / 2;

    double median = quickselect(ordered_orth, 0, subset_len - 1, mid_index + 1);

    new_pt = single_projection(median, pt_a, subset, b_minus_a_vec); // Pass along point

    if (subset_len % 2 == 0)
    {
        double *pt2;
        long mid_index2 = (subset_len / 2) - 1;
        double median2 = quickselect(ordered_orth, 0, subset_len - 1, mid_index2 + 1);
        pt2 = single_projection(median2, pt_a, subset, b_minus_a_vec);

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

void build_tree(long id, long l_id, long r_id)
{
    node_t *node = nodes[id];
    long subset_len = node->subset_len;
    long *subset = node->subset;

    if (subset_len == 0)
        return;

    if (subset_len > 1)
    {
        long *a_b;
        double *pt_a;
        double *pt_b;
        double *orth;
        double *median;
        double radius;

        long *subset_L = (long *)malloc(sizeof(long) * subset_len / 2);
        long *subset_R = (long *)malloc(sizeof(long) * (subset_len / 2 + (subset_len % 2)));

        if (subset_len == 2)
        {
            pt_a = pts[subset[0]];
            pt_b = pts[subset[1]];

            orth = (double *)malloc(sizeof(double) * subset_len);
            orth[0] = pt_a[0];
            orth[1] = pt_b[0];
            double total;
            total = 0;
            median = (double *)malloc(sizeof(double) * n_dims);
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
            pt_a = pts[subset[a_b[0]]];
            pt_b = pts[subset[a_b[1]]];

            double *b_minus_a_vec = (double *)malloc(sizeof(double) * n_dims);
            difference(pt_a, pt_b, b_minus_a_vec);

            // Orthogonal projection
            orth = orth_projection(subset, subset_len, a_b[0], a_b[1], pt_a, b_minus_a_vec);

            // Find median point
            median = find_median(orth, subset, subset_len, pt_a, b_minus_a_vec);

            // Find radius
            //#pragma omp task shared(radius)
            radius = find_radius(median, subset, subset_len);

            // Split between Left and Right branches
            //#pragma omp task
            split(orth, median, subset, subset_len, subset_L, subset_R);

            // Wait for tasks to finish
            //#pragma omp taskwait

            // Create node with id, center coords and radius
            fill_node(id, median, radius);

            free(b_minus_a_vec);
            free(a_b);
        }

        node_t *left = nodes[l_id];
        node->L = left;
        left->subset_len = subset_len / 2;
        left->subset = subset_L;

        node_t *right = nodes[r_id];
        node->R = right;
        right->subset_len = subset_len / 2 + (subset_len % 2);
        right->subset = subset_R;

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
    if (!(L == NULL && R == NULL))
        free(root->center_coord);

    // Recursive call
    dump_tree(L);
    dump_tree(R);
}

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        printf("Usage: %s <n_dims> <n_points> <seed>\n", argv[0]);
        exit(1);
    }

    n_dims = atoi(argv[1]);
    n_points = atol(argv[2]);
    seed = atoi(argv[3]);

    omp_set_nested(2);

    double exec_time;
    exec_time = -omp_get_wtime();

    pts = get_points(argc, argv, &n_dims, &n_points);
    long *full_set = (long *)malloc(sizeof(long) * n_points);

    for (long i = 0; i < n_points; i++)
    {
        full_set[i] = i;
    }

    long num_nodes = (1 << ((long)(log(2 * n_points - 1) / log(2)) + 1)) - 1; //2 * n_points - 1;
    long num_levels = (long)(log(num_nodes) / log(2)) + 1;

    node_t *_nodes = (node_t *)malloc(num_nodes * sizeof(*_nodes));
    nodes = (node_t **)malloc(num_nodes * sizeof(*nodes));
    for (long i = 0; i < num_nodes; i++)
        nodes[i] = _nodes + i;

    node_t *root = nodes[0];
    root->subset = full_set;
    root->subset_len = n_points;

    long base_id = 0;
    long level_size = 1;
    long child_base_id;
    // #pragma omp parrallel
    // #pragma omp single
    // {

    for (long l = 0; l < num_levels; l++) // niveis
    {
        child_base_id = base_id + level_size;
        for (long i = 0; i < level_size; i++)
        {
            build_tree(base_id + i, child_base_id + 2 * i, child_base_id + 2 * i + 1);
        }

        level_size *= 2;
        base_id = child_base_id;
    }

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1lf\n", exec_time);

    printf("%d %ld\n", n_dims, num_nodes);

    dump_tree(root); // to the stdout!

    free(nodes);
    free(_nodes);
    free(pts[0]);
    free(pts);
    return 0;
}
/* parallel single
for niveis= 0 ; niveis 
    pragma for
        for i= ponto
 */