#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include "gen_points.c"

#define LEFT 0
#define RIGHT 1

typedef struct node
{
    int id;
    double *center_coord;
    double radius;
    struct node *L;
    struct node *R;
} node_t;

int n_dims;
long n_points;
int seed;
double **pts;
node_t *tree;

node_t *create_node(int id, double *center_coord, double radius)
{
    // Allocate memory for new node
    node_t *node = (node_t *)malloc(sizeof(node_t));

    node->id = id;
    node->center_coord = center_coord;
    node->radius = radius;

    node->L = NULL;
    node->R = NULL;
    return node;
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

    for (int dim = 0; dim < n_dims; dim++)
    {
        res[dim] = a[dim] - b[dim];
    }
}

void sum(double *a, double *b, double *res)
{

    for (int dim = 0; dim < n_dims; dim++)
    {
        res[dim] = a[dim] + b[dim];
    }
}

void vector_scalar_product(double scalar, double *v, double *res)
{

    for (int dim = 0; dim < n_dims; dim++)
    {
        res[dim] = scalar * v[dim];
    }
}

double *orth_projection(long *subset, long subset_len, long a, long b, double *pt_a, double *b_minus_a_vec)
{

    double *proj = (double *)malloc(sizeof(double) * subset_len);

    double b_minus_a = b_minus_a_vec[0];

    double inner_prod_b_minus_a = inner_product(b_minus_a_vec, b_minus_a_vec);

    for (long p = 0; p < subset_len; p++)
    {

        if (p == a || p == b) // If p == a or p == b
        {
            proj[p] = pts[subset[p]][0];
            continue;
        }

        double p_minus_a = pts[subset[p]][0] - pt_a[0];

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

void quicksort(double *number, long first, long last)
{
    long i, j, pivot;
    double temp;

    if (first < last)
    {
        pivot = first;
        i = first;
        j = last;

        while (i < j)
        {
            while (number[i] <= number[pivot] && i < last)
                i++;
            while (number[j] > number[pivot])
                j--;
            if (i < j)
            {
                temp = number[i];
                number[i] = number[j];
                number[j] = temp;
            }
        }

        temp = number[pivot];
        number[pivot] = number[j];
        number[j] = temp;
        quicksort(number, first, j - 1);
        quicksort(number, j + 1, last);
    }
}

double *find_median(double *orth, long *subset, long subset_len, double *pt_a, double *b_minus_a_vec)
{
    double *new_pt;
    double *ordered_orth = (double *)malloc(sizeof(double) * subset_len);
    memcpy(ordered_orth, orth, sizeof(double) * subset_len);

    //qsort(ordered_orth, subset_len, sizeof(double), cmpfunc);
    quicksort(ordered_orth, 0, subset_len - 1);

    if (subset_len % 2 == 1)
    {
        // Take index, find corresponding element in orth, calculate orth_proj of that point, return it
        int mid_index = (int)subset_len / 2;

        for (long i = 0; i < subset_len; i++)
        {
            if (orth[i] == ordered_orth[mid_index])
            {
                new_pt = single_projection(orth[i], pt_a, subset, b_minus_a_vec); // Pass along point
                break;
            }
        }
    }
    else
    {
        double *pt1, *pt2;
        int idx1 = (int)(subset_len / 2) - 1;
        int idx2 = (int)subset_len / 2;

        for (long i = 0; i < subset_len; i++)
        {
            if (orth[i] == ordered_orth[idx1])
            {
                pt1 = single_projection(orth[i], pt_a, subset, b_minus_a_vec);
            }
            else if (orth[i] == ordered_orth[idx2])
            {
                pt2 = single_projection(orth[i], pt_a, subset, b_minus_a_vec);
            }
        }

        new_pt = (double *)malloc(sizeof(double) * n_dims);
        for (int i = 0; i < n_dims; i++)
        {
            new_pt[i] = (pt1[i] + pt2[i]) / 2.0;
        }

        free(pt1);
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
    // if (subset_len == 1 || subset_len == 0)
    // {
    //     return 0;
    // }

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

node_t *build_tree(long *subset, long subset_len, long id)
{
    node_t *root;

    if (subset_len > 1)
    {

        // Find A and B
        long *a_b = furthest_apart(subset, subset_len);
        double *pt_a = pts[subset[a_b[0]]];
        double *pt_b = pts[subset[a_b[1]]];

        double *b_minus_a_vec = (double *)malloc(sizeof(double) * n_dims);
        difference(pt_a, pt_b, b_minus_a_vec);

        // Orthogonal projection
        double *orth = orth_projection(subset, subset_len, a_b[0], a_b[1], pt_a, b_minus_a_vec);

        // Find median point
        double *median = find_median(orth, subset, subset_len, pt_a, b_minus_a_vec);

        // Find radius
        double radius = find_radius(median, subset, subset_len);

        // Create node with id, center coords and radius
        root = create_node(id, median, radius);

        // Split among L, R
        long *subset_L = (long *)malloc(sizeof(long) * subset_len / 2);
        long *subset_R = (long *)malloc(sizeof(long) * (subset_len / 2 + (subset_len % 2)));
        split(orth, median, subset, subset_len, subset_L, subset_R);

        root->L = build_tree(subset_L, subset_len / 2, id + 1);
        root->R = build_tree(subset_R, subset_len / 2 + (subset_len % 2), id + subset_len - (subset_len % 2));

        free(orth);
        free(b_minus_a_vec);
        free(a_b);
    }
    else
    {

        // Create leaf
        root = create_node(id, pts[subset[0]], 0);
    }

    free(subset);

    return root;
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
    free(root);

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

    double exec_time;
    exec_time = -omp_get_wtime();

    pts = get_points(argc, argv, &n_dims, &n_points);
    long *full_set = (long *)malloc(sizeof(long) * n_points);

    for (long i = 0; i < n_points; i++)
    {
        full_set[i] = i;
    }

    long num_nodes = 2 * n_points - 1;
    node_t *root = build_tree(full_set, n_points, 0);

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1lf\n", exec_time);

    printf("%d %ld\n", n_dims, num_nodes);
    dump_tree(root); // to the stdout!
    free(pts[0]);
    free(pts);
    return 0;
}