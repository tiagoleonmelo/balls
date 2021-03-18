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

    for (int i = 0; i < n_dims; i++)
    {
        total += ((pts[a][i] - pts[b][i]) * (pts[a][i] - pts[b][i]));
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

    for (long i = 0; i < subset_len; i++)
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

double **orth_projection(long *subset, long subset_len, long *a_b)
{

    long a = subset[a_b[0]];
    long b = subset[a_b[1]];

    double **proj = (double **)malloc(sizeof(double *) * subset_len);

    double *b_minus_a = (double *)malloc(sizeof(double) * n_dims);
    difference(pts[b], pts[a], b_minus_a);

    for (long p = 0; p < subset_len; p++)
    {
        proj[p] = (double *)malloc(sizeof(double) * n_dims);

        if (subset[p] == a || subset[p] == b)
        {
            memcpy(proj[p], pts[subset[p]], sizeof(double) * n_dims);
            continue;
        }

        double *p_minus_a = (double *)malloc(sizeof(double) * n_dims);
        difference(pts[subset[p]], pts[a], p_minus_a);

        double scalar = inner_product(p_minus_a, b_minus_a) / inner_product(b_minus_a, b_minus_a);

        vector_scalar_product(scalar, b_minus_a, proj[p]);
        sum(proj[p], pts[a], proj[p]);

        free(p_minus_a);
    }

    free(b_minus_a);

    return proj;
}

int cmpfunc(const void *pt_1, const void *pt_2)
{
    double const *p1 = *(double const **)pt_1;
    double const *p2 = *(double const **)pt_2;

    if (p1[0] < p2[0])
        return -1;
    if (p1[0] > p2[0])
        return 1;
    return 0;
}

double *find_median(double **orth, long subset_len)
{
    double **ordered_orth = (double **)malloc(sizeof(double *) * subset_len);
    memcpy(ordered_orth, orth, sizeof(double *) * subset_len);

    qsort(ordered_orth, subset_len, sizeof(orth[0]), cmpfunc);

    double *new_pt = (double *)malloc(sizeof(double) * n_dims);
    if (subset_len % 2 == 1)
    {
        int mid_index = (int)subset_len / 2;
        memcpy(new_pt, ordered_orth[mid_index], sizeof(double) * n_dims);
    }
    else
    {
        int idx1 = (int)(subset_len / 2) - 1;
        int idx2 = (int)subset_len / 2;

        for (int i = 0; i < n_dims; i++)
        {
            new_pt[i] = (ordered_orth[idx1][i] + ordered_orth[idx2][i]) / 2.0;
        }
    }

    free(ordered_orth);
    return new_pt;
}

void split(double **orth, double *median, long* subset, long subset_len, long *left, long *right)
{
    long ind_left = 0;
    long ind_right = 0;

    for (int i = 0; i < subset_len; i++)
    {
        if (median[0] > orth[i][0])
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
    if (subset_len == 1 || subset_len == 0)
    {
        return 0;
    }

    double current_dist;
    double max_dist = -1;
    double total;

    for (long i = 0; i < subset_len; i++)
    {
        current_dist = 0;
        total = 0;
        for (int d = 0; d < n_dims; d++)
        {
            total += ((center[d] - pts[subset[i]][d]) * (center[d] - pts[subset[i]][d]));
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

        // Orthogonal projection
        double **orth = orth_projection(subset, subset_len, a_b);

        // Find median point
        double *median = find_median(orth, subset_len);

        // Find radius
        double radius = find_radius(median, subset, subset_len);

        // printf("Radius: %f, subset_len %ld, A: %ld, B: %ld\n", radius, subset_len, a_b[0], a_b[1]);

        // Create node with id, center coords and radius
        root = create_node(id, median, radius);

        // Split among L, R
        long *subset_L = (long *)malloc(sizeof(long) * subset_len / 2);
        long *subset_R = (long *)malloc(sizeof(long) * (subset_len / 2 + (subset_len % 2)));
        split(orth, median, subset, subset_len, subset_L, subset_R);



        root->L = build_tree(subset_L, subset_len / 2, id + 1);
        root->R = build_tree(subset_R, subset_len / 2 + (subset_len % 2), id + subset_len - (subset_len % 2));

        for (long i = 0; i < subset_len; i++)
        {
            free(orth[i]);
        }
        free(orth);
        free(a_b);
    }
    else
    {

        // Create leaf
        root = create_node(id, pts[subset[0]], 0);

        // printf("Radius: %d\n", 0);
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

    // Recursive call?
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
    return 0;
}
