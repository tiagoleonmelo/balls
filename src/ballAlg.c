#include <stdio.h>
#include <omp.h>
#include <math.h>
#include "gen_points.c"

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
        curr_dist = distance(ret[0], subset[i]);

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


void vector_scalar_product(double scalar, double* v, double* res)
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

    double *b_minus_a = (double *)malloc(sizeof(double) * n_dims);
    difference(pts[b], pts[a], b_minus_a);

    double **proj = (double **) malloc(sizeof(double *) * subset_len);

    for (long p = 0; p < subset_len; p++)
    {
        double *p_minus_a = (double *)malloc(sizeof(double) * n_dims);
        difference(pts[subset[p]], pts[a], p_minus_a);

        double scalar = inner_product(p_minus_a, b_minus_a) / inner_product(b_minus_a, b_minus_a);
        
        proj[p] = (double *)malloc(sizeof(double) * n_dims);
        vector_scalar_product(scalar, b_minus_a, proj[p]);
        sum(proj[p], pts[a], proj[p]);

        free(p_minus_a);
    }

    free(b_minus_a);

    return proj;
}

node_t *build_tree(long *subset, long subset_len)
{

    // Find A and B
    long *a_b = furthest_apart(subset, subset_len);

    // Orthogonal projection
    double **orth = orth_projection(subset, subset_len, a_b);

    for (int i = 0; i < 5; i++)
    {
        printf("%f %f\n", orth[i][0], orth[i][1]);
    }
    

    double *center_coords = (double *)malloc(sizeof(double) * n_dims);

    center_coords[0] = 5.1;
    center_coords[1] = 219.0;

    node_t *root = create_node(1290, center_coords, 3.2);

    // if subset > 0
    // root->L = build_tree(subset_L);
    // root->R = build_tree(subset_R);

    free(a_b);

    return root;
}

void dump_tree(node_t *root)
{

    node_t *L = root->L;
    node_t *R = root->R;

    // printf("%d %d %d %f ", root->id, L->id, R->id, root->radius);

    for (size_t i = 0; i < n_dims; i++)
    {
        printf("%f ", root->center_coord[i]);
    }

    if (root->L == NULL && root->R == NULL)
    {
        return;
    }

    // Recursive call?
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

    node_t *root = build_tree(full_set, n_points);

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1lf\n", exec_time);

    free(full_set);
    dump_tree(root); // to the stdout!
    return 0;
}
