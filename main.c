#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>

double **getMtrx(int n)
{
    double **matrix = calloc(n, sizeof(*matrix));
    if (!matrix) {
        return NULL;
    }
    for (int i = 0; i < n; ++i) {
        matrix[i] = calloc(n, sizeof(*matrix[i]));
        if (!matrix[i]) {
            return NULL;
        }
    }
    return matrix;
}

void delMtrx(double **matrix, int n)
{
    for (int i = 0; i < n; ++i) {
        free(matrix[i]);
    }
    free(matrix);
}

double **clone(double **matrix, int n)
{
    double **res = getMtrx(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            res[i][j] = matrix[i][j];
        }
    }
    return res;
}

void printVect(double *vector, int n)
{
    for (int i = 0; i < n; ++i) {
        printf("%lf ", vector[i]);
    }
    printf("\n");
}

void scanMtrx(double **matrix, int n)
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            scanf("%lf", &matrix[i][j]);
        }
    }
}

void scanVect(double *vector, int n)
{
    for (int i = 0; i < n; ++i) {
        scanf("%lf", &vector[i]);
    }
}

double *tridiagMtrx(double **matrix, double *f, int n)
{
    // x0 = xn = 0
    double *solution = calloc(n, sizeof(*solution));
    double *a = calloc(n, sizeof(*a));
    double *b = calloc(n, sizeof(*b));
    double *c = calloc(n, sizeof(*c));
    double *alpha = calloc(n, sizeof(*alpha));
    double *beta = calloc(n, sizeof(*beta));

    for (int i = 1; i < n; ++i) {
        a[i] = matrix[i][i - 1];
    }

    for (int i = 0; i < n - 1; ++i) {
        b[i] = matrix[i][i + 1];
    }

    for (int i = 0; i < n; ++i) {
        c[i] = matrix[i][i];
    }

    if (fabs(c[0]) < DBL_EPSILON) {
        perror("Плохая матрица!\n");
        exit(1);
    }

    for (int i = 0; i < n; ++i) {
        if (fabs(c[i] + a[i] * (i == 0 ? 0.0 : alpha[i - 1])) < DBL_EPSILON) {
            perror("Плохая матрица\n!");
            exit(1);
        }
        alpha[i] = -b[i] / (c[i] + a[i] * (i == 0 ? 0.0 : alpha[i - 1]));
    }

    for (int i = 0; i < n; ++i) {
        beta[i] = (f[i] - a[i] * (i == 0 ? 0.0 : beta[i - 1])) /
                (c[i] + a[i] * (i == 0 ? 0.0 : alpha[i - 1]));
    }

    for (int i = n - 1; i >= 0; --i) {
        solution[i] = beta[i] + alpha[i] * (i == n - 1 ? 0.0 : solution[i + 1]);
    }

    free(a);
    free(b);
    free(c);
    free(alpha);
    free(beta);
    delMtrx(matrix, n);

    return solution;
}

int main(int argc, char *argv[])
{
    int n;
    scanf("%d", &n);
    double **matrix = getMtrx(n);
    double *f = calloc(n, sizeof(*f));
    scanMtrx(matrix, n);
    scanVect(f, n);

    double *sol = tridiagMtrx(matrix, f, n);
    printVect(sol, n);

    free(sol);

    return 0;
}