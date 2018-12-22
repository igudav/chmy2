#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>

double h = 0.01;

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
    // x_{-1} = x_n = 0
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

    return solution;
}

void double_runge_kutt(FILE *out, double a, double b,double init,
                       double(*f)(double x, double y))
{
    double y = init;
    for (double x = a; x <= b; x += h){
        fprintf(out, "%f ", y);
        y = y + h * (f(x, y) + f(x + h, y + f(x, y) * h)) / 2;
    }
}


void square_runge_kutt(FILE *out, double a, double b, double init,
                       double(*f)(double x, double y))
{
    double y = init;
    for (double x = a; x <= b; x += h){
        fprintf(out, "%f ", y);
        double k1 = f(x, y);
        double k2 = f(x + h/2, y + k1 * h /2);
        double k3 = f(x + h/2, y + k2 * h/ 2);
        double k4 = f(x + h, y + h * k3);
        y = y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    }
}

void double_runge_kutt_system(FILE *out1, FILE *out2, double a, double b,  double init1,  double init2,
                              double(*f1)(double x, double y1, double y2), double(*f2)(double x, double y1, double y2))
{
    double y1 = init1;
    double y2 = init2;
    for (double x = a; x <= b; x += h){
        fprintf(out1, "%f ", y1);
        fprintf(out2, "%f ", y2);
        double k1 = f1(x, y1, y2);
        double k21 = f2(x, y1, y2);
        double k2 = f1(x + h/2, y1 + k1 * h /2, y2 + k21 * h / 2);
        double k22 = f2(x + h/2, y1 + k1 * h /2, y2 + k21 * h / 2);

        y1 = y1 + h * k2;
        y2 = y2 + h * k22;
    }
}

void square_runge_kutt_system(FILE *out1, FILE *out2, double a, double b,  double init1,  double init2,
                              double(*f1)(double x, double y1, double y2), double(*f2)(double x, double y1, double y2))
{
    double y1 = init1;
    double y2 = init2;
    for (double x = a; x <= b; x += h){
        fprintf(out1, "%f ", y1);
        fprintf(out2, "%f ", y2);
        double k1 = f1(x, y1, y2);
        double k21 = f2(x, y1, y2);

        double k2 = f1(x + h/2, y1 + k1 * h /2, y2 + k21 * h / 2);
        double k22 = f2(x + h/2, y1 + k1 * h /2, y2 + k21 * h / 2);

        double k3 = f1(x + h/2, y1 + k2 * h /2, y2 + k22 * h / 2);
        double k23 = f2(x + h/2, y1 + k2 * h /2, y2 + k22 * h / 2);

        double k4 = f1(x + h, y1 + k3 * h , y2 + k23 * h);
        double k24 = f2(x + h, y1 + k3 * h , y2 + k23 * h);

        y1 = y1 + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        y2 = y2 + h * (k21 + 2 * k22 + 2 * k23 + k24) / 6;
    }
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