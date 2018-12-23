#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <assert.h>

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

void rk2(FILE *out, double a, double b,double y0, double(*f)(double x, double y))
{
    double y = y0;
    for (double x = a; x <= b; x += h){
        fprintf(out, "%f ", y);
        y = y + h * (f(x, y) + f(x + h, y + f(x, y) * h)) / 2;
    }
    fprintf(out, "\n");
}


void rk4(FILE *out, double a, double b, double y0, double(*f)(double x, double y))
{
    double y = y0;
    for (double x = a; x <= b; x += h){
        fprintf(out, "%f ", y);
        double k1 = f(x, y);
        double k2 = f(x + h/2, y + k1 * h /2);
        double k3 = f(x + h/2, y + k2 * h/ 2);
        double k4 = f(x + h, y + h * k3);
        y = y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    }
    fprintf(out, "\n");
}

void rk2sys(FILE *out1, FILE *out2, double a, double b, double y0_1,  double y0_2,
        double(*f1)(double x, double y1, double y2), double(*f2)(double x, double y1, double y2))
{
    double y1 = y0_1;
    double y2 = y0_2;
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
    fprintf(out1, "\n");
    fprintf(out2, "\n");
}

void rk4sys(FILE *out1, FILE *out2, double a, double b, double y0_1, double y0_2,
        double(*f1)(double x, double y1, double y2), double(*f2)(double x, double y1, double y2))
{
    double y1 = y0_1;
    double y2 = y0_2;
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
    fprintf(out1, "\n");
    fprintf(out2, "\n");
}

// a_0*y(a) + b_0*y'(a) = c_0
// a_1*y(b) + b_1*y'(b) = c_1
// y" + y'*p(x) + y*q(x) = -f(x)
void bound_prob(FILE *out, double a_0, double b_0, double c_0, double a_1, double b_1, double c_1,
        double (*p)(double), double (*q)(double), double (*f)(double), double a, double b)
{
    size_t sz = (b - a) / h;
    assert(sz > 1);
    double **matrix = getMtrx(sz);
    double *vec = calloc(sz, sizeof(*f));

    matrix[0][0] = a_0 * h - b_0;
    matrix[0][1] = b_0;
    vec[0] = c_0 * h;
    matrix[sz - 1][sz - 2] = a_1 * h - b_1;
    matrix[sz - 1][sz - 1] = b_1;
    vec[sz - 1] = c_1 * h;

    for (int i = 1; i < sz - 1; ++i) {
        matrix[i][i - 1] = 2 - h * p(a + i * h);
        matrix[i][i] = -4 + 2 * h * h * q(a + i * h);
        matrix[i][i + 1] = 2 + h * p(a + i * h);
        vec[i] = -2 * h * h * f(a + i * h);
    }

    double *sol = tridiagMtrx(matrix, vec, sz);

    for (int i = 0; i < sz; ++i) {
        fprintf(out, "%lf ", sol[i]);
    }
    fprintf(out, "\n");

    delMtrx(matrix, sz);
    free(vec);
    free(sol);
}

double const5(double x)
{
    return -5.0;
}

double const4(double x)
{
    return 4.0;
}

double const8(double x)
{
    return -8.0;
}

int main(int argc, char *argv[])
{

    FILE *out = fopen("test2", "a");

    bound_prob(out, 1, 0, 1, 1, 0, 7, const5, const4, const8, 0.0, log(2));

    freopen("test2", "a", out);

    for (int i = 0; i < 100; ++i) {
        fprintf(out, "%lf ", 2 - 1.5 * exp(i * h) + 0.5 * exp(4 * i * h));
    }

    fprintf(out, "\n");

    return 0;
}