#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>

double h;

void printMtrx(double **matrix, int n)
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%7lf ", matrix[i][j]);
        }
        printf("\n");
    }
}

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

void rk2(FILE *out, double a, double b,double y0, double (*f)(double, double))
{
    double y = y0;
    fprintf(out, "X;Y\n");
    for (double x = a; x <= b; x += h){
        fprintf(out, "%lf;%lf\n", x, y);
        y = y + h / 2 * (f(x, y) + f(x + h, y + f(x, y) * h));
    }
}


void rk4(FILE *out, double a, double b, double y0, double (*f)(double, double))
{
    double y = y0;
    fprintf(out, "X;Y\n");
    for (double x = a; x <= b; x += h){
        fprintf(out, "%lf;%lf\n", x, y);
        double k1 = f(x, y);
        double k2 = f(x + h/2, y + k1 * h /2);
        double k3 = f(x + h/2, y + k2 * h/ 2);
        double k4 = f(x + h, y + h * k3);
        y = y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    }
}

void rk2sys(FILE *out, double a, double b, double y0_1,  double y0_2,
        double (*f1)(double, double, double), double (*f2)(double, double, double))
{
    double y1 = y0_1;
    double y2 = y0_2;
    fprintf(out, "X;Y1;Y2\n");
    for (double x = a; x <= b; x += h){
        fprintf(out, "%lf;%lf;%lf\n", x, y1, y2);
        double k1 = f1(x, y1, y2);
        double k2 = f2(x, y1, y2);
        y1 = y1 + h / 2 * (f1(x, y1, y2) + f1(x + h, y1 + k1 * h, y2 + k2 * h));
        y2 = y2 + h / 2 * (f2(x, y1, y2) + f2(x + h, y1 + k1 * h, y2 + k2 * h));
    }
}

void rk4sys(FILE *out, double a, double b, double y0_1, double y0_2,
        double (*f1)(double, double, double), double (*f2)(double, double, double))
{
    double y1 = y0_1;
    double y2 = y0_2;
    fprintf(out, "X;Y1;Y2\n");
    for (double x = a; x <= b; x += h){
        fprintf(out, "%lf;%lf;%lf\n", x, y1, y2);
        double k11 = f1(x, y1, y2);
        double k21 = f2(x, y1, y2);

        double k12 = f1(x + h/2, y1 + k11 * h /2, y2 + k21 * h / 2);
        double k22 = f2(x + h/2, y1 + k11 * h /2, y2 + k21 * h / 2);

        double k13 = f1(x + h/2, y1 + k12 * h /2, y2 + k22 * h / 2);
        double k23 = f2(x + h/2, y1 + k12 * h /2, y2 + k22 * h / 2);

        double k14 = f1(x + h, y1 + k13 * h , y2 + k23 * h);
        double k24 = f2(x + h, y1 + k13 * h , y2 + k23 * h);

        y1 = y1 + h * (k11 + 2 * k12 + 2 * k13 + k14) / 6;
        y2 = y2 + h * (k21 + 2 * k22 + 2 * k23 + k24) / 6;
    }
}

// a_0*y(a) + b_0*y'(a) = c_0
// a_1*y(b) + b_1*y'(b) = c_1
// y" + y'*p(x) + y*q(x) = -f(x)
void bound_prob(FILE *out, double a_0, double b_0, double c_0, double a_1, double b_1, double c_1,
        double (*p)(double), double (*q)(double), double (*f)(double), double a, double b)
{
    int sz = round((b - a) / h) + 1;

    double **matrix = getMtrx(sz);
    double *vec = calloc(sz, sizeof(*vec));

    matrix[0][0] = a_0 * h - b_0;
    matrix[0][1] = b_0;
    vec[0] = c_0 * h;
    matrix[sz - 1][sz - 2] = -b_1;
    matrix[sz - 1][sz - 1] = a_1 * h + b_1;
    vec[sz - 1] = c_1 * h;

    for (int i = 1; i < sz - 1; ++i) {
        matrix[i][i - 1] = 2 - h * p(a + i * h);
        matrix[i][i] = 2 * h * h * q(a + i * h) - 4;
        matrix[i][i + 1] = 2 + h * p(a + i * h);
        vec[i] = -2 * h * h * f(a + i * h);
    }

    double *sol = tridiagMtrx(matrix, vec, sz);

    fprintf(out, "X;Y\n");
    for (int i = 0; i < sz; ++i) {
        fprintf(out, "%lf;%lf\n", a + i * h, sol[i]);
    }

    delMtrx(matrix, sz);
    free(vec);
    free(sol);
}

double smp1_f(double x, double y) { return sin(x) - y; }

double smp2_f(double x, double y) { return y - y * x; }

double smp3_f(double x, double y) { return (x - x * x) * y; }

double sys1_f1(double x, double u, double v) { return cos(u + 1.1 * v) + 2.1; }

double sys1_f2(double x, double u, double v) { return 1.1 / (x + 2.1 * u * u) + x + 1; }

// y1' = y1 - y2 + 2sinx
// y2' = 2y1 - y2
// #840
double sys2_f1(double x, double u, double v) { return u - v + 2 * sin(x); }

double sys2_f2(double x, double u, double v) { return 2 * u - v; }

// y1' = 2y1 - y2
// y2' = y1 + 2e^x
// #841
double sys3_f1(double x, double u, double v) { return 2 * u - v; }

double sys3_f2(double x, double u, double v) { return u + 2 * exp(x); }

double p_1(double x) { return 2; }

double q_1(double x) { return -1 / x; }

double f_1(double x) { return -3; }

// y" - 2y' + y = 6xe^x
double p_2(double x) { return -2; }

double q_2(double x) { return 1; }

double f_2(double x) { return -6 * x * exp(x); }

//y" + 1/x * y' + 1/x^2 * y = 1/x^3
double p_3(double x) { return 1 / x; }

double q_3(double x) { return 1 / (x * x); }

double f_3(double x) { return -1 / (x * x * x); }


int main(int argc, char *argv[])
{
    h = 0.05;

    // тесты для задачи Коши
    FILE *smpd1 = fopen("../../../Документы/ЧМЫ_02/res/simpled1.csv", "w");
    FILE *smpd2 = fopen("../../../Документы/ЧМЫ_02/res/simpled2.csv", "w");
    FILE *smpd3 = fopen("../../../Документы/ЧМЫ_02/res/simpled3.csv", "w");
    FILE *smpq1 = fopen("../../../Документы/ЧМЫ_02/res/simpleq1.csv", "w");
    FILE *smpq2 = fopen("../../../Документы/ЧМЫ_02/res/simpleq2.csv", "w");
    FILE *smpq3 = fopen("../../../Документы/ЧМЫ_02/res/simpleq3.csv", "w");

    // тесты для системы
    FILE *sysd1 = fopen("../../../Документы/ЧМЫ_02/res/systemd1.csv", "w");
    FILE *sysd2 = fopen("../../../Документы/ЧМЫ_02/res/systemd2.csv", "w");
    FILE *sysd3 = fopen("../../../Документы/ЧМЫ_02/res/systemd3.csv", "w");
    FILE *sysq1 = fopen("../../../Документы/ЧМЫ_02/res/systemq1.csv", "w");
    FILE *sysq2 = fopen("../../../Документы/ЧМЫ_02/res/systemq2.csv", "w");
    FILE *sysq3 = fopen("../../../Документы/ЧМЫ_02/res/systemq3.csv", "w");

    // тесты для краевой задачи
    FILE *bnd1 = fopen("../../../Документы/ЧМЫ_02/res/bound1.csv", "w");
    FILE *bnd2 = fopen("../../../Документы/ЧМЫ_02/res/bound2.csv", "w");
    FILE *bnd3 = fopen("../../../Документы/ЧМЫ_02/res/bound3.csv", "w");

    rk2(smpd1, 0, 1, 10.0, smp1_f);
    rk2(smpd2, 0, 1, 5.0, smp2_f);
    rk2(smpd3, 0, 1, 1.0, smp3_f);
    rk4(smpq1, 0, 1, 10.0, smp1_f);
    rk4(smpq2, 0, 1, 5.0, smp2_f);
    rk4(smpq3, 0, 1, 1.0, smp3_f);

    rk2sys(sysd1, 0, 1, 0, 0, sys1_f1, sys1_f2);
    rk2sys(sysd2, 0, 1, 0, 0, sys2_f1, sys2_f2);
    rk2sys(sysd3, 0, 1, 0.5, 1, sys3_f1, sys3_f2);
    rk4sys(sysq1, 0, 1, 0, 0, sys1_f1, sys1_f2);
    rk4sys(sysq2, 0, 1, 0, 0, sys2_f1, sys2_f2);
    rk4sys(sysq3, 0, 1, 0.5, 1, sys3_f1, sys3_f2);

    bound_prob(bnd1, 1, 0, 2, 0.5, -1, 1, p_1, q_1, f_1, 0.2, 0.5);
    bound_prob(bnd2, 1, 0, 1, 2, -1, 0.5, p_2, q_2, f_2, 1.2, 1.5);
    bound_prob(bnd3, 1, 0, 1.2, 1, 2, 1.4, p_3, q_3, f_3, 0.4, 0.7);

    return 0;
}