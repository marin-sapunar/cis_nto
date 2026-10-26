#ifndef CBLAS_H
#define CBLAS_H
typedef enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102} CBLAS_ORDER;
typedef enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112} CBLAS_TRANSPOSE;
void cblas_dgemm(const enum CBLAS_ORDER, const enum CBLAS_TRANSPOSE, const enum CBLAS_TRANSPOSE, const int, const int, const int, const double, const double*, const int, const double*, const int, const double, double*, const int);
void cblas_dgemv(const enum CBLAS_ORDER, const enum CBLAS_TRANSPOSE, const int, const int, const double, const double*, const int, const double*, const int, const double, double*, const int);
void cblas_dger(const enum CBLAS_ORDER, const int, const int, const double, const double*, const int, const double*, const int, double*, const int);
void cblas_daxpy(const int, const double, const double*, const int, double*, const int);
void cblas_dscal(const int, const double, double*, const int);
#endif
