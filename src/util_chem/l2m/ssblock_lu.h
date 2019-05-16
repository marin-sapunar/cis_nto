#ifndef SSBLOCK_LU_H
#define SSBLOCK_LU_H

void ssblock_lu(double *csc, const int no, const int nv1, const int nv2, const int ns1, const int ns2, double *wf1, double *wf2, double *l1minor, double *ss, const int print_level);

#endif
