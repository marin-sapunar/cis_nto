/*
 * FUNCTION: c_wrapper
 *
 * DESCRIPTION:
 * This is a file with the wrapper functions for the C functions
 */

#include <stdio.h>

#include "ssblock_lu.h"

void ssblock_(double *csc, const int *no, const int *nv1, const int *nv2, const int *ns1, const int *ns2, double *wf1, double *wf2, double *l1minor, double *ss, const int *print_level)
{
	int no_int = *no;
	int nv1_int = *nv1;
	int nv2_int = *nv2;
	int ns1_int = *ns1;
	int ns2_int = *ns2;
	int print_level_int = *print_level;

	ssblock_lu(csc, no_int, nv1_int, nv2_int, ns1_int, ns2_int, wf1, wf2, l1minor, ss, print_level_int);

}
