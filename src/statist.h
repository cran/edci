#ifndef _statisth_
#define _statisth_

double mean (const double *const data, const int n);
void qsort (double *const array, int n, int links, int rechts, int q);

double indexmean (const int *const index, const double *const data, 
		  const int n);
void indexsort (int *const index, const double *const array, 
		int n, int links, int rechts, int q);

#endif _statisth_
