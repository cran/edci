#include "statist.h"

double mean (const double *const data, const int n)
{
    int i;
    double m;

    m=0;
    for (i=0;i<n;i++)
	m += data[i];
    
    return m/n;
}

void qsort (double *const array, int n, int links, int rechts, int q)
{
    double hilf;
    int l, r, mitte;
    
    l=links;
    r=rechts;
    
    //Rprintf("1) n=%i, l=%i, r=%i\n",n,l,r);
    if (l<r){
	mitte=(l+r)/2;
	//Rprintf("2) Mitte=%i\n",mitte);
	
	while (l!=r) {
	    while (array[l]<=array[mitte] && l<mitte)
		l++;
	    while (array[r]>=array[mitte] && r>mitte)
		r--;
	    //Rprintf("3) l=%i, r=%i, mitte=%i, %f, %f, %f\n",l,r,mitte,
	    //    array[l],array[mitte],array[r]);

	    if (l!=r){
		hilf=array[l];
		array[l]=array[r];
		array[r]=hilf;
		
		if (l==mitte) 
		    mitte=r;
		else if (r==mitte)
		    mitte=l;
	    }
	}
	//Rprintf("Mitte=%i, n=%i, q=%i, n/q-1=%f\n",mitte,n,q,n/q-1);
	if (q==1 || mitte>=n/q-1)
	    qsort(array,n,links,mitte-1,q);
	if(q==1 || mitte<=(n/q)+1)
	    qsort(array,n,mitte+1,rechts,q);
    }
}

double indexmean (const int *const index, const double *const data, 
		  const int n)
{
    int i;
    double m;

    m=0;
    for (i=0;i<n;i++)
	m += data[index[i]];
    
    return m/n;
}

void indexsort (int *const index, const double *const array, 
		int n, int links, int rechts, int q)
{
    int hilf;
    int l, r, mitte;
    
    l=links;
    r=rechts;
    
    if (l<r){
	mitte=(l+r)/2;
	
	while (l!=r) {
	    while (array[index[l]]<=array[index[mitte]] && l<mitte)
		l++;
	    while (array[index[r]]>=array[index[mitte]] && r>mitte)
		r--;

	    if (l!=r){
		hilf=index[l];
		index[l]=index[r];
		index[r]=hilf;
		
		if (l==mitte) 
		    mitte=r;
		else if (r==mitte)
		    mitte=l;
	    }
	}
	//Rprintf("Mitte=%i, n=%i, q=%i, n/q-1=%f\n",mitte,n,q,n/q-1);
	if (q==1 || mitte>=n/q-1)
	    indexsort(index,array,n,links,mitte-1,q);
	if(q==1 || mitte<=(n/q)+1)
	    indexsort(index,array,n,mitte+1,rechts,q);
    }
}
