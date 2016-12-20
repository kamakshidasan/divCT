#ifndef PSORT_H_
#define PSORT_H_
#include <cstdio>
#include <omp.h>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include "const.h"

////////////////////////////////////////////////////////////
// Miscellaneous functions
////////////////////////////////////////////////////////////

/** Swap to value */
template <class NumType>
inline void Swap(NumType& value, NumType& other){
    NumType temp = value;
    value = other;
    other = temp;
}


namespace
{
    struct NewCompare
    {
        bool operator()(const int i, const int j ) const
        {
        	if(function_values[i] < function_values[j]) {
        			return true;
        		}
        		if(function_values[i] == function_values[j]) {
        			if(vertex_pos[i] < vertex_pos[j]) {
        				return true;
        			}
        			if(vertex_pos[i] == vertex_pos[j]) {
        				if(offset[i] < offset[j]) {
        					return true;
        				}
        				if(offset[i] == offset[j]) {
        					if(i < j) {
        						return true;
        					}
        				}
        			}
        		}
        		return false;
        }

    };

   }
////////////////////////////////////////////////////////////
// Quick sort
////////////////////////////////////////////////////////////

/* use in the sequential qs */
template <class SortType>
long QsPartition(SortType outputArray[], long left, long right){
    const long part = right;
    Swap(outputArray[part],outputArray[left + (right - left ) / 2]);
    const SortType partValue = outputArray[part];
    --right;

    while(true){
//        while(outputArray[left] < partValue){
    	while(lessVertex(outputArray[left],partValue)){
            ++left;
        }
//        while(right >= left && partValue <= outputArray[right]){
        while(right >= left && !lessVertex(outputArray[right],partValue)){
            --right;
        }
        if(right < left) break;

        Swap(outputArray[left],outputArray[right]);
        ++left;
        --right;
    }

    Swap(outputArray[part],outputArray[left]);

    return left;
}

/* a sequential qs */
template <class SortType>
void QsSequential(SortType array[], const long left, const long right){
    if(left < right){
        const long part = QsPartition(array, left, right);
        QsSequential(array,part + 1,right);
        QsSequential(array,left,part - 1);
    }
}

/** A task dispatcher */
template <class SortType>
void QuickSortOmpTask(SortType array[], const long left, const long right, const int deep){
    if(left < right){
        if( deep ){
            const long part = QsPartition(array, left, right);
            #pragma omp task
            QuickSortOmpTask(array,part + 1,right, deep - 1);
            #pragma omp task
            QuickSortOmpTask(array,left,part - 1, deep - 1);
        }
        else {
            const long part = QsPartition(array, left, right);
            QsSequential(array,part + 1,right);
            QsSequential(array,left,part - 1);
        }
    }
}

/** The openmp quick sort */
template <class SortType>
void QuickSortOmp(SortType array[], const long size, int deep){
	//qsort(array,size,sizeof(int),compVertex);
	//return;
	int l =0;
	while (deep >>= 1) { ++l; }
	omp_set_nested(1);
    #pragma omp parallel
    {
        #pragma omp single nowait
        {
            QuickSortOmpTask(array, 0, size - 1 , deep);
        }
    }
}

////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////

bool isSorted(int array[], const int size){
    for(int idx = 1; idx < size ; ++idx){
        if(!lessVertex(array[idx-1],array[idx])){
            return false;
        }
    }
    return true;
}

void print(long long array[], const int size){
    for(int idx = 0 ;idx < size; ++idx){
        printf("%lld\t",array[idx]);
    }
    printf("\n");
}

void seqmerge(int a[],int size, int temp[]){ //standard code to merge two sorted sequence of elements
    int i,m,k,l;
    int mid = size/2 -1 ;
    int low = 0;
    int high = size-1;
    l = low;
    i = low;
    m = mid+1;

    while((l<=mid)&&(m<=high)){

         // lessVertex compares the function values and then sorts the indices
         if(lessVertex(a[l],a[m])){
             temp[i]=a[l];
             l++;
         }
         else{
             temp[i]=a[m];
             m++;
         }
         i++;
    }

    if(l>mid){
         for(k=m;k<=high;k++){
             temp[i]=a[k];
             i++;
         }
    }
    else{
         for(k=l;k<=mid;k++){
             temp[i]=a[k];
             i++;
         }
    }

    for(k=low;k<=high;k++){
         a[k]=temp[k];
    }
}

void mergesort(int a[], int size, int temp[], int threads) { //parallel merge sort
//a is the array base address, size is the length of the array this instance of merge sort is acting upon
//temp is the temporary array required, threads is the available no. of total threads that can be created at this level
    if (threads <= 1) { //run out of max no. of threads
        if (size>1){ //still elements to sort
            qsort(a,size,sizeof(int),compVertex);
            //std::sort(a, a+size,NewCompare());
            //sort using the library sort function deterministic quick sort for larger values and insertion sort for smaller values
            //Note: merge sort is not used as base case as it would incur in high no. of function call overheads
        }
    }
    else if (threads > 1) {
        #pragma omp parallel sections num_threads(2) // Call the two parallel sections below with two children threads
        {
            #pragma omp section // parallel section 1 acting on the first half of the parent array
            {
                mergesort(a, size/2, temp, threads/2); // hence the size is halved as well as the no. of threads
            }

            #pragma omp section // parallel section 2 acting on the last half of the parent array
            {
                mergesort(a + size/2, size - size/2,temp + size/2, threads - threads/2);
                //address of a and temp are updated as well as the size and number of threads
            }
        }
    }

    seqmerge(a,size,temp); //finally merge the children array by a single thread
}

#endif
