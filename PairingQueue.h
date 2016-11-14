/*
 * PairingQueue.h
 *
 *  Created on: 11-May-2012
 *      Author: harish
 */

#ifndef PAIRINGQUEUE_H_
#define PAIRINGQUEUE_H_

#include "const.h"

struct HeapNode {
	int v;
	int next;
	int child;

	int prev;
};


/*
 * Irrespective of min or max heaps
 */

int getTop(HeapNode * heap, int heap1) {
	return heap[heap1].v;
}

void createHeapNode(HeapNode * heap, int heap1, int v) {
	heap[heap1].v = v;
	heap[heap1].next = -1;
	heap[heap1].child = -1;
}

/*
 * Merges two Min heaps. heap1 and heap2 are indices to roots of the two heaps to merge
 * returns index of new heap
 */
int mergeMinHeap(HeapNode * heap, int heap1, int heap2) {
	if(heap1 == -1) {
		return heap2;
	}
	if(heap2 == -1) {
		return heap1;
	}

	int v1 = heap[heap1].v;
	int v2 = heap[heap2].v;
	if(!lessVertex(v1, v2)) {
		int t = heap1;
		heap1 = heap2;
		heap2 = t;
	}

	heap[heap2].next = heap[heap1].child;
	heap[heap1].child = heap2;
	heap[heap1].next = -1;
	return heap1;
}

/*
 * deletes the minimum element of heap1 (assumes a min heap)
 * returns index of new heap
 */
int deleteMin(HeapNode * heap, int heap1) {
	int pos = heap[heap1].child;
	heap[heap1].child = -1;
	if(pos == -1) {
		return -1;
	}

	int h1 = pos;
	int h2 = heap[h1].next;
	int next = -1;
	if(h2 != -1) {
		next = heap[h2].next;
	}
	int p = -1;
	while(h1 != -1) {
		int h = mergeMinHeap(heap,h1,h2);
		p = mergeMinHeap(heap,p,h);
		h1 = next;
		h2 = -1;
		if(h1 != -1) {
			h2 = heap[h1].next;
		}
		next = -1;
		if(h2 != -1) {
			next = heap[h2].next;
		}
	}
	return p;
}


/*
 * Merges two Max heaps. heap1 and heap2 are indices to roots of the two heaps to merge
 * returns index of new heap
 */
int mergeMaxHeap(HeapNode * heap, int heap1, int heap2) {
	if(heap1 == -1) {
		return heap2;
	}
	if(heap2 == -1) {
		return heap1;
	}

	int v1 = heap[heap1].v;
	int v2 = heap[heap2].v;
	if(lessVertex(v1, v2)) {
		int t = heap1;
		heap1 = heap2;
		heap2 = t;
	}

	heap[heap2].next = heap[heap1].child;
	heap[heap1].child = heap2;
	heap[heap1].next = -1;
	return heap1;
}

/*
 * deletes the minimum element of heap1 (assumes a min heap)
 * returns index of new heap
 */
int deleteMax(HeapNode * heap, int heap1) {
	int pos = heap[heap1].child;
	heap[heap1].child = -1;
	if(pos == -1) {
		return -1;
	}

	int h1 = pos;
	int h2 = heap[h1].next;
	int next = -1;
	if(h2 != -1) {
		next = heap[h2].next;
	}
	int p = -1;
	while(h1 != -1) {
		int h = mergeMaxHeap(heap,h1,h2);
		p = mergeMaxHeap(heap,p,h);
		h1 = next;
		h2 = -1;
		if(h1 != -1) {
			h2 = heap[h1].next;
		}
		next = -1;
		if(h2 != -1) {
			next = heap[h2].next;
		}
	}
	return p;
}


#endif /* PAIRINGQUEUE_H_ */
