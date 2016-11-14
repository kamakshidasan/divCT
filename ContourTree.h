#ifndef CONTOUTTREE_H
#define CONTOUTTREE_H

#include <iostream>
#include "StoreReebGraph.h"
#include <stdio.h>
#include <string.h>

using namespace std;

typedef	struct JSTree {
	int * next;
	int * prev;
	int * nsize;
	int * psize;
	bool split;
	bool * present;
	
	JSTree(int no, bool split, int maxStar) {
		this->split = split;
		
		if(split) {
			next = new int[no * maxStar];
			prev = new int[no];
		} else {
			next = new int[no];
			prev = new int[no * maxStar];
		}
		nsize = new int[no];
		psize = new int[no];
		present = new bool[no];
		memset(present, false,no);
		memset(nsize,0,no*sizeof(int));
		memset(psize,0,no*sizeof(int));
	}
} JSTree;

class ContourTree {
	
public:
	JSTree * joinTree;
	JSTree * splitTree;
	
	int nv; 
	int maxStar;
	int * joinNodes;
	int * splitNodes;
	int jct;
	int sct;
	int * q;

	ContourTree(int noVertices, int maxStar) {
		jct = 0;
		sct = 0;
		nv = noVertices;
		this -> maxStar = maxStar;
		
		joinTree = new JSTree(nv, false, this -> maxStar);
		splitTree = new JSTree(nv, true, this -> maxStar);
		
		joinNodes = new int[nv];
		splitNodes = new int[nv];
		
		q = new int[nv];

		memset(joinNodes, 0, (nv) * sizeof(int));
		memset(splitNodes, 0, (nv) * sizeof(int));
		memset(q, 0, (nv) * sizeof(int));
	}

	void addJoinArc(int from, int to) {
		if(!joinTree->present[from]) {
			joinTree->present[from] = true;
			joinNodes[jct ++] = from;
		}
		if(!joinTree->present[to]) {
			joinTree->present[to] = true;
			joinNodes[jct ++] = to;
		}
		
		joinTree->next[from] = to;
		joinTree->nsize[from] ++;
		joinTree->prev[to * maxStar + joinTree->psize[to]] = from;
		joinTree->psize[to] ++;
	}
	
	void addSplitArc(int from, int to) {
		if(!splitTree->present[from]) {
			splitTree->present[from] = true;
			splitNodes[sct ++] = from;
		}
		if(!splitTree->present[to]) {
			splitTree->present[to] = true;
			splitNodes[sct ++] = to;
		}
		
		splitTree->next[from * maxStar + splitTree->nsize[from]] = to;
		splitTree->prev[to] = from;
		splitTree->nsize[from] ++;
		splitTree->psize[to] ++;
	}
	
	void remove(int xi, JSTree * tree) {
		if(tree->psize[xi] == 1 && tree->nsize[xi] == 1) {
			int p = tree->prev[xi];
			int n = tree->next[xi];
			int pmul = 1;
			int nmul = 1;
			if(tree->split) {
				n = tree->next[xi * maxStar];
				nmul = maxStar;
			} else {
				p = tree->prev[xi * maxStar];
				pmul = maxStar;
			}
			tree->psize[xi] = 0;
			tree->nsize[xi]  = 0;
			
			removeAndAdd(tree->next, p*nmul, tree->nsize[p], xi, n);
			removeAndAdd(tree->prev, n*pmul, tree->psize[n], xi, p);
		} else if(tree->psize[xi] == 0 && tree->nsize[xi] == 1) {
			int n = tree->next[xi];
			int pmul = 1;
			if(tree->split) {
				n = tree->next[xi * maxStar];
			} else {
				pmul = maxStar;
			}
			tree->nsize[xi] = 0;
			remove(tree->prev, n * pmul, tree->psize[n], xi);
			tree->psize[n] --;
		} else if(tree->psize[xi] == 1 && tree->nsize[xi] == 0) {
			int p = tree->prev[xi];
			int nmul = 1;
			if(!tree->split) {
				p = tree->prev[xi * maxStar];
			} else {
				nmul = maxStar;
			}
			tree->psize[xi] = 0;
			remove(tree->next, p * nmul, tree->nsize[p], xi);
			tree->nsize[p] --;
		} else {
			println("Can this too happen??????");
			exit(0);
		}
	}

	void remove(int* arr, int start, int arrSize, int xi) {
		for(int i = start;i < start + arrSize;i ++) {
			if(arr[i] == xi) {
				if(i != start + arrSize - 1) {
					arr[i] = arr[start + arrSize - 1];
				}
				return;
			}
		}
		println("Shouldn't happen 3");
		exit(0);
	}

	void removeAndAdd(int * arr, int start, int arrSize, int rem, int add) {
		for(int i = start;i < start + arrSize;i ++) {
			if(arr[i] == rem) {
				arr[i] = add;
				return;
			}
		}
		println("Shouldn't happen 2");
		exit(0);
	}

	void mergeTrees(StoreReebGraph * rg) {
		int front = 0;
		int back = 0;
		int added = 0;
		for(int x = 0;x < jct; x++) {
			int v = joinNodes[x];
			if(splitTree->nsize[v] + joinTree->psize[v] == 1) {
				q[back ++] = v;
			}
		}
		
		while(back > front + 1) {
			int xi = q[front ++];
			if(splitTree->nsize[xi] == 0 && splitTree->psize[xi] == 0) {
				if(!(joinTree->nsize[xi] == 0 && joinTree->psize[xi] == 0)) {
					println("Shoudn't happen 1!!!");
					cout << "added : " << added << endl;
					exit(0);
				}
				continue;
			}
			if(splitTree->nsize[xi] == 0) {
				if(splitTree->psize[xi] > 1) {
					println("Can this happen too???");
					exit(0);
				}
				int xj = splitTree->prev[xi];
				remove(xi, joinTree);
				remove(xi, splitTree);
//				int fr = xj;
//				int to = xi;
//				rg->addArc(fr, to);
				added ++;
				if(splitTree->nsize[xj] + joinTree->psize[xj] == 1) {
					q[back ++] = xj;
				}
			} else {
				if(joinTree->nsize[xi] > 1) {
					println("Can this happen too???");
					exit(0);
				}
				if(joinTree->nsize[xi] == 0) {
					println("Can this happen too again???");
					exit(0);
				}
				int xj = joinTree->next[xi];
				
				remove(xi, joinTree);
				remove(xi, splitTree);
				added ++;
//				int fr = xi;
//				int to = xj;
//				rg->addArc(fr, to);

				if(splitTree->nsize[xj] + joinTree->psize[xj] == 1) {
					q[back ++] = xj;
				}
			}
		}
	}


	void println(const char * s) {
		cout << s << endl;
	}

};



#endif
