#ifndef STOREREEBGRAPH_H
#define STOREREEBGRAPH_H

#include "MyIntList.h"
#include "CleanReebGraph.h"

class StoreReebGraph {

public :
	int * nodes;
	int * nodeType;
	float * nodeFn;
	MyIntList * xarc;
	MyIntList * yarc;

	StoreReebGraph(int no) {
		nodes = new int[no];
		nodeType = new int[no];
		nodeFn = new float[no];
		xarc = new MyIntList((int) (no*1.5));
		yarc = new MyIntList((int) (no*1.5));
		curNode = 0;
	}
	
	int curNode;
	void addNode(int v, float fn, int type) {
		nodes[curNode] = v;
		nodeType[curNode] = type;
		nodeFn[curNode] = fn;
		curNode ++;
	}
	
	CleanReebGraph rg;
	void setup() {
		 for(int i = 0;i < curNode;i ++) {
			 rg.addNode(nodes[i], nodeFn[i], nodeType[i]);
		 }
		 for(int i = 0;i < xarc->length;i ++) {
			 rg.addArc(xarc->array[i], yarc->array[i]);
		 }
		 rg.setup();
	}
	
	void addArc(int v1, int v2) {
		xarc->add(v1);
		yarc->add(v2);
	}
	void removeDeg2Nodes() {
		rg.removeDeg2Nodes();
	}
	
	
	void outputReebGraph(const char * file) {
		rg.outputReebGraph(file);
	}

};


#endif
