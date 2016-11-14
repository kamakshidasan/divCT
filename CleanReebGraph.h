#ifndef CLEANREEBGRAPH_H
#define CLEANREEBGRAPH_H

#include <deque>
#include <map>
#include "const.h"

using namespace std;

class CleanReebGraph {

public:

	class Arc {
	public:
		int v1;
		int v2;
		float fn;
		short icol;
		bool removed;

		bool equals(Arc a) {
			return (a.icol == icol);
		}
	};

	class Node {
	public:

		int v;
		bool removed;
		float fn;
		int type;
		deque<Arc*> prev;
		deque<Arc*> next;
	};

	deque<Node*> nodes;
	deque<Arc*> arcs;
	map<int, int> vmap;

	int ct;
	short ect;
	float persistence;

	CleanReebGraph() {
		ct = 0;
		ect = 0;
		persistence = 0;
	}

	void addNode(int v, float fn, int type) {
		Node * n = new Node();
		n->fn = fn;
		n->v = v;
		n->type = type;
		nodes.push_back(n);
		vmap[v] = ct;
		ct ++;
	}

	void setup() {
	}

	void addArc(int v1, int v2) {
		Arc * a = new Arc();
		a->v1 = vmap[v1];
		a->v2 = vmap[v2];
		a->fn = nodes[a->v2]->fn - nodes[a->v1]->fn;

		ect ++;
		a->icol = ect;
		arcs.push_back(a);

		nodes[a->v1]->next.push_back(a);
		nodes[a->v2]->prev.push_back(a);
	}

	void removeDeg2Nodes() {
		// remove degree 2 vertices
		for(unsigned int i = 0;i < nodes.size();i ++) {
			if(nodes[i]->prev.size() == 0 && nodes[i]->next.size() == 0) {
				nodes[i]->removed = true;
			}
			if(!nodes[i]->removed && nodes[i]->next.size() == 1 && nodes[i]->prev.size() == 1) {
				mergeNode(i);
			}
		}
	}

	void mergeNode(int i) {
		Arc * e1 = nodes[i]->prev[0];
		Arc * e2 = nodes[i]->next[0];

		nodes[i]->removed = true;
		e1->v2 = e2->v2;
		e2->removed = true;
		if(e1->icol > e2->icol) {
			e1->icol = e2->icol;
		}
		e1->fn += e2->fn;
		int pos = -1;
		for(unsigned int x = 0;x < nodes[e1->v2]->prev.size();x ++) {
			if(nodes[e1->v2]->prev[x] == e2) {
				pos = x;
				break;
			}
		}
		if(pos == -1) {
			cout << "pos cannot be -1" << endl;;
		}
		nodes[e1->v2]->prev.erase(nodes[e1->v2]->prev.begin() + pos);
		nodes[e1->v2]->prev.push_back(e1);
	}

	void outputReebGraph(const char * file) {
		const char * str [5];
		str[REGULAR] = "REGULAR";
		str[MINIMUM] = "MINIMA";
		str[MAXIMUM] = "MAXIMA";
		str[SADDLE] = "SADDLE";

		ofstream of(file);
		int nv = 0;
		int ne = 0;
		for(unsigned int i = 0;i < nodes.size();i ++) {
			if(!nodes[i]->removed) {
				nv ++;
			}
		}

		for(unsigned int i = 0;i < arcs.size();i ++) {
			if(!arcs[i]->removed) {
				ne ++;
			}
		}

		of << nv << " " << ne << endl;
		for(unsigned int i = 0;i < nodes.size();i ++) {
			if(!nodes[i]->removed) {
				of << nodes[i]->v << " " << nodes[i]->fn << " " << str[nodes[i]->type] << endl;
			}
		}

		for(unsigned int i = 0;i < arcs.size();i ++) {
			if(!arcs[i]->removed) {
				of << nodes[arcs[i]->v1]->v << " " << nodes[arcs[i]->v2]->v  << endl;
			}
		}
		of.close();
	}
};

#endif
