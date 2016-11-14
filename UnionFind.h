/*
 * UnionFind.h
 *
 *  Created on: 13-May-2012
 *      Author: harish
 */

#ifndef UNIONFIND_H_
#define UNIONFIND_H_


struct Elem {
	int parent;
	int rank;
};

void ElemSet(struct Elem* elems, int i) {
	elems[i].parent = i;
	elems[i].rank = 0;
}

struct Element {
	Element *parent;
	int rank;
	int vertex;

	Element(int v) {
		parent = this;
		rank = 0;
		vertex = v;
	}
};

Element* find(Element* elem) {
	if (elem->parent == elem) {
		return elem;
	} else {
		return elem->parent = find(elem->parent);
	}
}

void uni(Element *x, Element *y) {
	if (x == y)
		return;

	Element *x_root = find(x);
	Element *y_root = find(y);

	if (x_root == y_root)
		return;

	if (x_root->rank > y_root->rank) {
		y_root->parent = x_root;
	} else if (y_root->rank > x_root->rank) {
		x_root->parent = y_root;
	} else {
		y_root->parent = x_root;
		x_root->rank++;
	}

}


#endif /* UNIONFIND_H_ */
