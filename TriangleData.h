#ifndef TRIANGLEDATA_H
#define TRIANGLEDATA_H

#include <iostream>
#include <list>
#include <algorithm>
#include <time.h>
#include <deque>
#include <iostream>

#include "MyIntList.h"

using namespace std;

#define less(v1, v2) ((fnVertices[v1] < fnVertices[v2]) || (fnVertices[v1] == fnVertices[v2] && v1 < v2))
#define PREV 0
#define NEXT 1
#define max(a,b) ((a < b)?b:a)
#define getIndex(tin, e, n) ((tin * 6 + e * 2 + n))

typedef struct Vertex {
	MyIntList star;
} Vertex;

class TriangleData {
public:


	float * fnVertices;
	int noVertices;
	int vertexCt;
	Vertex * vertices;

	MyIntList tris;
	int * triangles;
	MyIntList adjList;
	int triCt;

	int maxStar;
	
	bool storeCoords;
	
	float *x, *y, *z;
	
	TriangleData(bool storeCoords) {
		this -> storeCoords = storeCoords;
	}
	
	void addVertex(float xx, float yy, float zz, float ff) {
		if(storeCoords) {
			x[vertexCt] = xx;
			y[vertexCt] = yy;
			z[vertexCt] = zz;
		}
		fnVertices[vertexCt++] = ff;
	}

	void addTriangle(int v1, int v2, int v3) {
		int v[3];
		if (less(v1, v2)) {
			if (less(v1, v3)) {
				if (less(v2, v3)) {
					v[0] = v1;
					v[1] = v2;
					v[2] = v3;
				} else {
					v[0] = v1;
					v[1] = v3;
					v[2] = v2;
				}
			} else {
				v[0] = v3;
				v[1] = v1;
				v[2] = v2;
			}
		} else {
			if (less(v2, v3)) {
				if (less(v1, v3)) {
					v[0] = v2;
					v[1] = v1;
					v[2] = v3;
				} else {
					v[0] = v2;
					v[1] = v3;
					v[2] = v1;
				}
			} else {
				v[0] = v3;
				v[1] = v2;
				v[2] = v1;
			}
		}
		v1 = v[0];
		v2 = v[1];
		v3 = v[2];

		if(isAdded(v1,v2,v3)) {
			return;
		}
		tris.add(v1);
		tris.add(v2);
		tris.add(v3);
		
		vertices[v1].star.add(triCt);
		vertices[v2].star.add(triCt);
		vertices[v3].star.add(triCt);
		maxStar = max(maxStar, vertices[v1].star.length);
		maxStar = max(maxStar, vertices[v2].star.length);
		maxStar = max(maxStar, vertices[v3].star.length);
		triCt++;
	}

	
	bool isAdded(int v1, int v2, int v3) {
		int s1 = vertices[v1].star.size();
		int s2 = vertices[v2].star.size();
		int s3 = vertices[v3].star.size();

		if (s1 < s2 && s1 < s3) {
			return isPresent(v1, v1, v2, v3);
		} else if (s2 < s3) {
			return isPresent(v2, v1, v2, v3);
		} else {
			return isPresent(v3, v1, v2, v3);
		}
	}

	bool isPresent(int v, int v1, int v2, int v3) {
		int s = vertices[v].star.size();
		for (int i = 0; i < s; i++) {
			int t = vertices[v].star.get(i);
			t *= 3;
			int tv1 = tris.get(t);
			int tv2 = tris.get(t + 1);
			int tv3 = tris.get(t + 2);
			if (tv1 == v1 && tv2 == v2 && tv3 == v3) {
				return true;
			}
		}
		return false;
	}

	int noTris;
	void setupTriangles() {
		noTris = tris.length / 3;
		cout << "No. of triangles : " << noTris << endl;
		triangles = new int [tris.length];
		memcpy(triangles, tris.array, tris.length * sizeof(int));
		tris.clear();
	}

	void setNoOfVertices(int nv) {
		noVertices = nv;
		fnVertices = new float[noVertices];
		vertices = new Vertex[noVertices];
		
		if(storeCoords) {
			x = new float[nv];
			y = new float[nv];
			z = new float[nv];
		}
		vertexCt = 0;
		triCt = 0;
		maxStar = 0;
	}
	
	MyIntList fan;
	MyIntList getTriangleFan(int v1, int v2) {
		fan.length = 0;
		int s1 = vertices[v1].star.length;
		int s2 = vertices[v2].star.length;
		if(s1 < s2) {
			return getFan(vertices[v1], v1, v2);
		} else {
			return getFan(vertices[v2], v1, v2);
		}
	}

	MyIntList getFan(Vertex v, int v1, int v2) {
		for(int i = 0;i < v.star.length;i ++) {
			int t = v.star.array[i];
			int tin = t * 3;
			int tv1 = triangles[tin];
			int tv2 = triangles[tin + 1];
			int tv3 = triangles[tin + 2];
			if(tv1 == v1) {
				if(tv2 == v2 || tv3 == v2) {
					fan.add(t);
				}
			} else if(tv2 == v1) {
				if(tv3 == v2) {
					fan.add(t);
				}
			}
		}
		return fan;
	}
};


#endif
