#ifndef READFILE_H
#define READFILE_H

#include<iostream>
#include<fstream>
#include "TriangleData.h"

using namespace std;

int readInt(char * memBlock, int &pos) {
	int * beg = (int *)memBlock;
	int p = beg[pos];
	pos ++;
	return p;
}

float readFloat(char * memBlock, int &pos) {
	float * beg = (float *)memBlock;
	float p = beg[pos];
	pos ++;
	return p;
}

TriangleData * readTSFile(const char * fileName, int fin, bool adj) {
	ifstream f (fileName);
	int nv, nt;
	f >> nv >> nt;
	cout << "No. of vertices : " << nv << endl;
	cout << "No. of Tetrahedron : " << nt << endl;
	TriangleData * data = new TriangleData(adj);
	data -> setNoOfVertices(nv);
	for(int i = 0;i < nv;i ++) {
		float x[4];
		f >> x[0] >> x[1] >> x[2] >> x[3];
		if(fin != 0) {
			x[3] = x[fin - 1];
		}
		data -> addVertex(x[0],x[1],x[2],x[3]);
	}

	for(int i = 0;i < nt;i ++) {
		int v1, v2, v3, v4;
		f >> v1 >> v2 >> v3 >> v4;
		data -> addTriangle(v1, v2, v3);
		data -> addTriangle(v1, v2, v4);
		data -> addTriangle(v1, v3, v4);
		data -> addTriangle(v2, v3, v4);
	}
	f.close();
	data -> setupTriangles(); 
	return data;
}


TriangleData * readOFFFile(const char * fileName, int fin, bool adj) {
	ifstream f (fileName);
	int nv, nt;
	f >> nv >> nt;
	cout << "No. of vertices : " << nv << endl;
	cout << "No. of Triangles : " << nt << endl;
	TriangleData * data = new TriangleData(adj);
	data -> setNoOfVertices(nv);
	for(int i = 0;i < nv;i ++) {
		float x[4];
		f >> x[0] >> x[1] >> x[2] >> x[3];
		if(fin != 0) {
			x[3] = x[fin - 1];
		}
		data -> addVertex(x[0],x[1],x[2],x[3]);
	}

	for(int i = 0;i < nt;i ++) {
		int v1, v2, v3;
		f >> v1 >> v2 >> v3;
		data -> addTriangle(v1, v2, v3);
	}
	f.close();
	data -> setupTriangles(); 
	return data;
}

TriangleData * readSimFile(const char * fileName, int fin, bool adj) {
	fstream f(fileName, ios::in|ios::binary|ios::ate);
	int size;
	size = (int) f.tellg();
	char * memblock = new char[size];
	f.seekg (0, ios::beg);
	f.read (memblock, size);
	f.close();

	int ptr = 0;
	int nv = readInt(memblock, ptr);
	int dim = readInt(memblock, ptr);

	cout << "No. of vertices : " << nv << endl;
	cout << "Dimension: " << dim << endl;
	TriangleData * data = new TriangleData(adj);
	data -> setNoOfVertices(nv);

	for(int i = 0;i < nv;i ++) {
		float v[dim];
		for(int j = 0;j < dim;j ++) {
			v[j] = readFloat(memblock, ptr);
		}
		float fn = readFloat(memblock, ptr);
		if(fin != 0) {
			fn = v[fin - 1];
		}
		data -> addVertex(v[0],v[1],v[2],fn);
	}
	cout << "Finished reading vertices" << endl;
	int ns;
	int ct = 0;
	do{
		ns = readInt(memblock, ptr);
		if(ns != -1) {
			int v[ns];
			for(int i = 0;i < ns;i ++) {
				v[i] = readInt(memblock,ptr);
			}
			if(ns == 3) {
				data -> addTriangle(v[0],v[1],v[2]);
			} else if(ns == 4) {
				data -> addTriangle(v[0],v[1],v[2]);
				data -> addTriangle(v[0],v[1],v[3]);
				data -> addTriangle(v[0],v[2],v[3]);
				data -> addTriangle(v[1],v[2],v[3]);
			} else {
				cout << ct << " " << ns << endl;
				cout << "invalid input" << endl;
				exit(0);
			}
			ct ++;
		}
	}while(ns != -1);
	delete memblock;
	data -> setupTriangles();
//	exit(0);
	return data;
}


#endif
