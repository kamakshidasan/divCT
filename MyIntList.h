#ifndef MYINTLIST_H
#define MYINTLIST_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

class MyIntList {
public:
	int * array;
	int length;
	int arraylength;
	
	MyIntList(int no) {
		array = new int[no];
		length = 0;
		arraylength = no;
	}
	
	MyIntList() {
		array = new int[10];
		length = 0;
		arraylength = 10;
	}
	
	void add(int n) {
		if(array == NULL) {
			array = new int[10];
			arraylength = 10;
		}
		if(length == arraylength) {
			int newl = int(length * 1.5);
			int * array1 = new int [newl];
			memcpy(array1, array, length*sizeof(int));
			delete array;
			array = array1;
			arraylength = newl;
		}
		array[length ++] = n;
	}
	
	int size() {
		return length;
	}
	
	int get(int i) {
		return array[i];
	}

	void clear() {
		delete array;
		array = NULL;
		length = 0;
	}

	void set(int pos, int val) {
		array[pos] = val;		
	}
};

#endif
