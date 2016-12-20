#ifndef CONST_H_
#define CONST_H_

enum NODE_TYPE {
	UNDEFINED = 0, REGULAR=1, MINIMUM=2, MAXIMUM=3, SADDLE=4, ERRORNODE=5
};

enum SAD {
	FS_XY = 1, FS_YZ=2, FS_XZ=3, BS_1=4, BS_2=5
};

int MAX_ADJ;
bool grid;

#define null NULL

float * function_values;
int * vertex_pos;
char * offset;

bool lessVertex(int i, int j) {
	if(function_values[i] < function_values[j]) {
		return true;
	}
	if(function_values[i] == function_values[j]) {
		if(vertex_pos[i] < vertex_pos[j]) {
			return true;
		}
		if(vertex_pos[i] == vertex_pos[j]) {
            // offset = 0 if critical point is a maximum/minimum/boundary saddle
            // offset = 1 if critical point is face/body saddle
			if(offset[i] < offset[j]) {
				return true;
			}
			if(offset[i] == offset[j]) {
                // if the offset are equal, then just check the indices
				if(i < j) {
					return true;
				}
			}
		}
	}
	return false;
}

bool greaterVertex(int i, int j) {
	return !(lessVertex(i, j));
}

inline int compVertex(const void *pi, const void *pj) {
	const int i=*((const int *)pi);
	const int j=*((const int *)pj);

	if(function_values[i] != function_values[j]) {
		return ((function_values[i]>function_values[j])?1:-1);
	}
	else {
		if(vertex_pos[i] != vertex_pos[j]) {
			return (vertex_pos[i]>vertex_pos[j]?1:-1);
		}
		else {
			if(offset[i] != offset[j]) {
				return (offset[i]>offset[j]?1:-1);
			}
			else {
				if(i != j) {
					return (i>j?1:-1);
				}
			}
		}
	return 0;
	}
return 0;
}

struct compare {
	bool operator()(int a, int b) {
		return lessVertex(a,b);
	}
};

#endif /* CONST_H_ */
