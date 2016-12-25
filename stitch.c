#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <omp.h>
#include "Graph.h"

int divx = 2; //no. of divisions along x-axis
int divy = 2;
int  divz =  2;
int sub_sizex = 128; //subcube size
int sub_sizey = 128;
int sub_sizez = 128;
int dimx = 256;
int dimy = 256;
int dimz = 256;
struct timeval tv1, tv2;
struct Graph* CTGraph;
int ext[6];
int dim;
long long int * vp;
float * fv;
char * off;
int * p;
int * c;
int* p_;
int* c_;
int * U;
int *U_;
int * own;
int * temp;
int * vi;
int * j_n, *j_c, *s_n, *s_c, *v_i;
long long int* v_p;
char *offset;
float *f_v;
int n1;
int n2;
int vn1, vn2, b1, b2;
int final;
int critical;

enum {
	false, true
};

void Union(int child, int parent, int * Uk) {
	if (child == parent)
		return;
	Uk[child] = parent;
	Uk[parent] = -1;
}


int Find(int elem, int* Uk) {
	if (Uk[elem] == -2) {
		return -2;
	}

	int t_elem = elem;
	int tt_elem;
	while (Uk[elem] != -1) {
		elem = Uk[elem];
	}

	//path compression

	while ((Uk[t_elem] != -1)) {
		tt_elem = t_elem;
		t_elem = Uk[t_elem];
		Uk[tt_elem] = elem;
	}
	//printf("find\n");
	return elem;
}

int lessVertex(int i, int j) {
	
	//printf("\nless vertex\n");
	if (fv[i] < fv[j]) {
		return true;
	}
	if (fv[i] == fv[j]) {
		if (vp[i] < vp[j]) {
			return true;
		}
		if (vp[i] == vp[j]) {
			if (off[i] < off[j]) {
				return true;
			}
			if (off[i] == off[j]) {
				if (i < j) {
					return true;
				}
			}
		}
	}
	return false;
}

// Standard code to merge two sorted sequence of elements
void seqmerge() {

	int * a = vi;
	int size = n1 + n2;

	int i, m, k, l;

	int mid = n1 - 1;
	int low = 0;
	int high = size - 1;
	l = low;
	i = low;
	m = mid + 1;

	while ((l <= mid) && (m <= high)) {
		if (lessVertex(a[l], a[m])) {
			temp[i] = a[l];
			own[i] = 0;
			l++;
		}
		else {
			temp[i] = a[m];
			own[i] = 1;
			m++;
		}
		i++;
	}

	if (l > mid) {
		for (k = m; k <= high; k++) {
			temp[i] = a[k];
			own[i] = 1;
			i++;
		}
	}
	else {
		for (k = l; k <= mid; k++) {
			temp[i] = a[k];
			own[i] = 0;
			i++;
		}
	}

	//free(a);
	a = temp;
}

void cleanup_trees() {

	int i, j, t;

	int num_critical = n1 + n2;
	int * x = (int*) malloc(num_critical * sizeof(int));
	char * is_critical = (char*) malloc(num_critical * sizeof(char));
	v_i = (int*) malloc(num_critical * sizeof(int));

	critical = 0;

	for (i = 0; i < num_critical; i++) {
		j = temp[i]; // That sorted array after seqmerge()
		is_critical[j] = 0;

		long long int vpos = vp[j];
		long long int dimxy = dimx * dimy;
		long long int posz = vpos / dimxy;
		long long int xy = vpos % dimxy;
		long long int posy = xy / dimx;
		long long int posx = xy % dimx;

            // If node is a leaf
		if ((c[2 * j] == -1 && c[2 * j + 1] == -1) || (c_[2 * j] == -1 && c_[2 * j + 1] == -1) ||
            // If both children in either trees exist
            (c[2 * j + 1] != -1 && c[2 * j] != -1) || (c_[2 * j + 1] != -1 && c_[2 * j] != -1) ||
            // Lies on the x boundary
            (posx == ext[0]) || (posx == ext[1]) ||
            // Lies on the y boundary
            (posy == ext[2]) || (posy == ext[3]) ||
            // Lies on the z boundary
            (posz == ext[4]) || (posz == ext[5])) {
                x[j] = critical;
                is_critical[j] = 1;
                v_i[critical] = critical;
                critical++;
		}
	}

    // Create arrays!
	j_n = (int*) malloc(critical * sizeof(int));
	s_n = (int*) malloc(critical * sizeof(int));
	j_c = (int*) malloc(2 * critical * sizeof(int));
	s_c = (int*) malloc(2 * critical * sizeof(int));
	v_i = (int*) realloc(v_i,critical * sizeof(int));
	v_p = (long long int*) malloc(critical * sizeof(long long int));
	f_v = (float*) malloc(critical * sizeof(float));
	offset = (char*) malloc(critical * sizeof(char));

    // Initialize the necessary arrays
	#pragma omp parallel for schedule(dynamic) private(i)
	for (i = 0; i < critical; i++) {
		j_c[2 * i] = j_c[2 * i + 1] = s_c[2 * i] = s_c[2 * i + 1] = j_n[i] = s_n[i] = -1;
	}

	int l = 0; //leaf count

	#pragma omp parallel for schedule(dynamic) private(i,j,t)
	for (i = 0; i < num_critical; i++) {
		j = temp[i];

		if (is_critical[j] == 1) {
            		// Update the function value
			f_v[x[j]] = fv[j];
			// Update the vertex position
			v_p[x[j]] = vp[j];
			// Update the offset
			offset[x[j]] = off[j];

            		/* Clean the Join Tree */
            		// Get the parent
           		t = p[j];

            		// If parent is not critical, then join with a grandparent
			while ((t != -1) && (is_critical[t] != 1)) {
				t = p[t];
			}

            		// Parent does not exist
			if (t == -1) {
				j_n[x[j]] = -1;
			}

			// Adhitya: Should this be parent?
			// The newly found parent is critical, and therefore is the new join neighbour
			else {
				j_n[x[j]] = x[t];
                		// Mark this node as the child of the new parent
                		// If the first child has not been assigned, then take it's place :)
				if (j_c[2 * x[t]] == -1) {
					j_c[2 * x[t]] = x[j];
                		}
                		// Otherwise, the second one is up for grabs
				else {
					j_c[2 * x[t] + 1] = x[j];
				}
			}

            		/* Clean the Split Tree */
           		// Time for Pruning the Split Tree!
			t = p_[j];

			// If Split Parent is not critical, then split with grandparent
			while ((t != -1) && (is_critical[t] != 1)) {
				t = p_[t];
            		}

            		// Split parent does not exist
			if (t == -1) {
				s_n[x[j]] = -1;

			}

			// The newly found parent is critical, and therefore is the new split neighbour
			else {
				s_n[x[j]] = x[t];
                		// Mark this node as the child of the new parent
                		// If the first child has not been assigned, then take it's place :)
				if (s_c[2 * x[t]] == -1) {
					s_c[2 * x[t]] = x[j];
				}
                		// Otherwise, the second one is up for grabs
				else {
					s_c[2 * x[t] + 1] = x[j];
                		}
			}
		}
	}

    	// Critical points after cleaning
	final = critical;

	// Cleaning complete - Free all the arrays already!
	free(p);
	free(p_);
	free(c);
	free(c_);

	free(vp);
	free(fv);
    	free(x);
    	free(is_critical);
	//printf("cleanup tree: critical: %d. critical_prev: %d\n", critical, num_critical);

}

int min(int x, int y) {
    return x < y ? x : y;
}

int max(int x, int y) {
    return x < y ?  y: x;
}

int main(int argc, char **argv) {
	//omp_set_num_threads(1);

	// Process the params file
	FILE* param = fopen("params.txt","r");
	fscanf(param,"%d %d %d %d %d %d",&divx,&divy,&divz,&sub_sizex,&sub_sizey,&sub_sizez);
	fclose(param);
	dimx = divx*sub_sizex;
	dimy = divy*sub_sizey;
	dimz = divz*sub_sizez;
	
	// Adhitya: For the final stitch, give the argument value as 1
	int tree;
	if (argc == 5) {
		tree = atoi(argv[4]);
	}

	// Process the files produced by the earlier program
	char file1[30], file2[30], file[30];
	strcpy(file1, argv[1]);
	strcat(file1, ".dat");
	strcpy(file2, argv[2]);
	strcat(file2, ".dat");
	strcpy(file, argv[3]);
	strcat(file, ".dat");

	// Store size of integer, char, float
	size_t isz = sizeof(int);
	size_t csz = sizeof(char);
	size_t fsz = sizeof(float);

	int ext1[6], ext2[6];
	int fp1, fp2;

	// Open both the input files
	fp1 = open(file1, O_RDONLY);
	fp2 = open(file2, O_RDONLY);

	gettimeofday(&tv1, NULL);

	// Look how data is stored in the files towards the end of the previous program

	// Store the number of vertices
	read(fp1, (void*) (&vn1), isz); //fscanf(fp1,"%d",&n1);
	read(fp2, (void*) (&vn2), isz); //fscanf(fp2,"%d",&n2);

	// Store the extents
	read(fp1, (void*) (ext1), 6 * isz);
	read(fp2, (void*) (ext2), 6 * isz);

	//extents
	ext[0] = min(ext1[0], ext2[0]); //find out the extents
	ext[1] = max(ext1[1], ext2[1]);
	ext[2] = min(ext1[2], ext2[2]);
	ext[4] = min(ext1[4], ext2[4]);
	ext[3] = max(ext1[3], ext2[3]);
	ext[5] = max(ext1[5], ext2[5]);

	// Store the number of critical points
	read(fp1, (void*) (&n1), isz);
	read(fp2, (void*) (&n2), isz);
	printf("n1:%d,n2:%d\n", n1, n2);
	
	// Store the function values
	fv = (float*) malloc((n1 + n2) * sizeof(float));
	read(fp1, (void*) fv, n1 * fsz);
	read(fp2, (void*) (fv + n1), n2 * fsz);

	// Store vertex_index
	vi = (int*) malloc((n1 + n2) * sizeof(int));
	read(fp1, (void*) vi, n1 * isz);
	read(fp2, (void*) (vi + n1), n2 * isz);
	
	// Store vertex_pos
	vp = (long long int*) malloc((n1 + n2) * sizeof(long long int));
	read(fp1, (void*) vp, n1 * sizeof(long long int));
	read(fp2, (void*) (vp + n1), n2 * sizeof(long long int));
	
	// Store the join tree parent array
	p = (int*) malloc((n1 + n2) * sizeof(int));
	read(fp1, (void*) p, n1 * isz);
	read(fp2, (void*) (p + n1), n2 * isz);
	
	// Store the join tree children
	c = (int*) malloc(2 * (n1 + n2) * sizeof(int));
	read(fp1, (void*) c, 2 * n1 * isz);
	read(fp2, (void*) (c + 2 * n1), 2 * n2 * isz);
	
	// Store the split tree parents
	p_ = (int*) malloc((n1 + n2) * sizeof(int));
	read(fp1, (void*) p_, n1 * isz);
	read(fp2, (void*) (p_ + n1), n2 * isz);
	
	// Store the split tree children
	c_ = (int*) malloc(2 * (n1 + n2) * sizeof(int));
	read(fp1, (void*) c_, 2 * n1 * isz);
	read(fp2, (void*) (c_ + 2 * n1), 2 * n2 * isz);
	
	// Store the offset - to indicate if a point is face or body saddle
	off = (char*) malloc((n1 + n2) * sizeof(char));
	read(fp1, (void*) off, n1 * csz);
	read(fp2, (void*) (off + n1), n2 * csz);

	close(fp1);
	close(fp2);
	
	// Adhitya: Uncomment the following lines, before release
	//remove(file1);
	//remove(file2);
	
	// Union Find for Join and Split Trees
	U = (int*) malloc((n1 + n2) * sizeof(int));
	U_ = (int*) malloc((n1 + n2) * sizeof(int));

	own = (int*) malloc((n1 + n2) * sizeof(int));
	temp = (int*) malloc((n1 + n2) * sizeof(int));

	gettimeofday(&tv2, NULL);

	printf("\n alloc and read time = %f miliseconds\n",(double) (tv2.tv_usec - tv1.tv_usec) / 1000+ (double) (tv2.tv_sec - tv1.tv_sec) * 1000);


	int i, j, k;

	gettimeofday(&tv1, NULL);

	// Merge both trees together

	#pragma omp parallel for schedule(dynamic) private(i)
	for (i = 0; i < (n1 + n2); i++) {
		U[i] = -2;
		U_[i] = -2;

		if (i >= n1) {

            		// Process vertex indices
			vi[i] = n1 + vi[i];

           		// Process parents
			p[i] = (p[i] == -1) ? -1 : n1 + p[i];
			p_[i] = (p_[i] == -1) ? -1 : n1 + p_[i];

            		// Process children
			c[2 * i] = (c[2 * i] == -1) ? -1 : (n1 + c[2 * i]);
			c_[2 * i] = (c_[2 * i] == -1) ? -1 : (n1 + c_[2 * i]);

			c[2 * i + 1] = (c[2 * i + 1] == -1) ? -1 : (n1 + c[2 * i + 1]);
			c_[2 * i + 1] = (c_[2 * i + 1] == -1) ? -1 : (n1 + c_[2 * i + 1]);
		}
	}

	// Phew, finally better names to work with!
	critical = n1 + n2;
	j_n = p;
	s_n = p_;
	j_c = c;
	s_c = c_;

	v_p = vp;
	f_v = fv;
	offset = off;

	// Merge the vertex_indices
	// Duplicate points will appear next to each other now
	seqmerge();
	v_i = temp;

	int * S = temp;

    #pragma omp parallel sections private(i,k)
    {
        // Join Tree

        #pragma omp section
        {
            for (i = 1; i < n1 + n2; i++) {
                int z = S[i];
                int zm = S[i-1];

                // Check if z and zm are boundary duplicate points
                if ((vp[z] == vp[zm]) && (off[z] == 0) && (off[zm] == 0) && (own[i] != own[i - 1])) {
                    // Make z as parent of zm
                    p[zm] = z;
                    // Add zm to z's children list
                    if (c[2 * z] == -1) {
                        c[2 * z] = zm;
                    }
                    else {
                        // Adhitya: What happens to the second child in this case?
                        if (c[2 * z] != zm) {
                            c[2 * z + 1] = zm;
                        }
                    }
                    //printf("dup:%d\n",i);
                    Union(zm, z, U);
                }

                // For the first child of z
                if (c[2 * z] != -1) {
                    // If child is present in U
                    if (U[c[2 * z]] != -2) {
                        // Find the parent from U
                        k = Find(c[2 * z], U);
                        // If new parent and current parent are not the same
                        if (k != z) {
                            Union(k, z, U);
                            // Make z as parent of k
                            p[k] = z;
                            // Make k as child of z
                            c[2 * z] = k;

                            // If the second child is k, then it should not be repeated
                            if (c[2 * z + 1] == k) {
                                c[2 * z + 1] = -1;
                            }
                        }
                    }
                }

                // For the first child of z
                if (c[2 * z + 1] != -1) {
                    // If child is present in U
                    if (U[c[2 * z + 1]] != -2) {
                        // Find the parent from U
                        k = Find(c[2 * z + 1], U);
                        // If new parent and current parent are not the same
                        if (k != z) {
                            Union(k, z, U);
                            // Make z as parent of k
                            p[k] = z;
                            // Make k as child of z
                            if (c[2 * z] != k) {
                                c[2 * z + 1] = k;
                            }
                            // If first child is already k, then avoid duplication
                            else {
                                c[2 * z + 1] = -1;
                            }
                        }
                    }
                }
            }
            free(U);
	} // omp section 1 ends

        // Split Tree

        #pragma omp section
        {
            // This time go down the sorted vertices because we are computing a split tree
            for (i = (n1 + n2 - 2); i >= 0; i--) {
                int z = S[i];
                int zp = S[i+1];

                // If first child of z exists
                if (c_[2 * z] != -1) {
                    // If the child is present in the union
                    if (U_[c_[2 * z]] != -2) {
                        // Find the parent from U_
                        k = Find(c_[2 * z], U_);
                        // If both new parent and current parent are not the same
                        if (k != z) {
                            Union(k, z, U_);
                            // Make parent of k as z
                            p_[k] = z;
                            // Make k as child of z
                            c_[2 * z] = k;
                            // If the second child is k, then avoid duplication
                            if (c_[2 * z + 1] == k) {
                                c_[2 * z + 1] = -1;
                            }
                        }
                    }
                }

                // If second child of z exists
                if (c_[2 * z + 1] != -1) {
                    // If child is present in U_
                    if (U_[c_[2 * z + 1]] != -2) {
                        // Find the parent from U_
                        k = Find(c_[2 * z + 1], U_);
                        // If both new parent and current parent are not the same
                        if (k != z) {
                            Union(k, z, U_);
                            // Make parent of k as z
                            p_[k] = z;
                            // Make k as child of z
                            if (c_[2 * z] != k) {
                                c_[2 * z + 1] = k;
                            }
                            // If first child is already k, then avoid duplication
                            else {
                                c_[2 * z + 1] = -1;
                            }
                        }
                    }
                }

                // Check if z and zm are boundary duplicate points
                if ((vp[z] == vp[zp]) && (off[z] == 0) && (off[zp] == 0) && (own[i] != own[i + 1])) {
                    // Make zp as parent of z
                    p_[z] = zp;
                    // Add z to zp's children list
                    if (c_[2 * zp] == -1) {
                        c_[2 * zp] = z;
                    }
                    else {
                        if (c_[2 * zp] == z) {
                            c_[2 * zp + 1] = z;
                        }
                    }
                    Union(z, zp, U_);
                }
            }
            free(U_);
        } // omp section 2 ends
    }

    free(own);

	gettimeofday(&tv2, NULL);

	printf("\n stitch time = %f miliseconds\n",	(double) (tv2.tv_usec - tv1.tv_usec) / 1000	+ (double) (tv2.tv_sec - tv1.tv_sec) * 1000);

	gettimeofday(&tv1, NULL);

	if (tree == 1) {
		cleanup_trees();
    	}

	gettimeofday(&tv2, NULL);

	printf("\n cleanup time = %f miliseconds\n",(double) (tv2.tv_usec - tv1.tv_usec) / 1000	+ (double) (tv2.tv_sec - tv1.tv_sec) * 1000);

	// Write everything stitched to a file
	int fp = open(file, O_RDWR | O_CREAT | O_TRUNC, 0777);
	int num_critical = critical;
	int num_vert = 0;//vn1 + vn2 - b1;
	printf("\n no. of critical points: %d\n",critical);
	gettimeofday(&tv1, NULL);
	write(fp, (void*) (&num_vert), isz);
	write(fp, (void*) (ext), 6*isz);
	write(fp, (void*) (&critical), isz);
	write(fp, (void*) (f_v), num_critical * fsz);
	write(fp, (void*) (v_i), num_critical * isz);
	write(fp, (void*) (v_p), num_critical * sizeof(long long int));
	write(fp, (void*) (j_n), num_critical * isz);
	write(fp, (void*) (j_c), 2 * num_critical * isz);
	write(fp, (void*) (s_n), num_critical * isz);
	write(fp, (void*) (s_c), 2 * num_critical * isz);
	write(fp, (void*) (offset), num_critical * csz);

	close(fp);

	gettimeofday(&tv2, NULL);
	//sync(fp);
	printf("\n writing time = %f miliseconds\n",(double) (tv2.tv_usec - tv1.tv_usec) / 1000	+ (double) (tv2.tv_sec - tv1.tv_sec) * 1000);

	if (tree != 1) {
		return 0;
	}

	int *upper_degree = (int*) malloc(num_critical * sizeof(int));
	int *lower_degree = (int*) malloc(num_critical * sizeof(int));
	int *leaves = (int*) malloc(num_critical * sizeof(int));
	int *is_processed = (int*) malloc(num_critical * sizeof(int));

	int b, num_leaves = 0;
	// Initialize a Contour Tree
	CTGraph = createGraph(num_critical);
	gettimeofday(&tv1, NULL);
	
	// Initialize the arrays
	#pragma omp parallel for schedule(dynamic) private(b)
	for (b = 0; b < num_critical; b++) {
		upper_degree[b] = 0;
		lower_degree[b] = 0;
		is_processed[b] = 0;
	}

	// Look at Hamish Carr's Algorithm 4.2
	// Find the number of leaves
	#pragma omp parallel for schedule(dynamic) private(b)
	for (b = 0; b < num_critical; b++) {
        // lower_degree[j_n[b]]++;
		if (j_n[b] > -1) {
            		__sync_fetch_and_add(lower_degree + j_n[b], 1);
		}
		// upper_degree[s_n[b]]++;
		if (s_n[b] > -1) {
            		__sync_fetch_and_add(upper_degree + s_n[b], 1);
		}

        	// Check if node is a leaf in either of the trees
		if (((j_c[2 * b] == -1) && (j_c[2 * b + 1] == -1)) || ((s_c[2 * b] == -1) && (s_c[2 * b + 1] == -1))){
			leaves[num_leaves] = b;
			__sync_fetch_and_add(&num_leaves,1);
		}
	}


	#pragma omp parallel for schedule(dynamic) private(i)
	for (i = 0; i < (num_leaves); i++) {
		int lcnt,ucnt,j,u,l;
		int tid = omp_get_thread_num();
		//int i = get_global_id(0);
		j = leaves[i]; //ith leaf
		while (!is_processed[j]) { //TODO check

			//is_processed[j]++;
			__sync_fetch_and_add(is_processed+j,1);

            		//upper degree in join Y tree and lower degree in split Y' tree
			if (lower_degree[j] == 0 && upper_degree[j] == 1) {
				u = j_n[j];
				while (is_processed[u] && j_n[u] != -1) {
					u = j_n[u];
                		}
				#pragma omp critical
				{
                    			addEdge(CTGraph,u,j);
                		}
				if (is_processed[u]) {
					break;
				}

				//lcnt = --lower_degree[u];
				lcnt = __sync_sub_and_fetch (lower_degree+u,1);//remove the vertex
				ucnt = upper_degree[u];

                		//add (u,j)
				if ((lcnt == 0 && ucnt == 1) || (lcnt == 1 && ucnt == 0)) {
					j = u;
				}
				else {
					break;
				}
			}

			// lower degree in join Y tree and upper degree in split Y' tree
            		else if (lower_degree[j] == 1 && upper_degree[j] == 0) {
				l = s_n[j];
				while (is_processed[l] && s_n[l] != -1) {
					l = s_n[l];
                		}
				#pragma omp critical
				{
                    			addEdge(CTGraph,l,j);
                		}

                		// Don't process it again!
				if (is_processed[l]) {
					break;
                		}
				lcnt = lower_degree[l];
				//ucnt = --upper_degree[l];
				ucnt = __sync_sub_and_fetch (upper_degree+l,1);

				//add (u,j)
				if ((lcnt == 0 && ucnt == 1) || (lcnt == 1 && ucnt == 0)) {
					j = l;
				}
				else {
					break;
				}
			}
			else {
				break;
			}
		}
	}

	printf("\nContour Tree computed\n");
    	gettimeofday(&tv2, NULL);

	printf("\n tree merge time = %f milliseconds\n",(double) (tv2.tv_usec - tv1.tv_usec) / 1000	+ (double) (tv2.tv_sec - tv1.tv_sec) * 1000);

	//-------------------------------------uncomment to write final contour tree to file-----------------
	printf("\nwriting final contour tree\n");
	printGraph(CTGraph,f_v);
    //--------------------------------------------------------------------------------
	return 0;
}


