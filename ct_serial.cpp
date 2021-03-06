#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <vector>
#include <utility>
#include <queue>
#include <deque>
#include <cassert>
#include <climits>
#include <cstdio>
#include <cstring>
#include <set>
#include <map>
#include <sys/time.h>
#include <omp.h>
#include <sched.h>
#include <fcntl.h>
#include <unistd.h>
#include "const.h"
#include "PSort.h"
using namespace std;
#define MAX_ADJ 2
int divx = 2; //no. of divisions along x-axis
int divy = 2;
int divz = 2;
int sub_sizex = 128; //subcube size
int sub_sizey = 128;
int sub_sizez = 128;

int F000 = 0;
int F001 = 4;
int F010 = 2;
int F011 = 6;
int F100 = 1;
int F101 = 5;
int F110 = 3;
int F111 = 7;

int XY = 1;     //00001
int YZ = 2;     //00010
int XZ = 4;     //00100
int XYZ = 8;    //01000
int XYZ2 = 16;  //10000

int global_dimx = divx*sub_sizex;
int global_dimy = divy*sub_sizey;
int global_dimz = divz*sub_sizez;


struct Graph* CTGraph;
char* offset_c;
int subcube;
long long int* g_vp;
int num_cores;
int final;
int num_vert, num_critical, num_minima, num_maxima, num_saddle, num_regular;
int func;
int dimx, dimy, dimz;
float * vertices;
int * vertex_index;
int * svi;
int * leaves;
int *lower_root_count;
int *upper_root_count;

NODE_TYPE* type;
int * bfType;
int * extrema_map;
int *steep_asc, *steep_desc, *min_root, *max_root;
int *lower_roots;
int *upper_roots;
int * critical_points;
int * join_neigh;
int * join_children;
int * split_neigh;
int * split_children;
int*j_n;
int*s_n;
int*j_c;
int*s_c;
int*v_i;
long long int *v_p;
float*f_v;
timeval start, end, result;
timeval start1, start2, end2, end1;
long mtime, seconds, useconds;
long cpTime, cleanuptime, initSetupTime, initParTime, rootTime, sortTime, initHeapTime,
buildHeapTime, joinTime, splitTime, completeTime, ctTime,
serialtreeTime;

struct extent {
    int x[2];
    int y[2];
    int z[2];
} ext;

struct Elems {
    int parent;
    int rank;
    int rep;
};

// The convention followed for representing neighbors
// If p has coordinates (x,y,z)
// 0 : (x-1,y,z) 1 : (x+1,y,z)
// 2 : (x,y-1,z) 3 : (x,y+1,z)
// 4 : (x,y,z-1) 5 : (x,y,z+1)

// Neighbors adjacent to each of the six vertices around a vertex
int neighbors[] = { 2, 3, 4, 5, // 0
                    2, 3, 4, 5, // 1
                    0, 1, 4, 5, // 2
                    0, 1, 4, 5, // 3
                    0, 1, 2, 3, // 4
                    0, 1, 2, 3  // 5
                  };

// The convention followed for representing vertices of a cube
// 0 : (0,0,0)  1 : (1,0,0)  2 : (0,1,0)  3 : (1,1,0)
// 4 : (0,0,1)  5 : (1,0,1)  6 : (0,1,1)  7 : (1,1,1)

// Adjacency List for each of the eight vertices in a cube
int bodyAdj[] = {   1, 2, 4, // 0
                    0, 3, 5, // 1
                    0, 3, 6, // 2
                    1, 2, 7, // 3
                    0, 5, 6, // 4
                    1, 4, 7, // 5
                    2, 4, 7, // 6
                    3, 5, 6  // 7
                };

// Insert the element
void ElemSet(struct Elems* elems, int i) {
    elems[i].parent = i;
    elems[i].rank = 0;
    elems[i].rep = i;
}

// Functions as advertised - Finds the parent
// Additionally, it also updates the parent of the nodes involved
int findit(struct Elems* elems, int i) {
    if (elems[i].parent == i){
        return i;
    }

    int k = i;

    // The node did not have it's parent as the same
    // Traverse all such nodes and find the node which has it's parent as same as it
    while (elems[k].parent != k) {
        k = elems[k].parent;
    }

    // Update the parent of entire traversal
    while (elems[i].parent != i) {
        int t = elems[i].parent;
        elems[i].parent = k;
        i = t;
    }

    return k;
}

// Given two nodes, unifies them - Updation of parents is done by findit()
// Rank based parent updation between the two nodes are done here
// Node with higher rank is returned
int unify(struct Elems* elems, int x, int y) {
    if (x == y) {
        return x;
    }

    // Find the parent of x
    int x_root = findit(elems, x);

    // Find the parent of y
    int y_root = findit(elems, y);

    // If both are equal then nothing to unify
    if (x_root == y_root) {
        return x_root;
    }

    // If x has a higher rank, update the parent of y
    if (elems[x_root].rank > elems[y_root].rank) {
        elems[y_root].parent = x_root;
        return x_root;
    }

    // If y has a higher rank, update the parent of x
    else if (elems[y_root].rank > elems[x_root].rank) {
        elems[x_root].parent = y_root;
        return y_root;
    }

    // If both of them have the same rank, update the parent of y
    else {
        elems[y_root].parent = x_root;
        elems[x_root].rank++;
        return x_root;
    }
}

void timeval_subtract(timeval *result, timeval *x, timeval *y) {
    // Perform the carry for the later subtraction by updating y.
    if (x->tv_usec < y->tv_usec) {
        int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
        y->tv_usec -= 1000000 * nsec;
        y->tv_sec += nsec;
    }
    if (x->tv_usec - y->tv_usec > 1000000) {
        int nsec = (x->tv_usec - y->tv_usec) / 1000000;
        y->tv_usec += 1000000 * nsec;
        y->tv_sec -= nsec;
    }

    // Compute the time remaining to wait.
    //   tv_usec is certainly positive.
    result->tv_sec = x->tv_sec - y->tv_sec;
    result->tv_usec = x->tv_usec - y->tv_usec;
}

// Is the function value at i lesser than that of at j?
bool vertex_compare(float* vertices, int i, int j) {
    return ((vertices[i] < vertices[j]) || (vertices[i] == vertices[j] && i < j));
}

// Is the function value at i greater than that of at j?
bool isGreater(float * vertices, int v1, int v2) {
    return !vertex_compare(vertices, v1, v2);
}

bool isBodyCP(int v, int dimx, int dimy, int dimz, int x, int y, int z,
        int dimxy, int nv, float * vertices) {
    bool mx[8];
    bool mn[8];

    // Make sure that the point does not lie on the boundary
    // Otherwise, you wouldn't be able to +1 or +dimx or +dimxy with confidence!
    if (x == dimx - 1 || y == dimy - 1 || z == dimz - 1) {
        return false;
    }

    int body[8];

    body[0] = v;                // (x,y,z)
    body[1] = v + 1;            // (x+1,y,z)
    body[2] = v + dimx;         // (x,y+1,z)
    body[3] = body[2] + 1;      // (x+1,y+1,z)

    body[4] = v + dimxy;        // (x,y,z+1)
    body[5] = body[4] + 1;      // (x+1,y,z+1)
    body[6] = body[4] + dimx;   // (x,y+1,z+1)
    body[7] = body[6] + 1;      // (x+1,y+1,z+1)

    int noMax = 0;
    int noMin = 0;

    // Iterate over the 8 vertices in the cube
    for (int i = 0; i < 8; i++) {

        // Intialize all variables
        int mxct = 0;
        int mnct = 0;
        mx[i] = false;
        mn[i] = false;

        // Go through the all 3 vertices adjacent to the vertex
        for (int j = 0; j < 3; j++) {
            int adj = bodyAdj[i * 3 + j];
            if (isGreater(vertices, body[i], body[adj])) {
                mxct++;
            }
            else {
                mnct++;
            }
        }

        // If all 3 vertices adjacent to it are lesser in value
        if (mxct == 3) {
            noMax++;
            mx[i] = true;
        }

        // If all 3 vertices adjacent to it are greater in value
        if (mnct == 3) {
            noMin++;
            mn[i] = true;
        }
    }

    // Dear comment reader, look at Figure 7 of  V.  Pascucci  and  K.  Cole-McLaughlin
    // "Parallel  computation  of  the topology of level sets", Algorithmica, vol. 38, no. 1, pp. 249–268, 2003

    // Figure 7(d)
    if (noMax == 4 || noMin == 4) {
        // Ask the real Aditya: According to the diagram, such a configuration has only one body saddle
        // Adhitya: I asked Vijay, such configurations are possible
        // TODO There are two body saddles. for now ignoring.
        return false;
    }

    // Figures 7(a), 7(c) have 0 body saddle(s)
    if (!((noMax == 2 && noMin == 1) || (noMax == 1 && noMin == 2))) {
        // no body saddle.
        return false;
    }

    // Figure 7(b)(ii) has 1 body saddle and Figure 7(b)(i) has 0 body saddle(s)
    if (noMax == 2) {
        if (mx[0] && mx[7]) {
            return true;
        }
        if (mx[1] && mx[6]) {
            return true;
        }
        if (mx[2] && mx[5]) {
            return true;
        }
        if (mx[3] && mx[4]) {
            return true;
        }
    }

    // Reverse of Figure 7(b)(ii)
    if (noMin == 2) {
        if (mn[0] && mn[7]) {
            return true;
        }
        if (mn[1] && mn[6]) {
            return true;
        }
        if (mn[2] && mn[5]) {
            return true;
        }
        if (mn[3] && mn[4]) {
            return true;
        }
    }

    return false;
}

bool isFaceCP(int v, int axis, int dimx, int dimy, int dimz, int x, int y,
        int z, int dimxy, int nv, float * vertices) {
    int face[4];

    // Make sure that there exists a diagnally opposite point
    // with higher coordinates on the same plane
    if (axis == XY && (x == dimx - 1 || y == dimy - 1)) {
        return false;
    }

    if (axis == YZ && (z == dimz - 1 || y == dimy - 1)) {
        return false;
    }

    if (axis == XZ && (x == dimx - 1 || z == dimz - 1)) {
        return false;
    }

    face[0] = v;
    if (axis == XY) {
        face[1] = v + 1;        // (x+1,y,z)
        face[2] = v + dimx;     // (x,y+1,z)
        face[3] = face[2] + 1;  // (x+1,y+1,z)
    }
    if (axis == XZ) {
        face[1] = v + 1;        // (x+1,y,z)
        face[2] = v + dimxy;    // (x,y,z+1)
        face[3] = face[2] + 1;  // (x+1,y,z+1)
    }

    if (axis == YZ) {
        face[1] = v + dimx;         // (x,y+1,z)
        face[2] = v + dimxy;        // (x,y,z+1)
        face[3] = face[2] + dimx;   // (x,y+1,z+1)
    }

    // Origin of arrow indicates which side is greater

    /*

    Both diagonally opposite points on the plane are maxima

    |2|---------->|3|
     |             ^
     |             |
     |             |
     |             |
     v             |
    |0|<--------- |1|

    */
    if (isGreater(vertices, face[0], face[1]) && isGreater(vertices, face[0], face[2]) &&
        isGreater(vertices, face[3], face[1]) && isGreater(vertices, face[3], face[2])) {
        return true;
    }
    /*

    Both diagonally opposite points on the plane are minima

    |2|<----------|3|
     ^             |
     |             |
     |             |
     |             |
     |             v
    |0|---------> |1|

    */

    if (isGreater(vertices, face[1], face[0]) && isGreater(vertices, face[2], face[0]) &&
        isGreater(vertices, face[1], face[3]) && isGreater(vertices, face[2], face[3])) {
        return true;
    }
    return false;
}

void classifyCriticalPoints() {

    gettimeofday(&start1, NULL);

    // create required memory
    type = new NODE_TYPE[num_vert];
    bfType = new int[num_vert];
    steep_asc = new int[num_vert];
    steep_desc = new int[num_vert];

    critical_points = new int[num_vert];
    leaves = new int[num_vert];
    extrema_map = new int[num_vert];

    int nv = num_vert;

    gettimeofday(&end1, NULL);
    timeval_subtract(&result, &end1, &start1);
    mtime = ((result.tv_sec) * 1000 + result.tv_usec / 1000.0) + 0.5;
    printf("Time to classify critical points: %ld milliseconds\n", mtime);
    cpTime = mtime;
    for (int v = 0; v < num_vert; v++) {
        int dimxy = dimx * dimy;
        int z = v / dimxy;
        int xy = v % dimxy;
        int y = xy / dimx;
        int x = xy % dimx;

        // 6 neighbors, v1 to v6
        // v0, v1 along x, v2, v3 along y and v4, v5 along z
        int verts[6] = {-1, -1, -1, -1, -1, -1};
        int bndc[6] = {-1, -1, -1, -1, -1, -1};
        int * bnd = bndc;

        // Used for classifying face extrema as critical
        int b[36] = {1, 2, 3, 4, 5, 0, // x = 0
                     0, 2, 3, 4, 5, 1, // x = dimx - 1
                     3, 0, 1, 4, 5, 2, // y = 0
                     2, 0, 1, 4, 5, 3, // y = dimy - 1
                     5, 0, 1, 2, 3, 4, // z = 0
                     4, 0, 1, 2, 3, 5  // z = dimz - 1
                    };

        if (x != 0) {
            // (x-1,y,z)
            verts[0] = v - 1;
        } else
            bnd = b;
        if (x != dimx - 1) {
            // (x+1,y,z)
            verts[1] = v + 1;
        } else
            bnd = b + 6;
        if (y != 0) {
            // (x,y-1,z)
            verts[2] = v - dimx;
        } else
            bnd = b + 12;
        if (y != dimy - 1) {
            // (x,y+1,z)
            verts[3] = v + dimx;
        } else
            bnd = b + 18;
        if (z != 0) {
            // (x,y,z-1)
            verts[4] = v - dimxy;
        } else
            bnd = b + 24;
        if (z != dimz - 1) {
            // (x,y,z+1)
            verts[5] = v + dimxy;
        } else
            bnd = b + 30;
        int done[6] = {-1, -1, -1, -1, -1, -1};
        // At the end of routine, done has 1, and 2 if it belongs to lower link, and 3 and 4 if it belongs to upper link
        // the numbers denoting the component

        // Can the count be more than 2 in both cases?
        // No. Everytime a vertex has been processed, all 4 of its neighbors are checked.
        // So, 5 vertices are checked in one single shot and and its neighbors are also checked using BFS in the queue
        // And all of them will be assigned the same count (if at all).
        // The remaining one vertex who wouldn't have been checked, will be checked during the iteration
        int lowerLinkCt = 0;
        int upperLinkCt = 0;

        // Represents a queue
        int q[6] = {0, 0, 0, 0, 0, 0};
        int begin = 0;
        int end = 0;

        steep_desc[v] = v;
        steep_asc[v] = v;
        type[v] = REGULAR;

        // count lower link components
        for (int i = 0; i < 6; i++) {
            // is valid vertex
            if (verts[i] != -1) {
                // is it part of lower link, and not considered in the component found so far
                if (done[i] == -1 && isGreater(vertices, v, verts[i])) {
                    // new lower link component
                    lowerLinkCt++;

                    // empty queue
                    begin = 0;
                    end = 0;

                    // add i to queue
                    q[end++] = i;
                    done[i] = lowerLinkCt;
                    // do bfs
                    while (end > begin) {
                        int vv = q[begin++];

                        for (int j = 0; j < 4; j++) {
                            int nbr = neighbors[vv * 4 + j];

                            if (done[nbr] == -1 && verts[nbr] != -1
                                    && isGreater(vertices, v, verts[nbr])) {
                                // add to queue
                                q[end++] = nbr;
                                done[nbr] = lowerLinkCt;
                            }
                        }
                    }
                }
            }
        }
        // count upper link components
        for (int i = 0; i < 6; i++) {
            // is valid vertex
            if (verts[i] != -1) {
                // is it part of upper link, and not considered in the component found so far
                if (done[i] == -1 && isGreater(vertices, verts[i], v)) {
                    // new upper link component
                    upperLinkCt++;

                    // empty queue
                    begin = 0;
                    end = 0;

                    // add i to queue
                    q[end++] = i;
                    done[i] = upperLinkCt + 2;
                    // do bfs
                    while (end > begin) {
                        int vv = q[begin++];

                        for (int j = 0; j < 4; j++) {
                            int nbr = neighbors[vv * 4 + j];

                            if (done[nbr] == -1 && verts[nbr] != -1
                                    && isGreater(vertices, verts[nbr], v)) {
                                // add to queue
                                q[end++] = nbr;
                                done[nbr] = upperLinkCt + 2;
                            }
                        }
                    }
                }
            }
        }

        // do whatever processing is required

        // A vertex is regular if its upper link and one lower link have exactly one component
        if (lowerLinkCt == 1 && upperLinkCt == 1) {
            // Since we have already defined the type of a vertex to be regular above,
            // we dont have to assign/check again in this condition

            // If the point is on the boundary
            if (bnd[0] != -1) {
                // The point outside the boundary is bnd[5]
                // Take the coordinate opposite to it (bnd[0])
                // If that coordinate is independently a link and
                // it's neighbours are present in a different link
                // then the boundary extrema is taken as a saddle
                if ((done[bnd[0]] != done[bnd[1]])
                        && (done[bnd[0]] != done[bnd[2]])
                        && (done[bnd[0]] != done[bnd[3]])
                        && (done[bnd[0]] != done[bnd[4]])) {
                    type[v] = SADDLE;
                }

                // Otherwise, it is just another regular point with exactly one different components
                else {
                    type[v] = REGULAR;
                }
            }
        }

        // A vertex is a minimum if its lower link is empty
        else if (lowerLinkCt == 0 && upperLinkCt == 1) {
            type[v] = MINIMUM;
        }

        // A vertex is a maximum if its upper link is empty
        else if (lowerLinkCt == 1 && upperLinkCt == 0) {
            type[v] = MAXIMUM;
        }

        else if (lowerLinkCt == 0 || upperLinkCt == 0) {
            type[v] = UNDEFINED;
        }

        // A vertex in the input structured grid is a saddle only when its
        // lower (upper) link lies on a plane normal to one of the axes and
        // its upper (lower) link consists of two isolated vertices.
        else {
            // This is clearly satisfied for (lowerLinkCt, upperLinkCt) pairs of (1,2) and (2,1)
            // Q: What about cases like (0,2), (2,0) and (2,2)?
            // A: Consider the case where points 5 and 4 are isolated and are upper links
            // Now, the rest (0,1,2,3) have to be either lower or upper links.
            // Even if one among them are part of an upper link, then points 5 and 4 are not isolated.
            // This implies the rest of the neighbors will be part of the lower link
            // Clearly (0,2) and (2,0) are not possible. Similarly, (2,2) is not possible

            type[v] = SADDLE;
        }

        // store face and body saddles.
        bfType[v] = 0;
        if (isFaceCP(v, XY, dimx, dimy, dimz, x, y, z, dimxy, nv, vertices)) {
            bfType[v] |= XY;
        }
        if (isFaceCP(v, YZ, dimx, dimy, dimz, x, y, z, dimxy, nv, vertices)) {
            bfType[v] |= YZ;
        }
        if (isFaceCP(v, XZ, dimx, dimy, dimz, x, y, z, dimxy, nv, vertices)) {
            bfType[v] |= XZ;
        }
        if (isBodyCP(v, dimx, dimy, dimz, x, y, z, dimxy, nv, vertices)) {
            bfType[v] |= XYZ;
        }

        // store indice if they are ascending/descending
        for (int i = 0; i < 6; i++) {
            if (done[i] == 1) {
                steep_desc[v] = verts[i];
            }
            if (done[i] == 2) {
                steep_desc[v] = verts[i];
            }
            if (done[i] == 3) {
                steep_asc[v] = verts[i];
            }
            if (done[i] == 4) {
                steep_asc[v] = verts[i];
            }
        }
    }
    gettimeofday(&end1, NULL);
    timeval_subtract(&result, &end1, &start1);
    mtime = ((result.tv_sec) * 1000 + result.tv_usec / 1000.0) + 0.5;
    //	printf("Time to classify critical points: %ld milliseconds\n", mtime);
    cpTime = mtime;
}

void updateFaceCP(int vIndex, int v, int axis, int dimx, int dimy, int dimz, int x, int y, int z, int dimxy, int nv, float * vertices, float * function_values, int* lower_roots, int* upper_roots, int* lower_root_count, int* upper_root_count, int* vertex_pos, char* offset) {
    int face[4];

    // Make sure that the point do not lie on the boundary of the plane
    if (axis == FS_XY && (x == dimx - 1 || y == dimy - 1)) {
        return;
    }

    if (axis == FS_YZ && (z == dimz - 1 || y == dimy - 1)) {
        return;
    }

    if (axis == FS_XZ && (x == dimx - 1 || z == dimz - 1)) {
        return;
    }

    face[0] = v;
    if (axis == FS_XY) {
        face[1] = v + 1;            // (x+1,y,z)
        face[2] = v + dimx;         // (x,y+1,z)
        face[3] = face[2] + 1;      // (x+1,y+1,z)
    }
    if (axis == FS_XZ) {
        face[1] = v + 1;            // (x+1,y,z)
        face[2] = v + dimxy;        // (x,y,z+1)
        face[3] = face[2] + 1;      // (x+1,y,z+1)
    }

    if (axis == FS_YZ) {
        face[1] = v + dimx;         // (x,y+1,z)
        face[2] = v + dimxy;        // (x,y,z+1)
        face[3] = face[2] + dimx;   // (x,y+1,z+1)
    }

    /*

    Both diagonally opposite points on the plane are maxima

    |2|---------->|3|
     |             ^
     |             |
     |             |
     |             |
     v             |
    |0|<--------- |1|

    */

    if (isGreater(vertices, face[0], face[1]) && isGreater(vertices, face[0], face[2]) &&
        isGreater(vertices, face[3], face[1]) && isGreater(vertices, face[3], face[2])) {

        lower_root_count[vIndex] = 2;
        upper_root_count[vIndex] = 2;

        // Obviously, 1 and 2 are lower roots!
        lower_roots[vIndex * MAX_ADJ + 0] = face[1];
        lower_roots[vIndex * MAX_ADJ + 1] = face[2];

        // Obviously, 0 and 3 are upper roots!
        upper_roots[vIndex * MAX_ADJ + 0] = face[0];
        upper_roots[vIndex * MAX_ADJ + 1] = face[3];

        // Update the vertex position to the greater of the 2 maximas!
        if (isGreater(vertices, face[1], face[2])) {
            vertex_pos[vIndex] = face[1];
        }
        else {
            vertex_pos[vIndex] = face[2];
        }
        offset[vIndex] = 1;
    }

    /*

    Both diagonally opposite points on the plane are minima

    |2|<----------|3|
     ^             |
     |             |
     |             |
     |             |
     |             v
    |0|---------> |1|

    */

    if (isGreater(vertices, face[1], face[0]) && isGreater(vertices, face[2], face[0]) &&
        isGreater(vertices, face[1], face[3]) && isGreater(vertices, face[2], face[3])) {

        lower_root_count[vIndex] = 2;
        upper_root_count[vIndex] = 2;

        // Obviously, 0 and 3 are lower roots!
        lower_roots[vIndex * MAX_ADJ + 0] = face[0];
        lower_roots[vIndex * MAX_ADJ + 1] = face[3];

        // Obviously, 1 and 2 are upper roots!
        upper_roots[vIndex * MAX_ADJ + 0] = face[1];
        upper_roots[vIndex * MAX_ADJ + 1] = face[2];

        // Update the vertex position to the greater of the 2 minimas!
        if (isGreater(vertices, face[0], face[3])) {
            vertex_pos[vIndex] = face[0];
        }
        else {
            vertex_pos[vIndex] = face[3];
        }
        offset[vIndex] = 1;
    }

    {

        // Okay, now you know the function values at all the four corners!
        // Now it should be essential to find the coordinates of the saddle using the function values
        // and then use the coordinates to find out the function value at the saddle!

        /*

                    2    B (x,1)      3

                    +----+------------+
                    |    |            |
                    |    | S (x,y)    |
                    |    |            |
            C (0,y) +-----------------+ D (1,y)
                    |    |            |
                    |    |            |
                    |    |            |
                    +----+------------+

                    0    A (x,0)      1

        f(A) = f(0) (1-x) + f(1) (x) = f(0) + x(f(1) - f(0))
        f(B) = f(2) (1-x) + f(3) (x) = f(2) + x(f(3) - f(2))

        f(S) = f(A) (1-y) + f(B) (y) = f(A) + y(f(B) - f(A))
             = axy + bx + cy + d

        Then, solving for the coefficients

        a = f(0) + f(3) - f(1) - f(2)
        b = f(1) - f(0)
        c = f(2) - f(0)
        d = f(0)

        By taking the partial derivatives of axy + bx + cy + d, we get (ay+b), (ax+c)
        Now, by equating those to 0, we get the coordinates of the saddle point which are,
        x = -c/a and y = -b/a

        */
        float corners[4];

        // Store the function values in corners
        for (int i = 0; i < 4; i++) {
            corners[i] = vertices[face[i]];
        }

        float a = corners[0] + corners[3] - corners[1] - corners[2];

        if (a == 0) {
            // Just use the value of d here!
            function_values[vIndex] = corners[0];
            return;
        }

        float c = corners[2] - corners[0];
        float x = -c / a;

        if (x > 1 || x < 0) {
            function_values[vIndex] = corners[0];
            return;
        }

        float b = corners[1] - corners[0];
        float y = -b / a;

        if (y > 1 || y < 0) {
            function_values[vIndex] = corners[0];
            return;
        }

        // Use the equation of the bilinear equation!
        function_values[vIndex] = a * x * y + b * x + c * y + corners[0];
    }
}

// Get index for the saddle in range lowMax to highMin
int getIndx(float* vertices, float current, int lowMax, int highMin) {
    // The face saddle will have the boundaries of lowMax and highMin
    float fmax = vertices[lowMax];
    if (current == fmax) {
        return lowMax - 1;
    }
    return highMin;
}

void updateBodyCP(int vIndex, int v, int dimx, int dimy, int dimz, int x, int y,
        int z, int dimxy, int nv, float * vertices, float * function_values,
        int* lower_roots, int* upper_roots, int* lower_root_count,
        int* upper_root_count, int* vertex_pos, char* offset) {

    bool mx[8];
    bool mn[8];

    // Store in a temporary variable
    int vIndex1 = vIndex;

    // Make sure that the point does not lie on the boundary
    // Otherwise, you wouldn't be able to +1 or +dimx or +dimxy with confidence!
    if (x == dimx - 1 || y == dimy - 1 || z == dimz - 1) {
        return;
    }

    int body[8];

    body[0] = v;                // (x,y,z)
    body[1] = v + 1;            // (x+1,y,z)
    body[2] = v + dimx;         // (x,y+1,z)
    body[3] = body[2] + 1;      // (x+1,y+1,z)

    body[4] = v + dimxy;        // (x,y,z+1)
    body[5] = body[4] + 1;      // (x+1,y,z+1)
    body[6] = body[4] + dimx;   // (x,y+1,z+1)
    body[7] = body[6] + 1;      // (x+1,y+1,z+1)

    int noMax = 0;
    int noMin = 0;

    // Highest amongst all minimas
    int highMin;

    // Lowest amongst all maximas
    int lowMax;

    // Iterate over the 8 vertices in the cube
    for (int i = 0; i < 8; i++) {
        int mxct = 0;
        int mnct = 0;
        mx[i] = false;
        mn[i] = false;

        // Go through the all 3 vertices adjacent to the vertex
        for (int j = 0; j < 3; j++) {
            int adj = bodyAdj[i * 3 + j];
            if (isGreater(vertices, body[i], body[adj])) {
                mxct++;
            }
            else {
                mnct++;
            }
        }

        // If all 3 vertices adjacent to it are lesser in value
        if (mxct == 3) {
            noMax++;
            mx[i] = true;
            // Adhitya: Is this correct? What happens when noMax is 2? -- Yes! It is OR!
            if (noMax == 1 || isGreater(vertices, lowMax, body[i])) {
                // Reduce the value of lowMax to current Maxima
                lowMax = body[i];
            }
        }

        // If all 3 vertices adjacent to it are greater in value
        if (mnct == 3) {
            noMin++;
            mn[i] = true;
            if (noMin == 1 || vertex_compare(vertices, highMin, body[i])) {
                // Increase the value of highMin to current Minima
                highMin = body[i];
            }
        }
    }

    // Dear comment reader, look at Figure 7 of  V.  Pascucci  and  K.  Cole-McLaughlin
    // "Parallel  computation  of  the topology of level sets", Algorithmica, vol. 38, no. 1, pp. 249–268, 2003

    // Figure 7(d)
    if (noMax == 4 || noMin == 4) {
        // TODO There are two body saddles. for now ignoring.
        return;
    }

    if (!((noMax == 2 && noMin == 1) || (noMax == 1 && noMin == 2))) {
        // no body saddle.
        return;
    }

    // Figure 7(b)(ii) has 1 body saddle and Figure 7(b)(i) has 0 body saddle(s)
    // In all these cases there will be two distinct upper roots!
    if (noMax == 2) {
        if (mx[0] && mx[7]) {
            upper_roots[vIndex1 * MAX_ADJ + 0] = body[0];
            upper_roots[vIndex1 * MAX_ADJ + 1] = body[7];
        }
        if (mx[1] && mx[6]) {
            upper_roots[vIndex1 * MAX_ADJ + 0] = body[1];
            upper_roots[vIndex1 * MAX_ADJ + 1] = body[6];
        }
        if (mx[2] && mx[5]) {
            upper_roots[vIndex1 * MAX_ADJ + 0] = body[2];
            upper_roots[vIndex1 * MAX_ADJ + 1] = body[5];
        }
        if (mx[3] && mx[4]) {
            upper_roots[vIndex1 * MAX_ADJ + 0] = body[3];
            upper_roots[vIndex1 * MAX_ADJ + 1] = body[4];
        }

        int min = 0;

        // Find the lowest amongst all roots
        for (int i = 1; i < 8; i++) {
            if (isGreater(vertices, body[min], body[i])) {
                min = i;
            }
        }
        lower_root_count[vIndex1] = 1;
        upper_root_count[vIndex1] = 2;
        lower_roots[vIndex1 * MAX_ADJ] = body[min];
        offset[vIndex] = 1;
    }

    // Reverse of Figure 7(b)(ii)
    // In all these cases there will be two distinct lower roots!
    if (noMin == 2) {
        if (mn[0] && mn[7]) {
            lower_roots[vIndex1 * MAX_ADJ + 0] = body[0];
            lower_roots[vIndex1 * MAX_ADJ + 1] = body[7];
        }
        if (mn[1] && mn[6]) {
            lower_roots[vIndex1 * MAX_ADJ + 0] = body[1];
            lower_roots[vIndex1 * MAX_ADJ + 1] = body[6];
        }
        if (mn[2] && mn[5]) {
            lower_roots[vIndex1 * MAX_ADJ + 0] = body[2];
            lower_roots[vIndex1 * MAX_ADJ + 1] = body[5];
        }
        if (mn[3] && mn[4]) {
            lower_roots[vIndex1 * MAX_ADJ + 0] = body[3];
            lower_roots[vIndex1 * MAX_ADJ + 1] = body[4];
        }

        // Find the greatest amongst all roots
        int max = 0;
        for (int i = 1; i < 8; i++) {
            if (isGreater(vertices, body[i], body[max])) {
                max = i;
            }
        }
        lower_root_count[vIndex1] = 2;
        upper_root_count[vIndex1] = 1;
        upper_roots[vIndex1 * MAX_ADJ] = body[max];
        offset[vIndex] = 1;
    }

    {

    /*

                      6 (0,1,1)      F         7 (1,1,1)
                       ________________________
                     /|             /         /|
                    / |            /         / |
                   /  |           /         /  |
                  /   |          /         /   |
               G /--------------/---------/    |
                /     |        / Q       /H    |
               /      |       /         /      |
              /       |      /         /       |
    4 (0,0,1)/________|_____/_________/ 5      |
             |        |     E         | (1,0,1)|
             |        |               |        |
             |        |__________B____|________| 3 (1,1,0)
             |       /             /  |       /
             |      / 2 (0,1,0)   /   |      /
             |     /             /    |     /
             |   C/-------------/-----|----/D
             |   /             / P    |   /
             |  /             /       |  /
             | /             /        | /
             |/_____________/_________|/
                         A
             0 (0,0,0)                1 (1,0,0)


        f(A) = f(0) (1-x) + f(1) (x)
        f(B) = f(2) (1-x) + f(3) (x)
        f(E) = f(4) (1-x) + f(5) (x)
        f(F) = f(6) (1-x) + f(7) (x)

        f(P) = f(A) (1-y) + f(B) (y)
        f(Q) = f(E) (1-y) + f(F) (y)

        f(S) = f(P) (1-z) + f(Q) (z)

        The Trilinear Interpolant equation is,

        f(x,y,z) = axyz + bxy + cxz + dyz + ex + fy + gz + h

        You may want to solve for the coefficients by yourself.

        For a better understanding, read the subtopic titled "Trilinear Interpolant on a Parallelopiped"
        V.  Pascucci  and  K.  Cole-McLaughlin, "Parallel  computation  of  the topology of level sets"

        The only small difference is that, our coder decided to divide r by (bc - ae) instead of multiplying it
        And then he used that to find x. And then he applied that formula of x to find y and z.

        Do not worry, the formulas are correct. Adhitya spent the entire day of 8-12-2016 in finding this fact.
    */

        float corners[8];

        for (int i = 0; i < 8; i++) {
            corners[i] = vertices[body[i]];
        }

        float a = -corners[F000] + corners[F001] + corners[F010] - corners[F011]
                + corners[F100] - corners[F101] - corners[F110] + corners[F111];
        float b = corners[F000] - corners[F010] - corners[F100] + corners[F110];
        float c = corners[F000] - corners[F001] - corners[F100] + corners[F101];
        float d = corners[F000] - corners[F001] - corners[F010] + corners[F011];
        float e = -corners[F000] + corners[F100];
        float g = -corners[F000] + corners[F010];
        float h = -corners[F000] + corners[F001];
        float k = corners[F000];

        if (noMax == 2) {
            function_values[vIndex] = vertices[highMin];
        }
        else {
            function_values[vIndex] = vertices[lowMax];
        }

        // Just check that the denominator does not hit the ceiling!
        if (a == 0) {
            if (b == 0 || c == 0 || d == 0) {
                vertex_pos[vIndex] = getIndx(vertices, function_values[vIndex1], lowMax, highMin);
                return;
            }
            // Adhitya: I am just guessing that the limits were taken and they ended up with this formula!
            // Adhitya: My guess could be wrong. Ask Vijay!
            float x = (d * e - c * g - b * h) / (2 * b * c);
            if (x < 0 || x > 1) {
                vertex_pos[vIndex] = getIndx(vertices, function_values[vIndex1], lowMax, highMin);
                return;
            }
            float y = (-d * e + c * g - b * h) / (2 * b * d);
            if (y < 0 || y > 1) {
                vertex_pos[vIndex] = getIndx(vertices, function_values[vIndex1], lowMax, highMin);
                return;
            }
            float z = (-d * e - c * g + b * h) / (2 * c * d);
            if (z < 0 || z > 1) {
                vertex_pos[vIndex] = getIndx(vertices, function_values[vIndex1], lowMax, highMin);
                return;
            }
            function_values[vIndex1] = a * x * y * z + b * x * y  + c * x * z + d * y * z + e * x + g * y + h * z + k;
            vertex_pos[vIndex] = getIndx(vertices, function_values[vIndex1], lowMax, highMin);
        }

        else {
            // Again, just check the denominator of x
            if ((a * e - b * c) == 0) {
                vertex_pos[vIndex] = getIndx(vertices, function_values[vIndex1], lowMax, highMin);
                return;
            }
            bool ret1 = true;
            bool ret2 = true;
            float r = (b * d - a * g) * (c * d - a * h) / (b * c - a * e);
            if (r < 0) {
                vertex_pos[vIndex] = getIndx(vertices, function_values[vIndex1], lowMax, highMin);
                return;
            }

            r = (float) sqrt(r);
            float x = (-d + r) / a;
            if (x < 0 || x > 1) {
                ret1 = false;
            }
            float y = -(c * x + h) / (a * x + d);
            if (y > 1 || y < 0 || (a * x + d) == 0) {
                ret1 = false;
            }
            float z = -(b * x + g) / (a * x + d);
            if (z > 1 || z < 0 || (a * x + d) == 0) {
                ret1 = false;
            }
            if (ret1 != false) {
                function_values[vIndex1] = a * x * y * z + b * x * y + c * x * z + d * y * z + e * x + g * y + h * z + k;
                vertex_pos[vIndex] = getIndx(vertices, function_values[vIndex1], lowMax, highMin);
                return;
            }

            // This comes from the paper itself - This can be derived
            x = (-d - r) / a;
            if (x < 0 || x > 1) {
                ret2 = false;
            }

            // Use the value of x for calculating y
            // ax + d is equivalent to sqrt(r)
            y = -(c * x + h) / (a * x + d);
            if (y > 1 || y < 0 || (a * x + d) == 0) {
                ret2 = false;
            }

            // Use the value of x for calculating z
            z = -(b * x + g) / (a * x + d);
            if (z > 1 || z < 0 || (a * x + d) == 0) {
                ret2 = false;
            }
            if (ret2 != false) {
                function_values[vIndex1] = a * x * y * z + b * x * y + c * x * z + d * y * z + e * x + g * y + h * z + k;
                vertex_pos[vIndex] = getIndx(vertices, function_values[vIndex1], lowMax, highMin);
                return;
            }
            vertex_pos[vIndex] = getIndx(vertices, function_values[vIndex1], lowMax, highMin);
            // TODO two saddles, then both ret1 and ret2 will be true, so use appropriate indices.
        }
    }

}

void fillCriticalPoints() {
    int cp = 0;

    for (cp = 0; cp < num_critical; cp++) {
        int v = critical_points[cp];

        lower_root_count[cp] = 0;
        upper_root_count[cp] = 0;

        // If the critical point is a maximum/minimum/boundary saddle
        if (v < num_vert) {
            vertex_pos[cp] = v;
            offset[cp] = 0;
            function_values[cp] = vertices[v];

            int dimxy = dimx * dimy;
            int z = v / dimxy;
            int xy = v % dimxy;
            int y = xy / dimx;
            int x = xy % dimx;

            // 6 neighbors, v1 to v6
            // v0, v1 along x, v2, v3 along y and v4, v5 along z
            int verts[6] = {-1, -1, -1, -1, -1, -1};

            if (x != 0) {
                // (x-1,y,z)
                verts[0] = v - 1;
            }
            if (x != dimx - 1) {
                // (x+1,y,z)
                verts[1] = v + 1;
            }
            if (y != 0) {
                // (x,y-1,z)
                verts[2] = v - dimx;
            }
            if (y != dimy - 1) {
                // (x,y+1,z)
                verts[3] = v + dimx;
            }
            if (z != 0) {
                // (x,y,z-1)
                verts[4] = v - dimxy;
            }
            if (z != dimz - 1) {
                // (x,y,z+1)
                verts[5] = v + dimxy;
            }

            int done[6] = {-1, -1, -1, -1, -1, -1};
            // At the end of routine, done has 1, and 2 if it belongs to lower link, and 3 and 4 if it belongs to upper link
            // the numbers denoting the component
            int lowerLinkCt = 0;
            int upperLinkCt = 0;

            int q[6] = {0, 0, 0, 0, 0, 0};
            int begin = 0;
            int end = 0;

            // count lower link components
            for (int i = 0; i < 6; i++) {
                // is valid vertex
                if (verts[i] != -1) {
                    // is it part of lower link, and not considered in the component found so far
                    if (done[i] == -1 && isGreater(vertices, v, verts[i])) {
                        // new lower link component
                        lowerLinkCt++;

                        // empty queue
                        begin = 0;
                        end = 0;

                        // add i to queue
                        q[end++] = i;
                        done[i] = lowerLinkCt;
                        // do bfs
                        while (end > begin) {
                            int vv = q[begin++];

                            for (int j = 0; j < 4; j++) {
                                int nbr = neighbors[vv * 4 + j];

                                if (done[nbr] == -1 && verts[nbr] != -1 && isGreater(vertices, v, verts[nbr])) {
                                    // add to queue
                                    q[end++] = nbr;
                                    done[nbr] = lowerLinkCt;
                                }
                            }
                        }
                    }
                }
            }

            // count upper link components
            for (int i = 0; i < 6; i++) {
                // is valid vertex
                if (verts[i] != -1) {
                    // is it part of upper link, and not considered in the component found so far
                    if (done[i] == -1 && isGreater(vertices, verts[i], v)) {
                        // new upper link component
                        upperLinkCt++;

                        // empty queue
                        begin = 0;
                        end = 0;

                        // add i to queue
                        q[end++] = i;
                        done[i] = upperLinkCt + 2;
                        // do bfs
                        while (end > begin) {
                            int vv = q[begin++];

                            for (int j = 0; j < 4; j++) {
                                int nbr = neighbors[vv * 4 + j];

                                if (done[nbr] == -1 && verts[nbr] != -1 && isGreater(vertices, verts[nbr], v)) {
                                    // add to queue
                                    q[end++] = nbr;
                                    done[nbr] = upperLinkCt + 2;
                                }
                            }
                        }
                    }
                }
            }

            lower_root_count[cp] = lowerLinkCt;
            upper_root_count[cp] = upperLinkCt;

            // Get the vertex from each component so that one can traverse
            for (int i = 0; i < 6; i++) {
                // Store the lower link component!
                if (done[i] == 1) {
                    lower_roots[cp * MAX_ADJ + 0] = verts[i];
                }
                // If there is a second lower link component, then it should be stored!
                if (done[i] == 2) {
                    lower_roots[cp * MAX_ADJ + 1] = verts[i];
                }
                // Store the upper link component!
                if (done[i] == 3) {
                    upper_roots[cp * MAX_ADJ + 0] = verts[i];
                }
                // If there is a second upper link component, then it should be stored!
                if (done[i] == 4) {
                    upper_roots[cp * MAX_ADJ + 1] = verts[i];
                }
            }
        }

        // its a body or face saddle
        else {
            for (int i = 1; i <= 3; i++) {
                // Adhitya: I feel that v is calculated wrongly. It should be v = v - (i * num_vert)
                // Answer: No. In initializeDataStructures(), whenever there is a saddle we increment i by
                // n * face and we keep overwriting it in critical_points[].
                // Over here, we check it each time over each face, so all cases will be satisfied
                v -= num_vert;
                if ((v < num_vert)&&(v >= 0)) {
                    // face saddle with axis i;
                    int dimxy = dimx * dimy;
                    int z = v / dimxy;
                    int xy = v % dimxy;
                    int y = xy / dimx;
                    int x = xy % dimx;

                    updateFaceCP(cp, v, i, dimx, dimy, dimz, x, y, z, dimxy, num_vert, vertices, function_values,
                            lower_roots, upper_roots, lower_root_count, upper_root_count, vertex_pos, offset);
                    //continue;
                }
            }

            // body saddle
            v -= num_vert;

            if ((v < num_vert)&&(v >= 0)){
                int dimxy = dimx * dimy;
                int z = v / dimxy;
                int xy = v % dimxy;
                int y = xy / dimx;
                int x = xy % dimx;

                updateBodyCP(cp, v, dimx, dimy, dimz, x, y, z, dimxy, num_vert, vertices, function_values,
                    lower_roots, upper_roots, lower_root_count, upper_root_count, vertex_pos, offset);
            }
        }
    }

    delete[] vertices;
}

int find_min_root(int i, int *steep_desc) {
    int k = i;

    // steep_desc was initialised with the index values in classifyCriticalPoints()
    // If the array contained the index value itself, then it was not processed at all
    // This is loop condition is just a fancy way of writing an infinite loop :P
    while (k != steep_desc[k]) {

        k = steep_desc[k];

        //if (critical_points[k]>num_vert) break;
        // The basic idea over here is that, keep on descending downwards till you hit a saddle
        if (type[k] == SADDLE) break;
    }
    return k;
}

int find_max_root(int i, int *steep_asc) {
    int k = i;

    // steep_asc was initialised with the index values in classifyCriticalPoints()
    // If the array contained the index value itself, then it was not processed at all
    // (I will not repeat my joke twice)
    while (k != steep_asc[k]) {

        k = steep_asc[k];
        //if (critical_points[k]>num_vert) break;
        // The basic idea over here is that, keep on ascending upwards till you hit a saddle
        if (type[k] == SADDLE) break;
    }
    return k;
}

void find_roots() {
    gettimeofday(&start1, NULL);
    for (int i = 0; i < num_critical; i++) {
        int c = i*MAX_ADJ;
        int c1 = c;
        // Use one vertex of the lower roots - Remember if there are two components,
        // you have *a* vertex from each of them
        for (int j = 0; j < lower_root_count[i]; j++) {
            //lower_roots[i*MAX_ADJ+j] = extrema_map[find_min_root(lower_roots[i*MAX_ADJ+j], steep_desc)];
            lower_roots[c] = find_min_root(lower_roots[c], steep_desc);
            // find_min_root returns the minima
            // We stored the critical point index of such extremas in extrema_map[] : initializeDataStructures()
            lower_roots[c] = extrema_map[lower_roots[c]];
            c++;
        }
        c = c1;
        for (int j = 0; j < upper_root_count[i]; j++) {
            // upper_roots[i*MAX_ADJ+j] = extrema_map[find_max_root(upper_roots[i*MAX_ADJ+j], steep_asc)];
            // find_max_root returns the maxima
            upper_roots[c] = find_max_root(upper_roots[c], steep_asc);
            upper_roots[c] = extrema_map[upper_roots[c]];
            c++;
        }
    }
    gettimeofday(&end1, NULL);
    timeval_subtract(&result, &end1, &start1);
    delete[] extrema_map;
    delete[] steep_asc;
    delete[] steep_desc;
    mtime = ((result.tv_sec) * 1000 + result.tv_usec / 1000.0) + 0.5;
    printf("Time to find roots: %ld milliseconds\n", mtime);
    rootTime = mtime;
}

void initarray() {
    function_values = new float[num_critical];
    vertex_index = new int[num_critical];

    join_neigh = new int[num_critical];
    join_children = new int[2 * num_critical];

    split_neigh = new int[num_critical];
    split_children = new int[2 * num_critical];

    vertex_pos = new int[num_critical];
    offset = new char[num_critical];

    lower_root_count = new int[num_critical];
    upper_root_count = new int[num_critical];

    lower_roots = new int[num_critical * MAX_ADJ];
    upper_roots = new int[num_critical * MAX_ADJ];

    for (int i = 0; i < num_critical; i++) {
        join_neigh[i] = -1;
        split_neigh[i] = -1;

        join_children[2 * i] = -1;
        join_children[2 * i + 1] = -1;

        split_children[2 * i] = -1;
        split_children[2 * i + 1] = -1;
        //is_processed[i] = 0;
        vertex_index[i] = i;
    }
}

void sortVertices() {
    gettimeofday(&start1, NULL);
    svi = new int[num_critical];
    //	sort(vertex_index, vertex_index + num_critical, compare());

    // vertex_index[i] = i - As simple as that for all num_critical
    // num_critical = num_minima + num_maxima + num_saddle

    mergesort(vertex_index, num_critical, svi, num_cores);
    //QuickSortOmp(vertex_index, num_critical,15);
    delete[] vertex_index;
    vertex_index = svi;

    gettimeofday(&end1, NULL);

    timeval_subtract(&result, &end1, &start1);
    mtime = ((result.tv_sec) * 1000 + result.tv_usec / 1000.0) + 0.5;
    printf("Time to sort vertices : %ld milliseconds\n", mtime);
    sortTime = mtime;
    //isSorted(vertex_index, num_critical);
    //printf("\n...%s....\n",isSorted(vertex_index, num_critical)?"true":"false");
}

void initializeDataStructures() {

    gettimeofday(&start1, NULL);
    int cpNo = 0;
    int leavesCt = 0;
    int xyIndex = num_vert * FS_XY;
    int yzIndex = num_vert * FS_YZ;
    int xzIndex = num_vert * FS_XZ;
    int bsindxex = num_vert * BS_1;

    for (int i = 0; i < num_vert; i++) {
        extrema_map[i] = 0;
        if (type[i] == MINIMUM) {
            num_minima++;
            leaves[leavesCt++] = cpNo;
            critical_points[cpNo] = i;
            extrema_map[i] = cpNo;
            cpNo++;
        }
        else if (type[i] == MAXIMUM) {
            num_maxima++;
            leaves[leavesCt++] = cpNo;
            critical_points[cpNo] = i;
            extrema_map[i] = cpNo;
            cpNo++;
        }
        // This takes Boundary saddles into consideration
        else if (type[i] == SADDLE) {
            num_saddle++;
            critical_points[cpNo] = i;
            extrema_map[i] = cpNo;
            cpNo++;
        }
        // There is a face saddle on the XY plane
        if ((bfType[i] & XY) != 0) {
            num_saddle++;
            // Adhitya: I have absolutely no clue why this array is stored in such a way
            // Answer: Just so that v > num_vert, its easier to detect saddles in fillCriticalPoints()
            // xyIndex = num_vert * FS_XY
            critical_points[cpNo] = (i + xyIndex);
            cpNo++;
        }
        // There is a face saddle on the YZ plane
        if ((bfType[i] & YZ) != 0) {
            num_saddle++;
            critical_points[cpNo] = (i + yzIndex);
            cpNo++;
        }
        // There is a face saddle on the XZ place
        if ((bfType[i] & XZ) != 0) {
            num_saddle++;
            critical_points[cpNo] = (i + xzIndex);
            cpNo++;
        }
        if ((bfType[i] & XYZ) != 0) {
            num_saddle++;
            critical_points[cpNo] = (i + bsindxex);
            cpNo++;
        }
    }

    num_critical = num_minima + num_maxima + num_saddle;
    num_regular = num_vert - num_critical;

    gettimeofday(&end1, NULL);
    timeval_subtract(&result, &end1, &start1);
    mtime = ((result.tv_sec) * 1000 + result.tv_usec / 1000.0) + 0.5;
    initSetupTime = mtime;

    gettimeofday(&start1, NULL);

    initarray();

    fillCriticalPoints();

    gettimeofday(&end1, NULL);
    timeval_subtract(&result, &end1, &start1);
    mtime = ((result.tv_sec) * 1000 + result.tv_usec / 1000.0) + 0.5;
    printf("Time to initialize data structures: %ld milliseconds\n", mtime);
    initParTime = mtime;
}

void findCriticalPointsGrid() {

    classifyCriticalPoints();

    initializeDataStructures();
}

void cleanup_trees() {
    gettimeofday(&start1, NULL);

    int critical = 0;
    int* x = new int[num_critical];
    int* is_critical = new int[num_critical];

    for (int i = 0; i < num_critical; i++) {
        int j = i; //check
        is_critical[j] = 0; // Initialize array
        x[j] = -1;
        long long int vpos = vertex_pos[j]; // Position in vertices

        // Familiar chant: we used this everywhere!
        long long int dimxy = dimx * dimy;
        long long int posz = vpos / dimxy;
        long long int xy = vpos % dimxy;
        long long int posy = xy / dimx;
        long long int posx = xy % dimx;

            // Check if critical point is a degree-2 node in Join Tree
        if ((join_children[2 * j] == -1)  || (join_children[2 * j + 1] != -1) ||
            // Check if critical point is a degree-2 node in Split Tree
            (split_children[2 * j] == -1) || (split_children[2 * j + 1] != -1) ||
            // Lies on the x boundary
            (posx == 0) || (posx == dimx - 1) ||
            // Lies on the y boundary
            (posy == 0) || (posy == dimy - 1) ||
            // Lies on the z boundary
            (posz == dimz - 1) || (posz == 0)) {
                // This point is the most critical that you can get :P
                x[j] = critical;
                is_critical[j] = 1;
                critical++;
        }

    }

    // Total number of critical points that are degree-2 in
    // both Join and Split Trees and do not lie on the boundary
    final = critical;
    printf("\n------ new critical: %d---------\n", final);

    /*
    join_neigh = j_n; split_neigh = s_n; join_children = j_c;
    split_children = s_c; vertex_index = v_i; function_values = f_v;
    vertex_pos = v_p
    */
    j_n = new int[critical];
    s_n = new int[critical];
    j_c = new int[2 * critical];
    s_c = new int[2 * critical];
    v_i = new int[critical];
    v_p = new long long int[critical];
    f_v = new float[critical];
    offset_c = new char[critical];

    for (int i = 0; i < 2 * critical; i++) {
        j_c[i] = s_c[i] = -1;
    }

    int c = 0;
    int leaf_count = 0;

    long long int scxy = divx * divy;
    long long int scz = subcube / scxy;
    long long int sxy = subcube % scxy;
    long long int scy = sxy / divx;
    long long int scx = sxy % divx;

    /*
        i --> index
        x[i] --> indice of critical point after pruning
        vertex_index[i] --> indice after being sorted
        x[vertex_index[i]] --> In these pruned set of critical points, What is it's indice in the sorted order?
        join_neigh[i] --> Join neighbour for the original indice
        join_neigh[vertex_index[i] --> In the pruned set of critical points, What is the sorted indice of its join neighbour?
    */
    for (int i = 0; i < num_critical; i++) {
        int j = vertex_index[i];
        int t = join_neigh[j];

        if ((join_children[2 * j] == -1) || (split_children[2 * j] == -1)){
            leaves[leaf_count++] = x[j];
        }

        // Check if point after pruning is critical
        if (is_critical[j] == 1) { // x[j] will exist in this case
            v_i[c] = x[j]; //check
            f_v[x[j]] = function_values[j];

            long long int vpos = vertex_pos[j];
            long long int dimxy = dimx * dimy;
            long long int posz = vpos / dimxy;
            long long int xy = vpos % dimxy;
            long long int posy = xy / dimx;
            long long int posx = xy % dimx;

            long long int final_pos = (scz * sub_sizez + posz) * global_dimx * global_dimy +
                                      (scy * sub_sizey + posy) * global_dimx +
                                      (scx * sub_sizex + posx);

            v_p[x[j]] = final_pos;
            offset_c[x[j]] = offset[j];

            // If Join neighbour is not critical, then join with parent
            while ((t != -1) && (is_critical[t] != 1)) {
                t = join_neigh[t];
            }

            // Join neighbour does not exist
            if (t == -1) {
                j_n[x[j]] = -1;
            }

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

            // Time for Pruning the Split Tree!
            t = split_neigh[j];

            // If Split neighbour is not critical, then join with parent
            while ((t != -1) && (is_critical[t] != 1)) {
                t = split_neigh[t];
            }

            // Split neighbour does not exist
            if (t == -1) {
                s_n[x[j]] = -1;
            }

            // The newly found parent is critical, and therefore is the new split neighbour
            else {
                s_n[x[j]] = x[t];
                // Mark this node as the child of the new parent
                // If the first child has not been assigned, then take it's place :)
                if (s_c[2 * x[t]] == -1){
                    s_c[2 * x[t]] = x[j];
                }
                // Otherwise, the second one is up for grabs
                else {
                    s_c[2 * x[t] + 1] = x[j];
                }
            }
            c++;
        }
    }

    // The temporary arrays can finally be removed and assigned to the much more fancier arrays
    int temp = num_critical;
    delete[] join_neigh;
    delete[] split_neigh;
    delete[] join_children;
    delete[] split_children;
    delete[] vertex_index;
    delete[] function_values;
    join_neigh = j_n;
    split_neigh = s_n;
    join_children = j_c;
    split_children = s_c;

    num_critical = critical;
    vertex_index = v_i;
    function_values = f_v;
    offset = offset_c;
    gettimeofday(&end1, NULL);
    timeval_subtract(&result, &end1, &start1);
    mtime = ((result.tv_sec) * 1000 + result.tv_usec / 1000.0) + 0.5;
    cleanuptime = mtime;
    printf(
            "Time to cleanup contour tree: %ld milliseconds, leaves: %d. critical_prev: %d\n",
            mtime, leaf_count, temp);
}

// Read the grid for input file
void readInputGrid(char * file) {
    ifstream input_file(file, ios::in | ios::binary);
    unsigned char * data = new unsigned char[num_vert];
    input_file.read((char*) data, sizeof (unsigned char) * num_vert);
    input_file.close();

    vertices = new float[num_vert * 6];

    for (int i = 0; i < num_vert; i++) {
        vertices[i] = (float) data[i];
    }
}

int max_lower_root_count = 0;

void serialtree() {
    //testRoots();
    sortVertices();

    gettimeofday(&start1, NULL);

    // Compute the Join Union
    struct Elems *jU = new struct Elems[num_critical];
    //int * Rep = new int[num_critical];

    for (int s = 0; s < num_critical; s++) {
        // vertex_index now has the list of sorted indexes
        int i = vertex_index[s];

        // Just for the parent node
        ElemSet(jU, i);

        // lower_root_count ranges between 0 and 2
        for (int j = 0; j < lower_root_count[i]; j++) {
            int k = i * MAX_ADJ + j;
            int rt = lower_roots[k];

            // Returns the newly found parent
            int f = findit(jU, rt);

            // What the node actually represents
            int r = jU[f].rep;

            // The newly found parent represents the current node, it's already been processed
            if (r == i) {
                // Break from inner loop
                continue;
            }

            // So the parent we found now, will be joined with our current parent
            // This means both of them will be neighbours
            join_neigh[r] = i;

            // For the first component
            if (join_children[2 * i] == -1) {
                join_children[2 * i] = r;
            }

            // For the second component
            else {
                join_children[2 * i + 1] = r;
            }

            // Better: unify(jU, newly_found_parent, current_parent)
            int g = unify(jU, f, i);
            jU[g].rep = i;
        }
    }
    delete[] jU;

    // Compute the Split Union
    struct Elems *sU = new struct Elems[num_critical];
    //int * Rep = new int[num_critical];

    for (int s = 0; s < num_critical; s++) {

        int i = vertex_index[num_critical - s - 1];

        // Just for the parent node
        ElemSet(sU, i);

        // Traverse the Upper Roots this time
        for (int j = 0; j < upper_root_count[i]; j++) {

            int k = i * MAX_ADJ + j;
            int rt = upper_roots[k];

            // Returns the newly found parent
            int f = findit(sU, rt);

            // What the node actually represents
            int r = sU[f].rep;

            // The newly found parent represents the current node, it's already been processed
            if (r == i){
                continue;
            }

            // So the parent we found now, and this will split with our current parent
            // This means both of them will be neighbours
            split_neigh[r] = i;

            // For the first component
            if (split_children[2 * i] == -1){
                split_children[2 * i] = r;
            }

            // For the second component
            else {
                split_children[2 * i + 1] = r;
            }

            // Better: unify(sU, newly_found_parent, current_parent)
            int g = unify(sU, f, i);
            sU[g].rep = i;

        }
    }
    delete[] sU;

    gettimeofday(&end1, NULL);
    timeval_subtract(&result, &end1, &start1);
    serialtreeTime = ((result.tv_sec) * 1000 + result.tv_usec / 1000.0) + 0.5;
    printf("\nxxxxxxxxxxxx serial tree time: %ld xxxxxxxxxxxxxxx\n", serialtreeTime);
}

void printTimings() {
    cout << "Time split up......." << endl;
    cout << "Classify Critical points              : " << cpTime << "\t ms" << endl;
    cout << "Init Data structures (Par step)       : " << initParTime << "\t ms"
            << endl;
    cout << "find roots                            : " << rootTime << "\t ms"
            << endl;
    cout << "Sort critical points                  : " << sortTime << "\t ms"
            << endl;
    cout << "Time to build pairing heaps           : " << buildHeapTime << "\t ms"
            << endl;
    cout << "Compute join tree                     : " << joinTime << "\t ms"
            << endl;
    cout << "Compute Split tree                    : " << splitTime << "\t ms"
            << endl;
    cout << "Combine join & split trees            : " << ctTime << "\t ms" << endl;
    long par = cpTime + rootTime + sortTime + joinTime + splitTime + ctTime
            + initParTime + buildHeapTime;
    cout << "----------------------------------------" << endl;
    cout << "Time taken by parallel steps          : " << par << "\t ms" << endl;
    cout << endl;

    long seq = initHeapTime + completeTime + initSetupTime;
    cout << "Init Data structures (setup)          : " << initSetupTime << "\t ms"
            << endl;
    cout << "Init pairing heaps                    : " << initHeapTime << "\t ms"
            << endl;
    cout << "Complete join & split trees           : " << completeTime << "\t ms"
            << endl;
    cout << "----------------------------------------" << endl;
    cout << "Time taken by sequential steps        : " << seq << "\t ms" << endl;
    cout << "----------------------------------------" << endl;
    cout << "----------------------------------------" << endl;
    long tot = par + seq;
    cout << "Time taken to compute contour tree    : " << tot << "\t ms" << endl
            << endl;
}

int main(int argc, char **argv) {
    // Check if the required number of arguments are present
    if (argc != 7) {
        cout << "Usage: ./ct input_file dimx dimy dimz output_file_id subcube_no." << endl;
        //exit(1);
    }

    // Open file with parameters specified
    FILE* param = fopen("params.txt", "r");
    fscanf(param, "%d %d %d %d %d %d", &divx, &divy, &divz, &sub_sizex, &sub_sizey, &sub_sizez);
    fclose(param);

    // Find dimensions of the complete data
    global_dimx = divx * sub_sizex;
    global_dimy = divy * sub_sizey;
    global_dimz = divz * sub_sizez;

    // Store starting of computation time
    gettimeofday(&start1, NULL);

    // Process the arguments passed to the program
    char* file = argv[1];
    dimx = atoi(argv[2]);
    dimy = atoi(argv[3]);
    dimz = atoi(argv[4]);

    //int pid = atoi(argv[6]);
    int pid = 1;

    // Create a string for the filename from the arguments provided
    char filename[20];
    strcpy(filename, argv[5]);
    strcat(filename, ".dat");

    // Adhitya: Can this be better?
    subcube = atoi(argv[6]) - 1;

    //global_dim = atoi(argv[6]);

    //num_cores = atoi(argv[7]);
    //omp_set_num_threads(num_cores);

    // Computer number of vertices
    num_vert = dimx * dimy * dimz;

    long long int scxy = divx * divy;
    long long int sxy = subcube % scxy;
    long long int scx = sxy % divz;
    long long int scy = sxy / divz;
    long long int scz = subcube / scxy;

    // Find the extents
    // Adhitya: Find out what happens if the extents are greater than global dimensions
    // Adhitya: More interestingly, how does splitting take place? - It happens in subdiv.m
    ext.x[0] = scx * sub_sizex;
    ext.x[1] =
            (((scx + 1) * sub_sizex) < global_dimx) ?
            (scx + 1) * sub_sizex : ((scx + 1) * sub_sizex - 1);

    ext.y[0] = scy * sub_sizey;
    ext.y[1] =
            (((scy + 1) * sub_sizey) < global_dimy) ?
            (scy + 1) * sub_sizey : ((scy + 1) * sub_sizey - 1);

    ext.z[0] = scz * sub_sizez;
    ext.z[1] =
            (((scz + 1) * sub_sizez) < global_dimz) ?
            (scz + 1) * sub_sizez : ((scz + 1) * sub_sizez - 1);

    // Read the Input Grid file
    readInputGrid(file);

    cout << "No.  of vertices : " << num_vert << endl;

    gettimeofday(&end1, NULL);
    timeval_subtract(&result, &end1, &start1);
    mtime = ((result.tv_sec) * 1000 + result.tv_usec / 1000.0) + 0.5;
    printf("Time for reading input: %ld millitseconds\n", mtime);

    gettimeofday(&start, NULL);
    findCriticalPointsGrid();

    printf("------------no. of critical points: %d----\n", num_critical);
    find_roots();
    serialtree();
    cleanup_trees();
    gettimeofday(&end1, NULL);
    timeval_subtract(&result, &end1, &start);
    mtime = ((result.tv_sec) * 1000 + result.tv_usec / 1000.0) + 0.5;
    //printTimings();
    printf("\n\n Total time: %ld milliseconds\n critical points:%d\n\n", mtime,
            num_critical);
    //printf("\n touched: %d minima: %d, maxima %d saddles %d final %d\n",
    //touched_critical_points, num_minima, num_maxima, num_saddle, final);
    size_t csz = sizeof (char);
    size_t isz = sizeof (int);
    size_t fsz = sizeof (float);
    gettimeofday(&start, NULL);
    int dummy;
    int fp = open(filename, O_RDWR | O_CREAT | O_TRUNC, 0777);
    dummy = write(fp, (void*) (&num_vert), isz);
    dummy = write(fp, (void*) (&(ext.x)), 2 * isz);
    dummy = write(fp, (void*) (&(ext.y)), 2 * isz);
    dummy = write(fp, (void*) (&(ext.z)), 2 * isz);
    dummy = write(fp, (void*) (&num_critical), isz);
    dummy = write(fp, (void*) (function_values), num_critical * fsz);
    dummy = write(fp, (void*) (vertex_index), num_critical * isz);
    dummy = write(fp, (void*) (v_p), num_critical * sizeof (long long int));
    dummy = write(fp, (void*) (join_neigh), num_critical * isz);
    dummy = write(fp, (void*) (join_children), 2 * num_critical * isz);
    dummy = write(fp, (void*) (split_neigh), num_critical * isz);
    dummy = write(fp, (void*) (split_children), 2 * num_critical * isz);
    dummy = write(fp, (void*) (offset), num_critical * csz);
    //fsync(fp);
    close(fp);

    gettimeofday(&end1, NULL);
    timeval_subtract(&result, &end1, &start);
    mtime = ((result.tv_sec) * 1000 + result.tv_usec / 1000.0) + 0.5;

    printf("\n---------writing time------: %ld\n\n", mtime);

}

