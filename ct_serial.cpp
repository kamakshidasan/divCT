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
#include "ContourTree.h"
#include "const.h"
#include "PairingQueue.h"
#include "UnionFind.h"
#include "PSort.h"
#include "TriangleData.h"
#include "ReadFile.h"
#include "Graph.h"
using namespace std;
#define MAX_ADJ 2
int divx = 2; //no. of divisions along x-axis
int divy = 2;
int divz = 2;
int sub_sizex = 128; //subcube size
int sub_sizey = 128;
int sub_sizez = 128;

int global_dimx = divx*sub_sizex;
int global_dimy = divy*sub_sizey;
int global_dimz = divz*sub_sizez;

struct Graph* CTGraph;
char* offset_c;
int subcube;

struct extent {
    int x[2];
    int y[2];
    int z[2];
} ext;
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

struct Elems {
    int parent;
    int rank;
    int rep;
};

void ElemSet(struct Elems* elems, int i) {
    elems[i].parent = i;
    elems[i].rank = 0;
    elems[i].rep = i;
}

int findit(struct Elems* elems, int i) {
    if (elems[i].parent == i)
        return i;
    int k = i;
    while (elems[k].parent != k) {
        k = elems[k].parent;
    }
    while (elems[i].parent != i) {
        int t = elems[i].parent;
        elems[i].parent = k;
        i = t;
    }
    return k;
}

int unify(struct Elems* elems, int x, int y) {
    if (x == y)
        return x;
    int x_root = findit(elems, x);
    int y_root = findit(elems, y);
    if (x_root == y_root)
        return x_root;
    if (elems[x_root].rank > elems[y_root].rank) {
        elems[y_root].parent = x_root;
        return x_root;
    } else if (elems[y_root].rank > elems[x_root].rank) {
        elems[x_root].parent = y_root;
        return y_root;
    } else {
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

bool vertex_compare(float* vertices, int i, int j) {
    return ((vertices[i] < vertices[j]) || (vertices[i] == vertices[j] && i < j));
}

bool isGreater(float * vertices, int v1, int v2) {
    return !vertex_compare(vertices, v1, v2);
}

int F000 = 0;
int F001 = 4;
int F010 = 2;
int F011 = 6;
int F100 = 1;
int F101 = 5;
int F110 = 3;
int F111 = 7;

int XY = 1;
int YZ = 2;
int XZ = 4;
int XYZ = 8;
int XYZ2 = 16;

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

bool isBodyCP(int v, int dimx, int dimy, int dimz, int x, int y, int z,
        int dimxy, int nv, float * vertices) {
    bool mx[8];
    bool mn[8];

    if (x == dimx - 1 || y == dimy - 1 || z == dimz - 1) {
        return false;
    }
    int body[8];
    body[0] = v;
    body[1] = v + 1;
    body[2] = v + dimx;
    body[3] = body[2] + 1;

    body[4] = v + dimxy;
    body[5] = body[4] + 1;
    body[6] = body[4] + dimx;
    body[7] = body[6] + 1;

    int noMax = 0;
    int noMin = 0;
    for (int i = 0; i < 8; i++) {
        int mxct = 0;
        int mnct = 0;
        mx[i] = false;
        mn[i] = false;
        for (int j = 0; j < 3; j++) {
            int adj = bodyAdj[i * 3 + j];
            if (isGreater(vertices, body[i], body[adj])) {
                mxct++;
            } else {
                mnct++;
            }
        }
        if (mxct == 3) {
            //max
            noMax++;
            mx[i] = true;
        }
        if (mnct == 3) {
            noMin++;
            mn[i] = true;
        }
    }
    if (noMax == 4 || noMin == 4) {
        // TODO There are two body saddles. for now ignoring.
        return false;
    }
    if (!(noMax == 2 && noMin == 1 || noMax == 1 && noMin == 2)) {
        // no body saddle.
        return false;
    }

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
        face[1] = v + 1;
        face[2] = v + dimx;
        face[3] = face[2] + 1;
    }
    if (axis == XZ) {
        face[1] = v + 1;
        face[2] = v + dimxy;
        face[3] = face[2] + 1;
    }

    if (axis == YZ) {
        face[1] = v + dimx;
        face[2] = v + dimxy;
        face[3] = face[2] + dimx;
    }

    if (isGreater(vertices, face[0], face[1])
            && isGreater(vertices, face[0], face[2])
            && isGreater(vertices, face[3], face[1])
            && isGreater(vertices, face[3], face[2])) {
        return true;
    }

    if (isGreater(vertices, face[1], face[0])
            && isGreater(vertices, face[2], face[0])
            && isGreater(vertices, face[1], face[3])
            && isGreater(vertices, face[2], face[3])) {
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
            // we dont have to do it again in this condition
            if (bnd[0] != -1) { // If the point is on the boundary
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
                else {
                    // Otherwise, it is just another regular point with exactly one different components
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
        face[1] = v + 1;
        face[2] = v + dimx;
        face[3] = face[2] + 1;
    }
    if (axis == FS_XZ) {
        face[1] = v + 1;
        face[2] = v + dimxy;
        face[3] = face[2] + 1;
    }

    if (axis == FS_YZ) {
        face[1] = v + dimx;
        face[2] = v + dimxy;
        face[3] = face[2] + dimx;
    }

    if (isGreater(vertices, face[0], face[1]) && isGreater(vertices, face[0], face[2]) && isGreater(vertices, face[3], face[1]) && isGreater(vertices, face[3], face[2])) {
        lower_root_count[vIndex] = 2;
        upper_root_count[vIndex] = 2;
        lower_roots[vIndex * MAX_ADJ + 0] = face[1];
        lower_roots[vIndex * MAX_ADJ + 1] = face[2];
        upper_roots[vIndex * MAX_ADJ + 0] = face[0];
        upper_roots[vIndex * MAX_ADJ + 1] = face[3];

        if (isGreater(vertices, face[1], face[2])) {
            vertex_pos[vIndex] = face[1];
        } else {
            vertex_pos[vIndex] = face[2];
        }
        offset[vIndex] = 1;
    }

    if (isGreater(vertices, face[1], face[0]) && isGreater(vertices, face[2], face[0]) && isGreater(vertices, face[1], face[3]) && isGreater(vertices, face[2], face[3])) {
        lower_root_count[vIndex] = 2;
        upper_root_count[vIndex] = 2;
        lower_roots[vIndex * MAX_ADJ + 0] = face[0];
        lower_roots[vIndex * MAX_ADJ + 1] = face[3];
        upper_roots[vIndex * MAX_ADJ + 0] = face[1];
        upper_roots[vIndex * MAX_ADJ + 1] = face[2];

        if (isGreater(vertices, face[0], face[3])) {
            vertex_pos[vIndex] = face[0];
        } else {
            vertex_pos[vIndex] = face[3];
        }
        offset[vIndex] = 1;
    }

    {
        float corners[4];
        for (int i = 0; i < 4; i++) {
            corners[i] = vertices[face[i]];
        }
        float a = corners[0] + corners[3] - corners[1] - corners[2];
        if (a == 0) {
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
        function_values[vIndex] = a * x * y + b * x + c * y + corners[0];
    }
}

int getIndx(float* vertices, float fsad, int lowMax, int highMin) {
    float fmax = vertices[lowMax];
    if (fsad == fmax) {
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

    int vIndex1 = vIndex;

    if (x == dimx - 1 || y == dimy - 1 || z == dimz - 1) {
        return;
    }
    int body[8];
    body[0] = v;
    body[1] = v + 1;
    body[2] = v + dimx;
    body[3] = body[2] + 1;

    body[4] = v + dimxy;
    body[5] = body[4] + 1;
    body[6] = body[4] + dimx;
    body[7] = body[6] + 1;

    int noMax = 0;
    int noMin = 0;
    int highMin, lowMax;
    for (int i = 0; i < 8; i++) {
        int mxct = 0;
        int mnct = 0;
        mx[i] = false;
        mn[i] = false;
        for (int j = 0; j < 3; j++) {
            int adj = bodyAdj[i * 3 + j];
            if (isGreater(vertices, body[i], body[adj])) {
                mxct++;
            } else {
                mnct++;
            }
        }
        if (mxct == 3) {
            //max
            noMax++;
            mx[i] = true;
            if (noMax == 1 || isGreater(vertices, lowMax, body[i])) {
                lowMax = body[i];
            }
        }
        if (mnct == 3) {
            noMin++;
            mn[i] = true;
            if (noMin == 1 || vertex_compare(vertices, highMin, body[i])) {
                highMin = body[i];
            }
        }
    }
    if (noMax == 4 || noMin == 4) {
        // TODO There are two body saddles. for now ignoring.
        return;
    }
    if (!(noMax == 2 && noMin == 1 || noMax == 1 && noMin == 2)) {
        // no body saddle.
        return;
    }

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

        float corners[8];

        for (int i = 0; i < 8; i++) {
            corners[i] = vertices[body[i]];
        }

        float a = -corners[F000] + corners[F001] + corners[F010] - corners[F011]
                + corners[F100] - corners[F101] - corners[F110] + corners[F111];
        float b = corners[F000] - corners[F010] - corners[F100] + corners[F110];
        float c = corners[F000] - corners[F001] - corners[F010] + corners[F011];
        float d = corners[F000] - corners[F001] - corners[F100] + corners[F101];
        float e = -corners[F000] + corners[F100];
        float f = -corners[F000] + corners[F010];
        float g = -corners[F000] + corners[F001];
        float h = corners[F000];

        if (noMax == 2) {
            function_values[vIndex] = vertices[highMin];
        } else {
            function_values[vIndex] = vertices[lowMax];
        }
        if (a == 0) {
            if (b == 0 || c == 0 || d == 0) {
                vertex_pos[vIndex] =
                        getIndx(vertices, function_values[vIndex1], lowMax,
                        highMin
                        );
                return;
            }
            float x = (e * c - d * f - b * g) / (2 * b * d);
            if (x < 0 || x > 1) {
                vertex_pos[vIndex] =
                        getIndx(vertices, function_values[vIndex1], lowMax,
                        highMin
                        );
                return;
            }
            float y = (-e * c + d * f - b * g) / (2 * b * c);
            if (y < 0 || y > 1) {
                vertex_pos[vIndex] =
                        getIndx(vertices, function_values[vIndex1], lowMax,
                        highMin
                        );
                return;
            }
            float z = (-e * c - d * f + b * g) / (2 * c * d);
            if (z < 0 || z > 1) {
                vertex_pos[vIndex] =
                        getIndx(vertices, function_values[vIndex1], lowMax,
                        highMin
                        );
                return;
            }
            function_values[vIndex1] = a * x * y * z + b * x * y + c * y * z
                    + d * z * x + e * x + f * y + g * z + h;
            vertex_pos[vIndex] =
                    getIndx(vertices, function_values[vIndex1], lowMax,
                    highMin
                    );
        } else {
            if ((a * e - b * d) == 0) {
                vertex_pos[vIndex] = getIndx(vertices, function_values[vIndex1], lowMax, highMin);
                return;
            }
            bool ret1 = true;
            bool ret2 = true;
            float r = (b * c - a * f) * (a * g - c * d) / (a * e - b * d);
            if (r < 0) {
                vertex_pos[vIndex] = getIndx(vertices, function_values[vIndex1], lowMax, highMin);
                return;
            }

            r = (float) sqrt(r);
            float x = (-c + r) / a;
            if (x < 0 || x > 1) {
                ret1 = false;
            }
            float y = -(d * x + g) / (a * x + c);
            if (y > 1 || y < 0 || (a * x + c) == 0) {
                ret1 = false;
            }
            float z = -(b * x + f) / (a * x + c);
            if (z > 1 || z < 0 || (a * x + c) == 0) {
                ret1 = false;
            }
            if (ret1 != false) {
                function_values[vIndex1] = a * x * y * z + b * x * y + c * y * z + d * z * x + e * x + f * y + g * z + h;
                vertex_pos[vIndex] = getIndx(vertices, function_values[vIndex1], lowMax, highMin);
                return;
            }

            x = (-c - r) / a;
            if (x < 0 || x > 1) {
                ret2 = false;
            }
            y = -(d * x + g) / (a * x + c);
            if (y > 1 || y < 0 || (a * x + c) == 0) {
                ret2 = false;
            }
            z = -(b * x + f) / (a * x + c);
            if (z > 1 || z < 0 || (a * x + c) == 0) {
                ret2 = false;
            }
            if (ret2 != false) {
                function_values[vIndex1] = a * x * y * z + b * x * y + c * y * z + d * z * x + e * x + f * y + g * z + h;
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

        //printf("\nnumcritical:%d, %d\n",cp,num_critical);
        lower_root_count[cp] = 0;
        upper_root_count[cp] = 0;
        if (v < num_vert) {
            // vertex critical
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
            // get the vertex from each component so that one can traverse
            for (int i = 0; i < 6; i++) {
                if (done[i] == 1) {
                    lower_roots[cp * MAX_ADJ + 0] = verts[i];
                }
                if (done[i] == 2) {
                    lower_roots[cp * MAX_ADJ + 1] = verts[i];
                }
                if (done[i] == 3) {
                    upper_roots[cp * MAX_ADJ + 0] = verts[i];
                }
                if (done[i] == 4) {
                    upper_roots[cp * MAX_ADJ + 1] = verts[i];
                }
            }
        } else {
            // its a body or face saddle

            for (int i = 1; i <= 3; i++) {
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

            int dimxy = dimx * dimy;
            int z = v / dimxy;
            int xy = v % dimxy;
            int y = xy / dimx;
            int x = xy % dimx;

            if ((v < num_vert)&&(v >= 0)) updateBodyCP(cp, v, dimx, dimy, dimz, x, y, z, dimxy, num_vert, vertices, function_values,
                    lower_roots, upper_roots, lower_root_count, upper_root_count, vertex_pos, offset);
        }
        //printf("helllllll\n");
        //printf("\nnumcritical:%d, %d\n",cp,num_critical);
    }

    delete[] vertices;
}

int find_min_root(int i, int *steep_desc) {
    int k = i;
    while (k != steep_desc[k]) {

        k = steep_desc[k];

        //if (critical_points[k]>num_vert) break;
        if (type[k] == SADDLE) break;
    }
    return k;
}

int find_max_root(int i, int *steep_asc) {
    int k = i;
    while (k != steep_asc[k]) {

        k = steep_asc[k];
        //if (critical_points[k]>num_vert) break;
        if (type[k] == SADDLE) break;
    }
    return k;
}

void find_roots() {
    gettimeofday(&start1, NULL);
    for (int i = 0; i < num_critical; i++) {
        int c = i*MAX_ADJ;
        int c1 = c;
        for (int j = 0; j < lower_root_count[i]; j++) {
            //lower_roots[i*MAX_ADJ+j] = extrema_map[find_min_root(lower_roots[i*MAX_ADJ+j], steep_desc)];
            lower_roots[c] = find_min_root(lower_roots[c], steep_desc);
            lower_roots[c] = extrema_map[lower_roots[c]];
            c++;
        }
        c = c1;
        for (int j = 0; j < upper_root_count[i]; j++) {
            // upper_roots[i*MAX_ADJ+j] = extrema_map[find_max_root(upper_roots[i*MAX_ADJ+j], steep_asc)];
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
        } else if (type[i] == MAXIMUM) {
            num_maxima++;
            leaves[leavesCt++] = cpNo;
            critical_points[cpNo] = i;
            extrema_map[i] = cpNo;
            cpNo++;
        } else if (type[i] == SADDLE) {
            num_saddle++;
            critical_points[cpNo] = i;
            extrema_map[i] = cpNo;
            cpNo++;
        }
        if ((bfType[i] & XY) != 0) {
            num_saddle++;
            critical_points[cpNo] = (i + xyIndex);
            cpNo++;
        }
        if ((bfType[i] & YZ) != 0) {
            num_saddle++;
            critical_points[cpNo] = (i + yzIndex);
            cpNo++;
        }
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
    //	printf("Time to initialize data structures: %ld milliseconds\n", mtime);
    initParTime = mtime;
}

void findCriticalPointsGrid() {

    classifyCriticalPoints();

    initializeDataStructures();
}

void cleanup_trees() {
    gettimeofday(&start1, NULL);
    int critical = 0;
    int b1 = dimx * dimy;
    int * x = new int[num_critical];
    int * is_critical = new int[num_critical];
    for (int i = 0; i < num_critical; i++) {
        int j = i; //check
        is_critical[j] = 0;
        long long int vpos = vertex_pos[j];
        long long int dimxy = dimx * dimy;
        long long int posz = vpos / dimxy;
        long long int xy = vpos % dimxy;
        long long int posy = xy / dimx;
        long long int posx = xy % dimx;
        if ((join_children[2 * j] == -1) || (split_children[2 * j] == -1)
                || (join_children[2 * j + 1] != -1)
                || (split_children[2 * j + 1] != -1) || (posx == 0)
                || (posx == dimx - 1) || (posy == 0) || (posy == dimy - 1)
                || (posz == dimz - 1) || (posz == 0)) {
            x[j] = critical;
            is_critical[j] = 1;
            critical++;
        }
    }
    final = critical;

    printf("\n------ new critical: %d---------\n", final);
    j_n = new int[critical];
    s_n = new int[critical];
    j_c = new int[2 * critical];
    s_c = new int[2 * critical];
    v_i = new int[critical];
    v_p = new long long int[critical];
    f_v = new float[critical];
    //u_d = new int[critical];
    //l_d = new int[critical];
    offset_c = new char[critical];
    for (int i = 0; i < 2 * critical; i++) {
        j_c[i] = s_c[i] = -1;
        //u_d[i / 2] = l_d[i / 2] = 0;
    }

    int c = 0;
    //int l=0;//leaf count

    //assuming cubic division
    /*long long int sub_size = dimx;
    if((dimx==dimy) && (dimy==dimz) && ((global_dim%dimx) == 0)) sub_size = dimx;
    else if ((dimx==dimy) && (dimy==dimz)) sub_size = dimx-1;
    else if (dimy < dimx)
            sub_size = dimy;
    else if (dimz < dimx)
            sub_size = dimz;*/
    /*long long int divx = (global_dimx / dimx);
    long long int divy = (global_dimy / dimy);
    long long int divz = (global_dimz / dimz);

     */
    long long int scxy = divx * divy;
    long long int scz = subcube / scxy;
    long long int sxy = subcube % scxy;
    long long int scy = sxy / divx;
    long long int scx = sxy % divx;
    int leaf_count = 0;
    for (int i = 0; i < num_critical; i++) {
        int j = vertex_index[i];
        int t = join_neigh[j];

        if ((join_children[2 * j] == -1) || (split_children[2 * j] == -1))
            leaves[leaf_count++] = x[j];

        if (is_critical[j] == 1) {
            v_i[c] = x[j]; //check
            f_v[x[j]] = function_values[j];

            long long int vpos = vertex_pos[j];
            long long int dimxy = dimx * dimy;
            long long int posz = vpos / dimxy;
            long long int xy = vpos % dimxy;
            long long int posy = xy / dimx;
            long long int posx = xy % dimx;

            long long int final_pos = (scz * sub_sizez + posz) * global_dimx
                    * global_dimy + (scy * sub_sizey + posy) * global_dimx
                    + (scx * sub_sizex + posx);

            v_p[x[j]] = final_pos;
            offset_c[x[j]] = offset[j];
            while ((t != -1) && (is_critical[t] != 1)) {
                t = join_neigh[t];
            }
            if (t == -1) {
                j_n[x[j]] = -1;

            } else {
                j_n[x[j]] = x[t];
                if (j_c[2 * x[t]] == -1)
                    j_c[2 * x[t]] = x[j];
                else
                    j_c[2 * x[t] + 1] = x[j];
                //l_d[x[t]]++;
            }
            t = split_neigh[j];
            while ((t != -1) && (is_critical[t] != 1))
                t = split_neigh[t];
            if (t == -1) {
                s_n[x[j]] = -1;

            } else {
                s_n[x[j]] = x[t];
                if (s_c[2 * x[t]] == -1)
                    s_c[2 * x[t]] = x[j];
                else
                    s_c[2 * x[t] + 1] = x[j];
                //u_d[x[t]]++;
            }
            c++;

        }
    }
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
    //vertex_pos = v_p;
    //upper_degree = u_d;
    //lower_degree = l_d;
    offset = offset_c;
    //isprocessed,uppdegree,lowerdegree
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

void serialtree() {

    //testRoots();
    sortVertices();
    gettimeofday(&start1, NULL);
    struct Elems *jU = new struct Elems[num_critical];
    //int * Rep = new int[num_critical];
    for (int s = 0; s < num_critical; s++) {
        int i = vertex_index[s];
        ElemSet(jU, i);
        for (int j = 0; j < lower_root_count[i]; j++) {

            int k = i * MAX_ADJ + j;
            int rt = lower_roots[k];
            int f = findit(jU, rt);
            int r = jU[f].rep;

            if (r == i) continue;
            join_neigh[r] = i;
            if (join_children[2 * i] == -1)
                join_children[2 * i] = r;
            else
                join_children[2 * i + 1] = r;
            int g = unify(jU, f, i);
            jU[g].rep = i;

        }
    }
    delete[] jU;
    struct Elems *sU = new struct Elems[num_critical];
    //int * Rep = new int[num_critical];

    for (int s = 0; s < num_critical; s++) {
        int i = vertex_index[num_critical - s - 1];
        ElemSet(sU, i);
        for (int j = 0; j < upper_root_count[i]; j++) {

            int k = i * MAX_ADJ + j;
            int rt = upper_roots[k];
            //if(upper_root_count[i]==2) printf("\nroot%d : %d\n",j+1,rt);

            int f = findit(sU, rt);
            int r = sU[f].rep;
            //if(upper_root_count[i]==2) printf("\nroot%d : %d r: %d\n",j+1,rt,r);
            if (r == i) continue;
            //if(f==findit(sU,i)) continue;
            split_neigh[r] = i;
            if (split_children[2 * i] == -1)
                split_children[2 * i] = r;
            else
                split_children[2 * i + 1] = r;
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
        //./serialct Data/Fuel/16/fo_1_33_33_17.raw 33 33 17 f1 1
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

    // Adhitya: The below two lines are definitely not correct going by the error display
    //num_cores = atoi(argv[5]);
    //int pid = atoi(argv[6]);
    int pid = 1;

    // Create a string for the filename from the arguments provided
    char filename[20];
    strcpy(filename, argv[5]);
    strcat(filename, ".dat");

    // Adhitya: Can this be better?
    subcube = atoi(argv[6]) - 1;

    //global_dim = atoi(argv[6]);
    //omp_set_num_threads(num_cores);

    // Computer number of vertices
    num_vert = dimx * dimy * dimz;

    /* long long int sub_size = dimx;
     if((dimx==dimy) && (dimy==dimz) && ((global_dim%dimx) == 0)) sub_size = dimx;
     else if ((dimx==dimy) && (dimy==dimz)) sub_size = dimx-1;
     else if (dimy < dimx)
            sub_size = dimy;
    else if (dimz < dimx)
            sub_size = dimz;
     */


    long long int scxy = divx * divy;
    long long int sxy = subcube % scxy;
    long long int scx = sxy % divz;
    long long int scy = sxy / divz;
    long long int scz = subcube / scxy;

    // Find the extents
    // Adhitya: Find out what happens if the extents are greater than global dimensions
    // Adhitya: More interestingly, how does splitting take place?
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
    int touched_critical_points = 0;
    printf("\n\n Total time: %ld milliseconds\n critical points:%d\n\n", mtime,
            num_critical);
    //printf("\n touched: %d minima: %d, maxima %d saddles %d final %d\n",
    //touched_critical_points, num_minima, num_maxima, num_saddle, final);
    int fp, fp2;
    size_t csz = sizeof (char);
    size_t isz = sizeof (int);
    size_t fsz = sizeof (float);
    gettimeofday(&start, NULL);
    int dummy;
    fp = open(filename, O_RDWR | O_CREAT | O_TRUNC, 0777);
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

    //std::cin.ignore();

}

