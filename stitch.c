#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <omp.h>
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
enum {
	false, true
};
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
void Union(int child, int parent, int * Uk) {
	if (child == parent)
		return;
	Uk[child] = parent;
	Uk[parent] = -1;
	//printf("union\n");
}
int Find(int elem, int* Uk) {
	if (Uk[elem] == -2)
		return -2;

	int t_elem = elem;
	int tt_elem;
	while (Uk[elem] != -1)
		elem = Uk[elem];

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
void seqmerge() { //standard code to merge two sorted sequence of elements

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
		} else {
			temp[i] = a[m];
			own[i] = 1;
			m++;
		}
		i++;
	}
	//  printf("\n.%d..%d...%d...\n",vp[a[1]],vp[a[n1]],lessVertex(a[n1],a[1]));
//printf("l: %d, m: %d i: %d",l,m,i);
	if (l > mid) {
		for (k = m; k <= high; k++) {
			temp[i] = a[k];
			own[i] = 1;
			i++;
		}
	} else {
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

	//printf("inside cleanup---- n1:%d,n2:%d,b1:%d,b2:%d\n",n1,n2,b1,b2);

	int i, j, t;

	int num_critical = n1 + n2;
	int * x = (int*) malloc(num_critical * sizeof(int));
	char * is_critical = (char*) malloc(num_critical * sizeof(char));
	v_i = (int*) malloc(num_critical * sizeof(int));

	critical = 0;

	//#pragma omp parallel for schedule(dynamic)
	for (i = 0; i < num_critical; i++) {
		j = temp[i]; //check
		is_critical[j] = 0;
		/*if (c[2 * j] == j) {
			printf("\nself-child\n");
			c[2 * j] = -1;
		}
		if (c[2 * j + 1] == j) {
			printf("\nself-child\n");
			c[2 * j + 1] = -1;
		}
		if (c_[2 * j] == j) {
			printf("\nself-child\n");
			c_[2 * j] = -1;
		}
		if (c_[2 * j + 1] == j) {
			printf("\nself-child\n");
			c_[2 * j + 1] = -1;
		}
		if ((c[2 * j] == c[2 * j + 1])&&(c[2 * j]!=-1)) {
			printf("\ndouble-child\n");
			c[2 * j + 1] = -1;
		}
		if ((c_[2 * j] == c_[2 * j + 1])&&(c_[2 * j]!=-1)) {
			printf("\ndouble-child\n");
			c_[2 * j + 1] = -1;
		}*/


		long long int vpos = vp[j];
		long long int dimxy = dimx * dimy;
				long long int posz = vpos / dimxy;
				long long int xy = vpos % dimxy;
				long long int posy = xy / dimx;
				long long int posx = xy % dimx;

		if ((c[2 * j] == -1 && c[2 * j + 1] == -1)
				|| (c_[2 * j] == -1 && c_[2 * j + 1] == -1)
				|| (c[2 * j + 1] != -1 && c[2 * j] != -1)
				|| (c_[2 * j + 1] != -1 && c_[2 * j] != -1)
				|| (posx == ext[0])
				|| (posx == ext[1])
								|| (posy == ext[2])
								|| (posy == ext[3])
								|| (posz == ext[4])
								|| (posz == ext[5])){
			x[j] = critical;
			is_critical[j] = 1;
			v_i[critical] = critical;

			critical++;
		}

	}

	j_n = (int*) malloc(critical * sizeof(int));
	s_n = (int*) malloc(critical * sizeof(int));
	j_c = (int*) malloc(2 * critical * sizeof(int));
	s_c = (int*) malloc(2 * critical * sizeof(int));
	v_i = (int*) realloc(v_i,critical * sizeof(int));
	v_p = (long long int*) malloc(critical * sizeof(long long int));
	f_v = (float*) malloc(critical * sizeof(float));
	offset = (char*) malloc(critical * sizeof(char));

//u_d=new int[critical];
//l_d=new int[critical];

	#pragma omp parallel for schedule(dynamic) private(i)
	for (i = 0; i < critical; i++) {
		j_c[2 * i] = j_c[2 * i + 1] = s_c[2 * i] = s_c[2 * i + 1] = j_n[i] =
				s_n[i] = -1;
		//	u_d[i/2]=l_d[i/2]=0;
	}

	int cl = 0;
	int l = 0; //leaf count

	#pragma omp parallel for schedule(dynamic) private(i,j,t)
	for (i = 0; i < num_critical; i++) {
		j = temp[i];
		t = p[j];

		//if ((join_children[2*j]==-1)||(split_children[2*j]==-1)) leaves[l++]=x[j];

		if (is_critical[j] == 1) {
			//v_i[cl] = x[j]; //check
			f_v[x[j]] = fv[j];
			v_p[x[j]] = vp[j];
			offset[x[j]] = off[j];
			while ((t != -1) && (is_critical[t] != 1)) {
				t = p[t];
			}
			if (t == -1) {
				j_n[x[j]] = -1;

			} else {
				j_n[x[j]] = x[t];
				if (j_c[2 * x[t]] == -1)
					j_c[2 * x[t]] = x[j];
				else
					j_c[2 * x[t] + 1] = x[j];
				//u_d[x[t]]++;
			}
			t = p_[j];
			while ((t != -1) && (is_critical[t] != 1))
				t = p_[t];
			if (t == -1) {
				s_n[x[j]] = -1;

			} else {
				s_n[x[j]] = x[t];
				if (s_c[2 * x[t]] == -1)
					s_c[2 * x[t]] = x[j];
				else
					s_c[2 * x[t] + 1] = x[j];
				//l_d[x[t]]++;
			}
			//cl++;

		}
	}
	//int temp = num_critical;
	/*join_neigh=j_n;
	 split_neigh=s_n;

	 num_critical=critical;
	 vertex_index = v_i;
	 function_values = f_v;
	 vertex_pos = v_p;
	 upper_degree = u_d;
	 lower_degree = l_d;*/
	//isprocessed,uppdegree,lowerdegree
	final = critical;
	free(p);
	free(p_);
	free(c);
	free(c_);

	free(vp);
	free(fv);
free(x);
free(is_critical);
	//printf("cleanup tree: critical: %d. critical_prev: %d\n", critical,			num_critical);

	/*for( i=0;i<critical;i++)
	 {if (j_n[i]>critical)
	 printf("\nnanana\n")
	 ;
	 }*/


}

int min(int x, int y) {
	if (x < y)
		return x;
	else
		return y;
}
int max(int x, int y) {
	if (x < y)
		return y;
	else
		return x;
}

int main(int argc, char **argv) {

	//to be read from file

//omp_set_num_threads(1);

	//vn1 = 389 * 133 * 66;
	
	FILE* param = fopen("params.txt","w");
fscanf(param,"%d %d %d %d %d %d",&divx,&divy,&divz,&sub_sizex,&sub_sizey,&sub_sizez);
fclose(param);
dimx = divx*sub_sizex;
dimy = divy*sub_sizey;
dimz = divz*sub_sizez;
	
int tree;
if (argc == 5)
		tree = atoi(argv[4]);

	char file1[30], file2[30],file[30];
	strcpy(file1, argv[1]);
	strcat(file1, ".dat");
	strcpy(file2, argv[2]);
	strcat(file2, ".dat");
strcpy(file, argv[3]);
strcat(file, ".dat");

	/*int vp1[6] = {100,200,2,3,400,500};
	 int vp2[6] = {10,20,2,3,40,50};
	 int p1[6] = {-1,0,1,1,2,3};
	 int p2[6] = {-1,0,1,1,2,3};
	 int c1[12] = {1,-1,2,3,4,-1,5,-1,-1,-1,-1,-1};
	 int c2[12] = {1,-1,2,3,4,-1,5,-1,-1,-1,-1,-1};*/

	//to be obtained after merging f1 and f2
	//int S[12] = {10,11,5,4,3,9,2,8,7,1,0,6};
	size_t isz = sizeof(int);
	size_t csz = sizeof(char);
	size_t fsz = sizeof(float);
	int ext1[6];
	int ext2[6];

	int fp1, fp2;
	fp1 = open(file1, O_RDONLY);
	fp2 = open(file2, O_RDONLY);

	gettimeofday(&tv1, NULL);

	read(fp1, (void*) (&vn1), isz); //fscanf(fp1,"%d",&n1);
	read(fp2, (void*) (&vn2), isz); //fscanf(fp2,"%d",&n2);
	read(fp1, (void*) (ext1), 6 * isz);
	read(fp2, (void*) (ext2), 6 * isz);

	//extents
	ext[0] = min(ext1[0], ext2[0]); //find out the extents
	ext[1] = max(ext1[1], ext2[1]);
	ext[2] = min(ext1[2], ext2[2]);
	ext[4] = min(ext1[4], ext2[4]);
	ext[3] = max(ext1[3], ext2[3]);
	ext[5] = max(ext1[5], ext2[5]);

	read(fp1, (void*) (&n1), isz);
	read(fp2, (void*) (&n2), isz);
	printf("n1:%d,n2:%d\n", n1, n2);

//exit(0);

	fv = (float*) malloc((n1 + n2) * sizeof(float)); //function_values
	read(fp1, (void*) fv, n1 * fsz); 
	read(fp2, (void*) (fv + n1), n2 * fsz);

	vi = (int*) malloc((n1 + n2) * sizeof(int)); //vertex_index
	read(fp1, (void*) vi, n1 * isz);
	read(fp2, (void*) (vi + n1), n2 * isz);

	vp = (long long int*) malloc((n1 + n2) * sizeof(long long int));//vertex_pos
	read(fp1, (void*) vp, n1 * sizeof(long long int));
	read(fp2, (void*) (vp + n1), n2 * sizeof(long long int));

	p = (int*) malloc((n1 + n2) * sizeof(int));//join tree parent array
	read(fp1, (void*) p, n1 * isz);
	read(fp2, (void*) (p + n1), n2 * isz);

	c = (int*) malloc(2 * (n1 + n2) * sizeof(int));//join tree children array j[2i] and j[2i+1] are the two children of ith critical point
	read(fp1, (void*) c, 2 * n1 * isz);
	read(fp2, (void*) (c + 2 * n1), 2 * n2 * isz);

	p_ = (int*) malloc((n1 + n2) * sizeof(int));//similar arrays for split tree
	read(fp1, (void*) p_, n1 * isz);
	read(fp2, (void*) (p_ + n1), n2 * isz);

	c_ = (int*) malloc(2 * (n1 + n2) * sizeof(int));
	read(fp1, (void*) c_, 2 * n1 * isz);
	read(fp2, (void*) (c_ + 2 * n1), 2 * n2 * isz);

	off = (char*) malloc((n1 + n2) * sizeof(char));//offset to indicate if a point is face or body saddle
	read(fp1, (void*) off, n1 * csz);
	read(fp2, (void*) (off + n1), n2 * csz);

	close(fp1);
	close(fp2);

	//remove(file1);
	//remove(file2);

	U = (int*) malloc((n1 + n2) * sizeof(int));//UF for join tree
	U_ = (int*) malloc((n1 + n2) * sizeof(int));//UF for split tree

	own = (int*) malloc((n1 + n2) * sizeof(int));
	temp = (int*) malloc((n1 + n2) * sizeof(int));

	gettimeofday(&tv2, NULL);

	printf("\n alloc and read time = %f miliseconds\n",(double) (tv2.tv_usec - tv1.tv_usec) / 1000+ (double) (tv2.tv_sec - tv1.tv_sec) * 1000);



	int i, j, k;

	gettimeofday(&tv1, NULL);

	#pragma omp parallel for schedule(dynamic) private(i)
	for (i = 0; i < (n1 + n2); i++) {
		U[i] = -2;
		U_[i] = -2;

		if (i < n1) {
			/*fscanf(fp1,"%f",&fv[i]);
			 fscanf(fp1,"%d",&vi[i]);
			 fscanf(fp1,"%d",&vp[i]);//vp[i] = vp1[i];
			 fscanf(fp1,"%d",&p[i]);//p[i] = p1[i];
			 fscanf(fp1,"%d",&c[2*i]);//c[2*i]=c1[2*i];
			 fscanf(fp1,"%d",&c[2*i+1]);//c[2*i+1]=c1[2*i+1];
			 fscanf(fp1,"%d",&off[i]);*/

		} else {
			j = i - n1;
			//fscanf(fp2,"%f",&fv[i]);
			//fscanf(fp2,"%d",&vi[i]);
			vi[i] = vi[i] + n1;
			//fscanf(fp2,"%d",&vp[i]);//vp[i] = vp2[j];
			//vp[i] = vp[i] + vn1 - b1;
			//fscanf(fp2,"%d",&p[i]);//p[i] = (p2[j]==-1)?-1:n1+p2[j];
			p[i] = (p[i] == -1) ? -1 : n1 + p[i];
			p_[i] = (p_[i] == -1) ? -1 : n1 + p_[i];

			//fscanf(fp2,"%d",&c[2*i]);//c[2*i] = (c2[2*j]==-1)?-1:(n1+c2[2*j]);
			c[2 * i] = (c[2 * i] == -1) ? -1 : (n1 + c[2 * i]);
			c_[2 * i] = (c_[2 * i] == -1) ? -1 : (n1 + c_[2 * i]);
			//fscanf(fp2,"%d",&c[2*i+1]);//c[2*i+1] = (c2[2*j+1]==-1)?-1:(n1+c2[2*j+1]);
			c[2 * i + 1] = (c[2 * i + 1] == -1) ? -1 : (n1 + c[2 * i + 1]);
			c_[2 * i + 1] = (c_[2 * i + 1] == -1) ? -1 : (n1 + c_[2 * i + 1]);
			//fscanf(fp2,"%d",&off[i]);
		}
		//printf("%d\n",U[i]);
	}

	//for (i=0;i<(n1+n2);i++){if (fv[i]==255) printf("\nvi>n1+n2:%f\n",fv[i]);}

	critical = n1 + n2;
	j_n = p;
	s_n = p_;
	j_c = c;
	s_c = c_;

	v_p = vp;
	f_v = fv;
	offset = off;

	/*for (i=0;i<critical;i++){
	 if (j_n[i]==-1) printf("\njoin root:%d,value:%f\n",i,f_v[i]);
	 if (s_n[i]==-1) printf("\nsplit root:%d,value:%f\n",i,f_v[i]);

	 }*/

	seqmerge();
	v_i = temp;
	//printf("\nhere\n");

	int * S = temp;
	int jmp = 0;

#pragma omp parallel sections private(i,k)
{
	#pragma omp section
	{
	for (i = 1; i < n1 + n2; i++) {
		int z = S[i];
		int zm = S[i-1];
		//if(vp[z]==33682195) printf("\n,,,,value:%f,,offset:%d,,own:%d,,,index:%d,\n",fv[z],off[z],own[i],temp[i]);
		//	printf("\n%d\n",i);
		//if(z>=n1) {printf("damn: %d, z: %d %d %d %f\n",i,z,S[1],S[2],fv[S[n1+1]]);}
		//if(own[i]==1){printf("\n own1 %d func: %f \n",i,fv[z]);break;};
		if ((vp[z] == vp[zm]) && (off[z] == 0)
				&& (off[zm] == 0) && (own[i] != own[i - 1])) {

			//printf("\nbrdr\n");
			//if(jmp==0){printf("hello: pos: %d,val: %f,index: %d\n\n",vp[z],fv[z],temp[i]);jmp = 2;}
			p[zm] = z;
			if (c[2 * z] == -1)
				c[2 * z] = zm;
			else
				if (zm != c[2 * z])
					c[2 * z + 1] = zm;
			//printf("dup:%d\n",i);
			Union(zm, z, U);
		}
		if (c[2 * z] != -1) {
			if (U[c[2 * z]] != -2) {
				//printf("find1:%d U: %d\n",i,U[c[2*z]]);
				k = Find(c[2 * z], U);
				if (k != z) {
					Union(k, z, U);
					p[k] = z;
					c[2 * z ] = k;
					if (c[2 * z + 1] == k)
						c[2 * z + 1] = -1;
				}
			}
		}
		if (c[2 * z + 1] != -1) {
			if (U[c[2 * z + 1]] != -2) {
				//printf("find2:%d\n",i);
				k = Find(c[2 * z + 1], U);
				if (k != z) {
					Union(k, z, U);
					p[k] = z;
					if (c[2 * z ] != k)
						c[2 * z + 1] = k;
					else
						c[2 * z + 1] = -1;
				}
			}
		}

	}
	free(U);
		}//omp section 1 ends



	#pragma omp section
{
	for (i = (n1 + n2 - 2); i >= 0; i--) {
		//	printf("\n%d\n",i);
		//if(S[i]>=n1) {printf("damn: %d, S[i]: %d %d %d %f\n",i,S[i],S[1],S[2],fv[S[n1+1]]);}
		//if(own[i]==1){printf("\n own1 %d func: %f \n",i,fv[S[i]]);break;};
		int z = S[i];
		int zp = S[i+1];

		if (c_[2 * z] != -1) {
			if (U_[c_[2 * z]] != -2) {
				//printf("find1:%d U: %d\n",i,U[c[2*z]]);
				k = Find(c_[2 * z], U_);
				if (k != z) {
					Union(k, z, U_);
					p_[k] = z;
					c_[2 * z] = k;
					if (c_[2 * z + 1] == k)
											c_[2 * z + 1] = -1;
				}
			}
		}
		if (c_[2 * z + 1] != -1) {
			if (U_[c_[2 * z + 1]] != -2) {
				//printf("find2:%d\n",i);
				k = Find(c_[2 * z + 1], U_);
				if (k != z) {
					Union(k, z, U_);
					p_[k] = z;
					if (c_[2 * z ] != k)
						c_[2 * z + 1] = k;
					else
						c_[2 * z + 1] = -1;
				}
			}
		}


		if ((vp[z] == vp[zp]) && (off[z] == 0)
				&& (off[zp] == 0) && (own[i] != own[i + 1])) {
			//printf("hello\n\n");
			p_[z] = zp;
			if (c_[2 * zp] == -1)
				c_[2 * zp] = z;
			else
				c_[2 * zp + 1] = z;
			//printf("dup:%d\n",i);
			Union(z, zp, U_);
		}

	}
free(U_);
	}


}
free(own);

	gettimeofday(&tv2, NULL);

	printf("\n stitch time = %f miliseconds\n",			(double) (tv2.tv_usec - tv1.tv_usec) / 1000					+ (double) (tv2.tv_sec - tv1.tv_sec) * 1000);
	//printf("\nvn1:%d,vn2:%d",vn1,vn2);
	//printf("\n final:%d\n",final);

	/*for( i = 0;i < n1+n2;i ++) {
	 if(p[i] != -1 && !lessVertex(i, p[i])) {
	 printf("\npartial damn\n");
	 }
	 }*/

	gettimeofday(&tv1, NULL);
	//printf("before cleanup--- n1:%d,n2:%d,b1:%d,b2:%d\n",n1,n2,b1,b2);



	if(tree==1) cleanup_trees();
	gettimeofday(&tv2, NULL);

	printf("\n cleanup time = %f miliseconds\n",			(double) (tv2.tv_usec - tv1.tv_usec) / 1000					+ (double) (tv2.tv_sec - tv1.tv_sec) * 1000);
	int edge;

	/*for (i=0;i<critical;i++){
	 if (j_n[i]==-1) printf("\njoin root:%d,value:%f, find : %d \n",i,f_v[i],Find(1580910,U));
	 if (s_n[i]==-1) printf("\nsplit root:%d,value:%f,find :%d\n",i,f_v[i],Find(9,U_));



	 }*/



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
	printf("\n writing time = %f miliseconds\n",			(double) (tv2.tv_usec - tv1.tv_usec) / 1000					+ (double) (tv2.tv_sec - tv1.tv_sec) * 1000);

//printf("\ntobecalculated: %d, string %s\n",tree,argv[4]);
	if (tree != 1)
		return 0;

	int *upper_degree = (int*) malloc(num_critical * sizeof(int));
	int *lower_degree = (int*) malloc(num_critical * sizeof(int));
	int *leaves = (int*) malloc(num_critical * sizeof(int));
	int *is_processed = (int*) malloc(num_critical * sizeof(int));

	int b, num_leaves = 0;
	CTGraph = createGraph(num_critical);
gettimeofday(&tv1, NULL);
	#pragma omp parallel for schedule(dynamic) private(b)
	for (b = 0; b < num_critical; b++) {
		upper_degree[b] = 0;
		lower_degree[b] = 0;
		is_processed[b] = 0;
	}
	#pragma omp parallel for schedule(dynamic) private(b)
	for (b = 0; b < num_critical; b++) {
		if (j_n[b] > -1)
		__sync_fetch_and_add(lower_degree+j_n[b],1);
			//lower_degree[j_n[b]]++;
		if (s_n[b] > -1)
		__sync_fetch_and_add(upper_degree+s_n[b],1);
			//upper_degree[s_n[b]]++;
		if (((j_c[2 * b] == -1) && (j_c[2 * b + 1] == -1))
				|| ((s_c[2 * b] == -1) && (s_c[2 * b + 1] == -1)))
		{

			leaves[num_leaves] = b;
			__sync_fetch_and_add(&num_leaves,1);
		}
	}
	//int lcnt, ucnt, u, l;

	#pragma omp parallel for schedule(dynamic) private(i)
	for (i = 0; i < (num_leaves); i++) {
		int lcnt,ucnt,j,u,l;
		int tid = omp_get_thread_num();
		//int i = get_global_id(0);
		j = leaves[i]; //ith leaf
		//printf("%d\n ",i);
		while (!is_processed[j]) { //TODO check

			//is_processed[j]++;
			__sync_fetch_and_add(is_processed+j,1);
			//#pragma omp atomic
			//(touched_critical_points++);
			if (lower_degree[j] == 0 && upper_degree[j] == 1) {
				//upper degree in join Y tree and lower degree in split Y' tree
				//local maxima
				u = j_n[j];
				while (is_processed[u] && j_n[u] != -1)
					u = j_n[u];
				#pragma omp critical
				 addEdge(CTGraph,u,j);
				if (is_processed[u])
					break;
				//lcnt = --lower_degree[u];
				lcnt = __sync_sub_and_fetch (lower_degree+u,1);//remove the vertex
				/*#pragma omp critical
				 lcnt = --lower_degree[u];*/

				ucnt = upper_degree[u];

				if ((lcnt == 0 && ucnt == 1) || (lcnt == 1 && ucnt == 0)) {
					j = u;
				} else {
					break;
				}
			} else if (lower_degree[j] == 1 && upper_degree[j] == 0) {
				l = s_n[j];
				while (is_processed[l] && s_n[l] != -1)
					l = s_n[l];
				#pragma omp critical
				 addEdge(CTGraph,l,j);
				if (is_processed[l])
					break;
				lcnt = lower_degree[l];

				//ucnt = --upper_degree[l];
				ucnt = __sync_sub_and_fetch (upper_degree+l,1);
				//add (u,j)

				if ((lcnt == 0 && ucnt == 1) || (lcnt == 1 && ucnt == 0)) {
					j = l;
				} else {
					break;
				}
			} else {
				break;
			}
		}
	}
	printf("\nct computed\n");
gettimeofday(&tv2, NULL);

	printf("\n tree merge time = %f miliseconds\n",			(double) (tv2.tv_usec - tv1.tv_usec) / 1000					+ (double) (tv2.tv_sec - tv1.tv_sec) * 1000);
	
	
	//-------------------------------------uncomment to write final contour tree to file-----------------
	printf("\nwriting final contour tree\n");
	printGraph(CTGraph,f_v);
//--------------------------------------------------------------------------------
	return 0;
}


