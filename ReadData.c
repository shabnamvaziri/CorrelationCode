#include "def.h"

extern int        S, iter, J, K, I, Iperim, * keepx, max_iter, * y;
extern double* x1, theta, theta0, B, * b, * c, * cperim, Q, *correlation;
extern double     MAX_DOUBLE;
extern double     MAX_DOUBLE2;
extern double* xi, * p, * d, * alpha, * alpha1, * delta, * epsilon, * beta, * gamma1, * Delta1, obj, *tau, **newbeta;
extern double* vPerim, * v, rho, * f;
extern int* jPerim, jCounter, * selectedfacility;
extern int* supplier;
extern double		MasterTime, SubTime;
extern clock_t		startTemp, endTemp;


// Function used to read the input parameters of the instance being solved

void read_instance(const char* name)
{
	int i, j, k, s, jj;
	FILE* in;

	in = Open_File(name, "r");

	if (fscanf(in, "%d %d %d", &J, &K, &I) != 3) {
		fprintf(stderr, "ERROR: Cannot read instance size \n");
		exit(1);
	}
	//printf("Instance size: N:%d M:%d \n", N, M);
	printf("read");
	Initialize_memory();
	rho = 0.5;
	for (k = 0; k < K; k++) {
		if (fscanf(in, "%lf", &d[k]) != 1) {
			fprintf(stderr, "ERROR: Cannot read capacities and set-up cots \n");
			exit(1);
		}
	}
	printf("read1-1");
	for (j = 0; j < J; j++) {
		if (fscanf(in, "%lf", &f[j]) != 1) {
			fprintf(stderr, "ERROR: Cannot read capacities and set-up cots \n");
			exit(1);
		}
	}
	for (i = 0; i < I; i++) {
		for (j = 0; j < J; j++) {
			fscanf(in, "%lf", &cperim[i * (J + 1) + j]);
		}
	}
	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			fscanf(in, "%lf", &c[j * K + k]);
		}
	}
	printf("read-2");
	for (i = 0; i < I; i++) {
		fscanf(in, "%lf", &alpha[i]);
	}
	for (j = 0; j < J; j++) {
		fscanf(in, "%lf", &alpha1[j]);
	}
	for (j = 0; j < J; j++) {
		b[j] = 1;
	}

	fscanf(in, "%lf", &B);

	for (j = 0; j < J; j++) {
		for (jj = 0; jj < J; jj++) {
			fscanf(in, "%lf", &correlation[j * (J + 1) + jj]);
		}
	}
	
	fclose(in);
	
}

//Function used to allocate memory to arrays based on the size of the input parameter of the specific instance being solved

void Initialize_memory(void)
{
	int s_d;
	s_d = pow(2, J);
	max_iter = 200;
	printf("memory");
	d = create_double_vector(K);
	cperim = create_double_vector(I * (J + 1));
	c = create_double_vector((J + 1) * K);
	alpha = create_double_vector(I);
	alpha1 = create_double_vector((J + 1));
	b = create_double_vector((J + 1));
	//jPerim = create_int_vector(10000);
	keepx = create_int_vector((J + 1));
	y = create_int_vector((J + 1));
	f = create_double_vector(J);
	vPerim = create_double_vector(I * J);
	v = create_double_vector(J * K);
	//S = create_int_vector(iteer);
	supplier = create_int_vector(I);
	p = create_double_vector(s_d);
	beta = create_double_vector(K * s_d * max_iter);
	delta = create_double_vector(I * s_d * max_iter);
	epsilon = create_double_vector((J + 1) * s_d * max_iter);
	gamma1 = create_double_vector((J + 1) * s_d * max_iter);
	Delta1 = create_double_vector((J + 1) * s_d * max_iter);
	xi = create_double_vector((J + 1) * s_d);
	//openFacility = create_int_matrix(100, max_iter);
	selectedfacility = create_int_vector(J);
	correlation = create_double_vector((J + 1) * (J + 1));
	tau = create_double_vector((J + 1) * (J + 1) * s_d);
	newbeta = create_double_matrix((J + 1) * (J + 1) * s_d, max_iter);
}


//Function used to release memory of the arrays based on the size of the input parameter of the specific instance being solved


void free_memory(void)
{
	int i;

	//for(i=0;i<N;i++)
	//  free(c[i]);

	//free(c);
	//free(b);
	//free(f);
	//free(d);
	//free(open_facility);
	//free(customer_assign);

}


// Function used to open data file

FILE* Open_File(const char* name, const char* mode)  // This function opens the file "name" under mode "mode" (reading, writing, etc)
{
	FILE* file;

	if ((file = fopen(name, mode)) == NULL) {
		printf("\nError: File cannot be opened \n");
		//OK=1;
		exit(8);
	}
	return file;
}

// Functions to allocate memory to one and two dimensional arrays


int** create_int_matrix(int rows, int Columns)
{
	int i;
	int** ptr;

	if ((ptr = (int**)calloc(rows, sizeof(int*))) == NULL) {
		printf("\nError: Insuficient memory \n");
		exit(8);
	}
	for (i = 0; i < rows; i++)
		ptr[i] = create_int_vector(Columns);
	return ptr;
}

double** create_double_matrix(int rows, int Columns)
{
	int i;
	double** ptr;

	if ((ptr = (double**)calloc(rows, sizeof(double*))) == NULL) {
		printf("\nError: Insuficient memory \n");
		exit(8);
	}
	for (i = 0; i < rows; i++) {
		ptr[i] = create_double_vector(Columns);
	}
	return ptr;
}

double*** create_double_matrix3D(int rows, int Columns, int Columns2)
{
	int i;
	double*** ptr;

	if ((ptr = (double***)calloc(rows, sizeof(double*))) == NULL) {
		printf("\nError: Insuficient memory \n");
		exit(8);
	}
	for (i = 0; i < rows; i++) {
		ptr[i] = create_double_matrix(Columns, Columns2);
	}
	return ptr;
}




int* create_int_vector(int dim)
{
	int* ptr;

	if ((ptr = (int*)calloc(dim, sizeof(int))) == NULL) {
		printf("\nError: Insuficient memory \n");
		exit(8);
	}
	return ptr;
}

double* create_double_vector(int dim)
{
	double* ptr;

	if ((ptr = (double*)calloc(dim, sizeof(double))) == NULL) {
		printf("\nError: Insuficient memory \n");
		exit(8);
	}
	return ptr;
}



// CPLEX functions to allocate memeory to arrays

void i_vector(int** vector, int n, char* s)
{
	if ((*vector = (int*)calloc(n, sizeof(int))) == NULL)
		//error(s);
		printf("Error: Insuficient memory \n");
	return;
}

void d_vector(double** vector, int n, char* s)
{
	if ((*vector = (double*)calloc(n, sizeof(double))) == NULL)
		// error(s);
		printf("Error: Insuficient memory \n");
	return;
}

void c_vector(char** vector, int n, char* s)
{
	if ((*vector = (char*)calloc(n, sizeof(char))) == NULL)
		//error(s);
		printf("Error: Insuficient memory \n");
	return;
}



