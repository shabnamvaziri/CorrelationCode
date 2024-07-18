#include "def.h"

extern int        S, iter, J, K, I, Iperim, * keepx, max_iter, * y;
extern double* x1, theta, theta0, B, * b, * c, * cperim, Q;
extern double     MAX_DOUBLE;
extern double     MAX_DOUBLE2;
extern double* xi, * p, * d, * alpha, * alpha1, * delta, * epsilon, * beta, * gamma1, * Delta1, obj, * tau, ** newbeta;
extern double* vPerim, * v, rho, * f;
extern int* jPerim, jCounter;
extern int* supplier;
extern double		MasterTime, SubTime;
extern clock_t		startTemp, endTemp;

double solve_MasterProblem(void)
{
	int i, j, jj, k, s;
	int index, index1;  // auxiliar indices to fill in the constraint matrix
	double best_upper_bound, best_lower_bound;
	int nodecount;     //Variables to call cplex
	CPXLPptr  lp;      // data strucutre to store a problem in cplex ...................
	CPXENVptr env;     // cplex environment.............................................
	int       numcols; // number of variables ..........................................
	int       numrows; // number of constraints.........................................
	int       numnz;   // number of non-zero elements in the matrix ....................
	int       objsen;  // optimization sense (min:1, max:-1 ) ..........................
	double* obj;    // objective function coefficients ..............................
	double* rhs;    // right and side of constraints ................................
	char* sense;  // constraints sense (<=: 'L', =:'E', >=:'G') ...................
	int* matbeg; // index of first non-zero element in each row...................
	int* matind; // associated column of each non-zelo element ...................
	double* matval; // coefficient values fo the non-zero elements of constraints....
	double* lb;     // lower bounds of variables.....................................
	double* ub;     // upper bounds of variables.....................................
	int       status;  // optimization status......................... .................
	double* x;      // solution vector (double, even if the problem is integer) .....
	char      probname[16]; // problem name for cplex .......................................
	char* ctype;  // variable type ('C', 'I', 'B') only if integer.................
	double    value;   // objevtive value of solution ..................................
	double    num_v_var, num_vPerim_var, num_y_var, num_theta_var;
	double    pos_theta;
	double* pos_v, * pos_vPerim;
	int* pos_y;
	double tempV = 0;
	int jcounter = 0;
	int ii = 0;
	double temp1 = 0;
	double temp2 = 0;
	double* temp4, * temp3;
	int* BigI;
	double* sumed;
	sumed = create_double_vector(J);
	BigI = create_int_vector(J * max_iter);
	pos_v = create_double_vector(J * K);
	pos_vPerim = create_double_vector(I * J);
	pos_y = create_int_vector(J);
	temp4 = create_double_vector(J * max_iter);
	temp3 = create_double_vector(J * max_iter);

	printf("master= %lf\n", 1, 1.);

	//Initialize CPLEX environment
	env = CPXopenCPLEX(&status);
	if (env == NULL) {
		char  errmsg[1024];
		printf("Could not open CPLEX. \n");
		CPXgeterrorstring(env, status, errmsg);
		printf("%s", errmsg);
	}

	// Create the problem in CPLEX 
	strcpy(probname, "UFLP");
	lp = CPXcreateprob(env, &status, probname);
	if (env == NULL) {
		char  errmsg[1024];
		printf("Could not create LP. \n");
		CPXgeterrorstring(env, status, errmsg);
		printf("%s", errmsg);
	}
	CPXchgobjsen(env, lp, CPX_MIN);

	printf("master= %lf\n", 1.1, 1.1);
	//Define v_jk variables
	index1 = 0;  // index of columns
	numcols = J * K;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			pos_v[j * K + k] = index1;
			obj[index1] = c[j * K + k] * rho;
			ctype[index1] = 'C';
			lb[index1] = 0;
			ub[index1] = CPX_INFBOUND;
			index1++;
		}
	}

	status = CPXnewcols(env, lp, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_v_var = index1;
	printf("master= %lf\n", 2, 2.2);

	//Define vPerim_ij variables
	index1 = 0;  // index of columns
	numcols = I * J;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (i = 0; i < I; i++) {
		for (j = 0; j < J; j++) {
			pos_vPerim[i * J + j] = index1 + num_v_var;
			obj[index1] = cperim[i * (J + 1) + j] * rho;
			ctype[index1] = 'C';
			lb[index1] = 0;
			ub[index1] = CPX_INFBOUND;
			index1++;
		}
	}
	status = CPXnewcols(env, lp, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_vPerim_var = index1;
	printf("master= %lf\n", 1, 1.3);

	//Define theta variables
	index1 = 0;  // index of columns
	numcols = 1;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");


	pos_theta = index1 + num_v_var + num_vPerim_var;
	obj[index1] = 1 - rho;
	ctype[index1] = 'C';
	lb[index1] = -CPX_INFBOUND;
	ub[index1] = CPX_INFBOUND;
	index1++;

	status = CPXnewcols(env, lp, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_theta_var = index1;


	//Define y variables
	index1 = 0;  // index of columns
	numcols = J;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");


	for (j = 0; j < J; j++) {
		pos_y[j] = index1 + num_v_var + num_vPerim_var + num_theta_var;
		obj[index1] = f[j];
		ctype[index1] = 'B';
		lb[index1] = 0;
		ub[index1] = 1;
		index1++;
	}
	status = CPXnewcols(env, lp, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_y_var = index1;
	printf("master= %lf\n", 2, 2.);

	//Add constraint 1 // demand
	numrows = K;
	numnz = J * K;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (k = 0; k < K; k++) {
		sense[index1] = 'G';
		rhs[index1] = d[k];
		matbeg[index1++] = index;
		for (j = 0; j < J; j++) {
			matind[index] = pos_v[j * K + k];
			matval[index++] = 1;
		}
	}

	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//Add constraint 2 // load balance
	numrows = J;
	numnz = K * J + I * J;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (j = 0; j < J; j++) {
		sense[index1] = 'G';
		rhs[index1] = 0;
		matbeg[index1++] = index;
		for (k = 0; k < K; k++) {
			matind[index] = pos_v[j * K + k];
			matval[index++] = -1;
		}
		for (i = 0; i < I; i++) {
			matind[index] = pos_vPerim[i * J + j];
			matval[index++] = 1;
		}
	}


	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//Add constraint 3 \\ capacity of i
	numrows = I;
	numnz = J * I;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;

	for (i = 0; i < I; i++) {
		sense[index1] = 'L';
		rhs[index1] = alpha[i];
		matbeg[index1++] = index;
		for (j = 0; j < J; j++) {
			matind[index] = pos_vPerim[i * J + j];
			matval[index++] = 1;
		}
	}
	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);


	//Add constraint 4 \\ capacity of warehouse j 
	numrows = J;
	numnz = J * I + J;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (j = 0; j < J; j++) {
		sense[index1] = 'L';
		rhs[index1] = 0;
		matbeg[index1++] = index;

		matind[index] = pos_y[j];
		matval[index++] = -1 * alpha1[j];
		for (i = 0; i < I; i++) {
			matind[index] = pos_vPerim[i * J + j];
			matval[index++] = 1;
		}
	}

	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);


	//Add constraint theta >= theta 0
	numrows = 1;
	numnz = 1;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	sense[index1] = 'G';
	rhs[index1] = theta0;
	matbeg[index1++] = index;
	matind[index] = pos_theta;
	matval[index++] = 1;

	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//Add constraint SVI2
	numrows = 1;
	numnz = 1 + J * K + I * J;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	sense[index1] = 'G';
	rhs[index1] = 0;
	matbeg[index1++] = index;
	matind[index] = pos_theta;
	matval[index++] = 1;
	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			matind[index] = pos_v[j * K + k];
			matval[index++] = -1 * c[j * K + k];
		}
	}
	for (i = 0; i < I; i++) {
		for (j = 0; j < J; j++) {
			matind[index] = pos_vPerim[i * J + j];
			matval[index++] = -1 * cperim[i * (J + 1) + j];
		}
	}

	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);


	if (iter > 0) {

		for (ii = 0; ii < iter; ii++) {
			temp1 = 0;
			temp2 = 0;

			for (j = 0; j < J; j++) {
				temp3[j * max_iter + ii] = 0;
				temp4[j * max_iter + ii] = 0;
				BigI[j * max_iter + ii] = 0;
				sumed[j] = 0;
			}
			int countj = 0;
			for (j = 0; j < J; j++) {
				for (jj = 0; jj < J + 1; jj++) {
					for (s = 0; s < S; s++) {
						sumed[j] = sumed[j] + newbeta[(j * (J + 1) + jj) + (s * (J + 1) * (J + 1))][ii];
					}
				}
				if (sumed[j] > 0) {
					BigI[j * max_iter + ii] = 1;
					countj = 1;
				}
				else
					BigI[j * max_iter + ii] = 0;
			}

			for (k = 0; k < K; k++) {
				for (s = 0; s < S; s++) {
					temp1 += beta[(k * S + s) + (ii * S * K)] * p[s] * d[k];
				}
			}

			for (i = 0; i < I; i++) {
				for (s = 0; s < S; s++) {
					temp2 += delta[(i * S + s) + (ii * I * S)] * p[s] * alpha[i];
				}
			}
		
			for (j = 0; j < J; j++) {
				for (jj = 0; jj < J + 1; jj++) {
					for (s = 0; s < S; s++) {
						temp4[j * max_iter + ii] += p[s] * newbeta[(j * (J + 1) + jj) + (s * (J + 1) * (J + 1))][ii] * alpha1[j];
					}
				}
			}
			//Add constraint optimality cut  
			numrows = iter;
			numnz = J + 1;
			d_vector(&rhs, numrows, "open_cplex:2");
			c_vector(&sense, numrows, "open_cplex:3");
			i_vector(&matbeg, numrows, "open_cplex:4");
			i_vector(&matind, numnz, "open_cplex:6");
			d_vector(&matval, numnz, "open_cplex:7");

			index = 0;
			index1 = 0;
			sense[index1] = 'G';
			rhs[index1] = (temp1 - temp2);
			matbeg[index1++] = index;
			matind[index] = pos_theta;
			matval[index++] = 1;
			for (j = 0; j < J; j++) {
				matind[index] = pos_y[j];
				matval[index++] = temp4[j * max_iter + ii];
			}
			status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
			if (status)
				fprintf(stderr, "CPXaddrows failed.\n");
			free(matbeg);
			free(matind);
			free(matval);
			free(sense);
			free(rhs);

			if (countj == 1) {
				/////// add SVI1///////
				numrows = iter;
				numnz = J + 1;
				d_vector(&rhs, numrows, "open_cplex:2");
				c_vector(&sense, numrows, "open_cplex:3");
				i_vector(&matbeg, numrows, "open_cplex:4");
				i_vector(&matind, numnz, "open_cplex:6");
				d_vector(&matval, numnz, "open_cplex:7");

				index = 0;
				index1 = 0;
				sense[index1] = 'G';
				rhs[index1] = 1;
				matbeg[index1++] = index;

				for (j = 0; j < J; j++) {
					matind[index] = pos_y[j];
					matval[index++] = BigI[j * max_iter + ii];
				}
				status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
				if (status)
					fprintf(stderr, "CPXaddrows failed.\n");
				free(matbeg);
				free(matind);
				free(matval);
				free(sense);
				free(rhs);
			}
		}
	}

	//CPXwriteprob(env, lp, "model2.lp", NULL);                          //write the model in .lp format if needed (to debug)

	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); //output display
	CPXsetintparam(env, CPX_PARAM_THREADS, 4);
	//CPXsetintparam(env,CPX_PARAM_INTSOLLIM,1);    //stops after finding first integer sol.
	CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 3); //different levels of output display
	//CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,1);//0:balanced; 1:feasibility; 2:optimality,3:bestbound, 4:hiddenfeas
	CPXsetdblparam(env, CPX_PARAM_TILIM, 3600); // time limit
	//CPXsetdblparam(env,CPX_PARAM_TRELIM, 14000); // B&B memory limit
	status = CPXsetintparam(env, CPX_PARAM_MEMORYEMPHASIS, 1);	//conserve memory where possible
	CPXsetintparam(env, CPX_PARAM_NODEFILEIND, 0);
	CPXsetdblparam(env, CPX_PARAM_EPGAP, 0.000001); // e-optimal solution (%gap)
	//CPXsetdblparam(env,CPX_PARAM_EPAGAP, 0.0000000001); // e-optimal solution (absolute value)
	//CPXsetdblparam(env,CPX_PARAM_EPINT, 0.0000000001); // integer precision
	//CPXsetintparam(env,CPX_PARAM_THREADS, 1); // Number of threads to use
	//CPXsetdblparam(env,CPX_PARAM_EPRHS, 0.0000001);
	//CPXsetintparam(env,CPX_PARAM_REDUCE, 0);  // only needed when adding lazy constraints
	//CPXsetintparam(env,CPX_PARAM_HEURFREQ, -1); //heuristic frequency and intensisty 
	//CPXsetdblparam(env, CPX_PARAM_CUTSFACTOR, 1.0);  //limit the number of cuts added by cplex 1.0002
	//CPXsetdblparam(env,CPX_PARAM_CUTUP,UpperBound+.01); // provide an initial upper bound
	//CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,CPX_MIPEMPHASIS_OPTIMALITY);  // MIP emphasis: optimality, feasibility, moving best bound
	//CPXsetintparam(env,CPX_PARAM_PARALLELMODE, 1); 
	//CPXsetintparam(env,CPX_PARAM_PREIND,0);
	//CPXsetintparam(env,CPX_PARAM_MIPORDIND,CPX_ON); // Turn on or off the use of priorities on bracnhing variables
	//CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,CPX_MIPEMPHASIS_BESTBOUND);  // MIP emphasis: optimality, feasibility, moving best bound

	gettimeofday(&start, NULL);
	//startTemp = clock();
	CPXmipopt(env, lp);  //solve the integer program
	gettimeofday(&stop, NULL);
	//endTemp = clock();
	//MasterTime = (double)(endTemp - startTemp) / (double)(CLOCKS_PER_SEC);
	MasterTime = ((double)(stop.tv_sec - start.tv_sec) * 1000 + (double)(stop.tv_usec - start.tv_usec) / 1000) / 1000;

	printf("master= %lf\n", 5, 5.);

	i = CPXgetstat(env, lp);
	if (i == 101)
		printf("Optimal solution found\n");
	else if (i == 102)
		printf("e-optimal solution found\n");
	else if (i == 103)
		printf(" infeasible solution\n");
	else if (i == 107)
		printf("Time limit reached\n");
	else
		printf("Unknown stopping criterion (%d)\n", i);

	// retrive solution values
	CPXgetmipobjval(env, lp, &value);
	printf("Upper bound: %.2f   ", value);
	best_upper_bound = value;
	//best_upper_bound += 2 * x1F + x2F;
	// If CPLEX was able to find the optimal solution, the previous function provides the optimal solution value
	//if not, it provides the best upper bound
	CPXgetbestobjval(env, lp, &value);  //best lower bound in case the problem was not solved to optimality
	best_lower_bound = value;
	printf("Lower bound: %.2f  \n", value);

	nodecount = CPXgetnodecnt(env, lp);
	printf("Number of BB nodes : %ld  \n", nodecount);

	numcols = CPXgetnumcols(env, lp);
	d_vector(&x, numcols, "open_cplex:0");
	CPXgetmipx(env, lp, x, 0, numcols - 1);  // obtain the values of the decision variables

	if (lp != NULL) {
		status = CPXfreeprob(env, &lp);
		if (status) {
			fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
		}
	}
	if (env != NULL) {
		status = CPXcloseCPLEX(&env);
		if (status) {
			char  errmsg[1024];
			fprintf(stderr, "Could not close CPLEX environment.\n");
			CPXgeterrorstring(env, status, errmsg);
			fprintf(stderr, "%s", errmsg);
		}
	}

	tempV = 0;
	index = 0;
	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			v[j * K + k] = x[index];
			index++;
		}
	}

	for (i = 0; i < I; i++) {
		for (j = 0; j < J; j++) {
			vPerim[i * J + j] = x[index];
			index++;
		}
	}
	theta = x[index];
	printf("theta = %lf\n", theta);
	index++;
	for (j = 0; j < J; j++) {
		if (x[index] > 0.5) {
			y[j] = 1;
			//printf("y[%d]=%d", j, y[j]);
		}
		else
			y[j] = 0;
		index++;
	}
	printf("master= %lf\n", 6, 6.);
	free(x);
	return  best_upper_bound;
}


/* This simple routine frees up the pointer *ptr, and sets *ptr
to NULL */

static void free_and_null(char** ptr)
{
	if (*ptr != NULL) {
		free(*ptr);
		*ptr = NULL;
	}
} /* END free_and_null */

