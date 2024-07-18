#include "def.h"

extern int        S, iter, J, K, I, Iperim, * keepx, max_iter, * y;
extern double* x1, theta, theta0, B, * b, * c, * cperim, Q, * correlation;
extern double     MAX_DOUBLE;
extern double     MAX_DOUBLE2;
extern double* xi, * p, * d, * alpha, * alpha1, * delta, * epsilon, * beta, * gamma1, * Delta1, obj, * tau, ** newbeta;
extern double* vPerim, * v, rho, * f;
extern int* jPerim, jCounter;
extern int* supplier, statusI;
extern double		MasterTime, SubTime;
extern clock_t		startTemp, endTemp;

double solve_SubProblem(void)
{
	int i, j, jj, k, s;
	int index, index1;  // auxiliar indices to fill in the constraint matrix
	double best_upper_bound, best_lower_bound;
	int nodecount;     //Variables to call cplex
	CPXLPptr  lp2;      // data strucutre to store a problem in cplex ...................
	CPXENVptr env2;     // cplex environment.............................................
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
	char		probname[16]; // problem name for cplex .......................................
	char* ctype;  // variable type ('C', 'I', 'B') only if integer.................
	double    value;   // objevtive value of solution ..................................
	double      num_x_var;
	double* pos_x;
	double		temp1 = 0, temp2 = 0, temp3 = 0, * temp4;
	double      num_beta_var, num_delta_var, num_tau_var, num_gamma_var, num_newbeta_var;
	double* pos_beta, * pos_delta, * pos_tau, * pos_gamma, * M, *pos_newbeta;
	int auxilary = 0;
	int auxilary1 = 0;

	pos_x = create_double_vector(J + 1);
	pos_beta = create_double_vector(K * S);
	pos_delta = create_double_vector(I * S);
	M = create_double_vector(J + 1);
	pos_tau = create_double_vector((J + 1) * (J + 1) * S);
	pos_newbeta = create_double_vector((J + 1) * (J + 1) * S);
	pos_gamma = create_double_vector((J + 1) * S);
	
	auxilary = 0;
	auxilary1 = 0;
	y[J] = 1;
	alpha1[J] = 0;
	for (s = 0; s < S; s++) {
		xi[J * S + s] = 0;
	}
	
	b[J] = 100;
	for (k = 0; k < K; k++) {
		c[J * K + k] = 150000;
		alpha1[J] += d[k];
	}
	for (i = 0; i < I; i++) {
		cperim[i * (J + 1) + J] = 150000;
	}

	for (j = 0; j < J + 1; j++) {
		M[j] = MAX_DOUBLE2;
		correlation[J * (J + 1) + j] = 0;
	}


	//Initialize CPLEX environment
	env2 = CPXopenCPLEX(&status);
	if (env2 == NULL) {
		char  errmsg[1024];
		printf("Could not open CPLEX. \n");
		CPXgeterrorstring(env2, status, errmsg);
		printf("%s", errmsg);
	}

	// Create the problem in CPLEX 
	strcpy(probname, "UFLP");
	lp2 = CPXcreateprob(env2, &status, probname);
	if (env2 == NULL) {
		char  errmsg[1024];
		printf("Could not create LP. \n");
		CPXgeterrorstring(env2, status, errmsg);
		printf("%s", errmsg);
	}

	CPXchgobjsen(env2, lp2, CPX_MAX);

	//Define x variables
	index1 = 0;  // index of columns
	numcols = J + 1;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");


	for (j = 0; j < J + 1; j++) {
		pos_x[j] = index1;
		obj[index1] = 0;
		ctype[index1] = 'B';
		lb[index1] = 0;
		ub[index1] = 1;
		index1++;
	}

	status = CPXnewcols(env2, lp2, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_x_var = index1;

	//Define beta variables
	index1 = 0;  // index of columns
	numcols = K * S;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (k = 0; k < K; k++) {
		for (s = 0; s < S; s++) {
			//sprintf(colname[(int)(index1 + num_x_var)], "beta%3d_%3d", k,s);
			pos_beta[k * S + s] = index1 + num_x_var;
			obj[index1] = p[s] * d[k];
			ctype[index1] = 'C';
			lb[index1] = 0;
			ub[index1] = CPX_INFBOUND;
			index1++;
		}
	}

	status = CPXnewcols(env2, lp2, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_beta_var = index1;

	//Define delta variables
	index1 = 0;  // index of columns
	numcols = I * S;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (i = 0; i < I; i++) {
		for (s = 0; s < S; s++) {
			//sprintf(colname[(int)(index1 + num_x_var + num_beta_var)], "delta%3d_%3d", i, s);
			pos_delta[i * S + s] = index1 + num_x_var + num_beta_var;
			obj[index1] = -1 * (p[s] * alpha[i]);
			ctype[index1] = 'C';
			lb[index1] = 0;
			ub[index1] = CPX_INFBOUND;
			index1++;
		}
	}

	status = CPXnewcols(env2, lp2, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_delta_var = index1;

	//Define tau variables
	index1 = 0;  // index of columns
	numcols = (J + 1) * (J + 1) * S;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (j = 0; j < J + 1; j++) {
		for (jj = 0; jj < J + 1; jj++) {
			for (s = 0; s < S; s++) {
				pos_tau[(j * (J + 1) + jj) + (s * (J + 1) * (J + 1))] = index1 + num_x_var + num_beta_var + num_delta_var;
				obj[index1] = 0;
				ctype[index1] = 'C';
				lb[index1] = 0;
				ub[index1] = CPX_INFBOUND;
				index1++;
			}
		}
	}

	status = CPXnewcols(env2, lp2, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_tau_var = index1;

	//Define gamma variables
	index1 = 0;  // index of columns
	numcols = (J + 1) * S;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (j = 0; j < J + 1; j++) {
		for (s = 0; s < S; s++) {
			pos_gamma[j * S + s] = index1 + num_x_var + num_beta_var + num_delta_var + num_tau_var;
			obj[index1] = 0;
			ctype[index1] = 'C';
			lb[index1] = 0;
			ub[index1] = CPX_INFBOUND;
			index1++;
		}
	}

	status = CPXnewcols(env2, lp2, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_gamma_var = index1;

	//Define newbeta variables
	index1 = 0;  // index of columns
	numcols = (J + 1) * (J+1) * S;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (j = 0; j < J + 1; j++) {
		for (jj = 0; jj < J + 1; jj++) {
			for (s = 0; s < S; s++) {
				pos_newbeta[(j * (J + 1) + jj) + (s * (J + 1) * (J + 1))] = index1 + num_x_var + num_beta_var + num_delta_var + num_tau_var + num_gamma_var;
				obj[index1] = p[s] * alpha1[j] * y[j] * -1;
				ctype[index1] = 'C';
				lb[index1] = 0;
				ub[index1] = CPX_INFBOUND;
				index1++;
			}
		}
	}

	status = CPXnewcols(env2, lp2, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_newbeta_var = index1;


	//Add budget constraint 
	numrows = 1;
	numnz = J + 1;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	sense[index1] = 'L';
	rhs[index1] = B;
	matbeg[index1++] = index;
	for (j = 0; j < J + 1; j++) {
		matind[index] = pos_x[j];
		matval[index++] = b[j];
	}

	status = CPXaddrows(env2, lp2, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//Add constraint beta-gamma
	numrows = (J + 1) * K * S;
	numnz = 2 * (J + 1) * K * S;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (j = 0; j < J + 1; j++) {
		for (k = 0; k < K; k++) {
			for (s = 0; s < S; s++) {
				sense[index1] = 'L';
				rhs[index1] = c[j * K + k];
				matbeg[index1++] = index;
				matind[index] = pos_beta[k * S + s];
				matval[index++] = 1;
				matind[index] = pos_gamma[j * S + s];
				matval[index++] = -1;
			}
		}
	}

	status = CPXaddrows(env2, lp2, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//Add constraint 2 gamma-delta-epsilon  
	numrows = S * (J + 1) * I;
	numnz = 3 * S * (J + 1) * I * (J + 1);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (i = 0; i < I; i++) {
		for (j = 0; j < J + 1; j++) {
			for (s = 0; s < S; s++) {
				sense[index1] = 'L';
				rhs[index1] = cperim[i * (J + 1) + j];
				matbeg[index1++] = index;
				matind[index] = pos_gamma[j * S + s];
				matval[index++] = 1;
				matind[index] = pos_delta[i * S + s];
				matval[index++] = -1;
				for (jj = 0; jj < J + 1; jj++) {
					matind[index] = pos_tau[(j * (J + 1) + jj) + (s * (J + 1) * (J + 1))];
					matval[index++] = -1;
				}
			}
		}
	}

	status = CPXaddrows(env2, lp2, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);



	//Add constraint 35  
	numrows = S * (J + 1) * (J + 1);
	numnz = 2 * S * (J + 1) * (J + 1);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (j = 0; j < J + 1; j++) {
		for (jj = 0; jj < J + 1; jj++) {
			for (s = 0; s < S; s++) {
				sense[index1] = 'G';
				rhs[index1] = 0;
				matbeg[index1++] = index;
				matind[index] = pos_newbeta[(j * (J + 1) + jj) + (s * (J + 1) * (J + 1))];
				matval[index++] = 1;
				matind[index] = pos_tau[(j * (J + 1) + jj) + (s * (J + 1) * (J + 1))];
				matval[index++] = (1 - correlation[j * (J + 1) + jj]) * -1;
			}
		}
	}

	status = CPXaddrows(env2, lp2, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//Add constraint 36  
	numrows = S * (J + 1) * (J + 1);
	numnz = (2 * S * (J + 1) * (J + 1) + (J + 1)) * S * 2;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (j = 0; j < J + 1; j++) {
		for (jj = 0; jj < J + 1; jj++) {
			for (s = 0; s < S; s++) {
				sense[index1] = 'G';
				rhs[index1] = 0;
				matbeg[index1++] = index;
				matind[index] = pos_newbeta[(j * (J + 1) + jj) + (s * (J + 1) * (J + 1))];
				matval[index++] = 1;
				matind[index] = pos_tau[(j * (J + 1) + jj) + (s * (J + 1) * (J + 1))];
				matval[index++] = -1;
				matind[index] = pos_x[jj];
				matval[index++] = M[j] * xi[jj * S + s];
			}
		}
	}

	status = CPXaddrows(env2, lp2, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	///////x<y
	numrows = (J + 1);
	numnz = 2 * (J + 1);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (j = 0; j < J + 1; j++) {
		sense[index1] = 'L';
		rhs[index1] = y[j];
		matbeg[index1++] = index;
		matind[index] = pos_x[j];
		matval[index++] = 1;
	}

	status = CPXaddrows(env2, lp2, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);



//	CPXwriteprob(env2, lp2, "model3.lp", NULL);                          //write the model in .lp format if needed (to debug)
	CPXsetintparam(env2, CPX_PARAM_THREADS, 4);
	CPXsetintparam(env2, CPX_PARAM_SCRIND, CPX_ON); //output display
	//CPXsetintparam(env,CPX_PARAM_INTSOLLIM,1);    //stops after finding first integer sol.
	//CPXsetintparam(env2, CPX_PARAM_MIPDISPLAY, 3); //different levels of output display
	//CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,1);//0:balanced; 1:feasibility; 2:optimality,3:bestbound, 4:hiddenfeas
	CPXsetdblparam(env2, CPX_PARAM_TILIM, 10800); // time limit
	CPXsetdblparam(env2, CPX_PARAM_EPGAP, 0.01); // e-optimal solution (%gap)
	CPXsetintparam(env2, CPX_PARAM_NODEFILEIND, 0);
	status = CPXsetintparam(env2, CPX_PARAM_MEMORYEMPHASIS, 2);	//conserve memory where possible

	//CPXsetdblparam(env,CPX_PARAM_TRELIM, 14000); // B&B memory limit
	//CPXsetdblparam(env,CPX_PARAM_EPGAP, 0.0000000001); // e-optimal solution (%gap)
	//CPXsetdblparam(env,CPX_PARAM_EPAGAP, 0.0000000001); // e-optimal solution (absolute value)
	//CPXsetdblparam(env,CPX_PARAM_EPINT, 0.0000000001); // integer precision
	//CPXsetintparam(env,CPX_PARAM_THREADS, 1); // Number of threads to use
	//CPXsetdblparam(env,CPX_PARAM_EPRHS, 0.0000001);
	//CPXsetintparam(env,CPX_PARAM_REDUCE, 0);  // only needed when adding lazy constraints
	//CPXsetintparam(env,CPX_PARAM_HEURFREQ, -1); //heuristic frequency and intensisty 
	//CPXsetdblparam(env, CPX_PARAM_CUTSFACTOR, 1.0);  //limit the number of cuts added by cplex 1.0002
	//CPXsetdblparam(env,CPX_PARAM_CUTUP,UpperBound+.01); // provide an initial upper bound
	//CPXsetintparam(env2, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_OPTIMALITY);  // MIP emphasis: optimality, feasibility, moving best bound
	//CPXsetintparam(env,CPX_PARAM_PARALLELMODE, 1); 
	//CPXsetintparam(env,CPX_PARAM_PREIND,0);
	//CPXsetintparam(env,CPX_PARAM_MIPORDIND,CPX_ON); // Turn on or off the use of priorities on bracnhing variables
	if (B > 1)
		CPXsetintparam(env2, CPX_PARAM_MIPEMPHASIS, 1);  // MIP emphasis: optimality, feasibility, moving best bound


	gettimeofday(&start, NULL);
	//startTemp = clock();
	CPXmipopt(env2, lp2);  //solve the integer program
	gettimeofday(&stop, NULL);
	//endTemp = clock();
	//SubTime = (double)(endTemp - startTemp) / (double)(CLOCKS_PER_SEC);
	SubTime = ((double)(stop.tv_sec - start.tv_sec) * 1000 + (double)(stop.tv_usec - start.tv_usec) / 1000) / 1000;
	
	i = CPXgetstat(env2, lp2);
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

	statusI = i;

	// retrive solution values
	CPXgetmipobjval(env2, lp2, &value);
	printf("Upper bound: %.2f   ", value);
	best_upper_bound = value;
	// If CPLEX was able to find the optimal solution, the previous function provides the optimal solution value
	//if not, it provides the best upper bound
	CPXgetbestobjval(env2, lp2, &value);  //best lower bound in case the problem was not solved to optimality
	best_lower_bound = value;
	printf("Lower bound: %.2f  \n", value);

	nodecount = CPXgetnodecnt(env2, lp2);
	printf("Number of BB nodes : %ld  \n", nodecount);

	numcols = CPXgetnumcols(env2, lp2);
	d_vector(&x, numcols, "open_cplex:0");
	CPXgetmipx(env2, lp2, x, 0, numcols - 1);  // obtain the values of the decision variables

	if (lp2 != NULL) {
		status = CPXfreeprob(env2, &lp2);
		if (status) {
			fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
		}
	}
	if (env2 != NULL) {
		status = CPXcloseCPLEX(&env2);
		if (status) {
			char  errmsg[1024];
			fprintf(stderr, "Could not close CPLEX environment.\n");
			CPXgeterrorstring(env2, status, errmsg);
			fprintf(stderr, "%s", errmsg);
		}
	}


	index = 0;
	for (j = 0; j < J + 1; j++) {
		if (x[index] > 0.5)
			keepx[j] = 1;
		else
			keepx[j] = 0;
		//printf("x[%d] = %d\n", j, keepx[j]);
		index++;
	}
	for (k = 0; k < K; k++) {
		for (s = 0; s < S; s++) {
			beta[(k * S + s) + (iter * S * K)] = x[index];
			index++;
		}
	}
	for (i = 0; i < I; i++) {
		for (s = 0; s < S; s++) {
			delta[(i * S + s) + (iter * I * S)] = x[index];
			index++;
		}
	}
	double sumepsi = 0;
	for (j = 0; j < J + 1; j++) {
		for (jj = 0; jj < J + 1; jj++) {
			for (s = 0; s < S; s++) {
				index++;
			}
		}
	}
	for (j = 0; j < J + 1; j++) {
		for (s = 0; s < S; s++) {
			gamma1[(j * S + s) + (iter * (J + 1) * S)] = x[index];
			index++;
		}
	}
	for (j = 0; j < J + 1; j++) {
		for (jj = 0; jj < J + 1; jj++) {
			for (s = 0; s < S; s++) {
				newbeta[(j * (J + 1) + jj) + (s * (J + 1) * (J + 1))][iter] = x[index];
				/*if (newbeta[(j * (J + 1) + jj) + (s * (J + 1) * (J + 1))][iter] > 0) {
					printf("hi");
				}*/
				index++;
			}
		}
	}
	printf("best_upper_bound = %f\n", best_upper_bound + 0.0000001 * sumepsi);
	//theta = x[index];
	free(x);


	//return best_upper_bound ;
	return best_upper_bound + 0.0000001 * sumepsi;
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

