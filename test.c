#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>

// FUNCTIONS FOR SETUP
int numberCells(double, double *);
void setVector(int, double, double *, double *);
void setObst(double *);
void setGoal(double *, double *, double *);
void setInitialValue(double *, double *, double *);
void conditionValue(double *, double *, double *, int, int, int);
void setInitialPolicy(double *, double *, char *);
void conditionPolicy(double *, double *, char *, int, int, int);

// FUNCTIONS FOR VALUE ITERATION
void valueIteration(double *, double *, double *, char *, double *);
void conditionR(int, int, int, double *, double *);
void conditionTheta(int, int, int, double *, double *);
void conditionPhi(int, int, int, double *, double *);
void computeTotalCost(double *, double *);
double computeNewValue(double *);
char computeNewPolicy(double *);

// FUNCTIONS FOR ANALYSIS
double cpuSecond(void);

// DEFINE GLOBAL PARAMETERS
int nr, ntheta, nphi;
double perr;
double gamma1;
double vGoal, vObst, vMove;
double vInitial;
int numActions = 7;

int main(int argc, char **argv){
	double iStart=cpuSecond();

	// DEFINE PARAMETERS
	double dr, dtheta, dphi;
	double rdim[2], thetadim[2], phidim[2];
	double *rVec, *thetaVec, *phiVec;
	// - minimum grid resolution for r, theta, phi
	dr = 5.0, dtheta = 180.0, dphi = 180;
	// - dimensions of the state space
	rdim[0] = 0.0, rdim[1] = 10.0;
	thetadim[0] = 0.0, thetadim[1] = 360.0;
	phidim[0] = 0.0, phidim[1] = 360.0;
	// - number of grid cells
	nr = numberCells(dr, rdim);
	ntheta = numberCells(dtheta, thetadim);
	nphi = numberCells(dphi, phidim);
	// - vectors for r, theta, phi
	rVec = (double *)malloc(sizeof(double)*nr);
	thetaVec = (double *)malloc(sizeof(double)*ntheta);
	phiVec = (double *)malloc(sizeof(double)*nphi);
	setVector(nr, dr, rdim, rVec);
	setVector(ntheta, dtheta, thetadim, thetaVec);
	setVector(nphi, dphi, phidim, phiVec);
	// - probability of going the wrong way
	perr = 0.0;
	// attenuation rate
	gamma1 = 1.0;
	// - value of goal, collision, movement
	vGoal = 100.0;
	vObst = -100.0;
	vMove = -1.0;
	// initial guess at all values
	vInitial = 0.0;

	// DEFINE OBSTACLE AND GOAL LOCATIONS
	double *isobst, *isgoal;
	isobst = (double *)calloc(nr*ntheta*nphi, sizeof(double));
	isgoal = (double *)calloc(nr*ntheta*nphi, sizeof(double));
	setObst(isobst);
	setGoal(thetaVec, phiVec, isgoal);

	// DEFINE INITIAL GUESS AT VALUE AND POLICY
	double *J;
	char *U;
	J = (double *)calloc(nr*ntheta*nphi, sizeof(double));
	U = (char *)calloc(nr*ntheta*nphi, sizeof(char));
	setInitialValue(isobst, isgoal, J);
	setInitialPolicy(isobst, isgoal, U);

	// DO VALUE ITERATION
	int T = 100;
	double *Jprev;
	char *Uprev;
	Jprev = (double *)calloc(nr*ntheta*nphi, sizeof(double));
	Uprev = (char *)calloc(nr*ntheta*nphi, sizeof(char));

	for(int t=0; t<T; t++){
		printf("Iteration %d\n", t+1);

		// Iterate over all states.
		memcpy(Jprev, J, sizeof(double)*nr*ntheta*nphi);
		memcpy(Uprev, U, sizeof(char)*nr*ntheta*nphi);

		valueIteration(isobst, isgoal, J, U, Jprev);

		for(int x=0; x<nr*ntheta*nphi; x++){
			printf("%2d J=%3.1f U=%2d\n", x, J[x], U[x]);
		}
		printf("\n");
	}

	// free used memory
	free(rVec);
	free(thetaVec);
	free(phiVec);
	free(isobst);
	free(isgoal);
	free(J);
	free(U);
	free(Jprev);
	free(Uprev);
	
	double iElaps=cpuSecond()-iStart;
	printf("Time elapsed on GPU = %f ms\n", iElaps*1000.0f);	
	return(0);
}

/*--------------- SETUP FUNCTIONS ----------------*/

int numberCells(double d, double *dim){
	int n = 0;
	double diff;
	diff = dim[1]-dim[0];

	if(d<1 || d>diff){
		printf("value of resolution or dimension is invalid.\n");
	}
	else{
		n = floorf(diff/d+1.0);
	}

	return n;
}

void setVector(int n, double d, double *dim, double *Vec){
	double value;
	value = dim[0];

	for(int i=0; i<n; i++){
		Vec[i] = value;
		value += d;
	}
}

void setObst(double *isobst){
	for(int j=0; j<ntheta; j++){
		for(int k=0;k<nphi; k++){
			isobst[nr*ntheta*k+(nr-1)*ntheta+j] = 1;
		}
	}
}

void setGoal(double *thetaVec, double *phiVec, double *isgoal){
	for(int j=0; j<ntheta; j++){
		for(int k=0; k<nphi; k++){
			if(thetaVec[j]==phiVec[k])
				isgoal[nr*ntheta*k+j] = 1;
		}
	}
}

void setInitialValue(double *isobst, double *isgoal, double *J){
	for(int i=0; i<nr; i++){
		for(int j=0; j<ntheta; j++){
			for(int k=0; k<nphi; k++){
				conditionValue(isobst, isgoal, J, i, j, k);
			}
		}
	}
}

void conditionValue(double *isobst, double *isgoal, double *J, int i, int j, int k){
	if(isobst[nr*ntheta*k+ntheta*i+j]){
		J[nr*ntheta*k+ntheta*i+j] = vObst;
	}
	else if(isgoal[nr*ntheta*k+ntheta*i+j]){
		J[nr*ntheta*k+ntheta*i+j] = vGoal;
	}
	else{
		J[nr*ntheta*k+ntheta*i+j] = vInitial;
	}
}

void setInitialPolicy(double *isobst, double *isgoal, char *U){
	srand((unsigned)time(NULL));

	for(int i=0; i<nr; i++){
		for(int j=0; j<ntheta; j++){
			for(int k=0; k<nphi; k++){
				conditionPolicy(isobst, isgoal, U, i, j, k);
			}
		}
	}
}

void conditionPolicy(double *isobst, double *isgoal, char *U, int i, int j, int k){
	if(isobst[nr*ntheta*k+ntheta*i+j]){
		U[nr*ntheta*k+ntheta*i+j] = -1;
	}
	else if(isgoal[nr*ntheta*k+ntheta*i+j]){
		U[nr*ntheta*k+ntheta*i+j] = -1;
	}
	else{
		char r = rand() % 7;
		U[nr*ntheta*k+ntheta*i+j] = r;
	}
}

/*--------------- VALUE ITERATION FUNCTIONS ----------------*/

void valueIteration(double *isobst, double *isgoal, double *J, char *U, double *Jprev){
	double *tempCost, *totalCost;
	tempCost = (double *)calloc(numActions, sizeof(double));
	totalCost = (double *)calloc(numActions, sizeof(double));

	for(int i=0; i<nr; i++){
		for(int j=0; j<ntheta; j++){
			for(int k=0; k<nphi; k++){
				if(!isobst[nr*ntheta*k+ntheta*i+j] && !isgoal[nr*ntheta*k+ntheta*i+j]){
					tempCost[0]=Jprev[nr*ntheta*k+ntheta*i+j];
					// condition of r
					conditionR(i, j, k, tempCost, Jprev);

					// Compute the total expected cost for each of the possible actions.
					computeTotalCost(tempCost, totalCost);

					// Compute the new exptected cost-to-go, by taking the maximum over
					// possible actions.
					J[nr*ntheta*k+ntheta*i+j] = computeNewValue(totalCost);
					U[nr*ntheta*k+ntheta*i+j] = computeNewPolicy(totalCost);

				}
			}
		}
	}

	free(tempCost);
	free(totalCost);
}

void conditionR(int i, int j, int k, double *tempCost, double *Jprev){
	if(i==0){
		tempCost[1] = Jprev[nr*ntheta*k+ntheta*(i+1)+j];
		tempCost[2] = Jprev[nr*ntheta*k+ntheta*i+j];
	}
	else{
		tempCost[1] = Jprev[nr*ntheta*k+ntheta*(i+1)+j];
		tempCost[2] = Jprev[nr*ntheta*k+ntheta*(i-1)+j];
	}
	conditionTheta(i, j, k, tempCost, Jprev);
}

void conditionTheta(int i, int j, int k, double *tempCost, double *Jprev){
	if(j==0){
		tempCost[3] = Jprev[nr*ntheta*k+ntheta*i+(j+1)];
		tempCost[4] = Jprev[nr*ntheta*k+ntheta*i+(ntheta-1)];
	}
	else if(j==ntheta-1){
		tempCost[3] = Jprev[nr*ntheta*k+ntheta*i];
		tempCost[4] = Jprev[nr*ntheta*k+ntheta*i+(j-1)];
	}
	else{
		tempCost[3] = Jprev[nr*ntheta*k+ntheta*i+(j+1)];
		tempCost[4] = Jprev[nr*ntheta*k+ntheta*i+(j-1)];
	}
	conditionPhi(i, j, k, tempCost, Jprev);
}

void conditionPhi(int i, int j, int k, double *tempCost, double *Jprev){
	if(k==0){
		tempCost[5] = Jprev[nr*ntheta*(k+1)+ntheta*i+j];
		tempCost[6] = Jprev[nr*ntheta*(nphi-1)+ntheta*i+j];
	}
	else if(k==nphi-1){
		tempCost[5] = Jprev[ntheta*i+j];
		tempCost[6] = Jprev[nr*ntheta*(k-1)+ntheta*i+j];
	}
	else{
		tempCost[5] = Jprev[nr*ntheta*(k+1)+ntheta*i+j];
		tempCost[6] = Jprev[nr*ntheta*(k-1)+ntheta*i+j];
	}
}

void computeTotalCost(double *tempCost, double *totalCost){
	double tempCostTotal=0;
	
	for(int n=0; n<numActions; n++){
		tempCostTotal+=tempCost[n];
	}
	for(int n=0; n<numActions; n++){
		totalCost[n]=vMove+gamma1*((1-perr)*tempCost[n]+(perr/6)*(tempCostTotal-tempCost[n]));
	}
}

double computeNewValue(double *totalCost){
	double max;
	max = totalCost[0];

	for(int n=0; n<numActions; n++){
		if(totalCost[n]>max)
			max=totalCost[n];
	}

	return max;
}

char computeNewPolicy(double *totalCost){
	double max;
	char idx;
	max = totalCost[0];

	for(int n=0; n<numActions; n++){
		if(totalCost[n]>max){
			max=totalCost[n];
			idx=n;
		}

	}

	return idx;
}

/*------------ FUNCTIONS FOR ANALYSIS --------------*/
double cpuSecond(void){
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return ((double)tp.tv_sec+(double)tp.tv_usec*1.e-6);
}
