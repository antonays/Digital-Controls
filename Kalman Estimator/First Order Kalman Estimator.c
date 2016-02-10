// First Order Kalman Estimator - No Pointers
// Given a set of input and output samples, we are to identify the transfer function that generates that data.
// We analyzed the step response of the samples and identified the function to be a 1st order function.

// made with no pointers to allow functionality on a NI MyRio Device through LabView - TODO: squeeze everything in one function 
// -- Rhaeg --

#include <stdio.h>
#include <stdlib.h>

#define N 20 // Number of samples

double ud_rid [N];
double yd_rid [N];

static int n = 1 ; // update the order of the system
static int m = 2; // update this for different order m = 2*n;

void Estimate(double lambda, double gamma, double t_result[m]);

// Linear Algebra Functions
void MatrixByMatrix(double a[m][m],double b[m][m], double result[m][m]);
void vectorByMatrix(double matrix[m][m], double vec[m], double result[m]);
void MatrixByVector(double vector[m], double matrix[m][m], double result[m]);
void MatrixByScalar(double scalar, double mat[m][m], double result[m][m]);
void vectorByScalar( double v[m], double scalar, double result[m]);
void VectorByVectorCR(double x[m], double y[m], double result[m][m]);
double VectorByVectorRC(double v1[m], double v2[m]);
void MatrixMatrixAddition(double mat1[m][m], double mat2[m][m], double result[m][m]);
void vectorVectorAddition(double v1[m], double v2[m], double result[m]);
void initiateIdentityMatrix(double result[m][m]);

// Utility Functions
void initializeData();
void randomize(int randomFactor, double result[m]);
void getPhi(int index, int size, double result[m]);
void printArray(double arr[], int size);
void printMatrix(double mat[m][m]);

void main(){
	// lambda is defined as Correction factor
	// gamma is the initializer of the P Matrix
	double lambda = 0.95, gamma = 100;
	double result[m];
	
	// initialize vector yd_rid (output samples) and ud_rid (input samples) [ Hard Coded ]
	initializeData();
	
	Estimate(lambda,gamma, result);
	
	printf("\n Theta: a0, b0 \n");
	printArray(result, m);
}

void Estimate(double lambda, double gamma, double t_result[m]){
	int j=0, i,k;
	double y_=0, y=0;
	double identity[m][m]; 
	double p[m][m], p_[m][m]; // p = Pi    ;    p_ = Pi+1
	// theta contains the coefficients of the resultant Transfer Function -- updated on each iteration
	// phi contains the given samples -- updated on each iteration
	double phi[m], t[m], t_[m];  //   t = theta     ;    t_ = theta (i+1)
	double secondTermNumerator[m][m], Q[m] , R[m][m], L[m], secondTermDenominator, C[m], SecondTerm [m][m], A[m][m], B[m];
	
	// create identity matrix m x m
	initiateIdentityMatrix(identity);
	
	// initial values theta = rand(max=100)  **pg. 43 Greco
	randomize(m, t);

	// initial values for p = gamma * eye(2)   ** pg. 43 Greco
	MatrixByScalar(gamma,identity, p);
		
	for (i=1;i<N;i++){
		
		// get samples
		// phi = [ -yd_rid [ i - x  ]    ud_rid[ i - x ] ] ** x depends on the order
		getPhi(i,n,phi);		
		
		// I P_ = 1/lambda^2 * (Pi - (Pi*phi*phi'*Pi)/(lambda^2+phi'*Pi*phi))
			/* Note: i call the division of matrices inside the parentheses in the function of Pi+1 : SecondTerm which is made of (SeconTermNumerator) / (SeconTermDenominator)   ** pg. 35  Greco  */
			
			// Prepare Terms Of SeconTermNumerator 
			// ( Pi * phi * phi' * Pi ) = ( Q * phi' * Pi) = ( R * Pi ) = secondTermNumerator
			// Q = Pi * phi
			MatrixByVector(phi, p, Q);
			// R = Pi * phi * phi'
			VectorByVectorCR( Q, phi, R);
			// secondTermNumerator = R * Pi;
			MatrixByMatrix( R, p, secondTermNumerator);

			//Prepare Terms Of SeconTermDenominator
			// ( lambda^2+phi'*Pi*phi ) =  (lambda^2 + L * phi) = SeconTermDenominator
			// L = phi' * Pi
			vectorByMatrix( p, phi , L);
			secondTermDenominator = VectorByVectorRC(L ,phi) + (lambda*lambda); // vectorByVectorRC returns scalar (row by column);
			
			//Prepare SecondTerm
			MatrixByScalar((-1)/secondTermDenominator, secondTermNumerator, SecondTerm);
			
			// prepare Subtraction inside Parentheses A = (Pi - SecondTerm);
			MatrixMatrixAddition( p, SecondTerm, A);
			
			// prepare final p_
			MatrixByScalar( 1/(lambda*lambda), A, p_);
		
		
		//II y(i+1) = phi'*t - approximate value for the next sample
		y_ = VectorByVectorRC(phi, t); 
		//printf("\n y_aproximate(i+1): %f\n", y_);
		
		
		//III t_ = t + p_*phi*(yd_rid(i)-y_)
			// Prepare B = p_ * phi;
			MatrixByVector(phi,p_,B);
			// Prepare C  = B * (yd_rid(i)-y_)
			vectorByScalar(B, (yd_rid[i]-y_), C);
			// Prepare final t_;
			vectorVectorAddition(t, C, t_);
		
		
		// swap vars  p = p_  ;    t = t_
		for (k=0;k<m;k++){
			for (j=0;j<m;j++){
				p[k][j]=p_[k][j];
			}
		}
		for (k=0;k<m;k++){
				t[k]=t_[k];
		}
	}
	// copy result to output
	for (k=0;k<m;k++){
		t_result[k] = t_[k];
	}
}

void getPhi(int index, int size, double result[m]){
	int i,j=0;	
	for (i=1;i<size+1;i++){
		result[j] = (-1)*yd_rid[index-i];
		j++;
	}
	for (i=1;i<size+1;i++){
		result[j] = ud_rid[index-i];
		j++;
	}
}

void MatrixByScalar(double scalar, double mat[m][m], double result[m][m]){
	int i,j;
	for (i=0;i<m;i++){
		for (j=0;j<m;j++){
			result[i][j]=mat[i][j] * scalar;
		}
	}
}

void printMatrix(double mat[m][m]){
	int i,j;
	for (i=0;i<m;i++){
		for (j=0;j<m;j++){
			printf("[%.6f],",mat[i][j]);
		}
		
	}
}

void printArray(double arr[], int size){
	int i=0;
	for(i=0;i<size;i++){
		printf(" [%.6f] ", arr[i]);
	}
	printf("\n");
}

void initiateIdentityMatrix(double result[m][m]){
	int i=0,j=0;
	for (j=0;j<m;j++){
		for (i=0;i<m;i++){
			result[i][j]=0;
		}
	}
	for (i=0;i<m;i++){
		result[i][i]=1; 	
	}
}

void MatrixByMatrix(double a[m][m],double b[m][m], double result[m][m]) {
   int i,j,k;
   double sum=0;
	for (i = 0; i < m; i++) {
      for (j = 0; j < m; j++) {
        for (k = 0; k < m; k++) {
          sum = sum + a[i][k]*b[k][j];
        }
        result[i][j] = sum;
        sum = 0;
      }
    }
}

void vectorByMatrix(double matrix[m][m], double vec[m], double result[m]){
	int i,j;
	double tmpresult=0;
	for (i=0;i<m;i++){
		tmpresult=0;
		for (j=0;j<m;j++){
			tmpresult += matrix[i][j]*vec[j];
		}
		result[i] = tmpresult;
	}
}

void MatrixByVector(double vector[m], double matrix[m][m], double result[m]){
	int i,j;
	double tmpResult;
	for (i=0;i<m;i++){
		tmpResult=0;
		for (j=0;j<m;j++){
			tmpResult += matrix[i][j]*vector[j];
		}
		result[i] = tmpResult;
	}
}

void VectorByVectorCR(double x[m], double y[m], double result[m][m]){
	int i,j;
	for (i=0;i<m;i++){
		for (j=0;j<m;j++){
			result[i][j] = x[i]*y[j];
		}
	}
}

void MatrixMatrixAddition(double mat1[m][m], double mat2[m][m], double result[m][m]){
	int i,j;
	for (i=0;i<m;i++){
		for (j=0;j<m;j++){
			result[i][j] = mat1[i][j]+mat2[i][j];
		}
	}
}

double VectorByVectorRC(double v1[m], double v2[m]){
	int i;
	double result=0;
	for (i=0;i<m;i++){	
		result += v1[i]*v2[i];
	}
	return result;
}

void randomize(int randomFactor, double result[]){
	int i=0;
	for (i=0;i<m;i++){
		result[i] = rand() % randomFactor;
	}
}

void vectorByScalar( double v[m], double scalar, double result[m]){
	int i=0;
	double v2[m];
	for (i=0;i<m;i++){
		result[i] = v[i]*scalar;
	}
}

void vectorVectorAddition(double v1[m], double v2[m], double result[m]){
	int i=0;
	double v3[m];
	for (i=0;i<m;i++){
		result[i] = v1[i] + v2[i];
	}
}

void initializeData(){
	ud_rid[0]=-2.00683593750000;
	ud_rid[1]=-2.00683593750000;
	ud_rid[2]=1.99707031250000;
	ud_rid[3]=1.99707031250000;
	ud_rid[4]=1.99707031250000;
	ud_rid[5]=1.99707031250000; 
	ud_rid[6]=1.99707031250000;  
	ud_rid[7]=2.00195312500000; 
	ud_rid[8]=1.99707031250000; 
	ud_rid[9]=1.99707031250000;
	ud_rid[10]=2.00195312500000; 
	ud_rid[11]=2.00195312500000;
	ud_rid[12]=-2.00195312500000;
	ud_rid[13]=-2.00683593750000;
	ud_rid[14]=-2.00683593750000;
	ud_rid[15]=-2.00683593750000;
	ud_rid[16]=-2.00195312500000;  
	ud_rid[17]=-2.00683593750000; 
	ud_rid[18]=-2.00195312500000; 
	ud_rid[19]=-2.00683593750000;

	yd_rid[0]=1.98730468750000;
	yd_rid[1]=1.99218750000000;
	yd_rid[2]= 1.25000000000000;
	yd_rid[3]=-0.72753906250000;
	yd_rid[4]=-1.49414062500000;
	yd_rid[5]=-1.79199218750000;
	yd_rid[6]=-1.90917968750000;
	yd_rid[7]=-1.95312500000000;
	yd_rid[8]=-1.97265625000000;
	yd_rid[9]=-1.97753906250000;
	yd_rid[10]=-1.98242187500000;
	yd_rid[11]=-1.98242187500000;
	yd_rid[12]=-1.24023437500000;
	yd_rid[13]=0.73730468750000;
	yd_rid[14]=1.50390625000000;
	yd_rid[15]=1.80175781250000;
	yd_rid[16]=1.91894531250000;
	yd_rid[17]=1.95800781250000;
	yd_rid[18]=1.97753906250000;
	yd_rid[19]=1.98730468750000;
}