// First Order Kalman Estimator - No Pointers
// built for labview function node
// Anton Ayzenberg

#include <stdio.h>
#include <stdlib.h>

#define N 30

float ud_rid [N];
float yd_rid [N];

static int n = 2 ; // update the order
static int m = 4; // update this for different order m = 2*n;

void Estimate(float lambda, float gamma, float t_result[m]);

// Linear Algebra Functions
void MatrixByMatrix(float a[m][m],float b[m][m], float result[m][m]);
void vectorByMatrix(float matrix[m][m], float vec[m], float result[m]);
void MatrixByVector(float vector[m], float matrix[m][m], float result[m]);
void MatrixByScalar(float scalar, float mat[m][m], float result[m][m]);
void vectorByScalar( float v[m], float scalar, float result[m]);
void VectorByVectorCR(float x[m], float y[m], float result[m][m]);
float VectorByVectorRC(float v1[m], float v2[m]);
void MatrixMatrixAddition(float mat1[m][m], float mat2[m][m], float result[m][m]);
void vectorVectorAddition(float v1[m], float v2[m], float result[m]);
void initiateIdentityMatrix(float result[m][m]);

// Utility Functions
void initializeData();
void randomize(int randomFactor, float result[m]);
void getPhi(int index, int size, float result[m]);
void printArray(float arr[], int size);
void printMatrix(float mat[m][m]);

void main(){
	float lambda = 0.95, gamma = 100;
	float t_final[m];

	initializeData();
	
	Estimate(lambda,gamma, t_final);
	
	printf("\n Theta: a0, b0 \n");
	printArray(t_final, m);
}

void Estimate(float lambda, float gamma, float t_result[m]){
	int j=0, i,k,l;
	float y_=0, y=0;
	float p[m][m], p_[m][m];
	float phi[m], t[m], t_[m];
	float secondTermNumerator[m][m], Q[m] , R[m][m], L[m], secondTermDenominator, C[m], SecondTerm [m][m], A[m][m], B[m];
	float tmpResult;
	float result=0;
	
	
	for (j=0;j<m;j++){
		for (i=0;i<m;i++){
			p[i][j]=0;
		}
	}
	for (i=0;i<m;i++){
		p[i][i]=1 * gamma; 	
	}
	
	
	t[0]=1;
	t[1]=1;
	t[2]=1;
	t[3]=1;
		
	for (i=2;i<N;i++){
		
		phi[0] = (-1)*yd_rid[i-1];
		phi[1] = (-1)*yd_rid[i-2];
		phi[2] = ud_rid[i-1];
		phi[3] = ud_rid[i-2];
		
		
		for (k=0;k<m;k++){
			tmpResult =0;
			for (j=0;j<m;j++){
				 tmpResult += p[k][j]*phi[j];
			}
			Q[k] = tmpResult;
		}
		
		for (k=0;k<m;k++){
			for (j=0;j<m;j++){
				R[k][j] = Q[k]*phi[j];
			}
		}
		
		tmpResult = 0;
		for (l = 0; l < m; l++) {
			for (j = 0; j < m; j++) {
				for (k = 0; k < m; k++) {
					tmpResult += R[l][k]*p[k][j];
				}
				secondTermNumerator[l][j] = tmpResult;
				tmpResult = 0;
			}
		}
		
		for (l=0;l<m;l++){
			tmpResult = 0;
			for (j=0;j<m;j++){
				 tmpResult += p[l][j]*phi[j];
			}
			L[l] = tmpResult;
		}
		
		result = 0;
		for (l=0;l<m;l++){	
			result += L[l]*phi[l];
		}
		
		secondTermDenominator = result + (lambda*lambda);
			
		for (l=0;l<m;l++){
			for (j=0;j<m;j++){
				SecondTerm[l][j]=secondTermNumerator[l][j] * ((-1) / secondTermDenominator);
			}
		}

		for (l=0;l<m;l++){
			for (j=0;j<m;j++){
				A[l][j] = p[l][j]+SecondTerm[l][j];
			}
		}
		
		for (l=0;l<m;l++){
			for (j=0;j<m;j++){
				p_[l][j]=A[l][j] * 1/(lambda*lambda);
			}
		}		

		y_=0;
		for (l=0;l<m;l++){	
			y_ += phi[l]*t[l];
		}
		
		for (l=0;l<m;l++){
			tmpResult = 0;
			for (j=0;j<m;j++){
				 tmpResult += p_[l][j]*phi[j];
			}
			B[l] = tmpResult;
		}
		
		for (l=0;l<m;l++){
			C[l] = B[l]*(yd_rid[i]-y_);
		}
		
		for (l=0;l<m;l++){
			t_[l] = t[l] + C[l];
		}
		
		
		// swapVars
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
void printMatrix(float mat[m][m]){
	int i,j;
	for (i=0;i<m;i++){
		for (j=0;j<m;j++){
			printf("[%.6f],",mat[i][j]);
		}
		printf("\n");
	}
}

void printArray(float arr[], int size){
	int i=0;
	for(i=0;i<size;i++){
		printf(" [%.6f] ", arr[i]);
	}
	printf("\n");
}

void initializeData(){
	ud_rid[0]=-2.99316406250000;
	ud_rid[1]=-2.99316406250000;
	ud_rid[2]=-2.99316406250000;
	ud_rid[3]=-2.99804687500000;
	ud_rid[4]= 3.02246093750000;
	ud_rid[5]= 3.01269531250000;
	ud_rid[6]= 3.01269531250000;
	ud_rid[7]= 3.01269531250000;
	ud_rid[8]= 3.00781250000000;
	ud_rid[9]= 3.00781250000000;
	ud_rid[10]= 3.00292968750000;
	ud_rid[11]= 3.00781250000000;
	ud_rid[12]= 3.00292968750000;
	ud_rid[13]= 3.00781250000000;
	ud_rid[14]= 3.00292968750000;
	ud_rid[15]= 3.00292968750000;
	ud_rid[16]= 3.00292968750000;
	ud_rid[17]= 3.00292968750000;
	ud_rid[18]= 2.99804687500000;
	ud_rid[19]= 3.00292968750000;
	ud_rid[20]= 3.00781250000000;
	ud_rid[21]= 3.00292968750000;
	ud_rid[22]= 3.00781250000000;
	ud_rid[23]= 3.00781250000000;
	ud_rid[24]= 3.00292968750000;
	ud_rid[25]= 3.00292968750000;
	ud_rid[26]= 3.00292968750000;
	ud_rid[27]= 3.00292968750000;
	ud_rid[28]= 3.00292968750000;
	ud_rid[29]= 3.00292968750000;
	ud_rid[30]= 3.00292968750000;
	
	yd_rid[0]=3.01757812500000;
	yd_rid[1]=2.98828125000000;
	yd_rid[2]= 2.97851562500000;
	yd_rid[3]= 2.98828125000000;
	yd_rid[4]= 2.73437500000000;
	yd_rid[5]=-1.16210937500000;
	yd_rid[6]=-5.54687500000000;
	yd_rid[7]=-6.96289062500000;
	yd_rid[8]=-4.84375000000000;
	yd_rid[9]=-1.68457031250000;
	yd_rid[10]=-0.28320312500000;
	yd_rid[11]=-1.49902343750000;
	yd_rid[12]=-3.73046875000000;
	yd_rid[13]=-4.84863281250000;
	yd_rid[14]=-4.16015625000000;
	yd_rid[15]=-2.61230468750000;
	yd_rid[16]=-1.72851562500000;
	yd_rid[17]=-2.13378906250000;
	yd_rid[18]=-3.21777343750000;
	yd_rid[19]=-3.87695312500000;
	yd_rid[20]=-3.63281250000000;
	yd_rid[21]=-2.88574218750000;
	yd_rid[22]=-2.39257812500000;
	yd_rid[23]=-2.53906250000000;
	yd_rid[24]=-3.07128906250000;
	yd_rid[25]=-3.42773437500000;
	yd_rid[26]=-3.33496093750000;
	yd_rid[27]=-2.95898437500000;
	yd_rid[28]=-2.70507812500000;
	yd_rid[29]=-2.76367187500000;
	yd_rid[30]=-3.02734375000000;
}