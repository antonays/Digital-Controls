// First Order Kalman Estimator - With Pointers
// Rhaeg

#include <stdio.h>
#include <stdlib.h>

#define N 20

double ud_rid [N];
double yd_rid [N];

void initializeData();
void printArray(double *arr, int size);
double **MatrixByMatrix(double **a,double **b, int dim);
double *vectorByMatrix(int dim, double **matrix, double *vec);
double **allocateMat(int rows, int columns);
double *allocateVec(int dim);
double **initiateIdentity(int dim);
double **scalarMatrix(double scalar, double **mat, int dim);
void printMatrix(double **mat, int dim);
double *MatrixByVector(int dim, double *vector, double **matrix);
double **VectorByVectorCR(int dim, double *x, double *y);
double **MatrixMatrixAddition(double **mat1, double **mat2, int dim);
double VectorByVectorRC(double *v1, double *v2,int dim);
double *vectorByScalar( double *v, double scalar, int dim);
double *vectorVectorAddition(double *v1, double *v2, int dim);
double *randomize(int dim, int randomFactor);
double *getPhi(int index, int n);

void main(){
	double y_=0, y=0, lambda = 0.9;
	int Gamma =100;
	int n=1; //order
	int m = 2*n;
	
	int i;
	double **identity, **p, **p_;
	double *phi, *t, *t_;
	
	initializeData();
	identity = initiateIdentity(m);
	t = randomize(m, 2);
	p = scalarMatrix(Gamma,identity,m);
	
	for (i=1;i<N;i++){
			phi = getPhi(i,n);

			double **secondTermNumerator = MatrixByMatrix( VectorByVectorCR( m,MatrixByVector(m,phi,p),phi ),p,m );

			double secondTermDenominator = VectorByVectorRC( vectorByMatrix( m, p, phi ),phi,m ) + (lambda*lambda);

			p_ = scalarMatrix( 1/(lambda*lambda), MatrixMatrixAddition( p, scalarMatrix((-1)/secondTermDenominator, secondTermNumerator,m),m ), m);

			y_ = VectorByVectorRC(phi, t , m); 

			t_ = vectorVectorAddition(t, vectorByScalar(MatrixByVector(m, phi,p_), (yd_rid[i]-y_), m),m );

		p = p_;
		t = t_;
	}
	
	printf("\n Theta: a0, b0 \n");
	printArray(t_, m);
	
	free(identity);
	free(p);
	free(p_);
	free(phi);
	free(t);
	free(t_);
}



double *getPhi(int index, int n){
	int i,j=0;
	double *phi;
	phi = allocateVec(2*n);
	
	for (i=1;i<n+1;i++){
		phi[j] = (-1)*yd_rid[index-i];
		j++;
	}
	for (i=1;i<n+1;i++){
		phi[j] = ud_rid[index-i];
		j++;
	}
	return phi;
}

double **scalarMatrix(double scalar, double **mat, int dim){
	int i,j;
	double **tmpMat = mat;
	for (i=0;i<dim;i++){
		for (j=0;j<dim;j++)
		{
			tmpMat[i][j]=tmpMat[i][j] * scalar;
		}
	}
	return tmpMat;
}

void printMatrix(double **mat, int dim){
	int i,j;
	for (i=0;i<dim;i++){
		for (j=0;j<dim;j++)
		{
			printf("[%f],",mat[i][j]);
		}
		printf("\n");
	}
}

void printArray(double *arr, int size){
	int i=0;
	for(i=0;i<size;i++){
		printf("%.8f\n", arr[i]);
	}
}

double **initiateIdentity(int dim){
	int i=0,j=0;
	double **mat = allocateMat(dim,dim);
	for (j=0;j<dim;j++){
		for (i=0;i<dim;i++){
			mat[i][j]=0;
		}
	}
	for (i=0;i<dim;i++){
		mat[i][i]=1; 	
	}
	return mat;
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

double **MatrixByMatrix(double **mat1,double **mat2, int dim) {
   int i,j,k;
   double sum=0;
   double **result;
   result = allocateMat(dim,dim);
	for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
        for (k = 0; k < dim; k++) {
          sum = sum + mat1[i][k]*mat2[k][j];
        }
        result[i][j] = sum;
        sum = 0;
      }
    }
	return result;
}

double **allocateMat(int rows, int columns){
	int i;
	double **matrix;
	matrix = (double **)malloc(rows*sizeof(double*));
    if (matrix== NULL)
    {
        printf("Allocation Error");
        return 0;
    }
    for (i=0; i<rows; i++)
    {
        matrix[i] = (double *)malloc(columns*sizeof(double));
        if (matrix[i]== NULL)
        {
            printf("Allocation Error");
            return 0;
        }
    }
	return matrix;
}

double *allocateVec(int dim){
	double *vec = (double *)malloc(dim*sizeof(double));
    if (vec == NULL)
    {
        printf("Allocation Error");
        return 0;
    }
	return vec;
}

double *vectorByMatrix(int dim, double **matrix, double *vec){
	int i,j;
	double *result;
	result = allocateVec(dim);
	double tmpresult=0;
	for (i=0;i<dim;i++){
		tmpresult=0;
		for (j=0;j<dim;j++)
		{
			tmpresult += matrix[i][j]*vec[j];
		}
		result[i] = tmpresult;
	}
	return result;
}

double **TransposeMatrix(double **mat, int dim){
	int i,j;
	double **tmpMat = allocateMat(dim,dim);
	for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++)
         tmpMat[i][j] = mat[j][i];
   }
   return tmpMat;
}

double *MatrixByVector(int dim, double *vector, double **matrix)
{
	int i,j;
	double *tmpResultVec, tmpResult;
	tmpResultVec = allocateVec(dim);
	for (i=0;i<dim;i++){
		tmpResult=0;
		for (j=0;j<dim;j++)
		{
			tmpResult += matrix[i][j]*vector[j];
		}
		tmpResultVec[i] = tmpResult;
	}
	return tmpResultVec;
}

double **VectorByVectorCR(int dim, double *x, double *y)
{
	double **tmpMatrix, result=0;
	tmpMatrix = allocateMat(dim,dim);
	int i,j;
	for (i=0;i<dim;i++){
		for (j=0;j<dim;j++){
			tmpMatrix[i][j] = x[i]*y[j];
		}
	}
	return tmpMatrix;
}

double **MatrixMatrixAddition(double **mat1, double **mat2, int dim){
	int i,j;
	double **tmpMat;
	tmpMat = allocateMat(dim,dim);
	for (i=0;i<dim;i++){
		for (j=0;j<dim;j++){
			tmpMat[i][j] = mat1[i][j]+mat2[i][j];
		}
	}
	return tmpMat;
}


double VectorByVectorRC(double *v1, double *v2,int dim){
	int i;
	double result=0;
	for (i=0;i<dim;i++){	
		result += v1[i]*v2[i];
	}
	return result;
}

double *randomize(int dim, int randomFactor){
	int i=0;
	double *v;
	v = allocateVec(dim);
	for (i=0;i<dim;i++){
		v[i] = rand() % randomFactor;
	}
	return v;
}

double *vectorByScalar( double *v, double scalar, int dim){
	int i=0;
	double *v2 = allocateVec(dim);
	for (i=0;i<dim;i++){
		v2[i] = v[i]*scalar;
	}
	return v2;
}

double *vectorVectorAddition(double *v1, double *v2, int dim){
	int i=0;
	double *v3;
	v3 = allocateVec(dim);
	for (i=0;i<dim;i++){
		v3[i] = v1[i] + v2[i];
	}
	return v3;	
}