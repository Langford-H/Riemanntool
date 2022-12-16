#include <string.h>

typedef struct
{
	int row, cols;
	double*data;
}Matrix;

//function guide
void Matrixfree(Matrix*matrix);
void Input(Matrix*matrix, double*arraydata);//input a array to assign data in matrix
void matrixCpy(Matrix*matrix1, Matrix*matrix2);// (source, destination) you need to input two matrix that has been initialized
void matrixPrint(Matrix*matrix);
int matrixSize(Matrix*matrix);
Matrix*initial(Matrix*matrix, int row, int cols);
Matrix*matrixAdd(Matrix*matrix1, Matrix*matrix2, double sign);//if sign > 0 then it is add else the function preform minus
Matrix*matrixTimes(Matrix*matrix1, Matrix*matrix2);
void matrixTrans2(Matrix*matrix);//matrix transverse
void matrixTrans2(Matrix*matrix);//matrix transverse
void matrixAbs(Matrix*matrix);// take absolute value of all the elements in the matrix 
Matrix*matrixConstimes(Matrix*matrix1, Matrix*matrix2);//matrix1 times matirx2 elemently return a matrix

//basic function
void Matrixfree(Matrix*matrix)
{
	free(matrix->data);
	matrix->data = NULL;
}

void Input(Matrix*matrix, double*arraydata)//input a array to assign data in matrix
{
	if(matrix->data != NULL)
	{
		memcpy(matrix->data, arraydata, sizeof(double)*matrix->row*matrix->cols);
	}
}

void matrixCpy(Matrix*matrix1, Matrix*matrix2)// (source, destination) you need to input two matrix that has been initialized
{
	//matrix2 = initial(matrix2, matrix1->row, matrix1->cols);
	matrix2 -> row = matrix1 -> row;
	matrix2 -> cols = matrix1 -> cols;
	memcpy(matrix2->data, matrix1->data, sizeof(double)*(matrix1->row)*(matrix2->cols));
}

void matrixPrint(Matrix*matrix)
{
	for(int i = 0; i < matrix->row; i++)
	{
		for(int j = 0; j < matrix->cols; j++)
		{
			printf("%lf \t", matrix->data[i*matrix->cols+j]);
		}
		printf("\n");
	}
}

int matrixSize(Matrix*matrix)
{
	if(matrix->row>0 && matrix->cols>0)
	{
		int size =matrix->row * matrix->cols; 
		return size;
	}
	else
	{
		printf("Wrong matrix!");
		return 0;
	}
}

Matrix*initial(Matrix*matrix, int row, int cols)
{
	if(row > 0 && cols > 0)
	{
		/*
		if(matrix->row != 0 && matrix->cols > 0)//prevent assign malloc on one pointer twice
		{
			matrix = (Matrix*)realloc(matrix,sizeof(Matrix));
			matrix -> row = row;
			matrix -> cols = cols;
			matrix -> data = (double*)realloc(matrix->data,sizeof(double)*row*cols);
			memset(matrix->data, 0, sizeof(double)*row*cols); //set every elements in the matrix to be 0
			return matrix;

		}
		need to be done
		*/
		matrix = (Matrix*)malloc(sizeof(Matrix));
		matrix -> row = row;
		matrix -> cols = cols;
		matrix -> data = (double*)malloc(sizeof(double)*row*cols);
		memset(matrix->data, 0, sizeof(double)*row*cols); //set every elements in the matrix to be 0
		return matrix;
	}
	else
	{
		return NULL;
	}
}

//calculation
Matrix*matrixAdd(Matrix*matrix1, Matrix*matrix2, double sign)//if sign > 0 then it is add else the function preform minus
{
	double S;
	if(sign == 0)
	{
		S = -1;
	}
	else
	{
		S = sign;
	}
	if(matrix1->row == matrix2->row && matrix1->cols == matrix2->cols)
	{
		Matrix*matrix3 = initial(matrix3, matrix1->row, matrix1->cols);
		for(int i = 0; i < matrix1 -> row; i++)
		{
			for(int j = 0; j< matrix1 -> cols; j++)
			{
				matrix3->data[i * matrix3 -> cols + j] = matrix1->data[i * matrix3 -> cols + j] + S*matrix2->data[i * matrix3 -> cols + j];
			}
		}

		return matrix3;
	}
	else
	{
		printf("wrong input!");
		return NULL;
	}
}

Matrix*matrixTimes(Matrix*matrix1, Matrix*matrix2)
{
	if(matrix1->cols == matrix2->row)
	{
		Matrix*matrix3 = initial(matrix3, matrix1->row, matrix2->cols);
		for(int i = 0; i < matrix3->row; i++)
		{
			for(int j = 0; j < matrix3->cols; j++)
			{
				for(int k = 0; k < matrix1->cols; k++)
				{
					matrix3->data[i*matrix3->cols + j] += matrix1->data[i*matrix1->cols + k]*matrix2->data[k*matrix2->cols + j];
				}
			}
		}
		
		return matrix3;
	}
	else
	{
		printf("wrong input!");
		return NULL;
	}
}


void matrixTrans2(Matrix*matrix)//matrix transverse
{
	if(matrix->row > 0 && matrix->cols > 0)
	{
		Matrix*matrixTemp = initial(matrixTemp, matrix->row, matrix->cols);
		matrixCpy(matrix, matrixTemp);
		matrix->row = matrixTemp->cols;
		matrix->cols = matrixTemp->row;

		for(int i = 0; i < matrix->row; i++)
		{
			for(int j = 0; j < matrix->cols; j++)
			{
				matrix->data[i*matrix->cols+j] = matrixTemp->data[j*matrixTemp->cols+i];
			}
		}
	}
	else
	{
		printf("Wrong input! One dimension is less than zero!");
	}
}

Matrix* matrixTrans(Matrix*matrix)//matrix transverse which will return a Matrix
{
	if(matrix->row > 0 && matrix->cols > 0)
	{
		Matrix*matrixTemp = initial(matrixTemp, matrix->cols, matrix->row);

		for(int i = 0; i < matrixTemp->row; i++)
		{
			for(int j = 0; j < matrixTemp->cols; j++)
			{
				matrixTemp->data[i*matrixTemp->cols+j] = matrix->data[j*matrix->cols+i];
			}
		}
		return matrixTemp;
	}
	else
	{
		printf("Wrong input! One dimension is less than zero!");
		return NULL;
	}
}

void matrixAbs(Matrix*matrix)
{
	int n = matrixSize(matrix);
	for(int i = 0; i < n; i++)
	{
		if(matrix->data[i] < 0)
		{
			matrix->data[i] = -matrix->data[i];
		}
	}
}

Matrix*matrixConstimes(Matrix*matrix1, Matrix*matrix2)//matrix1 times matirx2 elemently return a matrix
{
	if (matrix1->cols != matrix2->cols)
	{
		printf("Worng intput!");
		return NULL;
	}

	Matrix*out = initial(out , matrix1->row, matrix1->cols);
	if(matrix1->row == matrix2->row)
	{
		for(int i = 0; i < matrixSize(matrix1); i++)
		{
			out->data[i] = matrix1->data[i] * matrix2->data[i];
		}
	}
	else
	{
		for(int i = 0; i < matrix1->cols; i++)
		{
			for(int j = 0; j < matrix1->row; j++)
			{
				out->data[j*matrix1->cols+i] = matrix1->data[j*matrix1->cols+i] * matrix2->data[i];
			}
		}
	}

	return out;
}
