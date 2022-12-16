#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrixtool.h"
#include "Riemanntool.h"
#include "predictor.h"
#include "filter.h"
#include "corrector.h"
#include "analoop.h"

Matrix*loop(Matrix*Q, int n, double CFL, double beta, double t)
{
	Matrix*Qbar, *Qwave;
	double ct = 0;
	double dt;
	while(ct < t)
	{
		dt = makedt(CFL, n, Q);
		Qbar = MaCfilter(Q, CFL, beta);
		Qwave = MaCpredictor(Qbar, dt, n);
		Matrixfree(Q);
		Q = MaCcorrector(Qbar, Qwave, dt, n);

		Matrixfree(Qbar);
		Matrixfree(Qwave);
		ct += dt;
	}
	return Q;
}

int main(int argc, char** argv)
{
	//input section
	int n, l;
	double CFL, t, beta;
	char name[20];
	char inforout[100];
	char output[100];
	char anaout[100];
	printf("please input mesh number\n");
	scanf("%d", &n);
	printf("please input length\n");
	scanf("%d", &l);
	int N = 2*n*l;
	printf("please input CFL number\n");
	scanf("%lf", &CFL);
	printf("please input ending time\n");
	scanf("%lf", &t);
	printf("please input beta\n");
	scanf("%lf", &beta);
	//printf("please output your file name:\n");
	//scanf("%s", name);
	sprintf(name, "%s", argv[1]);
	sprintf(inforout, "%s/%s_informatiom.txt", name, name);
	sprintf(output, "%s/%s_data.txt", name,name);
	sprintf(anaout, "%s/%s_anadata.txt", name,name);
	FILE*fp1, *fp2, *fp3;
	//intput infromation txt file
	double infor[6];
	double input;
	Matrix*Infor = initial(Infor, 3, N+1);
	Matrix*Q = initial(Q, 3, N+1);


	char initialinf[6][15] = {"Left rho", "Left p", "Left u", "Right rho", "Right p", "Right u"};
	

	//information file output
	fp1 = fopen(inforout,"w");
	fprintf(fp1, "number of mesh perunit = %d\n", n);
	fprintf(fp1, "number of length = %d\n", l);
	fprintf(fp1, "CFL number = %lf\n", CFL);
	fprintf(fp1, "terminal time = %lf\n", t);
	fprintf(fp1, "Diffusion coefficient = %lf\n", beta);


	for(int i = 0; i < 6; i++)
	{
		printf("please input initial value %s\n", initialinf[i]);
		scanf("%lf", &input);
		infor[i] = input;
		fprintf(fp1, "initial information %s = %lf\n", initialinf[i], infor[i]);
	}
	fclose(fp1);


	//initial section
	Rinitial(Infor, infor);
	RQoutput(Infor, Q);//don't free infor!


	Q = loop(Q, n, CFL, beta, t);
	Rinforoutput(Q, Infor);//so we have infor right here.
	double x[N+1];//index
	fp2 = fopen(output,"w");
	for(int j = 0; j < N+1; j++)
	{
		x[j] = -l + ( double )j*1/n;
		fprintf(fp2, "%lf \t %lf \t %lf \t %lf \n", x[j], Infor->data[j], Infor->data[1*(N+1)+j], Infor->data[2*(N+1)+j]);
	}
	
	Matrix*anadata = analoop(infor, t, l);
	int Na = anadata->cols;
	fp3 = fopen(anaout,"w");
	for(int i = 0; i < Na; i++)
	{
		fprintf(fp3, "%lf \t %lf \t %lf \t %lf \n", anadata->data[3*Na+i], anadata->data[i], anadata->data[1*Na+i], anadata->data[2*Na+i]);
	}
	fclose(fp3);

	return 0;
}
