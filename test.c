#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrixtool.h"
#include "Riemanntool.h"
#include "filter.h"
#include "predictor.h"
#include "analoop.h"

int main()
{
	//Matrix*Infor = initial(Infor, 3, 10);
	//Matrix*Q = initial(Q, 3, 10);
	//Matrix*F = initial(F, 3, 10);
	double infor[6];
	double input;
	
	/*
	for(int i = 0; i < 6; i++)
	{
		infor[i] = i+1;
		//infor[i+3] = i + 1;
	}
	*/
	for(int i = 0; i < 6; i++)
	{
		printf("please input initial value %d\n", i);
		scanf("%lf", &input);
		infor[i] = input;
	}

	//Rinitial(Infor, infor);
	//int i = identification(infor);
	Matrix*out;
	/*
	double*inforwave = expaninformake(infor, 1);
	for(int i = 0; i < 5; i++)
	{
		printf("%lf\n",inforwave[i] );
	}
	*/
	out = analoop(infor, 3.5, 20);
	int N = out->cols;
	FILE*fp;
	fp = fopen("test.txt", "w");
	for(int i = 0; i < N; i++)
	{
		fprintf(fp, "%lf \t %lf \t %lf \t %lf \n", out->data[3*N+i], out->data[i], out->data[1*N+i], out->data[2*N+i]);
	}
	//RQoutput(Infor, Q);
	//RFoutput(Infor, F);
	//Rinforoutput(Q, Infor);
	//RFoutput2(Q, F);
	//matrixPrint(Q);
	//Matrix*out = translation(Infor, -1);
	//Matrix*theta = thetaoutput(Q);
	//Matrix*Qbar = MaCfilter(Q, 1, 1);
	//matrixPrint(Qbar);
	return 0;
}
