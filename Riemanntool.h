#define gamma 1.4
#define epsilon 0.000001
#define ERROR 0.0000000000001
#define sq(x) (x*x)

//initial
void Rinitial(Matrix*Infor, double*infor)//Infor, initial information which includes rho p u on both sides, and the amount of mesh unit = n*2l = n*L
{
	int N = Infor->cols;
	double infordata[3*N];
	for(int i = 0; i < 0.5*N; i++)
	{
		infordata[i] = infor[0];
		infordata[1*N + i] = infor[1];
		infordata[2*N + i] = infor[2];
	}
	for(int i = 0.5*N; i < 1*N; i++)
	{
		infordata[i] = infor[3];
		infordata[1*N + i] = infor[4];
		infordata[2*N + i] = infor[5];
	}

	Input(Infor, infordata);

}//test successfully

//basic tool for construction
double Energy(double rho, double p, double u)//compute energy
{
	double E = p/(rho * (gamma - 1)) + 0.5 * sq(u);
	return E;
}//test successfully

double dQ(double q1, double q2, double q3, int flag)//q1 = rho, q2 = rho*u, q3 = rho*E
{
	if(flag == 1)
	{
		return q1;
	}
	else if(flag == 2)
	{
		double p;
		p = (q3/q1 - sq(q2/q1)*0.5)*q1*(gamma - 1);
		return p;
	}
	else if(flag == 3)
	{
		return q2/q1;
	}
	else
	{
		printf("wrong flag!");
		return 0;
	}
}//test successfully

double makeQ(double rho, double p, double u, int flag)//q1 = rho, q2 = rho*u, q3 = rho*E
{
	if(flag == 1)
	{
		return rho;
	}
	else if(flag == 2)
	{
		return rho*u;
	}
	else if(flag == 3)
	{
		return rho*(Energy(rho, p, u));
	}
	else
	{
		printf("wrong flag!");
		return 0;
	}
}//test successfully

double makeF(double rho, double p, double u, int flag)//q1 = rho, q2 = rho*u, q3 = rho*E
{
	if(flag == 1)
	{
		return rho*u;
	}
	else if(flag == 2)
	{
		return rho*sq(u) + p;
	}
	else if(flag == 3)
	{
		return u*(rho*Energy(rho, p, u) + p);
	}
	else
	{
		printf("wrong flag!");
		return 0;
	}
}//test successfully

double makeQF(double q1, double q2, double q3, int flag)//q1 = rho, q2 = rho*u, q3 = rho*E
{
	double u = q2/q1;
	double p = (q3/q1 - sq(q2/q1)*0.5)*q1*(gamma - 1); 
	if(flag == 1)
	{
		return q2;
	}
	else if(flag == 2)
	{
		return q2*u + p;
	}
	else if(flag == 3)
	{
		return (q3 + p) * u;
	}
	else
	{
		printf("wrong flag!");
		return 0;
	}
}//successfully

double Csound(double a1, double a2, double a3, int flag)
{
	if(flag == 1)
	{
		double rho = a1;
		double p = a2;
		double c = sqrt(gamma*p/rho);
		return c;
	}
	else if(flag == 2)
	{
		double rho = a1;
		double p = (a3/a1 - sq(a2/a1)*0.5)*a1*(gamma - 1); 
		double c = sqrt(gamma*p/rho);
		return c;
	}
}

double sp_rad(Matrix*Ain, int flag)
{
	int N = Ain->cols;
	double rA = 0.;
	double rAtmp;
	if(flag == 1)//specially use for infor
	{
		for (int i = 0; i < 1*N; i++)
		{
			rAtmp = fabs(Ain->data[2*N+i]) + Csound(Ain->data[i], Ain->data[1*N+i], Ain->data[2*N+i], flag);
			if(rA < rAtmp)
			{
				rA = rAtmp;
			}
		}
	}
	else if(flag == 2)//use for Q
	{
		for (int i = 0; i < 1*N; i++)
		{
			rAtmp = fabs(Ain->data[1*N+i]/Ain->data[i]) + Csound(Ain->data[i], Ain->data[1*N+i], Ain->data[2*N+i], flag);
			if(rA < rAtmp)
			{
				rA = rAtmp;
			}
		}
	}

	return rA;
}


//construction tool
//three row of Infor is rho p u respectively
void RQoutput(Matrix*Infor, Matrix*Q)//acquire Q based on Infor
{
	int N = Infor->cols;
	//double Qdata[3*N];
	for(int i = 0; i < 1*N; i++)
	{
		Q->data[i] = makeQ(Infor->data[i], Infor->data[1*N+i], Infor->data[2*N+i],1);
		Q->data[1*N+i] = makeQ(Infor->data[i], Infor->data[1*N+i], Infor->data[2*N+i],2);//u is on the third row of the Infor
		Q->data[2*N+i] = makeQ(Infor->data[i], Infor->data[1*N+i], Infor->data[2*N+i],3);
	}

	//Input( Q, Qdata);//here I use Input
}//test successfully

void RFoutput(Matrix*Infor, Matrix*F)//i->rho, 100*n+i->p, 200*n+i->u
{
	int N = Infor->cols;
	//double Fdata[3*N];
	for(int i = 0; i < 1*N; i++)
	{
		F->data[i] = makeF(Infor->data[i], Infor->data[1*N+i], Infor->data[2*N+i],1);
		F->data[1*N+i] = makeF(Infor->data[i], Infor->data[1*N+i], Infor->data[2*N+i],2);//u is on the third row of the Infor
		F->data[2*N+i] = makeF(Infor->data[i], Infor->data[1*N+i], Infor->data[2*N+i],3);
	}

	//Input( F, Fdata);

}//test successfully

void Rinforoutput(Matrix*Q, Matrix*Infor)//acquire Infor from Q
{
	int N = Infor->cols;
	for(int i = 0; i < 1*N; i++)
	{
		Infor->data[i] = dQ(Q->data[i], Q->data[1*N+i], Q->data[2*N+i], 1);
		Infor->data[1*N + i] = dQ(Q->data[i], Q->data[1*N+i], Q->data[2*N+i], 2);
		Infor->data[2*N + i] = dQ(Q->data[i], Q->data[1*N+i], Q->data[2*N+i], 3);
	}
}//test successfully

void RFoutput2(Matrix*Q, Matrix*F)//construct F from Q
{
	int N = Q->cols;
	for(int i = 0; i < 1*N; i++)
	{
		F->data[i] = makeQF(Q->data[i], Q->data[1*N+i], Q->data[2*N+i], 1);
		F->data[1*N + i] = makeQF(Q->data[i], Q->data[1*N+i], Q->data[2*N+i], 2);
		F->data[2*N + i] = makeQF(Q->data[i], Q->data[1*N+i], Q->data[2*N+i], 3);
	}
}//test successfully

// a new matrix add funcition
Matrix*translation(Matrix*matrix, int flag)//used for Qi+1,Qi-1, make it a matrix add rather than for loop
{
	Matrix*out = initial(out, matrix->row, matrix->cols);
	matrixCpy(matrix, out);
	int N = matrix->cols;
	int R = matrix->row;
	if(flag > 0.5*N)
	{
		printf("warning! the translation may exceed acceptable region!");
	}
	if (flag > 0)
	{
		for(int i = 0; i < N - flag; i++)
		{
			for(int j = 0; j < R; j++)
			{
				out->data[i+j*N] = matrix->data[j*N+i+flag];
			}
		}
	}
	else if(flag < 0)
	{
		for(int i = N + flag; i > -flag; i--)//bear in mind that flag is a nagtive number
		{
			for(int j = 0; j < R; j++)
			{
				out->data[i+j*N] = matrix->data[j*N+i+flag];
			}
			//out->data[i] = matrix->data[i+flag];
			//out->data[N+i] = matrix->data[N+i+flag];
			//out->data[2*N+i] = matrix->data[2*N+i+flag];
		}
	}
	return out;
}//test successfully
//remember to free after using

Matrix*thetaoutput(Matrix*Q)//theta has been initialized!
{
	int n = Q->cols;
	Matrix*theta = initial(theta, 1, n);
	Matrix*rho = initial(rho, 1, n);
	Matrix*rho1, *rhom1, *temp0, *temp1, *temp2, *temp3;

	for(int i = 0; i < n; i++)//copy rho
	{
		rho->data[i] = Q->data[i];
	}

	rho1 = translation(rho, 1);
	rhom1 = translation(rho, -1);
	temp0 = matrixAdd(rho1,rho,-2);
	temp1 = matrixAdd(temp0, rhom1, 1);
	matrixAbs(temp1);
	temp2 = matrixAdd(rho1, rho, -1);
	matrixAbs(temp2);
	temp3 = matrixAdd(rho, rhom1, -1);
	matrixAbs(temp3);
	for(int i = 0; i < n; i++)
	{
		theta->data[i] = temp1->data[i]/(temp2->data[i] + temp3->data[i] + epsilon);
	}

	Matrixfree(rho);
	Matrixfree(rho1);
	Matrixfree(rhom1);
	Matrixfree(temp0);
	Matrixfree(temp1);
	Matrixfree(temp2);
	Matrixfree(temp3);
	return theta;
}//remember to free after using

double makedt(double CFL, int n, Matrix*Q)
{
	double spec = sp_rad(Q, 2);
	double dx = ( double )1/n;
	double dt = CFL*dx/spec;
	return dt;	
}

//function below used for anayltical resolution

double ffunction(double pbar,double p,double rho)
{
	double f, temp1, temp2, c;
	c = Csound(rho, p, 0, 1);//when flag = 1, there is no need for the third arguement
	if(pbar > p)
	{
		temp1 = pow((gamma+1)*pbar/(2*gamma*p) + (gamma - 1)/(2*gamma), 0.5);
		temp2 = pbar - p;
		temp1 = temp1*c*rho;
		f = temp2/temp1;
	}
	else if(pbar < p)
	{
		temp1 = pow((pbar/p), (gamma-1)/(2*gamma));
		temp2 = 2*c/(gamma - 1);
		temp1 = temp1 - 1;
		f = temp1*temp2;
	}
	else
	{
		f = 0;
	}

	return f;
}

double Ffunction(double p,double p1,double rho1,double p2,double rho2)//F function
{
	return ffunction(p, p1, rho1)+ffunction(p, p2, rho2);
}

int identification(double*infor)
{
	double rhoL = infor[0];
	double pL = infor[1];
	double uL = infor[2];
	double rhoR = infor[3];
	double pR = infor[4];
	double uR = infor[5];

	double Fp = uL - uR;
	
	double FL = Ffunction(pL,pL,rhoL,pR,rhoR);
	double FR = Ffunction(pR,pL,rhoL,pR,rhoR);
	double F0 = Ffunction(0,pL,rhoL,pR,rhoR);
	if(Fp > FL && Fp > FR)
	{
		return 1;//double shock wave
	}
	else if(Fp < FL && Fp < FR && Fp >= F0)
	{
		return 2;//double expansion
	}
	else if(pR > pL && Fp < FR && Fp > FL)
	{
		return 3;//right expansion, left shock
	}
	else if(pL > pR && Fp < FL && Fp > FR)
	{
		return 4;//right shock, left expansion
	}
	else if(uL == uR && pL == pR)
	{
		return 5;
	}
	else
	{
		return 0;
	}

}//test successfully!

double pbarmake(double*infor)//output p after wave
{
	double rhoL = infor[0];
	double pL = infor[1];
	double uL = infor[2];
	double rhoR = infor[3];
	double pR = infor[4];
	double uR = infor[5];
	double Fp = uL - uR;

	double p1 = 0;
	double p2 = 10000000;//thus we cannot solve p bigger than this.
	double pbar, F1, F2;
	while(fabs(p2 - p1) > ERROR)
	{
		F1 = Ffunction(p1,pL,rhoL,pR,rhoR);
		F2 = Ffunction(p2,pL,rhoL,pR,rhoR);
		pbar = 0.5*(p1 + p2);
		double Fbar = Ffunction(pbar,pL,rhoL,pR,rhoR); 
		if(Fbar > Fp)
		{
			p2 = pbar;
		}
		else if(Fbar < Fp)
		{
			p1 = pbar;
		}
		else if(Fbar == Fp)
		{
			break;
		}
	}

	return pbar;
}

double ubarmake(double pbar,double*infor)//output u after wave
{
	double rhoL = infor[0];
	double pL = infor[1];
	double uL = infor[2];
	double rhoR = infor[3];
	double pR = infor[4];
	double uR = infor[5];

	double ubar = 0.5*(uL + uR - ffunction(pbar, pL, rhoL) + ffunction(pbar, pR, rhoR)); 

	return ubar;
}

//shock
double shockAmake(double pbar, double rho, double p)
{
	double A = rho*Csound(rho, p, 0, 1)*sqrt((gamma + 1)*pbar/(2*gamma*p) + (gamma-1)/(2*gamma));

}

double shockumake(double pbar, double ubar, double rho, double p, double u, int flag)//output speed of shock
{
	double shocku;
	if(flag == 1)
	{
		shocku = u - shockAmake(pbar, rho, p)/rho;
	}
	else if(flag == 2)
	{
		shocku = u + shockAmake(pbar, rho, p)/rho;
	}
	//double shocku = u - (pbar - p)/(rho*(u - ubar)); 
	return shocku;
}


double shockrhomake(double pbar, double ubar, double rho, double p, double u, int flag)
{
	double shockrho;
	if(flag == 1)
	{
		shockrho = rho*shockAmake(pbar, rho, p)/(shockAmake(pbar, rho, p) - rho*(u - ubar));
	}
	else if(flag == 2)
	{
		shockrho = rho*shockAmake(pbar, rho, p)/(shockAmake(pbar, rho, p) + rho*(u - ubar));
	}
	//double shockrho = rho*(u - shocku)/(ubar - shocku);
	return shockrho;
}

double*shockinformake(double*infor, int flag)//output everything about shock, left1, right 2
{
	double rhoL = infor[0];
	double pL = infor[1];
	double uL = infor[2];
	double rhoR = infor[3];
	double pR = infor[4];
	double uR = infor[5];

	double*shockinfor = (double*)malloc(4*sizeof(double));

	shockinfor[1] = pbarmake(infor);//pbar
	shockinfor[2] = ubarmake(shockinfor[1],infor);//ubar
	if(flag == 1)
	{
		shockinfor[3] = shockumake(shockinfor[1], shockinfor[2], rhoL, pL, uL, flag);//speed of shock
		shockinfor[0] = shockrhomake(shockinfor[1], shockinfor[2], rhoL, pL, uL, flag);//rho after wave
	}
	else if(flag == 2)
	{
		shockinfor[3] = shockumake(shockinfor[1], shockinfor[2], rhoR, pR, uR, flag);
		shockinfor[0] = shockrhomake(shockinfor[1], shockinfor[2], rhoR, pR, uR, flag);
	}

	return shockinfor;
}

//expansion
/*
double expanrhomake(double pbar, double rho, double p)//this is only needed when expan and shock exist at the same time
{
	double expanrho = rho*pow(pbar/p, 1/gamma);
	return expanrho;
}
*/
double expanCmake(double*infor,double ubar, int flag)
{
	double rhoL = infor[0];
	double pL = infor[1];
	double uL = infor[2];
	double rhoR = infor[3];
	double pR = infor[4];
	double uR = infor[5];
	double Cstar;
	
	if(flag == 1)
	{
		Cstar = Csound(rhoL, pL, 0, 1) + 0.5*(gamma - 1)*(uL - ubar);
	}
	else if(flag == 2)
	{
		Cstar = Csound(rhoR, pR, 0, 1) - 0.5*(gamma - 1)*(uR - ubar);
	}

	return Cstar;
}

double expanrhomake(double pbar, double Cstar)
{
	double expanrho = gamma*pbar/sq(Cstar);
	return expanrho;
}

double expanumake(double rho, double p, double u, int flag)//left flag = 1, right = 2
{
	if(flag == 1)
	{
		return u - Csound(rho, p, 0, 1);
	}
	else if(flag == 2)
	{
		return u + Csound(rho, p, 0, 1);
	}
}

double expaninCmake(double*infor, double t, double x, int flag)//left 1, right 2, output sound speed inside expan wave
{
	double rhoL = infor[0];
	double pL = infor[1];
	double uL = infor[2];
	double rhoR = infor[3];
	double pR = infor[4];
	double uR = infor[5];
	double expanC;
	if(flag == 1)//left
	{
		expanC = (gamma - 1)/(gamma + 1)*(uL - x/t) + 2*Csound(rhoL, pL, 0, 1)/(gamma + 1);
	}
	else if(flag == 2)//right
	{
		expanC = -1*(gamma - 1)/(gamma + 1)*(uR - x/t) + 2*Csound(rhoR, pR, 0, 1)/(gamma + 1);
	}
	return expanC;
}

double*expaninsideinformake(double*infor, double t, double x, int flag)
{
	double rhoL = infor[0];
	double pL = infor[1];
	//double uL = infor[2];
	double rhoR = infor[3];
	double pR = infor[4];
	//double uR = infor[5];

	double*expaninsideinfor = (double*)malloc(3*sizeof(double));
	double expanC = expaninCmake(infor, t, x, flag);
	if(flag == 1)
	{
		expaninsideinfor[1] = pL*pow(expanC/Csound(rhoL,pL,0,1),2*gamma/(gamma - 1));//p
		expaninsideinfor[0] = gamma*expaninsideinfor[1]/sq(expanC);//rho
		expaninsideinfor[2] = x/t + expanC;//u
	}
	else if(flag == 2)
	{
		expaninsideinfor[1] = pR*pow(expanC/Csound(rhoR,pR,0,1),2*gamma/(gamma - 1));//p
		expaninsideinfor[0] = gamma*expaninsideinfor[1]/sq(expanC);
		expaninsideinfor[2] = x/t - expanC;
	}
	return expaninsideinfor;
}

double*expaninformake(double*infor,int flag)
{
	double rhoL = infor[0];
	double pL = infor[1];
	double uL = infor[2];
	double rhoR = infor[3];
	double pR = infor[4];
	double uR = infor[5];

	double*expaninfor = (double*)malloc(5*sizeof(double));

	expaninfor[1] = pbarmake(infor);//p after wave
	expaninfor[2] = ubarmake(expaninfor[1],infor);//u after wave
	double Cstar = expanCmake(infor, expaninfor[2], flag);

	if(flag == 1)
	{
		expaninfor[0] = expanrhomake(expaninfor[1],Cstar);//after wave rho
		expaninfor[3] = uL - Csound(rhoL,pL,0,1);//wavehead speed
		expaninfor[4] = expaninfor[2] - Cstar;//wavetail speed
	}
	else if(flag == 2)
	{
		expaninfor[0] = expanrhomake(expaninfor[1],Cstar);//after wave rho
		expaninfor[3] = uR + Csound(rhoR,pR,0,1);//wavehead speed
		expaninfor[4] = expaninfor[2] + Cstar;//wavetail speed
	}
	
	return expaninfor;
}
