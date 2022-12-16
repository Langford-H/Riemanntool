Matrix*MaCcorrector(Matrix* Qbar, Matrix* Qwave, double dt, int n)
{
	Matrix*Qn1;
	int N = Qbar->cols;
	Matrix*F1 = initial(F1 , 3, N);
	Matrix*F2 = initial(F2 , 3, N);
	Matrix*temp1, *temp2, *temp3;

	temp1 = translation(Qwave, 1);//temp1 = Qwave1
	
	RFoutput2(Qwave, F1);//F(Qi)
	RFoutput2(temp1, F2);//F(Qi+1)
	Matrixfree(temp1);
	temp2 = matrixAdd(F2, F1, -1);
	temp3 = matrixAdd(Qbar, Qwave, 1);//Qwave+Qbar 
	temp1 = matrixAdd(temp3, temp3, -0.5);//0.5(Qwave+Qbar)

	Qn1 = matrixAdd(temp1, temp2, -0.5*dt*n);

	Matrixfree(temp1);
	Matrixfree(temp2);
	Matrixfree(temp3);
	Matrixfree(F1);
	Matrixfree(F2);
	return Qn1;
}
