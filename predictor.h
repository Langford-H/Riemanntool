Matrix*MaCpredictor(Matrix*Qbar, double dt, int n)
{
	int N = Qbar->cols;
	Matrix*Qwave;
	Matrix*F1 = initial(F1 , 3, N);
	Matrix*F2 = initial(F2 , 3, N);
	Matrix*temp1, *temp2;
	
	temp1 = translation(Qbar, -1);//temp = Qbar-1
	RFoutput2(Qbar, F1);//F(Qi)
	RFoutput2(temp1, F2);//F(Qi-1)
	temp2 = matrixAdd(F1, F2, -1);
	Qwave = matrixAdd(Qbar, temp2, -dt*n);

	Matrixfree(temp1);
	Matrixfree(temp2);
	Matrixfree(F1);
	Matrixfree(F2);
	return Qwave;
}
