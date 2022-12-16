Matrix*MaCfilter(Matrix* Q, double CFL, double beta)
{
	int N = Q->cols;
	Matrix*out, *theta;
	Matrix*temp3 = thetaoutput(Q);//create theta
	Matrix*temp1, *temp2;
  theta = matrixAdd(temp3, temp3, CFL*(1-CFL)*beta - 1);//update theta to nu

	temp1 = translation(Q, 1);
	temp2 = translation(Q, -1);
	out = matrixAdd(temp1, Q, -2);
	Matrixfree(temp1);//temp1 as Qi+1 will be freed
	temp1 = matrixAdd(out, temp2,1);//now temp1 is Qi+1, Qi, Qi-1
	Matrixfree(out);
	Matrixfree(temp2);
	temp2 = matrixConstimes(temp1, theta);//now temp2 = nui*(Qi+1-2Qi+Qi-1)
	out = matrixAdd(Q, temp2, 0.5);


	Matrixfree(theta);
	Matrixfree(temp1);
	Matrixfree(temp2);
	Matrixfree(temp3);
	return out;
}
// remember Q need to be free, this return bar{Qi}
