Matrix* analoop(double*infor, double t, int l)//input initial condition and terminal time, output matrix
{
	int na = 50;
	int Na = 2*na*l;
	int id = identification(infor);
	Matrix*out = initial(out, 4, Na+1);

	double*leftwave, *rightwave;
	double x[Na + 1];

	if(id == 1)//double shock
	{
		leftwave = shockinformake(infor, 1);
		rightwave = shockinformake(infor, 2);
		double rightx, leftx, middlex;
		rightx = rightwave[3]*t;
		leftx = leftwave[3]*t;
		middlex = rightwave[2]*t;
		for(int i = 0; i < Na+1; i++)
		{
			x[i] = -l + (double)i*1/na;
			out->data[3*(Na+1)+i] = x[i];
			if(x[i] <= leftx)
			{
				out->data[i] = infor[0];
				out->data[i+Na+1] = infor[1];
				out->data[i+2*(Na+1)] = infor[2];
			}
			else if(x[i] > rightx)
			{
				out->data[i] = infor[3];
				out->data[i+Na+1] = infor[4];
				out->data[i+2*(Na+1)] = infor[5];
			}
			else if(x[i] > leftx && x[i] <= middlex)
			{
				out->data[i] = leftwave[0];
				out->data[i+Na+1] = leftwave[1];
				out->data[i+2*(Na+1)] = leftwave[2];
			}
			else if(x[i] > middlex && x[i] <= rightx)
			{
				out->data[i] = rightwave[0];
				out->data[i+Na+1] = rightwave[1];
				out->data[i+2*(Na+1)] = rightwave[2];
			}
		}
		
	}
	else if(id == 2)//double expan
	{
		leftwave = expaninformake(infor, 1);
		rightwave = expaninformake(infor, 2);
		double righthx, righttx, lefthx, lefttx, middlex;
		lefthx = leftwave[3]*t;
		lefttx = leftwave[4]*t;
		righthx = rightwave[3]*t;
		righttx = rightwave[4]*t;
		middlex = rightwave[2]*t;

		for(int i = 0; i < Na+1; i++)
		{
			x[i] = -l + (double)i*1/na;
			out->data[3*(Na+1)+i] = x[i];
			if(x[i] <= lefthx)
			{
				out->data[i] = infor[0];
				out->data[i+Na+1] = infor[1];
				out->data[i+2*(Na+1)] = infor[2];
			}
			else if(x[i] > righthx)
			{
				out->data[i] = infor[3];
				out->data[i+Na+1] = infor[4];
				out->data[i+2*(Na+1)] = infor[5];
			}
			else if(x[i] > lefttx && x[i] <= middlex)
			{
				out->data[i] = leftwave[0];
				out->data[i+Na+1] = leftwave[1];
				out->data[i+2*(Na+1)] = leftwave[2];
			}
			else if(x[i] > middlex && x[i] <= righttx)
			{
				out->data[i] = rightwave[0];
				out->data[i+Na+1] = rightwave[1];
				out->data[i+2*(Na+1)] = rightwave[2];
			}
			else if(x[i] > lefthx && x[i] <= lefttx)
			{
				double*expaninsideinfor = expaninsideinformake(infor, t, x[i], 1); 
				out->data[i] = expaninsideinfor[0];
				out->data[i+Na+1] = expaninsideinfor[1];
				out->data[i+2*(Na+1)] = expaninsideinfor[2];
				free(expaninsideinfor);
			}
			else if(x[i] > righttx && x[i] <= righthx)
			{
				double*expaninsideinfor = expaninsideinformake(infor, t, x[i], 2); 
				out->data[i] = expaninsideinfor[0];
				out->data[i+Na+1] = expaninsideinfor[1];
				out->data[i+2*(Na+1)] = expaninsideinfor[2];
				free(expaninsideinfor);
			}
		}
	}
	else if(id == 3)//right expan, left shock
	{
		leftwave = shockinformake(infor, 1);//shcok
		rightwave = expaninformake(infor, 2);//expan
		double righthx, righttx, leftx, middlex;
		leftx = leftwave[3]*t;
		righthx = rightwave[3]*t;
		righttx = rightwave[4]*t;
		middlex = rightwave[2]*t;

		for(int i = 0; i < Na+1; i++)
		{
			x[i] = -l + (double)i*1/na;
			out->data[3*(Na+1)+i] = x[i];
			if(x[i] <= leftx)
			{
				out->data[i] = infor[0];
				out->data[i+Na+1] = infor[1];
				out->data[i+2*(Na+1)] = infor[2];
			}
			else if(x[i] > righthx)
			{
				out->data[i] = infor[3];
				out->data[i+Na+1] = infor[4];
				out->data[i+2*(Na+1)] = infor[5];
			}
			else if(x[i] > leftx && x[i] <= middlex)
			{
				out->data[i] = leftwave[0];
				out->data[i+Na+1] = leftwave[1];
				out->data[i+2*(Na+1)] = leftwave[2];
			}
			else if(x[i] > middlex && x[i] <= righttx)
			{
				out->data[i] = rightwave[0];
				out->data[i+Na+1] = rightwave[1];
				out->data[i+2*(Na+1)] = rightwave[2];
			}
			else if(x[i] > righttx && x[i] <= righthx)
			{
				double*expaninsideinfor = expaninsideinformake(infor, t, x[i], 2); 
				out->data[i] = expaninsideinfor[0];
				out->data[i+Na+1] = expaninsideinfor[1];
				out->data[i+2*(Na+1)] = expaninsideinfor[2];
				free(expaninsideinfor);
			}
		}

	}
	else if(id == 4)//right shock, left expan
	{
		leftwave = expaninformake(infor, 1);
		rightwave = shockinformake(infor, 2);
		double rightx, lefthx, lefttx, middlex;
		lefthx = leftwave[3]*t;
		lefttx = leftwave[4]*t;
		rightx = rightwave[3]*t;
		middlex = rightwave[2]*t;

		for(int i = 0; i < Na+1; i++)
		{
			x[i] = -l + (double)i*1/na;
			out->data[3*(Na+1)+i] = x[i];
			if(x[i] <= lefthx)
			{
				out->data[i] = infor[0];
				out->data[i+Na+1] = infor[1];
				out->data[i+2*(Na+1)] = infor[2];
			}
			else if(x[i] > rightx)
			{
				out->data[i] = infor[3];
				out->data[i+Na+1] = infor[4];
				out->data[i+2*(Na+1)] = infor[5];
			}
			else if(x[i] > lefttx && x[i] <= middlex)
			{
				out->data[i] = leftwave[0];
				out->data[i+Na+1] = leftwave[1];
				out->data[i+2*(Na+1)] = leftwave[2];
			}
			else if(x[i] > middlex && x[i] <= rightx)
			{
				out->data[i] = rightwave[0];
				out->data[i+Na+1] = rightwave[1];
				out->data[i+2*(Na+1)] = rightwave[2];
			}
			else if(x[i] > lefthx && x[i] <= lefttx)
			{
				double*expaninsideinfor = expaninsideinformake(infor, t, x[i], 1); 
				out->data[i] = expaninsideinfor[0];
				out->data[i+Na+1] = expaninsideinfor[1];
				out->data[i+2*(Na+1)] = expaninsideinfor[2];
				free(expaninsideinfor);
			}
		}

	}

	else if(id == 5)
	{
		double middlex = infor[2]*t;
		for(int i = 0; i < Na+1; i++)
		{
			x[i] = -l + (double)i*1/na;
			out->data[3*(Na+1)+i] = x[i];
			out->data[i+Na+1] = infor[1];
			out->data[i+2*(Na+1)] = infor[2];
			if(x[i] > middlex)
			{
				out->data[i] = infor[3];
			}
			if(x[i] <= middlex)
			{
				out->data[i] = infor[0];
			}
		}
	}
	
	else if(id == 0)
	{
		leftwave = expaninformake(infor, 1);
		rightwave = expaninformake(infor, 2);
		double righthx, righttx, lefthx, lefttx;
		lefthx = leftwave[3]*t;
		lefttx = (infor[2] + 2*Csound(infor[0], infor[1], 0, 1)/0.4)*t;
		righthx = rightwave[3]*t;
		righttx = (infor[5] + 2*Csound(infor[3], infor[4], 0, 1)/0.4)*t;

		for(int i = 0; i < Na+1; i++)
		{
			x[i] = -l + (double)i*1/na;
			out->data[3*(Na+1)+i] = x[i];
			if(x[i] <= lefthx)
			{
				out->data[i] = infor[0];
				out->data[i+Na+1] = infor[1];
				out->data[i+2*(Na+1)] = infor[2];
			}
			else if(x[i] > righthx)
			{
				out->data[i] = infor[3];
				out->data[i+Na+1] = infor[4];
				out->data[i+2*(Na+1)] = infor[5];
			}
			else if(x[i] > lefttx && x[i] <= righttx)
			{
				out->data[i] = 0;
				out->data[i+Na+1] = 0;
				out->data[i+2*(Na+1)] = 0;
			}
			else if(x[i] > lefthx && x[i] <= lefttx)
			{
				double*expaninsideinfor = expaninsideinformake(infor, t, x[i], 1); 
				out->data[i] = expaninsideinfor[0];
				out->data[i+Na+1] = expaninsideinfor[1];
				out->data[i+2*(Na+1)] = expaninsideinfor[2];
				free(expaninsideinfor);
			}
			else if(x[i] > righttx && x[i] <= righthx)
			{
				double*expaninsideinfor = expaninsideinformake(infor, t, x[i], 2); 
				out->data[i] = expaninsideinfor[0];
				out->data[i+Na+1] = expaninsideinfor[1];
				out->data[i+2*(Na+1)] = expaninsideinfor[2];
				free(expaninsideinfor);
			}
		}

	}
	

	return out;
}
