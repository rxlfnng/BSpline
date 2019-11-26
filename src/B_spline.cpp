#include <math.h>
#include <malloc.h>
#include <stdio.h>
#include <memory.h>
#include "BSpline.h"

#include <vector>

void SVDslove(double **A, int m, int n, double **inverse)
{
	double eps = 0.000001f;
	double** V = (double**)malloc(m * sizeof(double*));
	double** M = (double**)malloc(n * sizeof(double*));
	double** U = (double**)malloc(n * sizeof(double*));
	double** tmp = (double**)malloc(n * sizeof(double*));
	double** eigen = (double**)malloc(n * sizeof(double*));
	double** p = (double**)malloc(n * sizeof(double*));
	double* sum = (double*)malloc(n * sizeof(double));
	int i, j, k, i_t, j_t, iter = 0;
	double max = 9999, tan2a, sina, cosa, alpha;
	for (i = 0; i < m; i++)
	{
		V[i] = (double*)malloc(n * sizeof(double));
	}
	for (i = 0; i < n; i++)
	{
		M[i] = (double*)malloc(n * sizeof(double));
		U[i] = (double*)malloc(n * sizeof(double));
		tmp[i] = (double*)malloc(n * sizeof(double));
		eigen[i] = (double*)malloc(n * sizeof(double));
		p[i] = (double*)malloc(n * sizeof(double));
	}





	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			M[i][j] = 0;
			eigen[i][j] = 0;
			for (k = 0; k < m; k++)
			{
				double a = A[k][i];
				double b = A[k][j];
				M[i][j] += A[k][i] * A[k][j];
				eigen[i][j] += A[k][i] * A[k][j];
			}
		}
	}

	iter = 0;
	while (max > eps)
	{
		max = eigen[0][1];
		i_t = 0;
		j_t = 1;
		for (i = 0; i < n; i++)
		{
			for (j = i; j < n; j++)
			{

				if (j > i&&fabs(eigen[i][j]) > max)
				{
					i_t = i;
					j_t = j;
					max = fabs(eigen[i][j]);
				}
			}
		}
		tan2a = 2 * eigen[i_t][j_t] / (eigen[i_t][i_t] - eigen[j_t][j_t]);
		alpha = atan(tan2a) / 2;
		sina = sin(alpha);
		cosa = cos(alpha);
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				p[i][j] = 0;

			}
			p[i][i] = 1;
		}
		p[i_t][i_t] = cosa;
		p[j_t][j_t] = cosa;
		p[i_t][j_t] = -sina;
		p[j_t][i_t] = sina;

		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				tmp[i][j] = 0;
				for (k = 0; k < n; k++)
				{
					tmp[i][j] += p[k][i] * eigen[k][j];
				}
			}
		}
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				eigen[i][j] = 0;
				for (k = 0; k < n; k++)
				{
					eigen[i][j] += tmp[i][k] * p[k][j];
				}
			}
		}
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				if (iter == 0)
				{
					if (i == j)
					{
						tmp[i][j] = 1;
					}
					else
					{
						tmp[i][j] = 0;
					}


				}
				else
				{
					tmp[i][j] = U[i][j];
				}

			}
		}
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				U[i][j] = 0;
				for (k = 0; k < n; k++)
				{
					U[i][j] += tmp[i][k] * p[k][j];
				}
			}
		}
		iter++;
	}


	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (i == j)
			{
				if (eigen[i][j] != 0)
				{
					eigen[i][j] = 1 / sqrt(eigen[i][j]);
				}
				else
				{
					eigen[i][j] = 0;
				}

			}
			else
			{
				eigen[i][j] = 0;
			}
		}
	}

	for (i = 0; i < n; i++)
	{
		sum[i] = 0;
		for (j = 0; j < n; j++)
		{
			sum[i] += U[j][i] * U[j][i];
		}
	}
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (sum[j] != 0)
			{
				U[i][j] /= sqrt(sum[j]);
			}
			else
			{
				U[i][j] = 0;
			}
		}
	}
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			V[i][j] = 0;
			for (k = 0; k < n; k++)
			{
				V[i][j] += A[i][k] * U[k][j];
			}

		}
	}
	for (i = 0; i < n; i++)
	{
		sum[i] = 0;
		for (j = 0; j < m; j++)
		{
			sum[i] += V[j][i] * V[j][i];
		}
	}
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (sum[j] != 0)
			{
				V[i][j] /= sqrt(sum[j]);
			}
			else
			{
				V[i][j] = 0;
			}
		}
	}

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			tmp[i][j] = 0;
			for (k = 0; k < n; k++)
			{
				tmp[i][j] += U[i][k] * eigen[k][j];
			}
		}
	}


	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{
			inverse[i][j] = 0;
			for (k = 0; k < n; k++)
			{
				inverse[i][j] += tmp[i][k] * V[j][k];
			}
		}
	}


	free(sum);

	for (i = 0; i < m; i++)
	{
		free(V[i]);
	}
	for (i = 0; i < n; i++)
	{
		free(M[i]);
		free(U[i]);
		free(tmp[i]);
		free(eigen[i]);
		free(p[i]);
	}
	free(V);
	free(M);
	free(U);
	free(tmp);
	free(eigen);
	free(p);
}


void Decompose_LU(double **M, int dim, double** L, double **U)
{
	int i, j, k;
	double sum;
	for (i = 0; i < dim; ++i)
	{
		L[i][i] = 1;
	}
	for (i = 0; i < dim; ++i)
	{
		for (j = i; j < dim; ++j)
		{
			sum = 0;
			for (k = 0; k < i; ++k)
			{
				sum += L[i][k] * U[k][j];
			}
			U[i][j] = (M[i][j] - sum) / L[i][i];
		}

		for (j = i + 1; j < dim; ++j)
		{
			sum = 0;
			for (k = 0; k < i; ++k)
			{
				sum += L[j][k] * U[k][i];
			}
			L[j][i] = (M[j][i] - sum) / U[i][i];
		}
	}
	sum = 0;
	for (i = 0; i < dim - 1; i++)
	{
		sum += L[dim - 1][i] * U[i][dim - 1];
	}
	U[dim - 1][dim - 1] = M[dim - 1][dim - 1] - sum;


	return;

}
void LU_inverse(double **data, int dim, double** inverse)
{
	int i, j, k;
	double **L = (double**)malloc(dim * sizeof(double*));
	double **U = (double**)malloc(dim * sizeof(double*));
	double **L_inverse = (double**)malloc(dim * sizeof(double*));
	double **U_inverse = (double**)malloc(dim * sizeof(double*));
	for (i = 0; i < dim; ++i)
	{
		L[i] = (double*)calloc(dim, sizeof(double));
		U[i] = (double*)calloc(dim, sizeof(double));

		L_inverse[i] = (double*)malloc(dim * sizeof(double));
		U_inverse[i] = (double*)malloc(dim * sizeof(double));
	}

	Decompose_LU(data, dim, L, U);


	for (i = 0; i < dim; ++i)
	{
		for (j = 0; j < dim; ++j)
		{
			if (i == j)
			{
				L_inverse[i][i] = 1.0 / L[i][i];
				U_inverse[i][i] = 1.0 / U[i][i];
			}
			else if (i > j)
			{
				U_inverse[i][j] = 0;
			}
			else
			{
				L_inverse[i][j] = 0;
			}



		}
	}
	for (i = 0; i < dim - 1; ++i)
	{
		for (j = i + 1; j < dim; ++j)
		{
			if (j > i)
			{
				U_inverse[i][j] = 0;
				for (k = i; k < j; ++k)
				{
					U_inverse[i][j] += U_inverse[i][k] * U[k][j];
				}
				U_inverse[i][j] *= -U_inverse[j][j];
			}
		}
	}

	for (i = dim - 1; i >= 1; --i)
	{
		for (j = i - 1; j >= 0; --j)
		{
			if (i > j)
			{
				L_inverse[i][j] = 0;
				for (k = j + 1; k <= i; ++k)
				{
					L_inverse[i][j] += L_inverse[i][k] * L[k][j];
				}
				L_inverse[i][j] *= -L_inverse[j][j];



			}
		}
	}

	for (i = 0; i < dim; ++i)
	{
		for (j = 0; j < dim; ++j)
		{
			inverse[i][j] = 0;
			for (k = 0; k < dim; ++k)
			{
				inverse[i][j] += U_inverse[i][k] * L_inverse[k][j];
			}
		}
	}

	for (i = 0; i < dim; ++i)
	{
		free(L[i]);
		free(U[i]);

		free(L_inverse[i]);
		free(U_inverse[i]);
	}
	free(L);
	free(U);

	free(L_inverse);
	free(U_inverse);
}


double caculateDistancef(double* a, double* b, int size)
{
	double dis = 0;
	for (int i = 0; i < size; ++i)
	{
		dis += (b[i] - a[i])*(b[i] - a[i]);
	}
	return sqrt(dis);
}


BSpline::BSpline(Vect* InputData,int NumInput, int Dim, int Degree)
{
	m_NumInput = NumInput;
	m_Dim = Dim;
	m_Degree = Degree;
	m_InputData = InputData;
	SplineParameterization();
	m_KnotData = NULL;//节点
	m_ControlData = NULL;//控制点
	m_NumControl = -1;
	m_NumKnot = -1;
}

BSpline::~BSpline()
{
	if (m_ControlData)
	{
		free(m_ControlData);
	}
	if (m_KnotData)
	{
		free(m_KnotData);
	}
	if (m_ParamInputData)
	{
		free(m_ParamInputData);
	}
}

/*
n=m-p-1

*/
int BSpline::FindSpan(double u)
{
	int mid;
	
	if (u == m_KnotData[m_NumControl])
	{
		return m_NumControl - 1;
	}
	if (u == m_KnotData[0])
	{
		return m_Degree;
	}
	for (mid = m_Degree; mid < m_NumControl; ++mid)
	{
		if (u <= m_KnotData[mid + 1] && u > m_KnotData[mid])
		{
			break;
		}
	}
	//low = m_Degree;
	//high = m_NumControl;
	//mid = (low + high)/2;
	//while (u < m_KnotData[mid]||u >= m_KnotData[mid+1])
	//{
	//	if (u < m_KnotData[mid])
	//	{
	//		high = mid;
	//	} 
	//	else
	//	{
	//		low = mid;
	//	}
	//	mid = (low + high)/2;
	//}
	return mid;
}

void BSpline::BasisFuns(int i, double u, double* N)
{
	N[0] = 1.0;
	double saved;
	double* left = (double*)malloc((m_Degree + 1) * sizeof(double));
	double* right = (double*)malloc((m_Degree + 1) * sizeof(double));
	for (int j = 1; j <= m_Degree; ++j)
	{
		left[j] = u - m_KnotData[i + 1 - j];
		right[j] = m_KnotData[i + j] - u;
		saved = 0.0;
		for (int r = 0; r < j; ++r)
		{
			double temp = N[r] / (right[r + 1] + left[j - r]);
			N[r] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		N[j] = saved;
	}
	free(left);
	free(right);
}



////弦长参数化
//double BSpline::SplineParameterization(Vect* data, int n, int p, double *uk)
//{
//	double d = 0;
//	double *dist = (double*)malloc(n * sizeof(double));
//	uk[0] = 0;
//	dist[0] = 0;
//	uk[n - 1] = 1.0;
//	for (int i = 1; i < n; ++i)
//	{
//		dist[i] = caculateDistancef(data[i], data[i - 1], 2);
//		d += dist[i];
//	}
//	for (int i = 1; i < n - 1; ++i)
//	{
//		uk[i] = uk[i - 1] + dist[i]/d;
//	}
//
//	free(dist);
//	return d;
//}
////向心参数化
//void BSpline::SplineParameterization(Vect* data, int n, int p, double *uk)
//{
//	double d = 0;
//	double *dist = (double*)malloc(n * sizeof(double));
//	uk[0] = 0;
//	dist[0] = 0;
//	uk[n - 1] = 1.0;
//	for (int i = 1; i < n; ++i)
//	{
//		dist[i] = sqrt(caculateDistancef(data[i], data[i - 1], 2));
//		d += dist[i];
//	}
//	for (int i = 1; i < n - 1; ++i)
//	{
//		uk[i] = uk[i - 1] + dist[i] / d;
//	}
//
//
//
//	
//	for (int i = 0; i < n; ++i)
//	{
//		uk[i] = 1.0*i / (double)(n-1);
//	}
//	free(dist);
//}

//均匀参数化
void BSpline::SplineParameterization()
{
	m_ParamInputData = (double*)malloc(m_NumInput * sizeof(double));
	for (int i = 0; i < m_NumInput; ++i)
	{
		m_ParamInputData[i] = 1.0*i / (double)(m_NumInput - 1);
	}

}
//计算节点矢量
int BSpline::calculateKnotVector()
{
	int i, j;
	if (m_NumControl < m_Degree + 1)
	{
		return -1;
	}
	for (i = 0; i <= m_Degree; ++i)
	{
		m_KnotData[i] = 0.0;
	}
	for (i = m_NumControl; i < m_NumKnot; ++i)
	{
		m_KnotData[i] = 1.0;
	}
	for (i = 1; i < m_NumControl - m_Degree; ++i)
	{
		int id = i + m_Degree;
		double temp = 0.0;
		for (j = i; j <= i + m_Degree - 1; ++j)
		{
			temp += m_ParamInputData[j];
		}
		m_KnotData[id] = temp / m_Degree;
	}
	return 0;
}



void BSpline::GlobalCurveInterp()
{
	int i, j , k;

	
	double **A = (double**)calloc(m_NumControl ,  sizeof(double*));
	double **inverse = (double**)malloc(m_NumControl * sizeof(double*));
	for (i = 0; i < m_NumControl; ++i)
	{
		A[i] = (double*)calloc(m_NumControl, sizeof(double));
		inverse[i] = (double*)malloc(m_NumControl * sizeof(double));
	}
	
	for (i = 0; i < m_NumControl; ++i)
	{
		int span = FindSpan(m_ParamInputData[i]);
		BasisFuns(span, m_ParamInputData[i], A[i]+span - m_Degree);
	}
	LU_inverse(A, m_NumControl, inverse);
	//SVDslove(A, n, n, inverse);
	for (i = 0; i < m_NumControl; ++i)
	{
		for (j = 0; j < m_Dim; ++j)
		{
			m_ControlData[i][j] = 0;
			for (k = 0; k < m_NumControl; ++k)
			{
				m_ControlData[i][j] += inverse[i][k] * m_InputData[k][j];
			}
			
		}
	}
	
	for (i = 0; i < m_NumControl; ++i)
	{
		free(A[i]);
		free(inverse[i]);
	}
	free(A);
	free(inverse);
}

void BSpline::CurvePoint(double s, Vect *Point)
{
	double *N = (double*)malloc((m_Degree + 1) * sizeof(double));
	int span = FindSpan(s);
	BasisFuns(span, s, N);
	for (int j = 0; j < m_Dim; j++)
	{
		(*Point)[j] = 0;
		for (int i = 0; i <= m_Degree; ++i)
		{
			(*Point)[j] += N[i] * m_ControlData[span - m_Degree + i][j];
		}
	}
	free(N);
	
}
void BSpline::BSpline_CurveInterp(int NumOut, Vect *out)
{
	m_NumControl = m_NumInput;
	m_NumKnot = m_NumControl + m_Degree + 1;
	m_ControlData = (Vect*)malloc(m_NumControl * sizeof(Vect));
	m_KnotData = (double*)malloc(m_NumKnot * sizeof(double));
	
	calculateKnotVector();
	GlobalCurveInterp();
	
	double delta = 1.0 / (double)NumOut;
	for (int i = 0; i < NumOut; ++i)
	{
		double s = delta * i;
		Vect Point;
		CurvePoint(s, &Point);
		memcpy(out[i] , Point, sizeof(Vect));
		printf("[%f, %f, %f],\n", Point[0], Point[1], Point[2]);
	}
	free(m_ControlData);
	m_ControlData = NULL;
	free(m_KnotData);
	m_KnotData = NULL;
	m_NumControl = -1;
	m_NumKnot = -1;

}


//int BSpline::CurveFitKnotVector(double *uk, int n, int p, double *u)
//{
//	int i;
//	int m = n + p + 1;
//	if (n < p + 1)
//	{
//		return -1;
//	}
//
//	double d = (double)m / (double)(n - p);
//
//	for (i = 0; i <= p; ++i)
//	{
//		u[i] = 0.0;
//	}
//	for (i = n; i < m; ++i)
//	{
//		u[i] = 1.0;
//	}
//	for (i = 1; i < n - p; ++i)
//	{
//		int id = i + p;
//		int j = (int)(i*d);
//		double alpha = i * d - j;
//		u[id] = (1-alpha)*uk[j-1] + alpha*uk[j];
//	}
//	return 0;
//}
int BSpline::CurveFitKnotVector()
{
	int i;
	

	for (i = 0; i <= m_Degree; ++i)
	{
		m_KnotData[i] = 0.0;
	}
	for (i = m_NumControl; i < m_NumKnot; ++i)
	{
		m_KnotData[i] = 1.0;
	}
	for (i = 1; i < m_NumControl - m_Degree; ++i)
	{
		int id = i + m_Degree;
		m_KnotData[id] = 1.0*i / (m_NumControl - m_Degree);
	}
	return 0;
}
void BSpline::GlobalCurveFit()
{
	int i, j, k;
	Vect* Rk = (Vect*)malloc(m_NumInput * sizeof(Vect));
	Vect* R = (Vect*)malloc(m_NumControl * sizeof(Vect));
	double** N = (double**)calloc(m_NumInput, sizeof(double*));
	double** NTN = (double**)calloc(m_NumControl, sizeof(double*));
	double** inverse = (double**)malloc(m_NumControl * sizeof(double*));
	for (i = 0; i < m_NumInput; ++i)
	{
		N[i] = (double*)calloc(m_NumControl, sizeof(double));

	}
	for (i = 0; i < m_NumControl; ++i)
	{
		NTN[i] = (double*)calloc(m_NumControl, sizeof(double));
		inverse[i] = (double*)malloc(m_NumControl * sizeof(double));
	}
	
	for (i = 0; i < m_NumInput; ++i)
	{
		int span = FindSpan(m_ParamInputData[i]);
		BasisFuns(span, m_ParamInputData[i], N[i] + span - m_Degree);
	}
	for (i = 0; i < m_NumControl; ++i)
	{
		for (k = 0; k < m_Dim; ++k)
		{
			R[i][k] = 0;
			for (j = 0; j < m_NumInput; ++j)
			{

				double Nip = N[j][i];
				//printf("%f\n", N[j][i]);
				R[i][k] += Nip * m_InputData[j][k];
			}
		}
	}
	
	free(Rk);

	for (i = 0; i < m_NumControl; ++i)
	{
		for (j = 0; j < m_NumControl; ++j)
		{
			NTN[i][j] = 0;
			for (k = 0; k < m_NumInput; ++k)
			{
				NTN[i][j] += N[k][i] * N[k][j];
			}
		}
	}


	for (i = 0; i < m_NumInput; ++i)
	{
		free(N[i]);

	}
	free(N);

	LU_inverse(NTN, m_NumControl, inverse);


	for (i = 0; i < m_NumControl; ++i)
	{
		free(NTN[i]);

	}
	free(NTN);
	
	for (i = 0; i < m_NumControl; ++i)
	{
		for (j = 0; j < m_Dim; ++j)
		{
			m_ControlData[i][j] = 0;
			for (k = 0; k < m_NumControl; ++k)
			{
				m_ControlData[i][j] += inverse[i][k] * R[k][j];

			}

		}
		printf("P = [%lf, %lf, %lf],\n", m_ControlData[i][0], m_ControlData[i][1], m_ControlData[i][2]);
	}
	
	for (i = 0; i < m_NumControl; ++i)
	{
		free(inverse[i]);

	}
	free(inverse);

	free(R);
}
int BSpline::selectControlPoints(double eps)
{
	int i;
	bool flag = false;
	if (m_NumInput - m_Degree - 1 < 2)
	{
		return -1;
	}
	for (m_NumControl = m_NumInput - m_Degree - 1; m_NumControl > m_Degree + 1; --m_NumControl)
	{
		m_NumKnot = m_NumControl + m_Degree + 1;
		m_ControlData = (Vect*)malloc(m_NumControl * sizeof(Vect));
		m_KnotData = (double*)malloc(m_NumKnot * sizeof(double));
		
		CurveFitKnotVector();
		GlobalCurveFit();

		for (i = 0; i < m_NumInput; ++i)
		{
			Vect Point;

			CurvePoint(m_ParamInputData[i], &Point);
			double error = caculateDistancef(Point, m_InputData[i], 2);
			if (isnormal(error)&&error > eps)
			{
				flag = true;
				break;
			}
		}
		free(m_ControlData);
		m_ControlData = NULL;
		free(m_KnotData);
		m_KnotData = NULL;
		if (flag)
		{
			break;
		}
	}
	m_NumControl = m_NumControl;
	m_NumKnot = m_NumControl + m_Degree + 1;


	return 0;
}
int BSpline::BSpline_CurveFit(int NumControl, int NumOut, Vect *OutData)
{
	m_NumControl = NumControl;
	m_NumKnot = m_NumControl + m_Degree + 1;
	if (m_NumInput <= m_NumKnot)
	{
		return -1;
	}
	m_ControlData = (Vect*)malloc(m_NumControl * sizeof(Vect));
	m_KnotData = (double*)malloc(m_NumKnot * sizeof(double));
	CurveFitKnotVector();
	GlobalCurveFit();

	double delta = 1.0 / (double)NumOut;
	for (int i = 0; i < NumOut; ++i)
	{
		double s = delta * i;
		Vect Point;
		CurvePoint(s, &Point);
		memcpy(OutData[i], Point, sizeof(Vect));
		printf("[%lf, %lf, %lf],\n", Point[0], Point[1], Point[2]);
	}
	free(m_ControlData);
	m_ControlData = NULL;
	free(m_KnotData);
	m_KnotData = NULL;
	m_NumControl = -1;
	m_NumKnot = -1;
	return 0;
}

int BSpline::BSpline_CurveFit(double eps, int NumOut, Vect *OutData)
{
	int flag = selectControlPoints(eps);
	if (flag)
	{
		printf("input data is too less!!!");
		return -1;
	}
	m_ControlData = (Vect*)malloc(m_NumControl * sizeof(Vect));
	m_KnotData = (double*)malloc(m_NumKnot * sizeof(double));
	CurveFitKnotVector();
	GlobalCurveFit();
	double delta = 1.0 / (double)NumOut;
	for (int i = 0; i < NumOut; ++i)
	{
		double s = delta * i;
		Vect Point;
		CurvePoint(s, &Point);
		memcpy(OutData[i], Point, sizeof(Vect));
		printf("[%lf, %lf, %lf],\n", Point[0], Point[1], Point[2]);
	}
	free(m_ControlData);
	m_ControlData = NULL;
	free(m_KnotData);
	m_KnotData = NULL;
	m_NumControl = -1;
	m_NumKnot = -1;
	return 0;
}

CubicSpline::CubicSpline(Vect* InputData, int NumInput, int Dim, double lambda)
{
	int i, j;
	m_NumInput = NumInput;
	m_Dim = Dim;
	m_Lambda = lambda;
	m_pInputData = (Vect*)malloc(m_NumInput * sizeof(Vect));
	memcpy(m_pInputData, InputData, m_NumInput * sizeof(Vect));
	m_pStepData = (double*)malloc((m_NumInput - 1) * sizeof(double));
	double *p = (double*)malloc((m_NumInput - 1) * sizeof(double));
	Vect* q = (Vect*)malloc(m_NumInput * sizeof(Vect));
	a = (Vect*)malloc((m_NumInput - 1) * sizeof(Vect));
	b = (Vect*)malloc((m_NumInput - 1) * sizeof(Vect));
	c = (Vect*)malloc((m_NumInput - 1) * sizeof(Vect));
	d = (Vect*)malloc((m_NumInput - 1) * sizeof(Vect));

	SplineParameterization();
	
	
	
	//m_pStepData[0] = m_pParamInputData[1] - m_pParamInputData[0];//caculateDistancef(m_pInputData[1], m_pInputData[0], dd);
	//for (i = 1; i < m_NumInput - 1; ++i)
	//{
	//	m_pStepData[i] = m_pParamInputData[i + 1] - m_pParamInputData[i]; //caculateDistancef(m_pInputData[i + 1], m_pInputData[i], dd);
	//	p[i] = 2 * (m_pParamInputData[i + 1] - m_pParamInputData[i - 1]);
	//	for (j = 0; j < m_Dim; ++j)
	//	{
	//		q[i][j] = 3 * (m_pInputData[i + 1][j] - m_pInputData[i][j]) / m_pStepData[i] - 3 * (m_pInputData[i][j] - m_pInputData[i - 1][j]) / m_pStepData[i - 1];
	//	}
	//}

	///*Gaussian Elimination*/
	//for (i = 2; i < m_NumInput - 1; ++i)
	//{
	//	p[i] = p[i] - m_pStepData[i - 1] * m_pStepData[i - 1] / p[i - 1];
	//	for (j = 0; j < m_Dim; ++j)
	//	{
	//		q[i][j] = q[i][j] - q[i - 1][j] * m_pStepData[i - 1] / p[i - 1];
	//	}
	//}
	///*Backsubtitution*/
	//for (j = 0; j < m_Dim; ++j)
	//{
	//	b[m_NumInput - 2][j] = q[m_NumInput - 2][j] / p[m_NumInput - 2];
	//}
	//for (i = 3; i < m_NumInput; ++i)
	//{
	//	for (j = 0; j < m_Dim; ++j)
	//	{
	//		b[m_NumInput - i][j] = (q[m_NumInput - i][j] - m_pStepData[m_NumInput - i] * b[m_NumInput - i + 1][j]) / p[m_NumInput - i];
	//	}
	//}

	///*spline parameters*/
	//double ddy = (3.0 * m_pStepData[0]);
	//for (j = 0; j < m_Dim; ++j)
	//{
	//	a[0][j] = b[1][j] / ddy;
	//	b[0][j] = 0;
	//	c[0][j] = (m_pInputData[1][j] - m_pInputData[0][j]) / m_pStepData[0] - b[1][j] * m_pStepData[0] / 3;
	//	d[0][j] = m_pInputData[0][j];
	//	//b[m_NumInput - 2][j] = 0;
	//}
	//printf("[%lf, %lf, %lf, %lf],\n", b[1][0], b[1][1], b[1][2], m_pStepData[0]);
	//printf("[%lf, %lf, %lf],\n", a[0][0], a[0][1], a[0][2]);
	//for (i = 1; i < m_NumInput - 1; ++i)
	//{
	//	for (j = 0; j < m_Dim; ++j)
	//	{
	//		a[i][j] = (b[i + 1][j] - b[i][j]) / (2 * m_pStepData[i]);
	//		c[i][j] = (m_pInputData[i + 1][j] - m_pInputData[i][j]) / m_pStepData[i] - a[i][j] * m_pStepData[i] * m_pStepData[i] - b[i][j] * m_pStepData[i];//(b[i][j] + b[i - 1][j])*m_pStepData[i - 1] + c[i - 1][j];
	//		d[i][j] = m_pInputData[i][j];
	//	}
	//}


	double *r = (double*)malloc((m_NumInput - 1) * sizeof(double));
	double *f = (double*)malloc((m_NumInput - 1) * sizeof(double));
	double *u = (double*)malloc((m_NumInput - 1) * sizeof(double));
	double *v = (double*)malloc((m_NumInput - 1) * sizeof(double));
	double *w = (double*)malloc((m_NumInput - 1) * sizeof(double));
	double sigma = 1.0;
	double mu = 2 * (1 - m_Lambda) / (3 * m_Lambda);
	m_pStepData[0] = m_pParamInputData[1] - m_pParamInputData[0];
	r[0] = 3 / m_pStepData[0];
	for (i = 1; i < m_NumInput - 1; ++i)
	{
		m_pStepData[i] = m_pParamInputData[i + 1] - m_pParamInputData[i];
		r[i] = 3 / m_pStepData[i];
		f[i] = -(r[i - 1] + r[i]);
		p[i] = 2 * (m_pParamInputData[i + 1] - m_pParamInputData[i - 1]);
		for (j = 0; j < m_Dim; ++j)
		{
			q[i][j] = 3 * (m_pInputData[i + 1][j] - m_pInputData[i][j]) / m_pStepData[i] - 3 * (m_pInputData[i][j] - m_pInputData[i - 1][j]) / m_pStepData[i - 1];
		}
	}
	u[0] = 0.0;
	v[0] = 0.0;
	w[0] = 0.0;
	for (i = 1; i < m_NumInput - 1; ++i)
	{
		u[i] = r[i - 1]*r[i - 1]*sigma + f[i]* f[i] *sigma + r[i] * r[i] *sigma;
		u[i] = mu * u[i] + p[i];
		v[i] = f[i] * r[i] * sigma + r[i] * f[i + 1] * sigma;
		v[i] = mu * v[i] + m_pStepData[i];
		w[i] = mu * r[i] * r[i + 1] * sigma;
		
	}

	Quincunx(m_NumInput, u, v, w, q);

	/*spline parameters*/
	for (j = 0; j < m_Dim; ++j)
	{
		d[0][j] = m_pInputData[0][j] - mu * r[0] * q[1][j] * sigma;
		double dx = m_pInputData[1][j] - mu * (f[1] * q[1][j] + r[1] * q[2][j])*sigma;
		
		a[0][j] = q[1][j] / (3.0 * m_pStepData[0]);
		b[0][j] = 0;
		c[0][j] = (dx - d[0][j]) / m_pStepData[0] - q[1][j] * m_pStepData[0] / 3;
		//r[0] = 0;
	}
	for (i = 1; i < m_NumInput - 1; ++i)
	{
		for (j = 0; j < m_Dim; ++j)
		{
			a[i][j] = (q[i + 1][j] - q[i][j]) / (3 * m_pStepData[i]);
			b[i][j] = q[i][j];
			c[i][j] = (q[i][j] + q[i - 1][j])*m_pStepData[i - 1] + c[i - 1][j];
			d[i][j] = r[i - 1] * q[i - 1][j] + f[i] * q[i][j] + r[i]*q[i + 1][j];
			d[i][j] = m_pInputData[i][j] - mu * d[i][j] * sigma;

		}
	}

	free(r);
	free(f);
	free(u);
	free(v);
	free(w);
	free(p);
	free(q);

}
CubicSpline::~CubicSpline()
{
	free(m_pParamInputData);
	free(m_pStepData);
	free(m_pInputData);
	free(a);
	free(b);
	free(c);
	free(d);
}
void CubicSpline::Quincunx(int n, double* u, double *v, double *w, Vect* q)
{
	int i,j;
	u[0] = 0;
	/*factorisation*/
	//u[1] = u[1] - u[0] * v[0] * v[0];
	v[1] = v[1] / u[1];
	w[1] = w[1] / u[1];
	for (i = 2; i < n - 1; ++i)
	{
		u[i] = u[i] - u[i - 2] * w[i - 2] * w[i - 2] - u[i - 1] * v[i - 1] * v[i - 1];
		v[i] = (v[i] - u[i - 1] * v[i - 1] * w[i - 1]) / u[i];
		w[i] = w[i] / u[i];
	}

	/*forward substitution*/

	for (i = 2; i < n - 1; ++i)
	{
		for (j = 0; j < m_Dim; ++j)
		{
			q[i][j] = q[i][j] - v[i - 1] * q[i - 1][j] - w[i - 2] * q[i - 2][j];
		}
	}

	for (i = 1; i < n - 1; ++i)
	{
		for (j = 0; j < m_Dim; ++j)
		{
			q[i][j] = q[i][j] / u[i];
		}
	}

	/*back substitution*/
	for (j = 0; j < m_Dim; ++j)
	{
		q[0][j] = 0.0;
		q[n - 1][j] = 0;
	}
	for (i = n - 3; i > 0; --i)
	{
		for (j = 0; j < m_Dim; ++j)
		{

			q[i][j] = q[i][j] - v[i] * q[i + 1][j] - w[i] * q[i + 2][j];
		}
	}

}

int CubicSpline::FindSpan(double u)
{
	int mid;

	for (mid = 0; mid < m_NumInput - 1; ++mid)
	{
		if (u < m_pParamInputData[mid + 1] && u >= m_pParamInputData[mid])
		{
			break;
		}
	}
	return mid;
}

void CubicSpline::CulicSpline_CurveInterp(int NumOut, Vect *OutData)
{
	int i, j;
	double delta = 1.0 / NumOut;

	for (i = 0; i < NumOut; ++i)
	{
		double u = i * delta*m_ArcLen;
		int span = FindSpan(u);
		double sx = m_pParamInputData[span];
		double cx = u - sx;

		for (j = 0; j < m_Dim; ++j)
		{
			
			OutData[i][j] = a[span][j] * cx*cx*cx + b[span][j] * cx*cx + c[span][j] * cx + d[span][j];
		}

	}


}

//弦长参数化
void CubicSpline::SplineParameterization()
{
	double dist = 0;
	m_pParamInputData = (double*)malloc(m_NumInput * sizeof(double));
	m_pParamInputData[0] = 0;
	for (int i = 1; i < m_NumInput; ++i)
	{
		dist = caculateDistancef(m_pInputData[i], m_pInputData[i - 1], dd);
		m_pParamInputData[i] = m_pParamInputData[i - 1] + dist;
	}
	m_ArcLen = m_pParamInputData[m_NumInput - 1];

}
int main()
{
	//3D
	int dim = dd;
	double wx[23] = { 1623.382813, 1287.131104, 1209.996582, 1065.779297, 983.905090, 892.094543, 835.077087, 796.471741, 751.769348,713.332092,669.832581,601.487244,556.960999,475.868347,418.927307,344.341064, 300.111877, 253.929565, 200.157654, 161.352356, 106.556396, 49.665283, -6.967102 };
	double wy[23] = { 704.965027, 689.627258, 694.901245, 700.669373, 704.640991,707.915222, 709.012573, 703.724487, 698.676758,692.176392,690.623779,690.779785,690.075928,691.352539,695.278809,699.948914, 702.718506, 688.662476,677.834961, 669.927490, 661.460449, 654.651489, 651.501038 };
	double wz[23] = { -12.510620, 251.928284, 217.869019, 193.214111, 111.798584,11.280273, 8.179810, 123.177490, 215.058594,273.975037,304.385864,324.083374,296.066772,270.943787,198.506470,117.120605, 41.464844, 157.157593, 251.385193, 294.979492, 327.523560, 331.231323, 307.992493 };
	
	int NumInput = 23;
	int NumOut = 200;
	int NumControl = 18;
	double eps = 50.0;
	Vect* data = (Vect*)malloc(NumInput * sizeof(Vect));
	Vect *out = (Vect*)malloc(NumOut * sizeof(Vect));


	for (int i = 0; i < NumInput; ++i)
	{
		data[i][0] = wx[i];
		data[i][1] = wy[i];
		data[i][2] = wz[i];
	}
	//BSpline *bspline = new BSpline(data, NumInput, dim, 3);
	//int ret = bspline->BSpline_CurveFit(eps, NumOut, out);
	//bspline->BSpline_CurveFit(NumControl, NumOut, out);
	////bspline->BSpline_CurveInterp(NumOut, out);
	//delete(bspline);
	
	CubicSpline* cspline = new CubicSpline(data, NumInput, dim, 0.00001);
	cspline->CulicSpline_CurveInterp(NumOut, out);
	for (int i = 0; i < NumOut; ++i)
	{
		printf("[%f, %f, %f],\n", (float)(out[i][0]), (float)(out[i][1]), (float)(out[i][2]));
	}



	free(data);
	free(out);
	getchar();
	return 1;
	
}






