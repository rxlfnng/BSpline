#pragma once
#define dd 3
typedef double Vect[dd];

void SVDslove(double **A, int m, int n, double **inverse);
void Decompose_LU(double **M, int dim, double** L, double **U);
void LU_inverse(double **data, int dim, double** inverse);
double caculateDistancef(double* a, double* b, int size);

class BSpline
{
public:
	BSpline(Vect* InputData, int NumImput, int Dim, int Degree);
	~BSpline();
	int BSpline_CurveFit(int NumControl, int NumOut, Vect *OutData);
	int BSpline_CurveFit(double eps, int NumOut, Vect *OutData);//fit 
	void BSpline_CurveInterp(int NumOut, Vect *out);
private:

	int FindSpan(double u);
	//double OneBasisFun(int p, int m, double* U, int i, double u);
	void BasisFuns(int i, double u, double* N);
	//void DersBasisFuns(int i, double u, int p, int n, double* U, double** ders);
	//void DersOneBasisFun(int p, int m, double *U, int i, double u, int n, double* ders);
	
	void SplineParameterization();
	int calculateKnotVector();
	int CurveFitKnotVector();


	int selectControlPoints(double eps);
	

	void GlobalCurveFit();
	void GlobalCurveInterp();
	void CurvePoint(double s, Vect *Point);

	int m_Dim;
	int m_NumControl;
	int m_NumInput;
	int m_NumKnot;
	int m_Degree;


	double *m_KnotData;//节点
	double *m_ParamInputData;//参数化输入
	Vect *m_ControlData;//控制点
	Vect *m_InputData;
};
