#pragma once
#include "SystemOfNonlinearEquations.h"
#include <iostream>

SystemParameters::SystemParameters(int n, int m, int maxiter, int maxiterBeta, real epsF, real epsBeta, VectorOfFunctions* Function, Vector x0)
{
	this->n = n;
	this->m = m;
	this->maxiter = maxiter;
	this->maxiterBeta = maxiterBeta;
	this->epsF = epsF;
	this->epsBeta = epsBeta;
	this->F = Function;
	this->x0 = x0;
}

Vector DisjointCircles::ComputeInPoint(Vector point)
{
	Vector value(size);

	value(0) = pow(point(0) + 2, 2) + pow(point(1) - 2, 2) - 4;
	value(1) = pow(point(0) - 2, 2) + pow(point(1) - 2, 2) - 4;

	return value;
}

Vector IntersectingCirclesAtPoint::ComputeInPoint(Vector point)
{
	Vector value(size);

	value(0) = pow(point(0) + 2, 2) + pow(point(1) - 2, 2) - 4;
	value(1) = pow(point(0) - 2, 2) + pow(point(1) - 2, 2) - 4;

	return value;
}

Vector IntersectingCircles::ComputeInPoint(Vector point)
{
	Vector value(size);

	value(0) = pow(point(0) + 1, 2) + pow(point(1) - 2, 2) - 4;
	value(1) = pow(point(0) - 1, 2) + pow(point(1) - 2, 2) - 4;

	return value;
}

Vector IntersectionOfCirclesWithLine::ComputeInPoint(Vector point)
{
	Vector value(size);

	value(0) = pow(point(0) + 2, 2) + pow(point(1) - 2, 2) - 4;
	value(1) = pow(point(0) - 2, 2) + pow(point(1) - 2, 2) - 4;
	value(2) = point(0) + point(1) - 2;

	return value;
}

Vector IntersectingLines::ComputeInPoint(Vector point)
{
	Vector value(size);

	value(0) = point(0);
	value(1) = point(1);
	value(2) = point(0) + point(1) - 4;

	return value;
}

Vector IntersectingLines_WithWeight::ComputeInPoint(Vector point)
{
	Vector value(size);

	value(0) = point(0);
	value(1) = 300 * (point(1));
	value(2) = (point(0) + point(1) - 4);

	return value;
}

Vector IntersectingLineWithSineWave::ComputeInPoint(Vector point)
{
	Vector value(size);

	value(0) = sin(point(0)) - point(1);
	value(1) = point(0) - point(1) + 1;

	return value;
}

void Symmetrization::LeadToSquare(Matrix& matrix, Vector& vector)
{
	Matrix temp = matrix.Transpose();
	matrix = temp * matrix;
	vector = temp * vector;
}

void ExcludingRows::LeadToSquare(Matrix& matrix, Vector& vector)
{
	int m = matrix.Rows();
	int n = matrix.Columns();
	int rowsToDelete = m - n;
	int row;
	real min;

	for (int i = 0; i < rowsToDelete; i++)
	{
		row = 0;
		min = fabs(vector(0));
		for (int j = 0; j < m - i; j++)
		{
			if (fabs(vector(j)) < min)
			{
				min = fabs(vector(j));
				row = j;
			}
		}

		swap(vector(row), vector(m - i - 1));
		for (int j = 0; j < n; j++)
		{
			swap(matrix(row, j), matrix(m - i - 1, j));
		}
	}
}

void Convolution::LeadToSquare(Matrix& matrix, Vector& vector)
{
	int m = matrix.Rows();
	int n = matrix.Columns();
	int rowsToDelete = m - n + 1;
	int row;
	real min;
	real sum = 0.0;
	Vector rowsSum(n);

	for (int i = 0; i < rowsToDelete; i++)
	{
		row = 0;
		min = fabs(vector(0));
		for (int j = 0; j < m - i; j++)
		{
			if (fabs(vector(j)) < min)
			{
				min = fabs(vector(j));
				row = j;
			}
		}

		sum += pow(vector(row), 2);
		for (int j = 0; j < n; j++)
		{
			rowsSum(j) += pow(matrix(row, j), 2);
		}

		swap(vector(row), vector(m - i - 1));
		for (int j = 0; j < n; j++)
		{
			swap(matrix(row, j), matrix(m - i - 1, j));
		}
	}

	vector(m - 1) = sum;
	for (int i = 0; i < n; i++)
	{
		matrix(m - 1, i) = 2 * rowsSum(i);
	}
}

SystemOfNonlinearEquations::SystemOfNonlinearEquations(struct SystemParameters parameters, Squaring* squaring)
{
	this->n = parameters.n;
	this->m = parameters.m;
	this->maxiter = parameters.maxiter;
	this->maxiterBeta = parameters.maxiterBeta;
	this->epsF = parameters.epsF;
	this->epsBeta = parameters.epsBeta;
	this->F = parameters.F;
	this->x0 = parameters.x0;

	this->squaring = squaring;
}

Matrix SystemOfNonlinearEquations::FormJacobiMatrix(Vector xk)
{
	Matrix Jacobi(m, n);
	Vector Fp = F->ComputeInPoint(xk);

	real h = 1e-10;

	Vector temp = xk;

	

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			temp(j) += h;
			Vector Fp_h = F->ComputeInPoint(temp);

			for (int k = 0; k < m; k++)
			{
				Fp_h(k) -= Fp(k);
			}

			Jacobi(i, j) = Fp_h(i) / h;
			temp(j) = xk(j);
		}
	}

	return Jacobi;
}

Vector SystemOfNonlinearEquations::ComputeDirectionByGauss(Matrix matrix, Vector vector)
{
	vector = -1 * vector;

	int i;

	for (i = 0; i < n; i++)
	{
		real mainElement = 0.0;
		int row = 0;

		for (int j = i; j < n; j++)
		{
			if (mainElement < fabs(matrix(j, i)))
			{
				mainElement = matrix(j, i);
				row = j;
			}
		}

		if (row != i)
		{
			swap(vector(i), vector(row));
			for (int j = 0; j < n; j++)
			{
				swap(matrix(i, j), matrix(row, j));
			}
		}

		vector(i) /= mainElement;
		for (int j = i + 1; j < n; j++)
		{
			matrix(i, j) /= mainElement;
		}

		for (int j = i + 1; j < n; j++)
		{
			mainElement = matrix(j, i);

			for (int k = i; k < n; k++)
			{
				matrix(j, k) -= mainElement * matrix(i, k);
			}

			vector(j) -= mainElement * vector(i);
		}
	}

	for (i -= 2; i >= 0; i--)
	{
		for (int j = i + 1; j < n; j++)
		{
			vector(i) -= vector(j) * matrix(i, j);
		}
	}
	
	return vector;
}

Vector SystemOfNonlinearEquations::ComputeXk(Vector xk, Vector dx)
{
	Vector xk1(n);
	real beta = 1.0;
	real FNorm = F->ComputeInPoint(xk).EuqlideanNorm();

	for (int v = 0; v < maxiterBeta && beta > epsBeta; v++)
	{
		xk1 = xk + beta * dx;

		real FkNorm = F->ComputeInPoint(xk1).EuqlideanNorm();
		if (FkNorm < FNorm)
		{
			break;
		}

		beta /= 2;
	}

	cout << "beta: " << beta << endl;

	return xk1;
}

Vector SystemOfNonlinearEquations::Solve()
{
	real F0Norm = F->ComputeInPoint(x0).EuqlideanNorm();
	Vector xk = x0;
	real discrepancy = 1.0;

	for (int k = 0; k < maxiter && discrepancy > epsF; k++)
	{
		Vector Fk = F->ComputeInPoint(xk);
		Matrix Jacobi = FormJacobiMatrix(xk);

		if (m != n)
		{
			squaring->LeadToSquare(Jacobi, Fk);
		}

		Vector dx = ComputeDirectionByGauss(Jacobi, Fk);

		xk = ComputeXk(xk, dx);

		discrepancy = F->ComputeInPoint(xk).EuqlideanNorm() / F0Norm;
	}

	return xk;
}

Grid SystemOfNonlinearEquations::SearchProcess()
{
	real F0Norm = F->ComputeInPoint(x0).EuqlideanNorm();
	Vector xk = x0;
	real discrepancy = 1.0;

	for (int k = 0; k < maxiter && discrepancy > epsF; k++)
	{
		Vector Fk = F->ComputeInPoint(xk);
		Matrix Jacobi = FormJacobiMatrix(xk);

		if (m != n)
		{
			squaring->LeadToSquare(Jacobi, Fk);
		}

		Vector dx = ComputeDirectionByGauss(Jacobi, Fk);

		xk = ComputeXk(xk, dx);

		discrepancy = F->ComputeInPoint(xk).EuqlideanNorm() / F0Norm;

		cout << "discrepancy: " << discrepancy << endl;

		co_yield xk;
	}
}

Matrix SystemOfNonlinearEquations::AnaliticalJacobiMatrix(Vector xk)
{
	Matrix Jacobi(m, n);

	Jacobi(0, 0) = cos(xk(0));
	Jacobi(0, 1) = -1;
	Jacobi(1, 0) = 1;
	Jacobi(1, 1) = -1;

	return Jacobi;
}

/*	Непересекающиеся окружности
	Jacobi(0, 0) = 2 * (xk(0) + 3);
	Jacobi(0, 1) = 2 * (xk(1) - 2);
	Jacobi(1, 0) = 2 * (xk(0) - 3);
	Jacobi(1, 1) = 2 * (xk(1) - 2);
	

	В одной точке
	Jacobi(0, 0) = 2 * (xk(0) + 2);
	Jacobi(0, 1) = 2 * (xk(1) - 2);
	Jacobi(1, 0) = 2 * (xk(0) - 2);
	Jacobi(1, 1) = 2 * (xk(1) - 2);


	В двух точках
	Jacobi(0, 0) = 2 * (xk(0) + 1);
	Jacobi(0, 1) = 2 * (xk(1) - 2);
	Jacobi(1, 0) = 2 * (xk(0) - 1);
	Jacobi(1, 1) = 2 * (xk(1) - 2);


	Окружности и прямая
	Jacobi(0, 0) = 2 * (xk(0) + 2);
	Jacobi(0, 1) = 2 * (xk(1) - 2);
	Jacobi(1, 0) = 2 * (xk(0) - 2);
	Jacobi(1, 1) = 2 * (xk(1) - 2);
	Jacobi(2, 0) = 1;
	Jacobi(2, 1) = 1;


	Три прямых
	Jacobi(0, 0) = 1;
	Jacobi(0, 1) = 0;
	Jacobi(1, 0) = 0;
	Jacobi(1, 1) = 1;
	Jacobi(2, 0) = 1;
	Jacobi(2, 1) = 1;


	Три прямых, одна взевешена
	Jacobi(0, 0) = 1;
	Jacobi(0, 1) = 0;
	Jacobi(1, 0) = 0;
	Jacobi(1, 1) = 300;
	Jacobi(2, 0) = 1;
	Jacobi(2, 1) = 1;


	Синусоида и прямая
	Jacobi(0, 0) = cos(xk(0));
	Jacobi(0, 1) = -1;
	Jacobi(1, 0) = 1;
	Jacobi(1, 1) = -1;
*/