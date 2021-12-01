#pragma once
#include "common.h"
#include "LinearAlgebra.h"
#include "Grid.h"

class VectorOfFunctions
{
public:
	virtual Vector ComputeInPoint(Vector point) = 0;
};

class DisjointCircles : public VectorOfFunctions
{
private:
	int size = 2;

public:
	Vector ComputeInPoint(Vector point) override;
};

class IntersectingCirclesAtPoint : public VectorOfFunctions
{
private:
	int size = 2;

public:
	Vector ComputeInPoint(Vector point) override;
};

class IntersectingCircles : public VectorOfFunctions
{
private:
	int size = 2;

public:
	Vector ComputeInPoint(Vector point) override;
};

class IntersectionOfCirclesWithLine : public VectorOfFunctions
{
private:
	int size = 3;

public:
	Vector ComputeInPoint(Vector point) override;
};

class IntersectingLines : public VectorOfFunctions
{
private:
	int size = 3;

public:
	Vector ComputeInPoint(Vector point) override;
};

class IntersectingLines_WithWeight : public VectorOfFunctions
{
private:
	int size = 3;

public:
	Vector ComputeInPoint(Vector point) override;
};

class IntersectingLineWithSineWave : public VectorOfFunctions
{
private:
	int size = 2;

public:
	Vector ComputeInPoint(Vector point) override;
};

struct SystemParameters
{
	int n;
	int m;
	int maxiter;
	int maxiterBeta;
	real epsF;
	real epsBeta;
	VectorOfFunctions* F;
	Vector x0 = Vector(1);

	SystemParameters(int n, int m, int maxiter, int maxiterBeta, real epsF, real epsBeta, VectorOfFunctions* Function, Vector x0);
};

class Squaring
{
public:
	virtual void LeadToSquare(Matrix& matrix, Vector& vector) = 0;
};

class Symmetrization : public Squaring
{
public:
	void LeadToSquare(Matrix& matrix, Vector& vector) override;
};

class ExcludingRows : public Squaring
{
public:
	void LeadToSquare(Matrix& matrix, Vector& vector) override;
};

class Convolution : public Squaring
{
public:
	void LeadToSquare(Matrix& matrix, Vector& vector) override;
};

class SystemOfNonlinearEquations
{
private:
	int n;
	int m;
	int maxiter;
	int maxiterBeta;
	real epsF;
	real epsBeta;
	VectorOfFunctions* F;
	Vector x0 = Vector(1);

	Squaring* squaring;

public:
	SystemOfNonlinearEquations(struct SystemParameters parameters, Squaring* squaring);
	Vector Solve();
	Grid SearchProcess();
	

private:
	Matrix FormJacobiMatrix(Vector xk);
	Matrix AnaliticalJacobiMatrix(Vector xk);
	Vector ComputeDirectionByGauss(Matrix matrix, Vector vector);
	Vector ComputeXk(Vector xk, Vector dx);
};