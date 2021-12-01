#include "LinearAlgebra.h"
#include "SystemOfNonlinearEquations.h"
#include <iostream>

int main()
{
	Vector x0(2);
	x0(0) = 1.0;
	x0(1) = 0.0;
	SystemOfNonlinearEquations s
	(
		SystemParameters(2, 2, 1000, 1000, 1e-7, 1e-7, new DisjointCircles(), x0),
		new ExcludingRows()
	);

	cout << scientific;

	Grid result = s.SearchProcess();

	int j = 1;
	for (Vector xk : result)
	{
		
		for (int i = 0; i < 2; i++)
		{
			cout << xk(i) << " ";
		}

		cout << j << endl;

		cout << endl << endl;

		j++;
	}
}

/*
* auto c = new Symmetrization();
	auto f = new IntersectionOfCirclesWithLine();
	auto m = s.AnaliticalJacobiMatrix(x0);

	auto fk = f->ComputeInPoint(x0);
	c->LeadToSquare(m, fk);

	Vector result = s.ComputeDirectionByGauss(m, fk);

	

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			//cout << m(i, j) << " ";
		}
		cout << result(i) << " ";
		//cout << endl;
	}
* 
* 
* 
*
*/