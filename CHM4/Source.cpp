#include "LinearAlgebra.h"
#include "SystemOfNonlinearEquations.h"
#include <iostream>

int main(int argc, char* argv[])
{
	Vector x0(2);
	x0(0) = 5.0;
	x0(1) = 0.0;
	SystemOfNonlinearEquations s
	(
		SystemParameters(2, 2, 1000, 1000, 1e-7, 1e-7, new DisjointCircles(), x0),
		new Convolution()
	);



	//Vector result = s.Solve();

	Grid result = s.SearchProcess();

	for (Vector xk : result)
	{
		for (int i = 0; i < 2; i++)
		{
			cout << xk(i) << " ";
		}
		cout << endl;
	}

	return 0;
}