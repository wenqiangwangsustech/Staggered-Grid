/***************************************************************************************/
/******************************************2019-01-04***********************************/
/************************************Author:Wenqiang Wang*******************************/
/**********************Southern University of Science and Technology********************/
/***************************************************************************************/
#include <iostream>
#include "ElasticWaveEquation.h"

int main()
{
	GRID grid = { 100, 100, 100 };
	Delta delta = { 0.0005f, 5.0f, 5.0f, 5.0f };
	SOURCELOCATION sourceLocation = { 0 };	
	sourceLocation.x = grid.xLength / 2;
	sourceLocation.y = grid.yLength / 2;
	sourceLocation.z = grid.zLength / 2;
	int timeLength = 1000;
	float mainFrequncy = 30.0f;
	int order = 10;
	float * diff_coef = new float[order / 2];
	
	diff_coef[0] = 1.21124f;
	diff_coef[1] = -8.97217e-2f;
	diff_coef[2] = 1.38428e-2f;
	diff_coef[3] = -1.76566e-3f;
	diff_coef[4] = 1.18680e-4f;

	float * buoyancy;
	float * lambda;
	float * mu;

	ALLOCATE_3D_pragma(float, grid.xLength, grid.yLength, grid.zLength);

	for ( int i = 0; i < grid.xLength * grid.yLength * grid.zLength; i++)
	{
		buoyancy[i] = 1.0f / 2500.0f;
		lambda[i] = 2500.0f * pow( 2000.0f, 2);
		mu[i] = 2500.0f * pow(2000.0f, 2) / 3.0f;
	}

	ElasticWaveEquation ewq =
		ElasticWaveEquation( 
			timeLength,
			grid, delta,
			sourceLocation, mainFrequncy, diff_coef,
			buoyancy, lambda, mu,
			order);

	time_t start, stop;
	start = time(NULL);
	ewq.run();
	//ewq.calculateDamping();
	stop = time(NULL);

	cout << (stop - start) << endl;

	delete[] diff_coef;
	DELETE_3D_pragma(float);
//	system("pause");
}
