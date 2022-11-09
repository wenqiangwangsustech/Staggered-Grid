/***************************************************************************************/
/******************************************2019-01-05***********************************/
/************************************Author:Wenqiang Wang*******************************/
/**********************Southern University of Science and Technology********************/
/***************************************************************************************/
#include <iostream>
#include "ElasticWaveEquation.h"

int main()
{
	GRID grid = { 300,200, 100 };
	Delta delta = { 0.0005f, 10.0f, 10.0f, 10.0f };
	//Delta delta = { 0.00175f, 100.0f, 100.0f, 100.0f };
	//Delta delta = { 0.0005f, 5.0f, 5.0f, 5.0f };
	int pml = 10;

	SOURCELOCATION sourceLocation = { 0 };
	STATIONLOCATION stationLocation;	

	
	sourceLocation.x = grid.xLength / 2;
	sourceLocation.y = grid.yLength / 2;
	sourceLocation.z = grid.zLength / 2;

	//stationLocation.x = grid.xLength - 2 * pml;
	//stationLocation.y = grid.yLength - 2 * pml;
	//stationLocation.z = grid.zLength - 2 * pml;


	stationLocation.x = grid.xLength / 2 + pml;
	stationLocation.y = grid.yLength / 2 + pml;
	stationLocation.z = grid.zLength / 2 + pml;

	int timeLength = 1500;
	//float mainFrequncy = 2.0f;
	float mainFrequncy = 20.0f;

	DifferenceCoefficient dc;

	dc.diff_coef[0] = 1.21124f;
	dc.diff_coef[1] = -8.97217e-2f;
	dc.diff_coef[2] = 1.38428e-2f;
	dc.diff_coef[3] = -1.76566e-3f;
	dc.diff_coef[4] = 1.18680e-4f;

	Medium medCPU;

	ALLOCATE_3D_Medium( medCPU );

	
	for ( int i = 0; i < grid.xLength * grid.yLength * grid.zLength; i++)
	{
		medCPU.buoyancy[i] = 1.0f / 2600.0f;
		medCPU.lambda[i] = 2600.0f * pow( 5000.0f, 2);
		medCPU.mu[i] = 2600.0f * pow(5000.0f, 2) / 3.0f;
	}


	dim3 blocksPerGrid( ( grid.xLength + 31 ) / 32, ( grid.yLength + 15 ) / 16, grid.zLength );//the num of blocks in every grid
	dim3 threadsPerBlock( 32, 16 ); //the num of threads in every block
	GPUDim gpuDim = { blocksPerGrid, threadsPerBlock };

	time_t start, stop;
	start = time(NULL);
	// clock_t start, stop;
	// start = clock();
	// cudaEvent_t start, stop;
	// cudaEventCreate( &start );
	// cudaEventCreate( &stop );
	// cudaEventRecord( start );

	ElasticWaveEquation ewq = ElasticWaveEquation(
		timeLength, gpuDim, pml,
		grid, delta, medCPU,
		sourceLocation, stationLocation,
		dc, 
		mainFrequncy );

	ewq.run();
	stop = time(NULL);
	// stop = clock();
	cout << "Time Loss:" << (stop - start) << endl;

	// float elapseTime;
	//cudaEventRecord( stop );
	// cudaEventElapsedTime( &elapseTime, start, stop );
	// cout << "Time Loss:" << elapseTime << endl;

	// cudaEventDestroy( start );
	// cudaEventDestroy( stop );
	
	DELETE_3D_Medium( medCPU );

	//system("pause");
}
