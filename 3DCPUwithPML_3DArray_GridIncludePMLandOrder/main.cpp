/***************************************************************************************/
/******************************************2018-12-21***********************************/
/************************************Author:Wenqiang Wang*******************************/
/**********************Southern University of Science and Technology********************/
/***************************************************************************************/

#include <iostream>
using namespace std;
#include "ElasticWaveEquation3D.h"
int main()
{
	GRID grid = { 101, 101, 101 };
	int pml = 10;
	SOURCELOCATION sourceLocation = { 0 };

	sourceLocation.x = grid.xLength / 2 + pml + ORDER / 2;
	sourceLocation.y = grid.yLength / 2 + pml + ORDER / 2;
	sourceLocation.z = grid.zLength / 2 + pml + ORDER / 2;

	int timeLength = 600;

	DEF_3D_POINTER( float, density );
	DEF_3D_POINTER( float, buoyancy );
	DEF_3D_POINTER( float, velocity );
	DEF_3D_POINTER( float, lambda );
	DEF_3D_POINTER( float, mu );

	ALLOCATE_3D(float, density,	 grid.xLength + 2 * pml + ORDER, grid.yLength + 2 * pml + ORDER, grid.zLength + 2 * pml  + ORDER);
	ALLOCATE_3D(float, buoyancy, grid.xLength + 2 * pml + ORDER, grid.yLength + 2 * pml + ORDER, grid.zLength + 2 * pml  + ORDER);
	ALLOCATE_3D(float, velocity, grid.xLength + 2 * pml + ORDER, grid.yLength + 2 * pml + ORDER, grid.zLength + 2 * pml  + ORDER);
	ALLOCATE_3D(float, lambda,	 grid.xLength + 2 * pml + ORDER, grid.yLength + 2 * pml + ORDER, grid.zLength + 2 * pml  + ORDER);
	ALLOCATE_3D(float, mu,		 grid.xLength + 2 * pml + ORDER, grid.yLength + 2 * pml + ORDER, grid.zLength + 2 * pml  + ORDER);

	FOR_LOOP_3D( i, 0, j, 0, k, 0, grid, pml, ORDER,density[i][j][k] = 2500.0f );
	FOR_LOOP_3D( i, 0, j, 0, k, 0, grid, pml, ORDER,velocity[i][j][k] = Vs );
				 					
	FOR_LOOP_3D( i, 0, j, 0, k, 0, grid, pml, ORDER,buoyancy[i][j][k] = 1.0f / density[i][j][k] );
	FOR_LOOP_3D( i, 0, j, 0, k, 0, grid, pml, ORDER,lambda[i][j][k] = density[i][j][k] * pow( velocity[i][j][k], 2 ) );
	FOR_LOOP_3D( i, 0, j, 0, k, 0, grid, pml, ORDER,mu[i][j][k] = density[i][j][k] * pow( velocity[i][j][k], 2 ) / 3.0f );
	
	DELETE_3D(float, density,	 grid.xLength + 2 * pml + ORDER, grid.yLength + 2 * pml + ORDER, grid.zLength + 2 * pml + ORDER);
	DELETE_3D(float, velocity,   grid.xLength + 2 * pml + ORDER, grid.yLength + 2 * pml + ORDER, grid.zLength + 2 * pml + ORDER);

	ElasticWaveEquation3D ewq = 
		ElasticWaveEquation3D( 
			grid, timeLength, pml, 
			buoyancy, lambda, mu, 
			sourceLocation );
	
	time_t start, stop;
	start = time( NULL );
	ewq.propagate( );
	//ewq.calculateDamping();
	stop = time(NULL);

	cout << ( stop - start ) << endl;
    
	DELETE_3D(float, buoyancy,   grid.xLength + 2 * pml + ORDER, grid.yLength + 2 * pml + ORDER, grid.zLength + 2 * pml + ORDER);
	DELETE_3D(float, lambda,	 grid.xLength + 2 * pml + ORDER, grid.yLength + 2 * pml + ORDER, grid.zLength + 2 * pml + ORDER);
	DELETE_3D(float, mu,		 grid.xLength + 2 * pml + ORDER, grid.yLength + 2 * pml + ORDER, grid.zLength + 2 * pml + ORDER);
	system("pause");								   
}