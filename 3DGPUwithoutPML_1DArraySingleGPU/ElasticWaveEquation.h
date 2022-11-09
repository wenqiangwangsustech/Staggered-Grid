/***************************************************************************************/
/******************************************2019-01-05***********************************/
/************************************Author:Wenqiang Wang*******************************/
/**********************Southern University of Science and Technology********************/
/***************************************************************************************/
#ifndef ELASTIC_WAVE_EQUATION_H
#define ELASTIC_WAVE_EQUATION_H
#include "ElasticMacro.h"
extern __global__ void locateSource( GRID grid, float deltaT, float sourceValue, VelocityStress vsGPU,  SOURCELOCATION sourceLocation );
extern __global__ void velocityUpdate( GRID grid, Delta delta, SOURCELOCATION sourceLocation, VelocityStress vsGPU, DifferenceCoefficient dc, float * buoyancy );
extern __global__ void stressUpdate( GRID grid, Delta delta, VelocityStress vsGPU, DifferenceCoefficient dc, Medium medGPU );
class ElasticWaveEquation
{
public:
	VelocityStress vsGPU;
	VelocityStress vsCPU;
	
	int timeLength;
	GRID grid; Delta delta; Medium medGPU; Medium medCPU;
	SOURCELOCATION sourceLocation;
	DifferenceCoefficient dc;
	float mainFrequncy;
	
	GPUDim gpuDim;

	ElasticWaveEquation(
		int timeLength, GPUDim gpuDim,
		GRID grid, Delta delta,Medium med,
		SOURCELOCATION sourceLocation, DifferenceCoefficient dc, 
		float mainFrequncy );
	~ElasticWaveEquation();
	
	void propagate( float sourceValue );
	float sourceFunction( int timeIndex );
	void snapShootVolum(int it);
	void snapShootSlice(int it);
	void run();

};
/* dim3 blocksPerGrid( ( grid.yLength + 32 ) / 32, ( grid.yLength + 15 ) / 16, grid.zLength );//the num of blocks in every grid
	dim3 threadsPerBlock( 32, 16 ); //the num of threads in every block
	*/
#endif //ELASTIC_WAVE_EQUATION_H

