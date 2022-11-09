/***************************************************************************************/
/******************************************2019-01-05***********************************/
/************************************Author:Wenqiang Wang*******************************/
/**********************Southern University of Science and Technology********************/
/***************************************************************************************/
#ifndef ELASTIC_WAVE_EQUATION_H
#define ELASTIC_WAVE_EQUATION_H
#include "ElasticMacro.h"
extern __global__ void locateSource( GRID grid, float deltaT, float sourceValue, VelocityStress vsGPU );
extern __global__ void velocityUpdate( GRID grid, Delta delta, VelocityStress vsGPU, DifferenceCoefficient dc, Medium medGPU, int pml );
extern __global__ void stressUpdate( GRID grid, Delta delta, VelocityStress vsGPU, DifferenceCoefficient dc, Medium medGPU, int pml );
extern __global__ void recordStationData( GRID grid, int timeIndex, float * stationDataGPU, VelocityStress vsGPU );
__device__ __constant__ VelocityStressPML vsPML[nAREA];
__device__ __constant__ SOURCELOCATION sourceLocationGPU;
__device__ __constant__ STATIONLOCATION stationLocationGPU;

class ElasticWaveEquation
{
public:
	VelocityStress vsGPU;
	VelocityStress vsCPU;
	
	int timeLength;
	GRID grid; Delta delta; Medium medGPU; Medium medCPU;
	DifferenceCoefficient dc;
	float mainFrequncy;
	
	GPUDim gpuDim;
	int pml;

	VelocityStressPML vsPMLCPU[26];
	


	SOURCELOCATION sourceLocationCPU;
	STATIONLOCATION stationLocationCPU;

	float * stationDataCPU;
	float * stationDataGPU;

	ElasticWaveEquation(
		int timeLength, GPUDim gpuDim, int pml,
		GRID grid, Delta delta,Medium medCPU,
		SOURCELOCATION sourceLocation, STATIONLOCATION stationLocationCPU,
		DifferenceCoefficient dc, 
		float mainFrequncy );
	~ElasticWaveEquation();
	
	void propagate( int timeIndex, float sourceValue );
	float sourceFunction( int timeIndex );
	void snapShootVolum(int it);
	void snapShootSlice(int it);
	void run();
	void calStartEndPoint( );
	void vsPMLCPUTovsPML( );
	void stationLocationCPUTostationLocationGPU();
	void sourceLocationCPUTosourceLocationGPU();
	void stationDataOutput();

};
/* dim3 blocksPerGrid( ( grid.yLength + 32 ) / 32, ( grid.yLength + 15 ) / 16, grid.zLength );//the num of blocks in every grid
	dim3 threadsPerBlock( 32, 16 ); //the num of threads in every block
	*/
#endif //ELASTIC_WAVE_EQUATION_H

