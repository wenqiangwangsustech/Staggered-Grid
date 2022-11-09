/***************************************************************************************/
/******************************************2019-01-05***********************************/
/************************************Author:Wenqiang Wang*******************************/
/**********************Southern University of Science and Technology********************/
/***************************************************************************************/
#ifndef ELASTIC_MACRO_H
#define ELASTIC_MACRO_H
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
using namespace std;
#include <memory.h>
#include <cstdio>
#include <cmath>
#include <ctime>

#define PI 3.1415926f
#define INDEX( ix, iy, iz ) (ix) + (iy) * grid.xLength + (iz) * grid.yLength * grid.xLength
#define ORDER 10

#define BUOYANCYBAR( buoyancy1, buoyancy2 ) ( ( buoyancy1 + buoyancy2 ) / 2 ) 
#define BUOYANCYBARX( buoyancy, ix, iy, iz ) BUOYANCYBAR( buoyancy[INDEX( ix, iy, iz )], buoyancy[INDEX( ix + 1, iy, iz )] ) 
#define BUOYANCYBARY( buoyancy, ix, iy, iz ) BUOYANCYBAR( buoyancy[INDEX( ix, iy, iz )], buoyancy[INDEX( ix, iy + 1, iz )] )
#define BUOYANCYBARZ( buoyancy, ix, iy, iz ) BUOYANCYBAR( buoyancy[INDEX( ix, iy, iz )], buoyancy[INDEX( ix, iy, iz + 1 )] )

#define MUBAR( mu1, mu2, mu3, mu4 ) \
		( 4.0f / ( 1 / mu1 + 1 / mu2 + 1 / mu3 + 1 / mu4 ) )

#define MUBARXY( mu, ix, iy, iz ) \
		MUBAR( mu[INDEX( ix, iy, iz )], mu[INDEX( ix + 1, iy, iz )], mu[INDEX( ix, iy + 1, iz )], mu[INDEX( ix + 1, iy + 1, iz )] )
#define MUBARXZ( mu, ix, iy, iz ) \
		MUBAR( mu[INDEX( ix, iy, iz )], mu[INDEX( ix + 1, iy, iz )], mu[INDEX( ix, iy, iz + 1 )], mu[INDEX( ix + 1, iy + 1, iz )] )
#define MUBARYZ( mu, ix, iy, iz ) \
		MUBAR( mu[INDEX( ix, iy, iz )], mu[INDEX( ix, iy + 1, iz )], mu[INDEX( ix, iy, iz + 1 )], mu[INDEX( ix, iy + 1, iz + 1 )] )


#define CUDA_MALLOC_3D_VelocityStress( vs, memsize ) \
		cudaMalloc( ( void ** ) &vs.stressXX , memsize );\
		cudaMalloc( ( void ** ) &vs.stressYY , memsize );\
		cudaMalloc( ( void ** ) &vs.stressZZ , memsize );\
		cudaMalloc( ( void ** ) &vs.stressXY , memsize );\
		cudaMalloc( ( void ** ) &vs.stressXZ , memsize );\
		cudaMalloc( ( void ** ) &vs.stressYZ , memsize );\
		cudaMalloc( ( void ** ) &vs.velocityX, memsize );\
		cudaMalloc( ( void ** ) &vs.velocityY, memsize );\
		cudaMalloc( ( void ** ) &vs.velocityZ, memsize );\
		cudaMemset( vs.stressXX , 0, memsize );\
		cudaMemset( vs.stressYY , 0, memsize );\
		cudaMemset( vs.stressZZ , 0, memsize );\
		cudaMemset( vs.stressXY , 0, memsize );\
		cudaMemset( vs.stressXZ , 0, memsize );\
		cudaMemset( vs.stressYZ , 0, memsize );\
		cudaMemset( vs.velocityX, 0, memsize );\
		cudaMemset( vs.velocityY, 0, memsize );\
		cudaMemset( vs.velocityZ, 0, memsize );

#define ALLOCATE_3D_VelocityStress( vs, memsize ) \
		vs.stressXX  = new float[ memsize ];\
		vs.stressYY  = new float[ memsize ];\
		vs.stressZZ  = new float[ memsize ];\
		vs.stressXY  = new float[ memsize ];\
		vs.stressXZ  = new float[ memsize ];\
		vs.stressYZ  = new float[ memsize ];\
		vs.velocityX = new float[ memsize ];\
		vs.velocityY = new float[ memsize ];\
		vs.velocityZ = new float[ memsize ];\
		memset( vs.stressXX , 0, memsize * sizeof( float ) );\
		memset( vs.stressYY , 0, memsize * sizeof( float ) );\
		memset( vs.stressZZ , 0, memsize * sizeof( float ) );\
		memset( vs.stressXY , 0, memsize * sizeof( float ) );\
		memset( vs.stressXZ , 0, memsize * sizeof( float ) );\
		memset( vs.stressYZ , 0, memsize * sizeof( float ) );\
		memset( vs.velocityX, 0, memsize * sizeof( float ) );\
		memset( vs.velocityY, 0, memsize * sizeof( float ) );\
		memset( vs.velocityZ, 0, memsize * sizeof( float ) );


#define CUDA_MALLOC_3D_Medium( med, memsize ) \
		cudaMalloc( ( void ** ) &medGPU.buoyancy  , memsize );\
		cudaMalloc( ( void ** ) &medGPU.lambda	  , memsize );\
		cudaMalloc( ( void ** ) &medGPU.mu 	 	  , memsize );


#define CUDA_FREE_3D_VelocityStress( vs ) \
		cudaFree( vs.stressXX  ); \
		cudaFree( vs.stressYY  ); \
		cudaFree( vs.stressZZ  ); \
		cudaFree( vs.stressXY  ); \
		cudaFree( vs.stressXZ  ); \
		cudaFree( vs.stressYZ  ); \
		cudaFree( vs.velocityX ); \
		cudaFree( vs.velocityY ); \
		cudaFree( vs.velocityZ );

#define DELETE_3D_VelocityStress( vs )\
		delete[] vs.stressXX; \
		delete[] vs.stressYY; \
		delete[] vs.stressZZ; \
		delete[] vs.stressXY; \
		delete[] vs.stressXZ; \
		delete[] vs.stressYZ; \
		delete[] vs.velocityX; \
		delete[] vs.velocityY; \
		delete[] vs.velocityZ;

#define CUDA_FREE_3D_Medium( med ) \
		cudaFree( med.buoyancy );\
		cudaFree( med.lambda );\
		cudaFree( med.mu );

#define ALLOCATE_3D_Medium( med ) \
		med.buoyancy 	= new float[ grid.xLength * grid.yLength * grid.zLength ];\
		med.lambda 		= new float[ grid.xLength * grid.yLength * grid.zLength ];\
		med.mu 			= new float[ grid.xLength * grid.yLength * grid.zLength ];

#define DELETE_3D_Medium( med ) \
		delete [] med.buoyancy;\
		delete [] med.lambda;\
		delete [] med.mu;

typedef struct SOURCELOCATION
{
	int x;
	int y;
	int z;
}SOURCELOCATION;

typedef struct GRID
{
	int xLength;
	int yLength;
	int zLength;
}GRID;

typedef struct Delta
{
	float t;
	float x;
	float y;
	float z;
}Delta;

typedef struct VelocityStress
{
	float * stressXX;
	float * stressYY;
	float * stressZZ;
	float * stressXY;
	float * stressXZ;
	float * stressYZ;
	float * velocityX;
	float * velocityY;
	float * velocityZ;
}VelocityStress;

typedef struct Medium
{
	float * buoyancy;
	float * lambda;
	float * mu;
}Medium;

typedef struct DifferenceCoefficient
{
	float diff_coef[ORDER / 2];
}DifferenceCoefficient;

typedef struct GPUDim {
	dim3 blocksPerGrid;
	dim3 threadsPerBlock;
}GPUDim;

#endif // !ELASTIC_MACRO_H

