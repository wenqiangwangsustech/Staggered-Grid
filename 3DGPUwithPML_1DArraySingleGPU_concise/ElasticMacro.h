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
#include "common.h"
#include "helper_cuda.h"

#define PI 3.1415926f
#define R 0.001f
#define tau 4.3f 
#define INDEX( ix, iy, iz ) (ix) + (iy) * grid.xLength + (iz) * grid.yLength * grid.xLength
#define ORDER 10
#define order 5
//#define DAMPING0( V, pml, delta ) log10( 1.0f / R ) * tau * V / ( pml * delta )
#define DAMPING0( V, pml, delta ) tau * V / ( delta ) * ( 8.0f / 15.0f - 0.03f * ( pml ) + 1.0f / 1500.0f * ( ( pml ) * ( pml ) ) )
#define DAMPING( V, i, pml, delta ) DAMPING0( V, pml, delta ) * powf( ( i ) * 1.0f / ( ( pml ) * 1.0f ), 2 )
#define nAREA 26

#define START order 
#define CalculateAreaStart order + pml

#define CalculateAreaEndX grid.xLength - order - pml
#define CalculateAreaEndY grid.yLength - order - pml
#define CalculateAreaEndZ grid.zLength - order - pml

#define EndX grid.xLength - order
#define EndY grid.yLength - order
#define EndZ grid.zLength - order

#define AreaCondition( ix, iy, iz, start_x, start_y, start_z, end_x, end_y, end_z ) ( start_x <= ix && ix < end_x ) && ( start_y <= iy && iy < end_y ) && ( start_z <= iz && iz < end_z )

/*
#define VERTEX0Condition( ix, iy, iz, grid, order, pml ) ( ( START <= ix && ix < CalculateAreaStart ) && ( START <= iy && iy < CalculateAreaStart ) && ( START <= iz && iz < CalculateAreaStart ) )
#define VERTEX1Condition( ix, iy, iz, grid, order, pml ) ( ( START <= ix && ix < CalculateAreaStart ) && ( START <= iy && iy < CalculateAreaStart ) && ( CalculateAreaEndZ <= iz && iz < EndZ ) )
#define VERTEX2Condition( ix, iy, iz, grid, order, pml ) ( ( START <= ix && ix < CalculateAreaStart ) && ( CalculateAreaEndY <= iy && iy < EndY ) && ( START <= iz && iz < CalculateAreaStart ) )
#define VERTEX3Condition( ix, iy, iz, grid, order, pml ) ( ( START <= ix && ix < CalculateAreaStart ) && ( CalculateAreaEndY <= iy && iy < EndY ) && ( CalculateAreaEndZ <= iz && iz < EndZ ) )
#define VERTEX4Condition( ix, iy, iz, grid, order, pml ) ( ( CalculateAreaEndX <= ix && ix < EndX ) && ( START <= iy && iy < CalculateAreaStart ) && ( START <= iz && iz < CalculateAreaStart ) )
#define VERTEX5Condition( ix, iy, iz, grid, order, pml ) ( ( CalculateAreaEndX <= ix && ix < EndX ) && ( START <= iy && iy < CalculateAreaStart ) && ( CalculateAreaEndZ <= iz && iz < EndZ ) )
#define VERTEX6Condition( ix, iy, iz, grid, order, pml ) ( ( CalculateAreaEndX <= ix && ix < EndX ) && ( CalculateAreaEndY <= iy && iy < EndY ) && ( START <= iz && iz < CalculateAreaStart ) )
#define VERTEX7Condition( ix, iy, iz, grid, order, pml ) ( ( CalculateAreaEndX <= ix && ix < EndX ) && ( CalculateAreaEndY <= iy && iy < EndY ) && ( CalculateAreaEndZ <= iz && iz < EndZ ) )

#define SURFACE0Condition( ix, iy, iz, grid, order, pml ) ( ( START <= ix && ix < CalculateAreaStart ) && ( CalculateAreaStart <= iy && iy < CalculateAreaEndY ) && ( CalculateAreaStart <= iz && iz < CalculateAreaEndZ ) )
#define SURFACE1Condition( ix, iy, iz, grid, order, pml ) ( ( CalculateAreaEndX <= ix && ix < EndX ) && ( CalculateAreaStart <= iy && iy < CalculateAreaEndY ) && ( CalculateAreaStart <= iz && iz < CalculateAreaEndZ ) )
#define SURFACE2Condition( ix, iy, iz, grid, order, pml ) ( ( CalculateAreaStart <= ix && ix < CalculateAreaEndX ) && ( START <= iy && iy < CalculateAreaStart ) && ( CalculateAreaStart <= iz && iz < CalculateAreaEndZ ) )
#define SURFACE3Condition( ix, iy, iz, grid, order, pml ) ( ( CalculateAreaStart <= ix && ix < CalculateAreaEndX ) && ( CalculateAreaEndY <= iy && iy < EndY ) && ( CalculateAreaStart <= iz && iz < CalculateAreaEndZ ) )
#define SURFACE4Condition( ix, iy, iz, grid, order, pml ) ( ( CalculateAreaStart <= ix && ix < CalculateAreaEndX ) && ( CalculateAreaStart <= iy && iy < CalculateAreaEndY ) && ( START <= iz && iz < CalculateAreaStart ) )
#define SURFACE5Condition( ix, iy, iz, grid, order, pml ) ( ( CalculateAreaStart <= ix && ix < CalculateAreaEndX ) && ( CalculateAreaStart <= iy && iy < CalculateAreaEndY ) && ( CalculateAreaEndZ <= iz && iz < EndZ ) )

#define EDGE0Condition( ix, iy, iz, grid, order, pml )  ( ( START <= ix && ix < CalculateAreaStart ) && ( START <= iy && iy < CalculateAreaStart ) && ( CalculateAreaStart <= iz && iz < CalculateAreaEndZ ) )
#define EDGE1Condition( ix, iy, iz, grid, order, pml )  ( ( START <= ix && ix < CalculateAreaStart ) && ( CalculateAreaEndY <= iy && iy < EndY ) && ( CalculateAreaStart <= iz && iz < CalculateAreaEndZ ) )
#define EDGE2Condition( ix, iy, iz, grid, order, pml )  ( ( START <= ix && ix < CalculateAreaStart ) && ( CalculateAreaStart <= iy && iy < CalculateAreaEndY ) && ( START <= iz && iz < CalculateAreaStart ) )
#define EDGE3Condition( ix, iy, iz, grid, order, pml )  ( ( START <= ix && ix < CalculateAreaStart ) && ( CalculateAreaStart <= iy && iy < CalculateAreaEndY ) && ( CalculateAreaEndZ <= iz && iz < EndZ ) )

#define EDGE4Condition( ix, iy, iz, grid, order, pml )  ( ( CalculateAreaEndX <= ix && ix < EndX ) && ( START <= iy && iy < CalculateAreaStart ) && ( CalculateAreaStart <= iz && iz < CalculateAreaEndZ ) )
#define EDGE5Condition( ix, iy, iz, grid, order, pml )  ( ( CalculateAreaEndX <= ix && ix < EndX ) && ( CalculateAreaEndY <= iy && iy < EndY ) && ( CalculateAreaStart <= iz && iz < CalculateAreaEndZ ) )
#define EDGE6Condition( ix, iy, iz, grid, order, pml )  ( ( CalculateAreaEndX <= ix && ix < EndX ) && ( CalculateAreaStart <= iy && iy < CalculateAreaEndY ) && ( START <= iz && iz < CalculateAreaStart ) )
#define EDGE7Condition( ix, iy, iz, grid, order, pml )  ( ( CalculateAreaEndX <= ix && ix < EndX ) && ( CalculateAreaStart <= iy && iy < CalculateAreaEndY ) && ( CalculateAreaEndZ <= iz && iz < EndZ ) )

#define EDGE8Condition( ix, iy, iz, grid, order, pml )  ( ( CalculateAreaStart <= ix && ix < CalculateAreaEndX ) && ( START <= iy && iy < CalculateAreaStart ) && ( START <= iz && iz < CalculateAreaStart ) )
#define EDGE9Condition( ix, iy, iz, grid, order, pml )  ( ( CalculateAreaStart <= ix && ix < CalculateAreaEndX ) && ( START <= iy && iy < CalculateAreaStart ) && ( CalculateAreaEndZ <= iz && iz < EndZ ) )
#define EDGE10Condition( ix, iy, iz, grid, order, pml ) ( ( CalculateAreaStart <= ix && ix < CalculateAreaEndX ) && ( CalculateAreaEndY <= iy && iy < EndY ) && ( START <= iz && iz < CalculateAreaStart ) )
#define EDGE11Condition( ix, iy, iz, grid, order, pml ) ( ( CalculateAreaStart <= ix && ix < CalculateAreaEndX ) && ( CalculateAreaEndY <= iy && iy < EndY ) && ( CalculateAreaEndZ <= iz && iz < EndZ ) )
*/


#define VERTEX0 0
#define VERTEX1 1
#define VERTEX2 2
#define VERTEX3 3
#define VERTEX4 4
#define VERTEX5 5
#define VERTEX6 6
#define VERTEX7 7
#define SURFACE0 8
#define SURFACE1 9
#define SURFACE2 10
#define SURFACE3 11
#define SURFACE4 12
#define SURFACE5 13
#define EDGE0 14
#define EDGE1 15
#define EDGE2 16
#define EDGE3 17
#define EDGE4 18
#define EDGE5 19
#define EDGE6 20
#define EDGE7 21
#define EDGE8 22
#define EDGE9 23
#define EDGE10 24
#define EDGE11 25



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

#define CUDA_MALLOC_3D_VelocityStressPML( vsPML, memsize ) \
		cudaMalloc( ( void ** ) &vsPML.stressXXx , memsize );\
		cudaMalloc( ( void ** ) &vsPML.stressXXy , memsize );\
		cudaMalloc( ( void ** ) &vsPML.stressXXz , memsize );\
		cudaMalloc( ( void ** ) &vsPML.stressYYx , memsize );\
		cudaMalloc( ( void ** ) &vsPML.stressYYy , memsize );\
		cudaMalloc( ( void ** ) &vsPML.stressYYz , memsize );\
		cudaMalloc( ( void ** ) &vsPML.stressZZx , memsize );\
		cudaMalloc( ( void ** ) &vsPML.stressZZy , memsize );\
		cudaMalloc( ( void ** ) &vsPML.stressZZz , memsize );\
		cudaMalloc( ( void ** ) &vsPML.stressXYx , memsize );\
		cudaMalloc( ( void ** ) &vsPML.stressXYy , memsize );\
		/*cudaMalloc( ( void ** ) &vsPML.stressXYz , memsize );*/\
		cudaMalloc( ( void ** ) &vsPML.stressXZx , memsize );\
		/*cudaMalloc( ( void ** ) &vsPML.stressXZy , memsize );*/\
		cudaMalloc( ( void ** ) &vsPML.stressXZz , memsize );\
		/*cudaMalloc( ( void ** ) &vsPML.stressYZx , memsize );*/\
		cudaMalloc( ( void ** ) &vsPML.stressYZy , memsize );\
		cudaMalloc( ( void ** ) &vsPML.stressYZz , memsize );\
		cudaMalloc( ( void ** ) &vsPML.velocityXx, memsize );\
		cudaMalloc( ( void ** ) &vsPML.velocityXy, memsize );\
		cudaMalloc( ( void ** ) &vsPML.velocityXz, memsize );\
		cudaMalloc( ( void ** ) &vsPML.velocityYx, memsize );\
		cudaMalloc( ( void ** ) &vsPML.velocityYy, memsize );\
		cudaMalloc( ( void ** ) &vsPML.velocityYz, memsize );\
		cudaMalloc( ( void ** ) &vsPML.velocityZx, memsize );\
		cudaMalloc( ( void ** ) &vsPML.velocityZy, memsize );\
		cudaMalloc( ( void ** ) &vsPML.velocityZz, memsize );\
		cudaMemset( vsPML.stressXXx , 0, memsize );\
		cudaMemset( vsPML.stressXXy , 0, memsize );\
		cudaMemset( vsPML.stressXXz , 0, memsize );\
		cudaMemset( vsPML.stressYYx , 0, memsize );\
		cudaMemset( vsPML.stressYYy , 0, memsize );\
		cudaMemset( vsPML.stressYYz , 0, memsize );\
		cudaMemset( vsPML.stressZZx , 0, memsize );\
		cudaMemset( vsPML.stressZZy , 0, memsize );\
		cudaMemset( vsPML.stressZZz , 0, memsize );\
		cudaMemset( vsPML.stressXYx , 0, memsize );\
		cudaMemset( vsPML.stressXYy , 0, memsize );\
		/*cudaMemset( vsPML.stressXYz , 0, memsize );*/\
		cudaMemset( vsPML.stressXZx , 0, memsize );\
		/*cudaMemset( vsPML.stressXZy , 0, memsize );*/\
		cudaMemset( vsPML.stressXZz , 0, memsize );\
		/*cudaMemset( vsPML.stressYZx , 0, memsize );*/\
		cudaMemset( vsPML.stressYZy , 0, memsize );\
		cudaMemset( vsPML.stressYZz , 0, memsize );\
		cudaMemset( vsPML.velocityXx, 0, memsize );\
		cudaMemset( vsPML.velocityXy, 0, memsize );\
		cudaMemset( vsPML.velocityXz, 0, memsize );\
		cudaMemset( vsPML.velocityYx, 0, memsize );\
		cudaMemset( vsPML.velocityYy, 0, memsize );\
		cudaMemset( vsPML.velocityYz, 0, memsize );\
		cudaMemset( vsPML.velocityZx, 0, memsize );\
		cudaMemset( vsPML.velocityZy, 0, memsize );\
		cudaMemset( vsPML.velocityZz, 0, memsize );



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
		memset( vs.stressXX , 0, memsize );\
		memset( vs.stressYY , 0, memsize );\
		memset( vs.stressZZ , 0, memsize );\
		memset( vs.stressXY , 0, memsize );\
		memset( vs.stressXZ , 0, memsize );\
		memset( vs.stressYZ , 0, memsize );\
		memset( vs.velocityX, 0, memsize );\
		memset( vs.velocityY, 0, memsize );\
		memset( vs.velocityZ, 0, memsize );





#define CUDA_MALLOC_3D_Medium( med, memsize ) \
		cudaMalloc( ( void ** ) &med.buoyancy  , memsize );\
		cudaMalloc( ( void ** ) &med.lambda	  , memsize );\
		cudaMalloc( ( void ** ) &med.mu 	, memsize );


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

#define DELETE_MALLOC_3D_VelocityStressPML( vsPML ) \
		cudaFree( vsPML.stressXXx  );\
		cudaFree( vsPML.stressXXy  );\
		cudaFree( vsPML.stressXXz  );\
		cudaFree( vsPML.stressYYx  );\
		cudaFree( vsPML.stressYYy  );\
		cudaFree( vsPML.stressYYz  );\
		cudaFree( vsPML.stressZZx  );\
		cudaFree( vsPML.stressZZy  );\
		cudaFree( vsPML.stressZZz  );\
		cudaFree( vsPML.stressXYx  );\
		cudaFree( vsPML.stressXYy  );\
		/*cudaFree( vsPML.stressXYz  );*/\
		cudaFree( vsPML.stressXZx  );\
		/*cudaFree( vsPML.stressXZy  );*/\
		cudaFree( vsPML.stressXZz  );\
		/*cudaFree( vsPML.stressYZx  );*/\
		cudaFree( vsPML.stressYZy  );\
		cudaFree( vsPML.stressYZz  );\
		cudaFree( vsPML.velocityXx );\
		cudaFree( vsPML.velocityXy );\
		cudaFree( vsPML.velocityXz );\
		cudaFree( vsPML.velocityYx );\
		cudaFree( vsPML.velocityYy );\
		cudaFree( vsPML.velocityYz );\
		cudaFree( vsPML.velocityZx );\
		cudaFree( vsPML.velocityZy );\
		cudaFree( vsPML.velocityZz );


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

typedef struct VelocityStressPML
{
	float * stressXXx;
	float * stressXXy;
	float * stressXXz;
	float * stressYYx;
	float * stressYYy;
	float * stressYYz;
	float * stressZZx;
	float * stressZZy;
	float * stressZZz;
	float * stressXYx;
	float * stressXYy;
	//float * stressXYz;
	float * stressXZx;
	//float * stressXZy;
	float * stressXZz;
	//float * stressYZx;
	float * stressYZy;
	float * stressYZz;
	float * velocityXx;
	float * velocityXy;
	float * velocityXz;
	float * velocityYx;
	float * velocityYy;
	float * velocityYz;
	float * velocityZx;
	float * velocityZy;
	float * velocityZz;

	int startX;
	int startY;
	int startZ;

	int endX;
	int endY;
	int endZ;

	int xLength;
	int yLength;
	int zLength;

}VelocityStressPML;

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

typedef struct STATIONLOCATION
{
	int x;
	int y;
	int z;
}STATIONLOCATION;
/*
if ( ix == vsPML[VERTEX0].startX && iy == vsPML[VERTEX0].startY && iz == vsPML[VERTEX0].startZ )\
		{\
			printf("dampingMinusX = %f, dampingMinusY = %f, dampingMinusZ = %f, ", dampingMinusX, dampingMinusY, dampingMinusZ );\
			printf("dampingPlusX = %f, dampingPlusY = %f, dampingPlusZ = %f\n", dampingPlusX, dampingPlusY, dampingPlusZ );\
		}\
*/
#define VelocityPMLUpdate( AREA ) \
		float dampingMinusX = 1.0f - 0.5f * delta.t * dampingX;\
		float dampingMinusY = 1.0f - 0.5f * delta.t * dampingY;\
		float dampingMinusZ = 1.0f - 0.5f * delta.t * dampingZ;\
		float dampingPlusX  = 1.0f + 0.5f * delta.t * dampingX;\
		float dampingPlusY  = 1.0f + 0.5f * delta.t * dampingY;\
		float dampingPlusZ  = 1.0f + 0.5f * delta.t * dampingZ;\
		float diffX = 0.0f, diffY = 0.0f, diffZ = 0.0f;\
		int IX = ix - vsPML[AREA].startX;\
		int IY = iy - vsPML[AREA].startY;\
		int IZ = iz - vsPML[AREA].startZ;\
		int index = IX + IY * vsPML[AREA].xLength + IZ * vsPML[AREA].yLength * vsPML[AREA].xLength;\
		for ( i = 0;  i < order;  i++)\
		{\
			diffX += dc.diff_coef[i] * (vsGPU.stressXX[INDEX(ix + 1 + i, iy, iz)] - vsGPU.stressXX[INDEX(ix - i, iy, iz)]);\
			diffY += dc.diff_coef[i] * (vsGPU.stressXY[INDEX(ix, iy + i, iz)] - vsGPU.stressXY[INDEX(ix, iy - 1 - i, iz)]);\
			diffZ += dc.diff_coef[i] * (vsGPU.stressXZ[INDEX(ix, iy, iz + i)] - vsGPU.stressXZ[INDEX(ix, iy, iz - 1 - i)]);\
		}\
		diffX /= delta.x;\
		diffY /= delta.y;\
		diffZ /= delta.z;\
		vsPML[AREA].velocityXx[index] = ( dampingMinusX * vsPML[AREA].velocityXx[index] + delta.t * BUOYANCYBARX(medGPU.buoyancy, ix, iy, iz) * diffX ) / dampingPlusX;\
		vsPML[AREA].velocityXy[index] = ( dampingMinusY * vsPML[AREA].velocityXy[index] + delta.t * BUOYANCYBARX(medGPU.buoyancy, ix, iy, iz) * diffY ) / dampingPlusY;\
		vsPML[AREA].velocityXz[index] = ( dampingMinusZ * vsPML[AREA].velocityXz[index] + delta.t * BUOYANCYBARX(medGPU.buoyancy, ix, iy, iz) * diffZ ) / dampingPlusZ;\
		vsGPU.velocityX[INDEX(ix, iy, iz)] = vsPML[AREA].velocityXx[index] + vsPML[AREA].velocityXy[index] + vsPML[AREA].velocityXz[index];\
		\
		diffX = 0.0f, diffY = 0.0f, diffZ = 0.0f;\
		for (i = 0; i < order; i++)\
		{\
			diffX += dc.diff_coef[i] * (vsGPU.stressXY[INDEX(ix + i, iy, iz)] - vsGPU.stressXY[INDEX(ix - 1 - i, iy, iz)]);\
			diffY += dc.diff_coef[i] * (vsGPU.stressYY[INDEX(ix, iy + 1 + i, iz)] - vsGPU.stressYY[INDEX(ix, iy - i, iz)]);\
			diffZ += dc.diff_coef[i] * (vsGPU.stressYZ[INDEX(ix, iy, iz + i)] - vsGPU.stressYZ[INDEX(ix, iy, iz - 1 - i)]);\
		}\
		diffX /= delta.x;\
		diffY /= delta.y;\
		diffZ /= delta.z;\
		vsPML[AREA].velocityYx[index] = ( dampingMinusX * vsPML[AREA].velocityYx[index] + delta.t * BUOYANCYBARX(medGPU.buoyancy, ix, iy, iz) * diffX ) / dampingPlusX;\
		vsPML[AREA].velocityYy[index] = ( dampingMinusY * vsPML[AREA].velocityYy[index] + delta.t * BUOYANCYBARX(medGPU.buoyancy, ix, iy, iz) * diffY ) / dampingPlusY;\
		vsPML[AREA].velocityYz[index] = ( dampingMinusZ * vsPML[AREA].velocityYz[index] + delta.t * BUOYANCYBARX(medGPU.buoyancy, ix, iy, iz) * diffZ ) / dampingPlusZ;\
		vsGPU.velocityY[INDEX(ix, iy, iz)] = vsPML[AREA].velocityYx[index] + vsPML[AREA].velocityYy[index] + vsPML[AREA].velocityYz[index];\
		\
		diffX = 0.0f, diffY = 0.0f, diffZ = 0.0f;\
		for (i = 0; i < order; i++)\
		{\
			diffX += dc.diff_coef[i] * (vsGPU.stressXZ[INDEX(ix + i, iy, iz)] - vsGPU.stressXZ[INDEX(ix - 1 - i, iy, iz)]);\
			diffY += dc.diff_coef[i] * (vsGPU.stressYZ[INDEX(ix, iy + i, iz)] - vsGPU.stressYZ[INDEX(ix, iy - 1 - i, iz)]);\
			diffZ += dc.diff_coef[i] * (vsGPU.stressZZ[INDEX(ix, iy, iz + 1 + i)] - vsGPU.stressZZ[INDEX(ix, iy, iz - i)]);\
		}\
		diffX /= delta.x;\
		diffY /= delta.y;\
		diffZ /= delta.z;\
		vsPML[AREA].velocityZx[index] = ( dampingMinusX * vsPML[AREA].velocityZx[index] + delta.t * BUOYANCYBARX(medGPU.buoyancy, ix, iy, iz) * diffX ) / dampingPlusX;\
		vsPML[AREA].velocityZy[index] = ( dampingMinusY * vsPML[AREA].velocityZy[index] + delta.t * BUOYANCYBARX(medGPU.buoyancy, ix, iy, iz) * diffY ) / dampingPlusY;\
		vsPML[AREA].velocityZz[index] = ( dampingMinusZ * vsPML[AREA].velocityZz[index] + delta.t * BUOYANCYBARX(medGPU.buoyancy, ix, iy, iz) * diffZ ) / dampingPlusZ;\
		vsGPU.velocityZ[INDEX(ix, iy, iz)] = vsPML[AREA].velocityZx[index] + vsPML[AREA].velocityZy[index] + vsPML[AREA].velocityZz[index];


#define StressPMLUpdate( AREA ) \
		float dampingMinusX = 1.0f - 0.5f * delta.t * dampingX;\
		float dampingMinusY = 1.0f - 0.5f * delta.t * dampingY;\
		float dampingMinusZ = 1.0f - 0.5f * delta.t * dampingZ;\
		float dampingPlusX  = 1.0f + 0.5f * delta.t * dampingX;\
		float dampingPlusY  = 1.0f + 0.5f * delta.t * dampingY;\
		float dampingPlusZ  = 1.0f + 0.5f * delta.t * dampingZ;\
		float muBarXY = 0.0f, muBarXZ = 0.0f, muBarYZ = 0.0f;\
		\
		float diffX = 0.0f, diffY = 0.0f, diffZ = 0.0f;\
		int IX = ix - vsPML[AREA].startX;\
		int IY = iy - vsPML[AREA].startY;\
		int IZ = iz - vsPML[AREA].startZ;\
		int index = IX + IY * vsPML[AREA].xLength + IZ * vsPML[AREA].yLength * vsPML[AREA].xLength;\
		for (i = 0; i < order; i++)\
		{\
			diffX += dc.diff_coef[i] * (vsGPU.velocityX[INDEX(ix + i, iy, iz)] - vsGPU.velocityX[INDEX(ix - 1 - i, iy, iz)]);\
			diffY += dc.diff_coef[i] * (vsGPU.velocityY[INDEX(ix, iy + i, iz)] - vsGPU.velocityY[INDEX(ix, iy - 1 - i, iz)]);\
			diffZ += dc.diff_coef[i] * (vsGPU.velocityZ[INDEX(ix, iy, iz + i)] - vsGPU.velocityZ[INDEX(ix, iy, iz - 1 - i)]);\
		}\
		diffX /= delta.x;\
		diffY /= delta.y;\
		diffZ /= delta.z;\
		\
		vsPML[AREA].stressXXx[index] = ( dampingMinusX * vsPML[AREA].stressXXx[index] + delta.t * ( medGPU.lambda[INDEX(ix, iy, iz)] + 2.0f * medGPU.mu[INDEX(ix, iy, iz)] ) * diffX ) / dampingPlusX;\
		vsPML[AREA].stressXXy[index] = ( dampingMinusY * vsPML[AREA].stressXXy[index] + delta.t * ( medGPU.lambda[INDEX(ix, iy, iz)] * diffY ) ) / dampingPlusY;\
		vsPML[AREA].stressXXz[index] = ( dampingMinusZ * vsPML[AREA].stressXXz[index] + delta.t * ( medGPU.lambda[INDEX(ix, iy, iz)] * diffZ ) ) / dampingPlusZ;\
		vsGPU.stressXX[INDEX(ix, iy, iz)] = vsPML[AREA].stressXXx[index] + vsPML[AREA].stressXXy[index] + vsPML[AREA].stressXXz[index];\
		\
		vsPML[AREA].stressYYx[index] = ( dampingMinusX * vsPML[AREA].stressYYx[index] + delta.t * ( medGPU.lambda[INDEX(ix, iy, iz)] * diffX ) ) / dampingPlusX;\
		vsPML[AREA].stressYYy[index] = ( dampingMinusY * vsPML[AREA].stressYYy[index] + delta.t * ( medGPU.lambda[INDEX(ix, iy, iz)] + 2.0f * medGPU.mu[INDEX(ix, iy, iz)] ) * diffY ) / dampingPlusY;\
		vsPML[AREA].stressYYz[index] = ( dampingMinusZ * vsPML[AREA].stressYYz[index] + delta.t * ( medGPU.lambda[INDEX(ix, iy, iz)] * diffZ ) ) / dampingPlusZ;\
		vsGPU.stressYY[INDEX(ix, iy, iz)] = vsPML[AREA].stressYYx[index] + vsPML[AREA].stressYYy[index] + vsPML[AREA].stressYYz[index];\
		\
		vsPML[AREA].stressZZx[index] = ( dampingMinusX * vsPML[AREA].stressZZx[index] + delta.t * ( medGPU.lambda[INDEX(ix, iy, iz)] * diffX ) ) / dampingPlusX;\
		vsPML[AREA].stressZZy[index] = ( dampingMinusY * vsPML[AREA].stressZZy[index] + delta.t * ( medGPU.lambda[INDEX(ix, iy, iz)] * diffY ) ) / dampingPlusY;\
		vsPML[AREA].stressZZz[index] = ( dampingMinusZ * vsPML[AREA].stressZZz[index] + delta.t * ( medGPU.lambda[INDEX(ix, iy, iz)] + 2.0f * medGPU.mu[INDEX(ix, iy, iz)] ) * diffZ ) / dampingPlusZ;\
		vsGPU.stressZZ[INDEX(ix, iy, iz)] = vsPML[AREA].stressZZx[index] + vsPML[AREA].stressZZy[index] + vsPML[AREA].stressZZz[index];\
		\
		muBarXY = MUBARXY(medGPU.mu, ix, iy, iz);\
		muBarXZ = MUBARXZ(medGPU.mu, ix, iy, iz);\
		muBarYZ = MUBARYZ(medGPU.mu, ix, iy, iz);\
		\
		diffX = 0.0f, diffY = 0.0f;\
		for (i = 0; i < order; i++)\
		{\
			diffY += dc.diff_coef[i] * (vsGPU.velocityX[INDEX(ix, iy + 1 + i, iz)] - vsGPU.velocityX[INDEX(ix, iy - i, iz)]);\
			diffX += dc.diff_coef[i] * (vsGPU.velocityY[INDEX(ix + 1 + i, iy, iz)] - vsGPU.velocityY[INDEX(ix - i, iy, iz)]);\
		}\
		diffX /= delta.x;\
		diffY /= delta.y;\
		vsPML[AREA].stressXYx[index] = ( dampingMinusX * vsPML[AREA].stressXYx[index] + delta.t * muBarXY * diffX ) / dampingPlusX;\
		vsPML[AREA].stressXYy[index] = ( dampingMinusY * vsPML[AREA].stressXYy[index] + delta.t * muBarXY * diffY ) / dampingPlusY;\
		vsGPU.stressXY[INDEX(ix, iy, iz)] = vsPML[AREA].stressXYx[index] + vsPML[AREA].stressXYy[index];\
		\
		diffX = 0.0f, diffZ = 0.0f;\
		for (i = 0; i < order; i++)\
		{\
			diffZ += dc.diff_coef[i] * (vsGPU.velocityX[INDEX(ix, iy, iz + 1 + i)] - vsGPU.velocityX[INDEX(ix, iy, iz - i)]);\
			diffX += dc.diff_coef[i] * (vsGPU.velocityZ[INDEX(ix + 1 + i, iy, iz)] - vsGPU.velocityZ[INDEX(ix - i, iy, iz)]);\
		}\
		diffX /= delta.x;\
		diffZ /= delta.z;\
		vsPML[AREA].stressXZx[index] = ( dampingMinusX * vsPML[AREA].stressXZx[index] + delta.t * muBarXZ * diffX ) / dampingPlusX;\
		vsPML[AREA].stressXZz[index] = ( dampingMinusZ * vsPML[AREA].stressXZz[index] + delta.t * muBarXZ * diffZ ) / dampingPlusZ;\
		vsGPU.stressXZ[INDEX(ix, iy, iz)] = vsPML[AREA].stressXZx[index] + vsPML[AREA].stressXZz[index];\
		\
		diffY = 0.0f, diffZ = 0.0f;\
		for (i = 0; i < order; i++)\
		{\
			diffZ += dc.diff_coef[i] * (vsGPU.velocityY[INDEX(ix, iy, iz + 1 + i)] - vsGPU.velocityY[INDEX(ix, iy, iz - i)]);\
			diffY += dc.diff_coef[i] * (vsGPU.velocityZ[INDEX(ix, iy + 1 + i, iz)] - vsGPU.velocityZ[INDEX(ix, iy - i, iz)]);\
		}\
		diffY /= delta.y;\
		diffZ /= delta.z;\
		vsPML[AREA].stressYZy[index] = ( dampingMinusY * vsPML[AREA].stressYZy[index] + delta.t * muBarYZ * diffY ) / dampingPlusY;\
		vsPML[AREA].stressYZz[index] = ( dampingMinusZ * vsPML[AREA].stressYZz[index] + delta.t * muBarYZ * diffZ ) / dampingPlusZ;\
		vsGPU.stressYZ[INDEX(ix, iy, iz)] = vsPML[AREA].stressYZy[index] + vsPML[AREA].stressYZz[index];

#endif // !ELASTIC_MACRO_H

