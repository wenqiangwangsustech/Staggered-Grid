/***************************************************************************************/
/******************************************2019-01-04***********************************/
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

#define DEF_3D_POINTER( type, pointer3D )  type * pointer3D;
#define DEF_3D_VelocityStress( type ) \
		type * stressXX; \
		type * stressYY; \
		type * stressZZ; \
		type * stressXY; \
		type * stressXZ; \
		type * stressYZ; \
		type * velocityX; \
		type * velocityY; \
		type * velocityZ;

#define DEF_3D_pragma( type ) \
		type * buoyancy; \
		type * lambda; \
		type * mu;


#define ALLOCATE_3D_VelocityStress( type, xLength, yLength, zLength ) \
		stressXX  = new type [ xLength * yLength * zLength ];\
		stressYY  = new type [ xLength * yLength * zLength ];\
		stressZZ  = new type [ xLength * yLength * zLength ];\
		stressXY  = new type [ xLength * yLength * zLength ];\
		stressXZ  = new type [ xLength * yLength * zLength ];\
		stressYZ  = new type [ xLength * yLength * zLength ];\
		velocityX = new type [ xLength * yLength * zLength ];\
		velocityY = new type [ xLength * yLength * zLength ];\
		velocityZ = new type [ xLength * yLength * zLength ];\
		memset(stressXX , 0, sizeof(float) * (xLength * yLength * zLength));\
		memset(stressYY , 0, sizeof(float) * (xLength * yLength * zLength));\
		memset(stressZZ , 0, sizeof(float) * (xLength * yLength * zLength));\
		memset(stressXY , 0, sizeof(float) * (xLength * yLength * zLength));\
		memset(stressXZ , 0, sizeof(float) * (xLength * yLength * zLength));\
		memset(stressYZ , 0, sizeof(float) * (xLength * yLength * zLength));\
		memset(velocityX, 0, sizeof(float) * (xLength * yLength * zLength));\
		memset(velocityY, 0, sizeof(float) * (xLength * yLength * zLength));\
		memset(velocityZ, 0, sizeof(float) * (xLength * yLength * zLength));

#define ALLOCATE_3D_pragma( type, xLength, yLength, zLength ) \
		lambda  = new type [ xLength * yLength * zLength ];\
		mu  = new type [ xLength * yLength * zLength ];\
		buoyancy  = new type [ xLength * yLength * zLength ];

#define DELETE_3D_VelocityStress( type )\
		delete[] stressXX; \
		delete[] stressYY; \
		delete[] stressZZ; \
		delete[] stressXY; \
		delete[] stressXZ; \
		delete[] stressYZ; \
		delete[] velocityX; \
		delete[] velocityY; \
		delete[] velocityZ;

#define DELETE_3D_pragma( type ) \
		delete[] lambda; \
		delete[] mu; \
		delete[] buoyancy;

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

#endif // !ELASTIC_MACRO_H

