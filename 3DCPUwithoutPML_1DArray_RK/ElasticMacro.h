/***************************************************************************************/
/******************************************2019-01-24***********************************/
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

#define a_1 -0.30874f
#define a0  -0.6326f
#define a1  1.2330f
#define a2  -0.3334f
#define a3  0.04168f

#define alpha4_1 0.0f
#define alpha4_2 0.5f
#define alpha4_3 0.5f
#define alpha4_4 1.0f

#define beta4_1 1.0f / 6.0f
#define beta4_2 1.0f / 3.0f
#define beta4_3 1.0f / 3.0f
#define beta4_4 1.0f / 6.0f

#define alpha6_1 0.0f
#define alpha6_2 0.353323f
#define alpha6_3 0.999597f
#define alpha6_4 0.152188f
#define alpha6_5 0.534216f
#define alpha6_6 0.603907f

#define beta6_1 0.0467621f
#define beta6_2 0.137286f
#define beta6_3 0.170975f
#define beta6_4 0.197572f
#define beta6_5 0.282263f
#define beta6_6 0.165142f

//Wave denotes velocity and stress

#define LxF( Wave, ix, iy, iz, delta )  ( 1.0f / delta ) * \
( a_1 * Wave[INDEX( ix - 1, iy, iz )] + a0 * Wave[INDEX( ix, iy, iz )] + a1 * Wave[INDEX( ix + 1, iy, iz )] + a2 * Wave[INDEX( ix + 2, iy, iz )] + a3 * Wave[INDEX( ix + 3, iy, iz )] )

#define LyF( Wave, ix, iy, iz, delta )  ( 1.0f / delta ) * \
( a_1 * Wave[INDEX( ix, iy - 1, iz )] + a0 * Wave[INDEX( ix, iy, iz )] + a1 * Wave[INDEX( ix, iy + 1, iz )] + a2 * Wave[INDEX( ix, iy + 2, iz )] + a3 * Wave[INDEX( ix, iy + 3, iz )] )

#define LzF( Wave, ix, iy, iz, delta )  ( 1.0f / delta ) * \
( a_1 * Wave[INDEX( ix, iy, iz - 1 )] + a0 * Wave[INDEX( ix, iy, iz )] + a1 * Wave[INDEX( ix, iy, iz + 1 )] + a2 * Wave[INDEX( ix, iy, iz + 2 )] + a3 * Wave[INDEX( ix, iy, iz + 3 )] )


#define LxB( Wave, ix, iy, iz, delta ) ( - 1.0f / delta ) * \
( a_1 * Wave[INDEX( ix + 1, iy, iz )] + a0 * Wave[INDEX( ix, iy, iz )] + a1 * Wave[INDEX( ix - 1, iy, iz )] + a2 * Wave[INDEX( ix - 2, iy, iz )] + a3 * Wave[INDEX( ix - 3, iy, iz )] )

#define LyB( Wave, ix, iy, iz, delta ) ( - 1.0f / delta ) * \
( a_1 * Wave[INDEX( ix, iy + 1, iz )] + a0 * Wave[INDEX( ix, iy, iz )] + a1 * Wave[INDEX( ix, iy - 1, iz )] + a2 * Wave[INDEX( ix, iy - 2, iz )] + a3 * Wave[INDEX( ix, iy - 3, iz )] )

#define LzB( Wave, ix, iy, iz, delta ) ( - 1.0f / delta ) * \
( a_1 * Wave[INDEX( ix, iy, iz + 1 )] + a0 * Wave[INDEX( ix, iy, iz )] + a1 * Wave[INDEX( ix, iy, iz - 1 )] + a2 * Wave[INDEX( ix, iy, iz - 2 )] + a3 * Wave[INDEX( ix, iy, iz - 3 )] )


// #define LxF( Wave, ix, iy, iz, delta, alpha, h ) 1.0f / delta * \
// ( a_1 * ( Wave[INDEX( ix - 1, iy, iz )] + alpha * h ) + a0 * ( Wave[INDEX( ix, iy, iz )] + alpha * h ) + a1 * ( Wave[INDEX( ix + 1, iy, iz )] + alpha * h ) + a2 * ( Wave[INDEX( ix + 2, iy, iz )] + alpha * h ) + a3 * ( Wave[INDEX( ix + 3, iy, iz )] + alpha * h ) )

// #define LyF( Wave, ix, iy, iz, delta, alpha, h ) 1.0f / delta * \
// ( a_1 * ( Wave[INDEX( ix, iy - 1, iz )] + alpha * h ) + a0 * ( Wave[INDEX( ix, iy, iz )] + alpha * h ) + a1 * ( Wave[INDEX( ix, iy + 1, iz )] + alpha * h ) + a2 * ( Wave[INDEX( ix, iy + 2, iz )] + alpha * h ) + a3 * ( Wave[INDEX( ix, iy + 3, iz )] + alpha * h ) )

// #define LzF( Wave, ix, iy, iz, delta, alpha, h ) 1.0f / delta * \
// ( a_1 * ( Wave[INDEX( ix, iy, iz - 1 )] + alpha * h ) + a0 * ( Wave[INDEX( ix, iy, iz )] + alpha * h ) + a1 * ( Wave[INDEX( ix, iy, iz + 1 )] + alpha * h ) + a2 * ( Wave[INDEX( ix, iy, iz + 2 )] + alpha * h ) + a3 * ( Wave[INDEX( ix, iy, iz + 3 )] + alpha * h ) )


// #define LxB( Wave, ix, iy, iz, delta, alpha, h ) - 1.0f / delta * \
// ( a_1 * ( Wave[INDEX( ix + 1, iy, iz )] + alpha * h ) + a0 * ( Wave[INDEX( ix, iy, iz )] + alpha * h ) + a1 * ( Wave[INDEX( ix - 1, iy, iz )] + alpha * h ) + a2 * ( Wave[INDEX( ix - 2, iy, iz )] + alpha * h ) + a3 * ( Wave[INDEX( ix - 3, iy, iz )] + alpha * h ) )

// #define LyB( Wave, ix, iy, iz, delta, alpha, h ) - 1.0f / delta * \
// ( a_1 * ( Wave[INDEX( ix, iy + 1, iz )] + alpha * h ) + a0 * ( Wave[INDEX( ix, iy, iz )] + alpha * h ) + a1 * ( Wave[INDEX( ix, iy - 1, iz )] + alpha * h ) + a2 * ( Wave[INDEX( ix, iy - 2, iz )] + alpha * h ) + a3 * ( Wave[INDEX( ix, iy - 3, iz )] + alpha * h ) )

// #define LzB( Wave, ix, iy, iz, delta, alpha, h ) - 1.0f / delta * \
// ( a_1 * ( Wave[INDEX( ix, iy, iz + 1 )] + alpha * h ) + a0 * ( Wave[INDEX( ix, iy, iz )] + alpha * h ) + a1 * ( Wave[INDEX( ix, iy, iz - 1 )] + alpha * h ) + a2 * ( Wave[INDEX( ix, iy, iz - 2 )] + alpha * h ) + a3 * ( Wave[INDEX( ix, iy, iz - 3 )] + alpha * h ) )
#define LxF( Wave, ix, iy, iz, delta )  ( 1.0f / delta ) * \
( a_1 * Wave[INDEX( ix - 1, iy, iz )] + a0 * Wave[INDEX( ix, iy, iz )] + a1 * Wave[INDEX( ix + 1, iy, iz )] + a2 * Wave[INDEX( ix + 2, iy, iz )] + a3 * Wave[INDEX( ix + 3, iy, iz )] ) + A * alpha * h

#define LxF( Wave, ix, iy, iz, delta, alpha, h ) 1.0f / delta * \
( a_1 * ( Wave[INDEX( ix - 1, iy, iz )] + alpha * h ) + a0 * ( Wave[INDEX( ix, iy, iz )] + alpha * h ) + a1 * ( Wave[INDEX( ix + 1, iy, iz )] + alpha * h ) + a2 * ( Wave[INDEX( ix + 2, iy, iz )] + alpha * h ) + a3 * ( Wave[INDEX( ix + 3, iy, iz )] + alpha * h ) )

#define LyF( Wave, ix, iy, iz, delta, alpha, h ) 1.0f / delta * \
( a_1 * ( Wave[INDEX( ix, iy - 1, iz )] + alpha * h ) + a0 * ( Wave[INDEX( ix, iy, iz )] + alpha * h ) + a1 * ( Wave[INDEX( ix, iy + 1, iz )] + alpha * h ) + a2 * ( Wave[INDEX( ix, iy + 2, iz )] + alpha * h ) + a3 * ( Wave[INDEX( ix, iy + 3, iz )] + alpha * h ) )

#define LzF( Wave, ix, iy, iz, delta, alpha, h ) 1.0f / delta * \
( a_1 * ( Wave[INDEX( ix, iy, iz - 1 )] + alpha * h ) + a0 * ( Wave[INDEX( ix, iy, iz )] + alpha * h ) + a1 * ( Wave[INDEX( ix, iy, iz + 1 )] + alpha * h ) + a2 * ( Wave[INDEX( ix, iy, iz + 2 )] + alpha * h ) + a3 * ( Wave[INDEX( ix, iy, iz + 3 )] + alpha * h ) )


#define LxB( Wave, ix, iy, iz, delta, alpha, h ) - 1.0f / delta * \
( a_1 * ( Wave[INDEX( ix + 1, iy, iz )] + alpha * h ) + a0 * ( Wave[INDEX( ix, iy, iz )] + alpha * h ) + a1 * ( Wave[INDEX( ix - 1, iy, iz )] + alpha * h ) + a2 * ( Wave[INDEX( ix - 2, iy, iz )] + alpha * h ) + a3 * ( Wave[INDEX( ix - 3, iy, iz )] + alpha * h ) )

#define LyB( Wave, ix, iy, iz, delta, alpha, h ) - 1.0f / delta * \
( a_1 * ( Wave[INDEX( ix, iy + 1, iz )] + alpha * h ) + a0 * ( Wave[INDEX( ix, iy, iz )] + alpha * h ) + a1 * ( Wave[INDEX( ix, iy - 1, iz )] + alpha * h ) + a2 * ( Wave[INDEX( ix, iy - 2, iz )] + alpha * h ) + a3 * ( Wave[INDEX( ix, iy - 3, iz )] + alpha * h ) )

#define LzB( Wave, ix, iy, iz, delta, alpha, h ) - 1.0f / delta * \
( a_1 * ( Wave[INDEX( ix, iy, iz + 1 )] + alpha * h ) + a0 * ( Wave[INDEX( ix, iy, iz )] + alpha * h ) + a1 * ( Wave[INDEX( ix, iy, iz - 1 )] + alpha * h ) + a2 * ( Wave[INDEX( ix, iy, iz - 2 )] + alpha * h ) + a3 * ( Wave[INDEX( ix, iy, iz - 3 )] + alpha * h ) )



//o means orientation that can express x,y,z 
#define L1( X, Y, Z, Wave_X, Wave_Y, Wave_Z )		Lx##X( Wave_X, ix, iy, iz, delta.x ) + Ly##Y( Wave_Y, ix, iy, iz, delta.y ) + Lz##Z( Wave_Z, ix, iy, iz, delta.z )
#define L2To4( X, Y, Z, Wave_X, Wave_Y, Wave_Z )  	
Lx##X( Wave_X, ix, iy, iz, delta.x ) + Ly##Y( Wave_Y, ix, iy, iz, delta.y ) + Lz##Z( Wave_Z, ix, iy, iz, delta.z )




#define RK4( vsX, vsY, vsZ ) \
				h1 = delta.t *  L( F, F, F, vsX, vsY, vsZ );
				h2 = 
				h1 = delta.t * L##xyz##F1( vs, ix, iy, iz, delta.xyz );\
				h2 = delta.t * L##xyz##B ( vs, ix, iy, iz, delta.xyz, alpha4_2, h1 );\
				h3 = delta.t * L##xyz##F ( vs, ix, iy, iz, delta.xyz, alpha4_3, h2 );\
				h4 = delta.t * L##xyz##F ( vs, ix, iy, iz, delta.xyz, alpha4_4, h3 );\
				diff##xyz = beta4_1 * h1 + beta4_2 * h2 + beta4_3 * h3 + beta4_4 * h4;
#define RK RK4

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

