#ifndef ELASTIC_MACRO_HPP
#define ELASTIC_MACRO_HPP
/***************************************************************************************/
/******************************************2018-12-21***********************************/
/************************************Author:Wenqiang Wang*******************************/
/**********************Southern University of Science and Technology********************/
/***************************************************************************************/
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
using namespace std;
#include <memory.h>
#include <cstdio>
#include <cmath>
#include <ctime>
#define PI 3.1415926f
#define R 0.001f
#define tau 3.4f
#define Vs 2000.0f
#define ORDER 10
#define DAMPING0( V, pml, delta ) log10( 1.0f / R ) * tau * V / ( pml * delta )
//#define DAMPING0( V, pml, delta ) tau * V / ( ( pml ) * delta ) * ( 8.0f / 15.0f - 0.03f * ( pml ) + 1.0f / 1500.0f * ( ( pml ) * ( pml ) ) )
#define DAMPING( V, i, pml, delta ) DAMPING0( V, pml, delta ) * powf( ( i ) * 1.0f / ( ( pml ) * 1.0f ), 2 )


#if ORDER == 12
#define DIFF_12( phiF, phiB1, phiF1, phiB2, phiF2, phiB3, phiF3, phiB4, phiF4, phiB5, phiF5, phiB6, delta ) \
				( (   1.22134f * ( phiF  - phiB1 ) \
				  -9.69315e-2f * ( phiF1 - phiB2 ) \
				   1.74477e-2f * ( phiF2 - phiB3 ) \
				  -2.96729e-3f * ( phiF3 - phiB4 ) \
				   3.59005e-4f * ( phiF4 - phiB5 ) \
				  -2.18478e-5f * ( phiF5 - phiB6 ) \
				) / delta )

#define DIFF_X(  phi, i, j, k, deltaX  ) \
		DIFF_12( phi[i    ][j][k], phi[i - 1][j][k],\
				 phi[i + 1][j][k], phi[i - 2][j][k],\
				 phi[i + 2][j][k], phi[i - 3][j][k],\
				 phi[i + 3][j][k], phi[i - 4][j][k],\
				 phi[i + 4][j][k], phi[i - 5][j][k],\
				 phi[i + 5][j][k], phi[i - 6][j][k], deltaX )

#define DIFF_Y(  phi, i, j, k, deltaY  ) \
		DIFF_12( phi[i][j    ][k], phi[i][j - 1][k],\
				 phi[i][j + 1][k], phi[i][j - 2][k],\
				 phi[i][j + 2][k], phi[i][j - 3][k],\
				 phi[i][j + 3][k], phi[i][j - 4][k],\
				 phi[i][j + 4][k], phi[i][j - 5][k],\
				 phi[i][j + 5][k], phi[i][j - 6][k], deltaY )

#define DIFF_Z(  phi, i, j, k, deltaZ  ) \
		DIFF_12( phi[i][j][k    ], phi[i][j][k - 1],\
				 phi[i][j][k + 1], phi[i][j][k - 2],\
				 phi[i][j][k + 2], phi[i][j][k - 3],\
				 phi[i][j][k + 3], phi[i][j][k - 4],\
				 phi[i][j][k + 4], phi[i][j][k - 5],\
				 phi[i][j][k + 5], phi[i][j][k - 6], deltaZ )


#endif // ORDER == 12

#if ORDER == 10
#define DIFF_10( phiF, phiB1, phiF1, phiB2, phiF2, phiB3, phiF3, phiB4, phiF4, phiB5, delta ) \
				( (   1.21124f * ( phiF  - phiB1 ) \
				  -8.97217e-2f * ( phiF1 - phiB2 ) \
				  +1.38428e-2f * ( phiF2 - phiB3 ) \
				  -1.76566e-3f * ( phiF3 - phiB4 ) \
				  +1.18680e-4f * ( phiF4 - phiB5 ) \
				) / delta )
#define DIFF_X(  phi, i, j, k, deltaX  ) \
		DIFF_10( phi[i    ][j][k], phi[i - 1][j][k],\
				 phi[i + 1][j][k], phi[i - 2][j][k],\
				 phi[i + 2][j][k], phi[i - 3][j][k],\
				 phi[i + 3][j][k], phi[i - 4][j][k],\
				 phi[i + 4][j][k], phi[i - 5][j][k], deltaX )

#define DIFF_Y(  phi, i, j, k, deltaY  ) \
		DIFF_10( phi[i][j    ][k], phi[i][j - 1][k],\
				 phi[i][j + 1][k], phi[i][j - 2][k],\
				 phi[i][j + 2][k], phi[i][j - 3][k],\
				 phi[i][j + 3][k], phi[i][j - 4][k],\
				 phi[i][j + 4][k], phi[i][j - 5][k], deltaY )

#define DIFF_Z(  phi, i, j, k, deltaZ  ) \
		DIFF_10( phi[i][j][k    ], phi[i][j][k - 1],\
				 phi[i][j][k + 1], phi[i][j][k - 2],\
				 phi[i][j][k + 2], phi[i][j][k - 3],\
				 phi[i][j][k + 3], phi[i][j][k - 4],\
				 phi[i][j][k + 4], phi[i][j][k - 5], deltaZ )
#endif // ORDER == 10

#define FOR_LOOP( i, start, N, expression )\
		for( int i = (start); i < (N); i ++ )\
		{ \
			expression;\
		}
#define FOR_LOOP_3D( i, startI, j, startJ, k, startK, grid, pml, order, expression ) \
		FOR_LOOP( i, startI, grid.xLength + 2 * pml + order, FOR_LOOP( j, startJ, grid.yLength + 2 * pml + order, FOR_LOOP( k, startK, grid.zLength + 2 * pml + order, expression ) ) )

#define BUOYANCYBAR( buoyancy1, buoyancy2 ) ( ( buoyancy1 + buoyancy2 ) / 2 ) 
#define BUOYANCYBARX( buoyancy, i, j, k ) BUOYANCYBAR( buoyancy[i][j][k], buoyancy[i + 1][j][k] ) 
#define BUOYANCYBARY( buoyancy, i, j, k ) BUOYANCYBAR( buoyancy[i][j][k], buoyancy[i][j + 1][k] ) 
#define BUOYANCYBARZ( buoyancy, i, j, k ) BUOYANCYBAR( buoyancy[i][j][k], buoyancy[i][j][k + 1] ) 

#define MUBAR( mu1, mu2, mu3, mu4 ) \
		( 4.0f / ( 1 / mu1 + 1 / mu2 + 1 / mu3 + 1 / mu4 ) )

#define MUBARXY( mu, i, j, k ) \
		MUBAR( mu[i][j][k], mu[i + 1][j][k], mu[i][j + 1][k], mu[i + 1][j + 1][k] )
#define MUBARXZ( mu, i, j, k ) \
		MUBAR( mu[i][j][k], mu[i + 1][j][k], mu[i][j][k + 1], mu[i + 1][j][k + 1] )
#define MUBARYZ( mu, i, j, k ) \
		MUBAR( mu[i][j][k], mu[i][j + 1][k], mu[i][j][k + 1], mu[i][j + 1][k + 1] )

#define DEF_3D_POINTER( type, pointer3D )  type *** pointer3D;
#define ALLOCATE_3D( type, pointer3D, xLength, yLength, zLength ) \
		pointer3D = new type ** [xLength]; \
		for (int i = 0; i < xLength; i++) \
		{ \
			pointer3D[i] = new type * [yLength]; \
			for (int j = 0; j < yLength; j++) \
			{ \
				pointer3D[i][j]  = new float [zLength]; \
				memset( pointer3D[i][j], 0, sizeof( float ) * (zLength));\
			} \
		}
#define DELETE_3D( type, pointer3D, xLength, yLength, zLength ) \
		for (int i = 0; i < xLength; i++) \
		{ \
			for (int j = 0; j < yLength; j++) \
			{ \
				delete[] pointer3D[i][j]; \
			} \
			delete[] pointer3D[i]; \
		}\
		delete[] pointer3D;


#define DEF_3D_VelocityStress( type, name ) \
		DEF_3D_POINTER(type, stressXX##name )\
		DEF_3D_POINTER(type, stressYY##name )\
		DEF_3D_POINTER(type, stressZZ##name )\
		DEF_3D_POINTER(type, stressXY##name )\
		DEF_3D_POINTER(type, stressXZ##name )\
		DEF_3D_POINTER(type, stressYZ##name )\
		DEF_3D_POINTER(type, velocityX##name)\
		DEF_3D_POINTER(type, velocityY##name)\
		DEF_3D_POINTER(type, velocityZ##name)

#define ALLOCATE_3D_VelocityStress( type, name, xLength, yLength, zLength )\
		stressXX##name  = new type **[xLength]; \
		stressYY##name  = new type **[xLength]; \
		stressZZ##name  = new type **[xLength]; \
		stressXY##name  = new type **[xLength]; \
		stressXZ##name  = new type **[xLength]; \
		stressYZ##name  = new type **[xLength]; \
		velocityX##name = new type **[xLength]; \
		velocityY##name = new type **[xLength]; \
		velocityZ##name = new type **[xLength]; \
		for (int i = 0; i < xLength; i++) \
		{ \
			stressXX##name [i] = new type *[yLength]; \
			stressYY##name [i] = new type *[yLength]; \
			stressZZ##name [i] = new type *[yLength]; \
			stressXY##name [i] = new type *[yLength]; \
			stressXZ##name [i] = new type *[yLength]; \
			stressYZ##name [i] = new type *[yLength]; \
			velocityX##name[i] = new type *[yLength]; \
			velocityY##name[i] = new type *[yLength]; \
			velocityZ##name[i] = new type *[yLength]; \
			for (int j = 0; j < yLength; j++) \
			{ \
				stressXX##name [i][j] = new float[zLength]; \
				stressYY##name [i][j] = new float[zLength]; \
				stressZZ##name [i][j] = new float[zLength]; \
				stressXY##name [i][j] = new float[zLength]; \
				stressXZ##name [i][j] = new float[zLength]; \
				stressYZ##name [i][j] = new float[zLength]; \
				velocityX##name[i][j] = new float[zLength]; \
				velocityY##name[i][j] = new float[zLength]; \
				velocityZ##name[i][j] = new float[zLength]; \
				memset( stressXX##name [i][j], 0, sizeof(float) * (zLength) );\
				memset( stressYY##name [i][j], 0, sizeof(float) * (zLength) );\
				memset( stressZZ##name [i][j], 0, sizeof(float) * (zLength) );\
				memset( stressXY##name [i][j], 0, sizeof(float) * (zLength) );\
				memset( stressXZ##name [i][j], 0, sizeof(float) * (zLength) );\
				memset( stressYZ##name [i][j], 0, sizeof(float) * (zLength) );\
				memset( velocityX##name[i][j], 0, sizeof(float) * (zLength) );\
				memset( velocityY##name[i][j], 0, sizeof(float) * (zLength) );\
				memset( velocityZ##name[i][j], 0, sizeof(float) * (zLength) );\
			} \
		}


#define DELETE_3D_VelocityStress( type, name, xLength, yLength )\
		for (int i = 0; i < xLength; i++) \
		{ \
			for (int j = 0; j < yLength; j++) \
			{ \
				delete[] stressXX##name [i][j]; \
				delete[] stressYY##name [i][j]; \
				delete[] stressZZ##name [i][j]; \
				delete[] stressXY##name [i][j]; \
				delete[] stressXZ##name [i][j]; \
				delete[] stressYZ##name [i][j]; \
				delete[] velocityX##name[i][j]; \
				delete[] velocityY##name[i][j]; \
				delete[] velocityZ##name[i][j]; \
			} \
			delete[] stressXX##name [i];\
			delete[] stressYY##name [i];\
			delete[] stressZZ##name [i];\
			delete[] stressXY##name [i];\
			delete[] stressXZ##name [i];\
			delete[] stressYZ##name [i];\
			delete[] velocityX##name[i];\
			delete[] velocityY##name[i];\
			delete[] velocityZ##name[i];\
		}\
		delete[] stressXX##name ; \
		delete[] stressYY##name ; \
		delete[] stressZZ##name ; \
		delete[] stressXY##name ; \
		delete[] stressXZ##name ; \
		delete[] stressYZ##name ; \
		delete[] velocityX##name; \
		delete[] velocityY##name; \
		delete[] velocityZ##name; 

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


#endif // !ELASTIC_MACRO_HPP

