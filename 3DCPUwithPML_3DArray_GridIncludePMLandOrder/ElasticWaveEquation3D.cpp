/***************************************************************************************/
/******************************************2018-12-21***********************************/
/************************************Author:Wenqiang Wang*******************************/
/**********************Southern University of Science and Technology********************/
/***************************************************************************************/

#ifndef ELASTIC_WAVE_EQUATION3D_HPP
#define ELASTIC_WAVE_EQUATION3D_HPP

#include "ElasticWaveEquation3D.h"

ElasticWaveEquation3D::
ElasticWaveEquation3D(
	GRID grid, int timeLength, int pml, 
	float *** buoyancy, float *** lambda, float *** mu, 
	SOURCELOCATION sourceLocation, 
	float deltaX, float deltaY, float deltaZ, float deltaT, 
	float mainFrequncy )
	: grid(grid), timeLength(timeLength), pml(pml),
	buoyancy(buoyancy), lambda(lambda), mu(mu),
	deltaX(deltaX), deltaY(deltaY), deltaZ(deltaZ), deltaT(deltaT),
	mainFrequncy(mainFrequncy), sourceLocation(sourceLocation)
{
	initialize();
}

void ElasticWaveEquation3D::initialize()
{
	ALLOCATE_3D_VelocityStress( float,   , grid.xLength + 2 * pml + ORDER, grid.yLength + 2 * pml + ORDER, grid.zLength + 2 * pml + ORDER);
	ALLOCATE_3D_VelocityStress( float, _x, grid.xLength + 2 * pml + ORDER, grid.yLength + 2 * pml + ORDER, grid.zLength + 2 * pml + ORDER);
	ALLOCATE_3D_VelocityStress( float, _y, grid.xLength + 2 * pml + ORDER, grid.yLength + 2 * pml + ORDER, grid.zLength + 2 * pml + ORDER);
	ALLOCATE_3D_VelocityStress( float, _z, grid.xLength + 2 * pml + ORDER, grid.yLength + 2 * pml + ORDER, grid.zLength + 2 * pml + ORDER);
																 
	ALLOCATE_3D(float, damping_x, grid.xLength + 2 * pml + ORDER, grid.yLength + 2 * pml + ORDER, grid.zLength + 2 * pml + ORDER);
	ALLOCATE_3D(float, damping_y, grid.xLength + 2 * pml + ORDER, grid.yLength + 2 * pml + ORDER, grid.zLength + 2 * pml + ORDER);
	ALLOCATE_3D(float, damping_z, grid.xLength + 2 * pml + ORDER, grid.yLength + 2 * pml + ORDER, grid.zLength + 2 * pml + ORDER);
}														



ElasticWaveEquation3D::~ElasticWaveEquation3D()
{
	DEF_3D_VelocityStress( float,   , grid.xLength + 2 * pml + ORDER, grid.yLength + 2 * pml + ORDER, grid.zLength + 2 * pml + ORDER);
	DEF_3D_VelocityStress( float, _x, grid.xLength + 2 * pml + ORDER, grid.yLength + 2 * pml + ORDER, grid.zLength + 2 * pml + ORDER);
	DEF_3D_VelocityStress( float, _y, grid.xLength + 2 * pml + ORDER, grid.yLength + 2 * pml + ORDER, grid.zLength + 2 * pml + ORDER);
	DEF_3D_VelocityStress( float, _z, grid.xLength + 2 * pml + ORDER, grid.yLength + 2 * pml + ORDER, grid.zLength + 2 * pml + ORDER);
																 
	DELETE_3D(float, damping_x, grid.xLength + 2 * pml + ORDER, grid.yLength + 2 * pml + ORDER, grid.zLength + 2 * pml + ORDER);
	DELETE_3D(float, damping_y, grid.xLength + 2 * pml + ORDER, grid.yLength + 2 * pml + ORDER, grid.zLength + 2 * pml + ORDER);
	DELETE_3D(float, damping_z, grid.xLength + 2 * pml + ORDER, grid.yLength + 2 * pml + ORDER, grid.zLength + 2 * pml + ORDER);
}

void ElasticWaveEquation3D::propagate()
{
	calculateDamping();
	ofstream curveLine1("curveLine1.txt");
	ofstream curveLine2("curveLine2.txt");
	ofstream curveLine3("curveLine3.txt");
	for (int it = 0; it < timeLength; it++)
	{
		cout << "Time Iterator: " << it << endl;
		loadSeismicSource(it);

		velocityUpdate();
		stressUpdate();

		curveLine1 << stressXX[pml + 10][sourceLocation.y][sourceLocation.z] << endl;
		curveLine2 << stressXX[sourceLocation.x][pml + 10][sourceLocation.z] << endl;
		curveLine3 << stressXX[sourceLocation.x][sourceLocation.y][pml + 10] << endl;
		if ((it + 1) >= 200)
		{
			if (it%50 == 0)
			{
				this->snapshoot(it + 1);
			}
		}
		if ( ( it + 1 ) == 300 )
			this->snapshoot(it + 1);
		/*if ((it + 1) == 400)
			this->snapshoot(it + 1);*/
	}
	curveLine1.close();
	curveLine2.close();
	curveLine3.close();
}

inline void ElasticWaveEquation3D::snapshoot(int it)
{
	char fileName[6][100] = { 0 };
	FOR_LOOP(i, 0, 6, sprintf(fileName[i], "snapshot_%d_%d.txt", i, it));
	FOR_LOOP(i, 0, 6, cout << fileName[i] << endl);
	ofstream * snapshot = new ofstream[6];
	FOR_LOOP(i, 0, 6, snapshot[i].open(fileName[i]));
	FOR_LOOP(i, 0, 6, if (!snapshot[i].is_open()) cout << "Can't open the file" << endl);

	FOR_LOOP_3D(i, 0, j, 0, k, 0, grid, pml, ORDER, snapshot[0] << stressXX[i][j][k] << endl;);
	/*FOR_LOOP_3D(i, 0, j, 0, k, 0, grid, pml + ORDER, snapshot[1] << stressYY[i][j][k]<< endl;);
	FOR_LOOP_3D(i, 0, j, 0, k, 0, grid, pml, snapshot[2] << stressZZ[i][j][k]<< endl;);
	FOR_LOOP_3D(i, 0, j, 0, k, 0, grid, pml, snapshot[3] << stressXY[i][j][k]<< endl;);
	FOR_LOOP_3D(i, 0, j, 0, k, 0, grid, pml, snapshot[4] << stressXZ[i][j][k]<< endl;);
	FOR_LOOP_3D(i, 0, j, 0, k, 0, grid, pml, snapshot[5] << stressYZ[i][j][k]<< endl;);
	*/
	//FOR_LOOP_3D(i, 0, j, 0, k, 0, grid, pml, snapshot[0] << velocityX[i][j][k] << endl;);
	//FOR_LOOP_3D(i, 0, j, 0, k, 0, grid, pml, snapshot[1] << velocityY[i][j][k] << endl;);
	//FOR_LOOP_3D(i, 0, j, 0, k, 0, grid, pml, snapshot[2] << velocityZ[i][j][k] << endl;);

	FOR_LOOP(i, 0, 6, snapshot[i].close(); );


	delete[] snapshot;

}

inline float ElasticWaveEquation3D::sourceFunction( int timeIndex )
{
	float t0 = ceil( 1.0f / ( mainFrequncy * deltaT ) );
	float tmp = pow( PI * mainFrequncy * ( ( timeIndex - t0 )* deltaT - 1 / mainFrequncy ), 2 );
	return ( 1 - 2 * tmp ) * exp( -tmp );
}

inline void ElasticWaveEquation3D::loadSeismicSource(int timeIterator)
{
	stressXX[sourceLocation.x][sourceLocation.y][sourceLocation.z] += pow(deltaT, 2) * sourceFunction(timeIterator);
	stressYY[sourceLocation.x][sourceLocation.y][sourceLocation.z] += pow(deltaT, 2) * sourceFunction(timeIterator);
	stressZZ[sourceLocation.x][sourceLocation.y][sourceLocation.z] += pow(deltaT, 2) * sourceFunction(timeIterator);
	//stressXY[sourceLocation.x][sourceLocation.y][sourceLocation.z] += pow( deltaT, 2 ) * sourceFunction( it );
	//stressXZ[sourceLocation.x][sourceLocation.y][sourceLocation.z] += pow( deltaT, 2 ) * sourceFunction( it );
	//stressYZ[sourceLocation.x][sourceLocation.y][sourceLocation.z] += pow( deltaT, 2 ) * sourceFunction( it );
	cout << sourceFunction(timeIterator) << endl;
}



inline void ElasticWaveEquation3D::velocityUpdate( )
{
	float dampingMinusX = 0.0f, dampingPlusX = 0.0f;
	float dampingMinusY = 0.0f, dampingPlusY = 0.0f;
	float dampingMinusZ = 0.0f, dampingPlusZ = 0.0f;
	float daltaT_buoyancyBar = 0.0f;
	for (int i = ORDER / 2; i < grid.xLength + 2 * pml + ORDER / 2; i++)
	{
		for (int j = ORDER / 2; j < grid.yLength + 2 * pml + ORDER / 2; j++)
		{
			for (int k = ORDER / 2; k < grid.zLength + 2 * pml + ORDER / 2; k++)
			{
				dampingMinusX = 1.0f - 0.5f * deltaT * damping_x[i][j][k];
				dampingMinusY = 1.0f - 0.5f * deltaT * damping_y[i][j][k];
				dampingMinusZ = 1.0f - 0.5f * deltaT * damping_z[i][j][k];
				dampingPlusX  = 1.0f + 0.5f * deltaT * damping_x[i][j][k];
				dampingPlusY  = 1.0f + 0.5f * deltaT * damping_y[i][j][k];
				dampingPlusZ  = 1.0f + 0.5f * deltaT * damping_z[i][j][k];
				daltaT_buoyancyBar = deltaT * BUOYANCYBARX( buoyancy, i, j, k );
				velocityX_x[i][j][k] = ( dampingMinusX * velocityX_x[i][j][k] +
										daltaT_buoyancyBar * DIFF_X( stressXX, i+1, j  , k  , deltaX ) ) / dampingPlusX;
				velocityX_y[i][j][k] = ( dampingMinusY * velocityX_y[i][j][k] +
										daltaT_buoyancyBar * DIFF_Y( stressXY, i  , j  , k  , deltaY ) ) / dampingPlusY;
				velocityX_z[i][j][k] = ( dampingMinusZ * velocityX_z[i][j][k] +
										daltaT_buoyancyBar * DIFF_Z( stressXZ, i  , j  , k  , deltaZ ) ) / dampingPlusZ;
				velocityX[i][j][k] = velocityX_x[i][j][k] + velocityX_y[i][j][k] + velocityX_z[i][j][k];
				
				daltaT_buoyancyBar = deltaT * BUOYANCYBARY( buoyancy, i, j, k );
				velocityY_x[i][j][k] = ( dampingMinusX * velocityY_x[i][j][k] +
										daltaT_buoyancyBar * DIFF_X( stressXY, i  , j  , k  , deltaX ) ) / dampingPlusX;
				velocityY_y[i][j][k] = ( dampingMinusY * velocityY_y[i][j][k] +
										daltaT_buoyancyBar * DIFF_Y( stressYY, i  , j+1, k  , deltaY ) ) / dampingPlusY;
				velocityY_z[i][j][k] = ( dampingMinusZ * velocityY_z[i][j][k] +
										daltaT_buoyancyBar * DIFF_Z( stressYZ, i  , j  , k  , deltaZ ) ) / dampingPlusZ;
				velocityY[i][j][k] = velocityY_x[i][j][k] + velocityY_y[i][j][k] + velocityY_z[i][j][k];
				
				daltaT_buoyancyBar = deltaT * BUOYANCYBARZ( buoyancy, i, j, k );
				velocityZ_x[i][j][k] = ( dampingMinusX * velocityZ_x[i][j][k] +
										daltaT_buoyancyBar * DIFF_X( stressXZ, i  , j  , k  , deltaX ) ) / dampingPlusX;
				velocityZ_y[i][j][k] = ( dampingMinusY * velocityZ_y[i][j][k] +
										daltaT_buoyancyBar * DIFF_Y( stressYZ, i  , j  , k  , deltaY ) ) / dampingPlusY;
				velocityZ_z[i][j][k] = ( dampingMinusZ * velocityZ_z[i][j][k] +
										daltaT_buoyancyBar * DIFF_Z( stressZZ, i  , j  , k+1, deltaZ ) ) / dampingPlusZ;
				velocityZ[i][j][k] = velocityZ_x[i][j][k] + velocityZ_y[i][j][k] + velocityZ_z[i][j][k];
				
			}
		}
	}

}

inline void ElasticWaveEquation3D::stressUpdate()
{
	float tmp = 0.0f, tmpDiffX = 0.0f, tmpDiffY = 0.0f, tmpDiffZ = 0.0f;
	float tmpDiff = 0.0f, tmpLambdaDiff = 0.0f, muBarXY = 0.0f, muBarXZ = 0.0f, muBarYZ = 0.0f;
	float dampingMinus = 0.0f, dampingPlus = 0.0f;

	float dampingMinusX = 0.0f, dampingPlusX = 0.0f;
	float dampingMinusY = 0.0f, dampingPlusY = 0.0f;
	float dampingMinusZ = 0.0f, dampingPlusZ = 0.0f;

	for (int i = ORDER / 2; i < grid.xLength + 2 * pml + ORDER / 2; i++)
	{
		for (int j = ORDER / 2; j < grid.yLength + 2 * pml + ORDER / 2; j++)
		{
			for (int k = ORDER / 2; k < grid.zLength + 2 * pml + ORDER / 2; k++)
			{
				dampingMinusX = 1.0f - 0.5f * deltaT * damping_x[i][j][k];
				dampingMinusY = 1.0f - 0.5f * deltaT * damping_y[i][j][k];
				dampingMinusZ = 1.0f - 0.5f * deltaT * damping_z[i][j][k];
				dampingPlusX  = 1.0f + 0.5f * deltaT * damping_x[i][j][k];
				dampingPlusY  = 1.0f + 0.5f * deltaT * damping_y[i][j][k];
				dampingPlusZ  = 1.0f + 0.5f * deltaT * damping_z[i][j][k];
				
				tmpDiffX = DIFF_X(velocityX, i, j, k, deltaX);
				tmpDiffY = DIFF_Y(velocityY, i, j, k, deltaY);
				tmpDiffZ = DIFF_Z(velocityZ, i, j, k, deltaZ);
				
				stressXX_x[i][j][k] = (dampingMinusX * stressXX_x[i][j][k] + deltaT * (lambda[i][j][k] + 2.0f * mu[i][j][k]) * tmpDiffX) / dampingPlusX;
				stressXX_y[i][j][k] = (dampingMinusY * stressXX_y[i][j][k] + deltaT * lambda[i][j][k] * tmpDiffY) / dampingPlusY;
				stressXX_z[i][j][k] = (dampingMinusZ * stressXX_z[i][j][k] + deltaT * lambda[i][j][k] * tmpDiffZ) / dampingPlusZ;
				stressXX[i][j][k] = stressXX_x[i][j][k] + stressXX_y[i][j][k] + stressXX_z[i][j][k];
				
				stressYY_x[i][j][k] = (dampingMinusX * stressYY_x[i][j][k] + deltaT * lambda[i][j][k] * tmpDiffX) / dampingPlusX;
				stressYY_y[i][j][k] = (dampingMinusY * stressYY_y[i][j][k] + deltaT * (lambda[i][j][k] + 2.0f * mu[i][j][k]) * tmpDiffY) / dampingPlusY;
				stressYY_z[i][j][k] = (dampingMinusZ * stressYY_z[i][j][k] + deltaT * lambda[i][j][k] * tmpDiffZ) / dampingPlusZ;
				stressYY[i][j][k] = stressYY_x[i][j][k] + stressYY_y[i][j][k] + stressYY_z[i][j][k];
				
				stressZZ_x[i][j][k] = (dampingMinusX * stressZZ_x[i][j][k] + deltaT * lambda[i][j][k] * tmpDiffX) / dampingPlusX;
				stressZZ_y[i][j][k] = (dampingMinusY * stressZZ_y[i][j][k] + deltaT * lambda[i][j][k] * tmpDiffY) / dampingPlusY;
				stressZZ_z[i][j][k] = (dampingMinusZ * stressZZ_z[i][j][k] + deltaT * (lambda[i][j][k] + 2.0f * mu[i][j][k]) * tmpDiffZ) / dampingPlusZ;
				stressZZ[i][j][k] = stressZZ_x[i][j][k] + stressZZ_y[i][j][k] + stressZZ_z[i][j][k];
				
				muBarXY = MUBARXY(mu, i, j, k);
				muBarXZ = MUBARXZ(mu, i, j, k);
				muBarYZ = MUBARYZ(mu, i, j, k);
				
				stressXY_x[i][j][k] = (dampingMinusX * stressXY_x[i][j][k] + deltaT * muBarXY * DIFF_X(velocityY, i + 1, j, k, deltaX)) / dampingPlusX;
				stressXY_y[i][j][k] = (dampingMinusY * stressXY_y[i][j][k] + deltaT * muBarXY * DIFF_Y(velocityX, i, j + 1, k, deltaY)) / dampingPlusY;
				stressXY_z[i][j][k] = dampingMinusZ * stressXY_z[i][j][k] / dampingPlusZ;
				stressXY[i][j][k] = stressXY_x[i][j][k] + stressXY_y[i][j][k] + stressXY_z[i][j][k];
				
				stressXZ_x[i][j][k] = (dampingMinusX * stressXZ_x[i][j][k] + deltaT * muBarXZ * DIFF_X(velocityZ, i + 1, j, k, deltaX)) / dampingPlusX;
				stressXZ_y[i][j][k] = dampingMinusY * stressXZ_y[i][j][k] / dampingPlusY;
				stressXZ_z[i][j][k] = (dampingMinusZ * stressXZ_z[i][j][k] + deltaT * muBarXZ * DIFF_Z(velocityX, i, j, k + 1, deltaZ)) / dampingPlusZ;
				stressXZ[i][j][k] = stressXZ_x[i][j][k] + stressXZ_y[i][j][k] + stressXZ_z[i][j][k];
				
				stressYZ_x[i][j][k] = dampingMinusX * stressYZ_x[i][j][k] / dampingPlusX;
				stressYZ_y[i][j][k] = (dampingMinusY * stressYZ_y[i][j][k] + deltaT * muBarYZ * DIFF_Y(velocityZ, i, j + 1, k, deltaY)) / dampingPlusY;
				stressYZ_z[i][j][k] = (dampingMinusZ * stressYZ_z[i][j][k] + deltaT * muBarYZ * DIFF_Z(velocityY, i, j, k + 1, deltaZ)) / dampingPlusZ;
				stressYZ[i][j][k] = stressYZ_x[i][j][k] + stressYZ_y[i][j][k] + stressYZ_z[i][j][k];
				
			}
		}
	}
}



inline void ElasticWaveEquation3D::calculateDamping( )
{
	float damping0X = 0.0f, damping0Y = 0.0f, damping0Z = 0.0f;
	damping0X = Vs / deltaX * ( 8.0f / 15.0f - 3.0f / 100.0f * pml + 1.0f / 1500.0f * pml * pml );
	damping0Y = Vs / deltaY * ( 8.0f / 15.0f - 3.0f / 100.0f * pml + 1.0f / 1500.0f * pml * pml );
	damping0Z = Vs / deltaZ * ( 8.0f / 15.0f - 3.0f / 100.0f * pml + 1.0f / 1500.0f * pml * pml );
	for (int i = ORDER / 2; i < grid.xLength + 2 * pml + ORDER / 2; i++)
	{
		for (int j = ORDER / 2; j < grid.yLength + 2 * pml + ORDER / 2; j++)
		{
			for (int k = ORDER / 2; k < grid.zLength + 2 * pml + ORDER / 2; k++)
			{
				if (i >= ORDER / 2 + pml && i < grid.xLength + pml + ORDER / 2 && 
					j >= ORDER / 2 + pml && j < grid.yLength + pml + ORDER / 2 && 
					k >= ORDER / 2 + pml && k < grid.zLength + pml + ORDER / 2 )
					;
				else
				{
					if ( i < pml + ORDER / 2 )
					{
						damping_x[i][j][k] = DAMPING(Vs, pml + ORDER / 2 - i, pml, deltaX);
					}
					else if (i >= grid.xLength + pml + ORDER / 2)
					{
						damping_x[i][j][k] = DAMPING(Vs, grid.xLength + pml + ORDER / 2 - i - 1, pml, deltaX);
					}

					if (j < pml + ORDER / 2 )
					{
						damping_y[i][j][k] = DAMPING(Vs, pml + ORDER / 2 - j, pml, deltaY);
					}
					else if (j >= grid.yLength + pml + ORDER / 2 )
					{
						damping_y[i][j][k] = DAMPING(Vs, grid.yLength + pml + ORDER / 2 - j - 1, pml, deltaY);
					}

					if (k < pml + ORDER / 2 )
					{
						damping_z[i][j][k] = DAMPING(Vs, pml + ORDER / 2 - k, pml, deltaZ);
					}
					else if (k >= grid.zLength + pml + ORDER / 2)
					{
						damping_z[i][j][k] = DAMPING(Vs, grid.zLength + pml + ORDER / 2 - k - 1, pml, deltaZ);
					}
				}
			}
		}
	}

	char fileName[3][100] = { "dampingX.txt", "dampingY.txt", "dampingZ.txt" };
	ofstream * dampingData = new ofstream[3];
	FOR_LOOP(i, 0, 3, dampingData[i].open(fileName[i]));
	FOR_LOOP(i, 0, 3, if (!dampingData[i].is_open()) cout << "Can't open the file" << endl);

	FOR_LOOP_3D(i, 0, j, 0, k, 0, grid, pml, ORDER, dampingData[0] << damping_x[i][j][k] << endl;);
	FOR_LOOP_3D(i, 0, j, 0, k, 0, grid, pml, ORDER, dampingData[1] << damping_y[i][j][k] << endl;);
	FOR_LOOP_3D(i, 0, j, 0, k, 0, grid, pml, ORDER, dampingData[2] << damping_z[i][j][k] << endl;);
	FOR_LOOP(i, 0, 3, dampingData[i].close());

}

#endif 