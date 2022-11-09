/***************************************************************************************/
/******************************************2019-01-04***********************************/
/************************************Author:Wenqiang Wang*******************************/
/**********************Southern University of Science and Technology********************/
/***************************************************************************************/
#include "ElasticWaveEquation.h"

ElasticWaveEquation::ElasticWaveEquation( 
	int timeLength, GRID grid, Delta delta, SOURCELOCATION sourceLocation, 
	float mainFrequncy, float * diff_coef, 
	float * buoyancy, float * lambda, float * mu, 
	int order)
	: timeLength(timeLength), 
	grid(grid), delta(delta), sourceLocation(sourceLocation), 
	mainFrequncy(mainFrequncy), diff_coef(diff_coef), order(order), 
	buoyancy(buoyancy), lambda(lambda), mu(mu)
{
	ALLOCATE_3D_VelocityStress(float, grid.xLength, grid.yLength, grid.zLength);
	
}

ElasticWaveEquation::~ElasticWaveEquation()
{
	DELETE_3D_VelocityStress(float);
}

void ElasticWaveEquation::run()
{
	for ( int it = 0; it < timeLength; it++)
	{
		cout << "Time Iterator: " << it << endl;
		loadSeismicSource(it);
		propagate();
		/*if ((it + 1) >= 200)
		{
			if (it % 50 == 0)
			{
				this->snapshoot(it + 1);
			}
		}*/
		if ( it % 5 == 0 )
		{
			snapShootSlice( it );
		}
	}
}

void ElasticWaveEquation::propagate()
{
	int ix, iy, iz, i;
	float diffX = 0.0f, diffY = 0.0f, diffZ = 0.0f;
	float muBarXY = 0.0f, muBarXZ = 0.0f, muBarYZ = 0.0f;

	for (iz = order / 2; iz < grid.zLength - order / 2; iz++)
	{
		for (iy = order / 2; iy < grid.yLength - order / 2; iy++)
		{
			for (ix = order / 2; ix < grid.xLength - order / 2; ix++)
			{
				diffX = 0.0f, diffY = 0.0f, diffZ = 0.0f;
				for ( i = 0;  i < order / 2;  i++)
				{
					diffX += diff_coef[i] * (stressXX[INDEX(ix + 1 + i, iy, iz)] - stressXX[INDEX(ix - i, iy, iz)]);
					diffY += diff_coef[i] * (stressXY[INDEX(ix, iy + i, iz)] - stressXY[INDEX(ix, iy - 1 - i, iz)]);
					diffZ += diff_coef[i] * (stressXZ[INDEX(ix, iy, iz + i)] - stressXZ[INDEX(ix, iy, iz - 1 - i)]);
				}
				diffX /= delta.x;
				diffY /= delta.y;
				diffZ /= delta.z;
				velocityX[INDEX(ix, iy, iz)] += delta.t * BUOYANCYBARX(buoyancy, ix, iy, iz) * (diffX + diffY + diffZ);
				
				diffX = 0.0f, diffY = 0.0f, diffZ = 0.0f;
				for (i = 0; i < order / 2; i++)
				{
					diffX += diff_coef[i] * (stressXY[INDEX(ix + i, iy, iz)] - stressXY[INDEX(ix - 1 - i, iy, iz)]);
					diffY += diff_coef[i] * (stressYY[INDEX(ix, iy + 1 + i, iz)] - stressYY[INDEX(ix, iy - i, iz)]);
					diffZ += diff_coef[i] * (stressYZ[INDEX(ix, iy, iz + i)] - stressYZ[INDEX(ix, iy, iz - 1 - i)]);
				}
				diffX /= delta.x;
				diffY /= delta.y;
				diffZ /= delta.z;
				velocityY[INDEX(ix, iy, iz)] += delta.t  * BUOYANCYBARY(buoyancy, ix, iy, iz) * (diffX + diffY + diffZ);
				
				diffX = 0.0f, diffY = 0.0f, diffZ = 0.0f;
				for (i = 0; i < order / 2; i++)
				{
					diffX += diff_coef[i] * (stressXZ[INDEX(ix + i, iy, iz)] - stressXZ[INDEX(ix - 1 - i, iy, iz)]);
					diffY += diff_coef[i] * (stressYZ[INDEX(ix, iy + i, iz)] - stressYZ[INDEX(ix, iy - 1 - i, iz)]);
					diffZ += diff_coef[i] * (stressZZ[INDEX(ix, iy, iz + 1 + i)] - stressZZ[INDEX(ix, iy, iz - i)]);
				}
				diffX /= delta.x;
				diffY /= delta.y;
				diffZ /= delta.z;
				velocityZ[INDEX(ix, iy, iz)] += delta.t  * BUOYANCYBARZ(buoyancy, ix, iy, iz) * (diffX + diffY + diffZ);

				diffX = 0.0f, diffY = 0.0f, diffZ = 0.0f;
				/*if (ix == sourceLocation.x && iy == sourceLocation.y && iz == sourceLocation.z)
				{
					cout << "velocityX:" << velocityX[INDEX(ix, iy, iz)] << endl;
					cout << "velocityY:" << velocityY[INDEX(ix, iy, iz)] << endl;
					cout << "velocityZ:" << velocityZ[INDEX(ix, iy, iz)] << endl;
				}*/
			}
		}
	}

	for (iz = order / 2; iz < grid.zLength - order / 2; iz++)
	{
		for (iy = order / 2; iy < grid.yLength - order / 2; iy++)
		{
			for (ix = order / 2; ix < grid.xLength - order / 2; ix++)
			{
				diffX = 0.0f, diffY = 0.0f, diffZ = 0.0f;
				for (i = 0; i < order / 2; i++)
				{
					diffX += diff_coef[i] * (velocityX[INDEX(ix + i, iy, iz)] - velocityX[INDEX(ix - 1 - i, iy, iz)]);
					diffY += diff_coef[i] * (velocityY[INDEX(ix, iy + i, iz)] - velocityY[INDEX(ix, iy - 1 - i, iz)]);
					diffZ += diff_coef[i] * (velocityZ[INDEX(ix, iy, iz + i)] - velocityZ[INDEX(ix, iy, iz - 1 - i)]);
				}
				diffX /= delta.x;
				diffY /= delta.y;
				diffZ /= delta.z;
				stressXX[INDEX(ix, iy, iz)] += delta.t * 
					(( lambda[INDEX(ix, iy, iz)] + 2.0f * mu[INDEX(ix, iy, iz)]) * diffX +
					   lambda[INDEX(ix, iy, iz)] * (diffY + diffZ) );				 
				stressYY[INDEX(ix, iy, iz)] += delta.t * 								 
					(( lambda[INDEX(ix, iy, iz)] + 2.0f * mu[INDEX(ix, iy, iz)]) * diffY +
					   lambda[INDEX(ix, iy, iz)] * (diffX + diffZ) );					 
				stressZZ[INDEX(ix, iy, iz)] += delta.t * 								 
					(( lambda[INDEX(ix, iy, iz)] + 2.0f * mu[INDEX(ix, iy, iz)]) * diffZ +
					   lambda[INDEX(ix, iy, iz)] * (diffX + diffY) );
				

				/*if (ix == sourceLocation.x && iy == sourceLocation.y && iz == sourceLocation.z)
				{
					cout << "stressXX:" << stressXX[INDEX(ix, iy, iz)] << endl;
					cout << "stressYY:" << stressYY[INDEX(ix, iy, iz)] << endl;
					cout << "stressZZ:" << stressZZ[INDEX(ix, iy, iz)] << endl;
				}*/


				muBarXY = MUBARXY(mu, ix, iy, iz);
				muBarXZ = MUBARXZ(mu, ix, iy, iz);
				muBarYZ = MUBARYZ(mu, ix, iy, iz);

				diffX = 0.0f, diffY = 0.0f;
				for (i = 0; i < order / 2; i++)
				{
					diffY += diff_coef[i] * (velocityX[INDEX(ix, iy + 1 + i, iz)] - velocityX[INDEX(ix, iy - i, iz)]);
					diffX += diff_coef[i] * (velocityY[INDEX(ix + 1 + i, iy, iz)] - velocityY[INDEX(ix - i, iy, iz)]);
				}
				diffX /= delta.x;
				diffY /= delta.y;
				stressXY[INDEX(ix, iy, iz)] += delta.t * muBarXY * (diffY + diffX);

				diffX = 0.0f, diffZ = 0.0f;
				for (i = 0; i < order / 2; i++)
				{
					diffZ += diff_coef[i] * (velocityX[INDEX(ix, iy, iz + 1 + i)] - velocityX[INDEX(ix, iy, iz - i)]);
					diffX += diff_coef[i] * (velocityZ[INDEX(ix + 1 + i, iy, iz)] - velocityZ[INDEX(ix - i, iy, iz)]);
				}
				diffX /= delta.x;
				diffZ /= delta.z;
				stressXZ[INDEX(ix, iy, iz)] += delta.t * muBarXZ * (diffZ + diffX);

				diffY = 0.0f, diffZ = 0.0f;
				for (i = 0; i < order / 2; i++)
				{
					diffZ += diff_coef[i] * (velocityY[INDEX(ix, iy, iz + 1 + i)] - velocityY[INDEX(ix, iy, iz - i)]);
					diffY += diff_coef[i] * (velocityZ[INDEX(ix, iy + 1 + i, iz)] - velocityZ[INDEX(ix, iy - i, iz)]);
				}
				diffY /= delta.y;
				diffZ /= delta.z;
				stressYZ[INDEX(ix, iy, iz)] += delta.t * muBarYZ * (diffZ + diffY);

				/*if (ix == sourceLocation.x && iy == sourceLocation.y && iz == sourceLocation.z)
				{
					cout << "stressXX:" << stressXY[INDEX(ix, iy, iz)] << endl;
					cout << "stressYY:" << stressXZ[INDEX(ix, iy, iz)] << endl;
					cout << "stressZZ:" << stressXY[INDEX(ix, iy, iz)] << endl;
				}*/

			}
		}
	}
}

float ElasticWaveEquation::sourceFunction(int timeIndex)
{
	float t0 = ceil(1.0f / (mainFrequncy * delta.t));
	float tmp = pow(PI * mainFrequncy * ((timeIndex - t0) * delta.t - 1 / mainFrequncy), 2);
	return (1 - 2 * tmp) * exp(-tmp);
}

void ElasticWaveEquation::loadSeismicSource(int timeIterator)
{
	stressXX[INDEX(sourceLocation.x, sourceLocation.y, sourceLocation.z)] += pow(delta.t, 2) * sourceFunction(timeIterator);
	stressYY[INDEX(sourceLocation.x, sourceLocation.y, sourceLocation.z)] += pow(delta.t, 2) * sourceFunction(timeIterator);
	stressZZ[INDEX(sourceLocation.x, sourceLocation.y, sourceLocation.z)] += pow(delta.t, 2) * sourceFunction(timeIterator);
	
	cout << sourceFunction(timeIterator) << endl;
}



void ElasticWaveEquation::snapShootSlice(int it)
{
	char fileName[9][100] = { 0 };
	ofstream * snapshot = new ofstream[9];

	//sprintf(fileName[0], "snapshotSliceXX_YoZ_%d.txt", it);
	//sprintf(fileName[1], "snapshotSliceXX_XoZ_%d.txt", it);
	//sprintf(fileName[2], "snapshotSliceXX_XOY_%d.txt", it);
	sprintf(fileName[0], "velocityXSliceXX_YoZ_%d.txt", it);
	sprintf(fileName[1], "velocityXSliceXX_XoZ_%d.txt", it);
	sprintf(fileName[2], "velocityXSliceXX_XoY_%d.txt", it);

	sprintf(fileName[3], "stressXXSliceXX_YoZ_%d.txt", it);
	sprintf(fileName[4], "stressXXSliceXX_XoZ_%d.txt", it);
	sprintf(fileName[5], "stressXXSliceXX_XoY_%d.txt", it);

	sprintf(fileName[6], "stressXYSliceXX_YoZ_%d.txt", it);
	sprintf(fileName[7], "stressXYSliceXX_XoZ_%d.txt", it);
	sprintf(fileName[8], "stressXYSliceXX_XoY_%d.txt", it);

	snapshot[0].open(fileName[0]);
	snapshot[1].open(fileName[1]);
	snapshot[2].open(fileName[2]);

	snapshot[3].open(fileName[3]);
	snapshot[4].open(fileName[4]);
	snapshot[5].open(fileName[5]);

	snapshot[6].open(fileName[6]);
	snapshot[7].open(fileName[7]);
	snapshot[8].open(fileName[8]);
	int ix, iy, iz;
	
	for( iz = 0; iz < grid.zLength; ++iz ) {
		
		for( iy = 0; iy < grid.yLength; ++iy) {
			
			//snapshot[0] << vsCPU.stressXX[INDEX( sourceLocation.x, iy, iz)] << endl;
			snapshot[0] << vsCPU.velocityX[INDEX(sourceLocation.x, iy, iz)] << endl;
			snapshot[3] << vsCPU.stressXX[INDEX(sourceLocation.x, iy, iz)] << endl;
			snapshot[6] << vsCPU.stressXY[INDEX(sourceLocation.x, iy, iz)] << endl;

		}

	}

	for( iz = 0; iz < grid.zLength; ++iz ) {
		
		for( ix = 0; ix < grid.xLength; ++ix) {
			
			snapshot[1] << vsCPU.velocityX[INDEX(ix, sourceLocation.y, iz)] << endl;
			snapshot[4] << vsCPU.stressXX[INDEX(ix, sourceLocation.y, iz)] << endl;
			snapshot[7] << vsCPU.stressXY[INDEX(ix, sourceLocation.y, iz)] << endl;

		}

	}
	

	for( iy = 0; iy < grid.yLength; ++iy ) {
		
		for( ix = 0; ix < grid.xLength; ++ix) {
			
			snapshot[2] << vsCPU.velocityX[INDEX(ix, iy, sourceLocation.z)] << endl;
			snapshot[5] << vsCPU.stressXX[INDEX(ix, iy, sourceLocation.z)] << endl;
			snapshot[8] << vsCPU.stressXY[INDEX(ix, iy, sourceLocation.z)] << endl;

		}

	}
	for(unsigned i = 0; i < 9; ++i) {
		snapshot[i].close();
	}

	delete[] snapshot;

}

void ElasticWaveEquation::snapshoot(int it)
{
	char fileName[6][100] = { 0 };
	ofstream * snapshot = new ofstream[6];

	for (size_t i = 0; i < 1; i++)
	{
		sprintf(fileName[i], "snapshot_%d_%d.txt", i, it);
		cout << fileName[i] << endl;
		snapshot[i].open(fileName[i]);
	}
	/*for (size_t i = 0; i < 6; i++)
	{
		if (!snapshot[i].is_open()) 
			cout << "Can't open the file" << endl;

	}*/
	int ix, iy, iz;
	for (iz = order / 2; iz < grid.zLength - order / 2; iz++)
	{
		for (iy = order / 2; iy < grid.yLength - order / 2; iy++)
		{
			for (ix = order / 2; ix < grid.xLength - order / 2; ix++)
			{
				snapshot[0] << stressXX[INDEX(ix, iy, iz)] << endl;
			}
		}
	}
	delete[] snapshot;

}
