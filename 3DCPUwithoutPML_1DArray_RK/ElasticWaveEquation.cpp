/***************************************************************************************/
/******************************************2019-01-24***********************************/
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
	float diffx = 0.0f, diffy = 0.0f, diffz = 0.0f;
	float muBarXY = 0.0f, muBarXZ = 0.0f, muBarYZ = 0.0f;
	float h1 = 0.0f, h2 = 0.0f, h3 = 0.0f, h4 = 0.0f;
	// float hx1 = 0.0f, hx2 = 0.0f,  hx3 = 0.0f,  hx4 = 0.0f;
	// float hy1 = 0.0f, hy2 = 0.0f,  hy3 = 0.0f,  hy4 = 0.0f;
	// float hz1 = 0.0f, hz2 = 0.0f,  hz3 = 0.0f,  hz4 = 0.0f;

	for (iz = order / 2; iz < grid.zLength - order / 2; iz++)
	{
		for (iy = order / 2; iy < grid.yLength - order / 2; iy++)
		{
			for (ix = order / 2; ix < grid.xLength - order / 2; ix++)
			{
				diffx = 0.0f, diffy = 0.0f, diffz = 0.0f;
				RK( stressXX, x );		RK( stressXY, y );		RK( stressXZ, z );//macro gives the solution diffx diffy diffz
				velocityX[INDEX(ix, iy, iz)] += BUOYANCYBARX(buoyancy, ix, iy, iz) * (diffx + diffy + diffz);

				diffx = 0.0f, diffy = 0.0f, diffz = 0.0f;
				RK( stressXY, x );		RK( stressYY, y );		RK( stressYZ, z );
				velocityY[INDEX(ix, iy, iz)] += BUOYANCYBARY(buoyancy, ix, iy, iz) * (diffx + diffy + diffz);
				
				diffx = 0.0f, diffy = 0.0f, diffz = 0.0f;
				RK( stressXZ, x );		RK( stressYZ, y );		RK( stressZZ, z );
				velocityZ[INDEX(ix, iy, iz)] += BUOYANCYBARZ(buoyancy, ix, iy, iz) * (diffx + diffy + diffz);

			}
		}
	}

	for (iz = order / 2; iz < grid.zLength - order / 2; iz++)
	{
		for (iy = order / 2; iy < grid.yLength - order / 2; iy++)
		{
			for (ix = order / 2; ix < grid.xLength - order / 2; ix++)
			{
				diffx = 0.0f, diffy = 0.0f, diffz = 0.0f;
				RK( velocityX, x );	RK( velocityY, y );	RK( velocityZ, z );
				stressXX[INDEX(ix, iy, iz)] += delta.t * 
					(( lambda[INDEX(ix, iy, iz)] + 2.0f * mu[INDEX(ix, iy, iz)]) * diffx +
					   lambda[INDEX(ix, iy, iz)] * (diffy + diffz) );				 
				stressYY[INDEX(ix, iy, iz)] += delta.t * 								 
					(( lambda[INDEX(ix, iy, iz)] + 2.0f * mu[INDEX(ix, iy, iz)]) * diffy +
					   lambda[INDEX(ix, iy, iz)] * (diffx + diffz) );					 
				stressZZ[INDEX(ix, iy, iz)] += delta.t * 								 
					(( lambda[INDEX(ix, iy, iz)] + 2.0f * mu[INDEX(ix, iy, iz)]) * diffz +
					   lambda[INDEX(ix, iy, iz)] * (diffx + diffy) );
				
				muBarXY = MUBARXY(mu, ix, iy, iz);
				muBarXZ = MUBARXZ(mu, ix, iy, iz);
				muBarYZ = MUBARYZ(mu, ix, iy, iz);

				diffx = 0.0f, diffy = 0.0f;
				RK( velocityX, y );	RK( velocityY, x );
				stressXY[INDEX(ix, iy, iz)] += delta.t * muBarXY * (diffy + diffx);

				diffx = 0.0f, diffz = 0.0f;
				RK( velocityX, z );	RK( velocityZ, x );
				stressXZ[INDEX(ix, iy, iz)] += delta.t * muBarXZ * (diffz + diffx);

				diffy = 0.0f, diffz = 0.0f;
				RK( velocityY, z );	RK( velocityZ, y );
				stressYZ[INDEX(ix, iy, iz)] += delta.t * muBarYZ * (diffz + diffy);


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
