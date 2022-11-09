/***************************************************************************************/
/******************************************2019-01-05***********************************/
/************************************Author:Wenqiang Wang*******************************/
/**********************Southern University of Science and Technology********************/
/***************************************************************************************/
#include "ElasticWaveEquation.h"

ElasticWaveEquation::ElasticWaveEquation(
	int timeLength, GPUDim gpuDim,
	GRID grid, Delta delta, Medium medCPU,
	SOURCELOCATION sourceLocation, DifferenceCoefficient dc, 
	float mainFrequncy )
:	timeLength( timeLength ), gpuDim( gpuDim ),
	grid( grid ), delta( delta ), medCPU( medCPU ),
	sourceLocation( sourceLocation ), dc( dc ), 
	mainFrequncy( mainFrequncy )
{
	int memsize = sizeof( float ) * grid.xLength * grid.yLength * grid.zLength;
	CUDA_MALLOC_3D_VelocityStress( vsGPU, memsize );
	CUDA_MALLOC_3D_Medium( medGPU, memsize );
	ALLOCATE_3D_VelocityStress( vsCPU, memsize );
}

ElasticWaveEquation::~ElasticWaveEquation( )
{
	CUDA_FREE_3D_VelocityStress( vsGPU );
	CUDA_FREE_3D_Medium( medGPU );
	DELETE_3D_VelocityStress( vsCPU );
}

void ElasticWaveEquation::run( )
{
	float sourceValue;
	int memsize = sizeof( float ) * grid.xLength * grid.yLength * grid.zLength;

	cudaMemcpy( medGPU.buoyancy, medCPU.buoyancy, 	memsize, cudaMemcpyHostToDevice );
	cudaMemcpy( medGPU.lambda, 	 medCPU.lambda, 	memsize, cudaMemcpyHostToDevice );
	cudaMemcpy( medGPU.mu,		 medCPU.mu, 		memsize, cudaMemcpyHostToDevice );

	for ( int it = 0; it < timeLength; it ++ )
	{
		cout << "Time Iterator: " << it << endl;
		sourceValue = sourceFunction( it );
		cout << it << " " << sourceValue << endl;

		propagate( sourceValue );
		// if( it == 800 ) {
		// 	//cudaDeviceSynchronize( );
		// 	cudaMemcpy( vsCPU.stressXX, vsGPU.stressXX, memsize, cudaMemcpyDeviceToHost );
		// 	snapshoot( it );
		// }

		if ( it % 10 == 0 )
		{
			//cudaMemcpy( vsCPU.stressXX, vsGPU.stressXX, memsize, cudaMemcpyDeviceToHost );
			cudaMemcpy( vsCPU.velocityX, vsGPU.velocityX, memsize, cudaMemcpyDeviceToHost );
			cudaMemcpy( vsCPU.stressXX, vsGPU.stressXX, memsize, cudaMemcpyDeviceToHost );
			cudaMemcpy( vsCPU.stressXY, vsGPU.stressXY, memsize, cudaMemcpyDeviceToHost );
			snapShootSlice( it );
		}
		
	}
}

float ElasticWaveEquation::sourceFunction( int timeIndex )
{
	float t0 = ceil(1.0f / (mainFrequncy * delta.t));
	float tmp = pow(PI * mainFrequncy * ((timeIndex - t0) * delta.t - 1 / mainFrequncy), 2);
	return (1 - 2 * tmp) * exp(-tmp) * 1.0e8f;
}

void ElasticWaveEquation::propagate( float sourceValue )
{
	locateSource 	<<< gpuDim.blocksPerGrid, gpuDim.threadsPerBlock >>> ( grid, delta.t, sourceValue, vsGPU, sourceLocation );	
	//cudaDeviceSynchronize();
	velocityUpdate 	<<< gpuDim.blocksPerGrid, gpuDim.threadsPerBlock >>> ( grid, delta, sourceLocation, vsGPU, dc, medGPU.buoyancy );
	//cudaDeviceSynchronize( );
	stressUpdate 	<<< gpuDim.blocksPerGrid, gpuDim.threadsPerBlock >>> ( grid, delta, vsGPU, dc, medGPU );
	//cudaDeviceSynchronize( );
}

__global__ void locateSource( GRID grid, float deltaT, float sourceValue, VelocityStress vsGPU,  SOURCELOCATION sourceLocation )
{
	int ix = threadIdx.x + blockIdx.x * blockDim.x; 
	int iy = threadIdx.y + blockIdx.y * blockDim.y;
	int iz = blockIdx.z * blockDim.z;
	if( ix == sourceLocation.x && iy == sourceLocation.y && iz == sourceLocation.z ) 
	{
		vsGPU.stressXX[INDEX( ix, iy, iz )] += pow(deltaT, 2) * sourceValue;
		vsGPU.stressYY[INDEX( ix, iy, iz )] += pow(deltaT, 2) * sourceValue;
		vsGPU.stressZZ[INDEX( ix, iy, iz )] += pow(deltaT, 2) * sourceValue;
 			//printf("sourceLocation\n");
 		//printf("ix = %d, iy = %d, iz = %d, sourceValue = %f\n", ix, iy, iz, sourceValue );
	}
}

__global__ void velocityUpdate( GRID grid, Delta delta, SOURCELOCATION sourceLocation, VelocityStress vsGPU, DifferenceCoefficient dc, float * buoyancy ) //__global__ 
{
	int i;
	int ix = threadIdx.x + blockIdx.x * blockDim.x; 
	int iy = threadIdx.y + blockIdx.y * blockDim.y;
	int iz = blockIdx.z * blockDim.z;
	
	float diffX = 0.0f, diffY = 0.0f, diffZ = 0.0f;
	//printf("velocityUpdate\n");
	// if( ix == sourceLocation.x && iy == sourceLocation.y && iz == sourceLocation.z ) 
	// {
	// 	vsGPU.stressXX[INDEX( ix, iy, iz )] += pow(delta.t, 2) * sourceValue;
 	// 		vsGPU.stressYY[INDEX( ix, iy, iz )] += pow(delta.t, 2) * sourceValue;
 	// 		vsGPU.stressZZ[INDEX( ix, iy, iz )] += pow(delta.t, 2) * sourceValue;
 	// 		//printf("sourceLocation\n");
 	// 		printf("ix = %d, iy = %d, iz = %d, sourceValue = %f\n", ix, iy, iz, sourceValue );
	// }

	if ( 
		( ( iz >= ORDER / 2 ) && ( iz < grid.zLength - ORDER / 2 ) ) && 
		( ( iy >= ORDER / 2 ) && ( iy < grid.yLength - ORDER / 2 ) ) && 
		( ( ix >= ORDER / 2 ) && ( ix < grid.xLength - ORDER / 2 ) ) )
	{
		diffX = 0.0f, diffY = 0.0f, diffZ = 0.0f;
		for ( i = 0;  i < ORDER/ 2;  i++)
		{
			diffX += dc.diff_coef[i] * (vsGPU.stressXX[INDEX(ix + 1 + i, iy, iz)] - vsGPU.stressXX[INDEX(ix - i, iy, iz)]);
			diffY += dc.diff_coef[i] * (vsGPU.stressXY[INDEX(ix, iy + i, iz)] - vsGPU.stressXY[INDEX(ix, iy - 1 - i, iz)]);
			diffZ += dc.diff_coef[i] * (vsGPU.stressXZ[INDEX(ix, iy, iz + i)] - vsGPU.stressXZ[INDEX(ix, iy, iz - 1 - i)]);
		}
		diffX /= delta.x;
		diffY /= delta.y;
		diffZ /= delta.z;
		vsGPU.velocityX[INDEX(ix, iy, iz)] += delta.t * BUOYANCYBARX(buoyancy, ix, iy, iz) * (diffX + diffY + diffZ);
		
		diffX = 0.0f, diffY = 0.0f, diffZ = 0.0f;
		for (i = 0; i < ORDER/ 2; i++)
		{
			diffX += dc.diff_coef[i] * (vsGPU.stressXY[INDEX(ix + i, iy, iz)] - vsGPU.stressXY[INDEX(ix - 1 - i, iy, iz)]);
			diffY += dc.diff_coef[i] * (vsGPU.stressYY[INDEX(ix, iy + 1 + i, iz)] - vsGPU.stressYY[INDEX(ix, iy - i, iz)]);
			diffZ += dc.diff_coef[i] * (vsGPU.stressYZ[INDEX(ix, iy, iz + i)] - vsGPU.stressYZ[INDEX(ix, iy, iz - 1 - i)]);
		}
		diffX /= delta.x;
		diffY /= delta.y;
		diffZ /= delta.z;
		vsGPU.velocityY[INDEX(ix, iy, iz)] += delta.t  * BUOYANCYBARY(buoyancy, ix, iy, iz) * (diffX + diffY + diffZ);
		
		diffX = 0.0f, diffY = 0.0f, diffZ = 0.0f;
		for (i = 0; i < ORDER / 2; i++)
		{
			diffX += dc.diff_coef[i] * (vsGPU.stressXZ[INDEX(ix + i, iy, iz)] - vsGPU.stressXZ[INDEX(ix - 1 - i, iy, iz)]);
			diffY += dc.diff_coef[i] * (vsGPU.stressYZ[INDEX(ix, iy + i, iz)] - vsGPU.stressYZ[INDEX(ix, iy - 1 - i, iz)]);
			diffZ += dc.diff_coef[i] * (vsGPU.stressZZ[INDEX(ix, iy, iz + 1 + i)] - vsGPU.stressZZ[INDEX(ix, iy, iz - i)]);
		}
		diffX /= delta.x;
		diffY /= delta.y;
		diffZ /= delta.z;
		vsGPU.velocityZ[INDEX(ix, iy, iz)] += delta.t  * BUOYANCYBARZ(buoyancy, ix, iy, iz) * (diffX + diffY + diffZ);

	}
}

__global__ void stressUpdate( GRID grid, Delta delta, VelocityStress vsGPU, DifferenceCoefficient dc, Medium medGPU )//__global__ 
{

	int i;
	int ix = threadIdx.x + blockIdx.x * blockDim.x; 
	int iy = threadIdx.y + blockIdx.y * blockDim.y;
	int iz = blockIdx.z * blockDim.z;

	float diffX = 0.0f, diffY = 0.0f, diffZ = 0.0f;
	float muBarXY = 0.0f, muBarXZ = 0.0f, muBarYZ = 0.0f;

	if ( 
		( ( iz >= ORDER / 2 ) && ( iz < grid.zLength - ORDER / 2 ) ) && 
		( ( iy >= ORDER / 2 ) && ( iy < grid.yLength - ORDER / 2 ) ) && 
		( ( ix >= ORDER / 2 ) && ( ix < grid.xLength - ORDER / 2 ) ) )
	{
		diffX = 0.0f, diffY = 0.0f, diffZ = 0.0f;
		for (i = 0; i < ORDER / 2; i++)
		{
			diffX += dc.diff_coef[i] * (vsGPU.velocityX[INDEX(ix + i, iy, iz)] - vsGPU.velocityX[INDEX(ix - 1 - i, iy, iz)]);
			diffY += dc.diff_coef[i] * (vsGPU.velocityY[INDEX(ix, iy + i, iz)] - vsGPU.velocityY[INDEX(ix, iy - 1 - i, iz)]);
			diffZ += dc.diff_coef[i] * (vsGPU.velocityZ[INDEX(ix, iy, iz + i)] - vsGPU.velocityZ[INDEX(ix, iy, iz - 1 - i)]);
		}
		diffX /= delta.x;
		diffY /= delta.y;
		diffZ /= delta.z;
		vsGPU.stressXX[INDEX(ix, iy, iz)] += delta.t * 
			(( medGPU.lambda[INDEX(ix, iy, iz)] + 2.0f * medGPU.mu[INDEX(ix, iy, iz)]) * diffX +
			   medGPU.lambda[INDEX(ix, iy, iz)] * (diffY + diffZ) );				 
		vsGPU.stressYY[INDEX(ix, iy, iz)] += delta.t * 								 
			(( medGPU.lambda[INDEX(ix, iy, iz)] + 2.0f * medGPU.mu[INDEX(ix, iy, iz)]) * diffY +
			   medGPU.lambda[INDEX(ix, iy, iz)] * (diffX + diffZ) );					 
		vsGPU.stressZZ[INDEX(ix, iy, iz)] += delta.t * 								 
			(( medGPU.lambda[INDEX(ix, iy, iz)] + 2.0f * medGPU.mu[INDEX(ix, iy, iz)]) * diffZ +
			   medGPU.lambda[INDEX(ix, iy, iz)] * (diffX + diffY) );
	


		muBarXY = MUBARXY(medGPU.mu, ix, iy, iz);
		muBarXZ = MUBARXZ(medGPU.mu, ix, iy, iz);
		muBarYZ = MUBARYZ(medGPU.mu, ix, iy, iz);

		diffX = 0.0f, diffY = 0.0f;
		for (i = 0; i < ORDER/ 2; i++)
		{
			diffY += dc.diff_coef[i] * (vsGPU.velocityX[INDEX(ix, iy + 1 + i, iz)] - vsGPU.velocityX[INDEX(ix, iy - i, iz)]);
			diffX += dc.diff_coef[i] * (vsGPU.velocityY[INDEX(ix + 1 + i, iy, iz)] - vsGPU.velocityY[INDEX(ix - i, iy, iz)]);
		}
		diffX /= delta.x;
		diffY /= delta.y;
		vsGPU.stressXY[INDEX(ix, iy, iz)] += delta.t * muBarXY * (diffY + diffX);

		diffX = 0.0f, diffZ = 0.0f;
		for (i = 0; i < ORDER / 2; i++)
		{
			diffZ += dc.diff_coef[i] * (vsGPU.velocityX[INDEX(ix, iy, iz + 1 + i)] - vsGPU.velocityX[INDEX(ix, iy, iz - i)]);
			diffX += dc.diff_coef[i] * (vsGPU.velocityZ[INDEX(ix + 1 + i, iy, iz)] - vsGPU.velocityZ[INDEX(ix - i, iy, iz)]);
		}
		diffX /= delta.x;
		diffZ /= delta.z;
		vsGPU.stressXZ[INDEX(ix, iy, iz)] += delta.t * muBarXZ * (diffZ + diffX);

		diffY = 0.0f, diffZ = 0.0f;
		for (i = 0; i < ORDER / 2; i++)
		{
			diffZ += dc.diff_coef[i] * (vsGPU.velocityY[INDEX(ix, iy, iz + 1 + i)] - vsGPU.velocityY[INDEX(ix, iy, iz - i)]);
			diffY += dc.diff_coef[i] * (vsGPU.velocityZ[INDEX(ix, iy + 1 + i, iz)] - vsGPU.velocityZ[INDEX(ix, iy - i, iz)]);
		}
		diffY /= delta.y;
		diffZ /= delta.z;
		vsGPU.stressYZ[INDEX(ix, iy, iz)] += delta.t * muBarYZ * (diffZ + diffY);

	}
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

void ElasticWaveEquation::snapShootVolum(int it)
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
	for (iz = ORDER / 2; iz < grid.zLength - ORDER / 2; iz++)
	{
		for (iy = ORDER / 2; iy < grid.yLength - ORDER / 2; iy++)
		{
			for (ix = ORDER / 2; ix < grid.xLength - ORDER / 2; ix++)
			{
				snapshot[0] << vsCPU.stressXX[INDEX(ix, iy, iz)] << endl;
			}
		}
	}
	delete[] snapshot;

}

