/***************************************************************************************/
/******************************************2019-01-05***********************************/
/************************************Author:Wenqiang Wang*******************************/
/**********************Southern University of Science and Technology********************/
/***************************************************************************************/
#include "ElasticWaveEquation.h"

ElasticWaveEquation::ElasticWaveEquation(
		int timeLength, GPUDim gpuDim, int pml,
		GRID grid, Delta delta,Medium medCPU,
		SOURCELOCATION sourceLocationCPU, STATIONLOCATION stationLocationCPU, 
		DifferenceCoefficient dc, 
		float mainFrequncy )
:	timeLength( timeLength ), gpuDim( gpuDim ), pml( pml ),
	grid( grid ), delta( delta ), medCPU( medCPU ),
	sourceLocationCPU( sourceLocationCPU ), stationLocationCPU( stationLocationCPU ),
	dc( dc ), 
	mainFrequncy( mainFrequncy )
{
	int memsize = sizeof( float ) * grid.xLength * grid.yLength * grid.zLength;

	stationDataCPU = new float[timeLength];
	cudaMalloc( &stationDataGPU, timeLength * sizeof( float ) );
	cudaMemset( stationDataGPU, 0, timeLength * sizeof( float ) );

	calStartEndPoint( );

	CUDA_MALLOC_3D_VelocityStress( vsGPU, memsize );
	CUDA_MALLOC_3D_Medium( medGPU, memsize );


	ALLOCATE_3D_VelocityStress( vsCPU, memsize );

	for( int i = 0; i < nAREA; ++ i ) 
	{
		vsPMLCPU[i].xLength = vsPMLCPU[i].endX - vsPMLCPU[i].startX + 1;
		vsPMLCPU[i].yLength = vsPMLCPU[i].endY - vsPMLCPU[i].startY + 1;
		vsPMLCPU[i].zLength = vsPMLCPU[i].endZ - vsPMLCPU[i].startZ + 1;
		memsize = sizeof( float ) * vsPMLCPU[i].xLength * vsPMLCPU[i].yLength * vsPMLCPU[i].zLength;
		CUDA_MALLOC_3D_VelocityStressPML( vsPMLCPU[i], memsize );
	}
	vsPMLCPUTovsPML();
	sourceLocationCPUTosourceLocationGPU();
	stationLocationCPUTostationLocationGPU();
}

void ElasticWaveEquation::calStartEndPoint( )
{
	vsPMLCPU[VERTEX0].startX  	= START; 			 	vsPMLCPU[VERTEX0].startY  	= START; 			  	vsPMLCPU[VERTEX0].startZ  	= START; 			   
	vsPMLCPU[VERTEX0].endX  	= CalculateAreaStart; 	vsPMLCPU[VERTEX0].endY  	= CalculateAreaStart; 	vsPMLCPU[VERTEX0].endZ  	= CalculateAreaStart;
	vsPMLCPU[VERTEX1].startX  	= START; 				vsPMLCPU[VERTEX1].startY  	= START; 				vsPMLCPU[VERTEX1].startZ  	= CalculateAreaEndZ; 
	vsPMLCPU[VERTEX1].endX  	= CalculateAreaStart; 	vsPMLCPU[VERTEX1].endY  	= CalculateAreaStart; 	vsPMLCPU[VERTEX1].endZ  	= EndZ; 
	vsPMLCPU[VERTEX2].startX  	= START; 				vsPMLCPU[VERTEX2].startY  	= CalculateAreaEndY; 	vsPMLCPU[VERTEX2].startZ  	= START; 
	vsPMLCPU[VERTEX2].endX  	= CalculateAreaStart; 	vsPMLCPU[VERTEX2].endY  	= EndY; 				vsPMLCPU[VERTEX2].endZ  	= CalculateAreaStart; 
	vsPMLCPU[VERTEX3].startX  	= START; 				vsPMLCPU[VERTEX3].startY  	= CalculateAreaEndY; 	vsPMLCPU[VERTEX3].startZ  	= CalculateAreaEndZ; 
	vsPMLCPU[VERTEX3].endX  	= CalculateAreaStart; 	vsPMLCPU[VERTEX3].endY  	= EndY; 				vsPMLCPU[VERTEX3].endZ  	= EndZ; 

	vsPMLCPU[VERTEX4].startX  	= CalculateAreaEndX; 	vsPMLCPU[VERTEX4].startY  	= START; 				vsPMLCPU[VERTEX4].startZ  	= START; 
	vsPMLCPU[VERTEX4].endX  	= EndX; 				vsPMLCPU[VERTEX4].endY  	= CalculateAreaStart; 	vsPMLCPU[VERTEX4].endZ  	= CalculateAreaStart;
	vsPMLCPU[VERTEX5].startX  	= CalculateAreaEndX; 	vsPMLCPU[VERTEX5].startY  	= START; 				vsPMLCPU[VERTEX5].startZ  	= CalculateAreaEndZ; 
	vsPMLCPU[VERTEX5].endX  	= EndX; 				vsPMLCPU[VERTEX5].endY  	= CalculateAreaStart; 	vsPMLCPU[VERTEX5].endZ  	= EndZ;
	vsPMLCPU[VERTEX6].startX  	= CalculateAreaEndX; 	vsPMLCPU[VERTEX6].startY  	= CalculateAreaEndY; 	vsPMLCPU[VERTEX6].startZ	= START; 
	vsPMLCPU[VERTEX6].endX  	= EndX; 				vsPMLCPU[VERTEX6].endY  	= EndY; 				vsPMLCPU[VERTEX6].endZ  	= CalculateAreaStart; 
	vsPMLCPU[VERTEX7].startX  	= CalculateAreaEndX;	vsPMLCPU[VERTEX7].startY  	= CalculateAreaEndY; 	vsPMLCPU[VERTEX7].startZ  	= CalculateAreaEndZ; 
	vsPMLCPU[VERTEX7].endX  	= EndX; 				vsPMLCPU[VERTEX7].endY  	= EndY; 				vsPMLCPU[VERTEX7].endZ  	= EndZ; 


	vsPMLCPU[SURFACE0].startX 	= START; 				vsPMLCPU[SURFACE0].startY 	= CalculateAreaStart; 	vsPMLCPU[SURFACE0].startZ 	= CalculateAreaStart; 
	vsPMLCPU[SURFACE0].endX 	= CalculateAreaStart; 	vsPMLCPU[SURFACE0].endY 	= CalculateAreaEndY; 	vsPMLCPU[SURFACE0].endZ 	= CalculateAreaEndZ; 
	vsPMLCPU[SURFACE1].startX 	= CalculateAreaEndX; 	vsPMLCPU[SURFACE1].startY 	= CalculateAreaStart; 	vsPMLCPU[SURFACE1].startZ 	= CalculateAreaStart; 
	vsPMLCPU[SURFACE1].endX 	= EndX; 				vsPMLCPU[SURFACE1].endY 	= CalculateAreaEndY; 	vsPMLCPU[SURFACE1].endZ 	= CalculateAreaEndZ;
	vsPMLCPU[SURFACE2].startX 	= CalculateAreaStart; 	vsPMLCPU[SURFACE2].startY 	= START; 				vsPMLCPU[SURFACE2].startZ 	= CalculateAreaStart; 
	vsPMLCPU[SURFACE2].endX 	= CalculateAreaEndX; 	vsPMLCPU[SURFACE2].endY 	= CalculateAreaStart; 	vsPMLCPU[SURFACE2].endZ 	= CalculateAreaEndZ; 
	vsPMLCPU[SURFACE3].startX 	= CalculateAreaStart; 	vsPMLCPU[SURFACE3].startY 	= CalculateAreaEndY; 	vsPMLCPU[SURFACE3].startZ 	= CalculateAreaStart; 
	vsPMLCPU[SURFACE3].endX 	= CalculateAreaEndX; 	vsPMLCPU[SURFACE3].endY 	= EndY; 				vsPMLCPU[SURFACE3].endZ 	= CalculateAreaEndZ;
	vsPMLCPU[SURFACE4].startX 	= CalculateAreaStart; 	vsPMLCPU[SURFACE4].startY 	= CalculateAreaStart; 	vsPMLCPU[SURFACE4].startZ 	= START; 
	vsPMLCPU[SURFACE4].endX 	= CalculateAreaEndX; 	vsPMLCPU[SURFACE4].endY 	= CalculateAreaEndY; 	vsPMLCPU[SURFACE4].endZ 	= CalculateAreaStart;
	vsPMLCPU[SURFACE5].startX 	= CalculateAreaStart; 	vsPMLCPU[SURFACE5].startY 	= CalculateAreaStart; 	vsPMLCPU[SURFACE5].startZ 	= CalculateAreaEndZ; 
	vsPMLCPU[SURFACE5].endX 	= CalculateAreaEndX; 	vsPMLCPU[SURFACE5].endY 	= CalculateAreaEndY; 	vsPMLCPU[SURFACE5].endZ 	= EndZ; 


	vsPMLCPU[EDGE0].startX    	= START; 				vsPMLCPU[EDGE0].startY    	= START; 				vsPMLCPU[EDGE0].startZ    	= CalculateAreaStart; 
	vsPMLCPU[EDGE0].endX    	= CalculateAreaStart; 	vsPMLCPU[EDGE0].endY    	= CalculateAreaStart; 	vsPMLCPU[EDGE0].endZ    	= CalculateAreaEndZ; 
	vsPMLCPU[EDGE1].startX    	= START; 				vsPMLCPU[EDGE1].startY    	= CalculateAreaEndY; 	vsPMLCPU[EDGE1].startZ    	= CalculateAreaStart; 
	vsPMLCPU[EDGE1].endX    	= CalculateAreaStart; 	vsPMLCPU[EDGE1].endY    	= EndY; 				vsPMLCPU[EDGE1].endZ    	= CalculateAreaEndZ;
	vsPMLCPU[EDGE2].startX    	= START; 				vsPMLCPU[EDGE2].startY    	= CalculateAreaStart; 	vsPMLCPU[EDGE2].startZ   	= START; 
	vsPMLCPU[EDGE2].endX    	= CalculateAreaStart; 	vsPMLCPU[EDGE2].endY    	= CalculateAreaEndY; 	vsPMLCPU[EDGE2].endZ    	= CalculateAreaStart;
	vsPMLCPU[EDGE3].startX    	= START; 				vsPMLCPU[EDGE3].startY    	= CalculateAreaStart; 	vsPMLCPU[EDGE3].startZ    	= CalculateAreaEndZ; 
	vsPMLCPU[EDGE3].endX    	= CalculateAreaStart; 	vsPMLCPU[EDGE3].endY    	= CalculateAreaEndY; 	vsPMLCPU[EDGE3].endZ    	= EndZ;

	vsPMLCPU[EDGE4].startX    	= CalculateAreaEndX; 	vsPMLCPU[EDGE4].startY    	= START; 				vsPMLCPU[EDGE4].startZ    	= CalculateAreaStart; 
	vsPMLCPU[EDGE4].endX    	= EndX; 				vsPMLCPU[EDGE4].endY    	= CalculateAreaStart; 	vsPMLCPU[EDGE4].endZ    	= CalculateAreaEndZ;
	vsPMLCPU[EDGE5].startX    	= CalculateAreaEndX; 	vsPMLCPU[EDGE5].startY    	= CalculateAreaEndY; 	vsPMLCPU[EDGE5].startZ    	= CalculateAreaStart; 
	vsPMLCPU[EDGE5].endX    	= EndX;					vsPMLCPU[EDGE5].endY    	= EndY; 				vsPMLCPU[EDGE5].endZ   	= CalculateAreaEndZ;
	vsPMLCPU[EDGE6].startX    	= CalculateAreaEndX; 	vsPMLCPU[EDGE6].startY    	= CalculateAreaStart; 	vsPMLCPU[EDGE6].startZ    	= START; 
	vsPMLCPU[EDGE6].endX    	= EndX; 				vsPMLCPU[EDGE6].endY    	= CalculateAreaEndY; 	vsPMLCPU[EDGE6].endZ    	= CalculateAreaStart;
	vsPMLCPU[EDGE7].startX    	= CalculateAreaEndX; 	vsPMLCPU[EDGE7].startY    	= CalculateAreaStart; 	vsPMLCPU[EDGE7].startZ    	= CalculateAreaEndZ; 
	vsPMLCPU[EDGE7].endX    	= EndX; 				vsPMLCPU[EDGE7].endY    	= CalculateAreaEndY; 	vsPMLCPU[EDGE7].endZ    	= EndZ;

	vsPMLCPU[EDGE8].startX    	= CalculateAreaStart; 	vsPMLCPU[EDGE8].startY    	= START; 				vsPMLCPU[EDGE8].startZ    	= START; 
	vsPMLCPU[EDGE8].endX    	= CalculateAreaEndX; 	vsPMLCPU[EDGE8].endY    	= CalculateAreaStart; 	vsPMLCPU[EDGE8].endZ    	= CalculateAreaStart; 
	vsPMLCPU[EDGE9].startX    	= CalculateAreaStart; 	vsPMLCPU[EDGE9].startY    	= START; 				vsPMLCPU[EDGE9].startZ    	= CalculateAreaEndZ; 
	vsPMLCPU[EDGE9].endX    	= CalculateAreaEndX; 	vsPMLCPU[EDGE9].endY    	= CalculateAreaStart; 	vsPMLCPU[EDGE9].endZ    	= EndZ;
	vsPMLCPU[EDGE10].startX   	= CalculateAreaStart; 	vsPMLCPU[EDGE10].startY   	= CalculateAreaEndY; 	vsPMLCPU[EDGE10].startZ   	= START; 
	vsPMLCPU[EDGE10].endX   	= CalculateAreaEndX; 	vsPMLCPU[EDGE10].endY   	= EndY; 				vsPMLCPU[EDGE10].endZ   	= CalculateAreaStart; 
	vsPMLCPU[EDGE11].startX   	= CalculateAreaStart; 	vsPMLCPU[EDGE11].startY   	= CalculateAreaEndY; 	vsPMLCPU[EDGE11].startZ   	= CalculateAreaEndZ; 
	vsPMLCPU[EDGE11].endX   	= CalculateAreaEndX; 	vsPMLCPU[EDGE11].endY   	= EndY; 				vsPMLCPU[EDGE11].endZ   	= EndZ; 
}

void ElasticWaveEquation::vsPMLCPUTovsPML( )
{
	cudaMemcpyToSymbol( vsPML, vsPMLCPU, sizeof( VelocityStressPML ) * nAREA );
}

void ElasticWaveEquation::sourceLocationCPUTosourceLocationGPU()
{
	cudaMemcpyToSymbol( sourceLocationGPU, &sourceLocationCPU, sizeof( STATIONLOCATION ) );
}

void ElasticWaveEquation::stationLocationCPUTostationLocationGPU()
{
	cudaMemcpyToSymbol( stationLocationGPU, &stationLocationCPU, sizeof( STATIONLOCATION ) );
}

ElasticWaveEquation::~ElasticWaveEquation( )
{
	CUDA_FREE_3D_VelocityStress( vsGPU );
	CUDA_FREE_3D_Medium( medGPU );
	DELETE_3D_VelocityStress( vsCPU );

	delete[] stationDataCPU;
	cudaFree( stationDataGPU );

	for( int i = 0; i < nAREA; ++ i ) 
	{
		 DELETE_MALLOC_3D_VelocityStressPML( vsPMLCPU[i] );
	}
}

void ElasticWaveEquation::run( )
{
	float sourceValue;
	int memsize = sizeof( float ) * grid.xLength * grid.yLength * grid.zLength;
	getLastCudaError( "CHECK" );
	CHECK( cudaMemcpy( medGPU.buoyancy,  medCPU.buoyancy, 	memsize, cudaMemcpyHostToDevice ) );
	CHECK( cudaMemcpy( medGPU.lambda, 	 medCPU.lambda, 	memsize, cudaMemcpyHostToDevice ) );
	CHECK( cudaMemcpy( medGPU.mu,		 medCPU.mu, 		memsize, cudaMemcpyHostToDevice ) );

	// ofstream curveLine1("curveLine1.txt");
	// ofstream curveLine2("curveLine2.txt");
	// ofstream curveLine3("curveLine3.txt");

	for ( int it = 0; it < timeLength; it ++ )
	{
		cout << "Time Iterator: " << it << endl;
		sourceValue = sourceFunction( it );
		//cout << it << " " << sourceValue << endl;

		propagate( it, sourceValue );
		if ( it % 100 == 0 && it >= 500 )
		{
			//cudaMemcpy( vsCPU.stressXX, vsGPU.stressXX, memsize, cudaMemcpyDeviceToHost );
			cudaMemcpy( vsCPU.velocityX, vsGPU.velocityX, memsize, cudaMemcpyDeviceToHost );
			//			cudaMemcpy( vsCPU.velocityX, vsGPU.velocityX, memsize, cudaMemcpyDeviceToHost );
			cudaMemcpy( vsCPU.stressXX, vsGPU.stressXX, memsize, cudaMemcpyDeviceToHost );
			cudaMemcpy( vsCPU.stressXY, vsGPU.stressXY, memsize, cudaMemcpyDeviceToHost );
			snapShootSlice( it );
		}


		// curveLine1 << vsPMLCPU.stressXX[INDEX( pml + 10, sourceLocationCPU.y, sourceLocation.z )] << endl;
		// curveLine2 << stressXX[sourceLocation.x][pml + 10][sourceLocation.z] << endl;
		// curveLine3 << stressXX[sourceLocation.x][sourceLocation.y][pml + 10] << endl;
		
	}

	cudaMemcpy( stationDataCPU, stationDataGPU, timeLength * sizeof( float), cudaMemcpyDeviceToHost );
	stationDataOutput();

	//curveLine1.close();
	//curveLine2.close();
	//curveLine3.close();


}

float ElasticWaveEquation::sourceFunction( int timeIndex )
{
	float t0 = ceil(1.0f / (mainFrequncy * delta.t));
	float tmp = pow(PI * mainFrequncy * ((timeIndex - t0) * delta.t - 1 / mainFrequncy), 2);
	return (1 - 2 * tmp) * exp(-tmp);
}

void ElasticWaveEquation::propagate( int timeIndex, float sourceValue )
{
	//cout << "bx = " << gpuDim.blocksPerGrid.x 	<< "by = " << gpuDim.blocksPerGrid.y 	<< "bz = " << gpuDim.blocksPerGrid.z 	<< endl;
	//cout << "tx = " << gpuDim.threadsPerBlock.x << "ty = " << gpuDim.threadsPerBlock.y 	<< "tz = " << gpuDim.threadsPerBlock.z 	<< endl;
	
	getLastCudaError( "CHECK" );
	locateSource 	<<< 1, 32 >>> ( grid, delta.t, sourceValue, vsGPU );	
	velocityUpdate 	<<< gpuDim.blocksPerGrid, gpuDim.threadsPerBlock >>> ( grid, delta, vsGPU, dc, medGPU, pml );
	stressUpdate 	<<< gpuDim.blocksPerGrid, gpuDim.threadsPerBlock >>> ( grid, delta, vsGPU, dc, medGPU, pml );
	recordStationData <<< 1, 32 >>> ( grid, timeIndex, stationDataGPU, vsGPU );
	cudaDeviceSynchronize( );
	
}


__global__ void locateSource( GRID grid, float deltaT, float sourceValue, VelocityStress vsGPU )
{
	int ix = threadIdx.x + blockIdx.x * blockDim.x; 
	int iy = threadIdx.y + blockIdx.y * blockDim.y;
	int iz = blockIdx.z * blockDim.z;
	if( ix == 0 && iy == 0 && iz == 0 ) 
	{
		vsGPU.stressXX[INDEX( sourceLocationGPU.x, sourceLocationGPU.y, sourceLocationGPU.z )] += deltaT * deltaT * sourceValue;//pow(deltaT, 2)
		vsGPU.stressYY[INDEX( sourceLocationGPU.x, sourceLocationGPU.y, sourceLocationGPU.z )] += deltaT * deltaT * sourceValue;//pow(deltaT, 2)
		vsGPU.stressZZ[INDEX( sourceLocationGPU.x, sourceLocationGPU.y, sourceLocationGPU.z )] += deltaT * deltaT * sourceValue;//pow(deltaT, 2)
 		//printf("vsPML[SURFACE0].velocityXx[0] = %f\n", vsPML[SURFACE0].velocityXx[0]);
 		//printf("vsPML[EDGE0].velocityXx[0] = %f\n", vsPML[SURFACE0].velocityXx[0]);
 		//	printf("sourceLocationGPU\n");
 		//	printf("ix = %d, iy = %d, iz = %d, sourceValue = %f\n", ix, iy, iz, sourceValue );
	}
}

__global__ void recordStationData( GRID grid, int timeIndex, float * stationDataGPU, VelocityStress vsGPU )
{
	int ix = threadIdx.x + blockIdx.x * blockDim.x; 
	int iy = threadIdx.y + blockIdx.y * blockDim.y;
	int iz = blockIdx.z * blockDim.z;
	if( ix == 0 && iy == 0 && iz == 0 ) 
	{
		stationDataGPU[timeIndex] = vsGPU.stressXX[INDEX( stationLocationGPU.x, stationLocationGPU.y, stationLocationGPU.z )];//pow(deltaT, 2)
	}
}

__global__ void velocityUpdate( GRID grid, Delta delta, VelocityStress vsGPU, DifferenceCoefficient dc, Medium medGPU, int pml ) //__global__ 
{
	int i;
	int ix = threadIdx.x + blockIdx.x * blockDim.x; 
	int iy = threadIdx.y + blockIdx.y * blockDim.y;
	int iz = blockIdx.z * blockDim.z;

	if ( AreaCondition( ix, iy, iz, CalculateAreaStart, CalculateAreaStart, CalculateAreaStart, CalculateAreaEndX, CalculateAreaEndY, CalculateAreaEndZ ) )
	{
		float diffX = 0.0f, diffY = 0.0f, diffZ = 0.0f;
		for ( i = 0;  i < order;  i++)
		{
			diffX += dc.diff_coef[i] * (vsGPU.stressXX[INDEX(ix + 1 + i, iy, iz)] - vsGPU.stressXX[INDEX(ix - i, iy, iz)]);
			diffY += dc.diff_coef[i] * (vsGPU.stressXY[INDEX(ix, iy + i, iz)] - vsGPU.stressXY[INDEX(ix, iy - 1 - i, iz)]);
			diffZ += dc.diff_coef[i] * (vsGPU.stressXZ[INDEX(ix, iy, iz + i)] - vsGPU.stressXZ[INDEX(ix, iy, iz - 1 - i)]);
		}
		diffX /= delta.x;
		diffY /= delta.y;
		diffZ /= delta.z;
		vsGPU.velocityX[INDEX(ix, iy, iz)] += delta.t * BUOYANCYBARX(medGPU.buoyancy, ix, iy, iz) * (diffX + diffY + diffZ);
		
		diffX = 0.0f, diffY = 0.0f, diffZ = 0.0f;
		for (i = 0; i < order; i++)
		{
			diffX += dc.diff_coef[i] * (vsGPU.stressXY[INDEX(ix + i, iy, iz)] - vsGPU.stressXY[INDEX(ix - 1 - i, iy, iz)]);
			diffY += dc.diff_coef[i] * (vsGPU.stressYY[INDEX(ix, iy + 1 + i, iz)] - vsGPU.stressYY[INDEX(ix, iy - i, iz)]);
			diffZ += dc.diff_coef[i] * (vsGPU.stressYZ[INDEX(ix, iy, iz + i)] - vsGPU.stressYZ[INDEX(ix, iy, iz - 1 - i)]);
		}
		diffX /= delta.x;
		diffY /= delta.y;
		diffZ /= delta.z;
		vsGPU.velocityY[INDEX(ix, iy, iz)] += delta.t  * BUOYANCYBARY(medGPU.buoyancy, ix, iy, iz) * (diffX + diffY + diffZ);
		
		diffX = 0.0f, diffY = 0.0f, diffZ = 0.0f;
		for (i = 0; i < order; i++)
		{
			diffX += dc.diff_coef[i] * (vsGPU.stressXZ[INDEX(ix + i, iy, iz)] - vsGPU.stressXZ[INDEX(ix - 1 - i, iy, iz)]);
			diffY += dc.diff_coef[i] * (vsGPU.stressYZ[INDEX(ix, iy + i, iz)] - vsGPU.stressYZ[INDEX(ix, iy - 1 - i, iz)]);
			diffZ += dc.diff_coef[i] * (vsGPU.stressZZ[INDEX(ix, iy, iz + 1 + i)] - vsGPU.stressZZ[INDEX(ix, iy, iz - i)]);
		}
		diffX /= delta.x;
		diffY /= delta.y;
		diffZ /= delta.z;
		vsGPU.velocityZ[INDEX(ix, iy, iz)] += delta.t  * BUOYANCYBARZ(medGPU.buoyancy, ix, iy, iz) * (diffX + diffY + diffZ);

	}
else if( AreaCondition( ix, iy, iz, vsPML[VERTEX0].startX, vsPML[VERTEX0].startY, vsPML[VERTEX0].startZ, vsPML[VERTEX0].endX, vsPML[VERTEX0].endY, vsPML[VERTEX0].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, vsPML[VERTEX0].endX - ix, 			pml, delta.x );
		float dampingY = DAMPING( Vs_wave, vsPML[VERTEX0].endY - iy, 			pml, delta.y );
		float dampingZ = DAMPING( Vs_wave, vsPML[VERTEX0].endZ - iz, 			pml, delta.z );

		VelocityPMLUpdate(  VERTEX0 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[VERTEX1].startX, vsPML[VERTEX1].startY, vsPML[VERTEX1].startZ, vsPML[VERTEX1].endX, vsPML[VERTEX1].endY, vsPML[VERTEX1].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, vsPML[VERTEX1].endX - ix, 			pml, delta.x );
		float dampingY = DAMPING( Vs_wave, vsPML[VERTEX1].endY - iy, 			pml, delta.y );
		float dampingZ = DAMPING( Vs_wave, iz - vsPML[VERTEX1].startZ + 1, 		pml, delta.z );

		VelocityPMLUpdate( VERTEX1 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[VERTEX2].startX, vsPML[VERTEX2].startY, vsPML[VERTEX2].startZ, vsPML[VERTEX2].endX, vsPML[VERTEX2].endY, vsPML[VERTEX2].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, vsPML[VERTEX2].endX - ix, 			pml, delta.x );
		float dampingY = DAMPING( Vs_wave, iy - vsPML[VERTEX2].startY + 1, 		pml, delta.y );
		float dampingZ = DAMPING( Vs_wave, vsPML[VERTEX2].endZ - iz, 			pml, delta.z );

		VelocityPMLUpdate( VERTEX2 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[VERTEX3].startX, vsPML[VERTEX3].startY, vsPML[VERTEX3].startZ, vsPML[VERTEX3].endX, vsPML[VERTEX3].endY, vsPML[VERTEX3].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, vsPML[VERTEX3].endX - ix, 			pml, delta.x );
		float dampingY = DAMPING( Vs_wave, iy - vsPML[VERTEX3].startY + 1, 		pml, delta.y );
		float dampingZ = DAMPING( Vs_wave, iz - vsPML[VERTEX3].startZ + 1, 		pml, delta.z );

		VelocityPMLUpdate( VERTEX3 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[VERTEX4].startX, vsPML[VERTEX4].startY, vsPML[VERTEX4].startZ, vsPML[VERTEX4].endX, vsPML[VERTEX4].endY, vsPML[VERTEX4].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, ix - vsPML[VERTEX4].startX + 1, 		pml, delta.x );
		float dampingY = DAMPING( Vs_wave, vsPML[VERTEX4].endY - iy, 			pml, delta.y );
		float dampingZ = DAMPING( Vs_wave, vsPML[VERTEX4].endZ - iz, 			pml, delta.z );

		VelocityPMLUpdate( VERTEX4 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[VERTEX5].startX, vsPML[VERTEX5].startY, vsPML[VERTEX5].startZ, vsPML[VERTEX5].endX, vsPML[VERTEX5].endY, vsPML[VERTEX5].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, ix - vsPML[VERTEX5].startX + 1, 		pml, delta.x );
		float dampingY = DAMPING( Vs_wave, vsPML[VERTEX5].endY - iy, 			pml, delta.y );
		float dampingZ = DAMPING( Vs_wave, iz - vsPML[VERTEX5].startZ + 1, 		pml, delta.z );

		VelocityPMLUpdate( VERTEX5 ) 
	}
	else if( AreaCondition( ix, iy, iz, vsPML[VERTEX6].startX, vsPML[VERTEX6].startY, vsPML[VERTEX6].startZ, vsPML[VERTEX6].endX, vsPML[VERTEX6].endY, vsPML[VERTEX6].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, ix - vsPML[VERTEX6].startX + 1, 		pml, delta.x );
		float dampingY = DAMPING( Vs_wave, iy - vsPML[VERTEX6].startY + 1, 		pml, delta.y );
		float dampingZ = DAMPING( Vs_wave, vsPML[VERTEX6].endZ - iz, 			pml, delta.z );

		VelocityPMLUpdate( VERTEX6 ) 
	}
	else if( AreaCondition( ix, iy, iz, vsPML[VERTEX7].startX, vsPML[VERTEX7].startY, vsPML[VERTEX7].startZ, vsPML[VERTEX7].endX, vsPML[VERTEX7].endY, vsPML[VERTEX7].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, ix - vsPML[VERTEX7].startX + 1, 		pml, delta.x );
		float dampingY = DAMPING( Vs_wave, iy - vsPML[VERTEX7].startY + 1, 		pml, delta.y );
		float dampingZ = DAMPING( Vs_wave, iz - vsPML[VERTEX7].startZ + 1, 		pml, delta.z );

		VelocityPMLUpdate( VERTEX7 ) 
	}
	else if( AreaCondition( ix, iy, iz, vsPML[SURFACE0].startX, vsPML[SURFACE0].startY, vsPML[SURFACE0].startZ, vsPML[SURFACE0].endX, vsPML[SURFACE0].endY, vsPML[SURFACE0].endZ ) )
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, vsPML[SURFACE0].endX - ix, 			pml, delta.x );
		float dampingY = 0.0f;
		float dampingZ = 0.0f;

		VelocityPMLUpdate( SURFACE0 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[SURFACE1].startX, vsPML[SURFACE1].startY, vsPML[SURFACE1].startZ, vsPML[SURFACE1].endX, vsPML[SURFACE1].endY, vsPML[SURFACE1].endZ ) )
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, ix - vsPML[SURFACE1].startX + 1, 	pml, delta.x );
		float dampingY = 0.0f;
		float dampingZ = 0.0f;

		VelocityPMLUpdate( SURFACE1 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[SURFACE2].startX, vsPML[SURFACE2].startY, vsPML[SURFACE2].startZ, vsPML[SURFACE2].endX, vsPML[SURFACE2].endY, vsPML[SURFACE2].endZ ) )
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingY = DAMPING( Vs_wave, vsPML[SURFACE2].endY - iy, 			pml, delta.y );
		float dampingX = 0.0f;
		float dampingZ = 0.0f;

		VelocityPMLUpdate( SURFACE2 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[SURFACE3].startX, vsPML[SURFACE3].startY, vsPML[SURFACE3].startZ, vsPML[SURFACE3].endX, vsPML[SURFACE3].endY, vsPML[SURFACE3].endZ ) )
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingY = DAMPING( Vs_wave, iy - vsPML[SURFACE3].startY + 1, 	pml, delta.y );
		float dampingX = 0.0f;
		float dampingZ = 0.0f;

		VelocityPMLUpdate( SURFACE3 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[SURFACE4].startX, vsPML[SURFACE4].startY, vsPML[SURFACE4].startZ, vsPML[SURFACE4].endX, vsPML[SURFACE4].endY, vsPML[SURFACE4].endZ ) )
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingZ = DAMPING( Vs_wave, vsPML[SURFACE4].endZ - iz, 			pml, delta.z );
		float dampingX = 0.0f;
		float dampingY = 0.0f;

		VelocityPMLUpdate( SURFACE4 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[SURFACE5].startX, vsPML[SURFACE5].startY, vsPML[SURFACE5].startZ, vsPML[SURFACE5].endX, vsPML[SURFACE5].endY, vsPML[SURFACE5].endZ ) )
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingZ = DAMPING( Vs_wave, iz - vsPML[SURFACE5].startZ + 1, 	pml, delta.z );
		float dampingX = 0.0f;
		float dampingY = 0.0f;

		VelocityPMLUpdate( SURFACE5 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[EDGE0].startX, vsPML[EDGE0].startY, vsPML[EDGE0].startZ, vsPML[EDGE0].endX, vsPML[EDGE0].endY, vsPML[EDGE0].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, vsPML[EDGE0].endX - ix, 				pml, delta.x );
		float dampingY = DAMPING( Vs_wave, vsPML[EDGE0].endY - iy, 				pml, delta.y );
		float dampingZ = 0.0f;

		VelocityPMLUpdate( EDGE0 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[EDGE1].startX, vsPML[EDGE1].startY, vsPML[EDGE1].startZ, vsPML[EDGE1].endX, vsPML[EDGE1].endY, vsPML[EDGE1].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, vsPML[EDGE1].endX - ix, 				pml, delta.x );
		float dampingY = DAMPING( Vs_wave, iy - vsPML[EDGE1].startY + 1, 		pml, delta.y );
		float dampingZ = 0.0f;

		VelocityPMLUpdate( EDGE1 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[EDGE2].startX, vsPML[EDGE2].startY, vsPML[EDGE2].startZ, vsPML[EDGE2].endX, vsPML[EDGE2].endY, vsPML[EDGE2].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, vsPML[EDGE2].endX - ix, 				pml, delta.x );
		float dampingZ = DAMPING( Vs_wave, vsPML[EDGE2].endZ - iz, 				pml, delta.z );
		float dampingY = 0.0f;

		VelocityPMLUpdate( EDGE2 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[EDGE3].startX, vsPML[EDGE3].startY, vsPML[EDGE3].startZ, vsPML[EDGE3].endX, vsPML[EDGE3].endY, vsPML[EDGE3].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, vsPML[EDGE3].endX - ix, 				pml, delta.x );
		float dampingZ = DAMPING( Vs_wave, iz - vsPML[EDGE3].startZ + 1, 		pml, delta.z );
		float dampingY = 0.0f;

		VelocityPMLUpdate( EDGE3 )
	}
	else if( AreaCondition( ix, iy, iz, vsPML[EDGE4].startX, vsPML[EDGE4].startY, vsPML[EDGE4].startZ, vsPML[EDGE4].endX, vsPML[EDGE4].endY, vsPML[EDGE4].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, ix - vsPML[EDGE4].startX + 1, 		pml, delta.x );
		float dampingY = DAMPING( Vs_wave, vsPML[EDGE4].endY - iy, 				pml, delta.y );
		float dampingZ = 0.0f;

		VelocityPMLUpdate( EDGE4 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[EDGE5].startX, vsPML[EDGE5].startY, vsPML[EDGE5].startZ, vsPML[EDGE5].endX, vsPML[EDGE5].endY, vsPML[EDGE5].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, ix - vsPML[EDGE5].startX + 1, 		pml, delta.x );
		float dampingY = DAMPING( Vs_wave, iy - vsPML[EDGE5].startY + 1, 		pml, delta.y );
		float dampingZ = 0.0f;

		VelocityPMLUpdate( EDGE5 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[EDGE6].startX, vsPML[EDGE6].startY, vsPML[EDGE6].startZ, vsPML[EDGE6].endX, vsPML[EDGE6].endY, vsPML[EDGE6].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, ix - vsPML[EDGE6].startX + 1, 		pml, delta.x );
		float dampingZ = DAMPING( Vs_wave, vsPML[EDGE6].endZ - iz, 				pml, delta.z );
		float dampingY = 0.0f;
		
		VelocityPMLUpdate( EDGE6 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[EDGE7].startX, vsPML[EDGE7].startY, vsPML[EDGE7].startZ, vsPML[EDGE7].endX, vsPML[EDGE7].endY, vsPML[EDGE7].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, ix - vsPML[EDGE7].startX + 1, 		pml, delta.x );
		float dampingZ = DAMPING( Vs_wave, iz - vsPML[EDGE7].startZ + 1, 		pml, delta.z );
		float dampingY = 0.0f;
		
		VelocityPMLUpdate( EDGE7 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[EDGE8].startX, vsPML[EDGE8].startY, vsPML[EDGE8].startZ, vsPML[EDGE8].endX, vsPML[EDGE8].endY, vsPML[EDGE8].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingY = DAMPING( Vs_wave, vsPML[EDGE8].endY - iy, 				pml, delta.y );
		float dampingZ = DAMPING( Vs_wave, vsPML[EDGE8].endZ - iz, 				pml, delta.z );
		float dampingX = 0.0f;
		
		VelocityPMLUpdate( EDGE8 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[EDGE9].startX, vsPML[EDGE9].startY, vsPML[EDGE9].startZ, vsPML[EDGE9].endX, vsPML[EDGE9].endY, vsPML[EDGE9].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingY = DAMPING( Vs_wave, vsPML[EDGE9].endY - iy, 				pml, delta.y );
		float dampingZ = DAMPING( Vs_wave, iz - vsPML[EDGE9].startZ + 1, 		pml, delta.z );
		float dampingX = 0.0f;
		
		VelocityPMLUpdate( EDGE9 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[EDGE10].startX, vsPML[EDGE10].startY, vsPML[EDGE10].startZ, vsPML[EDGE10].endX, vsPML[EDGE10].endY, vsPML[EDGE10].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingY = DAMPING( Vs_wave, iy - vsPML[EDGE10].startY + 1, 		pml, delta.y );
		float dampingZ = DAMPING( Vs_wave, vsPML[EDGE10].endZ - iz, 			pml, delta.z );
		float dampingX = 0.0f;

		VelocityPMLUpdate( EDGE10 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[EDGE11].startX, vsPML[EDGE11].startY, vsPML[EDGE11].startZ, vsPML[EDGE11].endX, vsPML[EDGE11].endY, vsPML[EDGE11].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingY = DAMPING( Vs_wave, iy - vsPML[EDGE11].startY + 1, 		pml, delta.y );
		float dampingZ = DAMPING( Vs_wave, iz - vsPML[EDGE11].startZ + 1, 		pml, delta.z );
		float dampingX = 0.0f;
		
		VelocityPMLUpdate( EDGE11 );
	}
}

 __global__ void stressUpdate( GRID grid, Delta delta, VelocityStress vsGPU, DifferenceCoefficient dc, Medium medGPU, int pml )//__global__ 
{
	int i;
	int ix = threadIdx.x + blockIdx.x * blockDim.x; 
	int iy = threadIdx.y + blockIdx.y * blockDim.y;
	int iz = blockIdx.z * blockDim.z;

	if ( AreaCondition( ix, iy, iz, CalculateAreaStart, CalculateAreaStart, CalculateAreaStart, CalculateAreaEndX, CalculateAreaEndY, CalculateAreaEndZ ) )
	{	
		
		float muBarXY = 0.0f, muBarXZ = 0.0f, muBarYZ = 0.0f;

		float diffX = 0.0f, diffY = 0.0f, diffZ = 0.0f;
		for (i = 0; i < order; i++)
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
		for (i = 0; i < order; i++)
		{
			diffY += dc.diff_coef[i] * (vsGPU.velocityX[INDEX(ix, iy + 1 + i, iz)] - vsGPU.velocityX[INDEX(ix, iy - i, iz)]);
			diffX += dc.diff_coef[i] * (vsGPU.velocityY[INDEX(ix + 1 + i, iy, iz)] - vsGPU.velocityY[INDEX(ix - i, iy, iz)]);
		}
		diffX /= delta.x;
		diffY /= delta.y;
		vsGPU.stressXY[INDEX(ix, iy, iz)] += delta.t * muBarXY * (diffY + diffX);

		diffX = 0.0f, diffZ = 0.0f;
		for (i = 0; i < order; i++)
		{
			diffZ += dc.diff_coef[i] * (vsGPU.velocityX[INDEX(ix, iy, iz + 1 + i)] - vsGPU.velocityX[INDEX(ix, iy, iz - i)]);
			diffX += dc.diff_coef[i] * (vsGPU.velocityZ[INDEX(ix + 1 + i, iy, iz)] - vsGPU.velocityZ[INDEX(ix - i, iy, iz)]);
		}
		diffX /= delta.x;
		diffZ /= delta.z;
		vsGPU.stressXZ[INDEX(ix, iy, iz)] += delta.t * muBarXZ * (diffZ + diffX);

		diffY = 0.0f, diffZ = 0.0f;
		for (i = 0; i < order; i++)
		{
			diffZ += dc.diff_coef[i] * (vsGPU.velocityY[INDEX(ix, iy, iz + 1 + i)] - vsGPU.velocityY[INDEX(ix, iy, iz - i)]);
			diffY += dc.diff_coef[i] * (vsGPU.velocityZ[INDEX(ix, iy + 1 + i, iz)] - vsGPU.velocityZ[INDEX(ix, iy - i, iz)]);
		}
		diffY /= delta.y;
		diffZ /= delta.z;
		vsGPU.stressYZ[INDEX(ix, iy, iz)] += delta.t * muBarYZ * (diffZ + diffY);

	}
	else if( AreaCondition( ix, iy, iz, vsPML[VERTEX0].startX, vsPML[VERTEX0].startY, vsPML[VERTEX0].startZ, vsPML[VERTEX0].endX, vsPML[VERTEX0].endY, vsPML[VERTEX0].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, vsPML[VERTEX0].endX - ix, 			pml, delta.x );
		float dampingY = DAMPING( Vs_wave, vsPML[VERTEX0].endY - iy, 			pml, delta.y );
		float dampingZ = DAMPING( Vs_wave, vsPML[VERTEX0].endZ - iz, 			pml, delta.z );

		StressPMLUpdate(  VERTEX0 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[VERTEX1].startX, vsPML[VERTEX1].startY, vsPML[VERTEX1].startZ, vsPML[VERTEX1].endX, vsPML[VERTEX1].endY, vsPML[VERTEX1].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, vsPML[VERTEX1].endX - ix, 			pml, delta.x );
		float dampingY = DAMPING( Vs_wave, vsPML[VERTEX1].endY - iy, 			pml, delta.y );
		float dampingZ = DAMPING( Vs_wave, iz - vsPML[VERTEX1].startZ + 1, 		pml, delta.z );

		StressPMLUpdate( VERTEX1 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[VERTEX2].startX, vsPML[VERTEX2].startY, vsPML[VERTEX2].startZ, vsPML[VERTEX2].endX, vsPML[VERTEX2].endY, vsPML[VERTEX2].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, vsPML[VERTEX2].endX - ix, 			pml, delta.x );
		float dampingY = DAMPING( Vs_wave, iy - vsPML[VERTEX2].startY + 1, 		pml, delta.y );
		float dampingZ = DAMPING( Vs_wave, vsPML[VERTEX2].endZ - iz, 			pml, delta.z );

		StressPMLUpdate( VERTEX2 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[VERTEX3].startX, vsPML[VERTEX3].startY, vsPML[VERTEX3].startZ, vsPML[VERTEX3].endX, vsPML[VERTEX3].endY, vsPML[VERTEX3].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, vsPML[VERTEX3].endX - ix, 			pml, delta.x );
		float dampingY = DAMPING( Vs_wave, iy - vsPML[VERTEX3].startY + 1, 		pml, delta.y );
		float dampingZ = DAMPING( Vs_wave, iz - vsPML[VERTEX3].startZ + 1, 		pml, delta.z );

		StressPMLUpdate( VERTEX3 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[VERTEX4].startX, vsPML[VERTEX4].startY, vsPML[VERTEX4].startZ, vsPML[VERTEX4].endX, vsPML[VERTEX4].endY, vsPML[VERTEX4].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, ix - vsPML[VERTEX4].startX + 1, 		pml, delta.x );
		float dampingY = DAMPING( Vs_wave, vsPML[VERTEX4].endY - iy, 			pml, delta.y );
		float dampingZ = DAMPING( Vs_wave, vsPML[VERTEX4].endZ - iz, 			pml, delta.z );

		StressPMLUpdate( VERTEX4 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[VERTEX5].startX, vsPML[VERTEX5].startY, vsPML[VERTEX5].startZ, vsPML[VERTEX5].endX, vsPML[VERTEX5].endY, vsPML[VERTEX5].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, ix - vsPML[VERTEX5].startX + 1, 		pml, delta.x );
		float dampingY = DAMPING( Vs_wave, vsPML[VERTEX5].endY - iy, 			pml, delta.y );
		float dampingZ = DAMPING( Vs_wave, iz - vsPML[VERTEX5].startZ + 1, 		pml, delta.z );

		StressPMLUpdate( VERTEX5 ) 
	}
	else if( AreaCondition( ix, iy, iz, vsPML[VERTEX6].startX, vsPML[VERTEX6].startY, vsPML[VERTEX6].startZ, vsPML[VERTEX6].endX, vsPML[VERTEX6].endY, vsPML[VERTEX6].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, ix - vsPML[VERTEX6].startX + 1, 		pml, delta.x );
		float dampingY = DAMPING( Vs_wave, iy - vsPML[VERTEX6].startY + 1, 		pml, delta.y );
		float dampingZ = DAMPING( Vs_wave, vsPML[VERTEX6].endZ - iz, 			pml, delta.z );

		StressPMLUpdate( VERTEX6 ) 
	}
	else if( AreaCondition( ix, iy, iz, vsPML[VERTEX7].startX, vsPML[VERTEX7].startY, vsPML[VERTEX7].startZ, vsPML[VERTEX7].endX, vsPML[VERTEX7].endY, vsPML[VERTEX7].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, ix - vsPML[VERTEX7].startX + 1, 		pml, delta.x );
		float dampingY = DAMPING( Vs_wave, iy - vsPML[VERTEX7].startY + 1, 		pml, delta.y );
		float dampingZ = DAMPING( Vs_wave, iz - vsPML[VERTEX7].startZ + 1, 		pml, delta.z );

		StressPMLUpdate( VERTEX7 ) 
	}
	else if( AreaCondition( ix, iy, iz, vsPML[SURFACE0].startX, vsPML[SURFACE0].startY, vsPML[SURFACE0].startZ, vsPML[SURFACE0].endX, vsPML[SURFACE0].endY, vsPML[SURFACE0].endZ ) )
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, vsPML[SURFACE0].endX - ix, 			pml, delta.x );
		float dampingY = 0.0f;
		float dampingZ = 0.0f;

		StressPMLUpdate( SURFACE0 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[SURFACE1].startX, vsPML[SURFACE1].startY, vsPML[SURFACE1].startZ, vsPML[SURFACE1].endX, vsPML[SURFACE1].endY, vsPML[SURFACE1].endZ ) )
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, ix - vsPML[SURFACE1].startX + 1, 	pml, delta.x );
		float dampingY = 0.0f;
		float dampingZ = 0.0f;

		StressPMLUpdate( SURFACE1 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[SURFACE2].startX, vsPML[SURFACE2].startY, vsPML[SURFACE2].startZ, vsPML[SURFACE2].endX, vsPML[SURFACE2].endY, vsPML[SURFACE2].endZ ) )
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingY = DAMPING( Vs_wave, vsPML[SURFACE2].endY - iy, 			pml, delta.y );
		float dampingX = 0.0f;
		float dampingZ = 0.0f;

		StressPMLUpdate( SURFACE2 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[SURFACE3].startX, vsPML[SURFACE3].startY, vsPML[SURFACE3].startZ, vsPML[SURFACE3].endX, vsPML[SURFACE3].endY, vsPML[SURFACE3].endZ ) )
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingY = DAMPING( Vs_wave, iy - vsPML[SURFACE3].startY + 1, 	pml, delta.y );
		float dampingX = 0.0f;
		float dampingZ = 0.0f;

		StressPMLUpdate( SURFACE3 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[SURFACE4].startX, vsPML[SURFACE4].startY, vsPML[SURFACE4].startZ, vsPML[SURFACE4].endX, vsPML[SURFACE4].endY, vsPML[SURFACE4].endZ ) )
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingZ = DAMPING( Vs_wave, vsPML[SURFACE4].endZ - iz, 			pml, delta.z );
		float dampingX = 0.0f;
		float dampingY = 0.0f;

		StressPMLUpdate( SURFACE4 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[SURFACE5].startX, vsPML[SURFACE5].startY, vsPML[SURFACE5].startZ, vsPML[SURFACE5].endX, vsPML[SURFACE5].endY, vsPML[SURFACE5].endZ ) )
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingZ = DAMPING( Vs_wave, iz - vsPML[SURFACE5].startZ + 1, 	pml, delta.z );
		float dampingX = 0.0f;
		float dampingY = 0.0f;

		StressPMLUpdate( SURFACE5 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[EDGE0].startX, vsPML[EDGE0].startY, vsPML[EDGE0].startZ, vsPML[EDGE0].endX, vsPML[EDGE0].endY, vsPML[EDGE0].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, vsPML[EDGE0].endX - ix, 				pml, delta.x );
		float dampingY = DAMPING( Vs_wave, vsPML[EDGE0].endY - iy, 				pml, delta.y );
		float dampingZ = 0.0f;

		StressPMLUpdate( EDGE0 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[EDGE1].startX, vsPML[EDGE1].startY, vsPML[EDGE1].startZ, vsPML[EDGE1].endX, vsPML[EDGE1].endY, vsPML[EDGE1].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, vsPML[EDGE1].endX - ix, 				pml, delta.x );
		float dampingY = DAMPING( Vs_wave, iy - vsPML[EDGE1].startY + 1, 		pml, delta.y );
		float dampingZ = 0.0f;

		StressPMLUpdate( EDGE1 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[EDGE2].startX, vsPML[EDGE2].startY, vsPML[EDGE2].startZ, vsPML[EDGE2].endX, vsPML[EDGE2].endY, vsPML[EDGE2].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, vsPML[EDGE2].endX - ix, 				pml, delta.x );
		float dampingZ = DAMPING( Vs_wave, vsPML[EDGE2].endZ - iz, 				pml, delta.z );
		float dampingY = 0.0f;

		StressPMLUpdate( EDGE2 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[EDGE3].startX, vsPML[EDGE3].startY, vsPML[EDGE3].startZ, vsPML[EDGE3].endX, vsPML[EDGE3].endY, vsPML[EDGE3].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, vsPML[EDGE3].endX - ix, 				pml, delta.x );
		float dampingZ = DAMPING( Vs_wave, iz - vsPML[EDGE3].startZ + 1, 		pml, delta.z );
		float dampingY = 0.0f;

		StressPMLUpdate( EDGE3 )
	}
	else if( AreaCondition( ix, iy, iz, vsPML[EDGE4].startX, vsPML[EDGE4].startY, vsPML[EDGE4].startZ, vsPML[EDGE4].endX, vsPML[EDGE4].endY, vsPML[EDGE4].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, ix - vsPML[EDGE4].startX + 1, 		pml, delta.x );
		float dampingY = DAMPING( Vs_wave, vsPML[EDGE4].endY - iy, 				pml, delta.y );
		float dampingZ = 0.0f;

		StressPMLUpdate( EDGE4 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[EDGE5].startX, vsPML[EDGE5].startY, vsPML[EDGE5].startZ, vsPML[EDGE5].endX, vsPML[EDGE5].endY, vsPML[EDGE5].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, ix - vsPML[EDGE5].startX + 1, 		pml, delta.x );
		float dampingY = DAMPING( Vs_wave, iy - vsPML[EDGE5].startY + 1, 		pml, delta.y );
		float dampingZ = 0.0f;

		StressPMLUpdate( EDGE5 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[EDGE6].startX, vsPML[EDGE6].startY, vsPML[EDGE6].startZ, vsPML[EDGE6].endX, vsPML[EDGE6].endY, vsPML[EDGE6].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, ix - vsPML[EDGE6].startX + 1, 		pml, delta.x );
		float dampingZ = DAMPING( Vs_wave, vsPML[EDGE6].endZ - iz, 				pml, delta.z );
		float dampingY = 0.0f;
		
		StressPMLUpdate( EDGE6 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[EDGE7].startX, vsPML[EDGE7].startY, vsPML[EDGE7].startZ, vsPML[EDGE7].endX, vsPML[EDGE7].endY, vsPML[EDGE7].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingX = DAMPING( Vs_wave, ix - vsPML[EDGE7].startX + 1, 		pml, delta.x );
		float dampingZ = DAMPING( Vs_wave, iz - vsPML[EDGE7].startZ + 1, 		pml, delta.z );
		float dampingY = 0.0f;
		
		StressPMLUpdate( EDGE7 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[EDGE8].startX, vsPML[EDGE8].startY, vsPML[EDGE8].startZ, vsPML[EDGE8].endX, vsPML[EDGE8].endY, vsPML[EDGE8].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingY = DAMPING( Vs_wave, vsPML[EDGE8].endY - iy, 				pml, delta.y );
		float dampingZ = DAMPING( Vs_wave, vsPML[EDGE8].endZ - iz, 				pml, delta.z );
		float dampingX = 0.0f;
		
		StressPMLUpdate( EDGE8 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[EDGE9].startX, vsPML[EDGE9].startY, vsPML[EDGE9].startZ, vsPML[EDGE9].endX, vsPML[EDGE9].endY, vsPML[EDGE9].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingY = DAMPING( Vs_wave, vsPML[EDGE9].endY - iy, 				pml, delta.y );
		float dampingZ = DAMPING( Vs_wave, iz - vsPML[EDGE9].startZ + 1, 		pml, delta.z );
		float dampingX = 0.0f;
		
		StressPMLUpdate( EDGE9 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[EDGE10].startX, vsPML[EDGE10].startY, vsPML[EDGE10].startZ, vsPML[EDGE10].endX, vsPML[EDGE10].endY, vsPML[EDGE10].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingY = DAMPING( Vs_wave, iy - vsPML[EDGE10].startY + 1, 		pml, delta.y );
		float dampingZ = DAMPING( Vs_wave, vsPML[EDGE10].endZ - iz, 			pml, delta.z );
		float dampingX = 0.0f;

		StressPMLUpdate( EDGE10 );
	}
	else if( AreaCondition( ix, iy, iz, vsPML[EDGE11].startX, vsPML[EDGE11].startY, vsPML[EDGE11].startZ, vsPML[EDGE11].endX, vsPML[EDGE11].endY, vsPML[EDGE11].endZ ) ) 
	{
		float Vs_wave = sqrt( medGPU.mu[INDEX(ix, iy, iz)] * medGPU.buoyancy[INDEX(ix, iy, iz)] );
		float dampingY = DAMPING( Vs_wave, iy - vsPML[EDGE11].startY + 1, 		pml, delta.y );
		float dampingZ = DAMPING( Vs_wave, iz - vsPML[EDGE11].startZ + 1, 		pml, delta.z );
		float dampingX = 0.0f;
		
		StressPMLUpdate( EDGE11 );
	}
}
/*
void ElasticWaveEquation::snapShootSlice(int it)
{
	char fileName[6][100] = { 0 };

	ofstream * snapshot = new ofstream[6];

	sprintf(fileName[0], "snapshotSliceXX_YoZ_%d.txt", it);
	sprintf(fileName[1], "snapshotSliceXX_XoZ_%d.txt", it);
	sprintf(fileName[2], "snapshotSliceXX_XOY_%d.txt", it);
	sprintf(fileName[3], "surfSliceXX_YoZ_%d.txt", it);
	sprintf(fileName[4], "surfSliceXX_XoZ_%d.txt", it);
	sprintf(fileName[5], "surfSliceXX_XOY_%d.txt", it);
	// sprintf(fileName[0], "velocitySliceXX_YoZ_%d.txt", it);
	// sprintf(fileName[1], "velocitySliceXX_XoZ_%d.txt", it);
	// sprintf(fileName[2], "velocitySliceXX_XoY_%d.txt", it);
	snapshot[0].open(fileName[0]);
	snapshot[1].open(fileName[1]);
	snapshot[2].open(fileName[2]);
	snapshot[3].open(fileName[3]);
	snapshot[4].open(fileName[4]);
	snapshot[5].open(fileName[5]);
	int ix, iy, iz;
	
	for( iz = 0; iz < grid.zLength; ++iz ) {
		
		for( iy = 0; iy < grid.yLength; ++iy) {
			
			//snapshot[0] << vsCPU.stressXX[INDEX( sourceLocationGPU.x, iy, iz)] << endl;
			snapshot[0] << vsCPU.stressXX[INDEX(sourceLocationCPU.x, iy, iz)] << endl;
			snapshot[3] << vsCPU.stressXX[INDEX(CalculateAreaEndX, iy, iz)] << endl;

		}

	}

	for( iz = 0; iz < grid.zLength; ++iz ) {
		
		for( ix = 0; ix < grid.xLength; ++ix) {
			
			snapshot[1] << vsCPU.stressXX[INDEX(ix, sourceLocationCPU.y, iz)] << endl;
			snapshot[4] << vsCPU.stressXX[INDEX(ix, CalculateAreaEndY, iz)] << endl;

		}

	}
	

	for( iy = 0; iy < grid.yLength; ++iy ) {
		
		for( ix = 0; ix < grid.xLength; ++ix) {
			
			snapshot[2] << vsCPU.stressXX[INDEX(ix, iy, sourceLocationCPU.z)] << endl;
			snapshot[5] << vsCPU.stressXX[INDEX(ix, CalculateAreaEndZ, iz)] << endl;

		}

	}
	
	snapshot[0].close( );
	snapshot[1].close( );
	snapshot[2].close( );
	snapshot[3].close( );
	snapshot[4].close( );
	snapshot[5].close( );
	delete[] snapshot;

}

void ElasticWaveEquation::snapShootVolum(int it)
{
	char fileName[6][100] = { 0 };
	ofstream * snapshot = new ofstream[6];
	int N = 1;
	for (size_t i = 0; i < N; i++)
	{
		sprintf(fileName[i], "snapshot_%d_%d.txt", i, it);
		cout << fileName[i] << endl;
		snapshot[i].open(fileName[i]);
	}
	for (size_t i = 0; i < 6; i++)
	{
		if (!snapshot[i].is_open()) 
			cout << "Can't open the file" << endl;

	}
	int ix, iy, iz;
	for (iz = order; iz < grid.zLength - order; iz++)
	{
		for (iy = order; iy < grid.yLength - order; iy++)
		{
			for (ix = order; ix < grid.xLength - order; ix++)
			{
				snapshot[0] << vsCPU.stressXX[INDEX(ix, iy, iz)] << endl;
			}
		}
	}

	//for (size_t i = 0; i < N; i++)
	//{
		//snapshot[i].close( );
	//}

	delete[] snapshot;

}
*/

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


void ElasticWaveEquation::stationDataOutput()
{
	char fileName[100] = "stationData.txt";
	ofstream stationFile;
	stationFile.open( fileName );
	for(unsigned i = 0; i < timeLength; ++i) {
		stationFile << stationDataCPU[i] << endl;
	}

	stationFile.close();
}

