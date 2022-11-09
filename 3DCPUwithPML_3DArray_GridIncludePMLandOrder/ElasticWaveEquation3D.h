/***************************************************************************************/
/******************************************2018-12-21***********************************/
/************************************Author:Wenqiang Wang*******************************/
/**********************Southern University of Science and Technology********************/
/***************************************************************************************/

#ifndef ELASTIC_WAVE_EQUATION3D_H
#define ELASTIC_WAVE_EQUATION3D_H
#include "ElasticMacro.h"
class ElasticWaveEquation3D
{
public:
	ElasticWaveEquation3D(
		GRID grid, int timeLength, int pml,
		float *** buoyancy, float *** lambda, float *** mu,
		SOURCELOCATION sourceLocation,
		float deltaX = 5.0f, float deltaY = 5.0f, float deltaZ = 5.0f, float deltaT = 0.0005f,
		float mainFrequncy = 30.0f );

	~ElasticWaveEquation3D();
	
	void initialize();
	void propagate();
	void snapshoot(int timeLength);
	float sourceFunction(int timeIndex);
	void loadSeismicSource(int timeIterator);
	void velocityUpdate();
	void stressUpdate();
	void calculateDamping();

private:
	int timeLength;
	float deltaX, deltaY, deltaZ, deltaT;
	float mainFrequncy;
	float c1, c2;

	GRID grid;
	int pml;
	SOURCELOCATION sourceLocation;

	DEF_3D_POINTER(float, buoyancy);
	DEF_3D_POINTER(float, lambda);
	DEF_3D_POINTER(float, mu);
	
	DEF_3D_VelocityStress(float,   );
	DEF_3D_VelocityStress(float, _x);
	DEF_3D_VelocityStress(float, _y);
	DEF_3D_VelocityStress(float, _z);

	DEF_3D_POINTER(float, damping_x);
	DEF_3D_POINTER(float, damping_y);
	DEF_3D_POINTER(float, damping_z);


};
#endif
