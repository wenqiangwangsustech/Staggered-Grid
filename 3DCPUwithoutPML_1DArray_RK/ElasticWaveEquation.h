/***************************************************************************************/
/******************************************2019-01-24***********************************/
/************************************Author:Wenqiang Wang*******************************/
/**********************Southern University of Science and Technology********************/
/***************************************************************************************/
#ifndef ELASTICWAVEEQUATION_H
#define ELASTICWAVEEQUATION_H
#include "ElasticMacro.h"
class ElasticWaveEquation
{
public:
	DEF_3D_VelocityStress(float)
	float * buoyancy;
	float * lambda;
	float * mu;
	GRID grid;
	Delta delta;
	float * diff_coef;
	float mainFrequncy;
	int order;
	int timeLength;
	SOURCELOCATION sourceLocation;
	ElasticWaveEquation(
		int timeLength,
		GRID grid, Delta delta, 
		SOURCELOCATION sourceLocation, float mainFrequncy, float * diff_coef,
		float * buoyancy, float * lambda, float * mu, 
		int order);
	~ElasticWaveEquation();

	float sourceFunction(int timeIndex);
	void loadSeismicSource(int timeIterator);
	void snapshoot(int it);
	void propagate();
	void run();
	void snapShootSlice(int it);

};
#endif //ELASTICWAVEEQUATION_H
