// dllmain.cpp : Defines the entry point for the DLL application.
#include "stdafx.h"
#include "newLattice.h"
#include "functions.hpp"
#include "alphaExpansion.h"

#define MRF_API extern "C" __declspec(dllexport) 

Lattice *lat = NULL;
double *rndField = NULL;
AlphaExpansion *alpha = NULL;
int *image = NULL;
int *priorImage = NULL;

MRF_API void CreateLattice(int width, int length, int height, bool periodic);
MRF_API void InitializeLattice(int* initLabel, int expandLabel, double **nbCost, double **priorCost);
MRF_API bool Optimize(double freq);
MRF_API void UpdateLabels(int *labels, int expandLabel);
MRF_API void CleanUp();
MRF_API void InitializeLatticeRFIM(long int seed, double delta, double hfield);
MRF_API double GetIsingEnergy();
MRF_API double GetIsingMagnetization();
MRF_API int* AlphaExpansionPoisson(int width, int length, int height, double prioStrength, int nlabel, double &oldEnergy, double &newEnergy, int &nMoves);
MRF_API int* GetImagePointer(int width, int height);
MRF_API int* GetPriorImagePointer(int width, int height);

BOOL APIENTRY DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
					 )
{
	switch (ul_reason_for_call)
	{
	case DLL_PROCESS_ATTACH:
	case DLL_THREAD_ATTACH:
	case DLL_THREAD_DETACH:
	case DLL_PROCESS_DETACH:
		break;
	}
	return TRUE;
}

void CreateLattice(int width, int length, int height, bool periodic) {
	if (width > 1 && length > 1 && height > 0) lat = new Lattice(width, length, height, periodic);
}

void InitializeLattice(int* initLabel, int expandLabel, double **nbCost, double **priorCost) {
	if (lat != NULL) {
		lat->Init(initLabel, expandLabel, nbCost, priorCost);
	}
}

bool Optimize(double freq) {
	if (lat != NULL) {
		return lat->GndState(freq);
	} else {
		return false;
	}
}

void UpdateLabels(int *labels, int expandLabel) {
	if (lat != NULL) lat->UpdateLabel(labels, expandLabel);
}

void CleanUp() {
	delete lat;
	lat = NULL;
	delete rndField;
	rndField = NULL;
	delete image;
	image = NULL;
	delete alpha;
	alpha = NULL;
	delete priorImage;
	priorImage = NULL;
}

void InitializeLatticeRFIM(long int seed, double delta, double hfield) {
	if (lat == NULL) return;
	if (rndField != NULL) {
		delete rndField;
		rndField = NULL;
	}
	rndField = new double[lat->SIZE];

	double equd = delta;
	double equh = hfield;
	if (hfield < 0) {
		equd = -equd;
		equh = -equh;
	}

	if (seed < 0) {
		seed = -seed;
	}
	if (seed == 0) {
		seed = 1;
	}

	long int idum = -seed;
	for (int curr=0; curr<lat->SIZE; curr++) {
		rndField[curr] = NormalRnd(&idum) * equd + equh;
	}

	lat->IsingInit(rndField);
}

double GetIsingEnergy() {
	if (lat == NULL || rndField == NULL) return 0;
	return lat->IsingEnergy(rndField);
}

double GetIsingMagnetization() {
	if (lat == NULL || rndField == NULL) return 0;
	return lat->IsingMagnetization(rndField);
}

int* AlphaExpansionPoisson(int width, int length, int height, double prioStrength, int nlabel, double &oldEnergy, double &newEnergy, int &nMoves) {
	if (alpha != NULL) delete alpha;

	double scale;
	int baselevel;

	AlphaExpansion::ScaleImage(image, width*length*height, nlabel, scale, baselevel);

	alpha = new AlphaExpansion(width, length, height, image, nlabel, baselevel, scale, priorImage, poissonProbEst, poissonProb, prioStrength);
	nMoves = alpha->Expansion(2, 10, 100, oldEnergy, newEnergy);
	return alpha->GetUpdate();
}

int* GetImagePointer(int width, int height) {
	if (image != NULL) delete image;
	image = new int[width * height];
	return image;
}

int* GetPriorImagePointer(int width, int height) {
	if (priorImage != NULL) delete priorImage;
	priorImage = new int[width * height];
	return priorImage;
}