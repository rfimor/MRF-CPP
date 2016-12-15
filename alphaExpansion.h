#ifndef ALPHAEXPANSION_H
#define ALPHAEXPANSION_H

#include "newLattice.h"

class AlphaExpansion 
{
	public:
		const int W;
		const int L;
		const int H;

		const double EPS;

		AlphaExpansion(int width, int length, int height, int *image, int nlabels, double **nbCost, double **priorCost);
		AlphaExpansion(int width, int length, int height, int *image, int nlabels, int *priorImage, double (*nbCost)(int, int), double (*priorCost)(int, int), double priorStrength);
		~AlphaExpansion();

		double OneMove(int expandLabel);
		int* GetUpdate();
		double GetEnergy();
		int Expansion(int minRounds, int maxRounds, long int seed, double &oldEnergy, double &energy);

	private:
		Lattice *lat;

		int *image;
		int *newImage;
		int nlabel;
		double **nbCost;
		double **priorCost;
		bool selfCreatedCost;

		int baseLevel;
		double scale;
		bool isScaled;

		void scaleImage(int *image, int length);
		int scaleBack(int x);
};

#endif