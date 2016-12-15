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
		AlphaExpansion(int width, int length, int height, int *image, int nlabels, int baseLevel, double scale, int *priorImage, double (*nbCost)(int, int), double (*priorCost)(int, int), double priorStrength);
		~AlphaExpansion();

		double OneMove(int expandLabel);
		int* GetUpdate();
		double GetEnergy();
		int Expansion(int minRounds, int maxRounds, long int seed, double &oldEnergy, double &energy);

		static void ScaleImage(int *image, int length, int &nlabel, double &scale, int &baselevel);

	private:
		Lattice *lat;

		int *image;
		int *newImage;
		int nlabel;
		double **nbCost;
		double **priorCost;
		bool selfCreatedCost;

		int scaleBack(int x, double scale, int base) {
			double s = x / scale + base;
			int r = (int)s;
			return r >= s - 0.5 ? (r + 1) :r;
		}
};

#endif