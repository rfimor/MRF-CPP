#include "stdafx.h"
#include "alphaExpansion.h"
#include <set>
#include "functions.hpp"
#include <math.h>

AlphaExpansion::AlphaExpansion(int width, int length, int height, int *image, int nlabels, double **nbCost, double **priorCost) 
	: W(width), L(length), H(height), EPS(1e-8) {
	lat = new Lattice(width, length, height, false);

	this->nbCost = nbCost;
	this->priorCost = priorCost;
	this->nlabel = nlabels;
	this->image = image;

	selfCreatedCost = false;

	newImage = new int[lat->SIZE];
	for (int i=0; i<lat->SIZE; i++) newImage[i] = image[i];
}

//values in "labels" must be from 0 to nlabel - 1, baseLevel should be higher than 1
//"Image" is scaled to [0, nlabel-1]. "priorImage" is not scaled
AlphaExpansion::AlphaExpansion(int width, int length, int height, int *image, int nlabels, int baseLevel, double scale, int *priorImage, double (*nbCost)(int, int), double (*priorCost)(int, int), double priorStrength) 
	: W(width), L(length), H(height), EPS(1e-8) {
	lat = new Lattice(width, length, height, false);

	nlabel = nlabels;
	this->image = image;

	make2DArray(this->nbCost, nlabel, nlabel);
	make2DArray(this->priorCost, nlabel, lat->SIZE);

	for (int i=0; i < nlabel; i++) {
		for (int j=0; j < nlabel; j++) {
			this->nbCost[i][j] = nbCost(scaleBack(i, scale, baseLevel),scaleBack(j, scale, baseLevel));
		}
		for (int j=0; j < lat->SIZE; j++) {
			this->priorCost[i][j] = priorStrength * priorCost(scaleBack(i, scale, baseLevel), priorImage[j]);
		}
	}

	newImage = new int[lat->SIZE];
	for (int i=0; i<lat->SIZE; i++) newImage[i] = image[i];
}

AlphaExpansion::~AlphaExpansion() {
	if (selfCreatedCost) {
		delete2DArray(this->nbCost, nlabel, nlabel);
		delete2DArray(this->priorCost, nlabel, lat->SIZE);
	}
	delete lat;
	delete [] newImage;
}

int* AlphaExpansion::GetUpdate() {
	return newImage;
}

double AlphaExpansion::GetEnergy() {
	double energy = 0;
	int curr = 0;
	for (int k=0;k<H;k++)
		for (int j=0;j<L;j++)
			for (int i=0;i<W;i++) {
				energy += priorCost[newImage[curr]][curr];
				
				int next = i < W-1 ? curr + 1 : curr + 1 - W;
				energy += nbCost[newImage[curr]][newImage[next]];
				
				next = j < L-1 ? curr + W : curr - W * (L-1);
				energy += nbCost[newImage[curr]][newImage[next]];

				next = k < H-1 ? curr + W * L : curr - W * L * (H-1);
				energy += nbCost[newImage[curr]][newImage[next]];

				curr++;
			}
	return energy;
}

double AlphaExpansion::OneMove(int expandLabel) {
	lat->Init(newImage, expandLabel, nbCost, priorCost);
	lat->GndState(0);
	lat->UpdateLabel(newImage, expandLabel);
	return GetEnergy();
}

int AlphaExpansion::Expansion(int minRounds, int maxRounds, long int seed, double &oldEnergy, double &energy) {
	std::set<int> labels;

	long int idum = -seed;
	if (idum > 0) idum = -idum;
	if (idum == 0) idum = -1;

	oldEnergy = GetEnergy();
	energy = oldEnergy;
	double dE;
	int cnt = 0;
	for (int rd = 0; rd < maxRounds; rd++) {
		for (int i=0; i<nlabel; i++) labels.insert(i);
		while (labels.size() > 0) {
			int pl;
			do {
				pl = (int)floor(ran1(&idum) * nlabel);
			} while (labels.find(pl) == labels.end());
			labels.erase(pl);

			dE = energy;
			energy = OneMove(pl);
			cnt++;

			if (energy - dE > -EPS && cnt > nlabel*minRounds) {
				return cnt;
			}
		}
	}
	return cnt;
}

void AlphaExpansion::ScaleImage(int *image, int length, int &nlabel, double &scale, int &baselevel) {
	int max = image[0];
	int min = max;
	for (int i=1; i<length; i++) {
		if (image[i] > max) max = image[i];
		else if (image[i] < min) min = image[i];
	}

	if ((max - min) < nlabel) {
		nlabel = max - min + 1;
		baselevel = min;
		scale = 1;
		for (int i=1; i<length; i++) image[i] -= min;
	} else {
		scale = (double)(nlabel - 1) / (max - min);
		for (int i=1; i<length; i++) image[i] = (int)((image[i] - min) * scale);
		baselevel = min;
	}
}
