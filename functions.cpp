#include "stdafx.h"
#include <math.h>
#include "functions.hpp"

void stat(const double *data, const int dim, double & meanx, double & stderrx)
{
    meanx = 0;
    stderrx = 0;

    for (int i=0; i<dim; i++)
    {
         meanx += data[i];
         stderrx += data[i] * data[i];
    }
    meanx = meanx/dim;
    if (dim == 1)
    {
        stderrx = meanx;
    }
    else
    {
         stderrx = sqrt((stderrx - dim * meanx * meanx) / (dim - 1));
    }
    stderrx = stderrx/ sqrt((double)dim);
}

double ran1(long int* idum)
{
  const double RNMX =(1.-REPS);
  
  int j;
  long int k;
  static long int iy = 0;
  static long int iv[NTAB];
  double temp;

  if (*idum<=0||!iy) {
      if (-(*idum)<1) *idum = 1;
      else *idum = -(*idum);
      for (j=NTAB+7;j>=0;j--) {
          k = (*idum)/IQ;
          *idum = IA*(*idum-k*IQ)-IR*k;
          if (*idum<0) *idum+=IM;
          if (j<NTAB) iv[j] = *idum;
      }
      iy=iv[0];
  }
  k = (*idum)/IQ;
  *idum = IA*(*idum-k*IQ)-IR*k;
  if (*idum<0) *idum+=IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j]=*idum;
  if ((temp=AM*iy)>RNMX) return RNMX;
  else return temp;
}

double NormalRnd(long int *idum)
{
  
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;

  if (*idum < 0) iset = 0;
  if (iset == 0) {
      do {v1 = 2.*ran1(idum)-1.0;
          v2 = 2.*ran1(idum)-1.0;
          rsq = v1*v1 + v2*v2;
      }
      while (rsq>=1. || rsq == 0);
      fac = sqrt(-2.*log(rsq)/rsq);
      gset = v1*fac;
      iset = 1;
      return v2*fac;
  }
  else {
      iset =0;
      return gset;
  }
}

int irbit2(unsigned long int *iseed) {

	
	if (*iseed & IB27) {
		*iseed = ((*iseed ^ MASK)<<1) | IB1;
		return -1;
	}
	else {
		*iseed <<= 1;
		return 1;
	}
}

void bootstrap(const double *dataset, const int nsample, const int nboot,
               double func(const double *,const int), double &meanb, double &stderrb)
{
     double *resample = new double[nboot];
     double *sampletemp = new double[nsample];

     long int seed = -(long int)(nsample+101);
     for (int i=0; i<nboot; i++)
     {
         for (int j=0; j<nsample; j++)
         {
             double rsam = ran1(&seed) * (nsample-1);
             int nsam = (int)rsam;
             if (rsam - nsam >= 0.5) nsam++;

             sampletemp[j] = dataset[nsam];
         }
         resample[i] = func(sampletemp,nsample);
     }
     stat(resample,nboot,meanb,stderrb);
     stderrb *= sqrt((double)nboot);
}

void bunching(const double *dataset, const int nsample, const int ndiv,
              double func(const double *,const int), double &meanb, double &stderrb)
{
     if (nsample <= ndiv) return;
     double *resample = new double[ndiv];

     int nbunch = (int)((double)nsample/ndiv);
     int ntail = nsample-(nbunch*(ndiv-1));

//     assert(ntail>=nbunch);

     double *temp = new double[ntail];
     int ndim;

     for (int i=0; i<ndiv; i++)
     {
         if (i != ndiv-1)
         {
            for (int j=0; j<nbunch; j++)
            {
                temp[j] = dataset[i*nbunch+j];
            }
            ndim = nbunch;
         }
         else
         {
            for (int j=0; j<ntail; j++)
            {
                temp[j] = dataset[nsample-ntail+j];
            }
            ndim = ntail;
         }

         resample[i] = func(temp,ndim);
     }

     stat(resample,ndiv,meanb,stderrb);
}

double poissonProb(int val, int lambda) {
	return logfac(val) - val * log((double)(lambda < 1 ? 1 : lambda)) + lambda;
}

double poissonProbEst(int val1, int val2) {
	double a = ((val1 + val2) / 2.0);
	int avg = (int)a;
	avg = a - avg >= 0.5 ? avg + 1 : avg;

	return 0.5 * (poissonProb(val1, avg) + poissonProb(val2, avg));
}

double logfac(int val) {
	double r = 1;
	if (val < 10) {
		for (int i=1; i<=val; i++) r = r*val;
		return log(r);
	} else {
		return val * (log((double)val) - 1) + 0.5 * log(2.0 * PI * val);
	}
}
