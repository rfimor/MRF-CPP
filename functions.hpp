#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#define REPS 1.2e-100

#define IA 16807
#define IM 2147483647
#define AM 4.6566128752457969241057508272e-10

#define IQ 127773
#define IR 2836
#define NTAB 32

#define NDIV 67108864
  
#define IB1 1
#define IB2 2
#define IB5 16
#define IB27 67108864 
#define MASK 19

#define PI 3.14159265359

void stat(const double *data, const int dim, double & meanx, double & stderrx);
	
template <class T>
void make2DArray(T ** &a,int d1, int d2)
{
        a = new T * [d1];
        for (int i = 0; i < d1;i ++) {
             a[i] = new T [d2];
        }
}

template <class T>
void delete2DArray(T ** &a, int d1, int d2)
{
        for (int i = 0; i < d1;i ++) {
             delete [] a[i];
        }
        delete [] a;
}

template <class T>
void make3DArray(T *** &a,int d1, int d2, int d3)
{
        a = new T ** [d1];
        for (int i = 0; i < d1;i ++) {
             a[i] = new T* [d2];
             for (int j = 0; j < d2; j++) {
                  a[i][j] = new T[d3];
             }
        }
}

template <class T>
void delete3DArray(T *** &a, int d1, int d2, int d3)
{
        for (int i = 0; i < d1;i ++) {
             for (int j = 0; j < d2; j++) {
                  delete [] a[i][j];
             }
             delete [] a[i];
        }
        delete [] a;
}


double ran1(long int* idum);
double NormalRnd(long int *idum);
int irbit2(unsigned long int *seed);

template <class T>
T getInput(const char * mes);

void bunching(const double *dataset, const int nsample, const int ndiv,
              double func(const double *,const int), double &meanb, double &stderrb);

void bootstrap(const double *dataset, const double *weight, const int nsample, const int nboot,
               double func(const double *,const double *,const int), double &meanb, double &stderrb);

double poissonProb(int val, int lambda);
double poissonProbEst(int val1, int val2);
double logfac(int val);

#endif
