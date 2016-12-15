#ifndef NEWLATTICE_H
#define NEWLATTICE_H

#include <stdlib.h>
#include "structure.h"

const int MAX_EDGE = 8;

//------------------------------------------------------------------------
class Vertex;

class Edge //flow & capacity coming into the vertex
{
  public:
    Vertex* end;
    double capacity;
    double residual;
    Edge* reverse;

    Edge(void);
    ~Edge(void) {}
    double Flow(void) {return capacity - residual;}
} ;

//-------------------------------------------------------------------------

class Vertex
{
  public:
     Vertex* next;
     Vertex* prev;

     double excess;
     int height;

     Edge* adj[MAX_EDGE];

     Vertex(void);
     ~Vertex(void);

} ;

class VertexList
{
  public:
     Vertex* head;
     Vertex* tail;

     VertexList(void);
     ~VertexList(void);
     void InsertFront(Vertex* v);
     void InsertRear(Vertex* v);
     void Remove(Vertex* v);
     void Delete(Vertex* v);
     void Empty(void);
} ;

class Lattice
{
  public:
     const int W;
	 const int L;
	 const int H;
     const int SIZE;

     Lattice(int width, int length, int height, bool periodic);
     ~Lattice(void);
	
	void IsingInit(double* mField);
     void Init(int initLabel[], int expandLabel, double **nbCost, double **priorCost);

  private:

     Vertex** site;
     Vertex* sink;

     double eps;
/*
     unsigned long int nPush;        //number of push operations
     unsigned long int nLift;        // number of lift operations
     unsigned long int nDischarge;   //number of discharges
     unsigned long int nGap;         //number of gap relabling
     unsigned long int nUpdate;
     unsigned long int nGapNodes;
*/

     const int HIGHEST;
     unsigned long int workSinceUpdate;
     const int ALPHA;
     const int BETA;
     unsigned long int updateLimit;
     int maxActive;
     int maxHeight;

     VertexList **active;
     VertexList **inactive;

     bool gndFound;

     void GlobalUpdate(void);
     void Discharge(Vertex *u);

	 void distributeWeight(int i, int j, Edge *adj, int initLabel[], int expandLabel, double **nbCost, bool bidirection); 

  public:
     bool GndState(double frequency);
	void UpdateLabel(int *label, int expandLabel);
	double IsingEnergy(double *prior);
	double IsingMagnetization(double *prior);
} ;

#endif

 
