#include "stdafx.h"
#include <assert.h>

#include "newLattice.h"
#include "functions.hpp"

Edge::Edge() {
	end = NULL;
}

Vertex::Vertex()
{
	for (int i=0; i<MAX_EDGE; i++)
	{
		adj[i] = NULL;
	}
}

Vertex::~Vertex()
{
	for (int i=0; i<MAX_EDGE; i++)
	{
		delete adj[i];
	}
}
//--------------------------------------------------------------------------


VertexList::VertexList()
{
	head = NULL;
	tail = NULL;
}

void VertexList::InsertFront(Vertex* v)
{
	if (head == NULL)
	{
		tail = head = v;
		v->prev = NULL;
		v->next = NULL;
	}
	else
	{
		v->next = head;
		head = head->prev = v;
		v->prev = NULL;
	}
}

void VertexList::InsertRear(Vertex* v)
{
	if (head == NULL)
	{
		InsertFront(v);
	}
	else
	{
		v->prev = tail;
		tail = tail->next = v;
		v->next = NULL;
	}
}

void VertexList::Remove(Vertex* v)
{
	if (v->prev != NULL)
	{
		(v->prev)->next = v->next;
	}
	else
	{
		head = v->next;
	}

	if (v->next != NULL)
	{
		(v->next)->prev = v->prev;
	}
	else
	{
		tail = v->prev;
	}
}

void VertexList::Delete(Vertex* v)
{
	if (v->prev != NULL)
	{
		(v->prev)->next = v->next;
	}
	else
	{
		head = v->next;
	}

	if (v->next != NULL)
	{
		(v->next)->prev = v->prev;
	}
	else
	{
		tail = v->prev;
	}
	delete v;
}

void VertexList::Empty()
{
	while (head != NULL)
	{
		Remove(head);
	}
}

VertexList::~VertexList()
{
	Empty();
}

//----------------------------------------------------------------------------


Lattice::Lattice(int width, int length, int height, bool periodic)
	: W(width), L(length), H(height), SIZE(L*W*H), HIGHEST(SIZE+1), ALPHA(6), BETA(12)
{
	updateLimit = ALPHA*(SIZE+2) + (MAX_EDGE*SIZE);

	site = new Vertex*[SIZE];

	for (int i=0; i<SIZE; i++)
	{
		site[i] = new Vertex();
	}

	sink = new Vertex();
	sink->height = 0;
	sink->excess = 0;

	eps = 1e-20;

	active = new VertexList*[HIGHEST];
	inactive = new VertexList*[HIGHEST];

	for (int i=1; i<HIGHEST; i++)
	{
		active[i] = new VertexList();
		inactive[i] = new VertexList();
	}

	for (int i=0;i<W;i++)
		for (int j=0;j<L;j++)
			for (int k=0;k<H;k++)
			{
				int curr = i + j * W + k * W * L;

				site[curr]->adj[0] = new Edge();
				site[curr]->adj[0]->end = NULL;

				site[curr]->adj[MAX_EDGE-1] = new Edge();
				site[curr]->adj[MAX_EDGE-1]->end = sink;
				site[curr]->adj[MAX_EDGE-1]->reverse = NULL;

				if (i < W-1 || periodic) {
					int next = (i+1)%W + j * W + k * W * L;
					site[curr]->adj[1] = new Edge();
					site[next]->adj[4] = new Edge();
					site[curr]->adj[1]->end = site[next];
					site[curr]->adj[1]->reverse =  site[next]->adj[4];
					site[next]->adj[4]->end = site[curr];
					site[next]->adj[4]->reverse =  site[curr]->adj[1];
				}

				if (j < L -1 || periodic) {
					int next = i + (j+1)%L * W + k * W * L;
					site[curr]->adj[2] = new Edge();
					site[next]->adj[5] = new Edge();
					site[curr]->adj[2]->end = site[next];
					site[curr]->adj[2]->reverse =  site[next]->adj[5];
					site[next]->adj[5]->end = site[curr];
					site[next]->adj[5]->reverse =  site[curr]->adj[2];
				}

				if (H > 1) {
					if (k < H-1 || periodic) {
						int next = i + j * W + (k+1)%H * W * L;
						site[curr]->adj[3] = new Edge();
						site[next]->adj[6] = new Edge();
						site[curr]->adj[3]->end = site[next];
						site[curr]->adj[3]->reverse =  site[next]->adj[6];
						site[next]->adj[6]->end = site[curr];
						site[next]->adj[6]->reverse =  site[curr]->adj[3];
					}
				} 
			}

			gndFound = false;
}

Lattice :: ~Lattice()
{
	for (int i=1; i<HIGHEST; i++)
	{
		delete active[i];
		delete inactive[i];
	}

	for (int i=0; i<SIZE; i++)
	{
		delete site[i];
		site[i] = NULL;
	}

	delete active;
	delete inactive;
	delete site;
	delete sink;
}

void Lattice :: distributeWeight(int curr, int next, Edge *adj, int initLabel[], int expandLabel, double **nbCost, bool bidirection) {
	if (adj == NULL || (initLabel[curr] == expandLabel && initLabel[next] == expandLabel)) return;

	double A = nbCost[initLabel[curr]][initLabel[next]]; 
	double B = nbCost[initLabel[curr]][expandLabel];
	double C = nbCost[expandLabel][initLabel[next]];
	double D = nbCost[expandLabel][expandLabel];

	adj->residual = adj->capacity = B + C - A - D;

	if (C > A) {
		site[curr]->excess += C - A;
	} else if (C < A) {
		site[curr]->adj[MAX_EDGE-1]->capacity += A - C;
	}

	if (C > D) {
		site[next]->adj[MAX_EDGE-1]->capacity += C - D;
	} else if (C < D) {
		site[next]->excess += D - C;
	}

	site[curr]->adj[MAX_EDGE-1]->residual = site[curr]->adj[MAX_EDGE-1]->capacity;
	site[next]->adj[MAX_EDGE-1]->residual = site[next]->adj[MAX_EDGE-1]->capacity;

	if (bidirection) {
		distributeWeight(next, curr, adj->reverse, initLabel, expandLabel, nbCost, false);
	}
}

void Lattice :: Init(int initLabel[], int expandLabel, double **nbCost, double **priorCost)
	// nbCost is a n by n matrix describing cost between labels (n is the total number of labels)
	// priorCost is a c by n matrix describing prior cost
{
	for (int i=1; i<HIGHEST; i++)
	{
		active[i]->Empty();
		inactive[i]->Empty();
	}

	gndFound = false;
	sink->excess = 0;

	int curr = 0;
	for (int k=0;k<H;k++)
		for (int j=0;j<L;j++)
			for (int i=0;i<W;i++) {
				site[curr]->excess = 0;
				site[curr]->height = 0;
				site[curr]->adj[0]->capacity = 0;
				site[curr]->adj[0]->residual = 0;
				site[curr]->adj[MAX_EDGE-1]->capacity = 0;
				site[curr]->adj[MAX_EDGE-1]->residual = 0;

				double temp = priorCost[initLabel[curr]][curr] - priorCost[expandLabel][curr];
				if (temp<0.) {
					site[curr]->excess = -temp;
				} else if (temp>0.)	{
					site[curr]->adj[MAX_EDGE-1]->residual = site[curr]->adj[MAX_EDGE-1]->capacity = temp;
				}
				curr++;
			}

	curr = 0;
	for (int k=0;k<H;k++)
		for (int j=0;j<L;j++)
			for (int i=0;i<W;i++) {
				
				int next = i < W-1 ? curr + 1 : curr + 1 - W;
				distributeWeight(curr, next, site[curr]->adj[1], initLabel, expandLabel, nbCost, true);
				
				next = j < L-1 ? curr + W : curr - W * (L-1);
				distributeWeight(curr, next, site[curr]->adj[2], initLabel, expandLabel, nbCost, true);

				next = k < H-1 ? curr + W * L : curr - W * L * (H-1);
				distributeWeight(curr, next, site[curr]->adj[3], initLabel, expandLabel, nbCost, true);

				curr++;
			}
}


bool Lattice::GndState(double frequency)
{
	if (!gndFound)
	{
		/*       nPush = 0;
		nLift = 0;
		nDischarge = 0;
		nGap = 0;
		nUpdate = 0;
		nGapNodes = 0;
		*/

		GlobalUpdate();

		while (maxActive > 0) {
			if (active[maxActive]->head == NULL) {
				maxActive--;
			}
			else {
				Vertex *u = active[maxActive]->head;
				Discharge(u);
			}

			if (workSinceUpdate * frequency > updateLimit) {
				GlobalUpdate();
			}
		}

		//     find cut
		QueuePtr<Vertex>* Qlist = new QueuePtr<Vertex>();

		for (int i=0; i<SIZE; i++) {
			Edge *curr = site[i]->adj[MAX_EDGE-1];
			if (curr->residual > eps) {
				site[i]->height = 1;
				Qlist->Enqueue(site[i]);
			} else {
				site[i]->height = HIGHEST;
			}
		}

		Vertex *curr;

		while (Qlist->Dequeue(curr))
		{

			int newHeight = curr->height + 1;

			for (int i=1; i<MAX_EDGE-1; i++) {
				if (curr->adj[i] == NULL) continue;
				Edge *rev = curr->adj[i]->reverse;
				Vertex *next = curr->adj[i]->end;

				if (rev->residual > eps && next->height == HIGHEST)
				{

					next->height = newHeight;
					Qlist->Enqueue(next);
				}
			}
		}
		delete Qlist;
		gndFound = true;
	}
	return gndFound;
}

void Lattice::UpdateLabel(int *label, int expandLabel) {
       for (int i=0; i<SIZE; i++) {
		   if (site[i]->height < HIGHEST) label[i] = expandLabel;
       }
}

void Lattice::GlobalUpdate()
{
	//     nUpdate++;

	QueuePtr<Vertex>* Qlist = new QueuePtr<Vertex>();

	workSinceUpdate = 0;
	maxActive = 0;
	maxHeight = 0;

	for (int i=1; i<HIGHEST; i++) {
		inactive[i]->Empty();
		active[i]->Empty();
	}

	for (int i=0; i<SIZE; i++) {
		Edge *curr = site[i]->adj[MAX_EDGE-1];
		if (curr->residual > eps) {
			site[i]->height = 1;
			Qlist->Enqueue(site[i]);
			if (site[i]->excess > eps) {
				active[1]->InsertFront(site[i]);
				maxActive = 1;
			}
			else {
				inactive[1]->InsertFront(site[i]);
			}
		} else {
			site[i]->height = HIGHEST;
		}
	}

	maxHeight = maxActive;

	Vertex *curr;

	while (Qlist->Dequeue(curr))
	{

		int newHeight = curr->height + 1;

		for (int i=1; i<MAX_EDGE-1; i++) {
			if (curr->adj[i] == NULL) continue;

			Edge *rev = curr->adj[i]->reverse;
			Vertex *next = curr->adj[i]->end;

			if (rev->residual > eps && next->height == HIGHEST)
			{

				next->height = newHeight;

				if (newHeight > maxHeight)
				{
					maxHeight = newHeight;
				}

				//maintenance of (in)active nodes
				if (next->excess > eps)
				{
					if (newHeight > maxActive)
					{
						maxActive = next->height;
					}

					active[newHeight]->InsertFront(next);
				}
				else
				{
					inactive[newHeight]->InsertFront(next);
				}
				//**********

				Qlist->Enqueue(next);
			}
		}
	}
	delete Qlist;
}

void Lattice::Discharge(Vertex *ver)
{
	//     nDischarge ++;

	active[ver->height]->Remove(ver);

	int start = 0;

	while (ver->excess > eps && ver->height < HIGHEST)
	{
		int admHeight = ver->height - 1;

		for (int i=start; i<MAX_EDGE; i++) {
			Edge *curr = ver->adj[i];
			if (curr == NULL) continue;
			Vertex *end = curr->end;

			if (end != NULL)
			{
				if (curr->residual > eps && end->height == admHeight)
				{

					//		  nPush++;
					double temp = (curr->residual > ver->excess) ?
						ver->excess : curr->residual;
					curr->residual -= temp;
					if (curr->reverse != NULL)
					{
						curr->reverse->residual += temp;
					} 

					if (end->excess < eps && end != sink)
					{
						inactive[admHeight]->Remove(end);
						active[admHeight]->InsertFront(end);
						if (maxActive < admHeight)
						{
							maxActive = admHeight;
						}
					}

					ver->excess -= temp;
					end->excess += temp;

					if (ver->excess < eps)
					{
						break;
					}
				}
			}
		}

		if (ver->excess < eps)
		{
			break;
		}

		int oldHeight = ver->height;
		if (inactive[oldHeight]->head == NULL && active[oldHeight]->head == NULL)
		{
			//               nGap++;
			for (int i=oldHeight+1; i<=maxHeight; i++)
			{
				//assert(active[i]->head == NULL);//because use HIGHEST lable

				while (inactive[i]->head != NULL)
				{
					//                       nGapNodes++;
					inactive[i]->head->height = HIGHEST;
					inactive[i]->Remove(inactive[i]->head);
				}
			}
			maxHeight = maxActive = oldHeight-1;
			ver->height = HIGHEST;
		}
		else
		{
			//               nLift++;
			workSinceUpdate += BETA;

			int minHeight = HIGHEST;

			for (int i=1; i<MAX_EDGE-1; i++) {
				Edge* curr = ver->adj[i];
				if (curr == NULL) continue;
				workSinceUpdate++;
				if ( curr->end->height < minHeight && curr->residual > eps )
				{
					minHeight = curr->end->height;
					//start = i - 1;
					//start++;      //have to do this stupid things to keep it from ignoring by optimizition
				}
			}

			ver->height = 1 + minHeight;
		}
	}

	if (ver->height < HIGHEST)
	{
		inactive[ver->height]->InsertFront(ver);
		if (ver->height > maxHeight)
		{
			maxHeight = ver->height;
		}
	}
}

void Lattice::IsingInit(double *mField) {
	double cost1[2] = {0.0, 1.0};
	double cost2[2] = {1.0, 0.0};
	double *nbCost[] = {cost1, cost2};

	int *initLabel = new int[SIZE];
	for (int i=0; i<SIZE; i++) initLabel[i] = 0;

	double *pField = new double[SIZE];
	for (int i=0; i<SIZE; i++) pField[i] = -mField[i];
	double **prior = new double*[2];
	prior[0] = mField;
	prior[1] = pField;

	Init(initLabel, 1, nbCost, prior);
	delete [] initLabel;
	delete [] pField;
	delete [] prior;
}

double Lattice::IsingEnergy(double *prior) {
       double energy = 0;
       for (int i=0; i<SIZE; i++) {
		short int spin = (site[i]->height < HIGHEST) ? 1: -1;
		short int nbSpin1 = (site[i]->adj[1]->end->height < HIGHEST) ? 1: -1;
		short int nbSpin2 = (site[i]->adj[2]->end->height < HIGHEST) ? 1: -1;
		energy += spin * (nbSpin1 + nbSpin2);
		if (site[i]->adj[3] != NULL) {
			short int nbSpin3 = (site[i]->adj[3]->end->height < HIGHEST) ? 1: -1;
			energy +=  spin * nbSpin3;
		}
		energy += spin * prior[i];
       }
	   return -energy;
}


double Lattice::IsingMagnetization(double *prior) {
       double M = 0;
       for (int i=0; i<SIZE; i++) {
		   short int spin = (site[i]->height < HIGHEST) ? 1: -1;
		   M += spin;
       }
	   return M;
}



