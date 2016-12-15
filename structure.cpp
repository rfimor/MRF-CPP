#include "stdafx.h"
#include <math.h>
#include "structure.h"

State :: State(double c, double r,double m)
{
       couple = c;
       rndmag = r;
       mag = m;
}

State :: State(const State & p)
{
       couple = p.couple;
       rndmag = p.rndmag;
       mag = p.mag;
}

bool State :: operator == (const State & p)
{
      return (fabs(couple - p.couple) < 1e-9 &&
              fabs(rndmag - p.rndmag) < 1e-9 &&
              fabs(mag - p.mag) < 1e-9);
}

void State :: operator = (const State & p)
{
           if ( !(this == &p) ) {
              couple = p.couple;
              rndmag = p.rndmag;
              mag = p.mag;
           }
}

double State :: Energy(const Point2D *test)
{
       return -(couple + rndmag * test->delta + mag * test->hfield);
}

/*
bool State :: OnState(const Point3D *p,double &diff,const double EPS)
{
       diff = fabs((couple + rndmag * p->delta + mag * p->hfield + p->energy)
                    /p->energy);
       return (diff < EPS);
}

bool State :: OnState(const Point3D *p,const double EPS)
{
       double diff = fabs((couple + rndmag * p->delta + mag * p->hfield + p->energy)
                    /p->energy);
       return (diff < EPS);
}
*/

Point2D :: Point2D()
{}

Point2D :: Point2D(double h,double d)
{
       delta = d;
       hfield = h;
}

bool Point2D :: operator == (const Point2D & p)
{
       return (fabs(delta - p.delta) < 1e-9 &&
               fabs(hfield - p.hfield) < 1e-9 );
}

void Point2D :: operator = (const Point2D & p)
{
           if ( !(this == &p) ) {
              delta = p.delta;
              hfield = p.hfield;
           }
}

void Point2D :: operator = (const Point3D & p)
{
            delta = p.delta;
            hfield = p.hfield;
}

Point3D :: Point3D()
{}

Point3D :: Point3D(double h, double d,double e)
{
       delta = d;
       hfield = h;
       energy = e;
}

bool Point3D :: operator == (const Point3D & p)
{
       return (fabs(delta - p.delta) < 1e-9 &&
               fabs(hfield - p.hfield) < 1e-9 &&
               fabs(energy -p.energy) < 1e-9 );
}

void Point3D :: operator = (const Point3D & p)
{
           if ( !(this == &p) ) {
              delta = p.delta;
              hfield = p.hfield;
              energy = p.energy;
           }
}

