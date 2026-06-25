/*
This file is part of Agave software.

    Agave is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    Agave is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Agave.  If not, see <http://www.gnu.org/licenses/>.

Cite as follows:
 
R. Flores-Moreno, Agave 1, Guadalajara Jal. (2025)
*/

#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <iomanip>
#include <fstream>

#include <Global.h>
//#include <Integrator.h>
#include <Math.h>
#include <Parameter.h>
#include <Shell.h>
#include <Vector.h>

using namespace std;

Shell::Shell(void)
{
  z.clear();
  d.clear();

  is_normalized = false;
  radius = 1.0/TOLNUM;
};

void Shell::Normalize(void)
{
  if (is_normalized) return;

  int i,j,n;
  double factor,norm;

  n = z.size();

  // Normalized GTOs
  norm = pow(2.0,double(l))*pow(2.0/M_PI,0.75);
  for ( i = 0 ; i < n ; i++ )
  {
    d[i] = norm*pow(z[i],(2.0*l+3.0)*0.25)*d[i];
  }

  factor = pow(M_PI,1.5)/pow(2.0,(double)l);
  norm = 0.0;
  for ( i = 0 ; i < n ; i++ )
  {
    for ( j = 0 ; j < n ; j++ )
    {
      norm += d[i]*d[j]*factor/pow(z[i]+z[j],l+1.5);
    }
  }

  // Normalized STOs
  int lx,ly,lz;
  nco = 0;
  for ( lx = l ; lx >= 0 ; lx-- )
  {
    for ( ly = l - lx ; ly >= 0 ; ly-- )
    {
      lz = l - lx - ly;
      ncsto.push_back( 1.0/sqrt(norm*DoubleFactorial(2*lx-1)*
                   DoubleFactorial(2*ly-1)*DoubleFactorial(2*lz-1)) );
      nco++;
    }
  }

  is_normalized = true;
}

void Shell::WriteFile(const char* filename)
{
  ofstream f(filename,ios::app);

  double norm = pow(2.0,double(l))*pow(2.0/M_PI,0.75);
  int k = z.size();
  f << setw(2) << fixed << n << " "
    << setw(2) << fixed << l << " "
    << setw(3) << fixed << k <<endl;
  for ( int i = 0 ; i < k ; i++ )
  {
     f << setw(20) << fixed << setprecision(10) << z[i] << "     ";
     f << setw(20) << fixed << setprecision(10)<< d[i]/(norm*pow(z[i],(2.0*l+3.0)*0.25))<<endl;
  }
  f.close();
}

void Shell::PrintOut()
{
  double norm = pow(2.0,double(l))*pow(2.0/M_PI,0.75);
  int k = z.size();
  cout << setw(2) << fixed << n << " "
       << setw(2) << fixed << l << " "
       << setw(3) << fixed << k <<endl;
  for ( int i = 0 ; i < k ; i++ )
  {
     cout << setw(20) << fixed << setprecision(10) << z[i] << "     ";
     cout << setw(20) << fixed << setprecision(10)<< d[i]/(norm*pow(z[i],(2.0*l+3.0)*0.25))<<endl;
  }
}


void Shell::EvaluateGTO(double *xc, double *yc, double *zc, double* res,int indx, int nr)
{
  Vector a;
  int i,ig,lx,ly,lz,n;
  double expf,r2;
  double xyzp;

  n = z.size();
  for ( ig = 0 ; ig < nr ; ig++ )
  {
    a = Vector(xc[ig],yc[ig],zc[ig]);
    r2 = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    if (sqrt(r2)<radius)
    {
      expf = 0.0;
      for ( i = 0 ; i < n ; i++ )
      {
        expf += d[i]*exp(-z[i]*r2);
      }
      int ico = 0;
      for ( lx = l ; lx >= 0 ; lx-- )
      {
        for ( ly = l - lx ; ly >= 0 ; ly-- )
        {
          lz = l - lx - ly;
          if (ico==indx)
          {
            xyzp = 1.0;
            for (i=1;i<=lx;i++) xyzp *= a.x;
            for (i=1;i<=ly;i++) xyzp *= a.y;
            for (i=1;i<=lz;i++) xyzp *= a.z;
            res[ig] = ncsto[ico]*xyzp*expf;
          }
          ico++;
        }
      }
    } else res[ig] = 0.0;
  }
};


void Shell::EvaluateGTODerivative(double *xc, double *yc, double *zc, double* dx, double* dy, double *dz,int indx, int nr)
{
  Vector a;
  int i,ig,lx,ly,lz,n;
  double expf,expfp,r2,t;
  double rxp[MAXL+2],ryp[MAXL+2],rzp[MAXL+2];

  n = z.size();
  for ( ig = 0 ; ig < nr ; ig++ )
  {
    a = Vector(xc[ig],yc[ig],zc[ig]);
    r2 = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    expf = 0.0;
    expfp = 0.0;
    for ( i = 0 ; i < n ; i++ )
    {
      t = d[i]*exp(-z[i]*r2);
      expf += t;
      expfp += t*2.0*z[i];
    }
    rxp[0] = 1.0;
    ryp[0] = 1.0;
    rzp[0] = 1.0;
    for (lx=1;lx<=l+1;lx++) rxp[lx] = rxp[lx-1]*a.x;
    for (ly=1;ly<=l+1;ly++) ryp[ly] = ryp[ly-1]*a.y;
    for (lz=1;lz<=l+1;lz++) rzp[lz] = rzp[lz-1]*a.z;
    int ico = 0;
    for ( lx = l ; lx >= 0 ; lx-- )
    {
      for ( ly = l - lx ; ly >= 0 ; ly-- )
      {
        lz = l - lx - ly;
        if (ico==indx)
        {
          t = ncsto[ico]*ryp[ly]*rzp[lz];
          dx[ig] = t*expfp*rxp[lx+1];
          if (lx>0) dx[ig] -= double(lx)*t*expf*rxp[lx-1];
          t = ncsto[ico]*rxp[lx]*rzp[lz];
          dy[ig] = t*expfp*ryp[ly+1];
          if (ly>0) dy[ig] -= double(ly)*t*expf*ryp[ly-1];
          t = ncsto[ico]*ryp[ly]*rxp[lx];
          dz[ig] = t*expfp*rzp[lz+1];
          if (lz>0) dz[ig] -= double(lz)*t*expf*rzp[lz-1];
        }
        ico++;
      }
    }
  }
};

/*
double Shell::Radius()
{
  return radius;
}

void Shell::EvaluateRadius()
{
  size_t i,imin;
  double delta,dg,g,ldt,zr;

  // Get smallest exponent
  imin = 0;
  for (i=1;i<z.size();i++)
  {
    if (z[i]<z[imin])
      imin = i;
  }

  ldt = log(ABS(d[imin])/ABS(1.0e-10));
  radius = MAX(ldt/z[imin],TOLNUM);

  // For large radial powers the guess should be corrected ***
  if (l>0)
  {
    g = MAX(sqrt(double(l)/(2.0*ABS(d[imin]))),TOLNUM);
    if (g>radius) radius = 0.5*(radius+g);
  }

  for (i=0;i<100;i++)
  {
    zr = z[imin]*radius;
    g = ldt + double(l)*log(radius) - zr*radius;
    dg = double(l)/radius - 2.0*zr;
    delta = g/dg;
    radius = MAX(radius-delta,TOLNUM);
    if (ABS(delta)<1.0E-10) break;
  }
  if (ABS(delta)>=1.0E-10) 
  {
    radius = 1.0/TOLNUM;
    cout << "Unable to find shell radius, infinity assumed" << endl;
  }
}

void Shell::NormalizeInteraction(Integrator *gi)
{
  if (is_normalized) return;

  int dim,i,lx,ly;
  double ra[3] = {0.0,0.0,0.0};
  double blk[MAXNCO][MAXNCO];

  nco = ((l+1)*(l+2))/2;
  ncsto.clear();
  for ( i = 0 ; i < nco ; i++ )
    ncsto.push_back( 1.0 );
  
  gi->Interaction2(ra,ra,this,this,blk,0);

  i = 0;
  for ( lx = l ; lx >= 0 ; lx-- )
  {
    for ( ly = l - lx ; ly >= 0 ; ly-- )
    {
      ncsto[i] = 1.0/sqrt(blk[i][i]);
      i++;
    }
  }

  is_normalized = true;
}

// Compute Radial potential (only exponential part)
double Shell::RadialPotential(double r)
{
  size_t i;
  double expf,r2;

  // Build potential
  r2 = r*r;

  // Contract the potential of the current shell 
  expf = 0.0;
  for (size_t i=0;i<z.size();i++)
    expf += d[i]*exp(-z[i]*r2);

  return expf;
}

*/
