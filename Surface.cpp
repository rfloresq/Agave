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

#include <iostream>

#include <QtWidgets>

#include <Agave.h>
#include <Atom.h>
#include <Basis.h>
#include <Evaluator.h>
#include <Global.h>
#include <Grid.h>
#include <Surface.h>
#include <System.h>
#include <Species.h>
#include <Structure.h>
#include <Set.h>
#include <Shell.h>
#include <Plotter.h>
#include <Math.h>

using namespace std;

Surface::Surface(System* s,Grid *g)
{
  sys = s;
  grid = g;
  scalar = new double[grid->npoint[0]];
};

void Surface::GetCubeValues(int i,int j,int k, double* val)
{
  int l,ll;
  ll = grid->npoint[1]*grid->npoint[2]*(k-1);

  l = ll + grid->npoint[1]*j + i-1;    
  val[0] = scalar[l];
  val[1] = scalar[l+1];
  l = ll + grid->npoint[1]*(j-1) + i;
  val[2] = scalar[l];
  val[3] = scalar[l-1];
  ll = grid->npoint[1]*grid->npoint[2]*k;
  l = ll + grid->npoint[1]*j + i-1;
  val[4] = scalar[l];
  val[5] = scalar[l+1];
  l = ll + grid->npoint[1]*(j-1) + i;
  val[6] = scalar[l];
  val[7] = scalar[l-1];
}

void Surface::Setup(void)
{
}

void Surface::Build(void)
{

  if (type==COLOURED_ISO)
  {
    //BuildParticleDensity(ELECTRON,1.0,false);
  }
  else
  {
    if ( strncasecmp(name,"ORBITAL",7) == 0 )  
      BuildOrbital();
    //else if ( strncasecmp(name,"CHARGE DENSITY",14) == 0 ) 
    //  BuildChargeDensity(false);
    else if ( strncasecmp(name,"PARTICLE DENSITY",16) == 0 ) 
      BuildParticleDensity(species,1.0,false);
    //else if ( strncasecmp(name,"SPIN DENSITY",16) == 0 ) 
    //  BuildParticleDensity(species,-1.0,false);
    //else if ( strncasecmp(name,"RESPONSE",8) == 0 ) 
    //  BuildParticleDensity(species,1.0,true);
    //else if ( strncasecmp(name,"FUKUI",5) == 0 ) 
    //  BuildFukui();
    else if ( strncasecmp(name,"GTO",3) == 0 ) 
      BuildGTO();
  }
};

void Surface::BuildGTO()
{
  int i,ib,is,io,j,k,ni,nj,nk;

  ni = grid->npoint[1];
  nj = grid->npoint[2];
  nk = grid->npoint[3];
  double x[ni];
  double y[ni];
  double z[ni];

  nscalar = 0;
  i = 0;
  for (ib=0;ib<(signed)sys->species[species]->basis->set.size();ib++)
  {
    for (is=0;is<(signed)sys->species[species]->basis->set[ib]->shell.size();is++)
    {
      Shell *shell = sys->species[species]->basis->set[ib]->shell[is];
      for (io=0;io<shell->nco;io++)
      {
        if (i==orbital_number) 
        {
          double x0,y0,z0;

          x0 = sys->mol->atom[ib]->x;
          y0 = sys->mol->atom[ib]->y;
          z0 = sys->mol->atom[ib]->z;
          for ( k = 0 ; k < nk ; k++ )
          {
            for ( j = 0 ; j < nj ; j++ )
            {
              grid->GetRowOfPoints(1,j,k,x,y,z);
              for (int ii=0;ii<ni;ii++) 
              {
                x[ii] -= x0;
                y[ii] -= y0;
                z[ii] -= z0;
              }
              shell->EvaluateGTO(x,y,z,&scalar[nscalar],io,ni);
              nscalar += ni;
            }
          }
          return;
        }
        i++;
      }
    }
  }
};

void Surface::BuildOrbital()
{
  if (species>=sys->species.size()) 
  {
    cout << "WARNING: Attempt to plot for non existing species"<<endl;
    return;
  } 
  if (spin>=sys->species[species]->nspin) 
  {
    cout << "WARNING: Attempt to plot for non existing spin"<<endl;
    return;
  }
  if (orbital_number>=sys->species[species]->nmo[spin])
  {
    cout << "WARNING: Attempt to non existing orbital"<<endl;
    cout << orbital_number << " from "
         << sys->species[species]->nmo[spin]
	 << endl;
    return;
  }
  
  int i,io,j,k,ni,nj,nk;
  int nco = sys->species[species]->basis->nco;
  double gto[nco];
  double moc[nco];
  sys->species[species]->C[spin]->GetColValues(orbital_number,moc);

  ni = grid->npoint[1];
  nj = grid->npoint[2];
  nk = grid->npoint[3];
  double x[ni];
  double y[ni];
  double z[ni];

  nscalar = 0;
  for ( k = 0 ; k < nk ; k++ )
  {
    for ( j = 0 ; j < nj ; j++ )
    {
      grid->GetRowOfPoints(1,j,k,x,y,z);
      for (i=0;i<ni;i++)
      {
        GetBasis(sys->species[species]->basis,x[i],y[i],z[i],gto); 
        scalar[nscalar] = 0.0;
        for (io=0;io<nco;io++)
          scalar[nscalar] += moc[io]*gto[io];
        nscalar++;
      }
    }
  }
};



void Surface::GetDerivatives(double* x, double *y, double *z, double *dx,
double *dy, double *dz, int n)
{
  if (type==COLOURED_ISO)
  {
    // GetParticleDensityDerivatives(sys->species[ELECTRON],x,y,z,dx,dy,dz,n,1.0,false);
  }
  else
  {
    if ( strncasecmp(name,"ORBITAL",7) == 0 ) 
      GetOrbitalDerivatives(sys->species[species],spin,
		            orbital_number,x,y,z,dx,dy,dz,n);
	  /*
    else if ( strncasecmp(name,"CHARGE DENSITY",14) == 0 ) 
      GetChargeDensityDerivatives(sys,x,y,z,dx,dy,dz,n,false);
          */
    else if ( strncasecmp(name,"PARTICLE DENSITY",16) == 0 ) 
      GetParticleDensityDerivatives(sys->species[species],
		                    x,y,z,dx,dy,dz,n,1.0,false);
    /*
    else if ( strncasecmp(name,"SPIN DENSITY",12) == 0 ) 
      GetParticleDensityDerivatives(sys->species[species],
		                x,y,z,dx,dy,dz,n,-1.0,false);
    else if ( strncasecmp(name,"RESPONSE",8) == 0 ) 
      GetParticleDensityDerivatives(sys->species[species],
		                x,y,z,dx,dy,dz,n,1.0,true);
    //else if ( strncasecmp(name,"FUKUI",5) == 0 ) 
    //  GetFukuiDerivatives(sys,x,y,z,dx,dy,dz,n,name);
    */ 
    else if ( strncasecmp(name,"GTO",3) == 0 ) 
      GetGTODerivatives(sys->species[species]->basis,orbital_number,x,y,z,dx,dy,dz,n);
  }
};


void Surface::GetValues(double* x, double *y, double *z, double *v, int n)
{
	/*
  if ( strncasecmp(name,"ORBITAL",7) == 0 ) 
    GetOrbital(sys->species[species]->spin[spin],
	       sys->species[species]->basis,orbital_number,x,y,z,v,n);
  else if ( strncasecmp(name,"CHARGE DENSITY",14) == 0 ) 
    GetChargeDensity(sys,x,y,z,v,n,false);
  else if ( strncasecmp(name,"PARTICLE DENSITY",16) == 0 ) 
    GetParticleDensity(sys->species[species],x,y,z,v,n,1.0,false);
  else if ( strncasecmp(name,"SPIN DENSITY",12) == 0 ) 
    GetParticleDensity(sys->species[species],x,y,z,v,n,-1.0,false);
  else if ( strncasecmp(name,"RESPONSE",8) == 0 ) 
    GetParticleDensity(sys->species[species],x,y,z,v,n,1.0,true);
  //else if ( strncasecmp(name,"FUKUI",5) == 0 ) 
  //  GetFukui(sys,x,y,z,v,n,name);
  else if ( strncasecmp(name,"ASC",3) == 0 ) 
    GetASC(sys,x,y,z,v,n);
  else if ( strncasecmp(name,"GTO",3) == 0 ) 
  if ( strncasecmp(name,"GTO",3) == 0 ) 
    GetGTO(sys->species[species]->basis,orbital_number,x,y,z,v,n);
    */
};

/*
void Surface::GetExtrema(double *vmin, double *vmax)
{
  int i;
  double val;

  *vmin = 99999999.0;
  *vmax = 0.0;
  for (i=0;i<nscalar;i++)
  {
    val = ABS(scalar[i]);
    *vmin = MIN(*vmin,val);
    *vmax = MAX(*vmax,val);
  }
}
*/

// Functions
// Density like functions
void Surface::BuildDensityLike(Matrix *P, Matrix *Q, Basis* basis,const char* option)
{
  int i,io,j,jo,k,ni,nj,nk;

  int nco = basis->nco;
  double gto[nco];
  double row[nco];
  double row2[nco];

  ni = grid->npoint[1];
  nj = grid->npoint[2];
  nk = grid->npoint[3];
  double x[ni];
  double y[ni];
  double z[ni];

  nscalar = 0;
  for ( k = 0 ; k < nk ; k++ )
  {
    for ( j = 0 ; j < nj ; j++ )
    {
      grid->GetRowOfPoints(1,j,k,x,y,z);
      for (i=0;i<ni;i++)
      {
        GetBasis(basis,x[i],y[i],z[i],gto); 
        scalar[nscalar] = 0.0;
        for (io=0;io<nco;io++)
        {
          P->GetRowValues( io , row );
	  if ( strncasecmp(option,"Average",7) == 0 ) 
          {
            Q->GetRowValues( io , row2 );
            for (jo=0;jo<nco;jo++)
	    {
              row[jo] += row2[jo];
              row[jo] *= 0.5;
	    }
	  }
	  else if ( strncasecmp(option,"Substract",9) == 0 ) 
          {
            Q->GetRowValues( io , row2 );
            for (jo=0;jo<nco;jo++)
              row[jo] -= row2[jo];
	  }
          for (jo=0;jo<nco;jo++)
            scalar[nscalar] += row[jo]*gto[io]*gto[jo];
        }
        nscalar++;
      }
    }
  }
};


/*
// Charge density
void Surface::BuildChargeDensity(bool response)
{
  int i,ispec,ispin,nbuf;
  
  nbuf = grid->npoint[0];
  double buf[nbuf];

  for (i=0;i<nbuf;i++)
    buf[i] = 0.0;

  for (ispec=0;ispec<sys->species.size();ispec++)
  {
    double charge = sys->species[ispec]->Charge();
    for (ispin=0;ispin<sys->species[ispec]->spin.size();ispin++)
    {
      if (response)
        BuildDensityLike(sys->species[ispec]->spin[ispin]->P1,
                         sys->species[ispec]->spin[ispin]->P1,
	                 sys->species[ispec]->basis,"None");
      else
        BuildDensityLike(sys->species[ispec]->spin[ispin]->P,
                         sys->species[ispec]->spin[ispin]->P,
	  	         sys->species[ispec]->basis,"None");
      for (i=0;i<nbuf;i++)
        buf[i] += scalar[i]*charge;
    }
  }

  for (i=0;i<nbuf;i++)
    scalar[i] = buf[i];
}
*/

// Particle density
void Surface::BuildParticleDensity(int species,double sf, bool response)
{
  int i,ispin,nbuf;
  
  nbuf = grid->npoint[0];
  double buf[nbuf];

  for (i=0;i<nbuf;i++)
    buf[i] = 0.0;

  for (ispin=0;ispin<sys->species[species]->nspin;ispin++)
  {
	  /*
    if (response)
      BuildDensityLike(sys->species[species]->spin[ispin]->P1,
                       sys->species[species]->spin[ispin]->P1,
                       sys->species[species]->basis,"None");
    else
    */
      BuildDensityLike(sys->species[species]->P[ispin],
                       sys->species[species]->P[ispin],
  	               sys->species[species]->basis,"None");
    if (spin==ALPHA)
      for (i=0;i<nbuf;i++)
        buf[i] += scalar[i];
    else if (spin==BETA)
      for (i=0;i<nbuf;i++)
        buf[i] += sf*scalar[i];
  }

  for (i=0;i<nbuf;i++)
    scalar[i] = buf[i];
}

/*
// Electronic Fukui functions
void Surface::BuildFukui()
{
  int ispec;
  
  ispec = 0; // Electrons
  if ( strncasecmp(name,"FUKUI REMOVAL",13) == 0 ) 
  {
    BuildDensityLike(sys->species[ispec]->FUKUIL,
                     sys->species[ispec]->FUKUIR,
                     sys->species[ispec]->basis,"None");
  }
  else if ( strncasecmp(name,"FUKUI ADDITION",14) == 0 ) 
  {
    BuildDensityLike(sys->species[ispec]->FUKUIR,
                     sys->species[ispec]->FUKUIL,
                     sys->species[ispec]->basis,"None");
  }
  else if ( strncasecmp(name,"FUKUI AVERAGE",13) == 0 ) 
  {
    BuildDensityLike(sys->species[ispec]->FUKUIL,
                     sys->species[ispec]->FUKUIR,
                     sys->species[ispec]->basis,"Average");
  }
  else if ( strncasecmp(name,"FUKUI DUAL",10) == 0 ) 
  {
    BuildDensityLike(sys->species[ispec]->FUKUIR,
                     sys->species[ispec]->FUKUIL,
                     sys->species[ispec]->basis,"Substract");
  }
}
*/
