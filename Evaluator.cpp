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

#include <Atom.h>
#include <Basis.h>
#include <Evaluator.h>
#include <Global.h>
#include <Integrator.h>
#include <Math.h>
#include <Matrix.h>
#include <System.h>
#include <Structure.h>
#include <Species.h>
#include <Set.h>
#include <Shell.h>

using namespace std;

void GetGTODerivatives(Basis *basis, int id, double* x, double *y, double *z, double *dx, double *dy, double *dz, int n)
{
  int i,ib,is,io;

  double wx[n];
  double wy[n];
  double wz[n];

  i = 0;
  for (ib=0;ib<(signed)basis->set.size();ib++)
  {
    for (is=0;is<(signed)basis->set[ib]->shell.size();is++)
    {
      Shell *shell = basis->set[ib]->shell[is];
      for (io=0;io<shell->nco;io++)
      {
        if (i==id) 
        {
          double x0,y0,z0;

          x0 = basis->mol->atom[ib]->x;
          y0 = basis->mol->atom[ib]->y;
          z0 = basis->mol->atom[ib]->z;
          for (int ii=0;ii<n;ii++) 
          {
            wx[ii] = x[ii] - x0;
            wy[ii] = y[ii] - y0;
            wz[ii] = z[ii] - z0;
          }
          shell->EvaluateGTODerivative(wx,wy,wz,dx,dy,dz,io,n);
          return;
        }
        i++;
      }
    }
  }
};

void GetBasisDerivatives(Basis *basis, double x, double y, double z,
 double *dx, double *dy, double *dz)
{
  int i,ib,is,io,lla;
  Shell *shell;

  double xx[1];
  double yy[1];
  double zz[1];
  double dxx[1];
  double dyy[1];
  double dzz[1];

  for (ib=0;ib<(signed)basis->set.size();ib++)
  {
    xx[0] = x - basis->mol->atom[ib]->x;
    yy[0] = y - basis->mol->atom[ib]->y;
    zz[0] = z - basis->mol->atom[ib]->z;
    for (is=0;is<(signed)basis->set[ib]->shell.size();is++)
    {
      shell = basis->set[ib]->shell[is];
      lla = shell->ll + basis->set[ib]->ll;
      for (io=0;io<shell->nco;io++)
      {
        i = lla + io;
        shell->EvaluateGTODerivative(xx,yy,zz,dxx,dyy,dzz,io,1);
        dx[i] = dxx[0];
        dy[i] = dyy[0];
        dz[i] = dzz[0];
      }
    }
  }
};


void GetOrbitalDerivatives(Species* spec, int spin,int id, double* x, double *y, double *z, double *dx, double *dy, double *dz, int n)
{
  int i,io;

  if (!spec->C[spin]) return; // No result for null MO matrix

  int nco = spec->basis->nco;
  double moc[nco];
  spec->C[spin]->GetColValues( id , moc );

  double gto_dx[nco];
  double gto_dy[nco];
  double gto_dz[nco];

  for (i=0;i<n;i++) 
  {
    GetBasisDerivatives(spec->basis,x[i],y[i],z[i],gto_dx,gto_dy,gto_dz);
    dx[i] = 0.0;
    dy[i] = 0.0;
    dz[i] = 0.0;
    for (io=0;io<nco;io++) 
    {
      dx[i] += gto_dx[io]*moc[io];
      dy[i] += gto_dy[io]*moc[io];
      dz[i] += gto_dz[io]*moc[io];
    }
  }
};

/*


void GetChargeDensityDerivatives(System* sys,double* x, double *y, double *z, double *dx, double *dy, double *dz, int n, bool response)
{
  int i,ispec,ispin;

  double dxbuf[n],dybuf[n],dzbuf[n];
  for (i=0;i<n;i++) 
  {
    dx[i] = 0.0;
    dy[i] = 0.0;
    dz[i] = 0.0;
  }
  for (ispec=0;ispec<sys->species.size();ispec++)
  {
    double charge = sys->species[ispec]->Charge();

    for (ispin=0;ispin<sys->species[ispec]->spin.size();ispin++)
    {
      if (response) 
	GetDensityLikeDerivatives(sys->species[ispec]->spin[ispin]->P1,
	                          sys->species[ispec]->spin[ispin]->P1,
                                  sys->species[ispec]->basis,
				  x,y,z,dxbuf,dybuf,dzbuf,n,"None");
      else 
	GetDensityLikeDerivatives(sys->species[ispec]->spin[ispin]->P,
	                          sys->species[ispec]->spin[ispin]->P,
		                  sys->species[ispec]->basis,
				  x,y,z,dxbuf,dybuf,dzbuf,n,"None");
      for (i=0;i<n;i++) 
      {
        dx[i] += dxbuf[i]*charge;
        dy[i] += dybuf[i]*charge;
        dz[i] += dzbuf[i]*charge;
      }
    }
  }
}
*/

void GetParticleDensityDerivatives(Species* species,double* x, double *y, double *z, double *dx, double *dy, double *dz, int n, double sf, bool response)
{
  int i,ispin;

  double dxbuf[n],dybuf[n],dzbuf[n];
  for (i=0;i<n;i++) 
  {
    dx[i] = 0.0;
    dy[i] = 0.0;
    dz[i] = 0.0;
  }
  for (ispin=0;ispin<species->nspin;ispin++)
  {
	  /*
    if (response) 
      GetDensityLikeDerivatives(species->spin[ispin]->P1,
	                        species->spin[ispin]->P1,
                                species->basis,
				x,y,z,dxbuf,dybuf,dzbuf,n,"None");
    else 
    */
      GetDensityLikeDerivatives(species->P[ispin],species->P[ispin],
		                species->basis,
				x,y,z,dxbuf,dybuf,dzbuf,n,"None");
    if (ispin==ALPHA)
      for (i=0;i<n;i++) 
      {
        dx[i] += dxbuf[i];
        dy[i] += dybuf[i];
        dz[i] += dzbuf[i];
      }
    else if (ispin==BETA)
      for (i=0;i<n;i++) 
      {
        dx[i] += sf*dxbuf[i];
        dy[i] += sf*dybuf[i];
        dz[i] += sf*dzbuf[i];
      }
  }
}

/*
void GetFukuiDerivatives(System* sys,double* x, double *y, double *z, double *dx, double *dy, double *dz, int n,char* name)
{
  int ispec;
  
  ispec = 0; // Electrons
  if ( strncasecmp(name,"FUKUI REMOVAL",13) == 0 ) 
  {
    GetDensityLikeDerivatives(sys->species[ispec]->FUKUIL,
                              sys->species[ispec]->FUKUIR,
                              sys->species[ispec]->basis,
			      x,y,z,dx,dy,dz,n,"None");
  }
  else if ( strncasecmp(name,"FUKUI ADDITION",14) == 0 ) 
  {
    GetDensityLikeDerivatives(sys->species[ispec]->FUKUIR,
                              sys->species[ispec]->FUKUIL,
                              sys->species[ispec]->basis,
			      x,y,z,dx,dy,dz,n,"None");
  }
  else if ( strncasecmp(name,"FUKUI AVERAGE",13) == 0 ) 
  {
    GetDensityLikeDerivatives(sys->species[ispec]->FUKUIL,
                              sys->species[ispec]->FUKUIR,
                              sys->species[ispec]->basis,
			      x,y,z,dx,dy,dz,n,"Average");
  }
  else if ( strncasecmp(name,"FUKUI DUAL",10) == 0 ) 
  {
    GetDensityLikeDerivatives(sys->species[ispec]->FUKUIR,
                              sys->species[ispec]->FUKUIL,
                              sys->species[ispec]->basis,
			      x,y,z,dx,dy,dz,n,"Substract");
  }
}

      
void GetGTO(Basis *basis, int id, double* x, double *y, double *z, double *v, int n)
{
  int i,ib,is,io;

  double wx[n];
  double wy[n];
  double wz[n];

  i = 0;
  for (ib=0;ib<(signed)basis->set.size();ib++)
  {
    for (is=0;is<(signed)basis->set[ib]->shell.size();is++)
    {
      Shell *shell = basis->set[ib]->shell[is];
      for (io=0;io<shell->nco;io++)
      {
        if (i==id) 
        {
          double x0,y0,z0;

          x0 = basis->mol->atom[ib]->x;
          y0 = basis->mol->atom[ib]->y;
          z0 = basis->mol->atom[ib]->z;
          for (int ii=0;ii<n;ii++) 
          {
            wx[ii] = x[ii] - x0;
            wy[ii] = y[ii] - y0;
            wz[ii] = z[ii] - z0;
          }
          shell->EvaluateGTO(wx,wy,wz,v,io,n);
          return;
        }
        i++;
      }
    }
  }
};
*/

void GetBasis(Basis* basis, double x, double y, double z,double *v)
{
  int i,ib,io,is,lla;
  double xx[1];
  double yy[1];
  double zz[1];

  for (ib=0;ib<(signed)basis->set.size();ib++)
  {
    for (is=0;is<(signed)basis->set[ib]->shell.size();is++)
    {
      Shell *shell = basis->set[ib]->shell[is];
      lla = shell->ll + basis->set[ib]->ll;
      for (io=0;io<shell->nco;io++)
      {
        i = lla + io;
        if (basis->set[ib]->is_neighbor)
        {
          xx[0] = x - basis->mol->atom[ib]->x;
          yy[0] = y - basis->mol->atom[ib]->y;
          zz[0] = z - basis->mol->atom[ib]->z;
          shell->EvaluateGTO(xx,yy,zz,&v[i],io,1);
        }
        else v[i] = 0.0;
      }
    }
  }
};


/*
void GetOrbital(Spin* spin, Basis* basis, int id, double* x, double *y, double *z, double *v, int n)
{
  int i,io;

  if (!spin->C) return; // No result for null MO matrix

  int nco = basis->nco;
  double gto[nco];
  double moc[nco];
  spin->C->GetColValues( id , moc );


  for (i=0;i<n;i++) 
  {
    GetBasis(basis,x[i],y[i],z[i],gto);
    v[i] = 0.0;
    for (io=0;io<nco;io++) 
      v[i] += gto[io]*moc[io];
  }
};




void GetChargeDensity(System* sys,double* x, double *y, double *z, double *v, int n, bool response)
{
  int i,ispec,ispin;
  
  double buf[n];

  for (i=0;i<n;i++)
    v[i] = 0.0;
  for (ispec=0;ispec<sys->species.size();ispec++)
  {
    double charge = sys->species[ispec]->Charge();
    for (ispin=0;ispin<sys->species[ispec]->spin.size();ispin++)
    {
      if (response) GetDensityLike(sys->species[ispec]->spin[ispin]->P1,
                                   sys->species[ispec]->spin[ispin]->P1,
                                   sys->species[ispec]->basis,x,y,z,buf,n,"None");
      else GetDensityLike(sys->species[ispec]->spin[ispin]->P,
                          sys->species[ispec]->spin[ispin]->P,
		          sys->species[ispec]->basis,x,y,z,buf,n,"None");
      for (i=0;i<n;i++)
        v[i] += buf[i]*charge;
    }
  }
}
void GetParticleDensity(Species* species,double* x, double *y, double *z, double *v, int n, double sf, bool response)
{
  int i,ispin;
  
  double buf[n];

  for (i=0;i<n;i++)
    v[i] = 0.0;
  for (ispin=0;ispin<species->spin.size();ispin++)
  {
    if (response) GetDensityLike(species->spin[ispin]->P1,
                                 species->spin[ispin]->P1,
                                 species->basis,x,y,z,buf,n,"None");
    else GetDensityLike(species->spin[ispin]->P,
                        species->spin[ispin]->P,
		        species->basis,x,y,z,buf,n,"None");
    if (ispin==ALPHA)
      for (i=0;i<n;i++)
        v[i] += buf[i];
    else if (ispin==BETA)
      for (i=0;i<n;i++)
        v[i] += sf*buf[i];
  }
}

void GetFukui(System* sys,double* x, double *y, double *z, double *v, int n,char* name)
{
  int ispec;
  
  ispec = 0; // Electrons
  if ( strncasecmp(name,"FUKUI REMOVAL",13) == 0 ) 
  {
    GetDensityLike(sys->species[ispec]->FUKUIL,
                   sys->species[ispec]->FUKUIR,
		   sys->species[ispec]->basis,x,y,z,v,n,"None");
  }
  else if ( strncasecmp(name,"FUKUI ADDITION",14) == 0 ) 
  {
    GetDensityLike(sys->species[ispec]->FUKUIR,
                   sys->species[ispec]->FUKUIL,
		   sys->species[ispec]->basis,x,y,z,v,n,"None");
  }
  else if ( strncasecmp(name,"FUKUI AVERAGE",13) == 0 ) 
  {
    GetDensityLike(sys->species[ispec]->FUKUIL,
                   sys->species[ispec]->FUKUIR,
		   sys->species[ispec]->basis,x,y,z,v,n,"Average");
  }
  else if ( strncasecmp(name,"FUKUI DUAL",10) == 0 ) 
  {
    GetDensityLike(sys->species[ispec]->FUKUIR,
                   sys->species[ispec]->FUKUIL,
		   sys->species[ispec]->basis,x,y,z,v,n,"Substract");
  }
}
*/

/*
void GetASC(System* sys,double* x, double *y, double *z, double *v, int n)
{
  double ax,ay,az;
  double dx,dy,dz;
  double px,py,pz;
  double d2,d2min;

  cout << "N points: "<<n<<endl;
  */
    /*
  for (int ip=0;ip<n;ip++)
  {
    px = x[ip];
    py = y[ip];
    pz = z[ip];
    d2min = 1.0e+10;
    for (size_t ia=0;ia<sys->fv.size();ia++)
    {
      ax = sys->fx[ia];
      ay = sys->fy[ia];
      az = sys->fz[ia];
      dx = ax - px;
      dy = ay - py;
      dz = az - pz;
      d2 = dx*dx + dy*dy + dz*dz;
      if (d2<d2min)
      {
	d2min = d2;
	v[ip] = sys->fv[ia];
      }
    }
  }
    */
/*
}
*/

/*
// Functions

// ==========================
// Density like functions
void GetDensityLike(Matrix *P, Matrix *Q, Basis* basis,
double* x, double *y, double *z, double *v, int n,const char* option)
{
  int i,io,jo;

  int nco = basis->nco;
  double s;
  double row[nco];
  double row2[nco];
  double gto[nco];

  for (int i=0;i<n;i++) 
  {
    GetBasis(basis,x[i],y[i],z[i],gto); 
    v[i] = 0.0;
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
      s = row[io]*gto[io];
      for (jo=0;jo<io;jo++) 
        s += 2.0*row[jo]*gto[jo];
      v[i] += s*gto[io];
    }
  }
};
*/

void GetDensityLikeDerivatives(Matrix *P,Matrix* Q,Basis* basis,
double* x, double *y, double *z, double *dx, double *dy, double *dz, 
int n, const char* option)
{
  int io,jo;

  int nco = basis->nco;
  double row[nco];
  double row2[nco];

  double gto[nco];
  double gto_dx[nco];
  double gto_dy[nco];
  double gto_dz[nco];

  for (int i=0;i<n;i++) 
  {
    GetBasis(basis,x[i],y[i],z[i],gto); 
    GetBasisDerivatives(basis,x[i],y[i],z[i],gto_dx,gto_dy,gto_dz);
    dx[i] = 0.0;
    dy[i] = 0.0;
    dz[i] = 0.0;
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
      {
        dx[i] += row[jo]*(gto_dx[io]*gto[jo]+gto[io]*gto_dx[jo]);
        dy[i] += row[jo]*(gto_dy[io]*gto[jo]+gto[io]*gto_dy[jo]);
        dz[i] += row[jo]*(gto_dz[io]*gto[jo]+gto[io]*gto_dz[jo]);
      }
    }
  }
};

/*
void GetParticleElectrostaticPotential(Species* spec,double* x,double *y,double *z,double *v, int n)
{
  size_t iset,ishl,jset,jshl;
  int ispin;
  int dao,i,iatom,igp,j,k,lla,llb,lmax,natom;
  double dist,xr,yr,zr;
  double ra[3],rb[3],rc[3];
  double **intblk;
  Shell *sa,*sb;


  // Dimensions
  lmax = spec->basis->LMax();
  dao = ((lmax+1)*(lmax+2))/2;
  natom = spec->sys->mol->Natom();

  // Memory allocation
  intblk = new double*[dao];
  for (i=0;i<dao;i++)
    intblk[i] = new double[dao];

  // Nuclear contribution
  for (igp=0;igp<n;igp++)
  {
    v[igp] = 0.0;
    for (iatom=0;iatom<natom;iatom++)
    {
      xr = spec->sys->mol->atom[iatom]->x - x[igp];
      yr = spec->sys->mol->atom[iatom]->y - y[igp];
      zr = spec->sys->mol->atom[iatom]->z - z[igp];
      dist = max(TOLNUM,sqrt(xr*xr + yr*yr + zr*zr));
      v[igp] += spec->sys->mol->atom[iatom]->charge/dist;
    }
  }

  for (iset=0;iset<spec->basis->set.size();iset++)
  {
    spec->basis->mol->atom[iset]->Double(ra); 
    for (jset=0;jset<spec->basis->set.size();jset++)
    {
      spec->basis->mol->atom[jset]->Double(rb); 
      for (ishl=0;ishl<spec->basis->set[iset]->shell.size();ishl++)
      {
        sa = spec->basis->set[iset]->shell[ishl];
        lla = spec->basis->set[iset]->ll + sa->ll;
        for (jshl=0;jshl<spec->basis->set[jset]->shell.size();jshl++)
        {
          sb = spec->basis->set[jset]->shell[jshl];
          llb = spec->basis->set[jset]->ll + sb->ll;
          for (igp=0;igp<n;igp++)
          {
            rc[0] = x[igp]; rc[1] = y[igp]; rc[2] = z[igp];
            spec->sys->gaussi->Core(ra,rb,rc,sa,sb,intblk,0,0,true);
            for (i=0;i<sa->nco;i++)
              for (j=0;j<sb->nco;j++)
                for (ispin=0;ispin<spec->nspin;ispin++)
                  v[igp] += intblk[i][j]*(*spec->P[ispin])(lla+i,llb+j);
          }
        }
      }
    }
  }

  // Free memory
  for (i=0;i<dao;i++)
    delete[] intblk[i];
  delete[] intblk;
}

*/
