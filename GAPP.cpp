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
#include <string.h>

#include <iostream>
#include <vector>

#include <Array.h>
#include <Basis.h>
#include <GAPP.h>
#include <Global.h>
/*
#include <Matrix.h>
*/
#include <Math.h>
#include <Species.h>
#include <System.h>
#include <Units.h>


using namespace std;

GAPP::GAPP(System* isys)
{
  sys = isys;

}

// Generalized one particle propagator in diagonal
// (quasi-particle) approximation. (29.03.25, RFM)
void GAPP::Diagonal(const char* imethod,int speca,int specb,int spina,int spinb,int orba,int orbb,double *omega)
{
  int nit;
  double seconst,ps;

  strcpy(method,imethod);

  cout << "GAPP::Diagonal "<<endl;
  cout << "  Species: "<<speca<<" and "<<specb<<endl;
  cout << "  Spin: "<<spina<<" and "<<spinb<<endl;
  cout << "  Orbital: "<<orba<<" and "<<orbb<<endl;

  // Adjust pole structure if G0W0
  // rfm g0w0e();

  // Initialize
  nit = 0;
  ps = 0.0;
  *omega = sys->species[speca]->moe[spina][orba];
  if (method[1]=='2')
    *omega += sys->species[specb]->moe[spinb][orbb];

  /*
  // C onstant self-energy
  ConstantSelfEnergy(&seconst,speca,specb,spina,spinb,orba,orbb);
  cout << "Constant Self-Energy (eV): "
       <<HartreeToeV(seconst)<<endl;
  *omega += seconst;
  */

  // Iterations
  NewtonRaphson(omega,&ps,speca,specb,spina,spinb,orba,orbb,&nit);

  // Report
  cout << method <<" "<<nit
                 <<" "<<-HartreeToeV(*omega)
	         <<" "<<ps<<endl;
}

// Performs Newton-Raphson root search.
void GAPP::NewtonRaphson(double* omega,double* ps,int speca,int specb,int spina,int spinb,int orba,int orbb,int *nit)
{
  int it,maxit;
  double corr,ds,omega0,pomega,s,tol;
  double se[2];

  cout << "Newton-Raphson"<<endl;
  // Initialize
  omega0 = *omega;
  tol = 0.001/HartreeToeV(1.0);
  maxit = 30;

  // Newton-Raphson iteration 
  *nit = 0;
  for (it=0;it<maxit;it++)
  {
    (*nit)++;
    pomega = *omega;
    SelfEnergy(*omega,se,speca,specb,spina,spinb,orba,orbb);
    s = se[0];
    ds = se[1];
    *ps = 1.0/(1.0-ds);
    *omega = (*ps)*(omega0 + s - ds*pomega);
    cout << "CYCLE "<<*nit<<" BE: "<<-HartreeToeV(*omega)<<endl;
    corr = *omega - pomega;
    if (ABS(corr)<tol) 
    {
      corr = *omega - omega0 - s;
      if (ABS(corr)>tol) *ps = -(*ps);
      return;
    }
  }
  *nit = 0;
}

void GAPP::SelfEnergy(double omega,double* se,int speca,int specb,int spina,int spinb,int orba,int orbb)
{
  if (strncmp(method,"G1PP",4) == 0 )
    SelfEnergy1PP(omega,se,speca,spina,orba);
/*
  else if (strncmp(method,"G2PP",4) == 0 )
    SelfEnergy2PP(omega,se,speca,specb,spina,spinb,orba,orbb);
    */
}

// Self-enery for reference MC propagator.
// S = v^tDVDVDv + [u^tDVDv] + u^tDu
// S' = v^tDVD'VDv + [v^tDVDVD'v] +
//      [u^tD'VDv+u^tDVD'v] + u^tD'u 
void GAPP::SelfEnergy1PP(double omega,double* se,int speca,int spina,int orba)
{
  int domos,dptr,dumos,d2hp,d2ph,dv;

  // Initialize 
  se[0] = 0.0;
  se[1] = 0.0;

  // No self-energy in required 
  if (strncmp(method,"G1PP-1",6) == 0 ) return;

  if (!((strncmp(method,"G1PP-2",6) == 0 )||
        (strncmp(method,"G1PP-P3",7) == 0 )||
        (strncmp(method,"G1PP-3",6) == 0 )|| 
        (strncmp(method,"G1PP-P4",7) == 0 ))) 
  {
    cout << "GAPP::SelfEnergy1PP, unknown method: "<<method<<endl;
    return;
  }

  // Dimensions
  domos = NumberOfHoles();
  dumos = NumberOfParticles();
  d2hp = domos*domos*dumos;
  d2ph = domos*dumos*dumos;
  dv = d2hp + d2ph;

  // Memory allocation 
  Array D(dv);
  Array FV(dv);
  Array FDV(dv);
  Array FDDV(dv);
  Array FDVDV(dv);
  Array FVDV(dv);
  Array FVDDV(dv);
  Array FU(dv);
  Array FDU(dv);
  Array FDDU(dv);

  // Eq. (25) of Ortiz JCP 1996 paper

  // v = (f|V|a_k)^(1),  u = (f|V|a_k)^(2)
  FVA(speca,spina,orba,&FV,&D,&dptr,true);
  if ((strncmp(method,"G1PP-P3",7) == 0 )||
      (strncmp(method,"G1PP-3",6) == 0 )|| 
      (strncmp(method,"G1PP-P4",7) == 0 )) 
  {
    if ((strncmp(method,"G1PP-P3",7) == 0 )&&
        (sys->OnlyElectrons())) FVA2(speca,spina,orba,&FU,false);
    else FVA2(speca,spina,orba,&FU,true);
  }

  // Apply energy denominators 
  FFF(omega,&D,&FV,&FDV,dptr);
  if ((strncmp(method,"G1PP-P3",7) == 0 )||
      (strncmp(method,"G1PP-3",6) == 0 )|| 
      (strncmp(method,"G1PP-P4",7) == 0 )) 
  {
    FFF(omega,&D,&FDV,&FDDV,dptr);
    if (strncmp(method,"G1PP-P3",7) == 0 )
    {
      FFF(omega,&D,&FU,&FDU,dptr);
      FFF(omega,&D,&FDU,&FDDU,dptr);
    }
  }

  // Apply V = (f|V|f)^(0,0) 
  if ((strncmp(method,"G1PP-P3",7) == 0 )||
      (strncmp(method,"G1PP-3",6) == 0 )|| 
      (strncmp(method,"G1PP-P4",7) == 0 )) 
  {
    if ((strncmp(method,"G1PP-P3",7) == 0 )&&
        (sys->OnlyElectrons())) FVF(&FDV,&FVDV,false);
    else FVF(&FDV,&FVDV,true);

    if (strncmp(method,"G1PP-P4",7) == 0 )
      FVF(&FDDV,&FVDDV,true);
  }

  // Apply energy denominators after V 
  if (strncmp(method,"G1PP-P4",7) == 0 )
    FFF(omega,&D,&FVDV,&FDVDV,dptr);

  // Second order 
  se[0] += FV.SumProduct(&FDV);
  se[1] -= FDV.SumProduct(&FDV);

  // P3 and 3
  if ((strncmp(method,"G1PP-P3",7) == 0 )||
      (strncmp(method,"G1PP-3",6) == 0 )|| 
      (strncmp(method,"G1PP-P4",7) == 0 )) 
  {
    int dsum;
    if ((strncmp(method,"G1PP-P3",7) == 0 )&&
        (sys->OnlyElectrons())) dsum = d2hp;
    else dsum = dptr;
    se[0] += FDV.SumProduct(&FU,dsum);
    se[1] -= FDDV.SumProduct(&FU,dsum);
    se[0] += FDV.SumProduct(&FVDV,dsum);
    se[1] -= 2.0*FDDV.SumProduct(&FVDV,dsum);
  }

  // Third order 
  if ((strncmp(method,"G1PP-3",6) == 0 )|| 
      (strncmp(method,"G1PP-P4",7) == 0 )) 
  {
    se[0] += FDV.SumProduct(&FU,dptr);
    se[1] -= FDDV.SumProduct(&FU,dptr);
  }

  // Fourth order. Eq. (5.1), Master's thesis of Tamayo
  if (strncmp(method,"G1PP-P4",7) == 0 )
  {
    se[0] += FDU.SumProduct(&FU,dptr);
    se[0] += 2.0*FDU.SumProduct(&FVDV,dptr);
    se[0] += FVDV.SumProduct(&FDVDV,dptr);
    se[1] -= FDDU.SumProduct(&FU,dptr);
    se[1] -= 2.0*FDDU.SumProduct(&FVDV,dptr);
    se[1] -= 2.0*FDVDV.SumProduct(&FVDDV,dptr);
    se[1] -= FDVDV.SumProduct(&FDVDV,dptr);
  }
}

/*
void GAPP::SelfEnergy2PP(double omega,double* se,int speca,int specb,int spina,int spinb,int orba,int orbb)
{
  int domos,dptr,dumos,d3hp,d3ph,dv;

  // Initialize 
  se[0] = 0.0;
  se[1] = 0.0;

  // No self-energy is required 
  if ((strncmp(method,"G2PP-0",6) == 0)||
      (strncmp(method,"G2PP-1",6) == 0)) return;

  if (!(strncmp(method,"G2PP-2",6) == 0 ))
  {
    cout << "GAPP::SelfEnergy2PP, unknown method: "<<method<<endl;
    return;
  }

  // Dimensions
  domos = NumberOfHoles();
  dumos = NumberOfParticles();
  d3hp = domos*domos*domos*dumos;
  d3ph = domos*dumos*dumos*dumos;
  dv = d3hp + d3ph;

  // Memory allocation 
  Array D(dv);
  Array GV(dv);

  // GVB
  GVB(speca,specb,spina,spinb,orba,orbb,&GV,&D);
}
*/

// Apply matrix [w 1-(f|F|f)]^{-1} (derivatives) to vector. 
void GAPP::FFF(double omega,Array* d,Array* v,Array* u,int n)
{
  int i;
  double den;

  u->SetZero();
  for (i=0;i<n;i++)
  {
    den = omega + d->Value(i);
    if (ABS(den)>TOLNUM) u->SetValue(i,v->Value(i)/den);
  }
}

int GAPP::NumberOfHoles()
{
  int dbas,domos;
  size_t ispec,ispin;

  domos = 0;
  for (ispec=0;ispec<sys->species.size();ispec++)
  {
    dbas = sys->species[ispec]->basis->nco;
    for (ispin=0;ispin<sys->species[ispec]->nspin;ispin++)
    {
      domos += sys->species[ispec]->npart[ispin];
    }
  }
  return domos;
}

int GAPP::NumberOfParticles()
{
  int dbas,dumos;
  size_t ispec,ispin;

  dumos = 0;
  for (ispec=0;ispec<sys->species.size();ispec++)
  {
    dbas = sys->species[ispec]->basis->nco;
    for (ispin=0;ispin<sys->species[ispec]->nspin;ispin++)
    {
      dumos += (dbas - sys->species[ispec]->npart[ispin]);
    }
  }
  return dumos;
}

// Load all eigenvalues ordered according to particle type.
void GAPP::GetEigenvalues(double *eigs)
{
  size_t a,i,ispec,ispin;
  int dbas,domo,domos,dorb,dumo,fo,fu;

  // Dimensions 
  domos = NumberOfHoles();

  fo = 0;
  fu = domos;
  for (ispec=0;ispec<sys->species.size();ispec++)
  {
    dbas = sys->species[ispec]->basis->nco;
    dorb = dbas;
    for (ispin=0;ispin<sys->species[ispec]->nspin;ispin++)
    {
      domo = sys->species[ispec]->npart[ispin];
      dumo = dorb - domo;

      for (i=0;i<domo;i++)
        eigs[fo+i] = sys->species[ispec]->moe[ispin][i];
      for (a=domo;a<dorb;a++)
        eigs[fu+a-domo] = sys->species[ispec]->moe[ispin][a];
      fo += domo;
      fu += dumo;
    }
  }
}

// Evaluates <P1P2|P3P4> for a given P1, spin-orbitals.
void GAPP::SOERI(int speca,int spina,int orba,vector<Matrix> ERIs)
{
  size_t ispec,ispin;
  int a,b,dbas,dbas2,fo,fo2,fu,fu2,i,j,t,tt;
  int domo,domos,domo2,dumo,dumo2,dorb,dorb2;

  // Initialize 
  for (i=0;i<ERIs.size();i++)
    ERIs[i].SetZero();

  // Dimensions
  domos = NumberOfHoles();
  domo = sys->species[speca]->npart[spina];
  dorb = sys->species[speca]->basis->nco;
  dbas = sys->species[speca]->basis->nco;
  dumo = dorb - domo;

  fo = 0; 
  fu = domos;
  for (ispec=0;ispec<sys->species.size();ispec++)
  {
    for (ispin=0;ispin<sys->species[ispec]->nspin;ispin++)
    {
      domo2 = sys->species[ispec]->npart[ispin];
      dorb2 = sys->species[ispec]->basis->nco;
      dumo2 = dorb2 - domo2;
      if ((speca==ispec)&&(spina==ispin)) break;
      fo += domo2;
      fu += dumo2;
    }
  }

  fo2 = 0; 
  fu2 = domos;
  for (ispec=0;ispec<sys->species.size();ispec++)
  {
    dbas2 = sys->species[ispec]->basis->nco;
    vector<Matrix> eri;
    eri.clear();
    for (int ieri=0;ieri<dbas2;ieri++)
      eri.push_back( Matrix(dbas,dbas2) );
    for (ispin=0;ispin<sys->species[ispec]->nspin;ispin++)
    {
      domo2 = sys->species[ispec]->npart[ispin];
      dorb2 = sys->species[ispec]->basis->nco;
      dumo2 = dorb2 - domo2;

      sys->species[speca]->TwoParticleIITensor(sys->species[ispec],
	 		                       spina,ispin,orba,eri);
      for (t=0;t<dorb2;t++)
      {
        if (t<domo2) tt = t + fo2;
        else tt = t - domo2 + fu2;

        // Copy: Mulliken notation => Dirac notation

        // Occupied-occupied 
        for (i=0;i<domo;i++)
        {
          for (j=0;j<domo2;j++)
          {
            ERIs[tt].ShiftValue(fo+i,fo2+j,eri[t](i,j));
            ERIs[tt].ShiftValue(fo2+j,fo+i,-eri[t](i,j));
          }
        }

        // Occupied-virtual 
        for (i=0;i<domo;i++)
        {
          for (a=domo2;a<dorb2;a++)
          {
            ERIs[tt].ShiftValue(fo+i,fu2+a-domo2,eri[t](i,a));
            ERIs[tt].ShiftValue(fu2+a-domo2,fo+i,-eri[t](i,a));
          }
        }

        // Virtual-occupied
	for (i=0;i<domo2;i++)
	{
	  for (a=domo;a<dorb;a++)
	  {
            ERIs[tt].ShiftValue(fu+a-domo,fo2+i,eri[t](a,i));
            ERIs[tt].ShiftValue(fo2+i,fu+a-domo,-eri[t](a,i));
	  }
	}

        // Virtual-virtual 
	for (a=domo;a<dorb;a++)
	{
	  for (b=domo2;b<dorb2;b++)
	  {
            ERIs[tt].ShiftValue(fu+a-domo,fu2+b-domo2,eri[t](a,b));
            ERIs[tt].ShiftValue(fu2+b-domo2,fu+a-domo,-eri[t](a,b));
	  }
	}
      }
      fo2 += domo2;
      fu2 += dumo2;
    }
    eri.clear();
  }
}

// Build (|f_{2ph}|V|a_p).
void GAPP::FVA(int speca,int spina,int orba,Array* v,Array* d,int *dptr,bool do2ph)
{
  int a,b,dh,dhp,dp,i,j;

  // Dimensions
  dh = NumberOfHoles();
  dp = NumberOfParticles();
  dhp = dh + dp;

  // Memory allocation 
  vector<Matrix> ERIs;
  ERIs.clear();
  for (i=0;i<dhp;i++)
    ERIs.push_back( Matrix(dhp,dhp) );

  // Get integrals 
  SOERI(speca,spina,orba,ERIs);

  // Load eigenvalues 
  double eig[dhp];
  GetEigenvalues(eig);

  // 2hp
  int k;
  *dptr = 0;
  for (i=0;i<dh;i++)
  {
    for (j=0;j<dh;j++)
    {
      for (a=dh;a<dhp;a++)
      {
	v->SetValue(*dptr,ERIs[a](i,j));
	d->SetValue(*dptr,eig[a] - eig[i] - eig[j]);
	(*dptr)++;
      }
    }
  }

  // 2ph 
  if (do2ph)
  {
    for (a=dh;a<dhp;a++)
    {
      for (b=dh;b<dhp;b++)
      {
        for (i=0;i<dh;i++)
	{
	  v->SetValue(*dptr,ERIs[i](a,b));
	  d->SetValue(*dptr,eig[i] - eig[a] - eig[b]);
	  (*dptr)++;
	}
      }
    }
  }
  
  // Scale
  v->Scale( 1.0/sqrt(2.0) );
  
  // Free allocation 
  ERIs.clear();
}

// Build (|f_{2ph}|V|a_p).
void GAPP::FVA2(int speca,int spina,int orba,Array* v,bool do2ph)
{
  int p,speco,spino,orbo;
  int a,b,c,dh,dhp,dp,i,iptr,j,k;
  double d,s;

  dh = NumberOfHoles();
  dp = NumberOfParticles();
  dhp = dh + dp;

  // Load eigenvalues 
  double eig[dhp];
  GetEigenvalues(eig);

  // Get integrals 
  vector<Matrix> ERIsP;
  vector<Matrix> ERIs;
  ERIsP.clear();
  ERIs .clear();
  for (i=0;i<dhp;i++)
  {
    ERIsP.push_back( Matrix(dhp,dhp) );
    ERIs.push_back( Matrix(dhp,dhp) );
  }
  SOERI(speca,spina,orba,ERIsP);

  for (i=0;i<dh;i++)
  {
    Decode(&i,&speco,&spino,&orbo,"occupied");
    SOERI(speco,spino,orbo,ERIs);
    for (j=0;j<dh;j++)
    {
      for (a=dh;a<dhp;a++)
      {
        s = 0.0;
        for (b=dh;b<dhp;b++)
        {
          // Sum_bc <pa||bc> c_ij^bc 
          for (c=b+1;c<dhp;c++)
          {
            d = eig[i] + eig[j] - eig[b] - eig[c];
            s = s + ERIsP[a](b,c)*ERIs[j](b,c)/d;
	  }
          // - Sum_bk <pk||jb> c_ik^ab 
          for (k=0;k<dh;k++)
          {
            d = eig[i] + eig[k] - eig[a] - eig[b];
            s = s - 2.0*ERIsP[k](j,b)*ERIs[k](a,b)/d;
	  }
	}
        iptr = Map2HP(i,j,a,dh,dp);
        v->SetValue(iptr,s/sqrt(2.0));
      }
    }
  }

  if (do2ph)
  {
    for (a=dh;a<dhp;a++)
    {
      Decode(&a,&speco,&spino,&orbo,"virtual");
      SOERI(speco,spino,orbo,ERIs);
      for (b=dh;b<dhp;b++)
      {
        for (i=0;i<dh;i++)
        {
          s = 0.0;
          for (j=0;j<dh;j++)
          {
            // Sum_jk <pi||jk> c_jk^ab 
            for (k=j+1;k<dh;k++)
            {
              d = eig[j] + eig[k] - eig[a] - eig[b];
              s = s + ERIsP[i](j,k)*ERIs[b](j,k)/d;
	    }
            // - Sum_ck <pc||bk> c_ik^ac 
            for (c=dh;c<dhp;c++)
            {
              d = eig[i] + eig[j] - eig[a] - eig[c];
              s =  s - 2.0*ERIsP[c](b,j)*ERIs[c](i,j)/d;
	    }
	  }
          iptr = Map2PH(a,b,i,dh,dp);
          v->SetValue(iptr,s/sqrt(2.0));
	}
      }
    }
  }
  
  // Free memory
  ERIsP.clear();
  ERIs.clear();
}

// Apply matrix (f_{2hp}|V|f_{2hp})^{(1)} to vector, V_{2hp}. 
void GAPP::FVF(Array* vi,Array* vf,bool do2ph)
{
  int speco,spino,orbo;
  int a,b,c,d,dh,dhp,dp,i,iptr,j,jptr,k,l;
  double s1,s2;

  dh = NumberOfHoles();
  dp = NumberOfParticles();
  dhp = dh + dp;

  // Allocate
  vector<Matrix> ERIs;
  ERIs .clear();
  for (i=0;i<dhp;i++)
    ERIs.push_back( Matrix(dhp,dhp) );

  // Multiply V_{ff} matrix by vector 
  for (i=0;i<dh;i++)
  {
    Decode(&i,&speco,&spino,&orbo,"occupied");
    SOERI(speco,spino,orbo,ERIs);
    for (j=0;j<dh;j++)
    {
      for (a=dh;a<dhp;a++)
      {
        s1 = 0.0;
        s2 = 0.0;
        for (l=0;l<dh;l++)
        {
          for (k=0;k<dh;k++)
          {
            jptr = Map2HP(k,l,a,dh,dp);
            s1 += vi->Value(jptr)*ERIs[j](k,l);
	  }
          for (b=dh;b<dhp;b++)
          {
            jptr = Map2HP(j,l,b,dh,dp);
            s2 += vi->Value(jptr)*ERIs[b](a,l);
	  }
        }
        iptr = Map2HP(i,j,a,dh,dp);
        vf->SetValue(iptr,-0.5*s1+2.0*s2); 
      }
    }
  }

  if (do2ph)
  {
    for (a=dh;a<dhp;a++)
    {
      Decode(&a,&speco,&spino,&orbo,"virtual");
      SOERI(speco,spino,orbo,ERIs);
      for (b=dh;b<dhp;b++)
      {
        for (i=0;i<dh;i++)
	{
          s1 = 0.0;
          s2 = 0.0;
          for (d=dh;d<dhp;d++)
          {
            for (c=dh;c<dhp;c++)
            {
              jptr = Map2PH(c,d,i,dh,dp);
              s1 += vi->Value(jptr)*ERIs[b](c,d);
	    }
            for (j=0;j<dh;j++)
	    {
              jptr = Map2PH(b,d,j,dh,dp);
              s2 += vi->Value(jptr)*ERIs[j](i,d);
	    }
          }
          iptr = Map2PH(a,b,i,dh,dp);
          vf->SetValue(iptr,0.5*s1-2.0*s2); 
	}
      }
    }
  }
  
  // Free memory
  ERIs.clear();
}

// Identify location of generic spin-orbital
void GAPP::Decode(int* so,int* speco,int* spino,int* orbo, const char* option)
{
  size_t ispec,ispin;
  int dh,dhp,dp,iorb;
  int k;

  if ( strncmp(option,"occup",5) == 0 ) 
  {
    k = 0; 
    for (ispec=0;ispec<sys->species.size();ispec++)
    {
      for (ispin=0;ispin<sys->species[ispec]->nspin;ispin++)
      {
        dh = sys->species[ispec]->npart[ispin];
        for (iorb=0;iorb<dh;iorb++)
        {
	  if (*so<0)
	  {
	    if ((*speco==ispec)&&(*spino==ispin)&&(*orbo==iorb))
	    {
              *so = k;
	      return;
	    }
	  }
	  else
	  {
            if (k==*so)
	    {
	      *speco = ispec;
	      *spino = ispin;
	      *orbo = iorb;
	      return;
	    }
	  }
          k++;
	}
      }
    }
  }
  else
  {
    k = NumberOfHoles();
    for (ispec=0;ispec<sys->species.size();ispec++)
    {
      for (ispin=0;ispin<sys->species[ispec]->nspin;ispin++)
      {
        dhp = sys->species[ispec]->basis->nco;
        dh = sys->species[ispec]->npart[ispin];
	dp = dhp - dh;
        for (iorb=dh;iorb<dhp;iorb++)
        {
	  if (*so<0)
	  {
	    if ((*speco==ispec)&&(*spino==ispin)&&(*orbo==iorb))
	    {
              *so = k;
	      return;
	    }
	  }
	  else
	  {
            if (k==*so)
	    {
	      *speco = ispec;
	      *spino = ispin;
	      *orbo = iorb;
	      return;
	    }
	  }
          k++;
	}
      }
    }
  }
}


int GAPP::Map2HP(int i,int j,int a,int dh, int dp)
{
   return  a - dh + dp*(j + dh*i);
}

int GAPP::Map2PH(int a,int b,int i,int dh, int dp)
{
   return  i + dh*(b - dh + dp*a);
}

/*
// Evaluation of constant self-energy for GAPP.
void GAPP::ConstantSelfEnergy(double* cse, int speca, int specb, int spina, int spinb, int orba, int orbb)
{
  if (strncmp(method,"G1PP",4) == 0 )
    ConstantSelfEnergy1PP(cse,speca,spina,orba);
  else if (strncmp(method,"G2PP",4) == 0 )
    ConstantSelfEnergy2PP(cse,speca,specb,spina,spinb,orba,orbb);
}
void GAPP::ConstantSelfEnergy1PP(double* cse, int speca, int spina, int orba)
{
  *cse = 0.0;
  if ((strncmp(method,"G1PP-1",6) == 0 )|| 
      (strncmp(method,"G1PP-2",6) == 0 )|| 
      (strncmp(method,"G1PP-P3",7) == 0 )) return;

  // Third order (a|V|a) terms 
  AVA3(cse,speca,spina,orba);
}
void GAPP::ConstantSelfEnergy2PP(double* cse, int speca, int specb, int spina, int spinb, int orba, int orbb)
{
  *cse = 0.0;
  if (strncmp(method,"G2PP-0",6) == 0 ) return;

  // First order: (aa|V|aa)^{(0,0)}
  AAVAA1(cse,speca,specb,spina,spinb,orba,orbb);
}

// Build (|f_{2ph}|V|a_p).
void GAPP::AVA3(double* ava3, int speca, int spina, int orba)
{
  int p,speco,spino,orbo;
  int a,b,c,dh,dhp,dp,i,iptr,j,k;
  double d,s,t1,t2;

  // Dimensions
  dh = NumberOfHoles();
  dp = NumberOfParticles();
  dhp = dh + dp;

  // Initialize 
  *ava3 = 0.0;

  // Load eigenvalues 
  double eig[dhp];
  GetEigenvalues(eig);

  // Get integrals 
  vector<Matrix> ERIsP;
  vector<Matrix> ERIs;
  ERIsP.clear();
  ERIs .clear();
  for (i=0;i<dhp;i++)
  {
    ERIsP.push_back( Matrix(dhp,dhp) );
    ERIs.push_back( Matrix(dhp,dhp) );
  }
  SOERI(speca,spina,orba,ERIsP);
  Decode(&p,&speca,&spina,&orba,"occupied");

  for (j=0;j<dh;j++)
  {
    Decode(&j,&speco,&spino,&orbo,"occupied");
    SOERI(speco,spino,orbo,ERIs);
    for (i=0;i<dh;i++)
    {
      // -\sum_{ijkab} <pk||pi> c_{jk}^{ab}(1)c_{ji}^{ba}(1) 
      for (k=0;k<dh;k++)
      {
        s = 0.0;
        for (a=dh;a<dhp;a++)
        {
          for (b=dh;b<dhp;b++)
          {
            // First order double excitation amplitudes 
            // c_{ij}^{ab}(1) = <ij||ab>/(e_i+e_j-e_a-e_b) 
            d = eig[i] + eig[j] - eig[a] - eig[b];
            t1 = ERIs[i](b,a)/d;
            d = eig[j] + eig[k] - eig[a] - eig[b];
            t2 = ERIs[k](a,b)/d;
            s = s + t1*t2;
          }
        }
        *ava3 += 0.5*ERIsP[k](p,i)*s;
      }
      // Second order single excitation amplitudes 
      // c_{i}^{a}(2)=\sum_{jkb} c_{jk}^{ba}(1)<jk||bi>/(e_i-e_a)
      //            -\sum_{jbc} c_{ji}^{bc}(1)<bc||ja>/(e_i-e_a)
      for (a=dh;a<dhp;a++)
      {
        s = 0.0;
        for (b=dh;b<dhp;b++)
        {
          for (k=0;k<dh;k++)
          {
            d = eig[j] + eig[k] - eig[a] - eig[b];
            t1 = ERIs[k](b,a)/d;
            d = eig[i] - eig[a];
            s -= t1*ERIs[k](b,i)/d;
          }
          for (c=dh;c<dhp;c++)
          {
            d = eig[j] + eig[i] - eig[c] - eig[b];
            t1 = ERIs[i](b,c)/d;
            d = eig[i] - eig[a];
            s += t1*ERIs[a](b,c)/d;
          }
        }
        // \sum_{ai} <pa||pi> c_{i}^{a}(2) 
        *ava3 += ERIsP[a](p,i)*s;
      }
    }
  }

  // \sum_{ijabc} <pa||pb> c_{ij}^{bc}(1)c_{ji}^{ac}(1) 
  for (c=dh;c<dhp;c++)
  {
    Decode(&c,&speco,&spino,&orbo,"virtual");
    SOERI(speco,spino,orbo,ERIs);
    for (a=dh;a<dhp;a++)
    {
      for (b=dh;b<dhp;b++)
      {
        s = 0.0;
        for (i=0;i<dh;i++)
        {
          for (j=0;j<dh;j++)
          {
            // First order double excitation amplitudes 
            //  c_{ij}^{ab}(1) = <ij||ab>/(e_i+e_j-e_a-e_b) 
            d = eig[i] + eig[j] - eig[b] - eig[c];
            t1 = ERIs[b](j,i)/d;
            d = eig[i] + eig[j] - eig[a] - eig[c];
            t2 = ERIs[a](i,j)/d;
            s += t1*t2;
	  }
        }
        *ava3 -= 0.5*ERIsP[a](p,b)*s;
      }
    }
  }

  // Free memory
  ERIsP.clear();
  ERIs .clear();
}

// Build (aa|Vaa)^{(0,0)} Two particle propagator
void GAPP::AAVAA1(double* aavaa1, int speca, int specb, int spina,int spinb, int orba, int orbb)
{
  size_t p;
  int a,b,c,dh,dhp,dp,i;

  // Dimensions
  dh = NumberOfHoles();
  dp = NumberOfParticles();
  dhp = dh + dp;

  // Get integrals 
  vector<Matrix> ERIsP;
  ERIsP.clear();
  for (i=0;i<dhp;i++)
  {
    ERIsP.push_back( Matrix(dhp,dhp) );
  }
  SOERI(speca,spina,orba,ERIsP); // FIXME: Too many transformations

  Decode(&i,&specb,&spinb,&orbb,"occupied");
  *aavaa1 = ERIsP[i](p,i);

  // Free memory
  ERIsP.clear();
}

// Build (g|V|b) for G2PP.
void GAPP::GVB(int speca,int specb,int spina,int spinb,int orba,int orbb,Array* v,Array* d)
{
  int a,b,c,dh,dhp,dp,i,j,k,p,ptr,q;

  // Dimensions
  dh = NumberOfHoles();
  dp = NumberOfParticles();
  dhp = dh + dp;

  // Memory allocation 
  vector<Matrix> ERIsA;
  vector<Matrix> ERIsB;
  ERIsA.clear();
  ERIsB.clear();
  for (i=0;i<dhp;i++)
  {
    ERIsA.push_back( Matrix(dhp,dhp) );
    ERIsB.push_back( Matrix(dhp,dhp) );
  }

  // Get integrals 
  SOERI(speca,spina,orba,ERIsA);
  SOERI(specb,spinb,orbb,ERIsB);

  // Load eigenvalues 
  double eig[dhp];
  GetEigenvalues(eig);

  // Identify active ones
  Decode(&p,&speca,&spina,&orba,"occupied");
  Decode(&q,&specb,&spinb,&orbb,"occupied");

  // 3hp
  ptr = 0;
  for (i=0;i<dh;i++)
  {
    for (j=0;j<dh;j++)
    {
      for (a=dh;a<dhp;a++)
      {
        v->SetValue(ptr,ERIsA[a](i,j));
        d->SetValue(ptr,eig[a] - eig[i] - eig[j]);
        ptr++;
      }
    }
  }

  // 3ph 
  for (a=dh;a<dhp;a++)
  {
    for (b=dh;b<dhp;b++)
    {
      for (i=0;i<dh;i++)
      {
        //RFM v->SetValue(*dptr,ERIs[i](a,b));
        d->SetValue(ptr,eig[i] - eig[a] - eig[b]);
        ptr++;
      }
    }
  }
  
  // Scale
  v->Scale( 1.0/sqrt(2.0) );
  
  // Free allocation 
  ERIsA.clear();
  ERIsB.clear();
}


*/
