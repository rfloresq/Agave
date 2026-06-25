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
#include <string.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <Global.h>
#include <Integrator.h>
#include <Species.h>
#include <Structure.h>
#include <System.h>

/*
#include <ADPT.h>
#include <Array.h>
#include <Basis.h>
#include <Matrix.h>
#include <Set.h>
#include <Shell.h>
#include <Spin.h>
#include <Symmetry.h>
#include <Atom.h>
#include <Units.h>
*/

using namespace std;

System::System()
{
  mol = new Structure();

  gaussi = new Integrator();

	/*
  G = 0;
  GI = 0;
  GIS = 0;
  */
};



/*
void System::ReadASC(const char *filename)
{
  ifstream f(filename);
  if ( f.fail() )
  {
    cout << "Error opening file "<<filename<<endl;
    return;
  };

  char line[MAXSTR];
  char str[MAXSTR];
  double x,y,z,v;

  f.getline(line,MAXSTR);
  f.getline(line,MAXSTR);
  f.getline(line,MAXSTR);
  f >> line;
  f >> line;
  f >> line;
  f >> fname;

  fx.clear();
  fy.clear();
  fz.clear();
  fv.clear();
  while (f.getline(line,MAXSTR))
  {
    stringstream ss(line);
    ss >> x;
    ss >> y;
    ss >> z;
    ss >> v;
    fx.push_back( x );
    fy.push_back( y );
    fz.push_back( z );
    fv.push_back( v );
  }
  f.close();
}
*/

void System::ReadFile(const char*filename,const char*fmt)
{
	/*
  if ( strcasecmp( fmt , "ASC" ) == 0 ) 
  {
    ReadASC(filename);
    return;
  }
  */

  // Load structure
  mol->ReadFile(filename,fmt);

  if ( strcasecmp( fmt , "XYZ" ) == 0 ) return;
  // if ( strcasecmp( fmt , "PIE" ) == 0 ) return;

  
  if ( strcasecmp( fmt , "AGV" ) == 0 ) 
  {
    ReadAGV(filename);
  } 
  else if ( strcasecmp( fmt , "D2K" ) == 0 ) 
  {
    ReadD2K(filename);
  } 

  for (int is=0;is<(signed)species.size();is++)
    species[is]->ReadFile(filename,fmt,is);
  
};

void System::WriteFile(const char*filename,const char*fmt)
{
  ofstream f(filename);
  if ( strcasecmp( fmt , "XYZ" ) != 0 ) f << " "<<endl;
  f.close();

  mol->WriteFile(filename);
  if ( strcasecmp( fmt , "XYZ" ) == 0 ) return;

  f = ofstream(filename,ios::app);
  f << "Agave input file "<<endl;
  f << "Number of species: "<<species.size()<<endl;
  f.close();

  for (int is=0;is<(signed)species.size();is++)
    species[is]->WriteFile(filename,is);
  
};

bool System::PrepareForSCF(vector<Species*> winspecies)
{
  size_t ispec;

  // Copy window species with particles
  species.clear();
  for (ispec=0;ispec<winspecies.size();ispec++)
    if (winspecies[ispec]->TotalNumberOfParticles()>0)
    {
      species.push_back( new Species(this,winspecies[ispec]->type) ); 
      species[species.size()-1]->CopyParameters( winspecies[ispec] );
    }
  if (species.size()<1)
  {
    cout << "TOTAL NUMBER OF PARTICLES IS ZERO!"<<endl;
    return false;
  }

  //symmetry = new Symmetry(this);
  //symmetry->PointGroup();
  //cout << "POINT GROUP: "<<symmetry->point_group<<endl;

  for (ispec=0;ispec<species.size();ispec++)
    species[ispec]->PrepareForSCF();

  /*
  for (ispec=0;ispec<species.size();ispec++)
  {
    species[ispec]->auxis = new Basis(mol);
    species[ispec]->auxis->GenerateAuxiliary( species[ispec]->basis, gaussi );
    if (species[ispec]->type==POSITRON)
      auxis_unique = species[ispec]->auxis;
  }

  // FIXME: Setting unique basis to electronic
  for (ispec=0;ispec<species.size();ispec++)
    species[ispec]->auxis = auxis_unique;

  for (ispec=0;ispec<species.size();ispec++)
    species[ispec]->PrepareForFitting();

  int dabas;

  dabas = auxis_unique->nco;

  if (dabas<1) cout << "NO AUXILIARY FUNCTIONS FOUND"<<endl;
  else
  {
    if (G) delete G;
    if (GI) delete GI;
    if (GIS) delete GIS;
    G = new Matrix(dabas,dabas);
    GI = new Matrix(dabas,dabas);
    GIS = new Matrix(dabas,dabas);
    FittingMatrix();
  }
  */
  return true;
};

/*
void System::PrepareForLinearResponse()
{
  for (int is=0;is<(signed)species.size();is++)
    species[is]->PrepareForLinearResponse();
};

bool System::PrepareForRT()
{
  bool res;
  for (size_t is=0;is<species.size();is++)
  {
    res = species[is]->PrepareForRT();
    if (!res) break;
  }
  return res;
};
*/

double System::Energy(bool electric)
{
  double energy = 0.0;
  for (int is=0;is<(signed)species.size();is++)
    energy += species[is]->Energy(electric);
  
  // Nuclear repulsion
  double nre = mol->NuclearRepulsionEnergy();
  energy += nre;

  return energy;
}

/*
void System::Dipole(double* dipole)
{
  int i;
  size_t ispec;
  double dt[3];

  for (i=0;i<3;i++)
    dipole[i] = 0.0;

  for (ispec=0;ispec<species.size();ispec++)
  {
    species[ispec]->Dipole(dt);
    for (i=0;i<3;i++)
      dipole[i] += dt[i];
  }
  
  // Nuclear dipole
  mol->NuclearDipole(dt);
  for (i=0;i<3;i++)
    dipole[i] += dt[i];
}
*/


bool System::OnlyElectrons()
{
  bool oe;
  size_t ispec;

  oe = true;
  for (ispec=0;ispec<species.size();ispec++)
    if (species[ispec]->type!=ELECTRON)
    {
      oe = false;
      break;
    }

  return oe;
}

/*
void System::FittingMatrix()
{
  BuildFittingMatrix(G,auxis_unique);

  // G^-1/2
  double *values;
  values = new double[auxis_unique->nco];
  
  GI->Copy(G);
  GI->SVDPower(values,double(-1.0),double(1.0e-8));

  GIS->Copy(G);
  GIS->SVDPower(values,double(-0.5),double(1.0e-7));

  delete values;
}

void System::BuildFittingMatrix(Matrix* G,Basis* auxis)
{
  int i,ib,is,j,jb,js,lla,llb;
  double ra[3],rb[3];
  double blk[MAXNCO][MAXNCO];
  Shell *sa,*sb;

  int lmax = auxis->LMax();
  int dblk = ((lmax+1)*(lmax+2))/2;

  G->SetZero();

  // Build G matrix
  for (ib=0;ib<(signed)auxis->set.size();ib++)
  {
    ra[0] = auxis->mol->atom[ib]->x; 
    ra[1] = auxis->mol->atom[ib]->y; 
    ra[2] = auxis->mol->atom[ib]->z; 
    for (is=0;is<(signed)auxis->set[ib]->shell.size();is++)
    {
      sa = auxis->set[ib]->shell[is];
      lla = auxis->set[ib]->ll+sa->ll;
      for (jb=0;jb<(signed)auxis->set.size();jb++)
      {
        rb[0] = auxis->mol->atom[jb]->x; 
        rb[1] = auxis->mol->atom[jb]->y; 
        rb[2] = auxis->mol->atom[jb]->z; 
        for (js=0;js<(signed)auxis->set[jb]->shell.size();js++)
        {
          sb = auxis->set[jb]->shell[js];
          llb = auxis->set[jb]->ll+sb->ll;

          gaussi->Interaction2(ra,rb,sa,sb,blk,0);
          for (i=0;i<sa->nco;i++)
            for (j=0;j<sb->nco;j++)
              G->SetValue(lla+i,llb+j,blk[i][j]);
        }
      }
    }
  }
}

void System::SolveGxb(Array* x)
{
  int dabas;
  Array *tmp;

  dabas = auxis_unique->nco;
  tmp = new Array(dabas);

  tmp->Copy(x);
  tmp->Transform(GI,x);

  delete tmp;
}
*/

/*
void System::EvaluateFittingCoefficients()
{
  // Big matrix
  int dabas,dabaso,dbig,i,j,lla,llb;
  size_t ispec,ispin,jspec,jspin,nspec,nspina,nspinb;
  double *bvec,*xvec;

  dbig = 0;
  for (ispec=0;ispec<species.size();ispec++)
  {
    dabas = species[ispec]->auxis->nco;
    //for (ispin=0;ispin<species[ispec]->spin.size();ispin++)
    dbig += dabas;
  }
  Matrix Gbig(dbig,dbig);
  bvec = new double[dbig];
  xvec = new double[dbig];

  Gbig.SetZero();
  for (i=0;i<dbig;i++)
  {
    bvec[i] = 0.0;
    xvec[i] = 0.0;
  }

  lla = 0;
  nspec = species.size();
  for (ispec=0;ispec<nspec;ispec++)
  {
    dabas = species[ispec]->auxis->nco;
    Array jvec(dabas);
    nspina = species[ispec]->spin.size();
    //for (ispin=0;ispin<nspina;ispin++)
    //{
      llb = 0;
      for (jspec=0;jspec<nspec;jspec++)
      {
        double factor = species[ispec]->Charge()*
	                species[jspec]->Charge();
        dabaso = species[jspec]->auxis->nco;
        Matrix Gmix(dabas,dabaso);
        species[ispec]->BuildFittingMatrix(species[jspec],&Gmix);
	//Gmix.Scale(factor);
        nspinb = species[jspec]->spin.size();
	for (i=0;i<dabas;i++)
	  for (j=0;j<dabaso;j++)
            Gbig.ShiftValue(lla+i,llb+j,Gmix(i,j) );
            //Gbig.ShiftedAdd( &Gmix, lla, llb );
        for (jspin=0;jspin<nspinb;jspin++)
        {
	    jvec.SetZero();
            species[jspec]->spin[jspin]->FittingVector(
	  	        species[jspec]->basis,species[ispec]->auxis,
	    	        species[jspec]->spin[jspin]->P,&jvec);
	    jvec.Scale(species[jspec]->Charge());
	    for (i=0;i<dabas;i++)
	      bvec[lla+i] += jvec.Value(i);
	}
          llb += dabaso;
	//}
      }
      lla += dabas;
    //}
  }

  //rfm Gbig.SVDPower(xvec,double(-1.0),double(1.0e-8));
  InvertSuper(Gbig,dbig);
  Gbig.VectorMultiply(bvec,xvec);

  lla = 0;
  for (ispec=0;ispec<species.size();ispec++)
  {
    dabas = species[ispec]->auxis->nco;
    for (ispin=0;ispin<species[ispec]->spin.size();ispin++)
    {
      species[ispec]->spin[ispin]->x->SetZero();
      for (i=0;i<dabas;i++)
      {
        if (species[ispec]->type==ELECTRON)
          species[ispec]->spin[ispin]->x->ShiftValue(i,0.5*xvec[lla+i]);
	else if (species[ispec]->type==POSITRON)
          if (ispin==0) 
            species[ispec]->spin[ispin]->x->ShiftValue(i,xvec[lla+i]);
	  else
            species[ispec]->spin[ispin]->x->ShiftValue(i,0.0);
      }
    }
    lla += dabas;
  }
  delete bvec;
  delete xvec;
}

void System::InvertSuper(Matrix r, int n)
{
  Matrix raux(n,n);
  Matrix tmp(n,n);

  raux.Copy( &r );
  raux.Transpose();
  raux.Multiply(&r,&tmp,false);

  double values[n];
  tmp.SVDPower(values,double(-1.0),double(1.0e-10));
  tmp.Multiply(&r,&raux,true);
  raux.Transpose();
  r.Copy( &raux );
}
*/

void System::ReadAGV(const char*filename)
{
  char line[MAXSTR];
  int nspec;

  ifstream f(filename);
  if ( f.fail() )
  {
    cout<<__FILE__<<__LINE__<<"Error opening file"<<endl;
    return;
  };

  nspec = 0;
  while (f.getline(line,MAXSTR))
  {
    int sindex = -1;
    string stringLine=line;
    sindex=stringLine.find("Number of species: ");
    if((sindex>=0)&&(sindex<(signed)stringLine.size()))
    {
      nspec = atoi(&line[19]);
      break;
    }
  }

  f.close();
};

void System::ReadD2K(const char*filename)
{
  char line[MAXSTR];
  int nspec;

  ifstream f(filename);
  if ( f.fail() )
  {
    cout<<__FILE__<<__LINE__<<"Error opening file"<<endl;
    return;
  };

  species.clear();
  nspec = 0;
  while (f.getline(line,MAXSTR))
  {
    int sindex = -1;
    string stringLine=line;
    sindex=stringLine.find(" NUMBER OF TYPES OF QUANTUM PARTICLES: ");
    if ((sindex>=0)&&(sindex<(signed)stringLine.size()))
    {
      nspec = atoi(&line[39]);

      species.push_back( new Species(this,ELECTRON) ); 
      for (int ispec=0;ispec<nspec;ispec++)
        species.push_back( new Species(this,PROTON) );  //FIXME: Not always H+
      break;
    }
  }

  f.close();
};
