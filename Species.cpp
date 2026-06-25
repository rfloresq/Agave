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
#include <fstream>

/*
#include <Array.h>
#include <ECP.h>
#include <Element.h>
#include <Parameter.h>
#include <Symmetry.h>
*/
#include <Atom.h>
#include <Basis.h>
#include <Global.h>
#include <Integrator.h>
#include <Math.h>
#include <Matrix.h>
#include <Set.h>
#include <Shell.h>
#include <Species.h>
#include <Structure.h>
#include <System.h>

Species::Species(System *isys,size_t itype)
{
  sys = isys;
  mol = sys->mol;
  type = itype;
  nspin = 2;
  X = 0;
  H = 0;
  for (int ispin=0;ispin<nspin;ispin++)
  {
    F[ispin] = 0;
    C[ispin] = 0;
    P[ispin] = 0;
  }

  /*
  D = 0;
  U = 0;
  FUKUIL = 0;
  FUKUIR = 0;
  */

  strcpy(basisname,"TZVP");
  basis = new Basis(mol);

  /*
  // Pseudopotentials
  strcpy(ppname,"NONE");
  ecp = 0;
  pps = new Basis(mol);
  */
};

void Species::CopyParameters(Species* other)
{
  strcpy(basisname,other->basisname);
  strcpy(ppname,other->ppname);
  strcpy(guess,other->guess);
  type = other->type;
  sys = other->sys;
  mol = other->mol;
  basis = other->basis;
  pps = other->pps;

  nspin = other->nspin;
  for (int ispin=0;ispin<nspin;ispin++)
    npart[ispin] = other->npart[ispin];
}

void Species::WriteFile(const char* filename,int id)
{
  ofstream f(filename,ios::app);
  f << endl << " Species ID: "<<id<<endl;
  for (int ispin=0;ispin<nspin;ispin++)
    f << endl << " Number of particles for spin "<<ispin
	      << " is " << npart[ispin] << endl;
  f.close();

  basis->WriteFile(filename);
};

void Species::ReadFile(const char*filename,const char*fmt, int q)
{ 
  cout << "Basis for species "<<q<<endl;
  if (basis) 
  {
    delete basis;
    basis = new Basis(mol);
  }
  basis->ReadFile(filename,fmt,q);

  if ( strcasecmp( fmt , "AGV" ) == 0 ) 
  {
    ReadAGV(filename,q);
  } 
  else if ( strcasecmp( fmt , "D2K" ) == 0 ) 
  {
    ReadD2K(filename,q);
  } 

};

double Species::Charge()
{
  double charge = -1.0;

  if (type==ELECTRON)
  {
    charge = -1.0;
  }
  else if ((type==HYDROGEN)||
	   (type==POSITRON)||
	   (type==PROTON))
  {
    charge = 1.0;
  }

  return charge;
}

double Species::Mass()
{
  double mass = 1.0;

  if ((type==ELECTRON)||(type==POSITRON))
  {
    mass = 1.0;
  }
  else if ((type==HYDROGEN)||(type==PROTON))
  {
    mass = 1836.2;
  }

  return mass;
}

void Species::ChangeGuess(const char* newguess)
{
  strcpy(guess,newguess);
}


bool Species::PrepareForSCF()
{
  cout << "Preparing for SCF species "<<endl;
  cout << TotalNumberOfParticles() << " of type "<<type<<endl;
  if (type==HYDROGEN) mol->NullifyHydrogens();

  if (type==ELECTRON) cout << " ELECTRON "<<endl;
  else if (type==POSITRON) cout << " POSITRON "<<endl;
  else if (type==HYDROGEN) cout << " HYDROGEN "<<endl;
  else if (type==PROTON) cout << " PROTON "<<endl;

  if (basis) delete basis;
  basis = new Basis(mol);
  basis->ReadLibrary("/usr/local/BASIS",basisname,type);
  basis->PrintOut();
  int dbas = basis->nco;
  cout << "Number of BASIS functions "<<dbas<<endl;

  if (X) delete X;
  X = new Matrix(dbas,dbas);
  OrthogonalizationMatrix();

  // Spin independent mono-particle matrix
  if (H) delete H;
  H = new Matrix(dbas,dbas);
  OneParticleMatrix(H,noshift,"kinetic");
  OneParticleMatrix(H,noshift,"core");
  
  //SALCs
  //if (U) delete U;
  //U = new Matrix(dbas,dbas);
  //sys->symmetry->BuildSALCs(this);

  // Initialize central force integrator
  /*
  if (!(strncmp(ppname,"NONE",4) == 0)) 
  {
    if (pps) delete pps;
    pps = new Basis(mol);
    ecp = 0;
    if ((strncmp(ppname,"ECP",3) == 0)|| 
        (strncmp(ppname,"RECP",4) == 0)|| 
        (strncmp(ppname,"QECP",4) == 0)) 
    {
      pps->ReadLibrary("/usr/local/ECPS",ppname,type);
      pps->PrintOut();
      ecp = new ECP(this);
      ecp->BuildMatrix(H);
    }
    else
    {
      pps->ReadLibrary("/usr/local/MCPS",ppname,type);
      cout << "PSEUDOPOTENTIAL: "<<ppname<<" not handled"<<endl;
    }
  }
  */

  int i;
  
  dbas = basis->nco;
  for (int ispin=0;ispin<nspin;ispin++)
  {
    nmo[ispin] = dbas;
    if (F[ispin]) delete F[ispin];
    if (C[ispin]) delete C[ispin];
    if (P[ispin]) delete P[ispin];
    F[ispin] = new Matrix(dbas,dbas);
    C[ispin] = new Matrix(dbas,dbas);
    P[ispin] = new Matrix(dbas,dbas);

    for (i=0;i<npart[ispin];i++)
      oo[ispin][i] = 1.0;
    for (i=npart[ispin];i<dbas;i++)
      oo[ispin][i] = 0.0;
  }
  return true;
};

/*
void Species::PrepareForFitting()
{
  for (size_t i=0;i<spin.size();i++)
    spin[i]->PrepareForFitting();

  int dabas;

  dabas = auxis->nco;

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
};

void Species::PrepareForLinearResponse()
{
  for (size_t i=0;i<spin.size();i++)
    spin[i]->PrepareForLinearResponse();
};

bool Species::PrepareForRT()
{
  bool res;

  for (size_t i=0;i<spin.size();i++)
  {
    res = spin[i]->PrepareForRT();
    if (!res) break;
  }

  return res;
};
*/

void Species::OrthogonalizationMatrix()
{
  int i,j;

  // Build overlap matrix
  X->SetZero();
  OneParticleMatrix(X,noshift,"overlap");
  //rfm X->PrintOut("OVERLAP",15,5);

  // Allocate auxiliary fields
  double *values;
  values = new double[basis->nco];

  X->SVDPower(values,-0.5,1.0e-8);
	  /*
  double *vec;
  vec = new double[basis->nco];
  X->Diagonalize(values); 
  
  for (i=0;i<basis->nco;i++)
  {
    if (values[i]<1.0e-14)
    {
      values[i] = 0.0;
      for (j=0;j<basis->nco;j++)
        X->SetValue(j,i,0.0);
    }
    else
    {
      X->GetColValues(i,vec);
      for (j=0;j<basis->nco;j++)
        X->SetValue(j,i,vec[j]/sqrt(values[i]));
    }
  }
  delete vec;
  */
  delete values;
}


void Species::OneParticleMatrix(Matrix *mat, double* shift, const char* option)
{
  int dao,i,iset,ishl,j,jset,jshl,k,lla,llb,lmax;
  double ra[3],rb[3],rc[3];
  double **intblk, **intblk2;
  Shell *sa,*sb;

  lmax = basis->LMax();
  dao = ((lmax+1)*(lmax+2))/2;

  intblk = new double*[dao];
  for (i=0;i<dao;i++)
    intblk[i] = new double[dao];
  if ( strncmp(option,"core",4) == 0 ) 
  {
    intblk2 = new double*[dao];
    for (i=0;i<dao;i++)
      intblk2[i] = new double[dao];
  };

  for (iset=0;iset<(signed)basis->set.size();iset++)
  {
    basis->mol->atom[iset]->Double(ra); 
    for (jset=0;jset<(signed)basis->set.size();jset++)
    {
      rb[0] = basis->mol->atom[jset]->x + shift[0]; 
      rb[1] = basis->mol->atom[jset]->y + shift[1]; 
      rb[2] = basis->mol->atom[jset]->z + shift[2]; 
      for (ishl=0;ishl<(signed)basis->set[iset]->shell.size();ishl++)
      {
        sa = basis->set[iset]->shell[ishl];
        lla = basis->set[iset]->ll + sa->ll;
        for (jshl=0;jshl<(signed)basis->set[jset]->shell.size();jshl++)
        {
          sb = basis->set[jset]->shell[jshl];
          llb = basis->set[jset]->ll + sb->ll;
          if ((jset>iset)||((jset==iset)&&(jshl>=ishl))) 
          {
            if ( strncmp(option,"kinet",5) == 0 )
            {
              sys->gaussi->Kinetic(ra,rb,sa,sb,intblk,0,true);
              // Apply mass 
              for (i=0;i<sa->nco;i++)
                for (j=0;j<sb->nco;j++)
                  intblk[i][j] /= Mass();
            }
            else if ( strncmp(option,"esp",3) == 0 )
            {
	      k = atoi(&option[3]);
              rc[0] = sys->mol->atom[k]->x; 
              rc[1] = sys->mol->atom[k]->y; 
              rc[2] = sys->mol->atom[k]->z; 
              sys->gaussi->Core(ra,rb,rc,sa,sb,intblk,0,0,true);
              for (i=0;i<sa->nco;i++)
                for (j=0;j<sb->nco;j++)
                  intblk[i][j] *= Charge();
            }
            else if ( strncmp(option,"core",4) == 0 )
            {
              for (i=0;i<sa->nco;i++)
                for (j=0;j<sb->nco;j++)
                  intblk[i][j] = 0.0;
              // Coulomb potential due to point charges
              for (k=0;k<sys->mol->Natom();k++)
              {
                if (sys->mol->atom[k]->charge>0.0)
                {
                  rc[0] = sys->mol->atom[k]->x; 
                  rc[1] = sys->mol->atom[k]->y; 
                  rc[2] = sys->mol->atom[k]->z; 
                  sys->gaussi->Core(ra,rb,rc,sa,sb,intblk2,0,0,true);
                  for (i=0;i<sa->nco;i++)
                    for (j=0;j<sb->nco;j++)
                      intblk[i][j] += intblk2[i][j]*
                      sys->mol->atom[k]->charge*Charge();
                }
              }
            }
            else if ( strncmp(option,"dipole",6) == 0 )
            {
              if ( strncmp(option,"dipole x",8) == 0 )
                sys->gaussi->Dipole(ra,rb,sa,sb,intblk,0);
              else if ( strncmp(option,"dipole y",8) == 0 )
                sys->gaussi->Dipole(ra,rb,sa,sb,intblk,1);
              else sys->gaussi->Dipole(ra,rb,sa,sb,intblk,2);
            }
            else 
            {
              sys->gaussi->Overlap(ra,rb,sa,sb,intblk,0,true);
            }
            for (i=0;i<sa->nco;i++)
              for (j=0;j<sb->nco;j++)
                mat->ShiftValue(lla+i,llb+j,intblk[i][j]);
          }
        }
      }
    }
  }

  for (i=0;i<dao;i++)
    delete[] intblk[i];
  delete[] intblk;
  if ( strncmp(option,"core",4) == 0 ) 
  {
    for (i=0;i<dao;i++)
      delete[] intblk2[i];
    delete[] intblk2;
  }

  mat->Symmetrize();
}

void Species::FockFromCore()
{
  for (int ispin=0;ispin<nspin;ispin++)
    F[ispin]->Copy(H);  
}

double Species::Energy(bool electric)
{
  double energy;

  energy = 0.0;
  for (int ispin=0;ispin<nspin;ispin++)
  {
    energy += 0.5*P[ispin]->QTrace(F[ispin]);
    energy += 0.5*P[ispin]->QTrace(H);
    //if (electric) energy += 0.5*P[ispin]->QTrace(D);
  }

  return energy;
}

/*
void Species::Dipole(double *dipole)
{
  int dbas,i;
  size_t ispin;
  double charge;
  double dt[3];

  for (i=0;i<3;i++)
    dipole[i] = 0.0;

  dbas = basis->nco;
  charge = Charge();
  Matrix M(dbas,dbas);

  M.SetZero();
  OneParticleMatrix(&M,noshift,"dipole x"); 
  for (ispin=0;ispin<spin.size();ispin++)
    dipole[0] += charge*spin[ispin]->P->QTrace(&M);

  M.SetZero();
  OneParticleMatrix(&M,noshift,"dipole y"); 
  for (ispin=0;ispin<spin.size();ispin++)
    dipole[1] += charge*spin[ispin]->P->QTrace(&M);

  M.SetZero();
  OneParticleMatrix(&M,noshift,"dipole z"); 
  for (ispin=0;ispin<spin.size();ispin++)
    dipole[2] += charge*spin[ispin]->P->QTrace(&M);

}
*/

void Species::TwoParticleMatrix(bool doex)
{
  bool skip[2];
  int ib,ibo,is,iso,jb,jbo,js,jso,lla,llao,llb,llbo;
  double energy;
  double ra[3],rao[3],rb[3],rbo[3];
  double blk[MAXNCO][MAXNCO][MAXNCO][MAXNCO];
  Shell *sa,*sao,*sb,*sbo;
  Species *other;
  Matrix *PO;

  for (size_t jspec=0;jspec<sys->species.size();jspec++)
  {
    other = sys->species[jspec];

    PO = new Matrix(other->basis->nco,other->basis->nco);
    other->GetDensityMatrix(PO);

    sys->gaussi->SetInteraction(&GammaCoulomb);

    // Effect of charges: repulsion or atraction
    PO->Scale( Charge()*other->Charge() );

    // Interaction/exchange contributions
    for (ib=0;ib<(signed)basis->set.size();ib++)
    {
      ra[0] = basis->mol->atom[ib]->x; 
      ra[1] = basis->mol->atom[ib]->y; 
      ra[2] = basis->mol->atom[ib]->z; 
      for (jb=0;jb<(signed)basis->set.size();jb++)
      {
        rb[0] = basis->mol->atom[jb]->x; 
        rb[1] = basis->mol->atom[jb]->y; 
        rb[2] = basis->mol->atom[jb]->z; 
        for (is=0;is<(signed)basis->set[ib]->shell.size();is++)
        {
          sa = basis->set[ib]->shell[is];
          lla = basis->set[ib]->ll + sa->ll;
          for (js=0;js<(signed)basis->set[jb]->shell.size();js++)
          {
            sb = basis->set[jb]->shell[js];
            llb = basis->set[jb]->ll + sb->ll;

            for (ibo=0;ibo<(signed)other->basis->set.size();ibo++)
            {
              rao[0] = other->basis->mol->atom[ibo]->x; 
              rao[1] = other->basis->mol->atom[ibo]->y; 
              rao[2] = other->basis->mol->atom[ibo]->z; 
              for (jbo=0;jbo<(signed)other->basis->set.size();jbo++)
              {
                rbo[0] = other->basis->mol->atom[jbo]->x; 
                rbo[1] = other->basis->mol->atom[jbo]->y; 
                rbo[2] = other->basis->mol->atom[jbo]->z; 
                for (iso=0;iso<(signed)other->basis->set[ibo]->shell.size();iso++)
                {
                  sao = other->basis->set[ibo]->shell[iso];
                  llao = other->basis->set[ibo]->ll + sao->ll;
                  for (jso=0;jso<(signed)other->basis->set[jbo]->shell.size();jso++)
                  {
                    sbo = other->basis->set[jbo]->shell[jso];
                    llbo = other->basis->set[jbo]->ll + sbo->ll;

                    sys->gaussi->Interaction(ra,rb,rao,rbo,sa,sb,sao,sbo,blk,TOLNUM,skip,0,0,0);
                    AddInteraction(lla,llb,llao,llbo,sa->nco,sb->nco,sao->nco,sbo->nco,blk,PO);
                    if (doex&&(other==this))
                      AddExchange(lla,llb,llao,llbo,sa->nco,sb->nco,sao->nco,sbo->nco,blk);
                  }
                }
              }
            }
          }
        }
      }
    }
    delete PO;
  }
}

/*
void Species::TwoParticleMatrixDF(bool response)
{
  int i,ib,ibo,is,iso,j,jb,js,k,lla,llao,llb;
  bool skip;
  double ra[3],rb[3],rao[3];
  double blk[MAXNCO][MAXNCO][MAXNCO];
  size_t ispin,jspec;
  Array *xo;
  Shell *sa,*sb,*sao;
  Species *other;

  double cblk[spin.size()][MAXNCO][MAXNCO];

  xo = new Array(auxis->nco);

  sys->gaussi->SetInteraction(&GammaCoulomb);

  for (jspec=0;jspec<sys->species.size();jspec++)
  {
    other = sys->species[jspec];
    other->GetAuxiliaryDensityVector(xo, Charge()*other->Charge(), response );

    // Effect of charges: repulsion or atraction
    for (ib=0;ib<(signed)basis->set.size();ib++)
    {
      basis->mol->atom[ib]->Double(ra); 
      for (is=0;is<(signed)basis->set[ib]->shell.size();is++)
      {
        sa = basis->set[ib]->shell[is];
        lla = basis->set[ib]->ll+sa->ll;
        for (jb=0;jb<(signed)basis->set.size();jb++)
        {
          basis->mol->atom[jb]->Double(rb); 
          for (js=0;js<(signed)basis->set[jb]->shell.size();js++)
          {
            sb = basis->set[jb]->shell[js];
            llb = basis->set[jb]->ll+sb->ll;

            for (ispin=0;ispin<spin.size();ispin++)
              for (i=0;i<sa->nco;i++)
                for (j=0;j<sb->nco;j++)
                  cblk[ispin][i][j] = 0.0;

            for (ibo=0;ibo<(signed)auxis->set.size();ibo++)
            {
              auxis->mol->atom[ibo]->Double(rao); 
              for (iso=0;iso<(signed)auxis->set[ibo]->shell.size();iso++)
              {
                sao = auxis->set[ibo]->shell[iso];
                llao = auxis->set[ibo]->ll+sao->ll;

                sys->gaussi->Interaction3(ra,rb,rao,sa,sb,sao,blk,skip,0,0);

                for (ispin=0;ispin<spin.size();ispin++)
                  for (i=0;i<sa->nco;i++)
                    for (j=0;j<sb->nco;j++)
                      for (k=0;k<sao->nco;k++)
		        if (other==this)
			{
		          if (response)
                            cblk[ispin][i][j] += (xo->Value(llao+k)+
		  			       spin[ispin]->z1->Value(llao+k))*
				               blk[i][j][k];
		          else
                            cblk[ispin][i][j] += (xo->Value(llao+k)+
				    	       spin[ispin]->z->Value(llao+k))*
				               blk[i][j][k];
			}
		        else
                          cblk[ispin][i][j] += xo->Value(llao+k)*blk[i][j][k];
              }
            }
            for (ispin=0;ispin<spin.size();ispin++)
              for (i=0;i<sa->nco;i++)
                for (j=0;j<sb->nco;j++)
		  if (response)
                    spin[ispin]->F1->ShiftValue(lla+i,llb+j,cblk[ispin][i][j]);
	          else
                    spin[ispin]->F->ShiftValue(lla+i,llb+j,cblk[ispin][i][j]);
          }
        }
      }
    }
  }
  delete xo;
}
*/

void Species::GetDensityMatrix(Matrix *PO)
{
  PO->SetZero();
  for (int ispin=0;ispin<nspin;ispin++)
    PO->Add(P[ispin]);
}

/*
void Species::GetPulayMatrix(Matrix *W)
{
  W->SetZero();
  for (size_t ispin=0;ispin<spin.size();ispin++)
    spin[ispin]->BuildPulayMatrix(W);
}

void Species::GetAuxiliaryDensityVector(Array *x, double sf, bool response)
{
  size_t ispin;

  x->SetZero();

  if (response)
    for (ispin=0;ispin<spin.size();ispin++)
      x->Add( spin[ispin]->x1 );
  else
    for (ispin=0;ispin<spin.size();ispin++)
      if (spin[ispin]->number_of_particles>0)
        x->Add( spin[ispin]->x );

  x->Scale( sf );
}
*/

void Species::NewMOs()
{
  int dbas = basis->nco;

  Matrix T(dbas,dbas);

  for (int ispin=0;ispin<nspin;ispin++)
  {
    // F' = X F X
    T.Copy(F[ispin]);
    TransformMatrixBasis(&T,true,false);

    // Diagonalization
    T.Diagonalize(&moe[ispin][0]); 

    // Back transformation of coefficients
    X->Transpose();
    T.Multiply(X,C[ispin],false);
    X->Transpose();
  }
}

void Species::NewDensity(double s)
{
  int i,mu,nu,dbas;
  double *moc;
  
  dbas = basis->nco;
  moc = new double[dbas];

  for (int ispin=0;ispin<nspin;ispin++)
  {
    if (type==ELECTRON) P[ispin]->Scale( 1.0 - s ); 
    else P[ispin]->Scale( 0.0 ); 

    for (i=0;i<=npart[ispin];i++)
    {
      C[ispin]->GetColValues( i , moc );
      for (mu=0;mu<dbas;mu++)
        for (nu=mu;nu<dbas;nu++)
          if (type==ELECTRON) P[ispin]->ShiftValue(mu,nu,oo[ispin][i]*moc[mu]*moc[nu]*s);
          else P[ispin]->ShiftValue(mu,nu,oo[ispin][i]*moc[mu]*moc[nu]);
    }

    // Symmetrize
    P[ispin]->Symmetrize();
  }

  delete moc;
}

/*
void Species::EvaluateFittingCoefficients()
{
  for (size_t ispin=0;ispin<spin.size();ispin++)
    spin[ispin]->EvaluateFittingCoefficients();
}

void Species::EvaluateFittingCoefficientsY()
{
  for (size_t ispin=0;ispin<spin.size();ispin++)
    spin[ispin]->EvaluateFittingCoefficientsY();
}
*/

void Species::ChangeNumber( double n, int spinid ) 
{
  npart[spinid] = n;
} 

void Species::ChangeBasis(char *newbas)
{
  strcpy(basisname,newbas);
}

/*
void Species::ChangePPs(char *newpp)
{
  strcpy(ppname,newpp);
}
*/

int Species::TotalNumberOfParticles()
{
  int np;
  np= 0;
  for (int ispin=0;ispin<nspin;ispin++)
    np += npart[ispin];

  return np;
}

/*
void Species::EvaluateElectricFieldMatrix(double* e)
{
  int dbas = basis->nco;

  Matrix DT(dbas,dbas);

  if (D) delete D;
  D = new Matrix(dbas,dbas);

  DT.SetZero();
  OneParticleMatrix(&DT,noshift,"dipole x"); 
  DT.Scale(e[0]);
  D->Add(&DT);

  DT.SetZero();
  OneParticleMatrix(&DT,noshift,"dipole y"); 
  DT.Scale(e[1]);
  D->Add(&DT);

  DT.SetZero();
  OneParticleMatrix(&DT,noshift,"dipole z"); 
  DT.Scale(e[2]);
  D->Add(&DT);
}

void Species::TwoParticleIIMOMatrix(Species* other,int spina,int spinao,int p,int q, Matrix eri)
{
  bool skip[2];
  int dbas,dbaso;
  int i,ib,ibo,is,iso,j,jb,jbo,js,jso,k,l,lla,llao,llb,llbo;
  double ss;
  double ra[3],rao[3],rb[3],rbo[3];
  double blk[MAXNCO][MAXNCO][MAXNCO][MAXNCO];
  double cblk[MAXNCO][MAXNCO][MAXNCO];
  Shell *sa,*sao,*sb,*sbo;

  dbas = basis->nco;
  dbaso = other->basis->nco;

  double wp[dbas];
  double wq[dbaso];

  spin[spina]->C->GetColValues(p,wp);
  other->spin[spinao]->C->GetColValues(q,wq);

  // Set as coulombic integrals
  sys->gaussi->SetInteraction(&GammaCoulomb);

  // Initialize
  eri.SetZero();

  for (ibo=0;ibo<(signed)other->basis->set.size();ibo++)
  {
    other->basis->mol->atom[ibo]->Double(rao); 
    for (jbo=0;jbo<(signed)other->basis->set.size();jbo++)
    {
      other->basis->mol->atom[jbo]->Double(rbo); 
      for (iso=0;iso<(signed)other->basis->set[ibo]->shell.size();iso++)
      {
        sao = other->basis->set[ibo]->shell[iso];
        llao = other->basis->set[ibo]->ll + sao->ll;
        for (jso=0;jso<(signed)other->basis->set[jbo]->shell.size();jso++)
        {
          sbo = other->basis->set[jbo]->shell[jso];
          llbo = other->basis->set[jbo]->ll + sbo->ll;

          for (ib=0;ib<(signed)basis->set.size();ib++)
          {
            basis->mol->atom[ib]->Double(ra); 
            for (is=0;is<(signed)basis->set[ib]->shell.size();is++)
            {
              sa = basis->set[ib]->shell[is];
              lla = basis->set[ib]->ll + sa->ll;

              // ( p mu | nu sigma )
              for (i=0;i<sao->nco;i++)
                for (j=0;j<sbo->nco;j++)
                  for (k=0;k<sa->nco;k++)
                    cblk[i][j][k] = 0.0;

              for (jb=0;jb<(signed)basis->set.size();jb++)
              {
                basis->mol->atom[jb]->Double(rb); 
                for (js=0;js<(signed)basis->set[jb]->shell.size();js++)
                {
                  sb = basis->set[jb]->shell[js];
                  llb = basis->set[jb]->ll + sb->ll;

                  sys->gaussi->Interaction(rao,rbo,ra,rb,sao,sbo,sa,sb,blk,TOLNUM,skip,0,0,0);
                  for (i=0;i<sao->nco;i++)
                    for (j=0;j<sbo->nco;j++)
                      for (k=0;k<sa->nco;k++)
                      {
                        ss = 0.0;
                        for (l=0;l<sb->nco;l++)
                          ss += wp[llb+l]*blk[i][j][k][l];
                          
                        cblk[i][j][k] += ss;
                      }
                }
              }

              // ( p mu | nu t )
              for (k=0;k<sa->nco;k++)
                for (j=0;j<sbo->nco;j++)
                {
                  ss = 0.0;
                  for (i=0;i<sao->nco;i++)
                    ss += wq[llao+i]*cblk[i][j][k];
                  eri(lla+k,llbo+j) += ss;
                }

            }
          }

        }
      }
    }
  }
  // ( p q | s t ) = < p s | q t >
  Matrix tmp(dbas,dbaso);
  eri.Multiply(other->spin[spinao]->C,&tmp,true);
  spin[spina]->C->Transpose();
  tmp.Multiply(spin[spina]->C,&eri,false);
  spin[spina]->C->Transpose();
  
  // Effect of charges: repulsion or atraction
  eri.Scale( Charge()*other->Charge());
}
*/

void Species::TwoParticleIITensor(Species* other,int spina,int spinao,int p,vector<Matrix> eri)
{
  bool skip[2];
  int dbas,dbaso;
  int i,ib,ibo,is,iso,j,jb,jbo,js,jso,k,l,lla,llao,llb,llbo,q;
  double ss;
  double ra[3],rao[3],rb[3],rbo[3];
  double blk[MAXNCO][MAXNCO][MAXNCO][MAXNCO];
  double cblk[MAXNCO][MAXNCO][MAXNCO];
  Shell *sa,*sao,*sb,*sbo;

  dbas = basis->nco;
  dbaso = other->basis->nco;

  double wp[dbas];
  C[spina]->GetColValues(p,wp);

  // Set as coulombic integrals
  sys->gaussi->SetInteraction(&GammaCoulomb);

  // Initialize
  for (size_t ieri=0;ieri<eri.size();ieri++)
    eri[ieri].SetZero();

  for (ibo=0;ibo<(signed)other->basis->set.size();ibo++)
  {
    other->basis->mol->atom[ibo]->Double(rao); 
    for (jbo=0;jbo<(signed)other->basis->set.size();jbo++)
    {
      other->basis->mol->atom[jbo]->Double(rbo); 
      for (iso=0;iso<(signed)other->basis->set[ibo]->shell.size();iso++)
      {
        sao = other->basis->set[ibo]->shell[iso];
        llao = other->basis->set[ibo]->ll + sao->ll;
        for (jso=0;jso<(signed)other->basis->set[jbo]->shell.size();jso++)
        {
          sbo = other->basis->set[jbo]->shell[jso];
          llbo = other->basis->set[jbo]->ll + sbo->ll;

          for (ib=0;ib<(signed)basis->set.size();ib++)
          {
            basis->mol->atom[ib]->Double(ra); 
            for (is=0;is<(signed)basis->set[ib]->shell.size();is++)
            {
              sa = basis->set[ib]->shell[is];
              lla = basis->set[ib]->ll + sa->ll;

              // ( p mu | nu sigma )
              for (i=0;i<sao->nco;i++)
                for (j=0;j<sbo->nco;j++)
                  for (k=0;k<sa->nco;k++)
                    cblk[i][j][k] = 0.0;

              for (jb=0;jb<(signed)basis->set.size();jb++)
              {
                basis->mol->atom[jb]->Double(rb); 
                for (js=0;js<(signed)basis->set[jb]->shell.size();js++)
                {
                  sb = basis->set[jb]->shell[js];
                  llb = basis->set[jb]->ll + sb->ll;

                  sys->gaussi->Interaction(rao,rbo,ra,rb,sao,sbo,sa,sb,blk,TOLNUM,skip,0,0,0);
                  for (i=0;i<sao->nco;i++)
                    for (j=0;j<sbo->nco;j++)
                      for (k=0;k<sa->nco;k++)
                      {
                        ss = 0.0;
                        for (l=0;l<sb->nco;l++)
                          ss += wp[llb+l]*blk[i][j][k][l];
                          
                        cblk[i][j][k] += ss;
                      }
                }
              }

              // ( p mu | nu t )
              for (q=0;q<dbaso;q++)
                for (k=0;k<sa->nco;k++)
                  for (j=0;j<sbo->nco;j++)
                  {
                    ss = 0.0;
                    for (i=0;i<sao->nco;i++)
                      ss += cblk[i][j][k]*
                            (*other->C[spinao])(llao+i,q);
                    eri[q](lla+k,llbo+j) += ss;
                  }

            }
          }

        }
      }
    }
  }
  // ( p q | s t ) = < p s | q t >
  Matrix tmp(dbas,dbaso);
  for (q=0;q<dbaso;q++)
  {
    eri[q].Multiply(other->C[spinao],&tmp,true);
    C[spina]->Transpose();
    tmp.Multiply(C[spina],&eri[q],false);
    C[spina]->Transpose();
    // Effect of charges: repulsion or atraction
    eri[q].Scale( Charge()*other->Charge());
  }
}

/*
void Species::FittingMatrix()
{
  BuildFittingMatrix(this,G);

  // G^-1/2
  double *values;
  values = new double[auxis->nco];
  
  GI->Copy(G);
  GI->SVDPower(values,double(-1.0),double(1.0e-8));

  GIS->Copy(G);
  GIS->SVDPower(values,double(-0.5),double(1.0e-7));

  delete values;
}

void Species::BuildFittingMatrix(Species* other,Matrix* G)
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
      for (jb=0;jb<(signed)other->auxis->set.size();jb++)
      {
        rb[0] = other->auxis->mol->atom[jb]->x; 
        rb[1] = other->auxis->mol->atom[jb]->y; 
        rb[2] = other->auxis->mol->atom[jb]->z; 
        for (js=0;js<(signed)other->auxis->set[jb]->shell.size();js++)
        {
          sb = other->auxis->set[jb]->shell[js];
          llb = other->auxis->set[jb]->ll+sb->ll;

          sys->gaussi->Interaction2(ra,rb,sa,sb,blk,0);
          for (i=0;i<sa->nco;i++)
            for (j=0;j<sb->nco;j++)
              G->SetValue(lla+i,llb+j,blk[i][j]);
        }
      }
    }
  }
}

void Species::SolveGxb(Array* x)
{
  int dabas;
  Array *tmp;

  dabas = auxis->nco;
  tmp = new Array(dabas);

  tmp->Copy(x);
  tmp->Transform(GI,x);

  delete tmp;
}
*/

void Species::TransformMatrixBasis(Matrix* M,bool ortho,bool density)
{
  int drow = M->NCol();
  int dcol = M->NRow();

  Matrix T(drow,dcol);
  Matrix W(drow,dcol);

  // Select transformation matrix
  if (ortho!=density)
  {
    W.Copy(X);
  }
  else
  {
    T.SetZero();
    double noshift[3] = {0.0,0.0,0.0};
    OneParticleMatrix(&T,noshift,"overlap");
    T.Multiply(X,&W,true);
  }

  // Apply transformation
  W.Transpose();
  M->Multiply(&W,&T,true);
  W.Transpose();
  T.Multiply(&W,M,false);
}

// Add up four center integrals to Fock matrix (Interaction)
void Species::AddInteraction(int ll1,int ll2,int ll3,int ll4,int n1,int n2,int n3,int n4, double b[MAXNCO][MAXNCO][MAXNCO][MAXNCO], Matrix *PO)
{
  int i,ispin,j,k,l;
  double p;

  for (k=0;k<n3;k++)
    for (l=0;l<n4;l++)
    {
      p = (*PO)(ll3+k,ll4+l);
      for (i=0;i<n1;i++)
        for (j=0;j<n2;j++)
          for (ispin=0;ispin<nspin;ispin++)
            F[ispin]->ShiftValue(ll1+i,ll2+j,b[i][j][k][l]*p);
    }
}

void Species::AddExchange(int ll1,int ll2,int ll3,int ll4,int n1,int n2,int n3,int n4, double b[MAXNCO][MAXNCO][MAXNCO][MAXNCO])
{
  int i,ispin,j,k,l;
  double p;

  for (k=0;k<n3;k++)
    for (l=0;l<n4;l++)
      for (i=0;i<n1;i++)
        for (j=0;j<n2;j++)
          for (ispin=0;ispin<nspin;ispin++)
          {
            p = -(*P[ispin])(ll3+k,ll2+j);
            F[ispin]->ShiftValue(ll1+i,ll4+l,b[i][j][k][l]*p);
          }
}

void Species::ReadAGV(const char*filename, int q)
{ 
}

void Species::ReadD2K(const char*filename, int q)
{ 
  ifstream f(filename);
  if ( f.fail() )
  {
    cout<<__FILE__<<__LINE__<<"Error opening file"<<endl;
    return;
  };

  char line[MAXSTR];
  char str[MAXSTR];
  double rval;
  int eindex,nblk,ll,qindex,ul;
  int dbas,ibas,ijump,jbas,n,njump;

  // Orbitals
  dbas = basis->nco;
  for (int ispin=0;ispin<nspin;ispin++)
  {
    if (C[ispin]) delete C[ispin];
    C[ispin] = new Matrix(dbas,dbas);
    nmo[ispin] = dbas;
  }
  
  while (f.getline(line,MAXSTR))
  {
    string stringLine=line;
    eindex=stringLine.find(" MO COEFFICIENTS");
    qindex=stringLine.find(" COEFFICIENTS FOR QUANTUM PARTICLE");
    nblk = 0;
    if(q==0&&(eindex>=0)&&(eindex<(signed)stringLine.size()))
    {
      nblk = 5;
      njump = 3;
    }
    else if(q>0&&(qindex>=0)&&(qindex<(signed)stringLine.size()))
    {
      nblk = 6;
      njump = 0;
    }
    if (nblk==0) continue;
    ll = 0;
    ul = MIN(ll + nblk-1,dbas-1);
    while(ul>=ll)
    {
      for (jbas = ll; jbas <= ul; jbas++)
      {
        f>>n;
        if (n-jbas!=1) break;
      }
      if (q==0)
      {
        // Energies
        for (jbas = ll; jbas <= ul; jbas++)
        {
          f>>str;
	}
        // Occupations
        for (jbas = ll; jbas <= ul; jbas++)
        {
          f>>str;
	}
      }
      //ul = jbas-1;
      for (ibas=0;ibas<dbas;ibas++)
      {
        f>>n;
	for (int ijump=0;ijump<njump;ijump++)
   	  f>>str;
	/*
        if (n-ibas!=1) 
	{
          cout << "UPS: something went wrong on reading"<<endl;
          cout << str << endl;
          cout << "n= "<<n<<" ibas = "<<ibas<<endl;
          exit(0);
        }
	*/
        for (jbas = ll; jbas <= ul; jbas++) 
        {
          f >> rval;
          C[0]->SetValue(ibas,jbas,rval); // FIXME: Loading only alphas
        }
      }
      ll = ul+1;
      ul = MIN(ll + nblk-1,dbas-1);
    }
  }

  // Electronic Fukui functions
  /*
  if (q==0)
  {
    if (FUKUIL) delete FUKUIL;
    if (FUKUIR) delete FUKUIR;
    FUKUIL = new Matrix(basis->nco,basis->nco);
    FUKUIR = new Matrix(basis->nco,basis->nco);
    int multiplicity;
    bool fukui;
    fukui = false;
    while (f.getline(line,MAXSTR))
    {
      bool left;
      int sindex = -1;
      int rindex = -1;
      string stringLine=line;
      sindex=stringLine.find(" NUMBER OF ELECTRONS:");
      if((sindex>=0)&&(sindex<(signed)stringLine.size()))
      {
        int npart = atoi(&line[21]);
	int nb =  (npart + 1 - multiplicity)/2;
	int na = npart - nb;
        spin[ALPHA]->ChangeNumber(na);
        spin[BETA]->ChangeNumber(nb);
      }
      sindex=stringLine.find(" MULTIPLICITY:");
      if((sindex>=0)&&(sindex<(signed)stringLine.size()))
      {
        multiplicity = atoi(&line[14]);
      }
      sindex=stringLine.find(" FUKUI^- MATRIX");
      if (sindex<0) sindex=stringLine.find(" FUKUI^+ MATRIX");
      if((sindex>=0)&&(sindex<(signed)stringLine.size()))
      {
        fukui = true;
        rindex=stringLine.find(" FUKUI^- MATRIX");
        if((rindex>=0)&&(rindex<(signed)stringLine.size()))
        { left = true; } else { left = false; }

        int nblk = 4;
        ll = 0;
        ul = MIN(ll + nblk-1,basis->nco-1);
        while(ul>=ll)
        {
          for (jbas = ll; jbas <= ul; jbas++)
          {
            f>>n;
            if (n-jbas!=1) break;
          }
          ul = jbas-1;
          for (ibas=0;ibas<basis->nco;ibas++)
          {
            f>>n;
   	    f>>str;
	    f>>str;
	    f>>str;
            if (n-ibas!=1) 
	    {
              cout << "UPS: something went wrong on reading"<<endl;
              cout << str << endl;
              cout << "n= "<<n<<" ibas = "<<ibas<<endl;
              exit(0);
            }
            for (jbas = ll; jbas <= ul; jbas++) 
            {
              f >> rval;
              if (left) FUKUIL->SetValue(ibas,jbas,rval);
	      else FUKUIR->SetValue(ibas,jbas,rval);
            }
          }
          ll = ul+1;
          ul = MIN(ll + nblk-1,basis->nco-1);
        }
      }
    }
    if (!fukui)
    {
      delete FUKUIL;
      delete FUKUIR;
    }
  }
*/


  //for (size_t i=0;i<spin.size();i++)
  //  spin[i]->ReadFile(basis,filename,fmt,i);

  f.close();
};

