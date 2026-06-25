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
// Master object for Quantum Chemistry calculations
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <fstream>

/*
#include <ADFT.h>
#include <ADPT.h>
#include <Basis.h>
#include <LFT.h>
*/
#include <Global.h>
#include <QM.h>
#include <SCF.h>
#include <Species.h>
#include <System.h>
/*
#include <RTADFT.h>
#include <RSGF.h>
#include <TDDFT.h>
#include <SSQ.h>
#include <Time.h>
*/

using namespace std;

QM::QM(System *isys)
{
  sys = isys;
  scf = new SCF(this);
  /*
  ssq = new SSQ(this);
  rsgf = new RSGF(this);
  lft = new LFT(this);
  tddft = new TDDFT(this);
  adft = 0;
  adpt = 0;
  rtadft = 0;
  setup = false;
  */

  species.clear();
  for (size_t itype=0;itype<NTYPES;itype++)
    species.push_back( new Species(sys,itype) ); 

  ChangeGuess("Core",ELECTRON);
  ChangeGuess("Core",POSITRON);
  ChangeGuess("Core",HYDROGEN);
}

void QM::PrepareForSCF()
{

  setup = false;

  if (sys->PrepareForSCF(species))
  {
	  /*
    if (( strcasecmp( scf->method , "Hartree-Fock" ) == 0 ) ||
        ( strcasecmp( scf->method , "Hartree-None" ) == 0 ))
    {
      if (adft) delete adft;
    }
    else 
    {
      if (!adft) adft = new ADFT(this);
      adft->SetFunctional(scf->method);
    }

    if (adft&&(sys->auxis_unique->nco==0))
    {
      cout << "Cannot run ADFT without auxiliary functions"<<endl;
      return;
    }
	  */

    scf->Setup();
    scf->Guess();
    setup = true;
  } 
}

/*
void QM::PrepareForADPT()
{
  Time *time = new Time();
  time->Start();

  if (!adft)
  {
    cout << "ADPT needs ADFT"<<endl;
    return;
  }

  sys->PrepareForLinearResponse();

  if (adpt) delete adpt;
  adpt = new ADPT(this);
  adpt->PrepareLinearResponse();

  cout << "TIME (QM::PrepareForADPT): "<<time->Seconds()<<endl;
  delete time;
}

void QM::PrepareForRT()
{
  if (!adft)
  {
    cout << "RT-ADFT needs ADFT"<<endl;
    return;
  }

  if (sys->PrepareForRT())
  {
    if (rtadft) delete rtadft;
    rtadft = new RTADFT(this);

    rtadft->Setup();
  }
  else
  {
    cout << "RT-ADFT setup failed"<<endl;
    return;
  }
}

void QM::PrepareForRSGF()
{
  setup = false;
  if (sys->PrepareForSCF(species)) setup = true;
}

void QM::PrepareForLFT()
{
  setup = false;
  if (sys->PrepareForSCF(species)) setup = true;
}

void QM::PrepareForTDDFT()
{
  setup = true;
}
*/

void QM::ChangeBasis(char *newbas, int type)
{
  for (size_t ispec=0;ispec<species.size();ispec++)
    if (species[ispec]->type == type ) 
      species[ispec]->ChangeBasis( newbas );
}

/*
void QM::ChangePPs(char *newpp, int type)
{
  for (size_t ispec=0;ispec<species.size();ispec++)
    if (species[ispec]->type == type ) 
      species[ispec]->ChangePPs( newpp );
}
*/

void QM::ChangeGuess(const char *newguess, int type)
{
  for (size_t ispec=0;ispec<species.size();ispec++)
    if (species[ispec]->type == type ) 
      species[ispec]->ChangeGuess( newguess );
}

/*
void QM::ReadFile(const char*filename)
{
} 

void QM::WriteFile(const char*filename)
{
} 

void QM::Polarizability(Matrix alpha)
{
  if (!adpt)
  {
    cout << "Prepare ADPT first"<<endl;
    return;
  }
  alpha.SetZero();
  adpt->Polarizability(&alpha);
};

void QM::ParticleBindingEnergy(int typea,int spina,int orbnuma,double *ie)
{
  if (!adpt)
  {
    cout << "Prepare ADPT first"<<endl;
    return;
  }

  Time *time = new Time();
  time->Start();

  adpt->ParticleBindingEnergy(typea,spina,orbnuma,ie);
  cout << "TIME (Binding energy from ADPT): "<<time->Seconds()<<endl;
  delete time;
};

void QM::NuclearFukuiFunction(int typea,int spina,int orbnuma)
{
  if (!adpt)
  {
    cout << "Prepare ADPT first"<<endl;
    return;
  }

  Time *time = new Time();
  time->Start();

  adpt->NuclearFukuiFunction(typea,spina,orbnuma);
  cout << "TIME (NNF from ADPT): "<<time->Seconds()<<endl;
  delete time;
};

void QM::HRemoval(double *he,int atomid,int n)
{
  if ((!adpt)&&(n>1))
  {
    cout << "Prepare ADPT first"<<endl;
    return;
  }

  Time *time = new Time();
  time->Start();

  adpt->HRemoval(he,atomid,n);
  cout << "TIME (H removal from ADPT): "<<time->Seconds()<<endl;
  delete time;
};
*/


void QM::ChangeNumber( double n, int type, int spin ) 
{
  for (size_t ispec=0;ispec<species.size();ispec++)
    if (species[ispec]->type==type)
      species[ispec]->ChangeNumber(n,spin); 
} 

