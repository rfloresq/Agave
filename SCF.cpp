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
#include <string.h>

#include <iostream>
#include <iomanip>

//#include <ADFT.h>   
#include <Math.h>
#include <QM.h>
#include <SCF.h>   
#include <Species.h>
#include <Structure.h>
#include <System.h>
//#include <Time.h>

using namespace std;

SCF::SCF(QM *iqm)
{
  qm = iqm;
  sys = qm->sys;

  Setup();
}

void SCF::Setup()
{
  cycle = 0;
  error = 99999.0;
  energy = 99999.0;
  tolerance = 1.0e-8;
  mixing = 1.0; //0.3;
}

void SCF::Guess()
{
  // H
  for (size_t ispec=0;ispec<sys->species.size();ispec++)
    sys->species[ispec]->FockFromCore();
  
  // C
  for (size_t ispec=0;ispec<sys->species.size();ispec++)
    sys->species[ispec]->NewMOs();
}

void SCF::Step()
{
  double penergy;

  // Ignore step request if QM is not properly set
  if (!qm->setup) return;

  //Time *time = new Time();
  //time->Start();

  cycle++;

  // P
  for (size_t ispec=0;ispec<sys->species.size();ispec++)
    if (cycle==1)
      sys->species[ispec]->NewDensity(1.0);
    else 
      sys->species[ispec]->NewDensity(mixing); 

/*
  // y
  for (size_t ispec=0;ispec<sys->species.size();ispec++)
    sys->species[ispec]->EvaluateFittingCoefficientsY(); 

  // x
  for (size_t ispec=0;ispec<sys->species.size();ispec++)
    sys->species[ispec]->EvaluateFittingCoefficients(); 

  // x
  //rfm sys->EvaluateFittingCoefficients(); 
  */

  // H
  for (size_t ispec=0;ispec<sys->species.size();ispec++)
    sys->species[ispec]->FockFromCore();

  // XC energy and potential
  //RFM if (qm->adft) qm->adft->PotentialAndEnergy(); 

  // F
  for (size_t ispec=0;ispec<sys->species.size();ispec++)
  {
    //if (qm->adft) sys->species[ispec]->TwoParticleMatrixDF(false); // Coulomb + xc
    //else if ( strcasecmp( method , "Hartree-None" ) == 0 )
    if ( strcasecmp( method , "Hartree-None" ) == 0 )
      sys->species[ispec]->TwoParticleMatrix(false); // Coulomb + exchange
    else sys->species[ispec]->TwoParticleMatrix(true); // Coulomb + exchange
  }

  penergy = energy;
  energy = Energy();
  error = penergy - energy;
  cout << "Energy: "<<cycle<<" "<<fixed<<setprecision(8)<<energy<<endl;
  
  // C
  for (size_t ispec=0;ispec<sys->species.size();ispec++)
    sys->species[ispec]->NewMOs();

  //cout << "TIME (SCF Step): "<<time->Seconds()<<endl;
  // delete time;
}

/*
void SCF::SetTolerance(double tol)
{
  tolerance = tol;
}

bool SCF::Converged()
{
  if ((ABS(error)<tolerance)&&(cycle>2)) return true;
  else return false;
}

void SCF::SeekConvergence(const char* filename)
{
  ofstream f(filename,ios::app);

  f <<endl<< " === SCF === "<<endl<<endl;

  Guess();

  int ncyc = 0;
  while (!Converged())
  {
    Step();
    f << setw(4)<< ncyc << " SCF energy = "
      <<fixed<<setw(15)<<setprecision(6)<<energy<<" Error = "
      <<fixed<<setw(15)<<setprecision(8)<<error<<endl;
    ncyc++;
    if (ncyc>1000) break;
  } 
  if (Converged()) f <<endl<< " === SCF CONVERGED === "<<endl<<endl;
  else f <<endl<< " === SCF NOT CONVERGED === "<<endl<<endl;

  f.close();
}
*/

double SCF::Energy()
{
  double energy = sys->Energy(false);
  //if (qm->adft) energy += qm->adft->exc;
  return energy;
}

double SCF::Error()
{
  return error;
}


void SCF::ChangeMethod(const char *new_method)
{
  strcpy(method,new_method);
}

