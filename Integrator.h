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
#ifndef AGV_INTEGRATOR_H
#define AGV_INTEGRATOR_H

#define PSDAOPointer(l) ((l)*(l)*(l)-(l))/6 + (l)*(l) + 2*(l) + 1;

#include <Parameter.h>

// +2 for kinetic energy integrals and avoid shifting full program
extern int gaop[MAXLI+3][MAXLI+3][MAXLI+3];
extern int raop[MAXLI+3][MAXLI+3][MAXLI+3];

class Shell;
class Set;

extern double ftab[2*MAXLI+6+1][MAXTABGAM+1];
extern void GammaCoulomb(int,double,double,double*);

class Integrator
{
  public:
    Integrator(void);

    void Overlap(double*,double*,Shell*,Shell*,double**,int,bool);
    void Kinetic(double*,double*,Shell*,Shell*,double**,int,bool);
    void Dipole(double*,double*,Shell*,Shell*,double**,int);
    void Core(double*,double*,double*,Shell*,Shell*,double**,int,int,bool);
    void Interaction(double*,double*,double*,double*,Shell*,Shell*,Shell*,Shell*,
              double [MAXNCO][MAXNCO][MAXNCO][MAXNCO], double,bool*,int,int,int);
    void Interaction3(double*,double*,double*,Shell*,Shell*,Shell*,
	      double [MAXNCO][MAXNCO][MAXNCO],bool,int,int);
    void Interaction2(double*,double*,Shell*,Shell*,double [MAXNCO][MAXNCO],int);
    void SetInteraction( void (*)(int,double,double,double*) );

  protected:

    void HRR(int,int,double*,double**);
    void NormalizeTwoIndex(Shell*,Shell*,double**,double**,bool);
    void (*Gamma)(int,double,double,double*);
};


#endif // AGV_INTEGRATOR_H

