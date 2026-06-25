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

#ifndef AGV_ATOM_H
#define AGV_ATOM_H

#include <Parameter.h>
#include <Vector.h>

//class AM;

class Atom : public Vector
{
  public:

    Atom(char*,double,double,double);
    Atom(int,double,double,double);

    char label[MAXSTR];
    int atomic_number;
    int ref_bond;
    int ref_angle;
    int ref_dihedral;
    double charge;
    double bond;
    double angle;
    double dihedral;
    /*
    double nuclear_magnetic_moment[3];
    AM *nuclear_spin;
    */

    double Distance(Atom*);
    double GetCovalentRadius(void);
    double GetVDWRadius(void);

    double* GetColor(void);

    void Setup(void);
    void SetSymbol(char*);
    void SetCharge(double);
    //void SetAtomicNumber(int);

  protected:

    int SymbolToAtomicNumber( char* );
};

#endif // AGV_ATOM_H
