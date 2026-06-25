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

#ifndef AGV_STRUCTURE_H
#define AGV_STRUCTURE_H

#include <vector>

#include <Vector.h>

using namespace std;

class Atom;

class Structure
{
  public:
    Structure(void);
    Structure(Structure*);

    vector<Atom*> atom;

    Atom* GetAtom(int);

    //double Mass(void);
    double NuclearRepulsionEnergy(void);

    int Natom(void);
    int NumberOfElectronsForNeutral(void);

    void AddAtom(int,int,int,int,double,double,double);
    //void CenterMass(void);
    void C2Z(int); 
    //void DeleteAtom(int);
    void DeleteLastAtom(void);
    //void NuclearDipole(double*);
    void NullifyHydrogens(void);
    void ReadFile(const char*,const char*);
    void ReadD2K(const char*);
    //void ReadPIE(const char*);
    void ReadXYZ(const char*);
    //void ReadZMatrix(const char*);
    void WriteFile(const char*);
    /*
    void WriteXYZ(const char*);
    void WriteZMatrix(const char*);
    */
    void Z2C(int); 

  protected:

    double Dihedral(int,int,int,int);

    void Center(void);
};

#endif // AGV_STRUCTURE_H
