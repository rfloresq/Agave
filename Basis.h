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

#ifndef AGV_BASIS_H
#define AGV_BASIS_H

#include <vector>

using namespace std;

//class Integrator;
class Set;
class Structure;

class Basis
{
  public:

    Basis(Structure*);

    Structure *mol;
    vector<Set*> set;           // One-center basis set
    int nco;

    //size_t NumberOfShells(void);

    //void GenerateAuxiliary(Basis*,Integrator*);
    void PrintOut(void);
    void ReadFile(const char*,const char*,int);
    void WriteFile(const char*);
    void LoadFromLibrary(const char*,const char*,int);
    void ReadLibrary(const char*,const char*,int);
    int LMax(void);
    //int NumberOfSphericalOrbitals(void);
};

#endif // AGV_BASIS_H
