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
#ifndef AGV_SET_H
#define AGV_SET_H

#include <Parameter.h>

#include <vector>

using namespace std;

//class Integrator;
class Shell;

class Set 
{
  public:

    Set(const char*,const char*);

    bool is_neighbor;
    char basname[MAXSTR];
    char elname[MAXSTR];
    int nco;
    int ll;
    int ul;
    vector<Shell*> shell; 

    //size_t NumberOfShells(void);

    //void GenerateAuxiliary(Set*,Integrator*);
    void PrintOut(void);
    void SetName(const char*);
    void ReadFile(const char*,const char*);
    void WriteFile(const char*);
    int LMax(void);
    //int NumberOfSphericalOrbitals(void);
};

#endif // AGV_SET_H
