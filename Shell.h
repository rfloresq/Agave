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

#ifndef AGV_SHELL_H
#define AGV_SHELL_H

#define PRIMARY      0     // Basis (electronic, protonic, ... any)
#define AUXILIARY    1     // Auxiliary basis set (electronic, ... )
#define LOCAL_PPS    2     // Effective/Model/Solvent/... potentials
#define SEMILOCAL    3     // Effective core potentials
#define NON_LOCAL    4     // Model core potential projector

#include <vector>

//class Integrator;

using namespace std;

class Shell 
{
  public:

    Shell(void);

    int n;
    int l;
    int nco;
    int ll;   // Lower limit of its CAOs in basis
    int ul;   // Upper limit of its CAOs in basis
    int type;
    vector<double> z;
    vector<double> d;
    vector<double> ncsto;

    //double RadialPotential(double);

    void EvaluateGTO(double*,double*,double*,double*,int,int);
    void EvaluateGTODerivative(double*,double*,double*,double*,double*,double*,int, int);
    //void EvaluateRadius(void);
    void Normalize(void);
    //void NormalizeInteraction(Integrator*);
    void PrintOut(void);
    void WriteFile(const char*);

    //double Radius(void);

  protected:

    bool is_normalized;
    double radius;
};

#endif // AGV_SHELL_H
