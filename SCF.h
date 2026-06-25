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
#ifndef AGV_SCF_H
#define AGV_SCF_H

#include <fstream>

#include <Parameter.h>

using namespace std;

class QM;
class System;

class SCF 
{
  public:
    SCF(QM*);

    QM *qm;
    System *sys;
    char method[MAXSTR];

    void ChangeMethod(const char*);
    //void Dipole(double*);
    void Guess(void);
    void Setup(void);
    //void SetTolerance(double);
    void Step(void);
    /*
    void NewDensity(void);
    void SeekConvergence(const char*);

    bool Converged(void);
    */

    double Energy(void);
    double Error(void);

  protected:

    double mixing;
    double energy;
    double error;
    double tolerance;
    int cycle;
};

#endif // AGV_SCF_H
