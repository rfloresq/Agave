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
#ifndef AGV_SYSTEM_H
#define AGV_SYSTEM_H

#include <vector>

//#include <Parameter.h>
//#include <Matrix.h>

using namespace std;

/*
class Array;
class Basis;
class Matrix;
class Symmetry;
*/

class Integrator;
class Species;
class Structure;

class System 
{
  public:

    System(void);

    Structure* mol;
    /*
    char fname[MAXSTR];
    Basis *auxis_unique;          
    Symmetry *symmetry;
    */
    Integrator *gaussi;
    /*
    Matrix *G;
    Matrix *GI;
    Matrix *GIS;
    vector<double> fx,fy,fz,fv;
    */
    vector<Species*> species;

    bool OnlyElectrons(void);
    //bool PrepareForRT(void);
    bool PrepareForSCF(vector<Species*>);

    double Energy(bool);

    /*
    //void EvaluateFittingCoefficients(void);
    void BuildFittingMatrix(Matrix*,Basis*);
    void FittingMatrix(void);
    void SolveGxb(Array*);
    void Dipole(double*);
    void PrepareForLinearResponse(void);
    */
    void ReadFile(const char*,const char*);
    void ReadAGV(const char*);
    void ReadD2K(const char*);
    //void ReadASC(const char *);
    void WriteFile(const char*,const char*);
    //void InvertSuper(Matrix,int);

};

#endif // AGV_SYSTEM_H
