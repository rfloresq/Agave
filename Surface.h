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

#ifndef AGV_SURFACE_H
#define AGV_SURFACE_H

#include <vector>

using namespace std;

class System;
class Basis;
class Grid;
class Set;
class Matrix;

class Surface 
{
  public:
    Surface(System*,Grid*);

    System* sys;
    Grid* grid;

    char name[MAXSTR];
    int type;
    int style;
    int shiny;
    int species;
    int spin;
    int orbital_number;
    int point_size;
    int line_width;
    int nscalar;
    double iso;
    double color[4];
    unsigned long list;
    double* scalar;

    void Build(void);
    void GetCubeValues(int,int,int,double*);
    void GetDerivatives(double*,double*,double*,double*,double*,double*,int);
    void GetValues(double*,double*,double*,double*,int);
    //void GetExtrema(double*,double*);
    void Setup(void);

  protected:

    //void BuildChargeDensity(bool);
    void BuildParticleDensity(int,double,bool);
    void BuildDensityLike(Matrix*,Matrix*,Basis*,const char*);
    //void BuildFukui(void);
    void BuildGTO(void);
    void BuildOrbital(void);
};

#endif  // AGV_SURFACE_H

