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
#ifndef AGV_GRID_H
#define AGV_GRID_H

#include <vector>

#include <Vector.h>
#include <Matrix.h>

using namespace std;

class Grid 
{
  public:
    Grid(void);

    void GetLattice(double*);
    void SetLattice(double*,double);
    void GetCubePoints(int,int,int,Vector r[8]);
    void GetRowOfPoints(int,int,int,double*,double*,double*);

    int npoint[4];
    Vector dstep[4];
    Vector vertex[4];
};

#endif  // AGV_GRID_H
