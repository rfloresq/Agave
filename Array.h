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
#ifndef AGV_ARRAY_H
#define AGV_ARRAY_H

class Matrix;

class Array
{
  public:
    Array(int);

    double *values;

    double Sum(void);
    double SumProduct(Array*);
    double SumProduct(Array*,int);
    double SumProduct(double*);
    double Value(int);

    void Add(Array*);
    void Add(double*);
    void Fill(double);
    void Scale(double);
    void SetZero(void);
    void Transform(Matrix*,Array*);
    void SetValue(int,double);
    void ShiftValue(int,double);
    void Copy(Array*);

  protected:

    int dim;
};

#endif // AGV_ARRAY_H
