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
#ifndef AGV_MATRIX_H
#define AGV_MATRIX_H

class Matrix
{
  public:
    Matrix(int,int);

    double *values;

    void operator += ( Matrix );
    double & operator () ( int, int );

    double AbsoluteSum(void);
    double Sum(void);
    double QTrace(Matrix*);
    double Trace(void);

    int NCol(void);
    int NRow(void);

    void Add(Matrix*);
    void Copy(Matrix*);
    void Copy(double**);
    void Diagonalize(double*);
    void ExtractFrom(Matrix*,int,int);
    void Fill(double);
    void FillRange(double,int,int,int,int);
    void GetColValues( int, double* );
    void GetRowValues( int, double* );
    void InsertTo(Matrix*,int,int);
    void Multiply(Matrix*,Matrix*,bool);
    void Print(const char*,const char*,int,int);
    void PrintOut(const char*,int,int);
    void Scale(double);
    void SetValue( int, int , double );
    void SetZero(void);
    void ShiftedAdd(Matrix*,int,int);
    void ShiftValue( int, int , double );
    void SVDPower(double*,double,double);
    void SwitchRows(int,int);
    void SwitchColumns(int,int);
    void Symmetrize(void);
    void Transform(Matrix*);
    void Transpose(void);
    void VectorMultiply(double*,double*);

  protected:

    int dim_col;
    int dim_row;
};

#endif // AGV_MATRIX_H
