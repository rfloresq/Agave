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

#ifndef AGV_VECTOR_H
#define AGV_VECTOR_H

class Matrix;

class Vector
{
public:
  Vector( double xx = 0.0 ,
          double yy = 0.0 ,
          double zz = 0.0 );
  Vector( double* );
  Vector(Vector*); 
  virtual ~Vector();

  double x;
  double y;
  double z; 

  bool Colineal(const Vector&);

  double Dot( const Vector & );
  double Norm(void);

  Vector Orthogonal(void); 
  Vector operator + ( const Vector & );
  Vector operator - ( const Vector & );
  Vector operator * ( double );
  Vector operator * ( int );
  Vector operator / ( double );
  Vector operator % ( double m[3][3] );
  Vector operator ^ ( double );
  Vector operator < ( const Vector & );
  Vector operator > ( const Vector & );

  virtual double operator [] ( short n );

  void Copy(Vector*); 
  void Double(double *);
  void GetReflectionMatrix(Matrix*);
  void GetRotationMatrix(double,Matrix*);
  void Rotate( const Vector & , double );
  void SetValue(int,double);
  void SetZero(void);
  void SphericalCoordinates(double*,double*,double*);

  void operator += ( const Vector & );
  void operator -= ( const Vector & );
  void operator >= ( const Vector & );
  void operator <= ( const Vector & );
  void operator *= ( double );
  void operator *= ( int );
  void operator /= ( double );
  void operator %= ( double m[3][3] );
  void operator ^= ( double );
};

#endif  // AGV_VECTOR_H

