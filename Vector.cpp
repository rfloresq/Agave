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
#include <math.h> 
#include <stdio.h>

#include <Global.h>
#include <Math.h>
#include <Matrix.h>
#include <Vector.h>

// Constructors
Vector::Vector( double *v )
{
  x = v[0];
  y = v[1];
  z = v[2];
};

Vector::Vector( double xx , double yy , double zz ) 
{
  x = xx;
  y = yy;
  z = zz;
};

Vector::Vector( Vector* v ) 
{
  x = (*v)[0];
  y = (*v)[1];
  z = (*v)[2];
};

Vector::~Vector()
{
}

void Vector::SetZero() 
{
  x = 0.0;
  y = 0.0;
  z = 0.0;
};

void Vector::Copy( Vector* v ) 
{
  x = (*v)[0];
  y = (*v)[1];
  z = (*v)[2];
};

void Vector::SetValue(int i,double v) 
{
  if (i==0) x = v;
  else if (i==1) y = v;
  else if (i==2) z = v;
};

// Access
double Vector::operator [] ( short n )
{
  if ( n == 0 )
  {
    return x;
  }
  else if ( n == 1 )
  {
    return y;
  }
  else if ( n == 2 )
  {
    return z;
  } 
  else
  {
    return 0.0;
  };
};

// this = result 
void Vector::operator += ( const Vector &cv )
{
  Vector v = cv;
  x += v[0];
  y += v[1];
  z += v[2];
};

void Vector::operator -= ( const Vector &cv )
{
  Vector v = cv; 
  x -= v[0];
  y -= v[1];
  z -= v[2];
};

void Vector::operator *= ( double s )
{
  x *= s;
  y *= s;
  z *= s;
};
void Vector::operator *= ( int n )
{
  x *= double(n);
  y *= double(n);
  z *= double(n);
};

void Vector::operator /= ( double s )
{
  x /= s;
  y /= s;
  z /= s;
};

void Vector::operator <= ( const Vector &cv )
{
  double m[3][3];
  Vector v = cv;

  m[0][0] = 0.0;   m[0][1] = -v[2]; m[0][2] = v[1];
  m[1][0] = v[2];  m[1][1] = 0.0;   m[1][2] = -v[0];
  m[2][0] = -v[1]; m[2][1] = v[0];  m[2][2] = 0.0;

  operator %= ( m );
};

void Vector::operator >= ( const Vector &v )
{
  operator <= ( v );
  operator *= ( -1.0 );
};

void Vector::operator %= ( double m[3][3] )
{
  double nx,ny,nz;

  nx = Dot( Vector( m[0][0] , m[0][1] , m[0][2] ) ); 
  ny = Dot( Vector( m[1][0] , m[1][1] , m[1][2] ) );
  nz = Dot( Vector( m[2][0] , m[2][1] , m[2][2] ) );

  x = nx;
  y = ny;
  z = nz;
};

void Vector::operator ^= ( double s ) 
{
  if ( Norm() > 0.0 )
  {
    operator *= ( s/Norm() );
  };
};


// this != result

Vector Vector::operator + ( const Vector &v )
{
  Vector u( x , y , z );
  u += v;
  return u;
};

Vector Vector::operator - ( const Vector &v )
{
  Vector u( x , y , z );
  u -= v;
  return u;
};

Vector Vector::operator < ( const Vector &v )
{
  Vector u( x , y , z );
  u <= v;
  return u;
};

Vector Vector::operator > ( const Vector &v )
{
  Vector u( x , y , z );
  u >= v;
  return u;
};

Vector Vector::operator % ( double m[3][3] )
{
  Vector u( x , y , z );
  u %= m;
  return u;
};

Vector Vector::operator ^ ( double s ) 
{
  Vector u( x , y , z );
  u ^= s;
  return u;
};

Vector Vector::operator * ( double s ) 
{
  Vector u( x , y , z );
  u *= s;
  return u;
};

Vector Vector::operator * ( int n ) 
{
  Vector u( x , y , z );
  u *= n;
  return u;
};

Vector Vector::operator / ( double s ) 
{
  Vector u( x , y , z );
  u /= s;
  return u;
};

// functions

double Vector::Dot( const Vector &cv )
{
  Vector v = cv;
  return( x*v[0] + y*v[1] + z*v[2] );
};

double Vector::Norm() 
{
  return  sqrt( Dot( Vector( x , y , z ) ) );
};

void Vector::Rotate( const Vector &v , double angle )
{
  Matrix m(3,3);
  double t[3][3];
  Vector u(v.x,v.y,v.z);
  u.GetRotationMatrix(angle,&m);
  t[0][0] = m(0,0);
  t[0][1] = m(0,1);
  t[0][2] = m(0,2);
  t[1][0] = m(1,0);
  t[1][1] = m(1,1);
  t[1][2] = m(1,2);
  t[2][0] = m(2,0);
  t[2][1] = m(2,1);
  t[2][2] = m(2,2);
  operator %= ( t );
};

void Vector::GetRotationMatrix(double angle, Matrix* m)
{
  double cosp,sinp,cost;

  Vector u(x,y,z);
  if ( u.Norm() == 0.0 ) 
  {
    return;
  }
  u ^= 1.0;

  angle = angle*M_PI/180.0; // Input angle in degrees

  cosp = cos(angle);
  sinp = sqrt(1.0-cosp*cosp);//sin(angle); 
  cost = 1.0 - cosp;

  m->SetValue(0,0,u[0]*u[0]*cost + cosp);
  m->SetValue(0,1,u[0]*u[1]*cost - u[2]*sinp);
  m->SetValue(0,2,u[0]*u[2]*cost + u[1]*sinp);
  m->SetValue(1,0,u[0]*u[1]*cost + u[2]*sinp);
  m->SetValue(1,1,u[1]*u[1]*cost + cosp);
  m->SetValue(1,2,u[1]*u[2]*cost - u[0]*sinp);
  m->SetValue(2,0,u[0]*u[2]*cost - u[1]*sinp);
  m->SetValue(2,1,u[1]*u[2]*cost + u[0]*sinp);
  m->SetValue(2,2,u[2]*u[2]*cost + cosp);
};

// Build 3D matrix representation for a single operation
void Vector::GetReflectionMatrix(Matrix* m)
{
  int i,j;

  Vector u(x,y,z);
  if ( u.Norm() == 0.0 ) 
  {
    m->SetZero();
    return;
  }
  u ^= 1.0;

  for (i=0;i<3;i++)
  {
    m->SetValue(i,i,1.0);
    for (j=0;j<3;j++)
    {
      m->ShiftValue(i,j,-2.0*u[i]*u[j]);
    }
  }
}

void Vector::Double(double *v)
{
  v[0] = x;
  v[1] = y;
  v[2] = z;
}

// Transform Cartesian to Spherical coordinates. 
// History: - Creation (19.10.05, RFM)
//          - Translated to C++ (21.07.15, RFM)
/*
void Vector::SphericalCoordinates(double *r, double *theta, double *phi)
{
  // Radius 
  *r = Norm();

  // Elevation angle 
  if (*r<TOLNUM) *theta = 0.0;
  else *theta = acos(z/(*r));

  // Azimuth angle 
  if (ABS(x)<TOLNUM) 
  {
    if (ABS(y)<TOLNUM) *phi = 0.0;
    else if (y<0.0) *phi = 1.5*M_PI;
    else *phi = 0.5*M_PI;
  }
  else if (x>0.0) *phi = atan(y/x);
  else *phi = atan(y/x) + M_PI;
}
*/

bool Vector::Colineal(const Vector &v)
{
  Vector u = v;

  if ( u.Norm() == 0.0 ) 
  {
    return true;
  }
  u ^= 1.0;
  u *= -Dot(u); 
  u.x += x;
  u.y += y;
  u.z += z;
  if (u.Norm()<TOLNUM) return true;
  else return false;
}

Vector Vector::Orthogonal() 
{
  double s;
  Vector u,v;

  if (Norm()<TOLNUM) return Vector( 1.0 , 0.0 , 0.0 );

  u = Vector( y , x , z );
  if (Colineal(u)) u = Vector( x , z , y );
  v = Vector( x , y , z );
  v ^= 1.0;
  s = u.Dot(v);
  v *= s;
  u -= v;
  u ^= 1.0;

  return u;
}

