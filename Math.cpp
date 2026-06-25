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

#include <Math.h>

double Factorial(int n)
{
  double r;

  if ( n < 0 )
  {
    r = 0.0;
  }
  else
  {
    r = 1.0;
    for ( int i = 1; i <= n; i++ )
    {
      r *= (double)i;
    }
  }
  return r;
}

double DoubleFactorial(int n)
{
  int i;
  double r;

  if ( n <= 1 )
  {
    r = 1.0;
  }
  else if ( MOD(n,2) == 1 )
  {
    r = 1.0;
    for ( i = 1 ; i <= n ; i += 2 )
    {
      r *= (double)i;
    }
  }
  else
  {
    r = 1.0;
    for ( i = 2 ; i <= n ; i += 2 )
    {
      r *= (double)i;
    }
  }

  return r;
}

double Noverk(int n,int k)
{
  double r;

  if (k>=0&&k<=n)
  {
    r = Factorial(n)/(Factorial(n-k)*Factorial(k));
  }
  else
  {
    r = 0.0;
  }
  return r;
}


