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
#ifndef AGV_MATH_H
#define AGV_MATH_H

#include <math.h>

#ifndef M_PI
  #define M_PI acos(-1.0)
#endif

#define ABS(a) ((a)>=0.0?(a):-(a))

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

#define RAD(a) ((a)*M_PI/180.0)
#define DEG(a) ((a)*180.0/M_PI)

#define MOD(a,b) (((int)(a))-((int)(b))*(((int)(a))/((int)(b)))) 

extern double Factorial(int);
extern double DoubleFactorial(int);
extern double Noverk(int,int);

#endif  // AGV_MATH_H
