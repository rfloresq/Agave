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
#include <Units.h>

// === Length ===

// Conversion factor taken from ?
double HartreeToeV(double a) { return ((a)*27.21138466); };
double eVToHartree(double a) { return ((a)/27.21138466); };

// Conversion factor 
double AngstromToBohr(double a) { return ((a)*1.8897261349309467); }; 
double BohrToAngstrom(double a) { return ((a)/1.8897261349309467); };

// === Angles ===
double DegreeToRadian(double deg) { return M_PI*deg/180.0; }; 
double RadianToDegree(double rad) { return 180.0*rad/M_PI; };

// === Time ===
//
double AUToAttosecond(double tau) 
{ 
  // Permeability of vacuum [N/A**2] 
  double c = 2.99792458e+8;
  double h = 6.626069311e-34;
  double e = 1.6021765314e-19;
  double m = 9.109382616e-31;
  double mu = 4.0*M_PI*1.0e-7;
  double a = 0.5*mu*c*e*e/h; 
  double r = 0.5*m*c*a*a/h;
  double f = 4.0*M_PI*c*r*1.0e-15;

  return (tau/f)*1.0e+3; 
}; 

double AttosecondToAU(double tatto) 
{ 
  double factor = 1.0/AUToAttosecond(1.0);
  return tatto*factor;
}

