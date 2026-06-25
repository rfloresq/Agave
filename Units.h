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
#ifndef AGV_UNITS_H
#define AGV_UNITS_H

// === Energy ===
extern double HartreeToeV(double);
extern double eVToHartree(double);

// === Length ===
extern double AngstromToBohr(double); 
extern double BohrToAngstrom(double); 

// === Angles ===
extern double DegreeToRadian(double);
extern double RadianToDegree(double);

// === Time ===
extern double AUToAttosecond(double); 
extern double AttosecondToAU(double); 

#endif  // AGV_UNITS_H
