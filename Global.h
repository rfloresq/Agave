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

#ifndef AGV_GLOBAL_H
#define AGV_GLOBAL_H

// Labels (species)
#define ELECTRON  0
#define POSITRON  1
#define HYDROGEN  2
#define PROTON    3
//#define METAL     3
//#define OXYGEN    4
#define NTYPES    4

// Labels (spins)
#define ALPHA 0
#define BETA  1

// Second quantization elementary operators
#define ANIHILATOR 0
#define CREATOR    1

// Numerical tolerance
#define MINRHO  1.0e-10
#define TOLNUM  1.0e-14

#endif // AGV_GLOBAL_H
