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

#ifndef AGV_ELEMENT_H
#define AGV_ELEMENT_H

#include <Parameter.h>

extern const char*  ELEMENT_SYMBOL[MAXATNUM+1];
extern const double ELEMENT_COV_R[MAXATNUM+1];
extern const double ELEMENT_VDW_R[MAXATNUM+1];
extern const double ELEMENT_BS_R[MAXATNUM+1];
extern const double ELEMENT_MASS[MAXATNUM+1];
extern const char*  ELEMENT_NAME[MAXATNUM+1];
extern double ELEMENT_COLOR[MAXATNUM+1][3];

int SymbolToAtomicNumber(char*);

#endif // AGV_ELEMENT_H
