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

#ifndef AGV_PARAMETER_H
#define AGV_PARAMETER_H

/*
#define MAXATOM  1000
*/
#define MAXATNUM  103
#define MAXLDER     0
#define MAXLBAS     4
/*
#define MAXLECPS    4
*/
#define MAXL     MAXLBAS + MAXLDER
#define MAXLI    MAXLBAS + MAXL
//#define MAXNSPEC    3
#define MAXSTO   1000
#define MAXSTR    256
#define MAXNCO   ((MAXL+1)*(MAXL+2))/2 
#define MAXSPIN       2
#define MAXTABGAM 120

#endif // AGV_PARAMETER_H
