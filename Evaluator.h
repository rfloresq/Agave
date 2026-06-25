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

#ifndef AGV_EVALUATOR_H
#define AGV_EVALUATOR_H

class Basis;
class Matrix;
class Set;
class Species;
class Spin;
class System;

/*
extern void GetGTO(Basis*,int,double*,double*,double*,double*, int);
*/
extern void GetBasis(Basis*,double,double,double,double*);
/*
extern void GetOrbital(Spin*,Basis*,int,double*,double*,double*,double*,int);
extern void GetDensityLike(Matrix*,Matrix*,Basis*,double*,double*,double*,
            double*,int,const char*);
extern void GetChargeDensity(System*,double*,double*,double*,double*,int,bool);
extern void GetParticleDensity(Species*,double*,double*,double*,double*,
	    int,double,bool);
extern void GetParticleElectrostaticPotential(Species*,double*,double*,double*,double*,int);
extern void GetFukui(System*,double*,double*,double*,double*,int,char*);
extern void GetASC(System*,double*,double*,double*,double*,int);
*/


extern void GetGTODerivatives(Basis*,int,double*,double*,double*,
            double*,double*,double*,int);
extern void GetBasisDerivatives(Basis*,double,double,double,
            double*,double*,double*);
extern void GetOrbitalDerivatives(Species*,int,int,double*,double*,double*, double*,double*,double*,int);
extern void GetDensityLikeDerivatives(Matrix*,Matrix*,Basis*,double*,
            double*,double*,double*,double*,double*,int,const char*);
/*
extern void GetChargeDensityDerivatives(System*,double*,double*,double*,double*,
            double*,double*,int,bool);
	    */
extern void GetParticleDensityDerivatives(Species*,double*,double*,double*, double*,double*,double*,int,double,bool);
/*
extern void GetFukuiDerivatives(System*,double*,double*,double*,double*,
            double*,double*,int,char*);
	    */

#endif  // AGV_EVALUATOR_H
