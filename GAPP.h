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
#ifndef AGV_GAPP_H
#define AGV_GAPP_H

#include <time.h>   

#include <vector>

#include <Matrix.h>   
#include <Parameter.h>   

using namespace std;

class Array;
class System;

class GAPP
{
  public:
    GAPP(System*);
    
    System* sys;

    int Map2HP(int,int,int,int,int);
    int Map2PH(int,int,int,int,int);
    int NumberOfHoles(void);
    int NumberOfParticles(void);

    void Decode(int*,int*,int*,int*,const char*);
    void Diagonal(const char*,int,int,int,int,int,int,double*);
    void GetEigenvalues(double*);
    void NewtonRaphson(double*,double*,int,int,int,int,int,int,int*);
    void SelfEnergy(double,double*,int,int,int,int,int,int);
    void SelfEnergy1PP(double,double*,int,int,int);
    //void SelfEnergy2PP(double,double*,int,int,int,int,int,int);
	 
  protected:

    char method[MAXSTR];

    /*
    void ConstantSelfEnergy(double*,int,int,int,int,int,int);
    void ConstantSelfEnergy1PP(double*,int,int,int);
    void ConstantSelfEnergy2PP(double*,int,int,int,int,int,int);
    */
    void SOERI(int,int,int,vector<Matrix>);
    void FFF(double,Array*,Array*,Array*,int);
    void FVA(int,int,int,Array*,Array*,int*,bool);
    void FVA2(int,int,int,Array*,bool);
    void FVF(Array*,Array*,bool);
    /*
    void AVA3(double*,int,int,int);
    void GVB(int,int,int,int,int,int,Array*,Array*);
    void AAVAA1(double*,int,int,int,int,int,int);
    */
};

#endif // AGV_GAPP_H
