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
#ifndef AGV_SPECIES_H
#define AGV_SPECIES_H

#include <vector>

#include <Parameter.h>

using namespace std;

/*
class Array;
class ECP;
class Integrator;
*/
class Basis;
class Matrix;
class Structure;
class System;

class Species
{
  public:

    Species(System*,size_t);

    char basisname[MAXSTR];
    char ppname[MAXSTR];
    char guess[MAXSTR];
    size_t type;
    System *sys;
    Structure *mol;
    Basis *basis,*auxis,*pps;          
    /*
    Matrix *G;
    Matrix *GI;
    Matrix *GIS;
    */
    Matrix *X,*H,*D,*U;
    Matrix *F[MAXSPIN],*C[MAXSPIN],*P[MAXSPIN];
    int nmo[MAXSPIN]; 
    double oo[MAXSPIN][MAXSTO];   // Orbital occupations
    /*
    Matrix *FUKUIL,*FUKUIR;
    */
    int nspin;
    int npart[MAXSPIN];

    //bool PrepareForRT(void);
    bool PrepareForSCF(void);

    /*
    void BuildFittingMatrix(Species*,Matrix*);
    void FittingMatrix(void);
    void SolveGxb(Array*);
    */

    void ChangeBasis(char*);
    /*
    void ChangePPs(char*);
    */
    void ChangeGuess(const char*);
    void ChangeNumber(double,int); 
    void CopyParameters(Species*);
    /*
    void Dipole(double*);
    void EvaluateElectricFieldMatrix(double*);
    void EvaluateFittingCoefficients(void);
    void EvaluateFittingCoefficientsY(void);
    */
    void FockFromCore(void);
    void GetDensityMatrix(Matrix*);
    /*
    void GetPulayMatrix(Matrix*);
    void GetAuxiliaryDensityVector(Array*,double,bool);
    */
    void NewDensity(double);
    void NewMOs(void);
    /*
    void PrepareForFitting(void);
    void PrepareForLinearResponse(void);
    */
    void ReadFile(const char*,const char*,int);
    void ReadAGV(const char*,int);
    void ReadD2K(const char*,int);
    void WriteFile(const char*,int);
    void OrthogonalizationMatrix(void);
    void OneParticleMatrix(Matrix*,double*,const char*);
    void TwoParticleMatrix(bool);
    /*
    void TwoParticleMatrixDF(bool);
    void TwoParticleIIMOMatrix(Species*,int,int,int,int,Matrix);
    */
    void TwoParticleIITensor(Species*,int,int,int,vector<Matrix>);

    double Charge(void);
    double Energy(bool);
    double Mass(void);

    int TotalNumberOfParticles(void);
    double moe[MAXSPIN][MAXSTO];      // Orbital energies

  protected:

    double noshift[3] = {0.0,0.0,0.0};
    // ECP *ecp;

    void AddInteraction(int,int,int,int,int,int,int,int,double b[MAXNCO][MAXNCO][MAXNCO][MAXNCO],Matrix*);
    void AddExchange(int,int,int,int,int,int,int,int,double b[MAXNCO][MAXNCO][MAXNCO][MAXNCO]);
    void TransformMatrixBasis(Matrix*,bool,bool);
};

#endif // AGV_SPECIES_H
