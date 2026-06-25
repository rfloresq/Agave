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
// Master object for Quantum Chemistry calculations
#ifndef AGV_QM_H
#define AGV_QM_H

#include <vector>

#include <Parameter.h>
#include <Matrix.h>

using namespace std;

/*
class ADFT;
class ADPT;
class LFT;
class Matrix;
class RTADFT;
class RSGF;
class TDDFT;
*/
class SCF;
class Species;
class System;
/*
class SSQ;
*/

class QM
{
  public:

    QM(System*);

    System* sys;
    /*
    ADFT* adft;
    ADPT* adpt;
    RSGF *rsgf;
    TDDFT *tddft;
    LFT *lft;
    RTADFT* rtadft;
    */
    SCF *scf;
    //SSQ *ssq;
    vector<Species*> species;
    bool setup;
    /*

    void HRemoval(double*,int,int);
    void NuclearFukuiFunction(int,int,int);
    void ParticleBindingEnergy(int,int,int,double*);
    void Polarizability(Matrix);
    */
    void PrepareForSCF(void);
    /*
    void PrepareForADPT(void);
    void PrepareForRT(void);
    void PrepareForRSGF(void);
    void PrepareForLFT(void);
    void PrepareForTDDFT(void);
    */
    void ChangeBasis(char*,int);
    //void ChangePPs(char*,int);
    void ChangeNumber(double,int,int);
    void ChangeGuess(const char*,int);
    /*
    void ReadFile(const char*);
    void WriteFile(const char*);
    */
};

#endif // AGV_QM_H
