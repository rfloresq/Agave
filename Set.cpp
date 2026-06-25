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
#include <stdlib.h>
#include <string.h>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include <Parameter.h>
#include <Math.h>
#include <Set.h>
#include <Shell.h>

using namespace std;

Set::Set(const char* ielname,const char* ibasname)
{
  shell.clear();
  nco = 0;
  is_neighbor = true;
  strcpy(elname,ielname);
  strcpy(basname,ibasname);
};

void Set::SetName(const char* ibasname)
{
  strcpy(basname,ibasname);
}

void Set::ReadFile(const char* basisfile,const char* fmt)
{
  char line[MAXSTR];

  bool labelecps;
  string stringLine=basisfile;
  int sindex=stringLine.find("ECPS");
  if(sindex>=0) labelecps = true;
  else labelecps = false;

  if ( strcasecmp( fmt , "D2K" ) == 0 )
  {
    cout << " LOADING BASIS "<<basname<<" FROM : "<<basisfile<<endl;
    ifstream f(basisfile);
    if ( f.fail() )
    {
      cout << "Error opening file [" << basisfile<< "]"<<endl;
      return;
    };
    string basisname = basname;
    string flag;
    flag = "-";
    flag += elname;
    while(f.getline(line,MAXSTR))
    {
      stringLine=line;
      sindex=stringLine.find(flag);
      if(sindex==1)
      {
        int tindex=stringLine.find(basisname);
        if (tindex>=0)
        {
          shell.clear();
          int k,n,iblk,nblk;
          double d,z;
          char jump[MAXSTR];
          f >> jump;
          if ( jump[0]=='#') 
          {
            f.getline(line,MAXSTR);
            f >> nblk;
          }
          else nblk = atoi(jump);
          int cp;
          for (iblk=0;iblk<nblk;iblk++)
          {
            Shell sb = Shell();
            f >> sb.n;
            f >> sb.l;
            f >> k;
            for (n=0;n<k;n++)
            {
              f >> z;
              f >> d;
              sb.z.push_back(z);
              sb.d.push_back(d);
            }
	    if (labelecps)
	    {
	      if (iblk+1==nblk) sb.type = LOCAL_PPS;
	      else sb.type = SEMILOCAL;
	    }
            shell.push_back( new Shell( sb ) );
	    shell[shell.size()-1]->Normalize();
            cp = shell.size() - 1;
            if (cp==0) shell[cp]->ll = 0;
            else shell[cp]->ll = shell[cp-1]->ul + 1;
            shell[cp]->ul = shell[cp]->ll + shell[cp]->nco - 1;
          }
          nco = shell[cp]->ul + 1;
          break;
        }
      }
    }
    f.close();
  }
}

void Set::WriteFile(const char* filename)
{
  ofstream f(filename,ios::app);
  f <<endl;
  f << "Atom type: "<<elname<<endl;
  f << "Basis name: "<<basname<<endl;
  f << "Number of shells: "<<shell.size()<<endl;
  f.close();
  for (size_t is=0;is<shell.size();is++)
    shell[is]->WriteFile(filename);
}

void Set::PrintOut()
{
  cout <<endl;
  cout << "Atom type: "<<elname<<endl;
  cout << "Basis name: "<<basname<<endl;
  cout << "Number of shells: "<<shell.size()<<endl;
  for (size_t is=0;is<shell.size();is++)
  {
    if (shell[is]->type==SEMILOCAL) cout << "SEMILOCAL "<<endl;
    shell[is]->PrintOut();
  }
}


int Set::LMax()
{
  int is,l,lmax;

  lmax = 0;
  for (is=0;is<(signed)shell.size();is++)
  {
    l = shell[is]->l;
    if (l>lmax) lmax = l;
  }
  return lmax;
}

/*
int Set::NumberOfSphericalOrbitals()
{
  int is,l;

  int nso = 0;
  for (is=0;is<(signed)shell.size();is++)
  {
    l = shell[is]->l;
    nso += 2*l + 1;
  }
  return nso;
}

size_t Set::NumberOfShells()
{
  return shell.size();
}


// Notice that auxiliary function sets are different to 
// original GEN-A2* for the case of hydrogen atoms
void Set::GenerateAuxiliary(Set *bas, Integrator* gi)
{
  int cp,i,iblock,is,l,lbmax,n,nauxis,nblock;
  double r,zbmax,zbmin,zmax,zmin;
  Shell *sa;

  // Max Z for basis
  lbmax = 0;
  zbmax = 0.0;
  for (is=0;is<(signed)bas->shell.size();is++)
  {
    sa = bas->shell[is];
    if (sa->l>lbmax) lbmax = sa->l;
    for (i=0;i<(signed)sa->z.size();i++)
      zbmax = MAX(zbmax,sa->z[i]);
  }
  zbmin = zbmax;
  for (is=0;is<(signed)bas->shell.size();is++)
  {
    sa = bas->shell[is];
    for (i=0;i<(signed)sa->z.size();i++)
      zbmin = MIN(zbmin,sa->z[i]);
  }

  r = 4.0;
  nblock = 3;
  if ( strncmp(bas->elname,"HYDROGEN",8) == 0 ) 
  {
    r = 5.0;
    nblock = 2;
  }

  // How many z?
  n = int(log(zbmax/zbmin)/log(r) + 0.5);

  // Largest z
  zmin = zbmin;
  zmax= zmin*pow(r,double(n-1));
  zmax= 2.0*zmax;  

  double z = zmax*r;
  Shell *sb;

  shell.clear();

  for (iblock=0;iblock<nblock;iblock++)
  {
    if (iblock==0) nauxis = MAX(1,n/nblock + MOD(n,nblock));
    else nauxis = MAX(1,n/nblock);
    for (i=0;i<nauxis;i++)
    {
      z = z/r;
      if (i==0) z *= (r+1.0)/r;
      for (l=0;l<=2*iblock;l++)
      {
        sb = new Shell();
        sb->d.push_back( 1.0 );
        sb->z.push_back( z );
        sb->n = 0;
        sb->l = l;
        sb->NormalizeInteraction(gi);
        sb->EvaluateRadius();
        shell.push_back ( sb );
        cp = shell.size() - 1;
        if (cp==0) shell[cp]->ll = 0;
        else shell[cp]->ll = shell[cp-1]->ul + 1;
        shell[cp]->ul = shell[cp]->ll + shell[cp]->nco - 1;
      }
      if (i==0) z *= r/(r+1.0);
    }
  }
  nco = shell[cp]->ul + 1;

  // Labeling
  // for (i=0;i<(signed)shell.size();i++)
  //  shell[i]->type = AUXILIARY_BASIS_SHELL;
}

*/
