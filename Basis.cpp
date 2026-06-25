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
#include <string.h>

#include <iostream>
#include <fstream>
#include <vector>

#include <Atom.h>
#include <Basis.h>
#include <Element.h>
#include <Global.h>
#include <Set.h>
#include <Structure.h>

using namespace std;

Basis::Basis(Structure* is)
{
  mol = is;

  set.clear();
};

void Basis::WriteFile(const char* filename)
{
  ofstream f(filename,ios::app);
  f << endl << " BASIS "<<endl;
  f.close();
  for (size_t i=0;i<set.size();i++)
    set[i]->WriteFile(filename);
};

void Basis::PrintOut()
{
  cout << endl << " BASIS "<<endl;
  for (size_t i=0;i<set.size();i++)
    set[i]->PrintOut();
};

int Basis::LMax()
{
  int l,lmax;

  lmax = 0;
  for (size_t ib=0;ib<set.size();ib++)
  {
    l = set[ib]->LMax();
    if (l>lmax) lmax = l;
  }
  return lmax;
}

// Load primary basis
void Basis::ReadFile(const char* filename,const char* fmt,int qp)
{
  int iatom,natom,iqp;
  char line[MAXSTR];

  natom = mol->Natom();

  if (natom<1)
    cout << "WARNING: Loading basis for ZERO atoms"<<endl;

  iqp = -1;
  if ( strcasecmp( fmt , "D2K" ) == 0 )
  {
    ifstream f(filename);
    if ( f.fail() )
    {
      cout <<__FILE__<<__LINE__<<"Error opening file"<<endl;
      return;
    };
    while(f.getline(line,MAXSTR))
    {
      string stringLine=line;
      int sindex=stringLine.find(" GENERAL BASIS SET");
      if((sindex>=0)&&(sindex<(signed)stringLine.size()))
      {
	iqp++;
	if (iqp!=qp) continue;
        int natom = mol->Natom();
        for (iatom = 0 ; iatom < natom ; iatom++)
          if (qp==0) set.push_back( new Set(ELEMENT_NAME[mol->atom[iatom]->atomic_number],&line[20]));
	  else set.push_back( new Set(ELEMENT_NAME[mol->atom[iatom]->atomic_number],&line[37]));
        f.getline(line,MAXSTR);
        stringLine=line;
        sindex=stringLine.find("BASIS");
        while((sindex>=0)&&(sindex<(signed)stringLine.size()))
        {
          int n;
          char str[3];
          n = 0;
          str[n++] = line[1];
          if (line[2]!=' ') str[n++] = line[2];
          str[n++] = '\0';
          n = SymbolToAtomicNumber(str);
          for (iatom = 0 ; iatom < natom ; iatom++)
          {
            if (mol->atom[iatom]->atomic_number==n)
            {
              if (line[12]==':') 
                set[iatom]->SetName(&line[14]);
              else
                set[iatom]->SetName(&line[15]);
            }
          }
          f.getline(line,MAXSTR);
          stringLine=line;
          sindex=stringLine.find("BASIS");
        }
        break;
      }
    }
    f.close();
  } 
  LoadFromLibrary("/usr/local/BASIS",fmt,ELECTRON);
}

void Basis::LoadFromLibrary(const char* libfile,const char* fmt,int type)
{
  int natom = mol->Natom();
  nco = 0;
  cout << " Loading basis for "<<natom<<" atoms "<<endl;
  for (int iatom=0;iatom<natom;iatom++)
  {
    cout << " Atom "<<iatom+1<<endl;
    if (type==HYDROGEN)
    {
      if (mol->atom[iatom]->atomic_number==1)
        set[iatom]->ReadFile(libfile,fmt);
    }
    else set[iatom]->ReadFile(libfile,fmt);
    if (iatom==0) set[iatom]->ll = 0;
    else set[iatom]->ll = set[iatom-1]->ul + 1;
    set[iatom]->ul = set[iatom]->ll + set[iatom]->nco - 1;
    nco += set[iatom]->nco;
  }
};


/*
int Basis::NumberOfSphericalOrbitals()
{
  int is;

  int nso = 0;
  for (is=0;is<(signed)set.size();is++)
  {
    nso = nso + set[is]->NumberOfSphericalOrbitals();
  }
  return nso;
}

size_t Basis::NumberOfShells()
{
  size_t nshl = 0;
  for (size_t is=0;is<set.size();is++)
  {
    nshl += set[is]->NumberOfShells();
  }
  return nshl;
}

// Notice that auxiliary function sets are different to 
// original GEN-A2* for the case of hydrogen atoms
void Basis::GenerateAuxiliary(Basis* basis, Integrator* gi)
{
  int iatom,natom;

  nco = 0;
  natom = mol->Natom();
  set.clear();
  for (iatom = 0 ; iatom < natom ; iatom++)
  {
    set.push_back( new Set(ELEMENT_NAME[mol->atom[iatom]->atomic_number],"GEN"));
    set[iatom]->GenerateAuxiliary(basis->set[iatom],gi);

    if (iatom==0) set[iatom]->ll = 0;
    else set[iatom]->ll = set[iatom-1]->ul + 1;
    set[iatom]->ul = set[iatom]->ll + set[iatom]->nco - 1;
    nco += set[iatom]->nco;
  }
}
*/

// Load primary basis
void Basis::ReadLibrary(const char* libfile,const char* basname, int type)
{
  int iatom,natom;
  char line[MAXSTR];

  natom = mol->Natom();

  if (natom<1)
  {
    cout << "WARNING: Loading basis for ZERO atoms"<<endl;
    return;
  }

  set.clear();
  for (iatom = 0 ; iatom < natom ; iatom++)
    set.push_back( new Set(ELEMENT_NAME[mol->atom[iatom]->atomic_number],basname));
  LoadFromLibrary(libfile,"D2K",type);
}

