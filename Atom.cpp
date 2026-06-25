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

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include <Atom.h>
#include <Element.h>
/*
#include <AM.h>
#include <Parameter.h>
#include <Units.h>
*/

Atom::Atom(char* il,double x,double y, double z)
    : Vector(x,y,z)
{
  strcpy(label,il);
  atomic_number = SymbolToAtomicNumber( il );
  Setup();
};

Atom::Atom(int an,double x,double y, double z)
    : Vector(x,y,z)
{
  strcpy(label,(char*)ELEMENT_SYMBOL[an]);
  atomic_number = an;
  Setup();
};

void Atom::Setup()
{
/*
  nuclear_magnetic_moment[0] = 0.0;
  nuclear_magnetic_moment[1] = 0.0;
  nuclear_magnetic_moment[2] = 0.0;
  */

  // notice that charge is reset
  SetCharge( double(atomic_number) );

  /*
  if (atomic_number==1)
    nuclear_spin = new AM(2,1);
  else
    nuclear_spin = new AM(1,1);
    */
};

/*
void Atom::SetAtomicNumber(int an)
{
  atomic_number = an;

  // notice that charge is reset
  SetCharge( double(atomic_number) );
}
*/

void Atom::SetSymbol(char* is)
{
  strcpy(label,is);
  atomic_number = SymbolToAtomicNumber( label );

  // notice that charge is reset
  SetCharge( double(atomic_number) );
}

void Atom::SetCharge(double c)
{
  charge = c;
}

double Atom::GetCovalentRadius(void)
{
  double r;

  r = ELEMENT_COV_R[atomic_number];

  return r;
};

double Atom::GetVDWRadius(void)
{
  double r;

  r = ELEMENT_VDW_R[atomic_number];

  return r;
};

static double color[3];
double* Atom::GetColor(void)
{
  color[0] = ELEMENT_COLOR[atomic_number][0];
  color[1] = ELEMENT_COLOR[atomic_number][1];
  color[2] = ELEMENT_COLOR[atomic_number][2];

  return color;
};

int Atom::SymbolToAtomicNumber( char *sym )
{
  int an = 0;
  char tst[3];

  if ( strlen(sym) == 0 )
  {
    return 0;
  }
  else
  {
    tst[0] = toupper(sym[0]);
    if ( strlen(sym) == 1 )
    {
      tst[1] = '\0';
    }
    else
    {
      if ( sym[1] >= '0' && sym[1] <= '9' )
      {
        tst[1] = '\0';
      }
      else
      {
        tst[1] = tolower(sym[1]);
        tst[2] = '\0';
      };
    };
  };

  for ( int i = 0 ; i <= MAXATNUM ; i++ )
  {
    if ( strcmp( tst , ELEMENT_SYMBOL[i] ) == 0 )
    {
      an = i;
      break;
    };
  };
  return an;

}

double Atom::Distance(Atom* other)
{
  double d;
  Vector v;

  v = Vector(this->x,this->y,this->z);
  v -= Vector(other->x,other->y,other->z);
  d = v.Norm(); 
  return d;
};

