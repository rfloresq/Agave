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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <QtGui>

#include <Atom.h>
#include <Structure.h>
#include <Units.h>

using namespace std;
using std::ifstream;

Structure::Structure(void)
{
  atom.clear();
};

Structure::Structure(Structure* other)
{
  atom.clear();
  for (int iatom = 0 ; iatom < (signed)other->atom.size() ; iatom++ )
    atom.push_back( new Atom( other->atom[iatom]->atomic_number,
			      other->atom[iatom]->x,
			      other->atom[iatom]->y,
			      other->atom[iatom]->z));
};

int Structure::Natom(void)
{
  return atom.size();
}


void Structure::Center(void)
{
  Vector ref(0.0,0.0,0.0);
  for (int iatom = 0 ; iatom < (signed)atom.size() ; iatom++ )
    ref += Vector(atom[iatom]->x,atom[iatom]->y,atom[iatom]->z);
  ref /= (double)(atom.size());

  for (int iatom = 0 ; iatom < (signed) atom.size() ; iatom++ )
  {
    atom[iatom]->x -= ref.x;
    atom[iatom]->y -= ref.y;
    atom[iatom]->z -= ref.z;
  }
};

/*
void Structure::DeleteAtom(int id)
{
  atom.erase(atom.begin()+size_t(id));
}
*/

Atom* Structure::GetAtom(int id)
{
  if (id<(signed)atom.size())
    return atom[id];
  else
  {
    exit(EXIT_FAILURE);
    return atom[0];
  }
}

void Structure::ReadXYZ(const char *filename)
{
  ifstream f(filename);
  if ( f.fail() )
  {
    cout << "Error opening file "<<filename<<endl;
    return;
  };

  int iatom,n;
  char str[MAXSTR];
  double x,y,z;
  Atom* c;

  f >> n;
  if (n<1) cout << "WARNING in ReadXYZ: Number of atoms non positive" <<endl;
  f.getline(str,MAXSTR);
  f.getline(str,MAXSTR);
  bool inbohr;
  if (strcasestr(str,"BOHR")!=NULL) inbohr = true;
  else inbohr = false;
  atom.clear();
  for ( iatom = 0 ; iatom < n ; iatom++ )
  {
    f >> str;
    f >> x;
    f >> y;
    f >> z;
    if (!inbohr)
    {
      x = AngstromToBohr(x);
      y = AngstromToBohr(y);
      z = AngstromToBohr(z);
    }
    c = new Atom(str,x,y,z);
    atom.push_back( c );
  }
  Center();

  f.close();
}


// Read deMon2k output file
void Structure::ReadD2K(const char*filename)
{
  ifstream f(filename);
  if ( f.fail() )
  {
    cout << __FILE__ << __LINE__ << " Error opening file"<<endl;
    void ErrorMessage(const char*, int, const char*, int);
    return;
  };
  
  Atom* c;
  int natom;
  char line[MAXSTR];
  char str[MAXSTR];
  double x,y,z;

  while (f.getline(line,MAXSTR))
  {
    int sindex;
    string stringLine=line;

    sindex=stringLine.find(" NUMBER OF ATOMS:");
    if((sindex>=0)&&(sindex<(signed)stringLine.size()))
    {
      natom = atoi(&line[17]);
      f >> str;
      f >> str;
      f >> str;
      f >> str;
      continue;
    }
    sindex=stringLine.find("INPUT ORIENTATION IN");
    //if(!((sindex>=0)&&(sindex<(signed)stringLine.size())))
    //  sindex=stringLine.find("STANDARD ORIENTATION IN");
    //if(!((sindex>=0)&&(sindex<(signed)stringLine.size())))
    //  sindex=stringLine.find("COORDINATES OF OPTIMIZATION STEP");
    if((sindex>=0)&&(sindex<(signed)stringLine.size()))
    {
      bool AA = false;
      if ( strstr( line, "ANGSTROM" ) != NULL ) AA = true;
      f.getline(line,MAXSTR);
      f.getline(line,MAXSTR);
      f.getline(line,MAXSTR); 
      atom.clear();
      for( int iatom = 0; iatom < natom ; iatom++ )
      {
        f >> str;
        f >> str;
        f >> x;
        f >> y;
        f >> z;
        if (AA) 
        {
          x = AngstromToBohr(x);
          y = AngstromToBohr(y);
          z = AngstromToBohr(z);
        }
        c = new Atom(str,x,y,z);
        atom.push_back( c );
        f.getline(line,MAXSTR);
      }
      Center();
    }
  }
  f.close();
}

/*
void Structure::ReadPIE(const char *filename)
{
  ifstream f(filename);
  if ( f.fail() )
  {
    cout << "Error opening file "<<filename<<endl;
    return;
  };

  char line[MAXSTR];
  char str[MAXSTR];
  int n,natom,sindex;
  bool inbohr;
  double x,y,z;
  Atom *c;

  natom = 0;
  inbohr = true;
  while (f.getline(line,MAXSTR))
  {
    string stringLine=line;

    sindex=stringLine.find("#define ANGSTROM");
    if((sindex>=0)&&(sindex<(signed)stringLine.size()))
    {
      inbohr = false;
    }

    sindex=stringLine.find("// Structure coordinates");
    if((sindex>=0)&&(sindex<(signed)stringLine.size()))
    {
      atom.clear();
      f.getline(line,MAXSTR);
      f.getline(line,MAXSTR);
      f.getline(line,MAXSTR);
      f.getline(line,MAXSTR);
      while ( line[0] != '}' )
      {
	stringstream ss(line);
	ss >> x;
	ss >> y;
	ss >> z;
        if (!inbohr)
        {
          x = AngstromToBohr(x);
          y = AngstromToBohr(y);
          z = AngstromToBohr(z);
        }
        c = new Atom("H",x,y,z);
        atom.push_back( c );
        f.getline(line,MAXSTR);
      }
    }
    sindex=stringLine.find("ATOM");
    if((sindex>=0)&&(sindex<(signed)stringLine.size()))
    {
      stringstream ss(line);
      ss >> str;
      ss >> str;
      ss >> str;
      n = strlen(str);
      str[n-1] = '\0';
      atom[natom]->SetSymbol( str );
      natom++;
    }
  }

  f.close();
}

*/


void Structure::ReadFile(const char*filename,const char*fmt)
{
  if ( strcasecmp( fmt , "AGV" ) == 0 )
  {
    ReadXYZ(filename);
  }
  else if ( strcasecmp( fmt , "D2K" ) == 0 )
  {
    ReadD2K(filename);
  }
  else if ( strcasecmp( fmt , "PIE" ) == 0 )
  {
    // ReadPIE(filename);
  }
  else if ( strcasecmp( fmt , "XYZ" ) == 0 )
  {
    ReadXYZ(filename);
  }
  else
  {
    cout << "Unknown format while reading Structure"<<endl;
  }
};

void Structure::WriteFile(const char* filename)
{
  ofstream f(filename,ios::app);
  f << "   " << Natom() << endl;
  f << "Geometry is given in angstroms "<<endl;
  for (int iatom=0;iatom<Natom();iatom++)
  {
    Atom* a = atom[iatom];
    f << a->label << "   " 
      << setw(15) << fixed << BohrToAngstrom(a->x) << "   " 
      << setw(15) << fixed << BohrToAngstrom(a->y) << "   " 
      << setw(15) << fixed << BohrToAngstrom(a->z) << "   " 
      << endl;
  }
  f << endl;
  f << endl;
  f.close();
};

// Comput nuclear repulsion energy in atomic units
double Structure::NuclearRepulsionEnergy()
{
  int iatom,jatom;
  double d,nr,x,y,z;

  nr = 0.0;
  for (iatom = 0 ; iatom < (signed)atom.size() ; iatom++ )
  {
    for (jatom = iatom+1 ; jatom < (signed)atom.size() ; jatom++ )
    {
      if (atom[iatom]->charge>0&&atom[jatom]->charge>0)
      {
        x = atom[iatom]->x - atom[jatom]->x;
        y = atom[iatom]->y - atom[jatom]->y;
        z = atom[iatom]->z - atom[jatom]->z;
        d = sqrt(x*x+y*y+z*z);
        nr += atom[iatom]->charge*atom[jatom]->charge/d;
      }
    }
  }
  return nr;
}

/*
void Structure::NuclearDipole(double *dipole)
{
  size_t iatom;
  int i;
  double charge;

  for (i=0;i<3;i++)
    dipole[i] = 0.0;

  for (iatom=0;iatom<atom.size();iatom++)
  {
    charge = atom[iatom]->charge;
    dipole[0] += charge*atom[iatom]->x;
    dipole[1] += charge*atom[iatom]->y;
    dipole[2] += charge*atom[iatom]->z;
  }
}
*/

int Structure::NumberOfElectronsForNeutral()
{
  int nelec;

  nelec = 0.0;
  for (size_t iatom = 0 ; iatom < atom.size() ; iatom++ )
    nelec += int(atom[iatom]->charge);

  return nelec;
}

void Structure::NullifyHydrogens()
{
  for (size_t iatom=0;iatom<atom.size();iatom++)
    if (atom[iatom]->atomic_number==1) 
      atom[iatom]->charge = 0.0;
};

/*
double Structure::Mass()
{
  int  an;   
  double mass;

  double total_mass = 0.0;
  for (int iatom = 0 ; iatom < (signed)atom.size() ; iatom++ )
  {
    an = atom[iatom]->atomic_number;
    mass = ELEMENT_MASS[an];
    total_mass += mass;
  }
  return total_mass;
};

void Structure::CenterMass(void)
{
  int an,iatom;
  double mass,total_mass;

  total_mass = 0.0;
  Vector ref(0.0,0.0,0.0);
  for (iatom=0;iatom<(signed)atom.size();iatom++ )
  {
    Vector a(atom[iatom]->x,atom[iatom]->y,atom[iatom]->z);
    an = atom[iatom]->atomic_number;
    mass = ELEMENT_MASS[an];
    total_mass += mass;
    a *= mass;
    ref += a;
  }
  ref /= total_mass;

  for (iatom=0;iatom<(signed)atom.size();iatom++ )
  {
    atom[iatom]->x -= ref.x;
    atom[iatom]->y -= ref.y;
    atom[iatom]->z -= ref.z;
  }
};

void Structure::ReadZMatrix(const char *filename)
{
  ifstream f(filename);
  if ( f.fail() ) return;

  char line[MAXSTR];
  char str[MAXSTR];
  int iatom,n;
  Atom* c;

  atom.clear();
  f >> n;
  f.getline(line,MAXSTR);
  f.getline(line,MAXSTR);
  for ( iatom = 0 ; iatom < n ; iatom++ )
  {
    f >> str;
    atom.push_back( new Atom(str,0.0,0.0,0.0) );
    c = atom[iatom];
    if (iatom>0) 
    {
      f >> c->ref_bond;
      f >> c->bond;
      c->bond = AngstromToBohr(c->bond);
    }
    if (iatom>1) 
    {
      f >> c->ref_angle;
      f >> c->angle;
    }
    if (iatom>2) 
    {
      f >> c->ref_dihedral;
      f >> c->dihedral;
    }
  }
  f.close();

  // Obtain cartesian representation
  Z2C(0); 
}
*/

void Structure::Z2C( int llatom ) 
{
  int iatom,bondref,angleref,dihedralref;

  for ( iatom = llatom ; iatom < (signed)atom.size() ; iatom++ )
  {
    bondref = atom[iatom]->ref_bond;

    if (iatom>0)
    {
      atom[iatom]->x = atom[bondref]->x;
      atom[iatom]->y = atom[bondref]->y;
      atom[iatom]->z = atom[bondref]->z;
    }
    else
    {
      atom[iatom]->x = 0.0;
      atom[iatom]->y = 0.0;
      atom[iatom]->z = 0.0;
    }

    if ( iatom > 0 )
    {
      Vector shift(0.0,0.0,1.0);
      if ( iatom > 1 )
      {
        angleref = atom[iatom]->ref_angle;
        shift.x = atom[angleref]->x - atom[bondref]->x;
        shift.y = atom[angleref]->y - atom[bondref]->y;
        shift.z = atom[angleref]->z - atom[bondref]->z;
        shift ^= (1.0);    
      };
      shift *= atom[iatom]->bond;
      atom[iatom]->x += shift.x;
      atom[iatom]->y += shift.y;
      atom[iatom]->z += shift.z;

      if ( iatom > 1 )
      {
        Vector axis;
        shift.x = atom[angleref]->x - atom[bondref]->x;
        shift.y = atom[angleref]->y - atom[bondref]->y;
        shift.z = atom[angleref]->z - atom[bondref]->z;
        if ( iatom == 2 )
        {
          axis = Vector(0.0,1.0,0.0) < shift;
        }
        else
        {
          dihedralref = atom[iatom]->ref_dihedral;
          axis.x = atom[dihedralref]->x - atom[bondref]->x;
          axis.y = atom[dihedralref]->y - atom[bondref]->y;
          axis.z = atom[dihedralref]->z - atom[bondref]->z;
          axis <= shift;
        };
        shift.x = atom[iatom]->x - atom[bondref]->x;
        shift.y = atom[iatom]->y - atom[bondref]->y;
        shift.z = atom[iatom]->z - atom[bondref]->z;
        if (shift.Colineal(axis)) axis = shift.Orthogonal();
        shift.Rotate(axis,atom[iatom]->angle);
        atom[iatom]->x = atom[bondref]->x + shift.x;
        atom[iatom]->y = atom[bondref]->y + shift.y;
        atom[iatom]->z = atom[bondref]->z + shift.z;

        if ( iatom > 2 )
        {
          axis.x = atom[angleref]->x - atom[bondref]->x;
          axis.y = atom[angleref]->y - atom[bondref]->y;
          axis.z = atom[angleref]->z - atom[bondref]->z;
          shift.x = atom[iatom]->x - atom[bondref]->x;
          shift.y = atom[iatom]->y - atom[bondref]->y;
          shift.z = atom[iatom]->z - atom[bondref]->z;
          shift.Rotate(axis,atom[iatom]->dihedral);
          atom[iatom]->x = atom[bondref]->x + shift.x;
          atom[iatom]->y = atom[bondref]->y + shift.y;
          atom[iatom]->z = atom[bondref]->z + shift.z;
        };
      };
    };
  };
  Center();
};


void Structure::AddAtom(int an, int r1, int r2,int r3,
                   double bond,double angle,double dihedral)
{
  Atom* a;

  a = new Atom(an,0.0,0.0,0.0);

  a->ref_bond = r1;
  a->ref_angle = r2;
  a->ref_dihedral = r3;

  a->bond = bond;
  a->angle = angle;
  a->dihedral = dihedral;

  atom.push_back( a );

  Z2C( atom.size()- 1 ); 
}

void Structure::DeleteLastAtom(void)
{
  int n = atom.size();
  if (n>0) atom.resize( n-1 );
  Center();
};


/*
// Point Group operations applier driver
int Structure::PGBDrv(char* pg,int* ref) 
{
  int n;

  // Aliasing
  if ( strcmp(pg,"Cs") == 0 ) { strcpy(pg,"C1h"); }
  else if ( strcmp(pg,"Ci") == 0 ) { strcpy(pg,"S2"); }
  else if ( strcmp(pg,"Vh") == 0 ) { strcpy(pg,"D2h"); } 
  else if ( strcmp(pg,"Vd") == 0 ) { strcpy(pg,"D2d"); }
  else if ( strcmp(pg,"S6v") == 0 ) { strcpy(pg,"D3d"); }
  else if ( strcmp(pg,"S8v") == 0 ) { strcpy(pg,"D4d"); };

  n = 0;
  if ((pg[0] != 'T')&&(pg[0] != 'O')&&(pg[0] != 'I')&&
      (pg[1] != '*')&&(pg[1] != 's'))
  {
    char s[10];
    sprintf(s,"%c",pg[1]);
    n = atoi(s);
  };

  Vector v(atom[ref[0]]->x,
	   atom[ref[0]]->y,
	   atom[ref[0]]->z);
  for (int iatom=0;iatom<(signed)atom.size();iatom++)
  {
    atom[iatom]->x -= v.x;
    atom[iatom]->y -= v.y;
    atom[iatom]->z -= v.z;
  }

  // Applying all operations to create missing atoms
  if ( pg[0] == 'C' ) { AppCGroup(pg,n,&ref[1]); }
  else if ( pg[0] == 'D' ) { AppDGroup(pg,n,&ref[1]); }
  else if ( pg[0] == 'S' ) { AppSGroup(n,&ref[1]); }
  else if ( pg[0] == 'T' ) { AppTGroup(&ref[1]); }
  else if ( pg[0] == 'O' ) { AppOGroup(&ref[1]); }
  else if ( pg[0] == 'I' ) { AppIGroup(&ref[1]); };

  Center();

  return 0;
};


// Purpose: Space Group operations applier driver
int Structure::SGBDrv(char* sg,int* ref) 
{
  int nr;

  nr = atom.size();

  if ( sg[1] == '1' )
  {
    cout << "PATIAL GROUPS NOT READY"<<endl;
    Vector axis(atom[ref[1]]->x,atom[ref[1]]->y,atom[ref[1]]->z);
    axis.x -= atom[ref[0]]->x;
    axis.y -= atom[ref[0]]->y;
    axis.z -= atom[ref[0]]->z;
    double sd = QInputDialog::getDouble(0,
                "Translation Unit Length","R = ",1.0,0.0,10000.0,3);
    double nt = QInputDialog::getInt(0,
                "Number on Monomers","NM = ",2,2,10000);
    nt = nt - 1;
    double s = 0.0;
    for ( int it = 0 ; it < nt ; it++ )
    {
      s = s + sd;
      AppTrans(axis,s,nr);
    }
  }
  else 
  {
    cout << "sgbdrv: Unsupported spatial group"<<endl;
    exit(EXIT_FAILURE);
  };

  Center();
  return 0;
};

// Apply I group operations 
int Structure::AppIGroup(int *ref)
{
  int i,iaxis,nr;

  nr = atom.size();

  Vector v,axisa,axisb,axis[12];
  Vector axis2[20],axis3[15],plane[15];

  axisa = Vector(atom[ref[0]]->x,atom[ref[0]]->y, atom[ref[0]]->z);
  axisa /= axisa.Norm();
  axisb = Vector(atom[ref[1]]->x,atom[ref[1]]->y, atom[ref[1]]->z);
  v = Vector( axisa );
  v *= axisb.Dot( axisa )/axisa.Dot( axisa );
  axisb -= v;
  axisb /= axisb.Norm();
  v = axisa > axisb;

  axis[0] = axisa;
  axis[1] = axisa;
  axis[1].Rotate( v , atan(2.0)*180.0/M_PI );
  for ( iaxis = 2 ; iaxis < 6 ; iaxis++ )
  {
    axis[iaxis] = axis[iaxis-1];
    axis[iaxis].Rotate( axisa , 72.0 );
  };
  axis[6] = axis[4];
  axis[6] *= -1.0; 
  axis[7] =  axis[5];
  axis[7] *= -1.0; 
  axis[8] = axis[1];
  axis[8] *= -1.0; 
  axis[9] = axis[2];
  axis[9] *= -1.0; 
  axis[10] = axis[3];
  axis[10] *= -1.0; 
  axis[11] = axis[0];
  axis[11] *= -1.0; 

  for ( iaxis = 0 ; iaxis < 12 ; iaxis++ )
  {
    for ( i = 0 ; i < 4 ; i++ )
    {
      AppSAxis( 10 , i + 1 , axis[iaxis] , nr );
    };
  };

  axis2[0] = axis[0] + axis[1] + axis[2];
  axis2[0] /= axis2[0].Norm();
  axis2[5] = axis[1] + axis[2] + axis[6];
  axis2[5] /= axis2[5].Norm();
  for ( iaxis = 1 ; iaxis < 5 ; iaxis++ )
  {
    axis2[iaxis] = axis2[iaxis-1];
    axis2[iaxis].Rotate( axisa , 72.0 );
    axis2[iaxis+5] = axis2[iaxis+4];
    axis2[iaxis+5].Rotate( axisa , 72.0 );
  };
  for ( iaxis = 10 ; iaxis < 20 ; iaxis++ )
  {
    axis2[iaxis] = axis2[iaxis-10];
    axis2[iaxis] *= -1.0; 
  };

  for ( iaxis = 0 ; iaxis < 20 ; iaxis++ )
  {
    for ( i = 0 ; i < 2 ; i++ ) 
    {
      AppSAxis( 6 , i + 1 , axis2[iaxis] , nr );
    };
  };

  for ( iaxis = 0 ; iaxis < 4 ; iaxis++ )
  {
    axis3[iaxis] = axis[0] + axis[iaxis+1];
    axis3[iaxis+5] = axis[iaxis+1] + axis[iaxis+2];
  };
  axis3[4] = axis[0] + axis[5];
  axis3[9] = axis[5] + axis[1];
  axis3[10] = axis[1] + axis[6];
  for ( iaxis = 11 ; iaxis < 15 ; iaxis++ )
  {
    axis3[iaxis] = axis3[iaxis-1];
    axis3[iaxis].Rotate( axisa , 72.0 ); 
  };

  for ( iaxis = 0 ; iaxis < 15 ; iaxis++ )
  {
    AppCAxis( 2 , 1 , axis3[iaxis] , nr );
  };

  for (iaxis = 0 ; iaxis < 5 ; iaxis++ )
  {
    plane[iaxis] = axisa > axis[iaxis+1];
    plane[iaxis+10] = axis[iaxis+1] > axis[iaxis+6];
  };
  for ( iaxis = 0 ; iaxis < 4 ; iaxis++ )
  {
    plane[iaxis+5] = axis[iaxis+1] > axis[iaxis+2];
  };
  plane[9] = axis[5] > axis[1];

  for ( iaxis = 0 ; iaxis < 15 ; iaxis++ )
  {
    AppSigma( plane[iaxis] , nr );
  };

  AppSAxis( 2 , 1 , axisa , nr );

  return 0;
};


// Apply O group operations 
int Structure::AppOGroup(int *ref)
{
  int i,j,iaxis;    
  int nr;            
  Vector axisa;    
  Vector c2axis;  
  Vector axis[5];

  nr = atom.size();

  axis[0] = Vector(atom[ref[0]]->x,atom[ref[0]]->y,atom[ref[0]]->z);
  axis[0] /= axis[0].Norm();
  axisa = axis[0];
  axis[1] = Vector(atom[ref[1]]->x,atom[ref[1]]->y,atom[ref[1]]->z);
  axis[2] = Vector( axis[0] );
  axis[2] *= axis[1].Dot( axis[0] )/axis[0].Dot( axis[0] );
  axis[1] -= axis[2];
  axis[1] /= axis[1].Norm();
  axis[2] = axis[0] > axis[1];
  axis[3] = axis[0];
  axis[4] = axis[1];

  for ( iaxis = 0 ; iaxis  < 3 ; iaxis++ )
  {
    c2axis = axis[iaxis+1] + axis[iaxis+2];
    c2axis /= c2axis.Norm(); 
    for ( i = 0 ; i < 3 ; i++ )
    {
      AppSAxis( 4 , i+1, axis[iaxis] , nr );
      AppCAxis( 4 , i+1, axis[iaxis] , nr );
      AppCAxis( 2 , 1, c2axis , nr );
      c2axis.Rotate( axis[iaxis] , 90.0 );
      AppCAxis( 2 , 1, c2axis , nr );
      for ( j = 0 ;  j < 3 ; j++ )
      {
        AppSigma( axis[iaxis+j] , nr );
      };
    };
  };

  AppSAxis( 2 , 1 , axisa , nr );

  axis[0] += axis[1];
  axis[0] += axis[2];
  axis[0] /= axis[0].Norm();
  for ( iaxis = 1 ; iaxis < 4 ; iaxis++ )
  {
    axis[iaxis] = axis[iaxis-1];
    axis[iaxis].Rotate( axisa , 90.0 );
  };
 
  for ( iaxis = 0 ; iaxis < 4 ; iaxis++ )
  {
    for ( i = 0 ; i < 5 ; i++ )
    {
      AppSAxis( 6 , i+1, axis[iaxis] , nr );
    };
  };

  return 0;
};

// Apply T group operations 
int Structure::AppTGroup(int *ref)
{
  int       i,iaxis,nr;
  double    tangle; 
  Vector axisb;  
  Vector v;
  Vector axis[4]; 

  nr = atom.size();

  axis[0] = Vector(atom[ref[0]]->x,atom[ref[0]]->y,atom[ref[0]]->z);
  axis[0] /= axis[0].Norm();
  axisb = Vector(atom[ref[1]]->x,atom[ref[1]]->y,atom[ref[1]]->z);
  axis[2] = Vector( axis[0] );
  axis[2] *= axisb.Dot( axis[0] )/axis[0].Dot( axis[0] );
  axisb -= axis[2];
  axisb /= axisb.Norm();

  tangle = acos(-1.0/3.0);
  axis[2] = axis[0] > axisb;
  axis[1] = axis[0];
  axis[1].Rotate( axis[2] , 180.0*tangle/M_PI );
  axis[2] = axis[1];
  axis[2].Rotate( axis[0] , 120.0 );
  axis[3] = axis[2];
  axis[3].Rotate( axis[0] , 120.0 );

  for ( iaxis = 0 ; iaxis  < 4 ; iaxis++ )
  {
    AppCAxis(3,1, axis[iaxis], nr );
    AppCAxis(3,2, axis[iaxis], nr );
    if ( iaxis == 0 )
    {
      v = axis[0] > axisb;
    }
    else 
    {
      v = axis[iaxis-1] > axis[iaxis];
    };
    for ( i = 0 ; i < 3 ; i++ )
    {
      AppSigma( v , nr );
      v.Rotate( axis[iaxis] , 120.0 );
    };
  };

  for ( iaxis = 0 ; iaxis  < 3 ; iaxis++ )
  {
    if ( iaxis == 2 )
    {
      v = axis[iaxis] + axis[1];
    }
    else
    {
      v = axis[iaxis] + axis[iaxis+1];
    };
    for ( i = 0 ; i < 3 ; i++ )
    {
      AppSAxis(4, i + 1, v , nr );
    };
  };

  return 0;
};

// Apply S group operations 
int Structure::AppSGroup(int n,int *ref)
{
  int  i;    
  int nr;    
  char s[10];
  
  if ( n%2 == 1 )
  {
    sprintf(s,"C%dh",n);
    return AppCGroup(s,n,ref);
  }
  else
  {
    nr = atom.size();

    Vector axisa(atom[ref[0]]->x,atom[ref[0]]->y,atom[ref[0]]->z);
    axisa /= axisa.Norm();

    if ( n >= 4 )
    {
      for ( i = 1 ; i <= n/2-1 ; i++ )
      {
        AppCAxis( n , i , axisa , nr );
      };
    };

    for ( i = 1 ; i < n ; i++ )
    {
      AppSAxis( n , i , axisa , nr );
    };

  };

  return 0;
};

// Apply D group operations 
int Structure::AppDGroup(char* pg,int n,int *ref)
{
  int i,top,step,nr; 

  nr = atom.size();

  Vector axisa(atom[ref[0]]->x,atom[ref[0]]->y,atom[ref[0]]->z);
  axisa /= axisa.Norm();

  if ( pg == (char*)"D*h" )
  {
    AppSAxis(2,1,axisa,nr);
    return 0;
  };
  
  for ( i = 1 ; i < n ; i++ )
  {
    AppCAxis( n , i , axisa , nr );
  };

  Vector axisb(atom[ref[1]]->x,atom[ref[1]]->y,atom[ref[1]]->z);
  Vector v( axisa );
  v *= axisb.Dot( axisa )/axisa.Dot( axisa );
  axisb -= v;
  axisb /= axisb.Norm();
  v = axisb;
  for ( i = 1 ; i <= n ; i++ )
  {
    AppCAxis( 2 , 1 , v , nr );
    v.Rotate( axisa , 180.0/((double)n) );
  };

  if ( pg[2] == 'h' )
  {
    if ( n%2 == 0 )
    {
      top = n - 1;
      step = 1;
    }
    else
    {
      top = 2*n-1;
      step = 2;
    };
    for ( i = 1 ; i <= top ; i += step )
    {
      AppSAxis( n , i , axisa , nr );
    };
  }
  else if ( pg[2] == 'd' )
  {
    for ( i = 1 ; i <= 2*n-1 ; i += 2 )
    {
      AppSAxis( n , i , axisa , nr );
    };
  };

  v = axisa > axisb;
  if ( pg[2] == 'd' )
  {
    v.Rotate( axisa , 180.0/((double)n) );
  };
  for ( i = 1 ; i <= n ; i++ )
  {
    AppSigma( v , nr );
    v.Rotate( axisa , 360.0/((double)n) );
  };

  if ( pg[2] == 'h' )
    AppSigma( axisa , nr );

  return 0;
};

// Apply C group operations 
int Structure::AppCGroup(char* pg,int n,int *ref)
{
  int i,nr;  

  if ( n < 2 ) return 0;

  Vector axisa(atom[ref[0]]->x,atom[ref[0]]->y,atom[ref[0]]->z);

  nr = atom.size();
  for ( i = 1 ; i < n ; i++ )
  {
    AppCAxis(n,i,axisa,nr);
  };
  
  if ( strlen(pg) == 2 ) return 0;
  
  if ( pg[2] == 'h' )
  {
    for ( i = 1 ; i < n ; i++ )
    {
      AppSAxis( n , i , axisa , nr );
    };
    AppSigma( axisa , nr );
  }
  else 
  {
    Vector axisb(atom[ref[1]]->x,atom[ref[1]]->y,atom[ref[1]]->z);
    Vector normal( axisa );
    normal *= axisb.Dot( axisa )/axisa.Dot( axisa );
    axisb -= normal;
    normal = axisa > axisb;
    for ( i = 1 ; i <= n ; i++ )
    {
      AppSigma(normal,nr);
      normal.Rotate(axisa,360.0/((double)n));
    };
  };
  return 0;
};


// Apply sigma plane reflexion 
void Structure::AppSigma(Vector normal,int nr)
{
  int iatom;
  double sp;
  Vector va,vb;

  normal ^= 1.0;

  for ( iatom = 0 ; iatom < nr ; iatom++ )
  {
    if ( atom[iatom]->atomic_number == 0 ) continue;

    va = Vector(atom[iatom]->x,atom[iatom]->y,atom[iatom]->z);
    sp = va.Dot(normal);
    vb = normal;
    vb *= 2.0*sp;
    va  -= vb;
    ClonAtom(iatom,va);
  };
};

// Apply Cn rotation operation 
void Structure::AppCAxis(int n,int i,Vector axis,int nr) 
{
  int iatom; 
  double phi; 
  Vector v;

  phi = 360.0*((double)i)/((double)n);
  axis ^= 1.0;

  for ( iatom = 0 ; iatom < nr ; iatom++ )
  {
    if ( atom[iatom]->atomic_number == 0 ) continue;

    v = Vector(atom[iatom]->x,atom[iatom]->y,atom[iatom]->z);
    v.Rotate(axis,phi); 

    ClonAtom(iatom,v);
  };
};

// Apply Sn rotation operation 
void Structure::AppSAxis(int n,int i,Vector axis,int nr)
{
  int iatom;
  double sp,phi;
  Vector v,va,vb;

  if ( i%2 == 0 )
  {
    AppCAxis(n,i,axis,nr);
  }
  else
  {
    phi = 360.0*((double)i)/((double)n);
    axis ^= 1.0;

    for ( iatom = 0 ; iatom < nr ; iatom++ )
    {
      if ( atom[iatom]->atomic_number == 0 ) continue;

      v = Vector(atom[iatom]->x,atom[iatom]->y,atom[iatom]->z);
      v.Rotate(axis,phi); 
      sp = Vector( v ).Dot( axis );
      va = Vector( v );
      vb = axis;
      vb *= 2.0*sp;
      va  -= vb;
      ClonAtom(iatom,va);
    };
  };
}; 

// Apply translation
void Structure::AppTrans(Vector axis,double shift, int nr)
{
  int iatom;
  Vector normal,v;

  normal = axis;
  normal ^= 1.0;
  normal *= shift;

  for ( iatom = 0 ; iatom < nr ; iatom++ )
  {
    if ( atom[iatom]->atomic_number == 0 ) continue;

    v = Vector(atom[iatom]->x,atom[iatom]->y,atom[iatom]->z);
    v += normal;
    ClonAtom(iatom,v);
  };
};

// Make a copy of a non-redundant atom if required
void Structure::ClonAtom(int from,Vector newpos)
{
  Vector v;  
  
  for ( int iatom = 0 ; iatom < (signed)atom.size() ; iatom++ )
  {
    if ( atom[iatom]->atomic_number == 0 ) continue;
    v = newpos;
    v.x -= atom[iatom]->x;
    v.y -= atom[iatom]->y;
    v.z -= atom[iatom]->z;
    if ( v.Norm() < 0.5 )
    {
      return;
    };
  };

  Atom *a = new Atom(atom[from]->atomic_number,newpos.x,newpos.y,newpos.z);
  atom.push_back( a );
};

void Structure::WriteXYZ(const char *filename)
{
  ofstream f(filename);
  f << atom.size() <<endl;
  f << endl;
  for (int i=0;i<(signed)atom.size();i++)
  {
    f << atom[i]->label << "  " 
      << BohrToAngstrom(atom[i]->x) << "   "
      << BohrToAngstrom(atom[i]->y) << "   "
      << BohrToAngstrom(atom[i]->z) << endl;
  }
  f.close();
};
*/

void Structure::C2Z( int llatom ) 
{
  int iatom;
  double val;
  Vector va,vb;

  for ( iatom = llatom ; iatom < (signed)atom.size() ; iatom++ )
  {
    atom[iatom]->ref_bond = 0;
    atom[iatom]->bond = 0.0;
    atom[iatom]->ref_angle = 0;
    atom[iatom]->angle = 0.0;
    atom[iatom]->ref_dihedral = 0;
    atom[iatom]->dihedral = 0.0;

    if (iatom>0)
    {
      va = Vector(atom[iatom]->x - atom[iatom-1]->x,
                  atom[iatom]->y - atom[iatom-1]->y,
                  atom[iatom]->z - atom[iatom-1]->z);
      atom[iatom]->ref_bond = iatom - 1;
      atom[iatom]->bond = va.Norm();
      if (iatom>1)
      {
        va ^= 1.0;
        vb = Vector(atom[iatom-2]->x - atom[iatom-1]->x,
                    atom[iatom-2]->y - atom[iatom-1]->y,
                    atom[iatom-2]->z - atom[iatom-1]->z);
        vb ^= 1.0;
        atom[iatom]->ref_angle = iatom - 2;
        atom[iatom]->angle = RadianToDegree(acos(va.Dot(vb)));
        if (iatom>2) 
        {
          atom[iatom]->ref_dihedral = iatom - 3;
          val = Dihedral(iatom-3,iatom-2,iatom-1,iatom);
          atom[iatom]->dihedral = val;
        }
      }
    }
    if (atom[iatom]->angle==180.0) 
    {
      va = Vector(atom[iatom]->x - atom[iatom-1]->x,
                  atom[iatom]->y - atom[iatom-1]->y,
                  atom[iatom]->z - atom[iatom-1]->z);
      vb = va.Orthogonal();
      vb *= va.Norm();
      va = Vector(atom[iatom-1]->x,atom[iatom-1]->y,atom[iatom-1]->z);
      va += vb;
      std::vector<Atom*>::iterator it;
      it = atom.begin();
      it = atom.insert(it+iatom,new Atom(0,va.x,va.y,va.z));
      iatom--;
    }
  }
}
  
double Structure::Dihedral(int iatom, int jatom, int katom, int latom)
{
  double val;
  double ra[3],rb[3],rc[3],rd[3];
  Vector auxa,auxb,auxc,va,vb,vc,vd,rva,rvc,rvd;

  ra[0] = atom[iatom]->x;
  ra[1] = atom[iatom]->y;
  ra[2] = atom[iatom]->z;
  va = Vector(ra);

  rb[0] = atom[jatom]->x;
  rb[1] = atom[jatom]->y;
  rb[2] = atom[jatom]->z;
  vb = Vector(rb);

  rc[0] = atom[katom]->x;
  rc[1] = atom[katom]->y;
  rc[2] = atom[katom]->z;
  vc = Vector(rc);

  rd[0] = atom[latom]->x;
  rd[1] = atom[latom]->y;
  rd[2] = atom[latom]->z;
  vd = Vector(rd);

  rva = va - vb;
  rvc = vc - vb;
  rvd = vd - vc;
  rva ^= 1.0;
  rvc ^= 1.0;
  rvd ^= 1.0;

  if (rva.Colineal(rvc)) 
  {
    if (rvd.Dot(rvc)<0.0) return 0.0; 
    else return 180.0; 
  }

  auxa = rvc;
  auxa *= -rva.Dot( rvc );
  auxa += rva;
  auxa ^= 1.0;

  auxb = rvc;
  auxb *= -rvd.Dot( rvc );
  auxb += rvd;
  auxb ^= 1.0;

  auxc = auxa;
  auxc.Rotate( rvc ,-90.0 );
  val = auxa.Dot( auxb );

  // Avoid running out of domain of acos()
  val = double(int(100000.0*val+0.5))/100000.0; 

  if ( auxb.Dot( auxc ) < 0.0 )
  {
    val = 360.0 - 180.0/M_PI*acos( val );
  }
  else
  {
    val =  180.0/M_PI*acos( val );
  }
  return val;
}

/*
void Structure::WriteZMatrix(const char *filename)
{
  ofstream f(filename);
  f << atom.size() <<endl;
  f << endl;
  for (int i=0;i<(signed)atom.size();i++)
  {
    f << atom[i]->label << "  "; 
    if (i>0) 
      f << fixed << setw(4) << atom[i]->ref_bond << "  " 
        << fixed << setw(10) << setprecision(3) 
        << BohrToAngstrom(atom[i]->bond) << "   ";
    if (i>1) 
      f << fixed << setw(4) << atom[i]->ref_angle << "  "
        << fixed << setw(10) << setprecision(2) << atom[i]->angle << "   ";
    if (i>2) 
      f << fixed << setw(4) << atom[i]->ref_dihedral << "  "
        << fixed << setw(10) << setprecision(2) << atom[i]->dihedral;
    f << endl;
  }
  f.close();
};

*/


