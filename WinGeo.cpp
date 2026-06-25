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
#include <math.h>

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <QtOpenGL/QGLWidget>
#include <QtWidgets>
#include <QStringList>

/*
#include <Vector.h>
*/
#include <Atom.h>
#include <Element.h>
#include <Agave.h>
#include <Structure.h>
#include <System.h>
#include <Units.h>
#include <Viewer.h>   
#include <WinElement.h>
#include <WinGeo.h>   

using namespace std;
using std::ifstream;


WinGeo::WinGeo( Agave *iagv )
              : QWidget( 0 )
{
  agv = iagv;
  mol = agv->sys->mol;
  setWindowTitle( "WinGeo Editor" );
  
  //QList<QString> elebas;
  //QList<QString> nucbas;

  QStringList Headers;
  Headers << "Atom" 
          << "Bond" << " " 
          << "Angle" << " " 
          << "Dihedral" << " " 
          << "e-Basis" << "Nuc-Basis";
  QGridLayout *mainLayout = new QGridLayout;

  tw = new QTableWidget( 100 , 9 );
  connect( tw , SIGNAL( itemChanged(QTableWidgetItem*) ), 
           this , SLOT( Edited( QTableWidgetItem*) ));
  tw->setHorizontalHeaderLabels(Headers);

  QToolButton *addButton =  new QToolButton( this );
  addButton->setText( "Add Atom" ); 
  connect(addButton,SIGNAL(clicked()),this,SLOT(AddAtomStart())); 

  QToolButton *delButton =  new QToolButton( this );
  delButton->setText( "Delete Last Atom" ); 
  connect(delButton,SIGNAL(clicked()),this,SLOT(DeleteAtom())); 

  /*
  QToolButton *fdButton =  new QToolButton( this );
  fdButton->setText( "Fit Distance" );
  connect(fdButton,SIGNAL(clicked()),this,SLOT(FitDistance()));
  */

  QComboBox *mComboBox = new QComboBox;
  mComboBox->addItem( QString("Measurement"));
  mComboBox->addItem( QString("Distance"));
  mComboBox->addItem( QString("Angle"));
  mComboBox->addItem( QString("Dihedral"));
  mComboBox->addItem( QString("Clean"));
  connect(mComboBox, SIGNAL(activated(const QString &)),
          this, SLOT( SetupMeasurement(const QString &)));

  /*
  QCheckBox* cb = new QCheckBox( this );
  cb->setText( "Def. Atom" );
  connect( cb , SIGNAL( clicked() ) , this , SLOT( UseDefault() ) );


  QCheckBox* cz = new QCheckBox( this );
  cz->setText( "Cartesians" ); 
  connect( cz , SIGNAL( clicked() ) , this , SLOT( ChangeCoordinates() ) );

  QCheckBox* cu = new QCheckBox( this );
  cu->setText( "Bohrs" ); 
  connect( cu , SIGNAL( clicked() ) , this , SLOT( ChangeUnits() ) );

  QComboBox *agComboBox = new QComboBox;
  agComboBox->addItem( QString("Add Group"));
  agComboBox->addItem( QString("--OH"));
  agComboBox->addItem( QString("--NH2"));
  agComboBox->addItem( QString("--CH3"));
  agComboBox->addItem( QString("File XYZ"));
  agComboBox->addItem( QString("File ZMT"));
  connect(agComboBox, SIGNAL(activated(const QString &)),
          this, SLOT( AddGroupStart(const QString &)));

	  */

  QComboBox *tagComboBox = new QComboBox;
  tagComboBox->addItem( "None");
  tagComboBox->addItem( "Symbols");
  tagComboBox->addItem( "Numbers");
  tagComboBox->addItem( "Sym+Num");
  tagComboBox->addItem( "Mulliken Charge");
  tagComboBox->addItem( "Mulliken Charge Response");
  tagComboBox->addItem( "Hirshfeld Charge");
  tagComboBox->addItem( "Hirshfeld Charge Response");
  connect(tagComboBox, SIGNAL(activated(const QString &)),
          agv->viewer, SLOT(ChangeTags(const QString &)));
  QLabel* lablab = new QLabel(tr("Labels:"));
  lablab->setBuddy(tagComboBox);

  /*
  QComboBox *pgComboBox = new QComboBox;
  pgComboBox->addItem( QString("Symmetry Group"));
  pgComboBox->addItem( QString("Td"));
  pgComboBox->addItem( QString("Oh"));
  pgComboBox->addItem( QString("Ih"));
  pgComboBox->addItem( QString("Cs"));
  pgComboBox->addItem( QString("Ci"));
  pgComboBox->addItem( QString("C1"));
  pgComboBox->addItem( QString("C2h"));
  pgComboBox->addItem( QString("C3h"));
  pgComboBox->addItem( QString("C4h"));
  pgComboBox->addItem( QString("C5h"));
  pgComboBox->addItem( QString("C2v"));
  pgComboBox->addItem( QString("C3v"));
  pgComboBox->addItem( QString("C4v"));
  pgComboBox->addItem( QString("C5v"));
  pgComboBox->addItem( QString("C6v"));
  pgComboBox->addItem( QString("C7v"));
  pgComboBox->addItem( QString("C8v"));
  pgComboBox->addItem( QString("C9v"));
  pgComboBox->addItem( QString("C*v"));
  pgComboBox->addItem( QString("S2"));
  pgComboBox->addItem( QString("S4"));
  pgComboBox->addItem( QString("S6"));
  pgComboBox->addItem( QString("S8"));
  pgComboBox->addItem( QString("D2h"));
  pgComboBox->addItem( QString("D3h"));
  pgComboBox->addItem( QString("D4h"));
  pgComboBox->addItem( QString("D5h"));
  pgComboBox->addItem( QString("D6h"));
  pgComboBox->addItem( QString("D7h"));
  pgComboBox->addItem( QString("D8h"));
  pgComboBox->addItem( QString("D9h"));
  pgComboBox->addItem( QString("D*h"));
  pgComboBox->addItem( QString("D3d"));
  pgComboBox->addItem( QString("D5d"));
  pgComboBox->addItem( QString("D7d"));
  pgComboBox->addItem( QString("D9d"));
  pgComboBox->addItem( QString("s1c1"));
  connect(pgComboBox, SIGNAL(activated(const QString &)),
          this, SLOT( BuildGroup(const QString &)));

  QSpinBox *risS = new QSpinBox( this );
  risS->setRange( 1 , 10000 );
  risS->setSingleStep( 1 );
  risS->setValue( 4 );
  QLabel* rislab = new QLabel(tr("RIS (alkane):"));
  rislab->setBuddy(risS);
  connect( risS , SIGNAL(valueChanged( int )) , this , SLOT( GoRIS( int ) ) );

  QSpinBox *tubeS = new QSpinBox( this );
  tubeS->setRange( 1 , 1000 );
  tubeS->setSingleStep( 1 );
  tubeS->setValue( 1 );
  QLabel* tubelab = new QLabel(tr("Nanotube:"));
  tubelab->setBuddy(tubeS);
  connect( tubeS , SIGNAL(valueChanged( int )) , this , SLOT( BuildTube(int)));
  */

  mainLayout->addWidget( addButton   , 0 , 0 , 1, 1 );
  mainLayout->addWidget( delButton   , 0 , 1 , 1, 1 );
  /*
  mainLayout->addWidget( cb          , 0 , 2 , 1, 1 );
  mainLayout->addWidget( rislab      , 0 , 3 , 1, 1 );
  mainLayout->addWidget( risS        , 0 , 4 , 1, 1 );
  mainLayout->addWidget( agComboBox  , 1 , 0 , 1, 1 );
  mainLayout->addWidget( pgComboBox  , 1 , 1 , 1, 1 );
  mainLayout->addWidget( cz          , 1 , 2 , 1, 1 );
  mainLayout->addWidget( tubelab     , 1 , 3 , 1, 1 );
  mainLayout->addWidget( tubeS       , 1 , 4 , 1, 1 );
  */
  //mainLayout->addWidget( fdButton    , 2 , 0 , 1, 1 );
  mainLayout->addWidget( mComboBox   , 2 , 1 , 1, 1 );
  /*
  mainLayout->addWidget( cu          , 2 , 2 , 1, 1 );
  */
  mainLayout->addWidget( lablab      , 2 , 3 , 1, 1 );
  mainLayout->addWidget( tagComboBox , 2 , 4 , 1, 1 );
  mainLayout->addWidget( tw          , 3 , 0 , 3, -1 );
  setLayout( mainLayout );

  CCSET = ZMatrix;
  default_element = 6;
  use_default = false;
  bohr = false;
  strcpy(group,"Atom");
}

void WinGeo::AddAtomStart( )
{
  if ( !use_default )
  {
    WinElement es(this,default_element);
    default_element = es.GetSelectedElement();
  };
  Enabled( false );
  if ( mol->Natom() < 2 )
  {
    AddAtom( mol->Natom(),0, 0 );
  }
  else 
  {
    if ( mol->Natom() < 3 )
    {
      agv->Report( "Select one atom" ); 
    }
    else if ( mol->Natom() < 4 )
    {
      agv->Report( "Select two atoms" );
    }
    else
    {
      agv->Report( "Select three atoms" );  
    };
    Enabled( false );
    agv->status |= SFB;
    agv->SetCursor( "Selection" );
  };
}

void WinGeo::Update()
{
  int iatom;
  Atom *a;

  disconnect( tw , SIGNAL( itemChanged(QTableWidgetItem*) ), 
             this , SLOT( Edited( QTableWidgetItem*) ));

  for ( iatom = 0 ; iatom < mol->Natom() ; iatom++ )
  {
    a = mol->GetAtom(iatom);
    if ( ! tw->item(iatom,0) )
    {
      tw->setItem(iatom,0, new QTableWidgetItem( " " ));
      tw->setItem(iatom,1, new QTableWidgetItem( " " ));
      tw->setItem(iatom,2, new QTableWidgetItem( " " ));
      tw->setItem(iatom,3, new QTableWidgetItem( " " ));
      tw->setItem(iatom,4, new QTableWidgetItem( " " ));
      tw->setItem(iatom,5, new QTableWidgetItem( " " ));
      tw->setItem(iatom,6, new QTableWidgetItem( " " ));
    }
    (tw->item(iatom,0))->setText(QString( a->label ));
    if (CCSET==ZMatrix) 
    {
      if ( iatom > 0 ) 
      {
        // Remember: users sees angstroms unless he asks otherwise, 
        // inside we have bohrs always
        (tw->item(iatom,1))->setText( QString::number( a->ref_bond ) );
        if (bohr) (tw->item(iatom,2))->setText( QString::number( a->bond ) );
        else (tw->item(iatom,2))->setText( QString::number( BohrToAngstrom(a->bond) ) );
        if ( iatom > 1 )
        {
          (tw->item(iatom,3))->setText( QString::number( a->ref_angle ) );
          (tw->item(iatom,4))->setText( QString::number( a->angle ) );
          if ( iatom > 2 )
          {
            (tw->item(iatom,5))->setText(QString::number( a->ref_dihedral ) );
            (tw->item(iatom,6))->setText(QString::number( a-> dihedral ) );
          }
          else
          {
            (tw->item( iatom , 5 ))->setText( " " ); 
            (tw->item( iatom , 6 ))->setText( " " );
          }
        }
        else
        {
          (tw->item( iatom , 3 ))->setText( " " ); 
          (tw->item( iatom , 4 ))->setText( " " );
          (tw->item( iatom , 5 ))->setText( " " );
          (tw->item( iatom , 6 ))->setText( " " );
        }
      }
      else
      {
        (tw->item( iatom , 1 ))->setText( " " ); 
        (tw->item( iatom , 2 ))->setText( " " );
        (tw->item( iatom , 3 ))->setText( " " );
        (tw->item( iatom , 4 ))->setText( " " );
        (tw->item( iatom , 5 ))->setText( " " );
        (tw->item( iatom , 6 ))->setText( " " );
      }
    }
    else
    { 
      (tw->item( iatom , 1 ))->setText( QString( "x" ) );
      (tw->item( iatom , 3 ))->setText( QString( "y" ) );
      (tw->item( iatom , 5 ))->setText( QString( "z" ) );
      // Remember: users sees angstroms unless he asks otherwise, 
      // inside we have bohrs always
      if (bohr) 
      {
        (tw->item( iatom , 2 ))->setText( QString::number( a->x ) );
        (tw->item( iatom , 4 ))->setText( QString::number( a->y ) );
        (tw->item( iatom , 6 ))->setText( QString::number( a->z ) );
      }
      else
      {
        (tw->item(iatom,2))->setText( QString::number(BohrToAngstrom(a->x)) );
        (tw->item(iatom,4))->setText( QString::number(BohrToAngstrom(a->y)) );
        (tw->item(iatom,6))->setText( QString::number(BohrToAngstrom(a->z)) );
      }
    }
  }

  connect( tw , SIGNAL( itemChanged(QTableWidgetItem*) ), 
           this , SLOT( Edited( QTableWidgetItem*) ));

  Enabled( true );
}

void WinGeo::DeleteAtom()
{
  int natom = mol->Natom();
  if ( natom > 0 )
  {
    mol->DeleteLastAtom();
    natom--;
    (tw->item( natom , 0 ))->setText( QString( " " ) );
    (tw->item( natom , 1 ))->setText( QString( " " ) );
    (tw->item( natom , 2 ))->setText( QString( " " ) );
    (tw->item( natom , 3 ))->setText( QString( " " ) );
    (tw->item( natom , 4 ))->setText( QString( " " ) );
    (tw->item( natom , 5 ))->setText( QString( " " ) );
    (tw->item( natom , 6 ))->setText( QString( " " ) );
    agv->viewer->Redraw();
    Enabled( true );
  };
}

void WinGeo::Enabled( bool enabled )
{
  setEnabled( enabled );
}

void WinGeo::Edited( QTableWidgetItem* item ) 
{
  int iatom = tw->currentRow();
  int col = tw->currentColumn();

  if ( col == 0 )
  {
    mol->atom[iatom]->SetSymbol( (item->text().toLatin1()).data() );
    item->setText( QString( ELEMENT_SYMBOL[mol->atom[iatom]->atomic_number] ) );
  }
  else if ( col == 1 )
  {
    item->setText( QString::number( mol->atom[iatom]->ref_bond ) );
  }
  else if ( col == 2 )
  {
    if (CCSET==ZMatrix) 
    {
      if (bohr) mol->atom[iatom]->bond = item->text().toDouble();
      else mol->atom[iatom]->bond = AngstromToBohr( item->text().toDouble() );
      mol->Z2C(iatom);
    }
    else 
    {
      if (bohr) mol->atom[iatom]->x = item->text().toDouble();
      else mol->atom[iatom]->x = AngstromToBohr( item->text().toDouble() );
      mol->C2Z(iatom);
    }
  }
  else if ( col == 3 )
  {
    item->setText( QString::number( mol->atom[iatom]->ref_angle ) );
  }
  else if ( col == 4 )
  {
    if (CCSET==ZMatrix) 
    {
      mol->atom[iatom]->angle = item->text().toDouble();
      mol->Z2C(iatom);
    }
    else 
    {
      if (bohr) mol->atom[iatom]->y = item->text().toDouble();
      else mol->atom[iatom]->y = AngstromToBohr( item->text().toDouble() );
      mol->C2Z(iatom);
    }
  }
  else if ( col == 5 )
  {
    item->setText( QString::number( mol->atom[iatom]->ref_dihedral ) );
  }
  else if ( col == 6 )
  {
    if (CCSET==ZMatrix) 
    {
      mol->atom[iatom]->dihedral = item->text().toDouble();
      mol->Z2C(iatom);
    }
    else 
    {
      if (bohr) mol->atom[iatom]->z = item->text().toDouble();
      else mol->atom[iatom]->z = AngstromToBohr( item->text().toDouble() );
      mol->C2Z(iatom);
    }
  } 
  agv->viewer->Redraw();
}

/*
void WinGeo::BuildGroup( const QString &s )
{
  if ( s == "Symmetry Group" ) return;
  sprintf( gname , "%s", (s.toLatin1()).data() );
  agv->status |= SFPG;
  Enabled ( false );
}

void WinGeo::AddGroupStart(const QString &s)
{
  bool hold;

  strcpy(group, (s.toLatin1()).data() );

  hold = use_default;
  use_default = true;

  AddAtomStart();

  use_default = hold;
}


void WinGeo::ChangeCoordinates()
{
  if (CCSET==Cartesians)
  {
    CCSET = ZMatrix;
  }
  else
  {
    CCSET = Cartesians;
  }
  Update();
}

void WinGeo::UseDefault()
{
  use_default =  (use_default ? false : true );
}

void WinGeo::FitDistance( void )
{
  agv->SetCursor( "Selection" );
  agv->status |= SFFV;
  agv->Report("Select two atoms");
}


void WinGeo::GoRIS(int n)
{
  //if (n!=100) return;

  double gr,r;

  // Initialize random number generator
  Random *random = new Random();

  gr = 0.0;
  int top = 0;
  while (top<1)
  {
    r = mol->GrowRIS( n , random);
    if (r>=0.0) 
    {
      cout << "Gyration radius = " << r <<endl;
      gr += r;
      top++;
    }
  }
  gr /= double(top);
  cout << "AVERAGE Gyration radius = " << gr <<endl;
  //return;
  Update();
  agv->viewer->Redraw();
}

void WinGeo::BuildTube(int n)
{
  int i,k,nc;
  double angle;
  double bond_length,h,l,r;

  nc = 9;
  bond_length = AngstromToBohr(1.3);

  h = bond_length*sin(M_PI/6.0);
  l = sqrt(4.0*(bond_length*bond_length-h*h));
  r = l/(2.0*sin(M_PI/nc));
  angle = 360.0/nc;
  
  Vector z(0.0,0.0,1.0);
  Vector pos(r,0.0,0.0);
  NATOM = 0;
  for (k=0;k<n;k++)
  {
    for (i=0;i<nc;i++)
    {
      ATOMNO[NATOM] = default_element;
      COORD[1][NATOM][0] = pos[0];
      COORD[1][NATOM][1] = pos[1];
      COORD[1][NATOM][2] = pos[2];
      pos.Rotate(z,angle);
      NATOM++;
    }
    angle /= 2.0;
    pos.Rotate(z,angle);
    angle *= 2.0;
    pos += Vector(0.0,0.0,h);
    for (i=0;i<nc;i++)
    {
      ATOMNO[NATOM] = default_element;
      COORD[1][NATOM][0] = pos[0];
      COORD[1][NATOM][1] = pos[1];
      COORD[1][NATOM][2] = pos[2];
      pos.Rotate(z,angle);
      NATOM++;
    }
    pos += Vector(0.0,0.0,bond_length);
  }
  center();
}

int WinGeo::SGDrv(int* ref) 
{
  if ( gname[0] == 's' )
  {
    return mol->SGBDrv(gname,ref);
  }
  else 
  {
    return mol->PGBDrv(gname,ref);
  };
};

void WinGeo::ChangeUnits()
{
  if (bohr)
  {
    bohr = false;
  }
  else
  {
    bohr = true;
  }
  Update();
}

// Reading from files is mediated by windows
// (File substitutes user, not window)
void WinGeo::Read(char* filename, char* fmt)
{
  if ( strcasecmp( fmt , "XYZ" ) == 0 ) mol->ReadXYZ(filename);
  else if ( strcasecmp( fmt , "ZMT" ) == 0 ) mol->ReadZMatrix(filename);
  Update(); 
  agv->viewer->Redraw();
}
*/


void WinGeo::SetupMeasurement( const QString &s )
{
  if ( s == "Measurement" ) return;
  if ( s == "Distance" && mol->Natom() > 1 ) 
  {
    agv->SetCursor( "Selection" );
    agv->status |= SFRV;
    agv->Report("Select two atoms");
  }
  else if ( s == "Angle" && mol->Natom() > 2)
  {
    agv->SetCursor( "Selection" );
    agv->status |= SFAV;
    agv->Report("Select three atoms");
  }
  else if ( s == "Dihedral" && mol->Natom() > 3 )
  {
    agv->SetCursor( "Selection" );
    agv->status |= SFDV;
    agv->Report("Select four atoms");
  }
  else if ( s == "Clean" && agv->viewer ) 
  {
    agv->viewer->monitored_bonds.clear();
    agv->viewer->monitored_angles.clear();
    agv->viewer->monitored_dihedrals.clear();
    agv->viewer->Redraw();
  };
}

void WinGeo::AddAtom( int na, int nb, int nc )
{
  double bond;

  if ( strncmp( group , "Atom", 4) == 0 )
  {
    if (na>0) bond = ELEMENT_COV_R[default_element] + 
                   ELEMENT_COV_R[mol->atom[na-1]->atomic_number];
    else bond = 0.0;
    mol->AddAtom( default_element, na - 1, nb - 1, nc - 1,
                         bond,109.467, 180.0);
  }
  else
  {
    //AddGroup(na-1,nb-1,nc-1);
    //strcpy( group , "Atom");
  }
  Update();
  agv->viewer->Redraw();
}

/*
void WinGeo::AddGroup( int n1, int n2, int n3 )
{
  int an,na,nb,nc;
  double bond;

  if ( strncmp(group,"--OH",4) == 0 ||
       strncmp(group,"--NH",4) == 0 ||
       strncmp(group,"--CH",4) == 0 ) 
  {   
    if ( strncmp(group,"--OH",4) == 0 ) an = 8;
    else if ( strncmp(group,"--NH",4) == 0 ) an = 7;
    else if ( strncmp(group,"--CH",4) == 0 ) an = 6;

    na = n1;
    nb = n2;
    nc = n3;
    mol->atom[na]->SetAtomicNumber( an );
    bond = ELEMENT_COV_R[an] + ELEMENT_COV_R[1];
    mol->AddAtom(1,na,nb,nc,bond,109.467,180.0);
    if ( strncmp(group,"--OH",4) != 0 ) 
    {
      mol->AddAtom(1,na,nb,nc,bond,109.467,300.0);
      if ( strncmp(group,"--CH",4) == 0 ) 
        mol->AddAtom(1,na,nb,nc,bond,109.467,60.0);
    }
  }
  else if ( strncmp(group,"File XYZ",8) == 0 || 
            strncmp(group,"File ZMT",8) == 0 )
  {
    cout << "LOADING OF FILE UNDER CONSTRUCTION"<<endl;
    char filename[MAXSTR];
    int iatom,patom;
    double angle,dihedral;

    patom = mol->Natom() - 1;

    Structure mol2 = Structure();

    if (!agv->GetFileName(filename))
      cout << "Unable to open file "<<endl;

    if ( strncmp(group,"File XYZ",8) == 0 ) mol2.ReadXYZ(filename);
    else if ( strncmp(group,"File ZMT",8) == 0 ) mol2.ReadZMatrix(filename);

    for (iatom=0;iatom<mol2.Natom();iatom++)
    {
      an = mol2.atom[iatom]->atomic_number;
      if (iatom==0)  
      {
        mol->atom[n1]->SetAtomicNumber( an );
      }
      else
      {
        na = patom + mol2.atom[iatom]->ref_bond;
        if (na==patom) na = n1;
        bond = mol2.atom[iatom]->bond;

        if (iatom<2)
        {
          nb = n2;
          angle = 109.467;
        }
        else
        {
          nb = patom + mol2.atom[iatom]->ref_angle;
          angle = mol2.atom[iatom]->angle;
        }
        if (nb==patom) nb = n1;

        if (iatom<3)
        {
          nc = n3;
          dihedral = 60.0 + 120.0*iatom;
        }
        else
        {
          nc = patom + mol2.atom[iatom]->ref_dihedral;
          dihedral = mol2.atom[iatom]->dihedral;
        }
        if (nc==patom) nc = n1;

        mol->AddAtom(an,na,nb,nc,bond,angle,dihedral);
      }
    }
  }
}
*/

