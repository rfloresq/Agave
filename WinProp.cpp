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
#include <iostream>

#include <QtWidgets>

/*
#include <ADFT.h>
#include <Curve.h>
#include <Math.h>
#include <Species.h>
#include <Spin.h>
#include <Structure.h>
#include <System.h>
#include <Agave.h>
#include <WinPlot.h>
*/
#include <GAPP.h>
#include <Global.h>
#include <QM.h>
#include <Units.h>
#include <WinProp.h>   

using namespace std;

WinProp::WinProp( Agave *ix , QM *iqm)
      : QWidget( 0 )
{
  agv = ix;
  qm = iqm;

  setWindowTitle( "Properties" );
  tabWidget = new QTabWidget;

  /*
  adptButton =  new QToolButton( this );  
  adptButton->setText( "ADPT" );
  connect( adptButton,SIGNAL(clicked()),this,SLOT(ADPT()));
  */

  propBox = new QComboBox;
  propBox->addItem(QString("Particle binding energies"));
  /*
  propBox->addItem(QString("H removal"));
  propBox->addItem(QString("Av. Polarizability"));
  propBox->addItem(QString("Nuclear Fukui function"));
  //propBox->addItem(QString("Energy"));
  */
  strcpy(property,(propBox->currentText().toLatin1()).data());
  connect(propBox, SIGNAL(activated(const QString &)),
          this, SLOT(ChangeProperty(const QString &)));

  theoryBox = new QComboBox;
  theoryBox->addItem(QString("G1PP-1"));
  theoryBox->addItem(QString("G1PP-2"));
  theoryBox->addItem(QString("G1PP-P3"));
  theoryBox->addItem(QString("G1PP-3"));
  theoryBox->addItem(QString("G1PP-P4"));
  /*
  theoryBox->addItem(QString("Linear Response"));
  theoryBox->addItem(QString("SCF"));
  //theoryBox->addItem(QString("MP2"));
  */
  strcpy(theory,(theoryBox->currentText().toLatin1()).data());
  connect(theoryBox, SIGNAL(activated(const QString &)),
          this, SLOT(ChangeTheory(const QString &)));

  QToolButton *evalButton =  new QToolButton( this );
  evalButton->setText( "Evaluate" ); 
  connect(evalButton,SIGNAL(clicked()), this,SLOT(Evaluate())); 

  typea = 0;
  QComboBox *specBox = new QComboBox;
  specBox->addItem(QString("Electron"));
  specBox->addItem(QString("Positron"));
  specBox->addItem(QString("Proton"));
  ChangeSpecies( specBox->currentText() );
  connect(specBox, SIGNAL(activated(const QString &)),
          this, SLOT(ChangeSpecies(const QString &)));
  typeb = 0;
  QComboBox *specBoxb = new QComboBox;
  specBoxb->addItem(QString("Electron"));
  specBoxb->addItem(QString("Positron"));
  specBoxb->addItem(QString("Proton"));
  ChangeSpeciesB( specBox->currentText() );
  connect(specBox, SIGNAL(activated(const QString &)),
          this, SLOT(ChangeSpeciesB(const QString &)));

  spina = 0;
  QComboBox *spinBox = new QComboBox;
  spinBox->addItem(QString("Alpha"));
  spinBox->addItem(QString("Beta"));
  ChangeSpin( spinBox->currentText() );
  connect(spinBox, SIGNAL(activated(const QString &)),
          this, SLOT(ChangeSpin(const QString &)));
  spinb = 0;
  QComboBox *spinBoxb = new QComboBox;
  spinBoxb->addItem(QString("Alpha"));
  spinBoxb->addItem(QString("Beta"));
  ChangeSpinB( spinBoxb->currentText() );
  connect(spinBox, SIGNAL(activated(const QString &)),
          this, SLOT(ChangeSpinB(const QString &)));

  ida = 0;
  QSpinBox *onS = new QSpinBox( this );
  onS->setRange( 1 , 1000 );
  onS->setSingleStep( 1 );
  onS->setValue( ida + 1 );
  connect( onS , SIGNAL(valueChanged( int )) ,
           this , SLOT( ChangeOrbitalNumber( int ) ) );
  idb = 0;
  QSpinBox *onSb = new QSpinBox( this );
  onSb->setRange( 1 , 1000 );
  onSb->setSingleStep( 1 );
  onSb->setValue( idb + 1 );
  connect( onSb , SIGNAL(valueChanged( int )) ,
           this , SLOT( ChangeOrbitalNumberB( int ) ) );

  QGridLayout *mainLayout = new QGridLayout;

  frame = new QWidget( 0 );
  tabWidget->addTab( frame , tr("Properties (a.u.)") );

  mainLayout->addWidget( tabWidget       , 1 , 0 , 6, 4 );
  mainLayout->addWidget( propBox         , 7 , 0 , 1, 1 );
  mainLayout->addWidget( theoryBox       , 7 , 1 , 1, 1 );
  mainLayout->addWidget( evalButton      , 7 , 2 , 1, 1 );
  mainLayout->addWidget( specBox         , 8 , 0 , 1, 1 );
  mainLayout->addWidget( spinBox         , 8 , 1 , 1, 1 );
  mainLayout->addWidget( onS             , 8 , 2 , 1, 1 );
  //mainLayout->addWidget( adptButton      , 8 , 3 , 1, 1 );
  mainLayout->addWidget( specBoxb        , 9 , 0 , 1, 1 );
  mainLayout->addWidget( spinBoxb        , 9 , 1 , 1, 1 );
  mainLayout->addWidget( onSb            , 9 , 2 , 1, 1 );

  setLayout( mainLayout );

  SCFListBox = new QListWidget;
  layout = new QVBoxLayout;
  layout->addWidget( SCFListBox );
  frame->setLayout( layout );
}

void WinProp::Evaluate()
{
  SCFListBox->clear();
  char str[MAXSTR];
  double value = 0.0;

  /*
  if ( strncmp(property,"Av. Polarizability",18) == 0 )
  {
    if ( strncmp(theory,"ADPT",4) != 0 ) 
    {
      cout << __FILE__<<__LINE__<<"Select ADPT for this property"<<endl;
      return;
    }
    Matrix alpha(3,3);
    qm->Polarizability(alpha);

    cout << "xx "<<alpha(0,0)<<endl;
    cout << "xy "<<alpha(0,1)<<endl;
    cout << "xz "<<alpha(0,2)<<endl;
    cout << "yx "<<alpha(1,0)<<endl;
    cout << "yy "<<alpha(1,1)<<endl;
    cout << "yz "<<alpha(1,2)<<endl;
    cout << "zx "<<alpha(2,0)<<endl;
    cout << "zy "<<alpha(2,1)<<endl;
    cout << "zz "<<alpha(2,2)<<endl;
    value = alpha.Trace()/3.0;
    sprintf( str,"%s(%s) = %20.8f (a.u.)\n", property,theory,value);
  }
  else if ( strncmp(property,"Particle binding",16) == 0 )
  */
  if ( strncmp(property,"Particle binding",16) == 0 )
  {
	  /*
    if ( strncmp(theory,"SCF",3) == 0 )
    {
      value = -qm->sys->species[typea]->spin[spina]->moe[ida];
    }
    else if ( strncmp(theory,"Linear Response",15) == 0 )  // ADPT
    {
      qm->ParticleBindingEnergy(typea,spina,ida,&value);
    }
    else if ((strncmp(theory,"G1PP",4) == 0 )|| 
             (strncmp(theory,"G2PP",4) == 0 )) // Propagator
	     */
    if ((strncmp(theory,"G1PP",4) == 0 )|| 
        (strncmp(theory,"G2PP",4) == 0 )) // Propagator
    {
      GAPP gapp(qm->sys);
      double corre;
      gapp.Diagonal(theory,typea,typeb,spina,spinb,ida,idb,&corre);
      value = -corre; // Sign corrected
    }
    value = HartreeToeV(value);
    sprintf( str,"%s(%s) = %20.3f eV\n", property,theory,value);
  }
  /*
  else if ( strncmp(property,"H removal",9) == 0 )
  {
    if ( strncmp(theory,"SCF",3) == 1 ) 
    {
      qm->HRemoval(&value,ida,1);
    }
    else if ( strncmp(theory,"Linear Response",15) == 0 ) 
    {
      qm->HRemoval(&value,ida,2);
    }
    value = HartreeToeV(value);
    sprintf( str,"%s(%s) = %20.3f eV\n", property,theory,value);
  }
  else if ( strncmp(property,"Nuclear Fukui function",22) == 0 )
  {
    if ( strncmp(theory,"Linear Response",15) == 0 ) 
    {
      qm->NuclearFukuiFunction(typea,spina,ida);
    } else cout << "This property only with linear response"<<endl;
    sprintf( str,"Computed, you can plot now\n");
  }
  */
  SCFListBox->insertItem( 0, QString(str) );

}

/*
void WinProp::ADPT()
{
  qm->PrepareForADPT();
}
*/

void WinProp::ChangeProperty( const QString &newp )
{
  strcpy(property,(newp.toLatin1()).data());
}

void WinProp::ChangeTheory( const QString &newt )
{
  strcpy(theory,(newt.toLatin1()).data());
}

void WinProp::ChangeSpecies( const QString &str )
{ SetSpecies(str,&typea); }
void WinProp::ChangeSpeciesB( const QString &str )
{ SetSpecies(str,&typeb); }
void WinProp::SetSpecies( const QString &str, int* specid )
{
  char name[MAXSTR];

  strcpy(name, (str.toLatin1()).data() );
  if ( strncasecmp(name,"Electron",8) == 0 ) *specid = ELECTRON;
  else if ( strncasecmp(name,"Positron",8) == 0 ) *specid = POSITRON;
  else if ( strncasecmp(name,"Proton",8) == 0 ) *specid = PROTON;
  else if ( strncasecmp(name,"Hydrogen",6) == 0 ) *specid = HYDROGEN;
}

void WinProp::ChangeSpin( const QString &str )
{ SetSpin(str,&spina); }
void WinProp::ChangeSpinB( const QString &str )
{ SetSpin(str,&spinb); }
void WinProp::SetSpin( const QString &str, int* spinid )
{
  char name[MAXSTR];

  strcpy(name, (str.toLatin1()).data() );
  if ( strncasecmp(name,"Alpha",5) == 0 ) *spinid = ALPHA;
  else if ( strncasecmp(name,"Beta",4) == 0 ) *spinid = BETA;
}

void WinProp::ChangeOrbitalNumber( int n )
{ SetOrbitalNumber(n-1,&ida); }
void WinProp::ChangeOrbitalNumberB( int n )
{ SetOrbitalNumber(n-1,&idb); }
void WinProp::SetOrbitalNumber( int n, int* orbid )
{
  *orbid = n;
}
