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
#include <fstream>
#include <iomanip>

/*
#include <Species.h>   
#include <SCF.h>   
#include <Spin.h>   
*/
#include <Global.h>   
#include <QM.h>   
#include <Structure.h>   
#include <System.h>   
#include <WinSpecies.h>   

using namespace std;

WinSpecies::WinSpecies(QM* iqm, int itype)
        : QWidget( 0 )
{
  qm = iqm;
  type = itype;

  guessBox = new QComboBox;
  guessBox->addItem(QString("Core"));
  if (type==HYDROGEN) guessBox->addItem(QString("Hydrogen"));
  ChangeGuess( guessBox->currentText() );
  connect(guessBox, SIGNAL(activated(const QString &)),
          this, SLOT(ChangeGuess(const QString &)));

  basisBox = new QComboBox;
  if (type==ELECTRON)
  {
    basisBox->addItem(QString("DZVP"));
    basisBox->addItem(QString("DZV"));
    basisBox->addItem(QString("STO-3G"));
    basisBox->addItem(QString("4-31G"));
    basisBox->addItem(QString("6-311G"));
    basisBox->addItem(QString("CC-PVDZ"));
    basisBox->addItem(QString("CC-PVTZ"));
    basisBox->addItem(QString("CC-PVQZ"));
    basisBox->addItem(QString("CC-PV5Z"));
    basisBox->addItem(QString("AUG-CC-PVTZ"));
  }
  else if (type==HYDROGEN)
  {
    basisBox->addItem(QString("H+DZSPDN"));
    basisBox->addItem(QString("CC-PVDZ"));
    basisBox->addItem(QString("CC-PVTZ"));
    basisBox->addItem(QString("CC-PVQZ"));
    basisBox->addItem(QString("CC-PV5Z"));
  }
  else if (type==PROTON)
  {
    basisBox->addItem(QString("H+DZSPDN"));
    basisBox->addItem(QString("CC-PVDZ"));
    basisBox->addItem(QString("CC-PVTZ"));
    basisBox->addItem(QString("CC-PVQZ"));
    basisBox->addItem(QString("CC-PV5Z"));
  }
  else if (type==POSITRON)
  {
    basisBox->addItem(QString("PSX-TZ"));
    basisBox->addItem(QString("ROMERO"));
    basisBox->addItem(QString("DZVP"));
  }

  QLabel* basislab = new QLabel(tr("Basis:"));
  basislab->setBuddy(basisBox);
  ChangeBasis( basisBox->currentText() );
  connect(basisBox, SIGNAL(activated(const QString &)),
          this, SLOT(ChangeBasis(const QString &)));
  /*

  QLabel* pplab;
  if (type==ELECTRON)
  {
    ppBox = new QComboBox;
    ppBox->addItem(QString("NONE"));
    ppBox->addItem(QString("ECP|SD"));
    ppBox->addItem(QString("QECP|SD"));
    ppBox->addItem(QString("RECP|SD"));
    ppBox->addItem(QString("ECP|LANL2DZ"));
    ppBox->addItem(QString("QECP|LANL2DZ"));
    ppBox->addItem(QString("MCP|LK"));
    ppBox->addItem(QString("RMCP|LK"));

    pplab = new QLabel(tr("Pseudopotentials:"));
    pplab->setBuddy(ppBox);
    ChangePPs( ppBox->currentText() );
    connect(ppBox, SIGNAL(activated(const QString &)),
            this, SLOT(ChangePPs(const QString &)));
  }
  */

  alphaSpin = new QDoubleSpinBox( this );
  alphaSpin->setRange( 0.0 , 100000.0 );
  alphaSpin->setSingleStep( 1.0 );
  alphaSpin->setDecimals( 0 );
  alphaSpin->setValue( 0.0 );
  QLabel* alphalab = new QLabel(tr("No. alpha:"));
  alphalab->setBuddy(alphaSpin);
  ChangeAlpha( alphaSpin->value() );
  connect( alphaSpin , SIGNAL(valueChanged( double )) ,
           this, SLOT( ChangeAlpha( double ) ) );

  betaSpin = new QDoubleSpinBox( this );
  betaSpin->setRange( 0.0 , 100000.0 );
  betaSpin->setSingleStep( 1.0 );
  betaSpin->setDecimals( 0 );
  betaSpin->setValue( 0.0 );
  QLabel* betalab = new QLabel(tr("No. beta:"));
  betalab->setBuddy(betaSpin);
  ChangeBeta( betaSpin->value() );
  connect( betaSpin , SIGNAL(valueChanged( double )) ,
           this, SLOT( ChangeBeta( double ) ) );

  QGridLayout *mainLayout = new QGridLayout;

  if (type==ELECTRON)
  {
    contButton =  new QToolButton( this );  
    contButton->setText( "Count" );
    connect( contButton,SIGNAL(clicked()),this,SLOT(Default()));
  }

  mainLayout->addWidget( basislab    ,      1 , 0 , 1 , 1);
  mainLayout->addWidget( basisBox    ,      1 , 1 , 1 , 1);
  /*
  if (type==ELECTRON)
  {
    mainLayout->addWidget( pplab     ,      1 , 2 , 1 , 1);
    mainLayout->addWidget( ppBox     ,      1 , 3 , 1 , 1);
  }
  */
  mainLayout->addWidget( alphalab    ,      2 , 0 , 1 , 1);
  mainLayout->addWidget( alphaSpin   ,      2 , 1 , 1 , 1);
  mainLayout->addWidget( betalab     ,      2 , 2 , 1 , 1);
  mainLayout->addWidget( betaSpin    ,      2 , 3 , 1 , 1);
  if (type==ELECTRON)
    mainLayout->addWidget( contButton   ,     4 , 0 , 1 , 1);
  mainLayout->addWidget( guessBox     ,     4 , 3 , 1 , 1);

  setLayout( mainLayout );
}

void WinSpecies::Default()
{
  // Count electrons
  int nalpha,nbeta,npart;
 
  if (type==ELECTRON) npart = qm->sys->mol->NumberOfElectronsForNeutral();
  else npart = 0;

  nbeta = npart/2;
  nalpha = npart - nbeta;

  qm->ChangeNumber( double(nalpha) , type, ALPHA); 
  qm->ChangeNumber( double(nbeta) , type, BETA);
  
  alphaSpin->setValue( double(nalpha) );
  betaSpin->setValue( double(nbeta) );
}

void WinSpecies::ChangeBasis( const QString &newbas )
{
  qm->ChangeBasis( (newbas.toLatin1()).data(), type );
}
/*
void WinSpecies::ChangePPs( const QString &newpp )
{
  qm->ChangePPs( (newpp.toLatin1()).data(), type );
}
*/

void WinSpecies::ChangeGuess( const QString &newguess )
{
  qm->ChangeGuess( (newguess.toLatin1()).data(), type );
}

void WinSpecies::ChangeAlpha( double n )
{
  qm->ChangeNumber( n , type , ALPHA);
}

void WinSpecies::ChangeBeta( double n )
{
  qm->ChangeNumber( n , type, BETA );
}

