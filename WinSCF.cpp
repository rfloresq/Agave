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

#include <Curve.h>
#include <Math.h>
#include <Parameter.h>   
#include <QM.h>   
#include <Agave.h>
#include <SCF.h>   
#include <Structure.h>
#include <System.h>
#include <Viewer.h>
#include <WinSCF.h>   
#include <WinPlot.h>

using namespace std;

WinSCF::WinSCF( Agave *ix , QM *iqm)
      : QWidget( 0 )
{
  agv = ix;
  qm = iqm;
  timer = new QTimer(this);

  setWindowTitle( "SCF" );
  tabWidget = new QTabWidget;

  setButton =  new QToolButton( this );  
  setButton->setText( "Setup" );
  connect( setButton,SIGNAL(clicked()),this,SLOT(Setup()));

  stepButton =  new QToolButton( this );  
  stepButton->setText( "Step" );
  connect( stepButton,SIGNAL(clicked()),this,SLOT(Step()));

  QGridLayout *mainLayout = new QGridLayout;

  QToolButton *startButton =  new QToolButton( this );
  startButton->setText( "Start" ); 
  connect(startButton,SIGNAL(clicked()), this,SLOT(Start())); 

  QToolButton *stopButton =  new QToolButton( this );
  stopButton->setText( "Stop" ); 
  connect(stopButton,SIGNAL(clicked()), this,SLOT(Stop())); 

  graphical = false;
  gbox = new QCheckBox( this );
  gbox->setText( QString("Graphical SCF") ); 
  connect( gbox , SIGNAL( clicked() ) , this , SLOT( SetGraphical() ) );

  methodBox = new QComboBox;
  methodBox->addItem(QString("Dirac-VWN"));
  methodBox->addItem(QString("Hartree-Fock"));
  methodBox->addItem(QString("Hartree-None"));
  methodBox->addItem(QString("Dirac-None"));
  methodBox->addItem(QString("None-None"));
  ChangeMethod( methodBox->currentText() );
  connect(methodBox, SIGNAL(activated(const QString &)),
          this, SLOT(ChangeMethod(const QString &)));

  frame = new QWidget( 0 );
  tabWidget->addTab( frame , tr("SCF Energy") );

  SCFPlotter = new Curve( this );
  SCFPlotter->resize( 300 , 200 ); 

  mainLayout->addWidget( tabWidget         , 1 , 0 , 5, 5 );
  mainLayout->addWidget( setButton         , 6 , 0 , 1, 1 );
  mainLayout->addWidget( stepButton        , 6 , 1 , 1, 1 );
  mainLayout->addWidget( startButton       , 6 , 2 , 1, 1 );
  mainLayout->addWidget( stopButton        , 6 , 3 , 1, 1 );
  mainLayout->addWidget( gbox              , 6 , 4 , 1, 1 );
  mainLayout->addWidget( methodBox         , 6 , 5 , 1, 1 );
  mainLayout->addWidget( SCFPlotter        , 9 , 0 , 6, 5 );
  setLayout( mainLayout );

  SCFListBox = new QListWidget;
  layout = new QVBoxLayout;
  layout->addWidget( SCFListBox );
  frame->setLayout( layout );
}

void WinSCF::Setup()
{
  qm->PrepareForSCF();
  energies.clear(); 
  SCFListBox->clear();
}

void WinSCF::Step()
{
  qm->scf->Step();

  energies.push_back( qm->scf->Energy() );
  int ncycle = energies.size();

  int i;
  char str[MAXSTR];
  i = ncycle - 1;
  sprintf( str,"%d: energy = %20.8f , error = %20.8f\n", 
      i, energies[i], qm->scf->Error());
  SCFListBox->insertItem( 0, QString(str) );

  int ll;
  if (ncycle>20) ll = ncycle - 20;
  else ll = 0;
  double *cval;
  cval = new double[20];
  for (i=ll;i<ncycle;i++)
    cval[i-ll] = double(i);
  double *eval;
  eval = new double[20];
  for (i=ll;i<ncycle;i++)
    eval[i-ll] = energies[i];
  SCFPlotter->SetData( cval , eval ,  ncycle-ll );
  delete[] cval;
  delete[] eval;
  if (graphical) 
  {
    agv->winplot->SetupPlot();
    agv->winplot->BuildPlot();
  }
  agv->viewer->Redraw();
  timer->start(40);
}

void WinSCF::Start()
{
  connect(timer,SIGNAL(timeout()),this,SLOT(Step()));
  Step();
}

void WinSCF::Stop()
{
  disconnect(timer,SIGNAL(timeout()), this,SLOT(Step()));
}

void WinSCF::SetGraphical( )
{ 
  if (graphical) graphical = false;
  else graphical = true;
}

void WinSCF::ChangeMethod( const QString &newmethod )
{
  qm->scf->ChangeMethod( (newmethod.toLatin1()).data() );
}

