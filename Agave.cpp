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

#include <QMainWindow>

#include <Model.h>
#include <System.h>   
#include <Parameter.h>   
#include <QM.h>   
#include <Agave.h>   
#include <Viewer.h>
#include <WinGeo.h>
#include <WinPlot.h>
#include <WinQM.h>


static const char *READ_FORMATS[] = {"XYZ","AGV","D2K",NULL};
static const char *WRITE_FORMATS[] = {"XYZ","AGV",NULL};

Agave::Agave( System* isys, QApplication *iapp, int argc, char** argv )
     : QMainWindow( 0 )
{
  sys = isys;
  qm = new QM(sys);

  status = 0;

  app = iapp;

  tabWidget = new QTabWidget;

  model = new Model( this );
  viewer = new Viewer( this, sys );
  wingeo = new WinGeo( this ); 
  winplot = new WinPlot( this , viewer );
  winqm = new WinQM(this,sys);

  statusbar = new QStatusBar( this );

  int i;

  // Configure reading
  i = 0;
  while ( READ_FORMATS[i] != NULL )
  {
    RFMTA.push_back( new QAction( tr( READ_FORMATS[i] ), this ) );
    RFMTA[i]->setStatusTip(tr("Open file in %1 format")
                           .arg(READ_FORMATS[i]));
    RFMTA[i]->setCheckable( true );
    connect( RFMTA[i], SIGNAL(triggered()), this, SLOT( ReadFile() ));
    i++;
  };

  // Configure writing
  i = 0;
  while ( WRITE_FORMATS[i] != NULL )
  {
    WFMTA.push_back( new QAction( tr( WRITE_FORMATS[i] ), this ) );
    WFMTA[i]->setStatusTip(tr("Open file in %1 format")
                           .arg(WRITE_FORMATS[i]));
    WFMTA[i]->setCheckable( true );
    connect( WFMTA[i], SIGNAL(triggered()), this, SLOT( WriteFile() ));
    i++;
  };

  /*
  imgAct = new QAction(tr("&Image"), this);
  imgAct->setStatusTip(tr("Write image file"));
  connect(imgAct, SIGNAL(triggered()), this, SLOT(WriteImage()));
  */

  quitAct = new QAction(tr("&Quit"), this);
  quitAct->setShortcut(tr("Ctrl+Q"));
  quitAct->setStatusTip(tr("Finish the application"));
  connect(quitAct, SIGNAL(triggered()),
          app, SLOT(closeAllWindows()));

  /*
  runAct = new QAction(tr("&Run"), this);
  runAct->setStatusTip(tr("Run calculation"));
  connect(runAct, SIGNAL(triggered()), this, SLOT(Run()));
  */

  aboutAct = new QAction(tr("&About"), this);
  aboutAct->setStatusTip(tr("Show the application's About box"));
  connect(aboutAct, SIGNAL(triggered()), this, SLOT(About()));

  fileMenu = new QMenu(tr("&File"),this);

  readMenu = new QMenu(tr("&Read"),this);
  writeMenu = new QMenu(tr("&Write"),this);


  for( size_t i = 0 ; i < RFMTA.size() ; i++ )
    readMenu->addAction( RFMTA[i] );
  fileMenu->addMenu(readMenu);

  for( size_t i = 0 ; i < WFMTA.size() ; i++ )
    writeMenu->addAction( WFMTA[i] );
  fileMenu->addMenu(writeMenu);

  fileMenu->addSeparator();

  /*
  fileMenu->addAction( imgAct );
  fileMenu->addSeparator();
  */

  fileMenu->addAction( quitAct );

  helpMenu = new QMenu(tr("&Help"),this);
  helpMenu->addAction(aboutAct);

  //calcMenu = new QMenu(tr("&Calculation"),this);
  //calcMenu->addAction(runAct);


  menuBar()->addMenu(fileMenu);
  menuBar()->addMenu(helpMenu);

  QWidget *visualizer = new QWidget(0);
  //QVBoxLayout *visuallayout = new QVBoxLayout;

  QGridLayout *mainLayout = new QGridLayout;
  visualizer->setLayout( mainLayout );


  //tabWidget->addTab( viewer , QString("Display") );
  tabWidget->addTab( wingeo , QString("Builder") );
  tabWidget->addTab( winplot , QString("Plotter") );
  tabWidget->addTab( winqm , QString("Calculator") );

  mainLayout->addWidget( viewer               , 1 , 1 , 10, 10 );
  mainLayout->addWidget( tabWidget            , 4 , 11 , 1, 10 );
  mainLayout->addWidget( statusbar            , 20, 11 ,1, 10 );

  setCentralWidget( visualizer );

  //timer = new QTimer(this);

  strcpy(program_name,"Agave 1");
  setWindowTitle(program_name);
  move( 0 , 0  );
  resize(2000 , 600);
}

/*
Agave::~Agave()
{
}
*/


void Agave::Report( const char *fmt, ... )
{
  va_list args;   
  char    s[MAXSTR];

  va_start( args , fmt );
  vsprintf(s , fmt , args );
  va_end( args );
 
  if (statusbar) statusbar->showMessage( QString( s ) );
}

void Agave::About()
{
  char str[MAXSTR];

  sprintf(str,"%s\n\nCopyright 2025 Roberto Flores-Moreno\nroberto.floresmoreno.qt@gmail.com",program_name);
  QMessageBox::about( this, "About Agave",str);
}

void Agave::ReadFile()
{
  char fileformat[MAXSTR];
  size_t i;

  for ( i = 0 ; i < RFMTA.size() ; i++ )
  {
    if ( RFMTA[i]->isChecked() ) 
    {
      RFMTA[i]->setChecked( false );
      sprintf(fileformat,"%s",((RFMTA[i]->text()).toStdString()).c_str());
      break;
    };
  };

  QString filter("All Files (*.*)");
  filter +=      "\nOnly    (*.xyz)";
  QString fn = QFileDialog::getOpenFileName(this,
                 QString("Read file"),
                 QDir::currentPath(),
                 filter);
  if ( ! ( fn.isEmpty() || fn.isNull() ) )
  {
    char filename[MAXSTR];

    sprintf(filename,"%s",(fn.toStdString()).c_str());
    Report( "Reading file %s ...", filename );
    string tit = program_name;
    tit += ":";
    tit += filename;
    setWindowTitle( tit.c_str() );
    sys->ReadFile(filename,fileformat);
    // qm->ReadFile(filename);
    viewer->Redraw();
    wingeo->Update();
    Report( "Reading file %s ... DONE", filename );
  };
}

void Agave::WriteFile()
{
  char sffmt[MAXSTR];
  char fileformat[MAXSTR];
  QString ffmt;
  size_t i;

  for ( i = 0 ; i < WFMTA.size() ; i++ )
  {
    if ( WFMTA[i]->isChecked() )
    {
      WFMTA[i]->setChecked( false );
      sprintf(sffmt,"%s",((WFMTA[i]->text()).toLatin1()).data());
      strcpy(fileformat,sffmt);
      // File extension comes only from first 3 characters
      sffmt[3] = '\0';
      ffmt = sffmt;
      break;
    };
  };
  QString initialPath = QDir::currentPath() + "/untitled.";
  initialPath = initialPath + ffmt.toLower();
  QString fn = QFileDialog::getSaveFileName( this,
               tr( "Write %1 file").arg(ffmt.toUpper()),
               initialPath, tr("All Files (*.%1)").arg(ffmt.toLower()) );
  if ( ! ( fn.isEmpty() || fn.isNull() ) )
  {
    char filename[MAXSTR];

    sprintf(filename,"%s",(fn.toLatin1()).data());
    Report( "Writing file %s ...", filename );
    sys->WriteFile(filename,fileformat);
    // qm->WriteFile(filename);
    Report( "Writing file %s ... DONE", filename );
  };
}

void Agave::SetCursor( const char *type )
{
  QString typestr = type;
  QCursor cursor;
  if ( typestr == "Normal" )
  {
    cursor = QCursor( Qt::ArrowCursor );
  }
  else if ( typestr == "Wait" )
  {
    cursor = QCursor( Qt::WaitCursor );
  }
  else if ( typestr == "Selection" )
  {
    cursor = QCursor( Qt::PointingHandCursor );
  };
  setCursor( cursor );
} 

/*
void Agave::Error( const char* filename , int linenumber, const char*msg , int el )
{
  QString title = "Attention";
  QString text;
  if ( el )
  {
    text = "Error occurred in ";
  }
  else
  {
    text = "Warning received from "; 
  };
  text += filename;
  text += ", line number ";
  text += QString::number( linenumber );
  text += "\n";
  text += msg;
  text += "\n";

  if ( el )
  {
    (void)QMessageBox::critical( 0, title, text, QMessageBox::Ok,
                                 QMessageBox::NoButton, QMessageBox::NoButton);
  }
  else
  {
    (void)QMessageBox::warning( 0, title, text, QMessageBox::Ok,
                                 QMessageBox::NoButton, QMessageBox::NoButton);
  };

  if ( el ) app->closeAllWindows();
}

void Agave::WriteImage()
{ 
  viewer->winimg->show(); 
}

void Agave::Run()
{
  sys->WriteFile("Agave.inp","AGV");
  qm->ReadFile("Agave.inp");
  app->closeAllWindows();
}
*/
