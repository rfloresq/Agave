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
#ifndef AGV_AGAVE_H
#define AGV_AGAVE_H

#define SFB    0x00000010L  // Selecting for bond
#define SFA    0x00000020L  // Selecting for angle
#define SFD    0x00000040L  // Selecting for dihedral
#define SFRV   0x00000080L  // Selecting for distance value
#define SFAV   0x00000100L  // Selecting for angle value
#define SFDV   0x00000200L  // Selecting for diheral value
#define SFSURF 0x00000400L  // Selecting for surface
#define SFPG   0x00000800L  // Selecting for point group
#define SFFV   0x00001000L  // Selecting for fitting value
#define SFPB   0x00002000L         
#define QUEUED 0x00004000L  // GUI has been switched off     

#include <QtWidgets>

#include <Parameter.h>

using namespace std;

class Model;
class QM;
class System;
class Viewer;
class WinGeo;
class WinPlot;
class WinQM;

class Agave : public QMainWindow
{
  Q_OBJECT

  public:

    Agave( System*, QApplication*, int argc, char** argv );
    //~Agave();

    void Report( const char*, ... );

    char program_name[MAXSTR];
    System* sys;
    QApplication* app;
    QM *qm;
    WinGeo *wingeo;
    WinQM *winqm;
    Model* model;
    Viewer *viewer;
    WinPlot *winplot;
    // QMenuBar* menu;

    long status;

    void SetCursor( const char* );
    //void Error(const char*, int, const char*, int);

  public slots:
 
    void About( void );

    void ReadFile( void );
    void WriteFile( void );
    //void Run(void);

  protected:

  /*
    void CreateActions();
    void CreateMenus();
    */

    QStatusBar *statusbar;

    QMenu *fileMenu;
    QMenu *readMenu;
    QMenu *writeMenu;
    QMenu *helpMenu;
    //QMenu *calcMenu;

    vector<QAction*> RFMTA;
    vector<QAction*> WFMTA;
    //QAction *imgAct;
    QAction *quitAct;
    //QAction *runAct;
    QAction *aboutAct;
    /*
    QAction *exitAction;

    QToolButton *inButton;
    QToolButton *outButton;
    QToolButton *numAtomButton;

    QTimer* timer;
    QSpinBox *lineS;
    */

    QTabWidget *tabWidget;
    /*
    QTabWidget *tabWidgetv;

  protected slots:

    void WriteImage( void );
    */
};

#endif // AGV_AGAVE_H

