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
#ifndef AGV_WINPLOT_H
#define AGV_WINPLOT_H

#include <vector>

#include <QWidget>

#include <Vector.h>

using namespace std;

class QListWidget;
class QTabWidget;
class QVBoxLayout;
class QDoubleSpinBox;
class QCheckBox;
class QRadioButton;

class Plotter;
class Agave;
class Viewer;
class Vector;

class WinPlot : public QWidget
{
  Q_OBJECT

  public:
    WinPlot( Agave* , Viewer* );

    Agave *agv;
    Viewer *viewer;
    Plotter *plotter;

    vector<double> temp;

  public slots:

    void Adjust(void);
    double AdjustBox(void);
    void ChangeIso( double );
    void ChangeLineWidth( int );
    void ChangePointSize( int );
    void ChangeBox( double );
    void ChangeMesh( double );
    void ChangeFlags( void );  
    void ChangeSpecies(const QString &);
    void ChangeSpin( const QString &);
    void ChangeOrbitalNumber( int );
    void SetupPlot(void);
    void BuildPlot(void);

  protected:

    void NewList( char** , int , char* ); 
    void SetupBox( QVBoxLayout* );

    int hold;
    int hide_surfaces;
    int drawref;
    int approximate_normals;
    int average_normals;

    int shiny;
    int line_width;
    int typea;
    int spina;
    int orba;
    int point_size;
    int lastplaneaxis;
    double isovalue;

    QDoubleSpinBox *cminS;
    QDoubleSpinBox *cmaxS;
    QDoubleSpinBox *cstepS;
    QDoubleSpinBox *isoS;
    QDoubleSpinBox *meshS;

    QDoubleSpinBox* axis[4][3];
    QCheckBox* flagbox[4][3];
    QTabWidget *tabWidget;
    vector<QWidget*> tabs;
    vector<QListWidget*> listados;
    vector<double> grdval;
    vector<double> grdvaldegen;
    vector<Vector> grdvec;
    vector<Vector> grdvecdegen;
};

#endif  // AGV_WINPLOT_H

