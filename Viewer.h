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
#ifndef AGV_VIEWER_H
#define AGV_VIEWER_H

#include <vector>
#include <QWheelEvent>
#include <QtOpenGL/QtOpenGL>

#include <Vector.h>
//Qt5 porting
#include <QtPrintSupport/QPrinter>
#include <QtPrintSupport/QPrintDialog>

using namespace std;

class wheelEvent;

class Plotter;
class Agave;
class System;
//class WinImage;

class Viewer : public QGLWidget
{
  Q_OBJECT

  public:
    Viewer( Agave*, System* );

    Agave *agv;
    System *sys;
    Plotter* plotter;
    //WinImage *winimg;
    double trans[3];
    double scale;
    GLdouble rot[16];
    vector<vector<size_t> > monitored_bonds;
    vector<vector<size_t> > monitored_angles;
    vector<vector<size_t> > monitored_dihedrals;
    int SelectedAtom(size_t);

    void DrawStrings( void );
    void DrawSurfaces( void );
    void PickAtom( void );
    void PickSurface( void );
    void Project( void );
    void Rotate( double , double , double , double );
    void ResetRotation(void);
    void Write(char*,Vector,size_t);
    void ZoomIn( void );
    void ZoomOut( void );

  public slots:
 
    void ChangeTags(const QString&);
    void Redraw(void);
    void SaveSelection(int);
    void wheelEvent(QWheelEvent *event);

  protected:

    double tag_color[4];
    QString tag_type;
    /*
    QString arrow_type;
    QPrinter *printer;
    */
    int mouse[2];
    vector<size_t> selected_atoms;
    float view_width;
    float view_height;
    float view_far;
    float view_near;

    void initializeGL( void );
    void resizeGL( int , int );
    void paintGL( void );

    double bg[4];

  protected slots:

    void keyPressEvent( QKeyEvent* );
    void mousePressEvent( QMouseEvent* );
    void mouseMoveEvent( QMouseEvent* );
    void mouseReleaseEvent( QMouseEvent* );
};

#endif  // AGV_VIEWER_H

