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
#ifndef AGV_CURVE_H
#define AGV_CURVE_H

#include <vector>

#include <QtGui>
//Qt5 port
#include <QtWidgets/QToolBox>
using namespace std;

class Curve : public QWidget
{
  Q_OBJECT

  public:

    enum Shape { Line, Points, Polyline, Polygon, 
                 Rect, RoundRect, Ellipse, Arc,
                 Chord, Pie, Path, Text };

    Curve(QWidget *parent = 0);

    QSize minimumSizeHint() const;
    QSize sizeHint() const;
    void SetData( double* , double* , int );
    void SetBounds( double , double , double, double );
    void Redraw(void);

  public slots:

    void setShape(Shape shape);
    void setPen(const QPen &pen);
    void setBrush(const QBrush &brush);
    void setAntialiased(bool antialiased);
    void setTransformed(bool transformed);
    void SaveImage( void );
    void HighlightPoint( int );

  protected:

    void paintEvent(QPaintEvent *event);
    void DrawStrings( void );

    vector<double> x;
    vector<double> y;

  private:

    Shape shape;
    QPen pen;
    QBrush brush;
    bool antialiased;
    bool transformed;
    double xmin;
    double ymin;
    double xmax;
    double ymax;
    int sel;
};

#endif  // AGV_CURVE_H

