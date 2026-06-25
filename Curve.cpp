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

#include <Curve.h>

static int offset = 2;

Curve::Curve( QWidget* wdg )
     : QWidget( wdg )
{
  shape = Polygon;
  antialiased = false;

  sel = -1;

  x.clear();
  y.clear();
}

QSize Curve::minimumSizeHint() const
{
  return QSize(100, 100);
}

QSize Curve::sizeHint() const
{
  return QSize(400, 200);
}

void Curve::setShape(Shape shape)
{
  this->shape = shape;
  update();
}

void Curve::setPen(const QPen &pen)
{
    this->pen = pen;
    update();
}

void Curve::setBrush(const QBrush &brush)
{
    this->brush = brush;
    update();
}

void Curve::setAntialiased(bool antialiased)
{
    this->antialiased = antialiased;
    update();
}

void Curve::setTransformed(bool transformed)
{
    this->transformed = transformed;
    update();
}

void Curve::Redraw()
{
  QPaintEvent *e;
  paintEvent(e);
}

void Curve::paintEvent(QPaintEvent *)
{
  int x1,x2,y1,y2,h,w;
  QPainter painter(this);
  painter.setPen( QPen(QColor( 255, 255, 255)));
  painter.setBrush(brush);
  if (antialiased)
    painter.setRenderHint(QPainter::Antialiasing);
  painter.setBackground( QBrush( QColor( 255, 255, 255, 0) ) );

  w = width();
  h = height();
  painter.fillRect( 0, 0 , w , h , QBrush(QColor( 85, 0, 127)));

  double scale = 0.9;
  //painter.drawLine( offset , h-offset , offset, offset );
  //painter.drawLine( offset , h-offset , w-offset, h-offset );
  for ( size_t i = 1 ; i < x.size() ; i++ )
  {
    x1 = (int)( ((x[i-1]-xmin)/(xmax-xmin))*w*scale);
    x2 = (int)( ((x[i]-xmin)/(xmax-xmin))*w*scale);
    y1 = (int)( ((y[i-1]-ymin)/(ymax-ymin))*h*scale);
    y2 = (int)( ((y[i]-ymin)/(ymax-ymin))*h*scale);
    painter.drawLine( x1+2*offset , h-y1-2*offset , x2+2*offset, h-y2-2*offset );
  }
  for ( size_t i = 0 ; i < x.size() ; i++ )
  {
    x2 = (int)( ((x[i]-xmin)/(xmax-xmin))*w*scale);
    y2 = (int)( ((y[i]-ymin)/(ymax-ymin))*h*scale) + 1;
    painter.drawRect( x2+2*offset, h-y2-2*offset , 2 , 2 );
  }
  if ( sel >= 0 ) 
  {
    x2 = (int)( ((x[sel]-xmin)/(xmax-xmin))*w*scale);
    y2 = (int)( ((y[sel]-ymin)/(ymax-ymin))*h*scale) + 1;
    painter.fillRect( x2+2*offset-4, h-y2-2*offset-4, 12 , 12 , 
                      QBrush(QColor( 255, 170, 128)));
  }
}

void Curve::SetData( double* xx, double *yy, int n)
{
  x.clear();
  y.clear();
  for ( int i = 0 ; i < n ; i++ )
  {
    x.push_back( xx[i] );
    y.push_back( yy[i] );
  }
  if ( x.size() == 0 || y.size() == 0) return;

  xmin = x[0];
  ymin = y[0];
  xmax = x[0];
  ymax = y[0];
  for ( int i = 1 ; i < n ; i++ )
  {
    if ( x[i] < xmin ) xmin = x[i];
    if ( x[i] > xmax ) xmax = x[i];
    if ( y[i] < ymin ) ymin = y[i];
    if ( y[i] > ymax ) ymax = y[i];
  }
}

void Curve::SetBounds( double xxmin, double xxmax ,
                              double yymin, double yymax )
{
  xmin = xxmin;
  xmax = xxmax;
  ymin = yymin;
  ymax = yymax;
}

void Curve::SaveImage()
{
}

void Curve::HighlightPoint(int n)
{
  if ( n < 0 || n >= x.size()) return;
  sel = n;
  update();
}
