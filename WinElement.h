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
#ifndef AGV_WINELEMENT_H
#define AGV_WINELEMENT_H

#include <Parameter.h>

#include "qdialog.h"

class QToolButton;

class WinElement : public QDialog 
{
  Q_OBJECT

  public:
    WinElement( QWidget *parent = 0 , int guess=0 );

  public slots:
    void SetDefaultElement( void );
    int GetSelectedElement( void );

  protected:
    QToolButton *element[120]; 
    int selected_element;

};

#endif // AGV_WINELEMENT_H

