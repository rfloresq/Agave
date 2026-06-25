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
#ifndef AGV_WINQM_H
#define AGV_WINQM_H

#include <QtGui>

#include <Parameter.h>
//Qt5 port
#include <QtWidgets/QToolButton>
/*
class QComboBox;
class QCheckBox;
class QToolButton;
class QDoubleSpinBox;
*/
class QTabWidget;

class System;
class Agave;
/*
class WinLFT;
class WinMS;
class WinRSGF;
class WinRMN;
class WinTDDFT;
class WinRT;
class WinSSQ;
*/
class QM;
class WinProp;
class WinSCF;
class WinSpecies;

class WinQM : public QWidget 
{
  Q_OBJECT

  public:

    WinQM(Agave*,System*);

    Agave *agv;
    /*
    WinLFT *winlft;
    WinMS *winms;
    WinRMN *winrmn;
    WinRSGF *winrsgf;
    WinTDDFT *wintddft;
    */
    WinSCF *winscf;
    WinProp *winprop;
    /*
    WinSSQ *winssq;
    WinRT *winrt;
    */
    WinSpecies *winelectron;
    WinSpecies *winhydrogen;
    WinSpecies *winpositron;
    WinSpecies *winproton;
    /*
    //WinSpecies *winm;
    //WinSpecies *wino;
    */
    QM *qm;

  public slots:
 
  protected:

    System* sys;

    // QToolButton *stepButton;

    QTabWidget *tabWidget;
    QTabWidget *tabSpecies;
};

#endif // AGV_WINQM_H

