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
#ifndef AGV_WINPROP_H
#define AGV_WINPROP_H

#include <vector>

#include <QtGui>

#include <Parameter.h>
//Qt5 port
#include <QtWidgets/QToolButton>
#include <QComboBox>
#include <QVBoxLayout>
#include <QListWidget>
#include <QTabWidget>

using namespace std;

class Agave;
class QM;
class Curve;

class WinProp : public QWidget
{
  Q_OBJECT

  public:
    WinProp( Agave*, QM* );

    Agave *agv;
    QM *qm;

  public slots:

    //void ADPT(void);
    void ChangeSpecies(const QString &);
    void ChangeSpin( const QString &);
    void ChangeOrbitalNumber( int );
    void ChangeSpeciesB(const QString &);
    void ChangeSpinB( const QString &);
    void ChangeOrbitalNumberB( int );
    void ChangeProperty( const QString & );
    void ChangeTheory( const QString & );
    void Evaluate(void);

  protected:

    void SetSpecies(const QString &, int*);
    void SetSpin( const QString &,int*);
    void SetOrbitalNumber( int,int* );

    char theory[MAXSTR];
    char property[MAXSTR];

    int ida,idb;
    int spina,spinb;
    int typea,typeb;

    //QToolButton *adptButton;
    QComboBox *theoryBox;
    QComboBox *propBox;
    QTabWidget *tabWidget;
    QVBoxLayout *layout;
    QWidget* frame;
    QListWidget *SCFListBox;

};

#endif // AGV_WINPROP_H
