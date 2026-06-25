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
#ifndef AGV_WINSPECIES_H
#define AGV_WINSPECIES_H

#include <QtGui>

#include <Parameter.h>
//Qt5 port
#include <QtWidgets/QToolButton>
#include <QComboBox>
#include <QLabel>
#include <QDoubleSpinBox>
#include <QGridLayout>
#include <QCheckBox>

class QComboBox;
class QCheckBox;
class QToolButton;
class QDoubleSpinBox;

class QM;

class WinSpecies : public QWidget 
{
  Q_OBJECT

  public:

    WinSpecies(QM*,int);

    QM *qm;

  public slots:
 
    void ChangeAlpha(double);
    void ChangeBeta(double);
    void Default(void);
    void ChangeBasis( const QString & );
    //void ChangePPs( const QString & );
    void ChangeGuess( const QString & );

  protected:

    int type;

    QComboBox *guessBox;
    QComboBox *basisBox;
    //QComboBox *ppBox;


    QToolButton *contButton;
    /*
    QToolButton *setButton;
    QToolButton *stepButton;
    */

    QDoubleSpinBox *alphaSpin;
    QDoubleSpinBox *betaSpin;

    /*
    QCheckBox* dirbox;
    QCheckBox* projbox;
    QCheckBox* resbox;
    */
};

#endif // AGV_WINSPECIES_H

