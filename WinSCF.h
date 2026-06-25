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
#ifndef AGV_WINSCF_H
#define AGV_WINSCF_H

#include <vector>

#include <QtGui>
//Qt5 port
#include <QListWidget>
#include <QVBoxLayout>
#include <QToolButton>
#include <QCheckBox>
using namespace std;

class QComboBox;

class Agave;
class QM;
class Curve;

class WinSCF : public QWidget
{
  Q_OBJECT

  public:
    WinSCF( Agave*, QM* );

    Agave *agv;
    QM *qm;

    void Enabled( bool );

  public slots:

    void Setup(void);
    void Step(void);
    void Start(void);
    void Stop(void);
    void SetGraphical(void);
    void ChangeMethod( const QString & );

  protected:

    bool graphical;
    double error;
    vector<double> energies;
    QTabWidget *tabWidget;
    QTimer* timer;
    QVBoxLayout *layout;
    QWidget* frame;
    QListWidget *SCFListBox;
    Curve* SCFPlotter;
    QComboBox *methodBox;

    QToolButton *setButton;
    QToolButton *stepButton;
    QCheckBox* gbox;
};

#endif // AGV_WINSCF_H
