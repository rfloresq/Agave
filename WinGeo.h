/*
This file is part of Agave software.

    AgavAgave is free software: you can redistribute it and/or modify
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
#ifndef AGV_WINGEO_H
#define AGV_WINGEO_H

#include <vector>

#include <QtGui>
#include <QStringList>

#include <Parameter.h>
//Qt5 port
#include <QTableWidgetItem>

#define Cartesians 0
#define ZMatrix    1

#include <Parameter.h>

using namespace std;

class Agave;
class Structure;

class WinGeo : public QWidget
{
  Q_OBJECT

  public:
    WinGeo( Agave* );

    Agave *agv;
    Structure *mol;
    int default_element;
    bool use_default;
    //char gname[80];

    void Enabled( bool );
    void AddAtom( int, int, int );
    /*
    void AddGroup(int,int,int);
    void Read(char*,char*);
    */

  public slots:

    void AddAtomStart( void );
    //void AddGroupStart(const QString &);
    void DeleteAtom( void );
    void Edited(QTableWidgetItem*);
    void SetupMeasurement(const QString &);
    /*
    void BuildGroup( const QString & );    
    void ChangeCoordinates( void );
    void ChangeUnits( void );
    */
    void Update( void );
    //void FitDistance( void );
    /*
    void UseDefault( void );
    void GoRIS(int);
    void BuildTube(int);
    int SGDrv(int*);
    */

  protected:

    bool bohr;
    char group[MAXSTR];
    int CCSET;
    QTableWidget *tw;
    QStringList *Headers;
};

#endif // AGV_WINGEO_H
