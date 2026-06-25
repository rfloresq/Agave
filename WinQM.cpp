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
#include <fstream>
#include <iomanip>

#include <Global.h>   
#include <QM.h>   
#include <Agave.h>   
#include <System.h>   
#include <WinProp.h>
#include <WinQM.h>   
#include <WinSCF.h>
#include <WinSpecies.h>
/*
#include <Species.h>   
#include <WinMS.h>
#include <WinLFT.h>
#include <WinRSGF.h>
#include <WinRMN.h>
#include <WinRT.h>
#include <WinSSQ.h>
#include <WinTDDFT.h>
*/

using namespace std;

WinQM::WinQM(Agave* ix,System* isys)
     : QWidget( 0 )
{
  agv = ix;
  sys = isys;
  qm = agv->qm;

  QGridLayout *mainLayout = new QGridLayout;

  tabSpecies = new QTabWidget;
  tabWidget = new QTabWidget;

  mainLayout->addWidget( tabSpecies     ,      1 , 0 , 1 , 4);
  mainLayout->addWidget( tabWidget    ,      6 , 0 , 1 , 4);

  setLayout( mainLayout );

  tabSpecies->clear();
  winelectron = new WinSpecies(qm,ELECTRON);
  tabSpecies->addTab( winelectron , QString("ELECTRON") );
  winpositron = new WinSpecies(qm,POSITRON);
  tabSpecies->addTab( winpositron , QString("POSITRON") );
  winhydrogen = new WinSpecies(qm,HYDROGEN);
  tabSpecies->addTab( winhydrogen , QString("HYDROGEN") );
  winproton = new WinSpecies(qm,PROTON);
  tabSpecies->addTab( winproton , QString("PROTON") );

  tabWidget->clear();
  winscf = new WinSCF( agv , qm );
  tabWidget->addTab( winscf , QString("SCF") );
  winprop = new WinProp( agv , qm );
  tabWidget->addTab( winprop , QString("Properties") );
  /*
  winssq = new WinSSQ( agv , qm );
  tabWidget->addTab( winssq , QString("SSQ") );
  winrt = new WinRT( agv , qm );
  tabWidget->addTab( winrt , QString("RT-DFT") );
  winrsgf = new WinRSGF( agv , qm );
  tabWidget->addTab( winrsgf , QString("RSGF") );
  winlft = new WinLFT( agv , qm );
  tabWidget->addTab( winlft , QString("LFT") );
  wintddft = new WinTDDFT( agv , qm );
  tabWidget->addTab( wintddft , QString("TD-DFT") );
  winms = new WinMS( agv , qm );
  tabWidget->addTab( winms , QString("Mass Spectrometry") );
  winrmn = new WinRMN( agv , qm );
  tabWidget->addTab( winrmn , QString("RMN") );
  */

  setWindowTitle( "Run calculation" );
  hide();
}



