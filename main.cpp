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

#include <QtWidgets>
#include <System.h>
#include <Agave.h>

using namespace std;

int main(int argc, char **argv)
{
  // System
  System* sys = new System();

  // GUI
  QApplication app(argc, argv);
  /* rfm
  QCoreApplication::setOrganizationName("Agave developer");
  QCoreApplication::setApplicationName("Agave");
  QCoreApplication::setApplicationVersion(QT_VERSION_STR);

  QCommandLineParser parser;
  parser.setApplicationDescription(QCoreApplication::applicationName());
  parser.addHelpOption();
  parser.addVersionOption();
  parser.addPositionalArgument("file", "The file to open.");
  parser.process(app);
  */

  Agave agv(sys,&app,argc,argv);
  agv.show();
  app.exec();

  return 0;
}
