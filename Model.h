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
#ifndef AGV_MODEL_H
#define AGV_MODEL_H

#include <GL/glu.h>

#include <QtCore>

class Agave;
class Atom;
class Structure;

class Model
{
  public:

    Model(Agave*);

    Agave* agv;
    Structure *mol;
    double srs;
    double crs;
    bool drawmol;
    GLUquadricObj* q;

    void DrawSphere( double* , double, int );
    void DrawCylinder( double* , double* , double, int );
    /*
    void DrawCone( double* , double* , double, int );
    */
    void DrawMolecule( int, bool );
    // void ChangeKind( const QString & );

  protected:

    void DrawAtom( Atom* , double, int );
    void DrawBond( Atom* , Atom* , double, int );
};

#endif // AGV_MODEL_H

