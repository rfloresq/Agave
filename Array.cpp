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
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <iomanip>
#include <fstream>

#include <Array.h>
#include <Matrix.h>

using namespace std;

Array::Array(int n)
{
  dim = n;

  values = new double[dim];
  SetZero();
}

double Array::Value( int i )
{
  if ((i<0)||(i>dim)) 
  {
    cout << "Attemp to extract non existing from Array"<<endl;
    cout << " i = "<<i<< " dim = "<<dim<<endl;
    return 0.0;
  }
  return values[i];
}

void Array::Transform(Matrix* m, Array *b)
{
  if (m->NCol()!=dim)
  {
    cout << "Array: Attemp to multiply with different dimension "<<endl;
    exit(EXIT_FAILURE);
  }
  m->VectorMultiply(values,b->values);
}

void Array::SetValue( int i, double val)
{
  values[i] = val;
}

void Array::ShiftValue( int i, double val)
{
  values[i] += val;
}

void Array::Fill(double val)
{
  int i;

  for (i=0;i<dim;i++)
    values[i] = val;
}

void Array::SetZero() { Fill(double(0.0)); }


double Array::Sum()
{
  int i;
  double t;

  t = 0.0;
  for (i=0;i<dim;i++)
    t += values[i];
  return t;
}

double Array::SumProduct(Array* other)
{
  int i;
  double t;

  t = 0.0;
  for (i=0;i<dim;i++)
    t += values[i]*other->Value(i);
  return t;
}

double Array::SumProduct(Array* other, int top)
{
  int i;
  double t;

  t = 0.0;
  for (i=0;i<top;i++)
    t += values[i]*other->Value(i);
  return t;
}

double Array::SumProduct(double* other)
{
  int i;
  double t;

  t = 0.0;
  for (i=0;i<dim;i++)
    t += values[i]*other[i];
  return t;
}

void Array::Copy(Array *other)
{
  int i;

  for (i=0;i<dim;i++)
    values[i] = other->Value(i);
}

void Array::Add(Array *other)
{
  for (int i=0;i<dim;i++)
    values[i] += other->values[i];
}

void Array::Add(double *other)
{
  for (int i=0;i<dim;i++)
    values[i] += other[i];
}

void Array::Scale(double s)
{
  for (int i=0;i<dim;i++)
    values[i] *= s;
}

