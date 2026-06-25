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
#include <lapacke.h>

#include <iostream>
#include <iomanip>
#include <fstream>

#include <Matrix.h>
#include <Math.h>

using namespace std;

Matrix::Matrix(int n,int m)
{
  dim_row = n;
  dim_col = m;

  values = new double[dim_row*dim_col];
  SetZero();
}

double & Matrix::operator () ( int i, int j )
{
  return values[i*dim_col+j];
}

void Matrix::operator += ( Matrix other )
{
  int i,j,k;

  k = 0;
  for (i=0;i<dim_row;i++)
  {
    for (j=0;j<dim_col;j++)
    {
      values[k] += other(i,j);
      k++;
    }
  }
}


void Matrix::Print(const char* filename, const char* header, int w, int nd)
{

  ofstream f(filename,ios::app);
  f << "=== "<<header<<" ==="<<endl;

  int nblk = 75/w;
  int ll,ul;
  ll = 0;
  while (ll<dim_col)
  {
    ul = MIN(dim_col-1,ll+nblk-1);
    f << "     ";
    for (int j=ll;j<=ul;j++)
      f << setw(w) << setprecision(nd) << fixed << j;
    f << endl;
    for (int i=0;i<dim_row;i++)
    {
      f << setw(4) << i << " ";
      for (int j=ll;j<=ul;j++)
        f << setw(w) << setprecision(nd) << fixed << values[i*dim_col+j];
      f << endl;
    }
    f << endl;
    ll = ul + 1;
  }
  f.close();
}
void Matrix::PrintOut(const char* header, int w, int nd)
{

  cout << "=== "<<header<<" ==="<<endl;

  int nblk = 75/w;
  int ll,ul;
  ll = 0;
  while (ll<dim_col)
  {
    ul = MIN(dim_col-1,ll+nblk-1);
    cout << "     ";
    for (int j=ll;j<=ul;j++)
      cout << setw(w) << setprecision(nd) << fixed << j;
    cout << endl;
    for (int i=0;i<dim_row;i++)
    {
      cout << setw(4) << i << " ";
      for (int j=ll;j<=ul;j++)
        cout << setw(w) << setprecision(nd) << fixed << values[i*dim_col+j];
      cout << endl;
    }
    cout << endl;
    ll = ul + 1;
  }
}

// Note: Eigenvectors are stored as rows to facilitate C handling
void Matrix::Diagonalize(double* eig)
{
  lapack_int info;

  info = LAPACKE_dsyev(LAPACK_COL_MAJOR,'V','U',dim_row,values,dim_col,eig);
  if (info>0) exit(EXIT_FAILURE);

  Transpose();
}

int Matrix::NRow()
{
  return dim_row;
}
int Matrix::NCol()
{
  return dim_col;
}

void Matrix::GetRowValues( int id, double *val )
{
  int ll = id*dim_col;
  for (int j=0;j<dim_col;j++)
    val[j] = values[ll+j];
}
void Matrix::GetColValues( int id, double *val )
{
  for (int j=0;j<dim_row;j++)
    val[j] = values[j*dim_col+id];
}

void Matrix::SetValue( int i, int j, double val)
{
  values[i*dim_col+j] = val;
}

void Matrix::ShiftValue( int i, int j, double val)
{
  values[i*dim_col+j] += val;
}

// Notice that transposition switches dimensionality
// Roberto Flores-Moreno, Ago 2015
void Matrix::Transpose(void)
{
  int i,j;
  int dim_row_h,dim_col_h;
  double *h;

  dim_row_h = dim_col;
  dim_col_h = dim_row;
  h = new double[dim_col_h*dim_row_h];

  for (i=0;i<dim_row;i++)
    for (j=0;j<dim_col;j++)
      h[j*dim_col_h+i] = values[i*dim_col+j];

  dim_row = dim_row_h;
  dim_col = dim_col_h;
  for (i=0;i<dim_row;i++)
    for (j=0;j<dim_col;j++)
      values[i*dim_col+j] = h[i*dim_col+j];

  delete[] h;
}

void Matrix::Multiply(Matrix *other,Matrix *res,bool right)
{
  int i,j,k;
  double s;

  if (right)
  { 
    if (dim_col != other->NRow() ) 
    {
      cout << "WARNING: Matrix dimensions are diferent" <<endl;
    }
    double others[dim_col];
    for (i=0;i<dim_row;i++)
    {
      for (j=0;j<other->NCol();j++)
      {
        other->GetColValues(j,others);
        s = 0.0;
        for (k=0;k<dim_col;k++)
        {
          s += values[i*dim_col+k]*others[k];
        } 
        res->SetValue(i,j,s);
      }
    }
  }
  else
  {
    if (other->NCol() != dim_row ) 
    {
      cout << "WARNING: Matrix dimensions are diferent" <<endl;
    }
    double others[dim_row];
    for (i=0;i<other->NRow();i++)
    {
      other->GetRowValues(i,others);
      for (j=0;j<dim_col;j++)
      {
        s = 0.0;
        for (k=0;k<dim_row;k++)
        {
          s += others[k]*values[k*dim_col+j];
        } 
        res->SetValue(i,j,s);
      }
    }
  }
}

void Matrix::Fill(double val)
{
  int i,j,k;

  k = 0;
  for (i=0;i<dim_row;i++)
  {
    for (j=0;j<dim_col;j++)
    {
      values[k] = val;
      k++;
    }
  }
}
void Matrix::FillRange(double val,int llrow,int llcol,int ulrow,int ulcol)
{
  int i,j,k;

  k = 0;
  for (i=0;i<dim_row;i++)
  {
    for (j=0;j<dim_col;j++)
    {
      if ((i>=llrow)&&(i<=ulrow)&&(j>=llcol)&&(j<=ulcol))
        values[k] = val;
      k++;
    }
  }
}

void Matrix::SetZero() { Fill(0.0); }

void Matrix::Symmetrize(void)
{
  int i,j;

  for (i=0;i<dim_row;i++)
    for (j=i+1;j<dim_col;j++)
      values[j*dim_col+i] = values[i*dim_col+j];
}

double Matrix::QTrace(Matrix *other)
{
  int i,j,ll;
  double t;
  double others[dim_col];

  if ((other->NRow()!=dim_row)||(other->NCol()!=dim_col))
  {
    cout << "Attemp to trace matrices of different dimensions "<<endl;
    exit(EXIT_FAILURE);
  }

  t = 0.0;
  for (i=0;i<dim_row;i++)
  {
    other->GetRowValues(i,others);
    ll = i*dim_col;
    for (j=0;j<dim_col;j++)
    {
      t += others[j]*values[ll+j];
    }
  }
  return t;
}

void Matrix::Copy(Matrix *other)
{
  int i,j,ll;
  double others[dim_col];

  for (i=0;i<dim_row;i++)
  {
    other->GetRowValues(i,others);
    ll = i*dim_col;
    for (j=0;j<dim_col;j++)
    {
      values[ll+j] = others[j];
    }
  }
}

void Matrix::Copy(double **other)
{
  int i,j,ll;

  for (i=0;i<dim_row;i++)
  {
    ll = i*dim_col;
    for (j=0;j<dim_col;j++)
    {
      values[ll+j] = other[i][j];
    }
  }
}

// Get negative power using singular value decomposition
void Matrix::SVDPower(double *eig,double x,double tol)
{
  int i,j;
  double w;

  Matrix tmp(dim_row,dim_col);
  tmp.Copy(this);

  tmp.Diagonalize(eig);

  for (i=0;i<dim_col;i++)
  {
    if (eig[i]<=tol)
    {
      eig[i] = 0.0;
      for (j=0;j<dim_row;j++)
        tmp(j,i) = 0.0;
    }
    else
    {
      w = pow(sqrt(eig[i]),ABS(x));
      if (x<0) w = 1.0/w;
      for (j=0;j<dim_row;j++)
        tmp(j,i) *= w;
    }
  }
  Matrix tmp2(dim_row,dim_col);
  tmp2.Copy(&tmp);
  tmp2.Transpose();
  tmp.Multiply(&tmp2,this,true); 
}

void Matrix::VectorMultiply(double *in,double *out)
{
  int i,ii,k;
  double s;

  for (i=0;i<dim_row;i++)
  {
    ii = i*dim_col;
    s = 0.0;
    for (k=0;k<dim_col;k++)
      s += values[ii+k]*in[k];
    out[i] = s;
  }
}

void Matrix::Add(Matrix *other)
{
  int i,j,ll;
  double others[other->NCol()];

  for (i=0;i<other->NRow();i++)
  {
    other->GetRowValues(i,others);
    ll = i*dim_col;
    for (j=0;j<other->NCol();j++)
    {
      values[ll+j] += others[j];
    }
  }
}

void Matrix::ExtractFrom(Matrix *other, int llrow, int llcol)
{
  int i,j,k,ll,llo;

  int odc = other->NCol();

  double others[odc];

  k = 0;
  for (i=0;i<dim_row;i++)
  {
    other->GetRowValues(i+llrow,others);
    for (j=0;j<dim_col;j++)
    {
      values[k] = others[j+llcol];
      k++;
    }
  }
}
void Matrix::InsertTo(Matrix *other, int llrow, int llcol)
{
  int i,j,k,ll,llo;

  int odc = other->NCol();

  double others[odc];

  k = 0;
  for (i=0;i<dim_row;i++)
  {
    other->GetRowValues(i+llrow,others);
    for (j=0;j<dim_col;j++)
    {
      other->SetValue(i+llrow,j+llcol,values[k]);
      k++;
    }
  }
}


void Matrix::Scale(double s)
{
  int i,j,ll;

  for (i=0;i<dim_row;i++)
  {
    ll = i*dim_col;
    for (j=0;j<dim_col;j++)
    {
      values[ll+j] *= s;
    }
  }
}


double Matrix::Trace()
{
  int i;
  double t;

  t = 0.0;
  for (i=0;i<dim_row;i++)
    t += values[i*dim_col+i];
  return t;
}

double Matrix::Sum()
{
  int i,j,ll;
  double suma;

  suma = 0.0;
  for (i=0;i<dim_row;i++)
  {
    ll = i*dim_col;
    for (j=0;j<dim_col;j++)
    {
      suma += values[ll+j];
    }
  }
  return suma;
}

double Matrix::AbsoluteSum()
{
  int i,j,ll;
  double suma;

  suma = 0.0;
  for (i=0;i<dim_row;i++)
  {
    ll = i*dim_col;
    for (j=0;j<dim_col;j++)
    {
      suma += ABS(values[ll+j]);
    }
  }
  return suma;
}

void Matrix::ShiftedAdd(Matrix *other, int llrow, int llcol)
{
  int i,j,ll;
  double others[other->NCol()];

  for (i=llrow;i<llrow+other->NRow();i++)
  {
    other->GetRowValues(i-llrow,others);
    ll = i*dim_col;
    for (j=llcol;j<llcol+other->NCol();j++)
    {
      values[ll+j] += others[j-llcol];
    }
  }
}

void Matrix::Transform(Matrix *U)
{
  Matrix tmp(dim_row,dim_col);
  Matrix tmp2(dim_row,dim_col);

  Multiply(U,&tmp,true);

  tmp2.Copy(U);
  tmp2.Transpose();
  tmp.Multiply(&tmp2,this,false);
}

void Matrix::SwitchRows(int a,int b)
{
  int j,lla,llb;
  double hold;

  lla = a*dim_col;
  llb = b*dim_col;

  for (j=0;j<dim_col;j++)
  {
    hold = values[lla+j];
    values[lla+j] = values[llb+j];
    values[llb+j] = hold;
  }
}

void Matrix::SwitchColumns(int a,int b)
{
  int i,ll;
  double hold;

  for (i=0;i<dim_row;i++)
  {
    ll = i*dim_col;
    hold = values[ll+a];
    values[ll+a] = values[ll+b];
    values[ll+b] = hold;
  }
}
