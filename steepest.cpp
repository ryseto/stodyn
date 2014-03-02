/* Steepest Descent -- Weiss' Algorithm 1
 * Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: steepest.c,v 2.4 2007/11/25 18:52:48 kichiki Exp $
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
#include "config.h"

#include <cstdio>
#include <cstdlib>
#include "libiter.h"
#include "memory-check.h" // CHECK_MALLOC


#ifdef HAVE_CBLAS_H
/* use ATLAS' CBLAS routines */

#ifdef OSX
#include <Accelerate/Accelerate.h>
#else
extern "C"{
#include <atlas/cblas.h>
}
#endif
#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
/* use Fortran BLAS routines */

double
ddot_(int* N, 
      double* X, int* incX, 
      double* Y, int* incY);
double
dnrm2_(int* N, 
       double* X, int* incX);
int
dcopy_(int* N,
       double* X, int* incX,
       double* Y, int* incY);
int
daxpy_(int* N,
       double* alpha,
       double* X, int* incX,
       double* Y, int* incY);
int
dscal_(int* N,
       double* alpha,
       double* X, int* incX);

# else // !HAVE_BLAS_H
/* use local BLAS routines */

#include "myblas.h"

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H


/* Steepest Descent -- Weiss' Algorithm 1
 * INPUT
 *   n : dimension of the problem
 *   b [n] : r-h-s vector
 *   atimes (int n, static double *x, double *b, void *param) :
 *        calc matrix-vector product A.x = b.
 *   atimes_param : parameters for atimes().
 *   it : struct iter. following entries are used
 *        it->max = kend : max of iteration
 *        it->eps = eps  : criteria for |r^2|/|b^2|
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x [n] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
int
steepest (int n, const double *b, double *x,
	  void (*atimes) (int, const double *, double *, void *),
	  void *atimes_param,
	  struct iter * it_param)
{
#ifndef HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */

  int i_1 = 1;
  double d_m1 = -1.0;

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  int ret = -1;
  double eps2 = it_param->eps * it_param->eps;
  int itmax = it_param->max;

  double *r  = (double *) malloc (sizeof (double) * n);
  double *ar = (double *) malloc (sizeof (double) * n);
  CHECK_MALLOC (r,  "steepest");
  CHECK_MALLOC (ar, "steepest");

  double res2 = 0.0; // for compiler warning

#ifdef HAVE_CBLAS_H
  /**
   * ATLAS version
   */

  double b2 = cblas_ddot (n, b, 1, b, 1); // (b,b)
  eps2 *= b2;

  atimes (n, x, r, atimes_param); // r = A.x
  cblas_daxpy (n, -1.0, b, 1, r, 1); // r = A.x - b
  
  int i;
  for (i = 0; i < itmax; i ++)
    {
      res2 = cblas_ddot (n, r, 1, r, 1); // res2 = (r, r)
      if (it_param->debug == 2)
	{
	  fprintf (it_param->out, "libiter-steepest %d %e\n", i, res2 / b2);
	}
      if (res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}
      
      atimes (n, r, ar, atimes_param); // ar = A.r
      double rar = cblas_ddot (n, r, 1, ar, 1); // rar = (r, A.r)

      double delta = - res2 / rar; // delta = - (r, r) / (r, A.r)
      
      cblas_daxpy (n, delta, r,  1, x, 1); // x += delta r
      cblas_daxpy (n, delta, ar, 1, r, 1); // r += delta A.r

    }

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /**
   * BLAS version
   */

  double b2 = ddot_ (&n, b, &i_1, b, &i_1); // (b,b)
  eps2 *= b2;

  atimes (n, x, r, atimes_param); // r = A.x
  daxpy_ (&n, &d_m1, b, &i_1, r, &i_1); // r = A.x - b
  
  int i;
  for (i = 0; i < itmax; i ++)
    {
      res2 = ddot_ (&n, r, &i_1, r, &i_1); // res2 = (r, r)
      if (it_param->debug == 2)
	{
	  fprintf (it_param->out, "libiter-steepest %d %e\n", i, res2 / b2);
	}
      if (res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}
      
      atimes (n, r, ar, atimes_param); // ar = A.r
      double rar = ddot_ (&n, r, &i_1, ar, &i_1); // rar = (r, A.r)

      double delta = - res2 / rar; // delta = - (r, r) / (r, A.r)
      
      daxpy_ (&n, &delta, r,  &i_1, x, &i_1); // x += delta r
      daxpy_ (&n, &delta, ar, &i_1, r, &i_1); // r += delta A.r

    }

# else // !HAVE_BLAS_H
  /**
   * local BLAS version
   */

  double b2 = my_ddot (n, b, 1, b, 1); // (b,b)
  eps2 *= b2;

  atimes (n, x, r, atimes_param); // r = A.x
  my_daxpy (n, -1.0, b, 1, r, 1); // r = A.x - b
  
  int i;
  for (i = 0; i < itmax; i ++)
    {
      res2 = my_ddot (n, r, 1, r, 1); // res2 = (r, r)
      if (it_param->debug == 2)
	{
	  fprintf (it_param->out, "libiter-steepest %d %e\n", i, res2 / b2);
	}
      if (res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}
      
      atimes (n, r, ar, atimes_param); // ar = A.r
      double rar = my_ddot (n, r, 1, ar, 1); // rar = (r, A.r)

      double delta = - res2 / rar; // delta = - (r, r) / (r, A.r)
      
      my_daxpy (n, delta, r,  1, x, 1); // x += delta r
      my_daxpy (n, delta, ar, 1, r, 1); // r += delta A.r

    }

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H


  free (r);
  free (ar);

  if (it_param->debug == 1)
    {
      fprintf (it_param->out, "libiter-steepest %d %e\n", i, res2 / b2);
    }

  it_param->niter = i;
  it_param->res2  = res2 / b2;
  return (ret);
}
