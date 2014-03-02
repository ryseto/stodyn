/* Excerpted from BLAS package
 * Copyright (C) 1999-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: myblas.c,v 1.5 2006/10/09 21:58:18 ichiki Exp $
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
#include <cmath>

#include "myblas.h"


/*    constant times a vector plus a vector.
 *    uses unrolled loops for increments equal to one.
 *    jack dongarra, linpack, 3/11/78.
 *    modified 12/3/93, array(1) declarations changed to array(*)
 * calc dy [i] = dy [i] + da * dx [i]
 * INPUT
 *    n : dimension
 *    da : constant
 *    dx [n * incx] : vector
 *    dy [n * incy] : vector
 * OUTPUT
 *    dy [n * incy] : vector
 */
void
my_daxpy (int n, double da, const double *dx, int incx, double *dy, int incy)
{
  int i, ix, iy;
  int m;

  if (n <= 0) return;
  if (da == 0.0) return;

  if (incx != 1
      || incy != 1)
    {
      /* code for unequal increments or equal increments
	 not equal to 1 */
      ix = 1;
      iy = 1;
      if (incx < 0) ix = (- n + 1) * incx + 1;
      if (incy < 0) iy = (- n + 1) * incy + 1;
      for (i = 0; i<n; i++)
	{
	  dy [iy] += da * dx [ix];
	  ix = ix + incx;
	  iy = iy + incy;
	}
    }
  else
    {
      /* code for both increments equal to 1 */

      /* clean-up loop */
      m = n % 4;
      if (m != 0)
	{
	  for (i = 0; i < m; i++)
	    dy [i] += da * dx [i];
	  if( n < 4 ) return;
	}

      for (i = m; i < n; i += 4)
	{
	  dy [i] += da * dx [i];
	  dy [i + 1] += da * dx [i + 1];
	  dy [i + 2] += da * dx [i + 2];
	  dy [i + 3] += da * dx [i + 3];
	}
    }
}

/* calc dz [i] = dy [i] + da * dx [i]
 * INPUT
 *    n : dimension
 *    da : constant
 *    dx [n * incx] : vector
 *    dy [n * incy] : vector
 * OUTPUT
 *    dz [n * incy] : vector
 */
void
my_daxpyz (int n, double da, const double *dx, int incx,
	   const double *dy, int incy,
	   double *dz, int incz)
{
  int i, ix, iy;
  int m;

  if (n <= 0) return;

  if (incx != 1
      || incy != 1)
    {
      /* code for unequal increments or equal increments
	 not equal to 1 */
      ix = 1;
      iy = 1;
      if (incx < 0) ix = (- n + 1) * incx + 1;
      if (incy < 0) iy = (- n + 1) * incy + 1;
      for (i = 0; i<n; i++)
	{
	  dz [iy] = dy [iy] + da * dx [ix];
	  ix = ix + incx;
	  iy = iy + incy;
	}
    }
  else
    {
      /* code for both increments equal to 1 */

      /* clean-up loop */
      m = n % 4;
      if (m != 0)
	{
	  for (i = 0; i < m; i++)
	    dz [i] = dy [i] + da * dx [i];
	  if( n < 4 ) return;
	}

      for (i = m; i < n; i += 4)
	{
	  dz [i] = dy [i] + da * dx [i];
	  dz [i + 1] = dy [i + 1] + da * dx [i + 1];
	  dz [i + 2] = dy [i + 2] + da * dx [i + 2];
	  dz [i + 3] = dy [i + 3] + da * dx [i + 3];
	}
    }
}


/*     copies a vector, x, to a vector, y.
 *     uses unrolled loops for increments equal to one.
 *     jack dongarra, linpack, 3/11/78.
 *     modified 12/3/93, array(1) declarations changed to array(*)
 * calc dy [i] = dx [i]
 * INPUT
 *    n : dimension
 *    dx [n * incx] : vector
 * OUTPUT
 *    dy [n * incy] : vector
 */
void
my_dcopy (int n, const double *dx, int incx,
	  double * dy, int incy)
{
  int i, ix, iy;
  int m;

  if (n <= 0) return;
  if (incx != 1
      || incy != 1)
    {
      /* code for unequal increments or equal increments
	 not equal to 1 */

      ix = 1;
      iy = 1;
      if (incx < 0) ix = (- n + 1) * incx + 1;
      if (incy < 0) iy = (- n + 1) * incy + 1;
      for (i = 0; i < n; i++)
	{
	  dy [iy] = dx [ix];
	  ix = ix + incx;
	  iy = iy + incy;
	}
    }
  else
    {
      /* code for both increments equal to 1 */

      /* clean-up loop */
      m = n % 7;
      if (m != 0 )
	{
	  for (i = 0; i < m; i ++)
	    dy [i] = dx [i];

	  if ( n  <  7 ) return;
	}

      for (i = m; i < n; i += 7)
	{
	  dy [i] = dx [i];
	  dy [i + 1] = dx [i + 1];
	  dy [i + 2] = dx [i + 2];
	  dy [i + 3] = dx [i + 3];
	  dy [i + 4] = dx [i + 4];
	  dy [i + 5] = dx [i + 5];
	  dy [i + 6] = dx [i + 6];
	}
    }
}

/*  DNRM2 returns the euclidean norm of a vector via the function
 *  name, so that
 *
 *     DNRM2 := sqrt( x'*x )
 *
 *  -- This version written on 25-October-1982.
 *     Modified on 14-October-1993 to inline the call to DLASSQ.
 *     Sven Hammarling, Nag Ltd.
 */
double
my_dnrm2_ (int n, const double *x, int incx)
{
  /* .. Parameters .. */
  static double one = 1.0;
  static double zero = 0.0;

  /* .. Local Scalars ..*/
  int ix;
  double absxi, norm, scale, ssq;


  if (n < 1 || incx < 1)
    norm = zero;
  else if (n == 1)
    norm = fabs (x [0]);
  else
    {
      scale = zero;
      ssq = one;
      /* The following loop is equivalent to this call to the LAPACK
       * auxiliary routine:
       * CALL DLASSQ( N, X, INCX, SCALE, SSQ )
       */
      for (ix = 0; ix < 1 + (n - 1) * incx; ix += incx)
	{
	  if (x [ix] != zero)
	    {
	      absxi = fabs (x [ix]);
	      if (scale < absxi)
		{
		  ssq = one + ssq * pow ((scale / absxi), 2.0);
		  scale = absxi;
		}
	      else
		ssq += pow ((absxi / scale), 2.0);
	    }
	}
      norm = scale * sqrt (ssq);
    }

  return (norm);
}

double
my_dnrm2 (int n, const double *x, int incx)
{
  int i;
  int ix;
  int m;
  double norm;


  if (n < 1 || incx < 1)
    return (0.0);
  else if (n == 1)
    return (fabs (x [0]));
  else if (incx != 1)
    {
      norm = 0.0;
      for (ix = 0; ix < 1 + (n - 1) * incx; ix += incx)
	if (x [ix] != 0.0)
	  norm += x [ix] * x [ix];
    }
  else /* incx == 1 */
    {
      norm = 0.0;
      m = n % 5;
      if (m != 0)
	{
	  for (i = 0; i < m; i ++)
	    norm += x [i] * x [i];
	}
      for (i = m; i < n; i += 5)
	{
	  norm +=
	    x [i] * x [i]
	    + x [i + 1] * x [i + 1]
	    + x [i + 2] * x [i + 2]
	    + x [i + 3] * x [i + 3]
	    + x [i + 4] * x [i + 4];
	}
    }
  return (sqrt (norm));
}


/*     forms the dot product of two vectors.
 *     uses unrolled loops for increments equal to one.
 *     jack dongarra, linpack, 3/11/78.
 *     modified 12/3/93, array(1) declarations changed to array(*)
 * calc (dy, dx)
 * INPUT
 *    n : dimension
 *    dx [n * incx] : vector
 *    dy [n * incy] : vector
 * OUTPUT (return value)
 */
double
my_ddot (int n, const double *dx, int incx, const double *dy, int incy)
{
  double dtemp;
  int i, ix, iy;
  int m;


  dtemp = 0.0;

  if (n <= 0) return (dtemp);
  if (incx != 1
      || incy != 1)
    {
      /* code for unequal increments or equal increments
	 not equal to 1 */
      ix = 1;
      iy = 1;
      if (incx < 0) ix = (- n + 1) * incx + 1;
      if (incy < 0) iy = (- n + 1) * incy + 1;
      for (i = 0; i < n; i ++)
	{
	  dtemp += dx [ix] * dy [iy];
	  ix = ix + incx;
	  iy = iy + incy;
	}
      return (dtemp);
    }
  else
    {
      /*code for both increments equal to 1*/

      /* clean-up loop */
      m = n % 5;
      if (m != 0)
	{
	  for (i = 0; i < m; i ++)
	    dtemp += dx [i] * dy [i];
	  if (n  <  5)
	    return (dtemp);
	}
      for (i = m; i < n; i += 5)
	{
	  dtemp +=
	    dx [i] * dy [i]
	    + dx [i + 1] * dy [i + 1]
	    + dx [i + 2] * dy [i + 2]
	    + dx [i + 3] * dy [i + 3]
	    + dx [i + 4] * dy [i + 4];
	}
      return (dtemp);
    }
}

/*     scales a vector by a constant.
 *     uses unrolled loops for increment equal to one.
 *     jack dongarra, linpack, 3/11/78.
 *     modified 3/93 to return if incx .le. 0.
 *     modified 12/3/93, array(1) declarations changed to array(*)
 * calc dx [i] = da * dx [i]
 */
void
my_dscal (int n, double da, double *dx, int incx)
{
  int i, m, nincx;


  if (n <= 0
      || incx <= 0)
    return;
  if (incx != 1)
    {
      /* code for increment not equal to 1 */

      nincx = n * incx;
      for (i = 0; i < nincx; i += incx)
        dx [i] *= da;
      return;
    }
  else
    {
      /* code for increment equal to 1 */

      /* clean-up loop */
      m = n % 5;
      if (m != 0)
	{
	  for (i = 0; i < m; i ++)
	    dx [i] *= da;
	  if (n < 5 ) return;
	}
      for (i = m; i < n; i += 5)
	{
	  dx [i] *= da;
	  dx [i + 1] *= da;
	  dx [i + 2] *= da;
	  dx [i + 3] *= da;
	  dx [i + 4] *= da;
	}
    }
}
/*     scales a vector by a constant.
 *     uses unrolled loops for increment equal to one.
 *     jack dongarra, linpack, 3/11/78.
 *     modified 3/93 to return if incx .le. 0.
 *     modified 12/3/93, array(1) declarations changed to array(*)
 * calc dz [i] = da * dx [i]
 */
void
my_dscalz (int n, double da, const double *dx, int incx,
	double *dz, int incz)
{
  int i, m, nincx;


  if (n <= 0
      || incx <= 0)
    return;
  if (incx != 1)
    {
      /* code for increment not equal to 1 */

      nincx = n * incx;
      for (i = 0; i < nincx; i += incx)
        dz [i] = da * dx [i];
      return;
    }
  else
    {
      /* code for increment equal to 1 */

      /* clean-up loop */
      m = n % 5;
      if (m != 0)
	{
	  for (i = 0; i < m; i ++)
	    dz [i] = da * dx [i];
	  if (n < 5 ) return;
	}
      for (i = m; i < n; i += 5)
	{
	  dz [i] = da * dx [i];
	  dz [i + 1] = da * dx [i + 1];
	  dz [i + 2] = da * dx [i + 2];
	  dz [i + 3] = da * dx [i + 3];
	  dz [i + 4] = da * dx [i + 4];
	}
    }
}
