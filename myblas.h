/* header file for blas.c --
 * Excerpted from BLAS package
 * Copyright (C) 1999-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: myblas.h,v 1.4 2006/10/09 21:58:40 ichiki Exp $
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
#ifndef	_MYBLAS_H_
#define	_MYBLAS_H_


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
my_daxpy (int n, double da, const double *dx, int incx, double *dy, int incy);

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
	   double *dz, int incz);

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
	  double * dy, int incy);

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
my_dnrm2_ (int n, const double *x, int incx);

double
my_dnrm2 (int n, const double *x, int incx);

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
my_ddot (int n, const double *dx, int incx, const double *dy, int incy);

/*     scales a vector by a constant.
 *     uses unrolled loops for increment equal to one.
 *     jack dongarra, linpack, 3/11/78.
 *     modified 3/93 to return if incx .le. 0.
 *     modified 12/3/93, array(1) declarations changed to array(*)
 * calc dx [i] = da * dx [i]
 */
void
my_dscal (int n, double da, double *dx, int incx);
/*     scales a vector by a constant.
 *     uses unrolled loops for increment equal to one.
 *     jack dongarra, linpack, 3/11/78.
 *     modified 3/93 to return if incx .le. 0.
 *     modified 12/3/93, array(1) declarations changed to array(*)
 * calc dz [i] = da * dx [i]
 */
void
my_dscalz (int n, double da, const double *dx, int incx,
	   double *dz, int incz);


#endif /* !_MYBLAS_H_ */
