/* orthomin scheme
 * Copyright (C) 1999-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: orthomin.c,v 2.9 2007/11/25 18:50:40 kichiki Exp $
 *
 * solver routines are translated into C by K.I. from fortran code
 * originally written by martin h. gutknecht
 *             ===============rog.f==============
 *             problem given by martin h. gutknecht
 *             numerical method: orthomin(k) method
 *             ver. 1.0 jul. 28 1995 by s. l. zhang
 *             ver. 1.1 aug. 31 1995 by s. l. zhang
 *    at slzhang.fort.iter.complex.orthomin.gutknecht-problem(gut.f)
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

#include <cstdio> /* fprintf() */
#include <cmath> /* log10() */
#include <cstdlib> /* malloc(), free() */
#include "memory-check.h" // CHECK_MALLOC
#include "libiter.h" // struct iter

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


#include "orthomin.h"


/* orthomin(k) method with BLAS/ATLAS
 * INPUT
 *   m : dimension of the problem
 *   b[m] : r-h-s vector
 *   atimes (int m, static double *x, double *b, void *param) :
 *        calc matrix-vector product A.x = b.
 *   atimes_param : parameters for atimes().
 *   it : struct iter. following entries are used
 *        it->max = kend : max of iteration
 *        it->eps = eps  : criteria for |r^2|/|b^2|
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x[m] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
int
otmk (int m, const double *b, double *x,
      void (*atimes) (int, const double *, double *, void *),
      void *atimes_param,
      struct iter *it)
{
#ifndef HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */

  int i_1 = 1;
  double d_m1 = -1.0;

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  int ret = -1;
  int kend = it->max;
  int kres = it->restart;
  double eps2 = it->eps * it->eps;

  /**
   * allocation of matrices
   * r   [m]
   * p   [(kres+1) * m]
   * ap  [(kres+1) * m]
   * beta[kres+1]
   * pap [kres+1]
   */
  double *r    = (double *)malloc (sizeof (double) * m);
  double *p    = (double *)malloc (sizeof (double) * m * (kres + 1));
  double *ap   = (double *)malloc (sizeof (double) * m * (kres + 1));
  double *beta = (double *)malloc (sizeof (double) * (kres + 1));
  double *pap  = (double *)malloc (sizeof (double) * (kres + 1));
  CHECK_MALLOC (r   , "otmk");
  CHECK_MALLOC (p   , "otmk");
  CHECK_MALLOC (ap  , "otmk");
  CHECK_MALLOC (beta, "otmk");
  CHECK_MALLOC (pap , "otmk");

  double res2 = 0.0; // for compiler warning.

#ifdef HAVE_CBLAS_H
  /**
   * ATLAS version
   */

  double b2 = cblas_ddot (m, b, 1, b, 1); // (b,b)
  eps2 *= b2;

  // r(0) = b - A.x
  atimes (m, x, r, atimes_param);
  cblas_daxpy (m, -1.0, b, 1, r, 1);
  cblas_dscal (m, -1.0, r, 1);

  // p(0) = r(0)
  cblas_dcopy (m, r, 1, p, 1);
  // ap(0) = A.p(0) (= A.r(0))
  atimes (m, p, ap, atimes_param);

  int iter;
  for (iter = 0; iter <= kend; iter ++)
    {
      res2 = cblas_ddot (m, r, 1, r, 1); // (r, r)
      if (it->debug == 2)
	{
	  fprintf (it->out, "otmk %d %e\n", iter, res2 / b2);
	}
      if(res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      int k1 = iter % (kres + 1);
      double rap = cblas_ddot (m, r, 1, ap + k1 * m, 1); // (r, ap(k1))
      // pap(k1) = (ap(k1), ap(k1))
      pap[k1] = cblas_ddot (m, ap + k1 * m, 1, ap + k1 * m, 1);

      double alpha = rap / pap[k1];
      cblas_daxpy (m, +alpha, p  + k1 * m, 1, x, 1); // x = x + alpha * p(k1)
      cblas_daxpy (m, -alpha, ap + k1 * m, 1, r, 1); // r = r - alpha * ap(k1)

      int k2 = (iter + 1) % (kres + 1);
      cblas_dcopy (m, r, 1, p + k2 * m, 1);     // p(k2) = r
      atimes (m, r, ap + k2 * m, atimes_param); // ap(k2) = A.r

      /*for (j = 0; j <= min0(kres-1,k); j ++)*/
      int jj = iter;
      if (jj > (kres - 1))
	{
	  jj = kres - 1;
	}
      int j;
      for (j = 0; j <= jj; j ++)
	{
	  int k3 = (iter - j) % (kres + 1);
	  // beta(k3) = - (ap(k2), ap(k3)) / (ap(k3), ap(k3))
	  // note that ap(k2) = A.r
	  beta[k3] = cblas_ddot (m, ap + k2 * m, 1, ap + k3 * m, 1);
	  beta[k3] = - beta[k3] / pap[k3];

	  // p(k2) += beta(k3) p(k3)
	  cblas_daxpy (m, beta [k3], p + k3 * m, 1, p + k2 * m, 1);
	  // ap(k2) += beta(k3) ap(k3)
	  cblas_daxpy (m, beta [k3], ap + k3 * m, 1, ap + k2 * m, 1);
	}
    }

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /**
   * BLAS version
   */

  double b2 = ddot_ (&m, b, &i_1, b, &i_1); // (b,b)
  eps2 *= b2;

  // r(0) = b - A.x
  atimes (m, x, r, atimes_param);
  daxpy_ (&m, &d_m1, b, &i_1, r, &i_1);
  dscal_ (&m, &d_m1, r, &i_1);

  // p(0) = r(0)
  dcopy_ (&m, r, &i_1, p, &i_1);
  // ap(0) = A.p(0) (= A.r(0))
  atimes (m, p, ap, atimes_param);

  int iter;
  for (iter = 0; iter <= kend; iter ++)
    {
      res2 = ddot_ (&m, r, &i_1, r, &i_1); // (r, r)
      if (it->debug == 2)
	{
	  fprintf (it->out, "otmk %d %e\n", iter, res2 / b2);
	}
      if(res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      int k1 = iter % (kres + 1);
      double rap = ddot_ (&m, r, &i_1, ap + k1 * m, &i_1); // (r, ap(k1))
      // pap(k1) = (ap(k1), ap(k1))
      pap[k1] = ddot_ (&m, ap + k1 * m, &i_1, ap + k1 * m, &i_1);

      double alpha = rap / pap[k1];
      double malpha = -alpha;
      daxpy_ (&m, &alpha,  p  + k1 * m, &i_1, x, &i_1); // x = x + alpha*p(k1)
      daxpy_ (&m, &malpha, ap + k1 * m, &i_1, r, &i_1); // r = r - alpha*ap(k1)

      int k2 = (iter + 1) % (kres + 1);
      dcopy_ (&m, r, &i_1, p + k2 * m, &i_1);   // p(k2) = r
      atimes (m, r, ap + k2 * m, atimes_param); // ap(k2) = A.r

      /*for (j = 0; j <= min0(kres-1,k); j ++)*/
      int jj = iter;
      if (jj > (kres - 1))
	{
	  jj = kres - 1;
	}
      int j;
      for (j = 0; j <= jj; j ++)
	{
	  int k3 = (iter - j) % (kres + 1);
	  // beta(k3) = - (ap(k2), ap(k3)) / (ap(k3), ap(k3))
	  // note that ap(k2) = A.r
	  beta[k3] = ddot_ (&m, ap + k2 * m, &i_1, ap + k3 * m, &i_1);
	  beta[k3] = - beta[k3] / pap[k3];

	  // p(k2) += beta(k3) p(k3)
	  daxpy_ (&m, beta + k3, p + k3 * m,  &i_1, p + k2 * m,  &i_1);
	  // ap(k2) += beta(k3) ap(k3)
	  daxpy_ (&m, beta + k3, ap + k3 * m, &i_1, ap + k2 * m, &i_1);
	}
    }

# else // !HAVE_BLAS_H
  /**
   * local BLAS version
   */

  double b2 = my_ddot (m, b, 1, b, 1); // (b,b)
  eps2 *= b2;

  // r(0) = b - A.x
  atimes (m, x, r, atimes_param);
  my_daxpy (m, -1.0, b, 1, r, 1);
  my_dscal (m, -1.0, r, 1);

  // p(0) = r(0)
  my_dcopy (m, r, 1, p, 1);
  // ap(0) = A.p(0) (= A.r(0))
  atimes (m, p, ap, atimes_param);

  int iter;
  for (iter = 0; iter <= kend; iter ++)
    {
      res2 = my_ddot (m, r, 1, r, 1); // (r, r)
      if (it->debug == 2)
	{
	  fprintf (it->out, "otmk %d %e\n", iter, res2 / b2);
	}
      if(res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      int k1 = iter % (kres + 1);
      double rap = my_ddot (m, r, 1, ap + k1 * m, 1); // (r, ap(k1))
      // pap(k1) = (ap(k1), ap(k1))
      pap[k1] = my_ddot (m, ap + k1 * m, 1, ap + k1 * m, 1);

      double alpha = rap / pap[k1];
      my_daxpy (m, +alpha, p  + k1 * m, 1, x, 1); // x = x + alpha * p(k1)
      my_daxpy (m, -alpha, ap + k1 * m, 1, r, 1); // r = r - alpha * ap(k1)

      int k2 = (iter + 1) % (kres + 1);
      my_dcopy (m, r, 1, p + k2 * m, 1);     // p(k2) = r
      atimes (m, r, ap + k2 * m, atimes_param); // ap(k2) = A.r

      /*for (j = 0; j <= min0(kres-1,k); j ++)*/
      int jj = iter;
      if (jj > (kres - 1))
	{
	  jj = kres - 1;
	}
      int j;
      for (j = 0; j <= jj; j ++)
	{
	  int k3 = (iter - j) % (kres + 1);
	  // beta(k3) = - (ap(k2), ap(k3)) / (ap(k3), ap(k3))
	  // note that ap(k2) = A.r
	  beta[k3] = my_ddot (m, ap + k2 * m, 1, ap + k3 * m, 1);
	  beta[k3] = - beta[k3] / pap[k3];

	  // p(k2) += beta(k3) p(k3)
	  my_daxpy (m, beta [k3], p + k3 * m, 1, p + k2 * m, 1);
	  // ap(k2) += beta(k3) ap(k3)
	  my_daxpy (m, beta [k3], ap + k3 * m, 1, ap + k2 * m, 1);
	}
    }

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  free (r);
  free (p);
  free (ap);
  free (beta);
  free (pap);

  if (it->debug == 1)
    {
      fprintf (it->out, "otmk %d %e\n", iter, res2 / b2);
    }

  it->niter = iter;
  it->res2  = res2 / b2;
  return (ret);
}

/* orthomin(k) method with preconditioning
 * INPUT
 *   m : dimension of the problem
 *   b[m] : r-h-s vector
 *   atimes (int m, static double *x, double *b, void *param) :
 *        calc matrix-vector product A.x = b.
 *   atimes_param : parameters for atimes().
 *   inv (int m, static double *b, double *x, void *param) :
 *        approx of A^{-1}.b = x for preconditioning.
 *   inv_param : parameters for the preconditioner inv().
 *   it : struct iter. following entries are used
 *        it->max = kend : max of iteration
 *        it->eps = eps  : criteria for |r^2|/|b^2|
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x[m] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
int
otmk_pc (int m, const double *b, double *x,
	 void (*atimes) (int, const double *, double *, void *),
	 void *atimes_param,
	 void (*inv) (int, const double *, double *, void *),
	 void *inv_param,
	 struct iter *it)
{
#ifndef HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */

  int i_1 = 1;
  double d_m1 = -1.0;

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  int ret = -1;
  int kend = it->max;
  int kres = it->restart;
  double eps2 = it->eps * it->eps;

  /**
   * allocation of matrices
   * r   [m]
   * p   [(kres+1) * m]
   * ap  [(kres+1) * m]
   * beta[kres+1]
   * pap [kres+1]
   * Kr  [m] for K^{-1}.r, where K is the preconditioner.
   */
  double *r    = (double *)malloc (sizeof (double) * m);
  double *p    = (double *)malloc (sizeof (double) * m * (kres + 1));
  double *ap   = (double *)malloc (sizeof (double) * m * (kres + 1));
  double *beta = (double *)malloc (sizeof (double) * (kres + 1));
  double *pap  = (double *)malloc (sizeof (double) * (kres + 1));
  CHECK_MALLOC (r   , "otmk_pc");
  CHECK_MALLOC (p   , "otmk_pc");
  CHECK_MALLOC (ap  , "otmk_pc");
  CHECK_MALLOC (beta, "otmk_pc");
  CHECK_MALLOC (pap , "otmk_pc");

  double res2 = 0.0; // for compiler warning.

#ifdef HAVE_CBLAS_H
  /**
   * ATLAS version
   */

  double b2 = cblas_ddot (m, b, 1, b, 1); // (b,b)
  eps2 *= b2;

  // r(0) = b - A.x
  atimes (m, x, r, atimes_param);
  cblas_daxpy (m, -1.0, b, 1, r, 1);
  cblas_dscal (m, -1.0, r, 1);

  // p(0) = K^{-1}.r(0)
  inv (m, r, p, inv_param);
  // ap(0) = A.p(0) (= A.r(0))
  atimes (m, p, ap, atimes_param);

  int iter;
  for (iter = 0; iter <= kend; iter ++)
    {
      res2 = cblas_ddot (m, r, 1, r, 1); // (r, r)
      if (it->debug == 2)
	{
	  fprintf (it->out, "otmk_pc %d %e\n", iter, res2 /b2);
	}
      if(res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      int k1 = iter % (kres + 1);
      double rap = cblas_ddot (m, r, 1, ap + k1 * m, 1); // (r, ap(k1))
      // pap(k1) = (ap(k1), ap(k1))
      pap[k1]    = cblas_ddot (m, ap + k1 * m, 1, ap + k1 * m, 1);

      // alpha = (r, ap(k1)) / (ap(k1), ap(k1))
      double alpha = rap / pap[k1];
      cblas_daxpy (m, +alpha, p  + k1 * m, 1, x, 1); // x = x + alpha * p(k1)
      cblas_daxpy (m, -alpha, ap + k1 * m, 1, r, 1); // r = r - alpha * ap(k1)
      
      int k2 = (iter + 1) % (kres + 1);
      inv (m, r, p + k2 * m, inv_param);                 // p(k2)  = K^{-1}.r
      atimes (m, p + k2 * m, ap + k2 * m, atimes_param); // ap(k2) = A.p(k2)

      /*for (j = 0; j <= min0(kres-1,k); j ++)*/
      int jj = iter;
      if (jj > (kres - 1))
	{
	  jj = kres - 1;
	}
      int j;
      for (j = 0; j <= jj; j ++)
	{
	  int k3 = (iter - j) % (kres + 1);
	  // beta(k3) = - (ap(k2), ap(k3)) / (ap(k3), ap(k3))
	  // note that ap(k2) = A.K^{-1}.r
	  beta[k3] = cblas_ddot (m, ap + k2 * m, 1, ap + k3 * m, 1);
	  beta[k3] = - beta[k3] / pap[k3];

	  // p(k2) += beta(k3) p(k3)
	  cblas_daxpy (m, beta [k3], p + k3 * m, 1, p + k2 * m, 1);
	  // ap(k2) += beta(k3) ap(k3)
	  cblas_daxpy (m, beta [k3], ap + k3 * m, 1, ap + k2 * m, 1);
	}
    }

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /**
   * BLAS version
   */

  double b2 = ddot_ (&m, b, &i_1, b, &i_1); // (b,b)
  eps2 *= b2;

  // r(0) = b - A.x
  atimes (m, x, r, atimes_param);
  daxpy_ (&m, &d_m1, b, &i_1, r, &i_1);
  dscal_ (&m, &d_m1, r, &i_1);

  // p(0) = K^{-1}.r(0)
  inv (m, r, p, inv_param);
  // ap(0) = A.p(0) (= A.r(0))
  atimes (m, p, ap, atimes_param);

  int iter;
  for (iter = 0; iter <= kend; iter ++)
    {
      res2 = ddot_ (&m, r, &i_1, r, &i_1); // (r, r)
      if (it->debug == 2)
	{
	  fprintf (it->out, "otmk_pc %d %e\n", iter, res2 / b2);
	}
      if(res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      int k1 = iter % (kres + 1);
      double rap = ddot_ (&m, r, &i_1, ap + k1 * m, &i_1); // (r, ap(k1))
      // pap(k1) = (ap(k1), ap(k1))
      pap[k1] = ddot_ (&m, ap + k1 * m, &i_1, ap + k1 * m, &i_1);

      // alpha = (r, ap(k1)) / (ap(k1), ap(k1))
      double alpha = rap / pap[k1];
      double malpha = -alpha;
      daxpy_ (&m, &alpha,  p  + k1 * m, &i_1, x, &i_1); // x = x + alpha*p(k1)
      daxpy_ (&m, &malpha, ap + k1 * m, &i_1, r, &i_1); // r = r - alpha*ap(k1)
      
      int k2 = (iter + 1) % (kres + 1);
      inv (m, r, p + k2 * m, inv_param);                 // p(k2)  = K^{-1}.r
      atimes (m, p + k2 * m, ap + k2 * m, atimes_param); // ap(k2) = A.p(k2)

      /*for (j = 0; j <= min0(kres-1,k); j ++)*/
      int jj = iter;
      if (jj > (kres - 1))
	{
	  jj = kres - 1;
	}
      int j;
      for (j = 0; j <= jj; j ++)
	{
	  int k3 = (iter - j) % (kres + 1);
	  // beta(k3) = - (ap(k2), ap(k3)) / (ap(k3), ap(k3))
	  // note that ap(k2) = A.K^{-1}.r
	  beta[k3] = ddot_ (&m, ap + k2 * m, &i_1, ap + k3 * m, &i_1);
	  beta[k3] = - beta[k3] / pap[k3];

	  // p(k2) += beta(k3) p(k3)
	  daxpy_ (&m, beta + k3, p + k3 * m,  &i_1, p + k2 * m,  &i_1);
	  // ap(k2) += beta(k3) ap(k3)
	  daxpy_ (&m, beta + k3, ap + k3 * m, &i_1, ap + k2 * m, &i_1);
	}
    }

# else // !HAVE_BLAS_H
  /**
   * local BLAS version
   */

  double b2 = my_ddot (m, b, 1, b, 1); // (b,b)
  eps2 *= b2;

  // r(0) = b - A.x
  atimes (m, x, r, atimes_param);
  my_daxpy (m, -1.0, b, 1, r, 1);
  my_dscal (m, -1.0, r, 1);

  // p(0) = K^{-1}.r(0)
  inv (m, r, p, inv_param);
  // ap(0) = A.p(0) (= A.r(0))
  atimes (m, p, ap, atimes_param);

  int iter;
  for (iter = 0; iter <= kend; iter ++)
    {
      res2 = my_ddot (m, r, 1, r, 1); // (r, r)
      if (it->debug == 2)
	{
	  fprintf (it->out, "otmk_pc %d %e\n", iter, res2 / b2);
	}
      if(res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      int k1 = iter % (kres + 1);
      double rap = my_ddot (m, r, 1, ap + k1 * m, 1); // (r, ap(k1))
      // pap(k1) = (ap(k1), ap(k1))
      pap[k1]    = my_ddot (m, ap + k1 * m, 1, ap + k1 * m, 1);

      // alpha = (r, ap(k1)) / (ap(k1), ap(k1))
      double alpha = rap / pap[k1];
      my_daxpy (m, +alpha, p  + k1 * m, 1, x, 1); // x = x + alpha * p(k1)
      my_daxpy (m, -alpha, ap + k1 * m, 1, r, 1); // r = r - alpha * ap(k1)
      
      int k2 = (iter + 1) % (kres + 1);
      inv (m, r, p + k2 * m, inv_param);                 // p(k2)  = K^{-1}.r
      atimes (m, p + k2 * m, ap + k2 * m, atimes_param); // ap(k2) = A.p(k2)

      /*for (j = 0; j <= min0(kres-1,k); j ++)*/
      int jj = iter;
      if (jj > (kres - 1))
	{
	  jj = kres - 1;
	}
      int j;
      for (j = 0; j <= jj; j ++)
	{
	  int k3 = (iter - j) % (kres + 1);
	  // beta(k3) = - (ap(k2), ap(k3)) / (ap(k3), ap(k3))
	  // note that ap(k2) = A.K^{-1}.r
	  beta[k3] = my_ddot (m, ap + k2 * m, 1, ap + k3 * m, 1);
	  beta[k3] = - beta[k3] / pap[k3];

	  // p(k2) += beta(k3) p(k3)
	  my_daxpy (m, beta [k3], p + k3 * m, 1, p + k2 * m, 1);
	  // ap(k2) += beta(k3) ap(k3)
	  my_daxpy (m, beta [k3], ap + k3 * m, 1, ap + k2 * m, 1);
	}
    }

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  free (r);
  free (p);
  free (ap);
  free (beta);
  free (pap);

  if (it->debug == 1)
    {
      fprintf (it->out, "otmk_pc %d %e\n", iter, res2);
    }

  it->niter = iter;
  it->res2  = res2 / b2;
  return (ret);
}
