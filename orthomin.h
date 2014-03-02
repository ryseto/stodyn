/* header file of orthomin.c --
 * orthomin scheme
 * Copyright (C) 1999-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: orthomin.h,v 2.8 2007/11/25 18:51:03 kichiki Exp $
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
#ifndef	_ORTHOMIN_H_
#define	_ORTHOMIN_H_


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
      struct iter *it);

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
	 struct iter *it);


#endif /* !_ORTHOMIN_H_ */
