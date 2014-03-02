/* header file of mygmres.c --
 * generalized minimum residual method
 * Copyright (C) 1998-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: gmres.h,v 2.9 2007/11/25 18:41:57 kichiki Exp $
 *
 * Reference :
 *   GMRES(m) : Y.Saad & M.H.Schultz, SIAM J.Sci.Stat.Comput.
 *   vol7 (1986) pp.856-869
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
#ifndef	_GMRES_H_
#define	_GMRES_H_


/* GMRES(m)
 * INPUT
 *   n : dimension of the problem
 *   b[n] : r-h-s vector
 *   atimes (int n, static double *x, double *b, void *param) :
 *        calc matrix-vector product A.x = b.
 *   atimes_param : parameters for atimes().
 *   it : struct iter. the following entries are used.
 *        it->max     : max of iteration
 *        it->restart : number of iteration at once
 *        it->eps     : criteria for |r^2|/|b^2|
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x[m] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
int
gmres_m (int n, const double *f, double *x,
	 void (*atimes) (int, const double *, double *, void *),
	 void *atimes_param,
	 struct iter *it);

/* GMRES(m) with preconditioning
 * INPUT
 *   n : dimension of the problem
 *   b[n] : r-h-s vector
 *   atimes (int n, static double *x, double *b, void *param) :
 *        calc matrix-vector product A.x = b.
 *   atimes_param : parameters for atimes().
 *   inv (int n, static double *b, double *x, void *param) :
 *        approx of A^{-1}.b = x for preconditioning.
 *   inv_param : parameters for the preconditioner inv().
 *   it : struct iter. the following entries are used.
 *        it->max     : max of iteration
 *        it->restart : number of iteration at once
 *        it->eps     : criteria for |r^2|/|b^2|
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x[n] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
int
gmres_m_pc (int n, const double *f, double *x,
	    void (*atimes) (int, const double *, double *, void *),
	    void *atimes_param,
	    void (*inv) (int, const double *, double *, void *),
	    void *inv_param,
	    struct iter *it);

/* plain GMRES (not restart version)
 * INPUT
 *   m : dimension of the problem
 *   b[m] : r-h-s vector
 *   atimes (int m, static double *x, double *b, void *param) :
 *        calc matrix-vector product A.x = b.
 *   atimes_param : parameters for atimes().
 *   it : struct iter. following entries are used
 *        it->max : max of iteration
 *        it->eps : criteria for |r^2|/|b^2|
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x[m] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
int
gmres (int n, const double *f, double *x,
       void (*atimes) (int, const double *, double *, void *),
       void *atimes_param,
       struct iter *it);


#endif /* !_GMRES_H_ */
