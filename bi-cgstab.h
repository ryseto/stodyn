/* header file of bi-cgstab.c --
 * wrapper for iterative solver routines
 * Copyright (C) 1999-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bi-cgstab.h,v 2.8 2007/11/25 18:42:45 kichiki Exp $
 *
 * solver routines are translated into C by K.I. from fortran code
 * originally written by martin h. gutknecht
 *             ===============rog.f==============
 *             problem given by martin h. gutknecht
 *             numerical method: bi-cgsta        method -- sta ()
 *             numerical method: bi-cgsta2       method -- st2 ()
 *             numerical method: gpbi-cg         method -- gpb ()
 *             ver. 1.0 aug. 08 1995  s. l. zhang
 *    at urus.slzhang.fort.iter.real.gpbcg.gutknecht-p(final.f)
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
#ifndef	_BI_CGSTAB_H_
#define	_BI_CGSTAB_H_


/* bi-cgstab method
 * INPUT
 *   m : dimension of the problem
 *   b [m] : r-h-s vector
 *   atimes (int m, static double *x, double *b, void *param) :
 *        calc matrix-vector product A.x = b.
 *   atimes_param : parameters for atimes().
 *   it : struct iter. following entries are used
 *        it->max = kend : max of iteration
 *        it->eps = eps  : criteria for |r^2|/|b^2|
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x [m] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
int
sta (int m, const double *b, double *x,
     void (*atimes) (int, const double *, double *, void *),
     void *atimes_param,
     struct iter *it);

/* bi-cgstab method with precondition
 * INPUT
 *   m : dimension of the problem
 *   b [m] : r-h-s vector
 *   atimes (int m, static double *x, double *b, void *param) :
 *        calc matrix-vector product A.x = b.
 *   atimes_param : parameters for atimes().
 *   inv (int m, static double *b, double *x, void *param) :
 *        approx of A^{-1}.b = x for preconditioning.
 *   inv_param : parameters for the preconditioner inv().
 *   it : struct iter. following entries are used
 *        it->max = kend : max of iteration
 *        it->eps = eps : log10 of cutoff
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x [m] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
int
sta_pc (int m, const double *b, double *x,
	void (*atimes) (int, const double *, double *, void *),
	void *atimes_param,
	void (*inv) (int, const double *, double *, void *),
	void *inv_param,
	struct iter *it);

/* bi-cgstab2 method
 * INPUT
 *   m : dimension of the problem
 *   b [m] : r-h-s vector
 *   atimes (int m, static double *x, double *b, void *param) :
 *        calc matrix-vector product A.x = b.
 *   atimes_param : parameters for atimes().
 *   it : struct iter. following entries are used
 *        it->max = kend : max of iteration
 *        it->eps = eps : log10 of cutoff
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x [m] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
int
sta2 (int m, const double *b, double *x,
      void (*atimes) (int, const double *, double *, void *),
      void *atimes_param,
      struct iter *it);

/* bi-cgstab2 method with precondition
 * INPUT
 *   m : dimension of the problem
 *   b [m] : r-h-s vector
 *   atimes (int m, static double *x, double *b, void *param) :
 *        calc matrix-vector product A.x = b.
 *   atimes_param : parameters for atimes().
 *   inv (int m, static double *b, double *x, void *param) :
 *        approx of A^{-1}.b = x for preconditioning.
 *   inv_param : parameters for the preconditioner inv().
 *   it : struct iter. following entries are used
 *        it->max = kend : max of iteration
 *        it->eps = eps : log10 of cutoff
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x [m] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
int
sta2_pc (int m, const double *b, double *x,
	 void (*atimes) (int, const double *, double *, void *),
	 void *atimes_param,
	 void (*inv) (int, const double *, double *, void *),
	 void *inv_param,
	 struct iter *it);

/* gpbi-cg method
 * ref: Zhang, SIAM J.Sci.Comput. 1997 vol.18 pp.537-551.
 * INPUT
 *   m : dimension of the problem
 *   b [m] : r-h-s vector
 *   atimes (int m, double *x, double *b) : calc matrix-vector product
 *   atimes_param : pointer to be passed to atimes routines
 *   it : struct iter. following entries are used
 *        it->max = kend : max of iteration
 *        it->eps = eps : log10 of cutoff
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x [m] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
int
gpb (int m, const double *b, double *x,
     void (*atimes) (int, const double *, double *, void *),
     void *atimes_param,
     struct iter *it);

/* gpbi-cg method with precondition
 * ref: Zhang, SIAM J.Sci.Comput. 1997 vol.18 pp.537-551.
 * INPUT
 *   m : dimension of the problem
 *   b [m] : r-h-s vector
 *   atimes (int m, static double *x, double *b, void *param) :
 *        calc matrix-vector product A.x = b.
 *   atimes_param : parameters for atimes().
 *   inv (int m, static double *b, double *x, void *param) :
 *        approx of A^{-1}.b = x for preconditioning.
 *   inv_param : parameters for the preconditioner inv().
 *   it : struct iter. following entries are used
 *        it->max = kend : max of iteration
 *        it->eps = eps : log10 of cutoff
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x [m] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
int
gpb_pc (int m, const double *b, double *x,
	void (*atimes) (int, const double *, double *, void *),
	void *atimes_param,
	void (*inv) (int, const double *, double *, void *),
	void *inv_param,
	struct iter *it);


#endif /* !_BI_CGSTAB_H_ */
