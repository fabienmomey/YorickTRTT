/*
 * gauss.c --
 *
 * Implement resolution of linear system, matrix inversion and determinant by
 * means of Gaussian elimination (LU decomposition).
 *
 * This code is largely customizable (see comments below) and is designed to
 * be embedded (#includ'ed) by another C code after proper definition of the
 * customization macros.  Note that these routines are targeted at small
 * linear systems (packages such as LAPACK would be more efficient for large
 * systems).
 *
 * This code is largely inspired by routines from "Numerical Recipes" (by
 * W.H. Press, S.A. Teukolsky, W.T. Vetterling & B.P. Flannery, Cambridge
 * University Press, 1992) with many improvements regarding the matrix storage
 * order and style, the use of zero-based indexes, workspace arrays, etc.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2011, the MiTiV Team.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can use, modify
 * and/or redistribute the software under the terms of the CeCILL-C license as
 * circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty and the software's author, the holder of the
 * economic rights, and the successive licensors have only limited liability.
 *
 * In this respect, the user's attention is drawn to the risks associated with
 * loading, using, modifying and/or developing or reproducing the software by
 * the user in light of its specific status of free software, that may mean
 * that it is complicated to manipulate, and that also therefore means that it
 * is reserved for developers and experienced professionals having in-depth
 * computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling
 * the security of their systems and/or data to be ensured and, more
 * generally, to use and operate it in the same conditions as regards
 * security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 *
 *-----------------------------------------------------------------------------
 *
 * $Id$
 * $Log$
 */

/*
 * Title: LU decomposition tools.
 *
 *   Version 1.0
 *
 *   Copyright (C) 2011, the MiTiV project.
 */

/*---------------------------------------------------------------------------*/
/* CUSTOMIZATION */

/*
 * This code is customizable to accommodate to most matrix storage order and
 * style.  Column-major and row-major storage are supported.  Matrices stored
 * in uni-dimensional arrays (possibly with strides, a.k.a. "leading
 * dimension" in LAPACK jargon) or as 2-D arrays are also supported.
 *
 * Assuming that i is the row number, and j is the column number (as in
 * conventinal matrix notation except that indices are zero-based) element
 * A(i,j) of matrix A is given by the following table:
 *
 *   +----------------------+ 
 *   | STYLE     A(i,j)     | 
 *   +----------------------+ 
 *   |   0     a[j*lda + i] |     Column-major storage, LAPACK style
 *   |   1     a[i*lda + j] |     Row-major storage
 *   |   2     a[j*n + i]   |     Column-major storage
 *   |   3     a[i*n + j]   |     Row-major storage
 *   |   4     a[j][i]      |     Column-major storage
 *   |   5     a[i][j]      |     Row-major storage, Numerical Recipes
 *   +-------------------- -+ 
 *
 * where N is the size of the matrix, LDA is its leading dimension.
 */
#ifndef STYLE
# define STYLE 0
#endif

/* The following macro defines the type for index values (must be a signed
   integer, like "int" or "long", the default is "int" as these routines are
   intended for small linear systems (LAPACK would be more efficient for large
   systems). */
#ifndef INTEGER
# define INTEGER int
#endif

/* The following macro defines the scope of the functions defined in this file
   (static, extern, etc.) */
#ifndef SCOPE
# define SCOPE extern
#endif

/* The following macro defines the type for real values (float or double). */
#ifndef REAL
# define REAL double
#endif

/* Macro TINY can be defined with a very small value (strictly greater than
   zero, e.g. 1E-40) to cope with near numerically singular matrices. */
/* #define TINY 1E-40 */

/* END OF CUSTOMIZATION */
/*---------------------------------------------------------------------------*/

/* TIMINGS
 *
 * Depending on the edition of Numerical Recipes, there are two versions of
 * Crout's algorithm for LU decomposition.  The one in the 2nd edition,
 * LUDEC_COL in this code, is more efficient for column-major storage order;
 * while the one given in the 3rd edition, LUDEC_ROW in this code, is more
 * efficient for row-major storage order.
 *
 * The following table gives the timings for solving linear problems of
 * different sizes (LU decomposition followed by forward and backward
 * substitutions) with LAPACK and Numerical Recipes algorithms.  LAPACK is
 * always the faster; in parenthesis is given the ratio between the time spent
 * by a given algorithm and by LAPACK version.  The computer is a Linux laptop
 * with Intel Core i7 (Q820) CPU at 1.73GHz and 8MB cache.

 * -----------------------------------------------------------------
 * Size     Order          LAPACK     LUDEC_COL       LUDEC_ROW
 * -----------------------------------------------------------------
 *  100     row-major      390 µs     560 µs          500 µs
 *          column-major   390 µs     560 µs          510 µs
 *  200     row-major      2.9 ms     4.0 ms          3.6 ms
 *          column-major   2.9 ms     4.1 ms          4.0 ms
 *  500     row-major       40 ms      61 ms           51 ms
 *          column-major    40 ms      62 ms          130 ms
 * 1000     row-major      295 ms    1400 ms          395 ms
 *          column-major   305 ms    1400 ms         4400 ms
 * 2000     row-major     2.53 s     21.2 s (8.4)     4.3 s ( 1.7)
 *          column-major  2.55 s     22.3 s (8.7)    48.3 s (18.9)
 * -----------------------------------------------------------------
 */


/* Note that we undefine macros prior to defining them to avoid clashes in
   case of multiple inclusion. */

/* Macro ABS is used to get a the absolute value of its argument. */
#undef ABS
#define ABS(x)  ((x) >= 0 ? (x) : -(x))

/* Macro FETCH is used to get a given matrix element. */
#undef FETCH
#if STYLE == 0
# define FETCH(a,i,j)       a[(j)*ld##a + (i)]
#elif STYLE == 1
# define FETCH(a,i,j)       a[(i)*ld##a + (j)]
#elif STYLE == 2
# define FETCH(a,i,j)       a[(j)*n + (i)]
#elif STYLE == 3
# define FETCH(a,i,j)       a[(i)*n + (j)]
#elif STYLE == 4
# define FETCH(a,i,j)       a[j][i]
#elif STYLE == 5
# define FETCH(a,i,j)       a[i][j]
#else
# error invalid value for STYLE macro
#endif

/* Macro VECTOR is used to declare a vector in the list of arguments of a
   function. */
#define VECTOR(v)    v[]

/* Macro MATRIX is used to declare a matrix in the list of arguments of a
   function. */
#undef MATRIX
#if STYLE <= 1
# define MATRIX(a)       a[], const INTEGER ld##a
#elif STYLE <= 3
# define MATRIX(a)       a[]
#else
# define MATRIX(a)       a[][]
#endif

/* Macro PUT_MATRIX is used to put a matrix argument in a function call. */
#undef PUT_MATRIX
#if STYLE <= 1
# define PUT_MATRIX(a)  a, ld##a
#else
# define PUT_MATRIX(a)  a
#endif

/* Definitions of functions. */

/* Choose the most efficient algorithm for LU decomposition according
   to the storage order. */ 
#ifdef LUDEC
#  if ((STYLE & 1) == 0) /* row-major storage order */
#    ifndef LUDEC_ROW
#      define LUDEC_ROW LUDEC
#    endif
#  else /* column-major storage orde */
#    ifndef LUDEC_COL
#      define LUDEC_COL LUDEC
#    endif
#  endif
#endif

#ifdef LUDEC
/*
 * Function: LUDEC
 *   Compute the LU decomposition a square matrix.
 *
 * Parameters:
 *   a   - on input, the _n_ by _n_ matrix _a_; on output, the LU
 *         decomposition of _a_;
 *   lda - the leading dimension of _a_;
 *   n   - the dimension of the matrix;
 *   pvt - an _n_ + 1 element vector used to store permutations of rows
 *         in the LU decomposition; on output, the last element of _pvt_
 *         is: 0 if the matrix is singular, or (-1)^m with m the
 *         number of permutations;
 *   ws  - an _n_ element workspace vector used to store row scaling.
 *
 * Returns:
 *   0 on success; -1 if the input matrix is singular.
 */
SCOPE int LUDEC(REAL MATRIX(a), const INTEGER n,
                INTEGER VECTOR(pvt), REAL VECTOR(ws));
#endif /* LUDEC */
#ifdef LUDEC_ROW
SCOPE int LUDEC_ROW(REAL MATRIX(a), const INTEGER n, INTEGER VECTOR(pvt),
                    REAL VECTOR(ws));
#endif /* LUDEC_COL */
#ifdef LUDEC_COL
SCOPE int LUDEC_COL(REAL MATRIX(a), const INTEGER n, INTEGER VECTOR(pvt),
                    REAL VECTOR(ws));
#endif /* LUDEC_COL */

#ifdef LUSLV
/*
 * Function: LUSLV
 *   Solve a linear system given the LU decomposition of the left-hand-side
 *   matrix.
 *
 * Parameters:
 *   a   - the LU decomposition of the _n_ by _n_ left-hand-side matrix;
 *   lda - the leading dimension of _a_;
 *   n   - the dimension of the matrix;
 *   pvt - the _n_ + 1 element vector with the row permutations of the
 *         LU decomposition.
 *   b   - an _n_ element vector with the right-hand-side vector on entry,
 *         the solution on return; note that nothing is done if the
 *         left-hand-side matrix is singular.
 *
 * Returns:
 *   0 on success; -1 if the left-hand-side matrix is singular.
 */
SCOPE int LUSLV(const REAL MATRIX(a), const INTEGER n,
                const INTEGER VECTOR(pvt), REAL VECTOR(b));
#endif /* LUSLV */

#ifdef LUINV
/*
 * Function: LUINV
 *   Compute inverse of a matrix from its LU decomposition.
 *
 * Parameters:
 *   a   - the LU decomposition of an _n_ by _n_ matrix;
 *   lda - the leading dimension of _a_;
 *   n   - the dimension of the matrix;
 *   pvt - the _n_ + 1 element vector with the row permutations of the
 *         LU decomposition;
 *   ws  - an _n_ element workspace vector.
 *
 * Returns:
 *   0 on success; -1 if the left-hand-side matrix is singular.
 */
SCOPE int LUINV(const REAL MATRIX(a), const INTEGER n,
                const INTEGER VECTOR(pvt), REAL MATRIX(b),
                REAL VECTOR(ws));
#endif /* LUINV */

#ifdef LUDET
/*
 * Function: LUDET
 *   Compute the determinant of a matrix from its LU decomposition.
 *
 * Parameters:
 *   a   - the LU decomposition of an _n_ by _n_ matrix;
 *   lda - the leading dimension of _a_;
 *   n   - the dimension of the matrix;
 *   pvt - the _n_ + 1 element vector with the row permutations of the
 *         LU decomposition.
 *
 * Returns:
 *   The determinant of the matrix.
 */
SCOPE REAL LUDET(const REAL MATRIX(a), const INTEGER n,
                 const INTEGER VECTOR(pvt));
#endif /* LUDET */

/*----------------------------------------------------------------------*/

#define A(i,j) FETCH(a,i,j)
#define B(i,j) FETCH(b,i,j)

#ifdef LUDEC_ROW
SCOPE int LUDEC_ROW(REAL MATRIX(a), const INTEGER n, INTEGER VECTOR(pvt),
                    REAL VECTOR(ws))
{
  const REAL ZERO = 0;
  const REAL ONE = 1;
  REAL big, tmp, scl;
  INTEGER i, imax, j, k, sgn, kp1;

  /* Assume that matrix is singular. */
  pvt[n] = 0;

  /* Store scaling of rows. */
  for (i = 0 ; i < n; ++i) {
    big = ZERO;
    for (j = 0; j < n; ++j) {
      if ((tmp = A(i,j)) < ZERO) {
        tmp = -tmp;
      }
      if (tmp > big) {
        big = tmp;
      }
    }
    if (big == ZERO) {
      /* singular matrix */
      return -1;
    }
    ws[i] = ONE/big;
  }

  /* Perform the LU decomposition with pivoting. */
  sgn = 1;
  for (k = 0; k < n; ++k) {
    /* search for the largest pivot element */
    big = ZERO;
    imax = k;
    for (i = k; i < n; ++i) {
      if ((tmp = ws[i]*A(i,k)) < ZERO) {
        tmp = -tmp;
      }
      if (tmp > big) {
        big = tmp;
        imax = i;
      }
    }
    if (imax != k) {
      /* apply row permutation */
      for (j = 0; j < n; ++j) {
        tmp = A(imax,j);
        A(imax,j) = A(k,j);
        A(k,j) = tmp;
      }
      ws[imax] = ws[k];
      sgn = -sgn;
    }
    pvt[k] = imax;
    if (A(k,k) == ZERO) {
      /* singular matrix */
#ifdef TINY
      if (TINY != ZERO) {
        A(k,k) = TINY;
      } else {
        return -1;
      }
#else
      return -1;
#endif
    }
    if ((kp1 = k + 1) < n) {
      scl = ONE/A(k,k);
      for (i = kp1; i < n; ++i) {
        tmp = (A(i,k) *= scl);
        for (j = kp1; j < n; ++j) {
          A(i,j) -= tmp*A(k,j);
        }
      }
    }
  }

  /* Store parity of row permutations and report success. */
  pvt[n] = sgn;
  return 0;
}
#endif /* LUDEC_ROW */

#ifdef LUDEC_COL /* LU decomposition by Crout's algorithm. */
SCOPE int LUDEC_COL(REAL MATRIX(a), const INTEGER n, INTEGER VECTOR(pvt),
                    REAL VECTOR(ws))
{
  const REAL ZERO = 0;
  const REAL ONE = 1;
  REAL big, tmp, sum;
  INTEGER i, imax, j, k, sgn;

  /* Assume that matrix is singular. */
  pvt[n] = 0;

  /* Store scaling of rows. */
  for (i = 0 ; i < n; ++i) {
    big = ZERO;
    for (j = 0; j < n; ++j) {
      if ((tmp = A(i,j)) < ZERO) {
        tmp = -tmp;
      }
      if (tmp > big) {
        big = tmp;
      }
    }
    if (big == ZERO) {
      /* singular matrix */
      return -1;
    }
    ws[i] = ONE/big;
  }

  /* Perform the LU decomposition with pivoting. */
  sgn = 1;
  for (j = 0; j < n; ++j) {
    for (i = 0; i < j; ++i) {
      sum = A(i,j);
      for (k = 0; k < i; ++k) {
        sum -= A(i,k)*A(k,j);
      }
      A(i,j) = sum;
    }
    /* search for the largest pivot element */
    big = ZERO;
    imax = j;
    for (i = j; i < n; ++i) {
      sum = A(i,j);
      for (k = 0; k < j; ++k) {
        sum -= A(i,k)*A(k,j);
      }
      A(i,j) = sum;
      if ((tmp = ws[i]*sum) < ZERO) {
        tmp = -tmp;
      }
      if (tmp > big) { /* FIXME: was "greater or equal" */
        big = tmp;
        imax = i;
      }
    }
    if (imax != j) {
      /* apply row permutation */
      for (k = 0; k < n; ++k) {
        tmp = A(imax,k);
        A(imax,k) = A(j,k);
        A(j,k) = tmp;
      }
      ws[imax] = ws[j];
      sgn *= -1;
    }
    pvt[j] = imax;
    if (A(j,j) == ZERO) {
      /* singular matrix */
#ifdef TINY
      if (TINY != ZERO) {
        A(j,j) = TINY;
      } else {
        return -1;
      }
#else
      return -1;
#endif
    }
    if (j + 1 < n) {
      tmp = ONE/A(j,j);
      for (i = j + 1; i < n; ++i) {
        A(i,j) *= tmp;
      }
    }
  }

  /* Store parity of row permutations and report success. */
  pvt[n] = sgn;
  return 0;
}
#endif /* LUDEC_COL */

#ifdef LUDET
SCOPE REAL LUDET(const REAL MATRIX(a), const INTEGER n,
                 const INTEGER VECTOR(pvt))
{
  REAL det;
  INTEGER i;
  if (pvt[n] != -1 || pvt[n] != 1) return (REAL)0;
  det = pvt[n]; /* initialize with the sign of the permutation */
  for (i = 0; i < n; ++i) {
    det *= A(i,i);
  }
  return det;
}
#endif /* LUDET */

#ifdef LUINV
# ifndef LUSLV
#  error LUSLV must be defined for LUINV
# endif
/* FIXME: no needs for a workspace if B has correct layout */
SCOPE int LUINV(const REAL MATRIX(a), const INTEGER n,
                const INTEGER VECTOR(pvt), REAL MATRIX(b),
                REAL VECTOR(ws))
{
  const REAL ZERO = 0;
  const REAL ONE = 1;
  INTEGER i, j;

  if (pvt[n] == 0) return -1; /* LHS matrix was singular */
  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      ws[i] = ZERO;
    }
    ws[j] = ONE;
    LUSLV(PUT_MATRIX(a), n, pvt, ws);
    for (i = 0; i < n; ++i) {
      B(i,j) = ws[i];
    }
  }
  return 0;
}
#endif /* LUINV */

#ifdef LUSLV
SCOPE int LUSLV(const REAL MATRIX(a), const INTEGER n,
                const INTEGER VECTOR(pvt), REAL VECTOR(b))
{
  REAL sum, tmp;
  INTEGER i, j, k, ip;

  /* Report failure if LHS matrix was singular. */
  if (pvt[n] == 0) return -1;

  /* Do the forward-substitution, unscrambling the row permutations and using
     variable K to store the index of the first non-vanishing element of B. */
  k = -1;
  for (i = 0; i < n; ++i) {
    ip = pvt[i];
    sum = b[ip];
    b[ip] = b[i];
    if (k >= 0) {
      for (j = k; j < i; ++j) {
        sum -= A(i,j)*b[j];
      }
    } else if (sum) {
      k = i;
    }
    b[i] = sum;
  }

  /* Do the back-substitution. */
  for (i = n - 1; i >= 0; --i) {
    sum = b[i];
    tmp = A(i,i);
    for (j = i + 1; j < n; ++j) {
      sum -= A(i,j)*b[j];
    }
    b[i] = sum/tmp;
  }

  /* Report success. */
  return 0;
}
#endif /* LUSLV */

/* Undefine private macros to allow for multiple inclusion. */
#undef A
#undef B
#if 0
# undef VECTOR
# undef MATRIX
# undef PUT_MATRIX
# undef FETCH
# undef ABS
#endif

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: nil
 * c-basic-offset: 2
 * fill-column: 78
 * coding: utf-8
 * End:
 */
