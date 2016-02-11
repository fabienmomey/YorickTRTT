/*
 * trtt_xform3d.c --
 *
 * Implementation of 3D transformations.
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

#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#include "trtt_xform3d.h"
#include "trtt_common.h"

/* It is possible to change the storage with the following macros.
   However note that the documentation of the code assumes row-major
   storage order. */
#define TRTT_XFORM3D_ROW_MAJOR
#undef  TRTT_XFORM3D_COLUMN_MAJOR

struct _trtt_xform3d {
  /*
   * A General 3D coordinate transform matrix is stored as 12
   * coefficients in row-major order:
   *
   *     / A11  A12  A13  A14 \   / A[0]  A[1]  A[2]  A[3]  \
   * A = | A21  A22  A23  A24 | = | A[4]  A[5]  A[6]  A[7]  |
   *     | A31  A32  A33  A34 |   | A[8]  A[9]  A[10] A[11] |
   *     \ A41  A42  A43  A44 /   \ 0     0     0     1     /
   */
  double a[12];
};

#define STATEMENT(code) do { code; } while (0)

#define RETURN_SUCCESS STATEMENT(return TRTT_SUCCESS)
#define RETURN_FAILURE STATEMENT(return TRTT_FAILURE)

#define TRTT_XFORM3D_COPY_COEFS(dst, src)           \
  STATEMENT(dst[ 0] = src[ 0];                  \
            dst[ 1] = src[ 1];                  \
            dst[ 2] = src[ 2];                  \
            dst[ 3] = src[ 3];                  \
            dst[ 4] = src[ 4];                  \
            dst[ 5] = src[ 5];                  \
            dst[ 6] = src[ 6];                  \
            dst[ 7] = src[ 7];                  \
            dst[ 8] = src[ 8];                  \
            dst[ 9] = src[ 9];                  \
            dst[10] = src[10];                  \
            dst[11] = src[11])

/* The macro _XFORM3D_FETCH is to access the element (i,j) of a matrix A with
   leading dimension LDA.  The macro _XFORM3D_FETCH_XFORM is to access the
   element of a 3D coordinate transform matrix. */
#ifdef TRTT_XFORM3D_ROW_MAJOR /* row-major storage order */
# define TRTT_XFORM3D_FETCH(a,lda,i,j)     ((a)[(i)*(lda) + (j)])
# define TRTT_XFORM3D_FETCH_XFORM(a,i,j)  TRTT_XFORM3D_FETCH(a,4,i,j)
#else /* column-major storage order */
# define TRTT_XFORM3D_FETCH(a,lda,i,j)     ((a)[(j)*(lda) + (i)])
# define TRTT_XFORM3D_FETCH_XFORM(a,i,j)  TRTT_XFORM3D_FETCH(a,3,i,j)
#endif

/* Stores the identity into a matrix A. */
#ifdef TRTT_XFORM3D_ROW_MAJOR /* row-major storage order */
# define TRTT_XFORM3D_RESET(a) STATEMENT(                   \
    a[0] = 1.0;  a[1] = 0.0;  a[2] = 0.0;  a[3] = 0.0;  \
    a[4] = 0.0;  a[5] = 1.0;  a[6] = 0.0;  a[7] = 0.0;  \
    a[8] = 0.0;  a[9] = 0.0; a[10] = 1.0; a[11] = 0.0)
#else /* column-major storage order */
# define TRTT_XFORM3D_RESET(a) STATEMENT(      \
    a[0] = 1.0;  a[1] = 0.0;  a[2] = 0.0;  \
    a[3] = 0.0;  a[4] = 1.0;  a[5] = 0.0;  \
    a[6] = 0.0;  a[7] = 0.0;  a[8] = 1.0;  \
    a[9] = 0.0; a[10] = 0.0; a[11] = 0.0)
#endif

/* Code for matrix dot product: C = A.B -- C must not be equal to
   A nor B */
#ifdef TRTT_XFORM3D_ROW_MAJOR /* row-major storage order */
# define TRTT_XFORM3D_XFORM_DOT_XFORM(c, a, b) STATEMENT(           \
     c[0] = a[0]*b[0] + a[1]*b[4] +  a[2]*b[8];                 \
     c[1] = a[0]*b[1] + a[1]*b[5] +  a[2]*b[9];                 \
     c[2] = a[0]*b[2] + a[1]*b[6] +  a[2]*b[10];                \
     c[3] = a[0]*b[3] + a[1]*b[7] +  a[2]*b[11] + a[3];         \
     c[4] = a[4]*b[0] + a[5]*b[4] +  a[6]*b[8];                 \
     c[5] = a[4]*b[1] + a[5]*b[5] +  a[6]*b[9];                 \
     c[6] = a[4]*b[2] + a[5]*b[6] +  a[6]*b[10];                \
     c[7] = a[4]*b[3] + a[5]*b[7] +  a[6]*b[11] + a[7];         \
     c[8] = a[8]*b[0] + a[9]*b[4] + a[10]*b[8];                 \
     c[9] = a[8]*b[1] + a[9]*b[5] + a[10]*b[9];                 \
    c[10] = a[8]*b[2] + a[9]*b[6] + a[10]*b[10];                \
    c[11] = a[8]*b[3] + a[9]*b[7] + a[10]*b[11] + a[11])
#else /* column-major storage order */
# define TRTT_XFORM3D_XFORM_DOT_XFORM(c, a, b) STATEMENT(           \
     c[0] = a[0]*b[0] + a[3]*b[1]  + a[6]*b[2];                 \
     c[1] = a[1]*b[0] + a[4]*b[1]  + a[7]*b[2];                 \
     c[2] = a[2]*b[0] + a[5]*b[1]  + a[8]*b[2];                 \
     c[3] = a[0]*b[3] + a[3]*b[4]  + a[6]*b[5];                 \
     c[4] = a[1]*b[3] + a[4]*b[4]  + a[7]*b[5];                 \
     c[5] = a[2]*b[3] + a[5]*b[4]  + a[8]*b[5];                 \
     c[6] = a[0]*b[6] + a[3]*b[7]  + a[6]*b[8];                 \
     c[7] = a[1]*b[6] + a[4]*b[7]  + a[7]*b[8];                 \
     c[8] = a[2]*b[6] + a[5]*b[7]  + a[8]*b[8];                 \
     c[9] = a[0]*b[9] + a[3]*b[10] + a[6]*b[11] + a[9];         \
    c[10] = a[1]*b[9] + a[4]*b[10] + a[7]*b[11] + a[10];        \
    c[11] = a[2]*b[9] + a[5]*b[10] + a[8]*b[11] + a[11])
#endif

/* Multiply matrix A and vector V to produce W = A.V -- W must
   not be equal to V if X=V[0], Y=V[1], and Z=V[2] */
#ifdef TRTT_XFORM3D_ROW_MAJOR /* row-major storage order */
# define TRTT_XFORM3D_XFORM_DOT_VECTOR(w, a, x, y, z) STATEMENT(    \
  w[0] = a[0]*x + a[1]*y +  a[2]*z +  a[3];                     \
  w[1] = a[4]*x + a[5]*y +  a[6]*z +  a[7];                     \
  w[3] = a[8]*x + a[9]*y + a[10]*z + a[11])
#else /* column-major storage order */
# define TRTT_XFORM3D_XFORM_DOT_VECTOR(w, a, x, y, z) STATEMENT(    \
  w[0] = a[0]*x + a[3]*y + a[6]*z +  a[9];                      \
  w[1] = a[1]*x + a[4]*y + a[7]*z + a[10];                      \
  w[3] = a[2]*x + a[5]*y + a[8]*z + a[11])
#endif

/* Code for dot product: B = R.A with R a 3x3 rotation matrix and A a general
   3D coordinate transform matrix -- operation must not be done in-place */
#ifdef TRTT_XFORM3D_ROW_MAJOR /* row-major storage order */
# define TRTT_XFORM3D_ROTATION_DOT_XFORM(b, r, a) STATEMENT(       \
     b[0] = r[0]*a[0] + r[1]*a[4] + r[2]*a[8];                  \
     b[1] = r[0]*a[1] + r[1]*a[5] + r[2]*a[9];                  \
     b[2] = r[0]*a[2] + r[1]*a[6] + r[2]*a[10];                 \
     b[3] = r[0]*a[3] + r[1]*a[7] + r[2]*a[11];                 \
     b[4] = r[3]*a[0] + r[4]*a[4] + r[5]*a[8];                  \
     b[5] = r[3]*a[1] + r[4]*a[5] + r[5]*a[9];                  \
     b[6] = r[3]*a[2] + r[4]*a[6] + r[5]*a[10];                 \
     b[7] = r[3]*a[3] + r[4]*a[7] + r[5]*a[11];                 \
     b[8] = r[6]*a[0] + r[7]*a[4] + r[8]*a[8];                  \
     b[9] = r[6]*a[1] + r[7]*a[5] + r[8]*a[9];                  \
    b[10] = r[6]*a[2] + r[7]*a[6] + r[8]*a[10];                 \
    b[11] = r[6]*a[3] + r[7]*a[7] + r[8]*a[11])
#else /* column-major storage order */
# define TRTT_XFORM3D_ROTATION_DOT_XFORM(b, r, a) STATEMENT(       \
     b[0] = r[0]*a[0] + r[3]*a[1]  + r[6]*a[2];                 \
     b[1] = r[1]*a[0] + r[4]*a[1]  + r[7]*a[2];                 \
     b[2] = r[2]*a[0] + r[5]*a[1]  + r[8]*a[2];                 \
     b[3] = r[0]*a[3] + r[3]*a[4]  + r[6]*a[5];                 \
     b[4] = r[1]*a[3] + r[4]*a[4]  + r[7]*a[5];                 \
     b[5] = r[2]*a[3] + r[5]*a[4]  + r[8]*a[5];                 \
     b[6] = r[0]*a[6] + r[3]*a[7]  + r[6]*a[8];                 \
     b[7] = r[1]*a[6] + r[4]*a[7]  + r[7]*a[8];                 \
     b[8] = r[2]*a[6] + r[5]*a[7]  + r[8]*a[8];                 \
     b[9] = r[0]*a[9] + r[3]*a[10] + r[6]*a[11];                \
    b[10] = r[1]*a[9] + r[4]*a[10] + r[7]*a[11];                \
    b[11] = r[2]*a[9] + r[5]*a[10] + r[8]*a[11])
#endif

trtt_xform3d_t *trtt_xform3d_new()
{
  trtt_xform3d_t *xform;
  double *a;

  xform = TRTT_NEW0(trtt_xform3d_t);
  if (xform == NULL) {
    return NULL;
  }
  a = xform->a;
  TRTT_XFORM3D_RESET(a);
  return xform;
}

int trtt_xform3d_init(trtt_xform3d_t *xform)
{
  double *a;
  a = xform->a;
  TRTT_XFORM3D_RESET(a);
  RETURN_SUCCESS;
}

void trtt_xform3d_delete(trtt_xform3d_t *xform)
{
  if (xform != NULL) {
    free(xform);
  }
}

int trtt_xform3d_reset(trtt_xform3d_t *xform)
{
  double *a;

  a = xform->a;
  TRTT_XFORM3D_RESET(a);
  RETURN_SUCCESS;
}

int trtt_xform3d_copy(trtt_xform3d_t *dst,
                      const trtt_xform3d_t *src)
{
  if (src != dst) {
    memcpy(dst, src, sizeof(*src));
  }
  RETURN_SUCCESS;
}

int trtt_xform3d_apply(double dst[3],
                       const trtt_xform3d_t *xform,
                       const double src[3])
{
  double src_x, src_y, src_z;
  const double *a;
  
  /* copy the input coordinates so that operation can be done "in-place" */
  src_x = src[0];
  src_y = src[1];
  src_z = src[2];
  a = xform->a;
  TRTT_XFORM3D_XFORM_DOT_VECTOR(dst, a, src_x, src_y, src_z);
  RETURN_SUCCESS;
}

int trtt_xform3d_combine(trtt_xform3d_t *dst,
                         const trtt_xform3d_t *left,
                         const trtt_xform3d_t *right)
{
  /* A temporary array is used in case operation is done "in-place" -- DST and
     LEFT or RIGHT are the same. */
  double t[12], *c;
  const double *a, *b;
  
  a = left->a;
  b = right->a;
  c = dst->a;
  if (c == a || c == b) {
    TRTT_XFORM3D_XFORM_DOT_XFORM(t, a, b);
    TRTT_XFORM3D_COPY_COEFS(c, t);
  } else {
    TRTT_XFORM3D_XFORM_DOT_XFORM(c, a, b);
  }

  RETURN_SUCCESS;
}

int trtt_xform3d_move(trtt_xform3d_t *dst,
                      double tx,
                      double ty,
                      double tz,
                      const trtt_xform3d_t *src)
{
  double *b;
  const double *a;
  
  a = src->a;
  b = dst->a;
  if (b != a) {
    TRTT_XFORM3D_COPY_COEFS(b, a);
  }
#ifdef TRTT_XFORM3D_ROW_MAJOR /* row-major storage order */
  b[ 3] += tx;
  b[ 7] += ty;
  b[11] += tz;
#else /* column-major storage order */
  b[ 9] += tx;
  b[10] += ty;
  b[11] += tz;
#endif

  RETURN_SUCCESS;
}

int trtt_xform3d_rotate(trtt_xform3d_t *dst,
                        double qw,
                        double qx,
                        double qy,
                        double qz,
                        const trtt_xform3d_t *src)
{
  double aww, axx, ayy, azz, nrm2, scl;
  double awx, awy, awz, axy, axz, ayz;
  double *b, r[9], t[12];
  const double *a;
  
  a = src->a;
  b = dst->a;

  /* Compute the elements of the rotation matrix (starting with the
     diagonal). */
#define R(i,j) TRTT_XFORM3D_FETCH(r,3,i,j)
  aww = qw*qw;
  axx = qx*qx;
  ayy = qy*qy;
  azz = qz*qz;
  nrm2 = (aww + axx) + (ayy + azz);
  if (nrm2 == 1.0) {
    /* Quaternion is normalized. */
    R(0,0) = (aww + axx) - (ayy + azz);
    R(1,1) = (aww + ayy) - (axx + azz);
    R(2,2) = (aww + azz) - (axx + ayy);
    scl = 2.0;
  } else if (nrm2 > 0.0) {
    /* Normalize the quaternion. */
    scl = 1.0/nrm2;
    R(0,0) = scl*((aww + axx) - (ayy + azz));
    R(1,1) = scl*((aww + ayy) - (axx + azz));
    R(2,2) = scl*((aww + azz) - (axx + ayy));
    scl += scl;
  } else {
    /* Quaternion cannot be normalized. */
    errno = EINVAL;
    RETURN_FAILURE; 
  }
  awx = scl*qw*qx;
  awy = scl*qw*qy;
  awz = scl*qw*qz;
  axy = scl*qx*qy;
  axz = scl*qx*qz;
  ayz = scl*qy*qz;
  R(0,1) = axy - awz;
  R(0,2) = axz + awy;
  R(1,0) = axy + awz;
  R(1,2) = ayz - awx;
  R(2,0) = axz - awy;
  R(2,1) = ayz + awx;
#undef R
  
  /* Multiply the rotation and the transform matrices. */
  if (a != b) {
    TRTT_XFORM3D_ROTATION_DOT_XFORM(b, r, a);
  } else {
    TRTT_XFORM3D_ROTATION_DOT_XFORM(t, r, a);
    TRTT_XFORM3D_COPY_COEFS(b, t);
  }

  RETURN_SUCCESS;
}

int trtt_xform3d_scale(trtt_xform3d_t *dst,
                       double x_scl,
                       double y_scl,
                       double z_scl,
                       const trtt_xform3d_t *src)
{
  const double *a;
  double *b;

  if (x_scl != 1.0 || y_scl != 1.0 || z_scl != 1.0 || dst != src) {
    a = src->a;
    b = dst->a;
#define XMAGNIFY(i) b[i] = x_scl*a[i]
    XMAGNIFY(0);
    XMAGNIFY(1);
    XMAGNIFY(2);
    XMAGNIFY(3);
#undef XMAGNIFY
#define YMAGNIFY(i) b[i] = y_scl*a[i]
    YMAGNIFY(4);
    YMAGNIFY(5);
    YMAGNIFY(6);
    YMAGNIFY(7);
#undef YMAGNIFY
#define ZMAGNIFY(i) b[i] = z_scl*a[i]
    ZMAGNIFY(8);
    ZMAGNIFY(9);
    ZMAGNIFY(10);
    ZMAGNIFY(11);
#undef ZMAGNIFY
  }

  RETURN_SUCCESS;
}

int trtt_xform3d_get_coefficients(double dst[12],
                                  const trtt_xform3d_t *src)
{
  const double *a;

  a = src->a;
  TRTT_XFORM3D_COPY_COEFS(dst, a);

  RETURN_SUCCESS;
}

int trtt_xform3d_set_coefficients(trtt_xform3d_t *dst, const double src[12])
{
  double *a;

  a = dst->a;
  TRTT_XFORM3D_COPY_COEFS(a, src);
  RETURN_SUCCESS;
}

/* Define macros and include source for LU decomposition. */
#ifdef TRTT_XFORM3D_ROW_MAJOR
# define STYLE   3       /* row-major storage order */
#else
# define STYLE   2       /* column-major storage order */
#endif
#define SCOPE   static
#define REAL    double
#define INTEGER int
#define LUDEC   ludec
#define LUSLV   luslv
#define LUINV   luinv
#undef  LUDET   /* no needs to compute the determinant */
#include "gauss.c"

int trtt_xform3d_invert(trtt_xform3d_t *dst,
                        const trtt_xform3d_t *src)
{
  const int n = 4;
  double lu[n*n], ainv[n*n], ws[n], *b;
  const double *a;
  int pvt[n+1];

  a = src->a;
  TRTT_XFORM3D_COPY_COEFS(lu, a);
  lu[12] = 0.0;
  lu[13] = 0.0;
  lu[14] = 0.0;
  lu[15] = 1.0; /* FIXME: scaling? */  
  if (LUDEC(lu, n, pvt, ws) != 0 || LUINV(lu, n, pvt, ainv, ws) != 0) {
    /* the matrix is singular */
    RETURN_FAILURE;
  }
  b = dst->a;
  TRTT_XFORM3D_COPY_COEFS(b, ainv);
  RETURN_SUCCESS;
}

void trtt_xform3d_print(FILE *output, const trtt_xform3d_t *xform)
{
  const double *a;
  int i, j;

  if (output == NULL) {
    output = stdout;
  }
  a = xform->a;
#define A(i,j) TRTT_XFORM3D_FETCH_XFORM(a,i,j)
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 4; ++j) {
      fprintf(output, "  %12.3e", A(i,j));
    }
    fputs("\n", output);
  }
#undef A
}

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
