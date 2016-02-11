/*
 * LongFessler2D.c -
 *
 * Implements Long & Fessler model for 2D tomography.
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <yapi.h>

#include "LongFessler2D.h"


#include <sys/types.h>
#include <unistd.h>

void Y_getpid(int argc)
{
  ypush_long(getpid());
}

#define MAX(a, b) ((a) >= (b) ? (a) : (b))
#define MIN(a, b) ((a) <= (b) ? (a) : (b))

static double min(double a, double b) { return MIN(a, b); }
static double max(double a, double b) { return MAX(a, b); }

static double integTriangle(double a, double b, double c)
{
  if (a >= b) return 0.0;
  return (b + a)*(b - a)/(c + c);
}

static double integRectangle(double a, double b)
{
  if (a >= b) return 0.0;
  return (b - a);
}

void
lf2d_direct(const lf2d_parameters_t* p, /* parameters */
            double beta,                /* direction of source w.r.t. y-axis */
            double detectorOffset,      /* offset of detector w.r.t. to its
                                           geometrical center */
            const double voxel[],       /* object voxels */
            double pixel[],             /* detector pixels */
            int clr)                    /* clear destination before
                                           integration */
{
  long i1, objectDimension1 = p->objectDimension1;
  long i2, objectDimension2 = p->objectDimension2;
  long k, kFirst, kLast, detectorDimension = p->detectorDimension;
  double voxelSize = p->voxelSize;
  double pixelSize = p->pixelSize;
  double pixelStep = p->pixelStep;
  double sourceDistance = p->sourceDistance;
  double detectorDistance = p->detectorDistance;

  double cosBeta = cos(beta);
  double sinBeta = sin(beta);
  double offset_y = cosBeta*(p->objectOffset1)+sinBeta*(p->objectOffset2);
  double offset_x = -sinBeta*(p->objectOffset1)+cosBeta*(p->objectOffset2);

  /* Position of the center of first voxel. */
  double x0 = offset_x - 0.5*(objectDimension1 - 1)*voxelSize;
  double y0 = offset_y - 0.5*(objectDimension2 - 1)*voxelSize;

  /* Parameters for computing the projection of a voxel on the detector.
     Projected coordinates are computed in pixel units. */
  double offset = detectorOffset/pixelStep + 0.5*(detectorDimension - 1);
  double scale = (sourceDistance + detectorDistance)/pixelStep;
  double size = pixelSize/pixelStep;

  /* Clear destination array (for integration of result). */
  if (clr) {
    memset(pixel, 0, sizeof(double)*detectorDimension);
  }

  /* Loop over all voxels. */
  for (i2 = 0; i2 < objectDimension2; ++i2) {
    double y = y0 + i2*voxelSize;  /* abscissa of voxel center */
    double y1 = y - 0.5*voxelSize;
    double y2 = y + 0.5*voxelSize;
    for (i1 = 0; i1 < objectDimension1; ++i1) {
      double x = x0 + i1*voxelSize;  /* ordinate of voxel center */
      double x1 = x - 0.5*voxelSize;
      double x2 = x + 0.5*voxelSize;
      double t0, t1, t2, t3;
      double alpha, q, qx, qy;

      /* Compute coordinates along detector axis of voxel at (x,y) in pixel
         units. */
#     define S(x,y) (scale*((y)*sinBeta + (x)*cosBeta) \
                     /(sourceDistance - ((y)*cosBeta - (x)*sinBeta)) + offset)
      t0 = S(x1,y1);
      t1 = S(x2,y1);
      t2 = S(x1,y2);
      t3 = S(x2,y2);
#     undef S

      /* Sort coordinates. */
#     define SORT(a,b) if (a > b) { double tmp = a; a = b; b = tmp; } else
      SORT(t0,t1);
      SORT(t2,t3);
      SORT(t0,t2);
      SORT(t1,t3);
      SORT(t1,t2);
#     undef SORT

      /* Normalization factor ALPHA is the trapezoid height times the pixel
         size (because the integral of the footprint is computed in pixel
         units).  (qx,qy) are the coordinates of the source-voxel vector. */
      qx = x + sourceDistance*sinBeta;
      qy = y - sourceDistance*cosBeta;
      q = (fabs(qy) > fabs(qx) ? qx/qy : qy/qx);
      alpha = (sqrt(1.0 + q*q)*(voxelSize*pixelSize)
               *voxel[i1 + i2*objectDimension1]);

      /* Index of first/last pixel impacted by this voxel. */
      kFirst = (long)floor(t0 - 0.5*size);
      if (kFirst < 0) kFirst = 0;
      kLast = (long)ceil(t3 + 0.5*size);
      if (kLast >= detectorDimension) kLast = detectorDimension - 1;

      /* Integrate projected voxel profile (approximated by a trapezoid) on the
         impacted pixels.  GAMMA is the integral of the trapezoid over the
         pixel (divided by the pixel size). */
      for (k = kFirst; k <= kLast; ++k) {
	double s1 = k - 0.5*size;
	double s2 = s1 + size;
	double gamma = (integTriangle(max(s1, t0) - t0,
                                      min(s2, t1) - t0, t1 - t0) +
                        integRectangle(max(s1, t1), min(s2, t2)) +
                        integTriangle(t3 - min(s2, t3),
                                      t3 - max(s1, t2), t3 - t2));
	pixel[k] += gamma*alpha;
      }
    }
  }
}

void
lf2d_line_integ_direct(const lf2d_parameters_t* p, /* parameters */
                       double beta,                /* direction of source w.r.t. y-axis */
                       double detectorOffset,      /* offset of detector w.r.t. to its
                                           geometrical center */
                       const double voxel[],       /* object voxels */
                       double pixel[],             /* detector pixels */
                       int clr)                    /* clear destination before
                                           integration */
{
  long i1, objectDimension1 = p->objectDimension1;
  long i2, objectDimension2 = p->objectDimension2;
  long k, kFirst, kLast, detectorDimension = p->detectorDimension;
  double voxelSize = p->voxelSize;
  double pixelSize = p->pixelSize;
  double pixelStep = p->pixelStep;
  double sourceDistance = p->sourceDistance;
  double detectorDistance = p->detectorDistance;

  double cosBeta = cos(beta);
  double sinBeta = sin(beta);
  double offset_y = cosBeta*(p->objectOffset1)+sinBeta*(p->objectOffset2);
  double offset_x = -sinBeta*(p->objectOffset1)+cosBeta*(p->objectOffset2);

  /* Position of the center of first voxel. */
  double x0 = offset_x - 0.5*(objectDimension1 - 1)*voxelSize;
  double y0 = offset_y - 0.5*(objectDimension2 - 1)*voxelSize;

  /* Parameters for computing the projection of a voxel on the detector.
     Projected coordinates are computed in pixel units. */
  double offset = detectorOffset/pixelStep + 0.5*(detectorDimension - 1);
  double scale = (sourceDistance + detectorDistance)/pixelStep;
  double size = pixelSize/pixelStep;

  /* Clear destination array (for integration of result). */
  if (clr) {
    memset(pixel, 0, sizeof(double)*detectorDimension);
  }

  /* Loop over all voxels. */
  for (i2 = 0; i2 < objectDimension2; ++i2) {
    double y = y0 + i2*voxelSize;  /* abscissa of voxel center */
    double y1 = y - 0.5*voxelSize;
    double y2 = y + 0.5*voxelSize;
    for (i1 = 0; i1 < objectDimension1; ++i1) {
      double x = x0 + i1*voxelSize;  /* ordinate of voxel center */
      double x1 = x - 0.5*voxelSize;
      double x2 = x + 0.5*voxelSize;
      double t0, t1, t2, t3;
      double alpha, q, qx, qy;

      /* Compute coordinates along detector axis of voxel at (x,y) in pixel
         units. */
#     define S(x,y) (scale*((y)*sinBeta + (x)*cosBeta) \
                     /(sourceDistance - ((y)*cosBeta - (x)*sinBeta)) + offset)
      t0 = S(x1,y1);
      t1 = S(x2,y1);
      t2 = S(x1,y2);
      t3 = S(x2,y2);
#     undef S

      /* Sort coordinates. */
#     define SORT(a,b) if (a > b) { double tmp = a; a = b; b = tmp; } else
      SORT(t0,t1);
      SORT(t2,t3);
      SORT(t0,t2);
      SORT(t1,t3);
      SORT(t1,t2);
#     undef SORT

      /* Normalization factor ALPHA is the trapezoid height times the pixel
         size (because the integral of the footprint is computed in pixel
         units).  (qx,qy) are the coordinates of the source-voxel vector. */
      qx = x + sourceDistance*sinBeta;
      qy = y - sourceDistance*cosBeta;
      q = (fabs(qy) > fabs(qx) ? qx/qy : qy/qx);
      alpha = (sqrt(1.0 + q*q)*voxelSize
               *voxel[i1 + i2*objectDimension1]);

      /* Index of first/last pixel impacted by this voxel. */
      kFirst = (long)floor(t0 - 0.5*size);
      if (kFirst < 0) kFirst = 0;
      kLast = (long)ceil(t3 + 0.5*size);
      if (kLast >= detectorDimension) kLast = detectorDimension - 1;

      /* Integrate projected voxel profile (approximated by a trapezoid) on the
         impacted pixels.  GAMMA is the integral of the trapezoid over the
         pixel (divided by the pixel size). */
      for (k = kFirst; k <= kLast; ++k) {
        if (k>=t0 && k<t1) {
          pixel[k] += alpha*(k-t0)/(t1-t0);
        } else if (k>=t1 && k<t2) {
          pixel[k] += alpha;     
        } else if (k>=t2 && k<t3) {
          pixel[k] += alpha*(k-t3)/(t2-t3);
        } else {
          pixel[k] += 0;
        }      
	/* double s1 = k - 0.5*size; */
	/* double s2 = s1 + size; */
	/* double gamma = (integTriangle(max(s1, t0) - t0, */
        /*                               min(s2, t1) - t0, t1 - t0) + */
        /*                 integRectangle(max(s1, t1), min(s2, t2)) + */
        /*                 integTriangle(t3 - min(s2, t3), */
        /*                               t3 - max(s1, t2), t3 - t2)); */
	/* pixel[k] += gamma*alpha */;
      }
    }
  }
}

void
lf2d_adjoint(const lf2d_parameters_t* p, /* parameters */
             double beta,                /* direction of source
                                            w.r.t. y-axis */
             double detectorOffset,      /* offset of detector w.r.t. to its
                                            geometrical center */
             const double pixel[],       /* detector pixels */
             double voxel[],             /* object voxels */
             int clr)                    /* clear destination before
                                            integration */
{
  long i1, objectDimension1 = p->objectDimension1;
  long i2, objectDimension2 = p->objectDimension2;
  long k, kFirst, kLast, detectorDimension = p->detectorDimension;
  double voxelSize = p->voxelSize;
  double pixelSize = p->pixelSize;
  double pixelStep = p->pixelStep;
  double sourceDistance = p->sourceDistance;
  double detectorDistance = p->detectorDistance;

  double cosBeta = cos(beta);
  double sinBeta = sin(beta);
  double offset_y = cosBeta*(p->objectOffset1)+sinBeta*(p->objectOffset2);
  double offset_x = -sinBeta*(p->objectOffset1)+cosBeta*(p->objectOffset2);

  /* Position of the center of first voxel. */
  double x0 = offset_x - 0.5*(objectDimension1 - 1)*voxelSize;
  double y0 = offset_y - 0.5*(objectDimension2 - 1)*voxelSize;

  /* Parameters for computing the projection of a voxel on the detector.
     Projected coordinates are computed in pixel units. */
  double offset = detectorOffset/pixelStep + 0.5*(detectorDimension - 1);
  double scale = (sourceDistance + detectorDistance)/pixelStep;
  double size = pixelSize/pixelStep;

  /* Clear destination array (for integration of result). */
  if (clr) {
    memset(voxel, 0, sizeof(double)*objectDimension1*objectDimension2);
  }

  /* Loop over all voxels. */
  for (i2 = 0; i2 < objectDimension2; ++i2) {
    double y = y0 + i2*voxelSize;  /* abscissa of voxel center */
    double y1 = y - 0.5*voxelSize;
    double y2 = y + 0.5*voxelSize;
    for (i1 = 0; i1 < objectDimension1; ++i1) {
      double x = x0 + i1*voxelSize;  /* ordinate of voxel center */
      double x1 = x - 0.5*voxelSize;
      double x2 = x + 0.5*voxelSize;
      double t0, t1, t2, t3, sum;
      double alpha, q, qx, qy;

      /* Compute coordinates along detector axis of voxel at (x,y) in pixel
         units. */
#     define S(x,y) (scale*((y)*sinBeta + (x)*cosBeta) \
                     /(sourceDistance - ((y)*cosBeta - (x)*sinBeta)) + offset)
      t0 = S(x1,y1);
      t1 = S(x2,y1);
      t2 = S(x1,y2);
      t3 = S(x2,y2);
#     undef S

      /* Sort coordinates. */
#     define SORT(a,b) if (a > b) { double tmp = a; a = b; b = tmp; } else
      SORT(t0,t1);
      SORT(t2,t3);
      SORT(t0,t2);
      SORT(t1,t3);
      SORT(t1,t2);
#     undef SORT

      /* Normalization factor.  (qx,qy) are the coordinates of the source-voxel
         vector. */
      qx = x + sinBeta*sourceDistance;
      qy = y - cosBeta*sourceDistance;
      q = (fabs(qy) > fabs(qx) ? qx/qy : qy/qx);
      alpha = sqrt(1.0 + q*q)*(voxelSize*pixelSize);

      /* Index of first/last pixel impacted by this voxel. */
      kFirst = (long)floor(t0 - 0.5*size);
      if (kFirst < 0) kFirst = 0;
      kLast = (long)ceil(t3 + 0.5*size);
      if (kLast >= detectorDimension) kLast = detectorDimension - 1;

      /* Integrate projected voxel profile (approximated by a trapezoid) on the
         impacted pixels.  GAMMA is the integral of the trapezoid over the
         pixel (divided by the pixel size). */
      sum = 0.0;
      for (k = kFirst; k <= kLast; ++k) {
	double s1 = k - 0.5*size;
	double s2 = s1 + size;
	double gamma = (integTriangle(max(s1, t0) - t0,
                                      min(s2, t1) - t0, t1 - t0) +
                        integRectangle(max(s1, t1), min(s2, t2)) +
                        integTriangle(t3 - min(s2, t3),
                                      t3 - max(s1, t2), t3 - t2));
        sum += alpha*gamma*pixel[k];
      }
      voxel[i1 + i2*objectDimension1] += sum;
    }
  }
}

void
lf2d_line_integ_adjoint(const lf2d_parameters_t* p, /* parameters */
                        double beta,                /* direction of source
                                            w.r.t. y-axis */
                        double detectorOffset,      /* offset of detector w.r.t. to its
                                            geometrical center */
                        const double pixel[],       /* detector pixels */
                        double voxel[],             /* object voxels */
                        int clr)                    /* clear destination before
                                            integration */
{
  long i1, objectDimension1 = p->objectDimension1;
  long i2, objectDimension2 = p->objectDimension2;
  long k, kFirst, kLast, detectorDimension = p->detectorDimension;
  double voxelSize = p->voxelSize;
  double pixelSize = p->pixelSize;
  double pixelStep = p->pixelStep;
  double sourceDistance = p->sourceDistance;
  double detectorDistance = p->detectorDistance;

  double cosBeta = cos(beta);
  double sinBeta = sin(beta);
  double offset_y = cosBeta*(p->objectOffset1)+sinBeta*(p->objectOffset2);
  double offset_x = -sinBeta*(p->objectOffset1)+cosBeta*(p->objectOffset2);

  /* Position of the center of first voxel. */
  double x0 = offset_x - 0.5*(objectDimension1 - 1)*voxelSize;
  double y0 = offset_y - 0.5*(objectDimension2 - 1)*voxelSize;

  /* Parameters for computing the projection of a voxel on the detector.
     Projected coordinates are computed in pixel units. */
  double offset = detectorOffset/pixelStep + 0.5*(detectorDimension - 1);
  double scale = (sourceDistance + detectorDistance)/pixelStep;
  double size = pixelSize/pixelStep;

  /* Clear destination array (for integration of result). */
  if (clr) {
    memset(voxel, 0, sizeof(double)*objectDimension1*objectDimension2);
  }

  /* Loop over all voxels. */
  for (i2 = 0; i2 < objectDimension2; ++i2) {
    double y = y0 + i2*voxelSize;  /* abscissa of voxel center */
    double y1 = y - 0.5*voxelSize;
    double y2 = y + 0.5*voxelSize;
    for (i1 = 0; i1 < objectDimension1; ++i1) {
      double x = x0 + i1*voxelSize;  /* ordinate of voxel center */
      double x1 = x - 0.5*voxelSize;
      double x2 = x + 0.5*voxelSize;
      double t0, t1, t2, t3, sum;
      double alpha, q, qx, qy;

      /* Compute coordinates along detector axis of voxel at (x,y) in pixel
         units. */
#     define S(x,y) (scale*((y)*sinBeta + (x)*cosBeta) \
                     /(sourceDistance - ((y)*cosBeta - (x)*sinBeta)) + offset)
      t0 = S(x1,y1);
      t1 = S(x2,y1);
      t2 = S(x1,y2);
      t3 = S(x2,y2);
#     undef S

      /* Sort coordinates. */
#     define SORT(a,b) if (a > b) { double tmp = a; a = b; b = tmp; } else
      SORT(t0,t1);
      SORT(t2,t3);
      SORT(t0,t2);
      SORT(t1,t3);
      SORT(t1,t2);
#     undef SORT

      /* Normalization factor.  (qx,qy) are the coordinates of the source-voxel
         vector. */
      qx = x + sinBeta*sourceDistance;
      qy = y - cosBeta*sourceDistance;
      q = (fabs(qy) > fabs(qx) ? qx/qy : qy/qx);
      alpha = sqrt(1.0 + q*q)*voxelSize;

      /* Index of first/last pixel impacted by this voxel. */
      kFirst = (long)floor(t0 - 0.5*size);
      if (kFirst < 0) kFirst = 0;
      kLast = (long)ceil(t3 + 0.5*size);
      if (kLast >= detectorDimension) kLast = detectorDimension - 1;

      /* Integrate projected voxel profile (approximated by a trapezoid) on the
         impacted pixels.  GAMMA is the integral of the trapezoid over the
         pixel (divided by the pixel size). */
      sum = 0.0;
      for (k = kFirst; k <= kLast; ++k) {
	/* double s1 = k - 0.5*size; */
	/* double s2 = s1 + size; */
	/* double gamma = (integTriangle(max(s1, t0) - t0, */
        /*                               min(s2, t1) - t0, t1 - t0) + */
        /*                 integRectangle(max(s1, t1), min(s2, t2)) + */
        /*                 integTriangle(t3 - min(s2, t3), */
        /*                               t3 - max(s1, t2), t3 - t2)); */
        /* sum += alpha*gamma*pixel[k]; */
        double gamma;
        if (k>=t0 && k<t1) {
          gamma = (k-t0)/(t1-t0);
        } else if (k>=t1 && k<t2) {
          gamma = 1;     
        } else if (k>=t2 && k<t3) {
          gamma = (k-t3)/(t2-t3);
        } else {
          gamma = 0;
        }   

        sum += alpha*gamma*pixel[k];
      }
      voxel[i1 + i2*objectDimension1] += sum;
    }
  }
}

void
lf2d_sparser(const lf2d_parameters_t* p, /* parameters */
            double beta,                /* direction of source w.r.t. y-axis */
            double detectorOffset,      /* offset of detector w.r.t. to its
                                           geometrical center */
            const long i1,              /* object voxels */
            const long i2,              /* object voxels */
            double pixel[],             /* detector pixels */
            int clr)                    /* clear destination before
                                           integration */
{
  long objectDimension1 = p->objectDimension1;
  long objectDimension2 = p->objectDimension2;
  long k, kFirst, kLast, detectorDimension = p->detectorDimension;
  double voxelSize = p->voxelSize;
  double pixelSize = p->pixelSize;
  double pixelStep = p->pixelStep;
  double sourceDistance = p->sourceDistance;
  double detectorDistance = p->detectorDistance;

  double cosBeta = cos(beta);
  double sinBeta = sin(beta);
  double offset_y = cosBeta*(p->objectOffset1)+sinBeta*(p->objectOffset2);
  double offset_x = -sinBeta*(p->objectOffset1)+cosBeta*(p->objectOffset2);

  /* Position of the center of first voxel. */
  double x0 = offset_x - 0.5*(objectDimension1 - 1)*voxelSize;
  double y0 = offset_y - 0.5*(objectDimension2 - 1)*voxelSize;
  
  /* Parameters for computing the projection of a voxel on the detector.
     Projected coordinates are computed in pixel units. */
  double offset = detectorOffset/pixelStep + 0.5*(detectorDimension - 1);
  double scale = (sourceDistance + detectorDistance)/pixelStep;
  double size = pixelSize/pixelStep;
  
  /* Clear destination array (for integration of result). */
  if (clr) {
    memset(pixel, 0, sizeof(double)*detectorDimension);
  }
  
  /* Loop over all voxels. */
  
  double y = y0 + i2*voxelSize;  /* abscissa of voxel center */
  double y1 = y - 0.5*voxelSize;
  double y2 = y + 0.5*voxelSize;
  
  double x = x0 + i1*voxelSize;  /* ordinate of voxel center */
  double x1 = x - 0.5*voxelSize;
  double x2 = x + 0.5*voxelSize;
  double t0, t1, t2, t3;
  double alpha, q, qx, qy;
  
  /* Compute coordinates along detector axis of voxel at (x,y) in pixel
     units. */
#     define S(x,y) (scale*((y)*sinBeta + (x)*cosBeta)                  \
                     /(sourceDistance - ((y)*cosBeta - (x)*sinBeta)) + offset)
  t0 = S(x1,y1);
  t1 = S(x2,y1);
  t2 = S(x1,y2);
  t3 = S(x2,y2);
#     undef S
  
  /* Sort coordinates. */
#     define SORT(a,b) if (a > b) { double tmp = a; a = b; b = tmp; } else
  SORT(t0,t1);
  SORT(t2,t3);
  SORT(t0,t2);
  SORT(t1,t3);
  SORT(t1,t2);
#     undef SORT
  
  /* Normalization factor ALPHA is the trapezoid height times the pixel
     size (because the integral of the footprint is computed in pixel
     units).  (qx,qy) are the coordinates of the source-voxel vector. */
  qx = x + sourceDistance*sinBeta;
  qy = y - sourceDistance*cosBeta;
  q = (fabs(qy) > fabs(qx) ? qx/qy : qy/qx);
  alpha = (sqrt(1.0 + q*q)*(voxelSize*pixelSize));
  
  /* Index of first/last pixel impacted by this voxel. */
  kFirst = (long)floor(t0 - 0.5*size);
  if (kFirst < 0) kFirst = 0;
  kLast = (long)ceil(t3 + 0.5*size);
  if (kLast >= detectorDimension) kLast = detectorDimension - 1;
  
  /* Integrate projected voxel profile (approximated by a trapezoid) on the
     impacted pixels.  GAMMA is the integral of the trapezoid over the
     pixel (divided by the pixel size). */
  for (k = kFirst; k <= kLast; ++k) {
    double s1 = k - 0.5*size;
    double s2 = s1 + size;
    double gamma = (integTriangle(max(s1, t0) - t0,
                                  min(s2, t1) - t0, t1 - t0) +
                    integRectangle(max(s1, t1), min(s2, t2)) +
                    integTriangle(t3 - min(s2, t3),
                                  t3 - max(s1, t2), t3 - t2));
    pixel[k] += gamma*alpha;
  }

}

void
lf2d_line_integ_sparser(const lf2d_parameters_t* p, /* parameters */
                        double beta,                /* direction of source w.r.t. y-axis */
                        double detectorOffset,      /* offset of detector w.r.t. to its
                                                       geometrical center */
                        const long i1,              /* object voxels */
                        const long i2,              /* object voxels */
                        double pixel[],             /* detector pixels */
                        int clr)                    /* clear destination before
                                                       integration */
{
  long objectDimension1 = p->objectDimension1;
  long objectDimension2 = p->objectDimension2;
  long k, kFirst, kLast, detectorDimension = p->detectorDimension;
  double voxelSize = p->voxelSize;
  double pixelSize = p->pixelSize;
  double pixelStep = p->pixelStep;
  double sourceDistance = p->sourceDistance;
  double detectorDistance = p->detectorDistance;
  
  double cosBeta = cos(beta);
  double sinBeta = sin(beta);
  double offset_y = cosBeta*(p->objectOffset1)+sinBeta*(p->objectOffset2);
  double offset_x = -sinBeta*(p->objectOffset1)+cosBeta*(p->objectOffset2);

  /* Position of the center of first voxel. */
  double x0 = offset_x - 0.5*(objectDimension1 - 1)*voxelSize;
  double y0 = offset_y - 0.5*(objectDimension2 - 1)*voxelSize;
  
  /* Parameters for computing the projection of a voxel on the detector.
     Projected coordinates are computed in pixel units. */
  double offset = detectorOffset/pixelStep + 0.5*(detectorDimension - 1);
  double scale = (sourceDistance + detectorDistance)/pixelStep;
  double size = pixelSize/pixelStep;
  
  /* Clear destination array (for integration of result). */
  if (clr) {
    memset(pixel, 0, sizeof(double)*detectorDimension);
  }
  
  /* Loop over all voxels. */
  
  double y = y0 + i2*voxelSize;  /* abscissa of voxel center */
  double y1 = y - 0.5*voxelSize;
  double y2 = y + 0.5*voxelSize;
  
  double x = x0 + i1*voxelSize;  /* ordinate of voxel center */
  double x1 = x - 0.5*voxelSize;
  double x2 = x + 0.5*voxelSize;
  double t0, t1, t2, t3;
  double alpha, q, qx, qy;
  
  /* Compute coordinates along detector axis of voxel at (x,y) in pixel
     units. */
#     define S(x,y) (scale*((y)*sinBeta + (x)*cosBeta)                  \
                     /(sourceDistance - ((y)*cosBeta - (x)*sinBeta)) + offset)
  t0 = S(x1,y1);
  t1 = S(x2,y1);
  t2 = S(x1,y2);
  t3 = S(x2,y2);
#     undef S
  
  /* Sort coordinates. */
#     define SORT(a,b) if (a > b) { double tmp = a; a = b; b = tmp; } else
  SORT(t0,t1);
  SORT(t2,t3);
  SORT(t0,t2);
  SORT(t1,t3);
  SORT(t1,t2);
#     undef SORT
  
  /* Normalization factor ALPHA is the trapezoid height times the pixel
     size (because the integral of the footprint is computed in pixel
     units).  (qx,qy) are the coordinates of the source-voxel vector. */
  qx = x + sourceDistance*sinBeta;
  qy = y - sourceDistance*cosBeta;
  q = (fabs(qy) > fabs(qx) ? qx/qy : qy/qx);
  alpha = (sqrt(1.0 + q*q)*voxelSize);
  
  /* Index of first/last pixel impacted by this voxel. */
  kFirst = (long)floor(t0 - 0.5*size);
  if (kFirst < 0) kFirst = 0;
  kLast = (long)ceil(t3 + 0.5*size);
  if (kLast >= detectorDimension) kLast = detectorDimension - 1;
  
  /* Integrate projected voxel profile (approximated by a trapezoid) on the
     impacted pixels.  GAMMA is the integral of the trapezoid over the
     pixel (divided by the pixel size). */
  for (k = kFirst; k <= kLast; ++k) {
    if (k>=t0 && k<t1) {
      pixel[k] += alpha*(k-t0)/(t1-t0);
    } else if (k>=t1 && k<t2) {
      pixel[k] += alpha;     
    } else if (k>=t2 && k<t3) {
      pixel[k] += alpha*(k-t3)/(t2-t3);
    } else {
      pixel[k] += 0;
    }   
    /* double s1 = k - 0.5*size; */
    /* double s2 = s1 + size; */
    /* double gamma = (integTriangle(max(s1, t0) - t0, */
    /*                               min(s2, t1) - t0, t1 - t0) + */
    /*                 integRectangle(max(s1, t1), min(s2, t2)) + */
    /*                 integTriangle(t3 - min(s2, t3), */
    /*                               t3 - max(s1, t2), t3 - t2)); */
    /* pixel[k] += gamma*alpha; */
  }

}

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * fill-column: 79
 * coding: utf-8
 * End:
 */
