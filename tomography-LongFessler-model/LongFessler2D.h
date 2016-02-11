/*
 * LongFessler2D.h -
 *
 * Definitions for Long & Fessler model for 2D tomography.
 */

#ifndef _LONGFESSLER2D_H
#define _LONGFESSLER2D_H 1

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _lf2d_parameters lf2d_parameters_t;

struct _lf2d_parameters {
  double pixelSize;                        /* detector pixel size */
  double pixelStep;                        /* detector pixel step (>=
                                              pixelSize) */
  double voxelSize;                        /* voxel size */
  double sourceDistance;                   /* distance origin-source */
  double detectorDistance;                 /* distance origin-detector */
  double objectOffset1, objectOffset2;     /* offsets of object center
                                              w.r.t. to its geometrical
                                              center */
  long objectDimension1, objectDimension2; /* dimensions of object */
  long detectorDimension;                  /* size of detector */
};

extern void
lf2d_direct(const lf2d_parameters_t* p, /* parameters */
            double beta,                /* direction of source w.r.t. y-axis */
            double detectorOffset,      /* offset of detector w.r.t. to its
                                           geometrical center */
            const double voxel[],       /* object voxels */
            double pixel[],             /* detector pixels */
            int clr);                   /* clear destination before
                                           integration */

extern void
lf2d_line_integ_direct(const lf2d_parameters_t* p, /* parameters */
                       double beta,                /* direction of source w.r.t. y-axis */
                       double detectorOffset,      /* offset of detector w.r.t. to its
                                                      geometrical center */
                       const double voxel[],       /* object voxels */
                       double pixel[],             /* detector pixels */
                       int clr);                   /* clear destination before
                                                      integration */
  
 extern void
 lf2d_adjoint(const lf2d_parameters_t* p, /* parameters */
              double beta,                /* direction of source
                                             w.r.t. y-axis */
              double detectorOffset,      /* offset of detector w.r.t. to its
                                             geometrical center */
              const double pixel[],       /* source (detector pixels) */
              double voxel[],             /* object voxels */
              int clr);                   /* clear destination before
                                             integration */
  
 extern void
 lf2d_line_integ_adjoint(const lf2d_parameters_t* p, /* parameters */
                         double beta,                /* direction of source
                                                        w.r.t. y-axis */
                         double detectorOffset,      /* offset of detector w.r.t. to its
                                                        geometrical center */
                         const double pixel[],       /* source (detector pixels) */
                         double voxel[],             /* object voxels */
                         int clr);                   /* clear destination before
                                                        integration */

#ifdef __cplusplus
}
#endif

#endif /* _LONGFESSLER2D_H */

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
