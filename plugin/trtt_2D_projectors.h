/*
 * TR-TT : iTeraTive Reconstruction for Tomography
 *
 * trtt_2D_projectors.h --
 *
 * TRTT 2D projectors.
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

/* PREPROCESSOR =============================================================== */

#ifndef TRTT_2D_PROJECTORS_H
#define TRTT_2D_PROJECTORS_H 1

/* PREPROCESSOR INCLUSIONS ==================================================== */

#include <stdio.h>
#include "trtt_common.h"
#include "trtt_xform3d.h"

/* PROTOYPES ================================================================== */

TRTT_BEGIN_C_DECLS

/*
 * Title: TR-TT Radon operators for tomography
 *
 *   Version 1.0
 *
 *   Copyright (C) 2011, the MiTiV project.
 *
 *   This package contains the operators performing the algebraic Radon
 *   transform at several dimensions and geometries.
 *
 *   The definitions are available by: : #include <trtt_operators.h>
 */


extern int trtt_PB2D(double *out, double *in, long nx, long ny, double s_scl, int s_deg, double *footprint, double *cum_footprint, long size_footprint, double step_footprint, double *coord_footprint, long nv, double v_scl, double *xos, double *xsd, int job);

extern int trtt_FB2D(double *out, double *in, long nx, long ny, double s_scl, int s_deg, double *footprint, double *cum_footprint, long size_footprint, double step_footprint, double *coord_footprint, long nv, double v_scl, double *xos, double *xsd, int job);

extern int trtt_PB2D_cubic(double *out, double *in, long nx, long ny, double s_scl, int s_deg, long nv, double v_scl, double *xos, double *xsd, int job);

extern int trtt_FB2D_cubic(double *out, double *in, long nx, long ny, double s_scl, int s_deg, long nv, double v_scl, double *xos, double *xsd, int job);

extern int trtt_PB2D_sparser(double *out, long i, long j, long nx, long ny, double s_scl, int s_deg, double *footprint, double *cum_footprint, long size_footprint, double step_footprint, double *coord_footprint, long nv, double v_scl, double *xos, double *xsd);

extern int trtt_FB2D_sparser(double *out, long i, long j, long nx, long ny, double s_scl, int s_deg, double *footprint, double *cum_footprint, long size_footprint, double step_footprint, double *coord_footprint, long nv, double v_scl, double *xos, double *xsd);

TRTT_END_C_DECLS

#endif /* TRTT_2D_PROJECTORS_H */

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
