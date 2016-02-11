/*
 * TRTT : iTeraTive Reconstruction for Tomography
 *
 * Ytrtt_2D_projectors.i --
 *
 * TRTT 2D Radon projectors.
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
 * the user in light of its specific status of free software, that may mean that
 * it is complicated to manipulate, and that also therefore means that it is
 * reserved for developers and experienced professionals having in-depth
 * computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling
 * the security of their systems and/or data to be ensured and, more generally,
 * to use and operate it in the same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 *
 *-----------------------------------------------------------------------------
 *
 * $Id$
 * $Log$
 */

/* REQUIREMENTS ============================================================== */

extern trtt_2D_projectors_info;
/* DOCUMENT trtt_2D_projectors_info

   Ytrtt_2D_projectors.i is a module of TRTT code for creating tomographic
   projectors, given TOMOGRAPHIC and DATA OBJECTS. These projectors are used as
   h_evaluator which carry references to the objects containing calibration
   parameters to perform the projection, following the Radon transform
   principle, of an image of voxels values with compatible dimensions. One
   projector is related on a single acquisition at a given instant, that is to
   say a given position of the Source associated with a given position of the
   Detector plane and the Object of interest. Two models are available for the
   projectors:

   - The Separable B-spline model : TRTT code is mainly developped to use this
     model. It consists in taking separable B-splines of desired degree as basis
     functions for the representation of the function associated with the object
     of interest. The projection of such a function on a detector plane
     (footprint) is approximated by a separable b-spline suitably scaled with a
     magnification and distorsion factor (depending on the geometry of
     projection and the incidence of the X rays on the detector plane). This
     approximation is justified by the quasi-sphericity of a b-spline of
     sufficiently high degree, yielding separability whatever the axis and
     isotropic behaviour. This projector has been implemented in C language and
     interfaced with Yorick, so as to be efficient for iterative reconstruction
     algorithms.
     
   - The Distance Driven model, developped by DeMan & Basu (2004). It is based
     on a cubic voxel representation of the object of interest. The projection
     of a cubic voxel on the detector plane is performed by determining the
     projection of the axial section of the cubic voxel, and approximating the
     quadrilateral form of the footprint by a rectangle. This model is also more
     precise than other standard models such Ray Driven, but is highly
     anisotropic and will suffer from high reconstruction errors when the number
     of projections is limited.

   These projectors are the static part of the tomographic data
   modelization. For applications of Dynamic tomography, the object of interest
   becomes a sequence of 2D images, and data corresponds to a given instant of
   acquisition. Therefore we need an operator which determines the temporal
   slice that we are projecting. This operator is a time interpolator of the
   2D+t image, using b-spline interpolation (SPL library) for a good
   modelization of the continuity of the sequence in the "temporal"
   direction. It is combined with a static projector to make a "dynamic"
   projector. If we want to perform static tomography, we do not need such
   projector.

   This module allows to create sets of projectors, given a set of DATA objects
   corresponding for example to a rotation of the source around the object of
   interest (standard framework for tomography).
   
   _
    |
    |_ COPS
       |
       |_ C01        => PROJECTOR
       |
       |_ C02        => PROJECTOR
       |
       |_ Ck         => PROJECTOR
       |
       :
       :
       |_ Cndata     => PROJECTOR
       |
       |_ Ck_list = ["C01","C02",...,"Ck",...,"Cndata"]

   Such a set of projectors is used to perform iterative reconstruction from a
   set of data sinograms, using the optimization tools developped at CRAL by
   Éric Thiébaut (OptimPack1 => SPL library).

   List of functions of the module (see the help of these functions):

   - trtt_single_2D_projector()
   - trtt_whole_2D_projector()
   - trtt_single_2D_time_interpolator()
   - trtt_single_distance_driven_2D_projector()
   - trtt_whole_distance_driven_2D_projector()
   - trtt_single_2D_combined_projector()
   - trtt_whole_2D_combined_projector()

   SEE ALSO trtt_tomdata_info, trtt_tomobj_info.

 */

/* GLOBALS =================================================================== */

/*|---------------|*/
/*| 2D PROJECTORS |*/
/*|---------------|*/

func trtt_single_Long_Fessler_line_integ_2D_projector(X, Yk, mode, status=)
/* DOCUMENT trtt_single_2D_projector, X, Yk, mode, status=

   Create a 2D projector with the Separable B-spline model, according to the
   geometry of projection MODE (TRTT_PARALLEL_BEAM or TRTT_FAN(CONE)_BEAM). X is
   the TOMOGRAPHIC OBJECT and Yk is the DATA object. The projector consists in a
   h_evaluator carrying references to X and Yk, associated with the adapted
   projection function. A STATUS can be given. If not it is created in the
   function and if an error occurs, it will be directly displayed.

   SEE ALSO: trtt_2D_projectors_info, _trtt_single_2D_projector, h_evaluator.
 */    
{
    if (!is_void(status)) {
        flag_disp_status = TRTT_FALSE;
        /* A status variable is given : in case of error, status is updated to
         * be passed to upper callers */
        if(trtt_error_query(status)) {
            return;
        }
    } else {
        flag_disp_status = TRTT_TRUE;
        /* We initialize a status variable for internal functions */
        status = trtt_error_init();
    }

    /* X must be a tomobj */
    if (!is_tomobj(X)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMOBJ, disp=flag_disp_status;
        return;
    }

    /* Yk must be a tomdata */
    if (!is_tomdata(Yk)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMDATA, disp=flag_disp_status;
        return;
    }

    /* Extract needed parameters */
    /* - X */
    obj_xi = X.xi;
    nt = X.nt;
    /* - Yk */
    data_xi = Yk.xi;
    xotilt = Yk.xotilt;
    xos = Yk.xos;
    xsd = Yk.xsd;
    t_index = Yk.t_index;
    
    if (is_void(Yk.xoff)) {
        xoff=0.0;
    } else {
        xoff=Yk.xoff;
    }

    if (is_void(Yk.yoff)) {
        yoff=0.0;
    } else {
        yoff=Yk.yoff;
    }

    if (is_void(Yk.voff)) {
        voff=0.0;
    } else {
        voff=Yk.voff;
    }
    
    /* xform combinations */
    xos_plus = trtt_xform3d_combine(xotilt, obj_xi);
    trtt_xform3d_combine, xos_plus, xos, xos_plus;
    
    xos_coefs = xos_plus.a;
    data_xi_inv = trtt_xform3d_invert(data_xi);
    xsd_plus = trtt_xform3d_combine(data_xi_inv, xsd);
    xsd_coefs = xsd_plus.a;

    /* Check t_index */
    if (!(t_index>=0 & t_index<=nt)) {
        trtt_error_set, status, TRTT_ERR_OUT_OF_RANGE_INDEX, msg="Incompatible TOMDATA 't_index' with dimension 'nt' of the TOMOBJ.", disp=flag_disp_status;
        return;
    }

    lf_param=lf2d_parameters(pixelSize=Yk.v_scl,
                             pixelStep=Yk.v_scl,
                             voxelSize=X.s_scl,
                             sourceDistance=Yk.Rsc,
                             detectorDistance=Yk.Rsd-Yk.Rsc,
                             objectOffset1=xoff,
                             objectOffset2=yoff,
                             objectDimension1=X.nx,
                             objectDimension2=X.ny,
                             detectorDimension=Yk.nv);
    
    /* Create the evaluator */
    p = h_new(X=X,
              Yk=Yk,
              xos_coefs = xos_coefs,
              xsd_coefs = xsd_coefs,
              offset=voff,
              beta=Yk.theta-(pi/2),
              lf_param = lf_param,
              sparser = symlink_to_name("lf2d_line_integ_sparser"));

    h_evaluator, p, "_trtt_single_Long_Fessler_line_integ_2D_projector";

    return p;
}
L_line_integ_single_2D = trtt_single_Long_Fessler_line_integ_2D_projector;

func _trtt_single_Long_Fessler_line_integ_2D_projector(p, x, job)
/* DOCUMENT _trtt_single_2D_projector, p, x, job

   Evaluator (direct and tranpose) of a 2D projector using the Separable
   B-spline model. According to the geometry of projection, a call to the suited
   C-implemented projector is done.

   SEE ALSO: trtt_single_2D_projector, h_evaluator.
 */ 
{
    lf_param=p.lf_param;
    nx=p.lf_param.objectDimension1;
    ny=p.lf_param.objectDimension2;
    nv=p.lf_param.detectorDimension;
    beta=p.beta;
    offset=p.offset;

    if (!job) {
        job=0n;
        out = array(double, nv);
    } else if (job==1) {
        out = array(double, nx, ny);
    } else {
        trtt_error_display_msg, code=TRTT_ERR_UNSUPPORTED_JOB;
        return;
    }

    if (!job) {
        lf2d_line_integ_direct, lf_param, beta, offset, x, out, 1n;
    } else if (job==1) {
        lf2d_line_integ_adjoint, lf_param, beta, offset, x, out, 1n;
    } else {
        trtt_error_display_msg, code=TRTT_ERR_UNSUPPORTED_JOB;
        return;
    }

    return out; 
}

func _trtt_single_Long_Fessler_line_integ_2D_sparser(p, i, j)
/* DOCUMENT _trtt_single_Long_Fessler_2D_sparser, p, x, job

   Evaluator (direct and tranpose) of a 2D projector using the Separable
   B-spline model. According to the geometry of projection, a call to the suited
   C-implemented projector is done.

   SEE ALSO: trtt_single_2D_projector, h_evaluator.
 */ 
{
    lf_param=p.lf_param;
    nx=p.lf_param.objectDimension1;
    ny=p.lf_param.objectDimension2;
    nv=p.lf_param.detectorDimension;
    beta=p.beta;
    offset=p.offset;
    __sparser = p.sparser;

    out = array(double, nv);
    __sparser, lf_param, beta, offset, i, j, out, 1n;
    
    return out; 
}

func trtt_whole_Long_Fessler_line_integ_2D_projector(X, Y, mode, status=)
/* DOCUMENT trtt_whole_2D_projector, X, Y, mode, status=

   Create a set of 2D projectors based on the Separable B-spline model, from a
   TOMOGRAPHIC OBJECT X and a set of DATA objects Y, according to the geometry
   of projection MODE. A STATUS can be given. If not it is created in the
   function and if an error occurs, it will be directly displayed.

   SEE ALSO: trtt_2D_projectors_info, trtt_single_2D_projector.
 */
{
    if (!is_void(status)) {
        flag_disp_status = TRTT_FALSE;
        /* A status variable is given : in case of error, status is updated to
         * be passed to upper callers */
        if(trtt_error_query(status)) {
            return;
        }
    } else {
        flag_disp_status = TRTT_TRUE;
        /* We initialize a status variable for internal functions */
        status = trtt_error_init();
    }

    /* X must be a tomobj */
    if (!is_tomobj(X)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMOBJ, disp=flag_disp_status;
        return;
    }

    // FIXME: Check if Y is a hash table containing TOMDATA objects
    if (!is_tomdata_set(Y)) {
        trtt_error_set, status, TRTT_ERR_NOT_TOMDATA_SET, disp=flag_disp_status;
        return;
    }

    Yk_list = Y.Yk_list;
    ndata = numberof(Yk_list);

    /* Create then hash table which contains the whole projectors */
    Cops = h_new();

    for (k=1; k<=ndata; ++k) {
        Yk = h_get(Y, Yk_list(k));

        Ck = L_line_integ_single_2D(X, Yk, mode, status=status);
        if(is_void(Ck)) {
            if (flag_disp_status) trtt_error_display_msg, status=status;
            return;
        }

        h_set, Cops, trtt_get_key_name(k, "C"), Ck;
    }

    Ck_list = h_keys(Cops);  
    h_set, Cops, Ck_list = Ck_list(sort(Ck_list));
    
    return Cops;
}
L_line_integ_whole_2D = trtt_whole_Long_Fessler_line_integ_2D_projector;
