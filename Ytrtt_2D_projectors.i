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

/* INTERNAL HIGH LEVEL FUNCTIONS ============================================= */

func _trtt_2D_projection_O_to_D (xn, yn, xos, xsd, mode)
/* DOCUMENT _trtt_2D_projection_O_to_D, xn, yn, xos, xsd, mode

   Project the 2D coordinates array (XN,YN) on a detector plane according to the
   spatial transforms coefficients XOS (Object to Source) and XSD (Source to
   Detector) and the geometry of projection MODE.

 */    
{
    /* O to S transform */
    xs = xos(1)*xn + xos(2)*yn(-,) + xos(4);
    ys = xos(5)*xn + xos(6)*yn(-,) + xos(8);

    /* S to D projection */
    if (anyof(xs==0.0) & mode == TRTT_FAN_BEAM)
        error, "IT WAS SUPPOSED TO BE HIGHLY IMPOSSIBLE !!!\n";

    a1 = xsd(1);
    if (a1 == 0)
        error, "IT WAS SUPPOSED TO BE HIGHLY IMPOSSIBLE !!!\n";
    
    a1inv = 1.0/a1;
    a2 = xsd(2);
    a4 = xsd(4);
    a5 = xsd(5);
    a6 = xsd(6);
    a8 = xsd(8);

    b1 = a6 + a1inv*a5*a2;
    b3 = a8 + a1inv*a5*a4;

    if (mode == TRTT_PARALLEL_BEAM) {
        return (b1*ys+b3);
    } else {
        tan_gamma = ys/xs;
        xeq = -(1.0/(a1+a2*tan_gamma))*a4; 
        return (b1*tan_gamma*xeq + b3);
    }
}


func _trtt_2D_distorsion_O_to_D (xn, yn, xos, xsd)
/* DOCUMENT _trtt_2D_magnification_O_to_D, xn, yn, xos, xsd

   Calculate the magnification factor on the detector of each voxel position
   (XN,YN) according to the spatial transforms coefficients XOS (Object to
   Source) and XSD (Source to Detector) in fan beam geometry. 
 */    
{
    /* O to S transform */
    xs = xos(1)*xn + xos(2)*yn(-,) + xos(4);
    ys = xos(5)*xn + xos(6)*yn(-,) + xos(8);

    /* S to D projection */
    if (anyof(xs==0.0))
        error, "IT WAS SUPPOSED TO BE HIGHLY IMPOSSIBLE !!!\n";

    return sqrt(1+(ys/xs)^2.0);
}

func _trtt_2D_magnification_O_to_D (xn, yn, xos, xsd)
/* DOCUMENT _trtt_2D_magnification_O_to_D, xn, yn, xos, xsd

   Calculate the magnification factor on the detector of each voxel position
   (XN,YN) according to the spatial transforms coefficients XOS (Object to
   Source) and XSD (Source to Detector) in fan beam geometry. 
 */    
{
    /* O to S transform */
    xs = xos(1)*xn + xos(2)*yn(-,) + xos(4);

    /* S to D projection */
    if (anyof(xs==0.0))
        error, "IT WAS SUPPOSED TO BE HIGHLY IMPOSSIBLE !!!\n";

    a4 = xsd(4);

    return abs(a4/xs);
}

func _trtt_2D_sparse_coefs(Cops, dyn=)
/* DOCUMENT _trtt_2D_sparse_coefs, Cops, dyn=

   Calculate the sparse projector matrix coefficients from each projector Ck of
   the set COPS. Sparse matrix components are added to each Ck in order not have
   to re-calculate it everytime.
 */
{
    type = Cops.type;
    Ck_list = Cops.Ck_list;
    ndata = numberof(Ck_list);
    if (!dyn) {
        X = h_get(Cops,Ck_list(1)).X;
    } else {
        X = h_get(Cops,Ck_list(1)).Ak.X;
    }
    nx = X.nx;
    ny = X.ny;

    if (!Cops.exist_sparse_coefs) {
        /* Get the sparse matrix parameters for each projector */
        for (k=1; k<=ndata; ++k) {
            /* Get the current projector */
            if (!dyn) {
                Ck = h_get(Cops,Ck_list(k));
            } else {
                Ck = h_get(Cops,Ck_list(k)).Ak;
            }
            /* Get the current detector size */
            nv = Ck.Yk.nv;

            // MCk_coefs = [];
            // MCk_row_indices = [];
            // MCk_col_indices = [];
            MCk_row_dimlist = [1,nv];
            MCk_col_dimlist = [1,nx*ny];

            f1=create("temp_trtt_fichier1");
            f2=create("temp_trtt_fichier2");
            f3=create("temp_trtt_fichier3");
            f1=open("temp_trtt_fichier1","r+");
            f2=open("temp_trtt_fichier2","r+");
            f3=open("temp_trtt_fichier3","r+");
            nsparsetot=0;
            
            for (j=0; j<ny; ++j) {
                for (i=0; i<nx; ++i) {
                    /* Calculate the matrix coefficients */
                    if (type==TRTT_SPLINE_DRIVEN) {
                        matrix_coefs = _trtt_single_2D_projector_sparser(Ck, i, j)(*);
                    } else if (type==TRTT_LONG_FESSLER) {
                        matrix_coefs = _trtt_single_Long_Fessler_2D_sparser(Ck, i, j)(*);
                    } else if  (type==TRTT_LONG_FESSLER_LINE_INTEG) {
                        matrix_coefs = _trtt_single_Long_Fessler_line_integ_2D_sparser(Ck, i, j)(*);
                    }
                    
                    /* Get the sparse indices */
                    sparse_indices = where(matrix_coefs != 0.0);
                    if(is_array(sparse_indices)) {
                        nsparse = numberof(sparse_indices);
                        nsparsetot+=nsparse;
                        
                        write, f1, format="%f\n", matrix_coefs(sparse_indices);
                        write, f2, format="%i\n", sparse_indices;
                        write, f3, format="%i\n", array(long(j*nx+i+1), nsparse);
                        /* Grow sparse matrix coefficients */
                        // grow, MCk_coefs, matrix_coefs(sparse_indices);
                        /* Grow sparse matrix column indices */
                        // grow, MCk_row_indices, sparse_indices;
                        /* Grow sparse matrix row indices */
                        // grow, MCk_col_indices, array(long(j*nx+i+1), nsparse);
                    }
                }
            }
            
            MCk_coefs = array(double,nsparsetot);
            MCk_row_indices = array(long,nsparsetot);
            MCk_col_indices = array(long,nsparsetot);
            read, f1, format="%f\n", MCk_coefs;
            read, f2, format="%i\n", MCk_row_indices;
            read, f3, format="%i\n", MCk_col_indices;
            close, f1;
            close, f2;
            close, f3;
            system, "rm -fr temp_trtt_fichier1 temp_trtt_fichier2 temp_trtt_fichier3";
            
            /* Store sparse matrix parameters */
            h_set, Ck, MCk_coefs = MCk_coefs;
            h_set, Ck, MCk_row_dimlist = MCk_row_dimlist;
            h_set, Ck, MCk_row_indices = MCk_row_indices;
            h_set, Ck, MCk_col_dimlist = MCk_col_dimlist;
            h_set, Ck, MCk_col_indices = MCk_col_indices;
        }
        h_set, Cops, exist_sparse_coefs = 1n;
    }

    return TRTT_SUCCESS;
}

func _trtt_2D_sparse_matrix(Cops, dyn=)
/* DOCUMENT _trtt_2D_sparse_matrix, Cops, dyn=

   Calculate the sparse projector matrix coefficients from the set of projectors
   COPS and return the sparse matrix. Sparse matrix components are added to COPS
   in order not have to re-calculate it everytime.
 */
{
    Ck_list = Cops.Ck_list;
    ndata = numberof(Ck_list);
    X = h_get(Cops,Ck_list(1)).X;

    _trtt_2D_sparse_coefs, Cops, dyn=dyn;
    
    /* Creating the whole sparse matrix */
    nx = X.nx;
    ny = X.ny;
    nt = X.nt;
    n_col = nx*ny*nt;
    n_row = 0.0;
    /* Initialize sparse matrix parameters */
    MCops_coefs = [];
    MCops_row_indices = [];
    MCops_col_indices = [];
    /* Initialize whole sinogram data */
    y = [];
    n_row_old = long(0);
    for (k=1; k<=ndata; ++k) {
        /* Get the current projector */
        if (!dyn) {
            Ck = h_get(Cops,Ck_list(k));
        } else {
            Ck = h_get(Cops,Ck_list(k)).Ak;
        }
        /* Get and store the current projection data */
        Yk = Ck.Yk;
        /* Get the current detector size */
        nv = Yk.nv;
        n_row += nv;
        grow, y, Yk.pixels;
        /* Get sparse matrix parameters */
        grow, MCops_coefs, Ck.MCk_coefs;
        grow, MCops_row_indices, n_row_old+Ck.MCk_row_indices;
        grow, MCops_col_indices, Ck.MCk_col_indices;
        n_row_old = n_row;
    }
    /* Last sparse matrix parameters */
    MCops_row_dimlist = [1,long(n_row)];
    MCops_col_dimlist = [1,long(n_col)];

    /* Create the sparse matrix */
    sparse_MCops = sparse_matrix(MCops_coefs,
                                 MCops_row_dimlist, MCops_row_indices,
                                 MCops_col_dimlist, MCops_col_indices);

    h_set, Cops, matrix = 1n;
    
    return sparse_MCops;
}

/* USER FUNCTIONS ============================================================ */

/*|---------------|*/
/*| 2D PROJECTORS |*/
/*|---------------|*/

func trtt_single_2D_projector(X, Yk, mode, status=)
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

    /* Create the evaluator */
    p = h_new(X=X,
              Yk=Yk,
              xos_coefs = xos_coefs,
              xsd_coefs = xsd_coefs);

    if (mode == TRTT_PARALLEL_BEAM) {
        h_set, p, projector = symlink_to_name("trtt_PB2D");
        h_set, p, sparser = symlink_to_name("trtt_PB2D_sparser");
    }
    else if (mode == TRTT_FAN_BEAM || mode == TRTT_CONE_BEAM) {
        h_set, p, projector = symlink_to_name("trtt_FB2D");
        h_set, p, sparser = symlink_to_name("trtt_FB2D_sparser");    
    }
    else {
        trtt_error_set, status, TRTT_ERR_BAD_GEOMETRY, disp=flag_disp_status;
        return;
    }

    h_evaluator, p, "_trtt_single_2D_projector";

    return p;
}
A_single_2D = trtt_single_2D_projector;

func _trtt_single_2D_projector(p, x, job)
/* DOCUMENT _trtt_single_2D_projector, p, x, job

   Evaluator (direct and tranpose) of a 2D projector using the Separable
   B-spline model. According to the geometry of projection, a call to the suited
   C-implemented projector is done.

   SEE ALSO: trtt_single_2D_projector, h_evaluator.
 */ 
{
    nx = p.X.nx;
    ny = p.X.ny;
    s_scl = p.X.s_scl;
    s_deg = p.X.s_deg;
    local footprint; eq_nocopy, footprint, p.X.footprint;
    local cum_footprint; eq_nocopy, cum_footprint, p.X.cum_footprint;
    size_footprint = p.X.size_footprint;
    step_footprint = p.X.step_footprint;
    local coord_footprint; eq_nocopy, coord_footprint, p.X.coord_footprint;
    nv = p.Yk.nv;
    v_scl = p.Yk.v_scl;
    local xos_coefs; eq_nocopy, xos_coefs, p.xos_coefs;
    local xsd_coefs; eq_nocopy, xsd_coefs, p.xsd_coefs;
    __projector = p.projector;

    if (!job) {
        job=0n;
        out = array(double, nv);
    } else if (job==1) {
        out = array(double, nx, ny);
    } else {
        trtt_error_display_msg, code=TRTT_ERR_UNSUPPORTED_JOB;
        return;
    }

    if(__projector(out, x, nx, ny, s_scl, s_deg,
                   footprint, cum_footprint,
                   size_footprint, step_footprint, coord_footprint,
                   nv, v_scl,
                   xos_coefs, xsd_coefs, job) == TRTT_FAILURE) {
        trtt_error_display_msg, code=TRTT_ERR_PROJECTOR_FAILURE;
        return;
    } else {
        return out;
    }             
}

func _trtt_single_2D_projector_sparser(Ck, i, j)
/* DOCUMENT _trtt_single_2D_projector_sparser, Ck, i, j

   Evaluator (direct and tranpose) of a 2D projector using the Separable
   B-spline model. According to the geometry of projection, a call to the suited
   C-implemented projector is done.

   SEE ALSO: trtt_single_2D_projector, h_evaluator.
 */ 
{
    nx = Ck.X.nx;
    ny = Ck.X.ny;
    s_scl = Ck.X.s_scl;
    s_deg = Ck.X.s_deg;
    local footprint; eq_nocopy, footprint, Ck.X.footprint;
    local cum_footprint; eq_nocopy, cum_footprint, Ck.X.cum_footprint;
    size_footprint = Ck.X.size_footprint;
    step_footprint = Ck.X.step_footprint;
    local coord_footprint; eq_nocopy, coord_footprint, Ck.X.coord_footprint;
    nv = Ck.Yk.nv;
    v_scl = Ck.Yk.v_scl;
    local xos_coefs; eq_nocopy, xos_coefs, Ck.xos_coefs;
    local xsd_coefs; eq_nocopy, xsd_coefs, Ck.xsd_coefs;
    __sparser = Ck.sparser;

    out = array(double, nv);

    if(__sparser(out, i, j, nx, ny, s_scl, s_deg,
                 footprint, cum_footprint,
                 size_footprint, step_footprint, coord_footprint,
                 nv, v_scl,
                 xos_coefs, xsd_coefs) == TRTT_FAILURE) {
        trtt_error_display_msg, code=TRTT_ERR_PROJECTOR_FAILURE;
        return;
    } else {
        return out;
    }             
}

func trtt_whole_2D_projector(X, Y, mode, status=)
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

        Ck = A_single_2D(X, Yk, mode, status=status);
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
A_whole_2D = trtt_whole_2D_projector;

func _trtt_mpy_2D_projector(x, Cops)
/* DOCUMENT _trtt_mpy_2D_projector, x, Cops

   MPY implementation of the 2D projection function of image X. The calculation
   of the projections is sequentially distributed in several instances of
   MP_SIZE operations, each one applying a given projector Ck of COPS, until
   every one has been applied. Then the resulting sinogram is returned.

   SEE ALSO: _trtt_single_2D_projector.
 */
{
    Ck_list = Cops.Ck_list;
    ndata = numberof(Ck_list);

    /* Extraction des paramètres communs */
    C1 = h_get(Cops,Ck_list(1));
    X = C1.X;
    nx = X.nx; ny = X.ny;
    s_scl = X.s_scl;
    s_deg = X.s_deg;
    nv = C1.Yk.nv;
    v_scl = C1.Yk.v_scl;

    /* Extraction de la fonction de projection */
    funcproj = name_of_symlink(C1.projector);

    /* Envoi des paramètres communs */
    mp_exec, "mp_handout, x, nx, ny, s_scl, s_deg, nv, v_scl, funcproj;";
    /* Récupération par tous les rangs de la fonction de projection */
    mp_exec, "__projector = symlink_to_name(funcproj);";

    /* Calcul du nombre d'occurrences */
    Nmult = ndata/mp_size;
    Nrest = ndata%mp_size;

    /* Création du tableau de "récupération" des données */
    out_final = array(double, nv, ndata);

    /* Calcul des occurrences */
    for (i=0; i<=Nmult; ++i) {

        if (i<Nmult) nproc = mp_size;
        else if (Nrest > 0) nproc = Nrest;
        else break;
        
        k_0 = i*mp_size + 1;
        Ck_0 = h_get(Cops,Ck_list(k_0));

        for (k=1; k<nproc; ++k) {
            Ck = h_get(Cops,Ck_list(k_0+k));
            xos_coefs = Ck.xos_coefs;
            xsd_coefs = Ck.xsd_coefs;
            nv = Ck.Yk.nv;
            v_scl = Ck.Yk.v_scl;

            /* Envoi du numéro de rang concerné + */
            mp_exec, "mp_handout, k, nproc";
            /* Envoi des paramètres propres */
            mp_exec, "if (!mp_rank) {mp_send, k, vpack(xos_coefs);} else if (mp_rank==k) {vunpack, mp_recv(0), xos_coefs;} ";
            mp_exec, "if (!mp_rank) {mp_send, k, vpack(xsd_coefs);} else if (mp_rank==k) {vunpack, mp_recv(0), xsd_coefs;} ";
        }

        /* Paramètres propres du rang 0 */
        xos_coefs = Ck_0.xos_coefs;
        xsd_coefs = Ck_0.xsd_coefs;

        /* Commande de projection */
        mp_exec, "out = array(double, nv); __projector, out, x, nx, ny, s_scl, s_deg, nv, v_scl, xos_coefs, xsd_coefs, 0;";
        mp_exec, "out_rank = array(double, nv, nproc); out_rank(,mp_rank)=out;";
        mp_exec, "sino = mp_handin(out_rank);";
        out_final(,k_0:k_0+nproc-1) = sino;
    }
    
    return out_final;
}

/*|-----------------------|*/
/*| 2D TIME INTERPOLATORS |*/
/*|-----------------------|*/

func trtt_single_2D_time_interpolator(X, Yk, status=)
/* DOCUMENT trtt_single_2D_time_interpolator, X, Yk, status=

   Create a 2D time interpolator with b-spline interpolation framework (SPL
   library). X is the TOMOGRAPHIC OBJECT and Yk is the DATA object. The
   interpolator consists in a h_evaluator carrying the time interpolator WT,
   created with the parameters contained in X and Yk, associated with a function
   which uses this interpolator on an image X. A STATUS can be given. If not it
   is created in the function and if an error occurs, it will be directly
   displayed.

   SEE ALSO: trtt_2D_projectors_info, _trtt_single_2D_time_interpolator, h_evaluator.
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
    nt = X.nt;
    nx = X.nx;
    ny = X.ny;
    t_scl = X.t_scl;
    t_deg = X.t_deg;
    /* - Yk */
    t_index = Yk.t_index;

    /* Check t_index */
    if (!(t_index>=0 & t_index<=nt)) {
        trtt_error_set, status, TRTT_ERR_OUT_OF_RANGE_INDEX, msg="Incompatible TOMDATA 't_index' with dimension 'nt' of the TOMOBJ.", disp=flag_disp_status;
        return;
    }

    /* Create the interpolator */
    wt = spl_interp_coefs_new_spline(t_index, nt, t_deg);

    p = h_new(wt = wt,
              nx = nx,
              ny = ny,
              interpolator = symlink_to_name("_trtt_2D_time_interpolator"));

    h_evaluator, p, "_trtt_single_2D_time_interpolator";

    return p;
}
B_single_2D = trtt_single_2D_time_interpolator;

func _trtt_single_2D_time_interpolator(p, x, job)
/* DOCUMENT _trtt_single_2D_time_interpolator, p, x, job

   Evaluator (direct and tranpose) of a 2D time interpolator based on b-spline
   interpolation.

   SEE ALSO: trtt_single_2D_time_interpolator, h_evaluator.
 */ 
{
    wt = p.wt;
    __interpolator = p.interpolator;
    nx = p.nx;
    ny = p.ny;
    m = wt.m;
    n = wt.n; // = nt
    p = wt.p;
    j = wt.j + 1;
    w = wt.w;

    if (!job) {
        job=0n;
        out = array(double, nx, ny);
    } else if (job==1) {
        out = array(double, nx, ny, n);
    } else {
        trtt_error_display_msg, code=TRTT_ERR_UNSUPPORTED_JOB;
        return;
    }
    
    __interpolator, out, x, p, j, w, job;
   
    return out;
}

func _trtt_2D_time_interpolator(&out, x, p, j, w, job)
{
    // Allocate array for the result.
    if (!job) {
        job=0n;
    } else if (job==1) {
    } else {
        error, "unsupported JOB";
    }

    if (!job) {
        out = (w(-,-,)*x(..,j))(..,sum);
    } else if (job==1) {
        out(..,j) = w(-,-,)*x(..,-:1:p);
    } else {
        error, "unsupported JOB";
    }

    return TRTT_SUCCESS;
}

/*|----------------------------|*/
/*| DISTANCE DRIVEN PROJECTORS |*/
/*|----------------------------|*/

func trtt_single_distance_driven_2D_projector(X, Yk, mode, status=)
/* DOCUMENT trtt_single_distance_driven_2D_projector, X, Yk, mode, status=
   
   This function computes a linear operator of algebraic Radon transform,
   performing the "distance driven" method of DeMan & Basu (2004). It is created
   from a tomobj X and a tomdata Yk. It can be performed both in "parallel" and
   "fan beam" MODE of projection. The "distance driven" method is a combination
   of the "voxel" and "ray driven" methods, considering voxels and detectors,
   not as points anymore, but with the thickness of their support. The interest
   is to prevent the artifacts that are generated in the "pixel driven"
   projection and in the "ray driven" backprojection. The projector consists in
   a h_evaluator carrying references to X and Yk, associated with the adapted
   projection function. A STATUS can be given. If not it is created in the
   function and if an error occurs, it will be directly displayed.

   SEE ALSO: trtt_2D_projectors_info,
             _trtt_single_distance_driven_PB2D_projector,
             _trtt_single_distance_driven_FB2D_projector,
             h_evaluator.
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
    xotilt = Yk.xotilt;
    xos = Yk.xos;
    xsd = Yk.xsd;
    t_index = Yk.t_index;
    /* xform combinations */
    // xos_plus = trtt_xform3d_combine(xotilt, obj_xi);
    // trtt_xform3d_combine, xos_plus, xos, xos_plus;
    xos_plus = trtt_xform3d_combine(xos, xotilt);
    xos_coefs = xos_plus.a;
    xsd_coefs = xsd.a;

    /* Check t_index */
    if (!(t_index>=0 & t_index<=nt)) {
        trtt_error_set, status, TRTT_ERR_OUT_OF_RANGE_INDEX, msg="Incompatible TOMDATA 't_index' with dimension 'nt' of the TOMOBJ.", disp=flag_disp_status;
        return;
    }

    p = h_new(X=X,
              Yk=Yk,
              xos_coefs = xos_coefs,
              xsd_coefs = xsd_coefs);
        
    if (mode == TRTT_PARALLEL_BEAM) {
        h_set, p, projector = symlink_to_name("_trtt_single_distance_driven_PB2D_projector");
        h_set, p, sparser = symlink_to_name("_trtt_single_distance_driven_PB2D_projector_sparser");
        // h_evaluator, p, "_trtt_single_distance_driven_PB2D_projector";
    }      
    else if (mode == TRTT_FAN_BEAM || mode == TRTT_CONE_BEAM){
        h_set, p, projector = symlink_to_name("_trtt_single_distance_driven_FB2D_projector");
        h_set, p, sparser = symlink_to_name("_trtt_single_distance_driven_FB2D_projector_sparser");
        // h_evaluator, p, "_trtt_single_distance_driven_FB2D_projector";
    }
    else {
        trtt_error_set, status, TRTT_ERR_BAD_GEOMETRY, disp=flag_disp_status;
        return;
    }

    h_evaluator, p, "_trtt_single_2D_projector";

    return p;
}
D_single_2D = trtt_single_distance_driven_2D_projector;
// Define an easier alias


func _trtt_single_distance_driven_PB2D_projector(&out, x, nx, ny, s_scl, s_deg, footprint, cum_footprint, size_footprint, step_footprint, coord_footprint, nv, v_scl, xos_coefs, xsd_coefs, job)    
/* DOCUMENT _trtt_single_distance_driven_PB2D_projector, &out, x, nx, ny, s_scl, s_deg, nv, v_scl, xos_coefs, xsd_coefs, job
   
   Evaluator of the projector (direct and transpose) based on the "distance
   driven" method of DeMan & Basu (2002), in "parallel beam" mode of projection.
   
   SEE ALSO: trtt_single_distance_driven_2D_projector.
*/        
{
    /* - Calculated parameters */
    xn = (indgen(nx)-0.5*(nx+1))*s_scl;
    yn = (indgen(ny)-0.5*(ny+1))*s_scl;
    v = (indgen(nv)-0.5*(nv+1))*v_scl;

    /* Pre-declaration of local variables used as aliases. */
    local objX;

    // Allocate array for the result.
    if (!job) {
        job=0n;
        xt = transpose(x);
    } else if (job==1) {
    } else {
        error, "unsupported JOB";
    }

    // We need the boundaries of the pixels detector
    v_bounds = (indgen(nv+1)-0.5*(nv+2))*v_scl;
    border_down = v_bounds(1);
    border_up = v_bounds(0);
    // We need the boundaries of the voxels
    xn_bounds = (indgen(nx+1)-0.5*(nx+2))*s_scl;
    yn_bounds = (indgen(ny+1)-0.5*(ny+2))*s_scl;

    /* The mapping direction is simply the detector direction. Hence we only
     * have to project the voxels on it, and the pixels detector are already
     * "placed".
     */
                
    /* We have to choose the axis direction to project on the
     * detector direction (mapping direction). The axis direction to
     * be projected will involve the projection of the boundaries of
     * the voxel support in this direction (active variable). The
     * positions of the boundaries are projected from the center of
     * the voxel in the other axis direction (passive variable).
     */
    alpha5 = abs(xos_coefs(5));
    alpha6 = abs(xos_coefs(6));
    // FIXME: the comparison of the 5th and the 6th coefficients of XOS
    // determine the mapping direction.
    
    if (alpha5 < alpha6) {
        axis = 1;
        // Boundaries of the support of the voxel
        n_active = ny;
        n_passive = nx;

        /* Table of projected variables */
        proj_var = _trtt_2D_projection_O_to_D(xn,yn_bounds,xos_coefs,xsd_coefs,TRTT_PARALLEL_BEAM);
        proj_var = transpose(proj_var);
        
        if (!job) eq_nocopy, objX, xt;
    } else {
        axis = 2;
        // Boundaries of the support of the voxel
        n_active = nx;
        n_passive = ny;

        /* Table of projected variables */
        proj_var = _trtt_2D_projection_O_to_D(xn_bounds,yn,xos_coefs,xsd_coefs,TRTT_PARALLEL_BEAM);
        
        if (!job) eq_nocopy, objX, x;
    }
    
    /* Weights due to orientation : normalized height of the voxels */
    h_norm = (s_scl^2.0)/avg(abs(proj_var(dif,))(*));
        
    for (q=1; q<=n_passive; ++q) {
        
        if (!job) {
            // The object is sliced over the passive axis : objX(,j)
            vox_values = objX(,q,1);
        }
        
        /* Projection of the voxels and re-ordering
         * voxels values */
        projections = proj_var(,q);
        test_vox_order = anyof(projections(dif)<0.0); 
        if (test_vox_order) {
            i_order = sort(projections);
            projections = projections(i_order);
            if (!job) vox_values = vox_values(i_order(2:0));
        }
        voxpos_down = min(projections);
        voxpos_up = max(projections);
        
        /* Constraint the voxel projection to the
         * detector boundaries if necessary */
        idproj = where(projections>border_down & projections<border_up);
        /* Find the voxel neighbour indexes */
        iddet_out = where(v_bounds>=voxpos_up);
        ids_tmp = digitize(v_bounds, projections)-1;
        if (is_array(iddet_out))
            ids_tmp(iddet_out) = n_active+1;
        id_neigh_vox = ids_tmp;
        if (is_array(idproj))
            grow, id_neigh_vox, idproj;
        id_neigh_vox = id_neigh_vox(sort(id_neigh_vox));
        
        /* Find the detector neighbour indexes */
        if (is_array(idproj)) {
            projections = projections(idproj);
        }
        id_neigh_det = grow(indgen(nv+1), digitize(projections, v_bounds)-1);
        id_neigh_det = id_neigh_det(sort(id_neigh_det));
        
        /* Create the whole positions (voxels + pixels
         * detector) array */                               
        grow, projections, v_bounds;
        projections = projections(sort(projections));
        /* Differential coefficients */
        coefs = h_norm*projections(dif);
        Nco = numberof(coefs);

        idvox = where(id_neigh_vox>0 & id_neigh_vox<(n_active+1));

        if (is_array(idvox)) {

            if (idvox(0)>Nco)
                idvox = idvox(1:-1);
        
            if (!job) {
                /* Find the voxels values associated with the coefficients */
                voxels = array(double, numberof(id_neigh_vox)-1);
                voxels(idvox) = vox_values(id_neigh_vox(idvox));
                
                /* Update the Sinogram slice */
                out += histogram(id_neigh_det(1:-1), coefs*voxels, top=nv);
            } else {
                out_q = array(double, n_active);
                /* Find the detector values associated with the coefficients */
                sino = x(id_neigh_det(1:-1));
                /* Update the Object slice */
                
                out_q = histogram(id_neigh_vox(idvox), (coefs*sino)(idvox), top=n_active);
                if (test_vox_order) {
                    out_q = out_q(i_order(2:0));
                }
                if (axis == 1) {
                    out(q,,1) += out_q;
                } else {
                    out(,q,1) += out_q;
                }
            }
        }
    }
    
    return TRTT_SUCCESS;
}

func _trtt_single_distance_driven_FB2D_projector(&out, x, nx, ny, s_scl, s_deg, footprint, cum_footprint, size_footprint, step_footprint, coord_footprint, nv, v_scl, xos_coefs, xsd_coefs, job)
/* DOCUMENT _trtt_single_distance_driven_FB2D_projector, &out, x, nx, ny, s_scl, s_deg, nv, v_scl, xos_coefs, xsd_coefs, job

   Evaluator of the projector (direct and transpose) based on the "distance
   driven" method of DeMan & Basu (2002), in "fan beam" mode of projection.
   
   SEE ALSO: trtt_single_distance_driven_2D_projector.
*/        
{
    /* - Calculated parameters */
    xn = (indgen(nx)-0.5*(nx+1))*s_scl;
    yn = (indgen(ny)-0.5*(ny+1))*s_scl;
    v = (indgen(nv)-0.5*(nv+1))*v_scl;

    /* Pre-declaration of local variables used as aliases. */
    local objX;

    // Allocate array for the result.
    if (!job) {
        job=0n;
        xt = transpose(x);
    } else if (job==1) {
    } else {
        error, "unsupported JOB";
    }

    // We need the boundaries of the pixels detector
    v_bounds = (indgen(nv+1)-0.5*(nv+2))*v_scl;
    border_down = v_bounds(1);
    border_up = v_bounds(0);
    // We need the boundaries of the voxels
    xn_bounds = (indgen(nx+1)-0.5*(nx+2))*s_scl;
    yn_bounds = (indgen(ny+1)-0.5*(ny+2))*s_scl;

    /* The mapping direction is simply the detector direction. Hence we only
     * have to project the voxels on it, and the pixels detector are already
     * "placed".
     */
                
    /* We have to choose the axis direction to project on the
     * detector direction (mapping direction). The axis direction to
     * be projected will involve the projection of the boundaries of
     * the voxel support in this direction (active variable). The
     * positions of the boundaries are projected from the center of
     * the voxel in the other axis direction (passive variable).
     */

    alpha5 = abs(xos_coefs(5));
    alpha6 = abs(xos_coefs(6));
    // FIXME: the comparison of the 5th and the 6th coefficients of XOS
    // determine the mapping direction.

    /* Calculate the magnification factor of each voxel position */
    magnification = _trtt_2D_magnification_O_to_D (xn, yn, xos_coefs, xsd_coefs);
    distorsion = _trtt_2D_distorsion_O_to_D (xn, yn, xos_coefs, xsd_coefs);
    
    if (alpha5 < alpha6) {
        axis = 1;
        // Boundaries of the support of the voxel
        n_active = ny;
        n_passive = nx;
        /* Calculate the magnification factor of each voxel position */
        magnification = transpose(magnification);
        distorsion = transpose(distorsion);
        /* Table of projected variables */
        proj_var = _trtt_2D_projection_O_to_D(xn,yn_bounds,xos_coefs,xsd_coefs,TRTT_FAN_BEAM);
        proj_var = transpose(proj_var);
        
        if (!job) eq_nocopy, objX, xt;
    } else {
        axis = 2;
        // Boundaries of the support of the voxel
        n_active = nx;
        n_passive = ny;
        /* Table of projected variables */
        proj_var = _trtt_2D_projection_O_to_D(xn_bounds,yn,xos_coefs,xsd_coefs,TRTT_FAN_BEAM);
        
        if (!job) eq_nocopy, objX, x;
    }
    
    for (q=1; q<=n_passive; ++q) {
        
        if (!job) {
            // The object is sliced over the passive axis : objX(,j)
            vox_values = objX(,q,1);
        }
        
        /* Projection of the voxels and re-ordering
         * voxels values */
        projections = proj_var(,q);
        magn = magnification(,q);
        delta = distorsion(,q);
        test_vox_order = anyof(projections(dif)<0.0);
        if (test_vox_order) {
            i_order = sort(projections);
            projections = projections(i_order);
            magn = magn(i_order(2:0));
            delta = delta(i_order(2:0));
            if (!job) vox_values = vox_values(i_order(2:0));
        }
        voxpos_down = min(projections);
        voxpos_up = max(projections);
        /* Weights due to orientation (depend on each beam
         * direction) : normalized height of the voxels */
        h_norm = delta*magn*(s_scl^2.0)*(1.0/projections(dif));
        
        /* Constraint the voxel projection to the
         * detector boundaries if necessary */
        idproj = where(projections>border_down & projections<border_up);
        /* Find the voxel neighbour indexes */
        iddet_out = where(v_bounds>=voxpos_up);
        ids_tmp = digitize(v_bounds, projections)-1;
        if (is_array(iddet_out))
            ids_tmp(iddet_out) = n_active+1;
        id_neigh_vox = ids_tmp;
        if (is_array(idproj))
            grow, id_neigh_vox, idproj;
        id_neigh_vox = id_neigh_vox(sort(id_neigh_vox));
        
        /* Find the detector neighbour indexes */
        if (is_array(idproj)) {
            projections = projections(idproj);
        }
        id_neigh_det = grow(indgen(nv+1), digitize(projections, v_bounds)-1);
        id_neigh_det = id_neigh_det(sort(id_neigh_det));
        
        /* Create the whole positions (voxels + pixels
         * detector) array */
        grow, projections, v_bounds;
        projections = projections(sort(projections));
        /* Differential coefficients */
        coefs = projections(dif);
        Nco = numberof(coefs);

        idvox = where(id_neigh_vox>0 & id_neigh_vox<(n_active+1));
        
        if (is_array(idvox)) {
            
            if (idvox(0)>Nco)
                idvox = idvox(1:-1);
            
            if (!job) {
                /* Find the voxels values associated with the coefficients */
                voxels = array(double, numberof(id_neigh_vox)-1);
                voxels(idvox) = (h_norm*vox_values)(id_neigh_vox(idvox));
                
                /* Update the Sinogram slice */
                out += histogram(id_neigh_det(1:-1), coefs*voxels, top=nv);
            } else {
                out_q = array(double, n_active);
                /* Find the detector values associated with the coefficients */
                sino = x(id_neigh_det(1:-1));
                /* Update the Object slice */
                out_q = histogram(id_neigh_vox(idvox),
                                  h_norm(id_neigh_vox(idvox))*(coefs*sino)(idvox),
                                  top=n_active);
                if (test_vox_order) {
                    out_q = out_q(i_order(2:0));
                }
                if (axis == 1) {
                    out(q,,1) += out_q;
                } else {
                    out(,q,1) += out_q;
                }
            }
        }
    }
            
    return TRTT_SUCCESS;
}

func trtt_whole_distance_driven_2D_projector(X, Y, mode, status=)
/* DOCUMENT trtt_whole_distance_driven_2D_projector, X, Y, mode, status=

   Create a set of 2D projectors based on the Distance Driven model, from a
   TOMOGRAPHIC OBJECT X and a set of DATA objects Y, according to the geometry
   of projection MODE. A STATUS can be given. If not it is created in the
   function and if an error occurs, it will be directly displayed.

   SEE ALSO: trtt_2D_projectors_info, trtt_single_distance_driven_2D_projector.
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

        Ck = D_single_2D(X, Yk, mode, status=status);
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
D_whole_2D = trtt_whole_distance_driven_2D_projector;

/*|------------------------|*/
/*| 2D COMBINED PROJECTORS |*/
/*|------------------------|*/

func trtt_single_2D_combined_projector(X, Yk, mode, type, status=)
/* DOCUMENT trtt_single_2D_combined_projector, X, Yk, mode, type, status=

   Create a 2D+t projector, combining a 2D spatial projector with the model
   given by TYPE (TRTT_SPLINE_DRIVEN or TRTT_DISTANCE_DRIVEN), and a time
   interpolator, according to the geometry of projection MODE
   (TRTT_PARALLEL_BEAM or TRTT_FAN(CONE)_BEAM), to perform Dynamic tomography. X
   is the TOMOGRAPHIC OBJECT and Yk is the DATA object. The projector consists in
   a h_evaluator carrying references to the projector Ck (carrying references to
   X and Yk) and and the time interpolator Bk, associated with the projection
   function. A STATUS can be given. If not it is created in the function and if
   an error occurs, it will be directly displayed.

   SEE ALSO: trtt_2D_projectors_info, _trtt_single_2D_combined_projector, h_evaluator.
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

    if (type == TRTT_SPLINE_DRIVEN) {
        Ak = A_single_2D(X, Yk, mode, status=status);
        if(is_void(Ak)) {
            if (flag_disp_status) trtt_error_display_msg, status=status;
            return;
        }
    } else if (type == TRTT_DISTANCE_DRIVEN) {
        Ak = D_single_2D(X, Yk, mode, status=status);
        if(is_void(Ak)) {
            if (flag_disp_status) trtt_error_display_msg, status=status;
            return;
        }
    }

    Bk = B_single_2D(X, Yk, status=status);
    if(is_void(Bk)) {
        if (flag_disp_status) trtt_error_display_msg, status=status;
        return;
    }

    p = h_new(Ak = Ak,
              Bk = Bk);

    h_evaluator, p, "_trtt_single_2D_combined_projector";

    return p;
}
C_single_2D = trtt_single_2D_combined_projector;

func _trtt_single_2D_combined_projector(p, x, job)
/* DOCUMENT _trtt_single_2D_combined_projector, p, x, job

   Evaluator (direct and tranpose) of a 2D combined projector. The function
   applies a combination of a spatial 2D projector Ck and a time interpolator
   Bk.

   JOB=0 (direct)

        out = Ck[Bk(x)]

   JOB=1 (transpose)

        out = Bk[Ck(x,1),1]

   SEE ALSO: trtt_single_2D_combined_projector, h_evaluator.
 */ 
{
    Ak = p.Ak;
    Bk = p.Bk;

    if (!job) {
        out = Ak(Bk(x));
    } else if (job==1) {
        out = Bk(Ak(x,1),1);
    } else {
        trtt_error_display_msg, code=TRTT_ERR_UNSUPPORTED_JOB;
        return;
    }

    return out;
}

func trtt_whole_2D_combined_projector(X, Y, mode, type, status=)
/* DOCUMENT trtt_whole_2D_combined_projector, X, Y, mode, type, status=

   Create a set of 2D combined projectors, from a TOMOGRAPHIC OBJECT X and a set
   of DATA objects Y, according to the geometry of projection MODE. A STATUS can
   be given. If not it is created in the function and if an error occurs, it
   will be directly displayed.

   SEE ALSO: trtt_2D_projectors_info, trtt_single_2D_combined_projector.
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

        Ck = C_single_2D(X, Yk, mode, type, status=status);
        if(is_void(Ck)) {
            trtt_error_display_msg, status=status;
            return;
        }

        h_set, Cops, trtt_get_key_name(k, "C"), Ck;
    }

    Ck_list = h_keys(Cops);
    h_set, Cops, Ck_list = Ck_list(sort(Ck_list));
    
    return Cops;
}
C_whole_2D = trtt_whole_2D_combined_projector;

/* DISTANCE DRIVEN SPARSER FUNCTIONS  ======================================== */

func _trtt_single_distance_driven_PB2D_projector_sparser(&out, i, j, nx, ny, s_scl, s_deg, footprint, cum_footprint, size_footprint, step_footprint, coord_footprint, nv, v_scl, xos_coefs, xsd_coefs)    
/* DOCUMENT _trtt_single_distance_driven_PB2D_projector_sparser, &out, i, j, nx, ny, s_scl, s_deg, nv, v_scl, xos_coefs, xsd_coefs
   
   Evaluator of the projector (direct and transpose) based on the "distance
   driven" method of DeMan & Basu (2002), in "parallel beam" mode of projection.
   
   SEE ALSO: trtt_single_distance_driven_2D_projector.
*/        
{
    /* - Calculated parameters */
    xn = ((indgen(nx)-0.5*(nx+1))*s_scl)(i+1);
    yn = ((indgen(ny)-0.5*(ny+1))*s_scl)(j+1);
    xn_bounds = ((indgen(nx+1)-0.5*(nx+2))*s_scl)(i+1:i+2);
    yn_bounds = ((indgen(ny+1)-0.5*(ny+2))*s_scl)(j+1:j+2);
    
    // We need the boundaries of the pixels detector
    v = (indgen(nv)-0.5*(nv+1))*v_scl;
    v_bounds = (indgen(nv+1)-0.5*(nv+2))*v_scl;
    border_down = v_bounds(1);
    border_up = v_bounds(0);
    
    /* The mapping direction is simply the detector direction. Hence we only
     * have to project the voxels on it, and the pixels detector are already
     * "placed".
     */
    
    /* We have to choose the axis direction to project on the
     * detector direction (mapping direction). The axis direction to
     * be projected will involve the projection of the boundaries of
     * the voxel support in this direction (active variable). The
     * positions of the boundaries are projected from the center of
     * the voxel in the other axis direction (passive variable).
     */
    alpha5 = abs(xos_coefs(5));
    alpha6 = abs(xos_coefs(6));
    // FIXME: the comparison of the 5th and the 6th coefficients of XOS
    // determine the mapping direction.
    
    if (alpha5 < alpha6) {
        axis = 1;
        // Boundaries of the support of the voxel
        n_active = ny;
        n_passive = nx;
        /* Table of projected variables */
        projections = _trtt_2D_projection_O_to_D(xn,yn_bounds,xos_coefs,xsd_coefs,TRTT_PARALLEL_BEAM)(1,);
    } else {
        axis = 2;
        // Boundaries of the support of the voxel
        n_active = nx;
        n_passive = ny;
        /* Table of projected variables */
        projections = _trtt_2D_projection_O_to_D(xn_bounds,yn,xos_coefs,xsd_coefs,TRTT_PARALLEL_BEAM);
    }
    
    /* Weights due to orientation : normalized height of the voxels */
    h_norm = (s_scl^2.0)/avg(abs(projections(dif))(*));
    
    projections = projections(sort(projections));
    projections(1) = max(projections(1),border_down);
    projections(2) = min(projections(2),border_up);
    
    if (projections(1)<projections(2)) {
        idpixfirst = (digitize(projections,v_bounds)-1)(1);
        
        idpixdet = where(v_bounds>projections(1) & v_bounds<projections(2));
        if (is_array(idpixdet)) {
            grow, projections, v_bounds(idpixdet);
            projections = projections(sort(projections));
        }
        
        coefs = h_norm*projections(dif);
        nco = numberof(coefs);
        /* Update the Sinogram slice */
        out(idpixfirst:idpixfirst+nco-1) += coefs;
    }
    
    return TRTT_SUCCESS;
}

func _trtt_single_distance_driven_FB2D_projector_sparser(&out, i, j, nx, ny, s_scl, s_deg, footprint, cum_footprint, size_footprint, step_footprint, coord_footprint, nv, v_scl, xos_coefs, xsd_coefs)    
/* DOCUMENT _trtt_single_distance_driven_FB2D_projector_sparser, &out, i, j, nx, ny, s_scl, s_deg, nv, v_scl, xos_coefs, xsd_coefs
   
   Evaluator of the projector (direct and transpose) based on the "distance
   driven" method of DeMan & Basu (2002), in "parallel beam" mode of projection.
   
   SEE ALSO: trtt_single_distance_driven_2D_projector.
*/        
{
    /* - Calculated parameters */
    xn = ((indgen(nx)-0.5*(nx+1))*s_scl)(i+1);
    yn = ((indgen(ny)-0.5*(ny+1))*s_scl)(j+1);
    xn_bounds = ((indgen(nx+1)-0.5*(nx+2))*s_scl)(i+1:i+2);
    yn_bounds = ((indgen(ny+1)-0.5*(ny+2))*s_scl)(j+1:j+2);
    
    // We need the boundaries of the pixels detector
    v = (indgen(nv)-0.5*(nv+1))*v_scl;
    v_bounds = (indgen(nv+1)-0.5*(nv+2))*v_scl;
    border_down = v_bounds(1);
    border_up = v_bounds(0);

    /* The mapping direction is simply the detector direction. Hence we only
     * have to project the voxels on it, and the pixels detector are already
     * "placed".
     */
                
    /* We have to choose the axis direction to project on the
     * detector direction (mapping direction). The axis direction to
     * be projected will involve the projection of the boundaries of
     * the voxel support in this direction (active variable). The
     * positions of the boundaries are projected from the center of
     * the voxel in the other axis direction (passive variable).
     */
    alpha5 = abs(xos_coefs(5));
    alpha6 = abs(xos_coefs(6));
    // FIXME: the comparison of the 5th and the 6th coefficients of XOS
    // determine the mapping direction.

    /* Calculate the magnification factor of each voxel position */
    magn = _trtt_2D_magnification_O_to_D (xn, yn, xos_coefs, xsd_coefs)(1);
    delta = _trtt_2D_distorsion_O_to_D (xn, yn, xos_coefs, xsd_coefs)(1);
    
    if (alpha5 < alpha6) {
        axis = 1;
        // Boundaries of the support of the voxel
        n_active = ny;
        n_passive = nx;
        /* Table of projected variables */
        projections = _trtt_2D_projection_O_to_D(xn,yn_bounds,xos_coefs,xsd_coefs,TRTT_FAN_BEAM)(1,);
    } else {
        axis = 2;
        // Boundaries of the support of the voxel
        n_active = nx;
        n_passive = ny;
        /* Table of projected variables */
        projections = _trtt_2D_projection_O_to_D(xn_bounds,yn,xos_coefs,xsd_coefs,TRTT_FAN_BEAM);
    }
    
    /* Projection of the voxels and re-ordering
     * voxels values */
    projections = projections(sort(projections));
    projections(1) = max(projections(1),border_down);
    projections(2) = min(projections(2),border_up);
        
    /* Weights due to orientation (depend on each beam
         * direction) : normalized height of the voxels */
    h_norm = delta*magn*(s_scl^2.0)*(1.0/projections(dif));
    
    if (projections(1)<projections(2)) {
        idpixfirst = (digitize(projections,v_bounds)-1)(1);
        
        idpixdet = where(v_bounds>projections(1) & v_bounds<projections(2));
        if (is_array(idpixdet)) {
            grow, projections, v_bounds(idpixdet);
            projections = projections(sort(projections));
        }
        
        coefs = h_norm*projections(dif);
        nco = numberof(coefs);
        /* Update the Sinogram slice */
        out(idpixfirst:idpixfirst+nco-1) += coefs;
    }
    
    return TRTT_SUCCESS;
}

/* OBSOLETE FUNCTIONS  ======================================================= */


func _trtt_2D_sparse_coefs_old(Cops, dyn=)
/* DOCUMENT _trtt_2D_sparse_coefs, Cops, dyn=

   Calculate the sparse projector matrix coefficients from each projector Ck of
   the set COPS. Sparse matrix components are added to each Ck in order not have
   to re-calculate it everytime.
 */
{
    Ck_list = Cops.Ck_list;
    ndata = numberof(Ck_list);
    if (!dyn) {
        X = h_get(Cops,Ck_list(1)).X;
    } else {
        X = h_get(Cops,Ck_list(1)).Ak.X;
    }
    nx = X.nx;
    ny = X.ny;

    if (!Cops.exist_sparse_coefs) {
        /* Get the sparse matrix parameters for each projector */
        for (k=1; k<=ndata; ++k) {
            /* Get the current projector */
            if (!dyn) {
                Ck = h_get(Cops,Ck_list(k));
            } else {
                Ck = h_get(Cops,Ck_list(k)).Ak;
            }
            /* Get the current detector size */
            nv = Ck.Yk.nv;

            MCk_coefs = [];
            MCk_row_indices = [];
            MCk_col_indices = [];
            MCk_row_dimlist = [1,nv];
            MCk_col_dimlist = [1,nx*ny];
            
            for (q=1; q<=nv; ++q) {
                /* Get the current canonical vector */
                sino = array(double,nv);
                sino(q) = 1.0;
                /* Calculate the matrix coefficients */
                matrix_coefs = Ck(sino,1)(*);
                /* Get the sparse indices */
                sparse_indices = where(matrix_coefs != 0.0);
                if(is_array(sparse_indices)) {
                    nsparse = numberof(sparse_indices);
                    /* Grow sparse matrix coefficients */
                    grow, MCk_coefs, matrix_coefs(sparse_indices);
                    /* Grow sparse matrix column indices */
                    grow, MCk_col_indices, sparse_indices;
                    /* Grow sparse matrix row indices */
                    grow, MCk_row_indices, array(long(q), nsparse);
                }
            }
            /* Store sparse matrix parameters */
            h_set, Ck, MCk_coefs = MCk_coefs;
            h_set, Ck, MCk_row_dimlist = MCk_row_dimlist;
            h_set, Ck, MCk_row_indices = MCk_row_indices;
            h_set, Ck, MCk_col_dimlist = MCk_col_dimlist;
            h_set, Ck, MCk_col_indices = MCk_col_indices;
        }
        h_set, Cops, exist_sparse_coefs = 1n;
    }

    return TRTT_SUCCESS;
}

func _trtt_2D_sparse_matrix_old(Cops)
/* DOCUMENT _trtt_2D_sparse_matrix, Cops

   Calculate the sparse projector matrix coefficients from the set of projectors
   COPS and return the sparse matrix. Sparse matrix components are added to COPS
   in order not have to re-calculate it everytime.
 */
{
    Ck_list = Cops.Ck_list;
    ndata = numberof(Ck_list);
    nv = h_get(Cops,Ck_list(1)).Yk.nv;
    if (!Cops.matrix) {
        /* FIXME: First method creating the whole matrix */
        /*-----------------------------------------------*/
        /* Calculate projector matrix and data sinogram */
        /* Initialize whole projection matrix */
        // MCops = array(double,nv*ndata,nx*ny);
        /* Initialize whole sinogram data */
        // y = array(double,nv,ndata);
        // for (k=1; k<=ndata; ++k) {
        /* Get the current projector */
        //     Ck = h_get(Cops,Ck_list(k));
        /* Get and store the current projection data */
        //     Yk = Ck.Yk;
        //     y(,k)=Yk.pixels;
        //     for (q=1; q<=nv; ++q) {
        /* Get the current canonical vector */
        //         sino = array(double,nv);
        //         sino(q) = 1.0;
        /* Calculate and store the current column of the matrix */
        //         MCops((k-1)*nv+q,) = Ck(sino,1)(*);
        //     }
        // }

        /* Create the sparse matrix */
        // sparse_MCops = sparse_squeeze(MCops);

        /* Store sparse matrix parameters */
        // h_set, Cops, MCops_coefs = sparse_MCops.coefs;
        // h_set, Cops, MCops_row_dimlist = sparse_MCops.row_dimlist;
        // h_set, Cops, MCops_row_indices = sparse_MCops.row_indices;
        // h_set, Cops, MCops_col_dimlist = sparse_MCops.col_dimlist;
        // h_set, Cops, MCops_col_indices = sparse_MCops.col_indices;
        // h_set, Cops, y = y;
        // h_set, Cops, matrix = 1n;


        /* FIXME: Other method which allow higher dimensioned systems */
        /* FIXME: BUT VERY LONGER !!!                                 */
        /*------------------------------------------------------------*/
        /* Initialize sparse matrix parameters */
        MCops_coefs = [];
        MCops_row_indices = [];
        MCops_col_indices = [];
        MCops_row_dimlist = [1,nv*ndata];
        MCops_col_dimlist = [1,nx*ny];
        /* Initialize whole sinogram data */
        y = array(double,nv,ndata);
        for (k=1; k<=ndata; ++k) {
            /* Get the current projector */
            Ck = h_get(Cops,Ck_list(k));
            /* Get and store the current projection data */
            Yk = Ck.Yk;
            y(,k)=Yk.pixels;
            for (q=1; q<=nv; ++q) {
                /* Get the current canonical vector */
                sino = array(double,nv);
                sino(q) = 1.0;
                /* Calculate the matrix coefficients */
                matrix_coefs = Ck(sino,1)(*);
                /* Get the sparse indices */
                sparse_indices = where(matrix_coefs != 0.0);
                nsparse = numberof(sparse_indices);
                /* Grow sparse matrix coefficients */
                grow, MCops_coefs, matrix_coefs(sparse_indices);
                /* Grow sparse matrix column indices */
                grow, MCops_col_indices, sparse_indices;
                /* Grow sparse matrix row indices */
                grow, MCops_row_indices, array(long((k-1)*nv+q), nsparse);
            }
        }
        /* Store sparse matrix parameters */
        h_set, Cops, MCops_coefs = MCops_coefs;
        h_set, Cops, MCops_row_dimlist = MCops_row_dimlist;
        h_set, Cops, MCops_row_indices = MCops_row_indices;
        h_set, Cops, MCops_col_dimlist = MCops_col_dimlist;
        h_set, Cops, MCops_col_indices = MCops_col_indices;
        h_set, Cops, matrix = 1n;
        /* Store sinogram */
        h_set, Cops, y = y;

        /* Create the sparse matrix */
        sparse_MCops = sparse_matrix(MCops_coefs,
                                     MCops_row_dimlist, MCops_row_indices,
                                     MCops_col_dimlist, MCops_col_indices);

    } else {
        /* Get sparse matrix parameters */
        local coefs; eq_nocopy, coefs, Cops.MCops_coefs;
        local row_dimlist; eq_nocopy, row_dimlist, Cops.MCops_row_dimlist;
        local row_indices; eq_nocopy, row_indices, Cops.MCops_row_indices;
        local col_dimlist; eq_nocopy, col_dimlist, Cops.MCops_col_dimlist;
        local col_indices; eq_nocopy, col_indices, Cops.MCops_col_indices;

        /* Create the sparse matrix */
        sparse_MCops = sparse_matrix(coefs,
                                     row_dimlist, row_indices,
                                     col_dimlist, col_indices);
    }

    return sparse_MCops;
}
