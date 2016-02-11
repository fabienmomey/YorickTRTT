/*
 * TRTT : iTeraTive Reconstruction for Tomography
 *
 * Ytrtt.i --
 *
 * TRTT modules loader.
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

extern TRTT_PATH;
TRTT_PATH = get_env("TRTT_PATH");

if (TRTT_PATH!=string()) {
    TRTT_PATH += "/";

    include, "yeti.i", 1;
    include, "yeti_yhdf.i", 1;
    // include, "OptimPack1_special.i", 1; //FIXME: modification de op_mnb
    include, "OptimPack1.i", 1;
    include, "opky.i", 1; // NOUVEL OPTIMPACK !!!
    // include, "lbfgsb.i", 1; // L-BFGS-B reverse communication
    include, "totvar.i", 1;
    include, "tomo_model.i", 1; // Modèle Long & Fessler

    include, "cmap.i", 1;
    include, "array.i", 1;
    include, "plot.i", 1;
    include, "utils.i";
    include, "linop.i", 1;
    // include, "rgl.i", 1;
    include, "img.i", 1;
    include, "gfx.i", 1;
    include, "yspl.i", 1;
    include, "fft_utils.i", 1;
    include, "romberg.i", 1;
    include, "Fabien/myYorick_routines.i", 1;
    include, "Michel/myplot.i", 1;
    include, "colormap.i", 1;
        
    include, TRTT_PATH+"plugin/Ytrtt_plugin.i", 1;
    include, TRTT_PATH+"Ytrtt_common.i", 1;
    include, TRTT_PATH+"Ytrtt_error.i", 1;
    include, TRTT_PATH+"Ytrtt_tomobj.i", 1;
    include, TRTT_PATH+"Ytrtt_tomdata.i", 1;
    include, TRTT_PATH+"Ytrtt_2D_projectors.i", 1;
    include, TRTT_PATH+"Ytrtt_2D_projectors_Long_Fessler.i", 1;
    include, TRTT_PATH+"Ytrtt_2D_projectors_Long_Fessler_line_integ.i", 1;
    include, TRTT_PATH+"Ytrtt_2D_projectors_comparator.i", 1;
    include, TRTT_PATH+"Ytrtt_3D_projectors.i", 1;
} else {
    error, "No path defined for YoricKTRTT code => << $ export TRTT_PATH=... >>";
}
