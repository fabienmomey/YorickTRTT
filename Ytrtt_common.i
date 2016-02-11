/*
 * TRTT : iTeraTive Reconstruction for Tomography
 *
 * Ytrtt_common.i --
 *
 * Common features for TRTT.
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

/* GLOBALS =================================================================== */

TRTT_WATER_ABSORPTION = 1.0; //1.928;
/* DOCUMENT Absorption of water (in mm^-1) for a photon energy of 70keV. */

TRTT_FALSE = 0n;
/* DOCUMENT Boolean value symbolizing a FALSE answer to a test. */

TRTT_TRUE = 1n;
/* DOCUMENT Boolean value symbolizing a TRUE answer to a test. */

TRTT_SUCCESS = 0;
/* DOCUMENT Returned value upon success. */

TRTT_FAILURE = -1;
/* DOCUMENT Returned value upon failure. */

/* Types of Projectors */
TRTT_SPLINE_DRIVEN = 1n;
/* DOCUMENT Spline Driven projector type. */
TRTT_DISTANCE_DRIVEN = 2n;
/* DOCUMENT Distance Driven projector type. */
TRTT_LONG_FESSLER = 3n;
/* DOCUMENT Long Fessler projector type. */
TRTT_LONG_FESSLER_LINE_INTEG = 4n;
/* DOCUMENT Long Fessler projector type. */

/* Geometries for Projectors*/
TRTT_PARALLEL_BEAM = 3n;
/* DOCUMENT Parallel projection mode. */
TRTT_CONE_BEAM = 4n;
/* DOCUMENT Conic projection mode. */
TRTT_FAN_BEAM = 5n;
/* DOCUMENT Fan projection mode. */

/* FUNCTIONS ================================================================= */

func _trtt_is_1dim(n)
/* DOCUMENT _trtt_is_1dim, n

   Checks if N is compatible with a dimension (scalar integer and strictly
   positive). It returns the result of the test (0 for FALSE or 1 for TRUE).

   SEE ALSO: _trtt_is_1scale, _trtt_is_1deg.
*/
{
    return (is_integer(n) & is_integer(n) & n>=1);
}
is_1dim = _trtt_is_1dim;

func _trtt_is_1scale(scl)
/* DOCUMENT _trtt_is_1scale, n

   Checks if SCL is compatible with a sampling step (scalar real and
   positive). It returns the result of the test (0 for FALSE or 1 for TRUE).

   SEE ALSO: _trtt_is_1dim, _trtt_is_1deg.
*/
{
    return (is_real(scl) & is_scalar(scl) & scl>=0);
}
is_1scale = _trtt_is_1scale;

func _trtt_is_1deg(deg)
/* DOCUMENT _trtt_is_1deg, n

   Checks if DEG is compatible with a spline degree (scalar integer between 0
   and 9). It returns the result of the test (0 for FALSE or 1 for TRUE).

   SEE ALSO: _trtt_is_1dim, _trtt_is_1scale.
*/
{
    return (is_integer(deg) & is_scalar(deg) & (deg>=0 && deg<=9));
}
is_1deg = _trtt_is_1deg;

func trtt_get_key_name(idx, ident_string)
/* DOCUMENT trtt_get_key_name, idx, ident_string
   
   Returns a name of key made of the IDENT_STRING and the given index IDX.
   Example: trtt_get_key_name(2, "A") could return string "A3"; */ 
{
    return swrite(format="%s%04d", ident_string, idx);
}

func trtt_recursive_h_delete(h, keyname)
/* DOCUMENT recursive_h_delete, h, keyname
   
   Deletes the key KEYNAME recursively from hash table H. Therefore if H
   contains hash tables, KEYNAME will be deleted from these ones too.  */
{
    keylist = h_keys(h);
    n = numberof(keylist);
    for (i=1;i<=n;i++)
    {
        if (keylist(i)==keyname) h_delete, h, keyname;
                
        if (is_hash (h_get(h,keylist(i)))) {
            h_rec = h_get(h,keylist(i));
            trtt_recursive_h_delete, h_rec, keyname;
        }
    }
}

func trtt_span(x1, xn, n)
/* DOCUMENT trtt_span, x1, xn, n
   
   Returns array of N doubles equally spaced from X1 to XN. It uses the function
   INDGEN of Yorick. If n=1, it returns [X1].  */
{
    if (n > 1)
        return ((xn-x1)/(n-1.0))*(indgen(n)-(1.0+n)/2.0)+(x1+xn)/2.0;
    else
        return [x1];
}
