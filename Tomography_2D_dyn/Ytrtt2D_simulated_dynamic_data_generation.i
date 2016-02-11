/*
 * TRTT : iTeraTive Reconstruction for Tomography
 *
 * Ytrtt_simulated_data_generation.i --
 *
 * Tools pour generating simulated data sets from ellispoids for TRTT.
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

/* INTERNAL FUNCTIONS ======================================================== */

// func _trtt_periodic_spline(t, a, f, d)
// {
//     T = 1/f;
//     if (d>-1) {
//         n = (long(2*t*f))%2;
//         id1 = where(n==0);
//         id2 = where(n==1);
//         t = t%(0.5*T);
//         sig = array(double, dimsof(t));
//         sig(id1) = a * (4.*f)^d * ( (0.25*T) - abs(t(id1)-0.25*T) )^d;
//         sig(id2) = a * (4.*f)^d * ( abs(t(id2)-0.25*T)^d - (0.25*T)^d );
//     } else {
//         sig = a * sin(2.*pi*f*t);
//     }
    
//     return sig;
// }

func _trtt_periodic_spline(t, a, f, d)
{
    T = 1/f;

    if (d>-1) {
        n = (long(2*t*f))%2;
        t = t%(0.5*T);
        return a * (-1)^(1-n) * (4.*f)^d * (abs(t-0.25*T)^d - (0.25*T)^d);
    } else if (d==-1) {
        return a * sin(2.*pi*f*t);
    } else {
        return a * cos(2.*pi*f*t);
    }
}

/* USER FUNCTIONS ============================================================ */

func trtt_create_ellipses_dynamic_phantom(Ellipses, nx, ny, s_scl, xoff, yoff, tacq, ndata) 
/* DOCUMENT trtt_create_ellipses_phantom_rigid_motion, Ellipses, nx, ny, s_scl, xoff, yoff, tacq, ndata
   This function creates the 2D phantom image using the geometrical parameters
   of a 2D tomobj. A set of ellipses ELLIPSES is passed to the
   function. ELLIPSES is a hash table of the form:

   Ellipses           
   |
   |
   |----- E1
   |       |----- x0     // x-coordinate of the center in meters
   |       |----- y0     // y-coordinate of the center in meters
   |       |----- a      // length of the semimajor axis in meters
   |       |----- b      // length of the semiminor axis in meters
   |       |----- alpha  // angle of rotation of the major axis from
                         // the x-axis in radians
   |       |----- mu     // attenuation (constant on the surface of
                         // the ellipse)
   |
   |
   |----- E2
   |       |----- x0        
   |       |----- y0       
   |       |----- a         
   |       |----- b         
   |       |----- alpha   
   |       |----- mu      
   |
   |
   |----- E3
   :
   :  
   :
              
   This hash table requires a new entry per ellipse.

   SEE ALSO: CT_system_define.
*/
{ 
    E_list = Ellipses.Ellipses_list;
    nE = numberof(E_list);

    /* Time coordinates */
    t = (indgen(ndata)-1.0)*(tacq/(ndata-1.0));
    t = t(-,-,);

    /* Calculate coordinates */
    x = (indgen(nx)-((nx+1.0)/2.0))*s_scl;
    y = (indgen(ny)-((ny+1.0)/2.0))*s_scl;

    x = x(,-,-);
    y = y(-,,-);

    // Create the phantom image
    Im = array(double, nx, ny, ndata);


    for (u=1; u<=nE; ++u)
    {
        // Get the current ellipse
        E = h_get(Ellipses, E_list(u));
                
        // Extract its geometrical parameters on all the frames (3rd
        // dimension)
        x0 = E.x0 + xoff;
        y0 = E.y0 + yoff;
        a = E.a /*- 0.5*dx*/;
        b = E.b /*- 0.5*dy*/;
        alpha = E.alpha;
        mu = E.mu*TRTT_WATER_ABSORPTION; /* FIXME: conversion to true absorption */
        tr_x = E.tr_x;
        tr_y = E.tr_y;
        dis_a = E.dis_a;
        dis_b = E.dis_b;
        
        if (is_array(tr_y)) {
            y0 += _trtt_periodic_spline(t, tr_y(1), tr_y(2), tr_y(3));
        }
        if (is_array(tr_x)) {
            x0 += _trtt_periodic_spline(t, tr_x(1), tr_x(2), tr_x(3));
        }
        if (is_array(dis_a)) {
            a += _trtt_periodic_spline(t, dis_a(1), dis_a(2), dis_a(3));
        }
        if (is_array(dis_b)) {
            b += _trtt_periodic_spline(t, dis_b(1), dis_b(2), dis_b(3));
        }
                
        // Calculate the coefficients of the equation of the ellipse
        A = (cos(alpha) / a)^2 + (sin(alpha) / b)^2;
        B = (sin(alpha) / a)^2 + (cos(alpha) / b)^2;
        C = 2 * cos(alpha) * sin(alpha) * ((1/a^2) - (1/b^2));
                
        // Resolve the equation by identifying the pixels included in the ellipse
        eq = (A*(x-x0)^2 + B*(y-y0)^2 + C*(x-x0)*(y-y0) <= 1);
        // Add the result to the phantom image with the attenuation value
        Im += mu * eq;
    }
    
    return Im;
}
trtt_create_dyn_phantom = trtt_create_ellipses_dynamic_phantom;

func trtt_ellipses_dynamic_parallel_beam_analytic_Radon_transform(Ellipses, v, nv, v_scl, theta, t, xoff, yoff, epsilon=)
/* DOCUMENT trtt_ellipses_dynamic_parallel_beam_analytic_Radon_transform, Ellipses, v, nv, v_scl, theta, t, xoff, yoff, epsilon=
   
   This function calculates analytically the Radon transform, in PARALLEL BEAM
   geometry, of a 2D tomobj defined as a set of ellipses (for example a
   Shepp-Logan phantom). All the parameters are extracted from the set of
   ellipses ELLIPSES. ELLIPSES is a hash table of the form:

   Ellipses           
   |
   |
   |----- E1
   |       |----- x0     // x-coordinate of the center in meters
   |       |----- y0     // y-coordinate of the center in meters
   |       |----- a      // length of the semimajor axis in meters
   |       |----- b      // length of the semiminor axis in meters
   |       |----- alpha  // angle of rotation of the major axis from
                         // the x-axis in radians
   |       |----- mu     // attenuation (constant on the surface of
                         // the ellipse)
   |
   |
   |----- E2
   |       |----- x0        
   |       |----- y0       
   |       |----- a         
   |       |----- b         
   |       |----- alpha   
   |       |----- mu      
   |
   |
   |----- E3
   :
   :  
   :

   The Radon transform is determined for NV detector pixels at positions V
   sampled at step V_SCL, for an angle of projection THETA. It is calculated on
   band integrals (thick beam). It is necessary to give data compatible with the
   algebraic projectors distance or spline driven.
              
   The function calculates the Radon transform of this set of ellipses which is
   the sum of the Radon transform of each single ellipse Ei:
   
   * Parallel beam : P (theta, v) = SUM_i [ P_Ei (theta, v) ]
   = SUM_i [ L_i (theta, v) * mu_i ]

   where L_i (theta, v) is the length of intersection of the ray in (theta, v)
   and the ellipse Ei, and mu_i is the uniform attenuation on Ei. THETA is the
   projection angle and V is the position on the detector.

   The function returns the sinogram of the Radon transform.

   SEE ALSO: trtt_ellipses_fan_beam_analytic_Radon_transform.
*/        
{
    E_list = Ellipses.Ellipses_list;
    nE = numberof(E_list);
    
    if (is_void(epsilon)) epsilon=1.e-6;
    max_doublings = 50;

    // Create the sinogram
    Sinogram = array(double, nv);
              
    for (i=1; i<=nE; i++)
    {
        // Get the current ellipse
        E = h_get(Ellipses, E_list(i));

        tr_x = E.tr_x;
        tr_y = E.tr_y;
        tr_alpha = E.tr_alpha;
        dis_a = E.dis_a;
        dis_b = E.dis_b;

        Enew = h_copy(E);

        if (is_array(tr_y)) {
            newy0 = E.y0 + _trtt_periodic_spline(t, tr_y(1), tr_y(2), tr_y(3));
            h_set, Enew, y0=newy0;
        }
        if (is_array(tr_x)) {
            newx0 = E.x0 + _trtt_periodic_spline(t, tr_x(1), tr_x(2), tr_x(3));
            h_set, Enew, x0=newx0;
        }
        if (is_array(tr_alpha)) {
            newalpha = E.alpha + _trtt_periodic_spline(t, tr_alpha(1), tr_alpha(2), tr_alpha(3));
            h_set, Enew, alpha=newalpha;
        }
        if (is_array(dis_a)) {
            newa = E.a + _trtt_periodic_spline(t, dis_a(1), dis_a(2), dis_a(3));
            h_set, Enew, a=newa;
        }
        if (is_array(dis_b)) {
            newb = E.b + _trtt_periodic_spline(t, dis_b(1), dis_b(2), dis_b(3));
            h_set, Enew, b=newb;
        }
        
        // Calculate the sinogram for the given ellipse       
        Sinogram += R_ellipse_single_PB(Enew, v, nv, v_scl, theta, epsilon, xoff, yoff);      
    }
        
    return Sinogram;
}
R_ellipses_dyn_PB = trtt_ellipses_dynamic_parallel_beam_analytic_Radon_transform; // Define an easier alias

func trtt_ellipses_dynamic_fan_beam_analytic_Radon_transform(Ellipses, v, nv, v_scl, theta, t, Rsc, Rsd, xoff, yoff, epsilon=)
/* DOCUMENT trtt_ellipses_dynamic_fan_beam_analytic_Radon_transform, Ellipses, v, nv, v_scl, theta, t, Rsc, Rsd, xoff, yoff, epsilon=
   
   This function calculates analytically the Radon transform, in a FAN BEAM
   GEOMETRY, of a 2D tomobj defined as a set of ellipses (for example a
   Shepp-Logan phantom). All the parameters are extracted from the set of
   ellipses ELLIPSES. ELLIPSES is a hash table of the form:

   Ellipses           
   |
   |
   |----- E1
   |       |----- x0     // x-coordinate of the center in meters
   |       |----- y0     // y-coordinate of the center in meters
   |       |----- a      // length of the semimajor axis in meters
   |       |----- b      // length of the semiminor axis in meters
   |       |----- alpha  // angle of rotation of the major axis from
                         // the x-axis in radians
   |       |----- mu     // attenuation (constant on the surface of
                         // the ellipse)
   |
   |
   |----- E2
   |       |----- x0        
   |       |----- y0       
   |       |----- a         
   |       |----- b         
   |       |----- alpha   
   |       |----- mu      
   |
   |
   |----- E3
   :
   :  
   :

   The Radon transform is determined for NV detector pixels at positions V
   sampled at step V_SCL, for an angle of projection THETA. It is calculated on
   band integrals (thick beam). It is necessary to give data compatible with the
   algebraic projectors distance or spline driven.
              
   The function calculates the Radon transform of this set of ellipses which is
   the sum of the Radon transform of each single ellipse Ei:
   
   * Parallel beam : P (theta, v) = SUM_i [ P_Ei (theta, v) ]
   = SUM_i [ L_i (theta, v) * mu_i ]

   where L_i (theta, v) is the length of intersection of the ray in (theta, v)
   and the ellipse Ei, and mu_i is the uniform attenuation on Ei. THETA is the
   projection angle and V is the position on the detector.

   * Fan beam : P (beta, v)
   = SUM_i [ P_parallel_Ei (beta + gamma_v , Rsc * sin(gamma_v)) ]
   
   where RSC is the distance Source - Object Center, gamma_v = tan(v/Rsd), and
   RSD is the distance Source - Detector. BETA is the projection angle and V is
   the position on the detector corresponding to the fan angle gamma_v.

   The function returns the sinogram of the Radon transform.

   SEE ALSO: trtt_ellipses_parallel_beam_analytic_Radon_transform.
*/        
{
    E_list = Ellipses.Ellipses_list;
    nE = numberof(E_list);

    if (is_void(epsilon)) epsilon=1.e-6;
    max_doublings = 50;

    // Create the sinogram
    Sinogram = array(double, nv);
              
    for (i=1; i<=nE; i++)
    {
        // Get the current ellipse
        E = h_get(Ellipses, E_list(i));
        
        tr_x = E.tr_x;
        tr_y = E.tr_y;
        tr_alpha = E.tr_alpha;
        dis_a = E.dis_a;
        dis_b = E.dis_b;

        Enew = h_copy(E);

        if (is_array(tr_y)) {
            newy0 = E.y0 + _trtt_periodic_spline(t, tr_y(1), tr_y(2), tr_y(3));
            h_set, Enew, y0=newy0;
        }
        if (is_array(tr_x)) {
            newx0 = E.x0 + _trtt_periodic_spline(t, tr_x(1), tr_x(2), tr_x(3));
            h_set, Enew, x0=newx0;
        }
        if (is_array(tr_alpha)) {
            newalpha = E.alpha + _trtt_periodic_spline(t, tr_alpha(1), tr_alpha(2), tr_alpha(3));
            h_set, Enew, alpha=newalpha;
        }
        if (is_array(dis_a)) {
            newa = E.a + _trtt_periodic_spline(t, dis_a(1), dis_a(2), dis_a(3));
            h_set, Enew, a=newa;
        }
        if (is_array(dis_b)) {
            newb = E.b + _trtt_periodic_spline(t, dis_b(1), dis_b(2), dis_b(3));
            h_set, Enew, b=newb;
        }
        
        // Calculate the sinogram for the given ellipse       
        Sinogram += R_ellipse_single_FB(Enew, v, nv, v_scl, theta, Rsc, Rsd, epsilon, xoff, yoff);      
    }

    return Sinogram;
}
R_ellipses_dyn_FB = trtt_ellipses_dynamic_fan_beam_analytic_Radon_transform; // Define an easier alias

func disk_dynamic_2D(Nnorm)
{
    DEG2RAD = pi/180;
    Ellipses = h_new(
        E1 = h_new( // bone
            x0 = 0.*Nnorm,  
            y0 = 0.*Nnorm, 
            a = 0.5*Nnorm, 
            b = 0.5*Nnorm,
            alpha = 0.,   
            mu = 1.0,
            dis_a=[0.1*Nnorm, 1./5., -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            dis_b=[0.1*Nnorm, 1./5., -1])
        );
    Ellipses_list = h_keys(Ellipses);
    h_set, Ellipses, Ellipses_list = Ellipses_list;
    
    return Ellipses;
}

func shepp_sparrow_dynamic_2D(Nnorm)
{
    DEG2RAD = pi/180;
    Ellipses = h_new(
        E1 = h_new( // bone
            x0 = 0.*Nnorm,  
            y0 = 0.*Nnorm, 
            a = 0.92*Nnorm, 
            b = 0.69*Nnorm, 
            alpha = pi/2,   
            mu = 1.0,
            dis_a=[0.004*Nnorm, 1./5., -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            dis_b=[0.004*Nnorm, 1./5., -1]), //  [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
        E2 = h_new( // muscle
            x0 = 0.*Nnorm,
            y0 = -0.0184*Nnorm,
            a = 0.874*Nnorm,
            b = 0.6624*Nnorm,
            alpha = pi/2,
            mu = -0.94,
            dis_a=[0.004*Nnorm, 1./5., -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            dis_b=[0.004*Nnorm, 1./5., -1]), //  [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
        E3 = h_new( // fat
            x0 = 0.22*Nnorm,
            y0 = 0.*Nnorm,
            a = 0.31*Nnorm,
            b = 0.11*Nnorm,
            alpha = 72.*DEG2RAD,
            mu = -0.16,
            dis_a=[0.04*Nnorm, 1./5., -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            dis_b=[0.04*Nnorm, 1./5., -1]), //  [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]   
        E4 = h_new( // fat
            x0 = -0.22*Nnorm,
            y0 = 0.*Nnorm,
            a = 0.41*Nnorm,
            b = 0.16*Nnorm,
            alpha = 108.*DEG2RAD,
            mu = -0.16,
            dis_a=[0.04*Nnorm, 1./5., -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            dis_b=[0.04*Nnorm, 1./5., -1]), // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
        E5 = h_new( // insert 1
            x0 = -0.22*Nnorm,
            y0 = 0.*Nnorm,
            a = 0.035*Nnorm,
            b = 0.035*Nnorm,
            alpha = 72.*DEG2RAD,
            mu = 0.1,
            tr_x=[0.016*Nnorm, 1./5., -2], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            tr_y=[0.06*Nnorm, 1./5., -1]), // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
        E6 = h_new( // insert 2
            x0 = 0.22*Nnorm,
            y0 = 0.*Nnorm,
            a = 0.015*Nnorm,
            b = 0.015*Nnorm,
            alpha = 72.*DEG2RAD,
            mu = 1.1,
            tr_x=[-0.008*Nnorm, 1./5., -2], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            tr_y=[0.03*Nnorm, 1./5., -1]), // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
        E7 = h_new( // blood/muscle
            x0 = 0.*Nnorm,
            y0 = 0.35*Nnorm,
            a = 0.25*Nnorm,
            b = 0.21*Nnorm,
            alpha = pi/2,
            mu = -0.01, // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            tr_y=[0.016*Nnorm, 1./5., -1]),
        E8 = h_new( // bone
            x0 = 0.*Nnorm,
            y0 = -0.605*Nnorm,
            a = 0.05*Nnorm,
            b = 0.05*Nnorm,
            alpha = 0.,
            mu = 0.94,
            tr_y=[-0.016*Nnorm, 1./5., -1]), // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
        E9 = h_new( /* FIXME: "Support" ellipse : for absorption scaling */
            x0 = 0.*Nnorm,  
            y0 = 0.*Nnorm, 
            a = 0.92*Nnorm, 
            b = 0.69*Nnorm, 
            alpha = pi/2,   
            mu = 1.0,
            dis_a=[0.004*Nnorm, 1./5., -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            dis_b=[0.004*Nnorm, 1./5., -1]) //  [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
        );
    Ellipses_list = h_keys(Ellipses);
    h_set, Ellipses, Ellipses_list = Ellipses_list;
    
    return Ellipses;
}

func shepp_sparrow_semi_dynamic_2D(Nnorm, tcycl=)
{
    if (is_void(tcycl)) tcycl=5.;
    DEG2RAD = pi/180;
    Ellipses = h_new(
        E1 = h_new( // bone
            x0 = 0.*Nnorm,  
            y0 = 0.*Nnorm, 
            a = 0.92*Nnorm, 
            b = 0.69*Nnorm, 
            alpha = pi/2,   
            mu = 2.0),
        E2 = h_new( // muscle
            x0 = 0.*Nnorm,
            y0 = -0.0184*Nnorm,
            a = 0.874*Nnorm,
            b = 0.6624*Nnorm,
            alpha = pi/2,
            mu = -0.94),
        E3 = h_new( // fat
            x0 = 0.22*Nnorm,
            y0 = 0.*Nnorm,
            a = 0.31*Nnorm,
            b = 0.11*Nnorm,
            alpha = 72.*DEG2RAD,
            mu = -0.16,
            dis_a=[0.04*Nnorm, 1./tcycl, -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            dis_b=[0.04*Nnorm, 1./tcycl, -1]), //  [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]   
        E4 = h_new( // fat
            x0 = -0.22*Nnorm,
            y0 = 0.*Nnorm,
            a = 0.41*Nnorm,
            b = 0.16*Nnorm,
            alpha = 108.*DEG2RAD,
            mu = -0.16,
            dis_a=[0.04*Nnorm, 1./tcycl, -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            dis_b=[0.04*Nnorm, 1./tcycl, -1]), // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
        E5 = h_new( // insert 1
            x0 = -0.22*Nnorm,
            y0 = 0.*Nnorm,
            a = 0.035*Nnorm,
            b = 0.035*Nnorm,
            alpha = 72.*DEG2RAD,
            mu = 0.1,
            tr_x=[0.016*Nnorm, 1./tcycl, -2], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            tr_y=[0.06*Nnorm, 1./tcycl, -1]), // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
        E6 = h_new( // insert 2 a
            x0 = 0.21*Nnorm,
            y0 = 0.*Nnorm,
            a = 0.015*Nnorm,
            b = 0.015*Nnorm,
            alpha = 72.*DEG2RAD,
            mu = 0.1,
            tr_x=[-0.008*Nnorm, 1./tcycl, -2], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            tr_y=[0.03*Nnorm, 1./tcycl, -1]), // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
        E7 = h_new( // insert 2 b
            x0 = 0.23*Nnorm,
            y0 = 0.03*Nnorm,
            a = 0.015*Nnorm,
            b = 0.015*Nnorm,
            alpha = 72.*DEG2RAD,
            mu = 0.15,
            tr_x=[-0.008*Nnorm, 1./tcycl, -2], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            tr_y=[0.03*Nnorm, 1./tcycl, -1]), // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
        E8 = h_new( // blood/muscle
            x0 = 0.*Nnorm,
            y0 = 0.35*Nnorm,
            a = 0.25*Nnorm,
            b = 0.21*Nnorm,
            alpha = pi/2,
            mu = -0.01,
            tr_y=[0.016*Nnorm, 1./tcycl, -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            dis_a=[0.01*Nnorm, 1./1., -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            dis_b=[0.01*Nnorm, 1./1., -1]), // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
        E9 = h_new( // bone
            x0 = 0.*Nnorm,
            y0 = -0.605*Nnorm,
            a = 0.05*Nnorm,
            b = 0.05*Nnorm,
            alpha = 0.,
            mu = 0.94),
        E10 = h_new( // insert 3
            x0 = 0.35*Nnorm,
            y0 = -0.65*Nnorm,
            a = 0.06*Nnorm,
            b = 0.02*Nnorm,
            alpha = pi/4,
            mu = 0.01,
            tr_x=[0.04*Nnorm, 1./tcycl, -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            tr_y=[0.04*Nnorm, 1./tcycl, -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]), // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            tr_alpha=[pi/12, 1./tcycl, -1]),
        E11 = h_new( // insert 4
            x0 = -0.35*Nnorm,
            y0 = -0.65*Nnorm,
            a = 0.06*Nnorm,
            b = 0.02*Nnorm,
            alpha = 3.*pi/4,
            mu = 0.01,
            tr_x=[-0.02*Nnorm, 1./tcycl, -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            tr_y=[0.02*Nnorm, 1./tcycl, -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]), // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            tr_alpha=[-pi/12, 1./tcycl, -1]),
        E12 = h_new( // insert 5
            x0 = 0.3*Nnorm,
            y0 = 0.6*Nnorm,
            a = 0.08*Nnorm,
            b = 0.04*Nnorm,
            alpha = -pi/4,
            mu = 0.025, // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            dis_a=[0.02*Nnorm, 1./tcycl, -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            dis_b=[0.01*Nnorm, 1./tcycl, -1]), // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
        E13 = h_new( // insert 6 a
            x0 = -0.55*Nnorm,
            y0 = 0.*Nnorm,
            a = 0.06*Nnorm,
            b = 0.02*Nnorm,
            alpha = pi/3,
            mu = 0.025, // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            tr_x=[0.02*Nnorm, 1./tcycl, -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            tr_y=[0.04*Nnorm, 1./tcycl, -1]), // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
        E14 = h_new( // insert 6 b
            x0 = -0.5*Nnorm,
            y0 = -0.05*Nnorm,
            a = 0.06*Nnorm,
            b = 0.02*Nnorm,
            alpha = pi/3,
            mu = -0.025, // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            tr_x=[-0.02*Nnorm, 1./tcycl, -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            tr_y=[-0.04*Nnorm, 1./tcycl, -1]), // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
        E15 = h_new( // insert 7
            x0 = 0.45*Nnorm,
            y0 = -0.25*Nnorm,
            a = 0.1*Nnorm,
            b = 0.06*Nnorm,
            alpha = pi/6,
            mu = -0.01, // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            dis_a=[-0.02*Nnorm, 1./tcycl, -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            dis_b=[-0.01*Nnorm, 1./tcycl, -1]) // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
        // E16 = h_new( /* FIXME: "Support" ellipse : for absorption scaling */
        //     x0 = 0.*Nnorm,  
        //     y0 = 0.*Nnorm, 
        //     a = 0.92*Nnorm, 
        //     b = 0.69*Nnorm, 
        //     alpha = pi/2,   
        //     mu = 1.0)
            );
    Ellipses_list = h_keys(Ellipses);
    h_set, Ellipses, Ellipses_list = Ellipses_list;
    
    return Ellipses;
}
