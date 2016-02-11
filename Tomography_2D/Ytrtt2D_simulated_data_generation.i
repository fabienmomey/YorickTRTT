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

/* USER FUNCTIONS ============================================================ */

func trtt_create_ellipses_phantom_evolve(Ellipses, nx, ny, s_scl, xoff, yoff, shape=, get_masks=) 
{ 
    if (is_void(shape)) shape = 0n;
    h = abs(shape)/2.0; /* half border thickness */
    
    E_list = Ellipses.Ellipses_list;
    nE = numberof(E_list);
    
    // Create the phantom image
    Im = array(double, nx, ny);

    for (k=1; k<=nE; ++k)
    {
        // Get the current ellipse
        E = h_get(Ellipses, E_list(k));          
        // Extract its geometrical parameters on all the frames (3rd
        // dimension)
        x0 = E.x0 + xoff;
        y0 = E.y0 + yoff;
        a = E.a /*- 0.5*dx*/;
        b = E.b /*- 0.5*dy*/;
        alpha = E.alpha;
        mu = E.mu*TRTT_WATER_ABSORPTION; /* FIXME: conversion to true absorption */

        sqa=a*a;
        sqb=b*b;
        
        /* distance map from ellipse center */
        x = (indgen(nx)-((nx+1.0)*0.5))*s_scl-x0;
        y = (indgen(ny)-((ny+1.0)*0.5))*s_scl-y0;
        u = x*cos(alpha)+y(-,)*sin(alpha);
        v = -1.0*x*sin(alpha)+y(-,)*cos(alpha);
        r = abs(u,v);
        theta = atan(v,u);
        rho=sqrt((sqa*sqb)/(sqb*(cos(theta)^2.0)+sqa*(sin(theta)^2.0)));
        id_rho=where(r==0.0);
        if(is_array(id_rho)) rho(id_rho)=0.0;

        /* All points strictly inside the pupil are set to 1. */
        tmax = (r <= rho - h); /* strictly inside */
        p = double(tmax);
        local outer;
        outer = where((! tmax) & (r < rho + h)); /* at the outer edge */
        a = 1.0/(h + h);
        if (is_array(outer)) {
            p(outer) = a*((rho(outer) + h) - r(outer));
        }

        Im += mu*p;

        /* Masks asked */
        if (get_masks=1n) {
            h_set, E, mask_inside=tmax;
            h_set, E, mask=(r < rho+h);
        }
    }

    h_set, Ellipses, phantom=Im;
    
    return Im;
}
trtt_create_phantom_evolve = trtt_create_ellipses_phantom_evolve;

func trtt_create_ellipses_phantom(Ellipses, nx, ny, s_scl, xoff, yoff, shape=, oversampling=) 
/* DOCUMENT trtt_create_ellipses_phantom_rigid_motion, Ellipses, nx, ny, s_scl, xoff, yoff, Ellipses, oversampling=
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
    // Oversampling in one direction
    if (is_void(oversampling)) oversampling = 0n; 
    if (oversampling%2) {
        msg_to_user, "The oversampling value must be even.";
    }
    ns = oversampling;

    if (is_void(shape)) shape = 0n;
    h = abs(shape)/2.0; /* half border thickness */
    
    E_list = Ellipses.Ellipses_list;
    nE = numberof(E_list);

    /* Calculate coordinates */
    x = (indgen(nx)-((nx+1.0)*0.5))*s_scl;
    y = (indgen(ny)-((ny+1.0)*0.5))*s_scl;

    // Oversampling step
    if (ns > 0) {
        x_spl = span(x + 0.5*((1.0-ns)/ns), x + 0.5*((ns-1.0)/ns), ns);
        y_spl = span(y + 0.5*((1.0-ns)/ns), y + 0.5*((ns-1.0)/ns), ns);
                
        x_spl = x_spl(sort(x_spl(*)));
        y_spl = y_spl(sort(y_spl(*)));

        x_spl = x_spl(,-);
        y_spl = y_spl(-,);

        // Create the phantom image
        Im = array(double, ns*nx, ns*ny);
    } else {
        x_spl = x(,-);
        y_spl = y(-,);

        // Create the phantom image
        Im = array(double, nx, ny);
    }

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

        // if (!shape) {
        // Calculate the coefficients of the equation of the ellipse
        A = (cos(alpha) / a)^2 + (sin(alpha) / b)^2;
        B = (sin(alpha) / a)^2 + (cos(alpha) / b)^2;
        C = 2 * cos(alpha) * sin(alpha) * ((1/a^2) - (1/b^2));
        // Resolve the equation by identifying the pixels included in the ellipse
        eq = (A*(x_spl - x0)^2 + B*(y_spl - y0)^2 + C*(x_spl-x0)*(y_spl-y0) <= 1);
        // Add the result to the phantom image with the attenuation value
        Im += mu * eq;
        // } else {
        //     A = (cos(alpha) / a)^2 + (sin(alpha) / b)^2;
        //     B = (sin(alpha) / a)^2 + (cos(alpha) / b)^2;
        //     C = 2 * cos(alpha) * sin(alpha) * ((1/a^2) - (1/b^2));
        //     r = sqrt(A*(x_spl - x0)^2 + B*(y_spl - y0)^2 + C*(x_spl-x0)*(y_spl-y0));
        //     /* All points strictly inside the pupil are set to 1. */
        //     ah=a-h;
        //     bh=b-h;
        //     A_in = (cos(alpha) / ah)^2 + (sin(alpha) / bh)^2;
        //     B_in = (sin(alpha) / ah)^2 + (cos(alpha) / bh)^2;
        //     C_in = 2 * cos(alpha) * sin(alpha) * ((1/ah^2) - (1/bh^2));
        //     rh_moins = sqrt(A_in*(x_spl - x0)^2 + B_in*(y_spl - y0)^2 + C_in*(x_spl-x0)*(y_spl-y0));
        //     eq_in = (rh_moins*rh_moins < 1);
        //     /* All points strictly outside the pupil are set to 1. */
        //     ah=a+h;
        //     bh=b+h;
        //     A_out = (cos(alpha) / ah)^2 + (sin(alpha) / bh)^2;
        //     B_out = (sin(alpha) / ah)^2 + (cos(alpha) / bh)^2;
        //     C_out = 2 * cos(alpha) * sin(alpha) * ((1/ah^2) - (1/bh^2));
        //     rh_plus = sqrt(A_out*(x_spl - x0)^2 + B_out*(y_spl - y0)^2 + C_out*(x_spl-x0)*(y_spl-y0));
        //     eq_outer = (rh_plus*rh_plus < 1);
           
        //     local outer;
            
        //     outer = where((! eq_in) & (eq_outer)); /* at the outer edge */
        //     eq_in=double(eq_in);
        //     /* linear B-spline */
        //     h_eq = avg(0.5*(rh_moins(outer)-rh_plus(outer)));
        //     a = 1.0/(h_eq+h_eq);
        //     eq_in(outer) = a*((1+h_eq) - r(outer));
        // }

        // Im += mu*eq_in;
    }
    
    if (ns > 0) {
        Phantom = array(double, nx, ny);

        for (u=1; u<=nx; ++u) {
            i = (u-1) * ns + 1;
            for (v=1; v<=ny; ++v) {
                j = (v-1) * ns + 1;
                Phantom(u,v) = avg(Im(i:i+ns-1,j:j+ns-1));
            }
        }
        return Phantom;
    } else {
        return Im;
    }
}
trtt_create_phantom = trtt_create_ellipses_phantom;

func trtt_get_ellipse_parallel_beam_analytic_Radon_profile(param, v)
/* DOCUMENT trtt_get_ellipse_parallel_beam_analytic_Radon_profile, param, v

   This function calculates the Radon profile of an ellipse for a parallel beam
   geometry, given its geometrical parameters in the hash table PARAM, of the
   form:
   
   PARAM
    |----- theta  // Angle of projection
    |----- x0     // x-coordinate of the center in meters
    |----- y0     // y-coordinate of the center in meters
    |----- a      // length of the semimajor axis in meters
    |----- b      // length of the semiminor axis in meters
    |----- alpha  // angle of rotation of the major axis from
                  // the x-axis in radians
    |----- mu     // attenuation (constant on the surface of
                  // the ellipse)

   The analytic expression of the Radon profile of an ellipse with parameters
   given above is :

                                    sqrt(rho²-(v-v0)²) 
   P (theta, v)] = 2 * mu * ab *  ------------------
                                           rho²

   Where rho² = a² cos²(theta-alpha) + b² sin²(theta-alpha)
   
   and v0 = sqrt(x0²+y0²) * cos(tan⁻¹(y0/x0)-theta)
    
 */
{
    theta = param.theta;
    a = param.a;
    b = param.b;
    x0 = param.x0;
    y0 = param.y0;
    alpha = param.alpha;

    gam = atan(y0,x0);
    
    theta_prime = theta-alpha;
    rho = (a*cos(theta_prime))^2+(b*sin(theta_prime))^2;
    v0 = sqrt(x0^2+y0^2) * cos(gam-theta);
    v_square = min(rho, (v-v0)^2);

    dim = numberof(v);
    profile = (2.0/rho)*a*b*sqrt(rho-v_square);
    if (dim==1) profile = profile(1);
    
    return profile;
}

func trtt_get_ellipse_fan_beam_analytic_Radon_profile(param, v)
/* DOCUMENT trtt_get_ellipse_fan_beam_analytic_Radon_profile, param, v

   This function calculates the Radon profile of an ellipse for a fan beam
   geometry, given its geometrical parameters in the hash table PARAM, of the
   form:
   
   PARAM
    |----- theta  // Angle of projection
    |----- x0     // x-coordinate of the center in meters
    |----- y0     // y-coordinate of the center in meters
    |----- a      // length of the semimajor axis in meters
    |----- b      // length of the semiminor axis in meters
    |----- alpha  // angle of rotation of the major axis from
                  // the x-axis in radians
    |----- mu     // attenuation (constant on the surface of
                  // the ellipse)
    |----- Rsc    // Source Center distance
    |----- Rsd    // Source Detector distance

   The analytic expression of the Radon profile of an ellipse with parameters
   given above is :

                                    sqrt(rho²-(v_eq-v0)²) 
   P (theta, v)] = 2 * mu * ab *  ------------------
                                           rho²

   Where rho² = a² cos²(theta-alpha) + b² sin²(theta-alpha)
   
   and v0 = sqrt(x0²+y0²) * cos(tan⁻¹(y0/x0)-theta)

   For the fan beam geometry, the following change of variable is done:
   
   /
   | gamma_v = atan(v/Rsd);
   | theta = beta + gamma_s
   | v_eq = Rsc * sin(gamma_v))
   \
    
 */
{
    theta = param.theta;
    a = param.a;
    b = param.b;
    x0 = param.x0;
    y0 = param.y0;
    alpha = param.alpha;
    Rsc = param.Rsc;
    Rsd = param.Rsd;

    gam = atan(y0,x0);

    gamma_v = atan(v/Rsd);
    beta = theta + gamma_v;
    
    theta_prime = beta-alpha;
    v_eq = Rsc * sin(gamma_v);
    
    rho = (a*cos(theta_prime))^2+(b*sin(theta_prime))^2;
    v0 = sqrt(x0^2+y0^2) * cos(gam-beta);
    v_square = min(rho, (v_eq-v0)^2);

    dim = numberof(v);
    profile = (2.0/rho)*a*b*sqrt(rho-v_square);
    if (dim==1) profile = profile(1);

    // tan_gamma = (v/Rsd);
    // cos_gamma = cos(atan(tan_gamma));
    // cos_th = cos(theta);
    // sin_th = sin(theta);
    // cos_alpha = cos(alpha);
    // sin_alpha = sin(alpha);

    // A = (cos_alpha/a)^2.0 + (sin_alpha/b)^2.0;
    // B = (sin_alpha/a)^2.0 + (cos_alpha/b)^2.0;
    // C = 2.0*cos_alpha*sin_alpha*((1.0/(a^2.0)) - (1.0/(b^2.0)));
    // D = 2.0*A*x0 + C*y0;
    // E = 2.0*B*x0 + C*y0;
    // X = (cos_th+tan_gamma*sin_th)*cos_gamma;
    // Y = (-1.0*sin_th+tan_gamma*cos_th)*cos_gamma;

    // d = A*X^2.0 + B*Y^2.0 + C*X*Y;
    // e = -1.0 * (D*X + E*Y);
    // f = A*x0^2.0 + B*y0^2.0 + C*x0*y0 -1.0;

    // delta = e^2.0 - 4.0*d*f;

    // dim = numberof(v);
    // profile = array(double, dim);
    // idelta = where(delta>0.0);
    // if (is_array(idelta)) profile(idelta) = abs(sqrt(delta(idelta))/d(idelta));
    // if (dim==1) profile = profile(1);
    
    return profile;
}

func trtt_ellipse_parallel_beam_analytic_Radon_profile_evaluator(theta,a,b,x0,y0,alpha)
/* DOCUMENT trtt_ellipse_parallel_beam_analytic_Radon_profile_evaluator theta, a, b, x0, y0, alpha

   Function creating the EVALUATOR of the Radon profile of an ellipse for a
   parallel beam geometry, given its geometrical parameters :
   
   PARAM
    |----- theta  // Angle of projection
    |----- x0     // x-coordinate of the center in meters
    |----- y0     // y-coordinate of the center in meters
    |----- a      // length of the semimajor axis in meters
    |----- b      // length of the semiminor axis in meters
    |----- alpha  // angle of rotation of the major axis from
                  // the x-axis in radians
    |----- mu     // attenuation (constant on the surface of
                  // the ellipse)
 */
{
    obj = h_new(theta = theta,
                a = a,
                b = b,
                x0 = x0,
                y0 = y0,
                alpha = alpha
        );
    h_evaluator, obj, trtt_get_ellipse_parallel_beam_analytic_Radon_profile;
    return obj;
}

func trtt_ellipse_fan_beam_analytic_Radon_profile_evaluator(theta,a,b,x0,y0,alpha,Rsc,Rsd)
/* DOCUMENT trtt_ellipse_parallel_beam_analytic_Radon_profile_evaluator theta, a, b, x0, y0, alpha, Rsc, Rsd

   Function creating the EVALUATOR of the Radon profile of an ellipse for a fan
   beam geometry, given its geometrical parameters :
   
   PARAM
    |----- theta  // Angle of projection
    |----- x0     // x-coordinate of the center in meters
    |----- y0     // y-coordinate of the center in meters
    |----- a      // length of the semimajor axis in meters
    |----- b      // length of the semiminor axis in meters
    |----- alpha  // angle of rotation of the major axis from
                  // the x-axis in radians
    |----- mu     // attenuation (constant on the surface of
                  // the ellipse)
    |----- Rsc    // Source Center distance
    |----- Rsd    // Source Detector distance
 */
{
    obj = h_new(theta = theta,
                a = a,
                b = b,
                x0 = x0,
                y0 = y0,
                alpha = alpha,
                Rsc = Rsc,
                Rsd = Rsd
        );
    h_evaluator, obj, trtt_get_ellipse_fan_beam_analytic_Radon_profile;
    return obj;
}

func trtt_single_ellipse_parallel_beam_analytic_Radon_transform(Ellipse, v, nv, v_scl, theta, epsilon, xoff, yoff)
/* DOCUMENT trtt_single_ellipse_parallel_beam_analytic_Radon_transform, Ellipse, v, nv, v_scl, theta, xoff, yoff

   This function calculates analytically the Radon transform, in PARALLEL BEAM
   geometry, of an ellipse given a set of parameters defined in an ELLIPSE hash
   table :

   Ellipse
   |----- x0      // x-coordinate of the center in meters
   |----- y0      // y-coordinate of the center in meters
   |----- a       // length of the semimajor axis in meters
   |----- b       // length of the semiminor axis in meters
   |----- alpha   // angle of rotation of the major axis from the x-axis
   // in radians
   |----- mu      // attenuation (constant on the surface of the ellipse)

   The Radon transform is calculated on NV detector pixels at positions V, for
   an angle of projection THETA.

   The beam passing through the object is considered as having a thickness the
   width of which is the same as a pixel detector size V_SCL. It gives data
   compatible with a distance or spline driven projector.

   The algorithm looks for the length of intersection between the ellipse and the
   beam passing through the object. The analytic expression of the Radon profile
   of an ellipse with parameters given above is :

                                 sqrt(rho²-(v-v0)²) 
   R (theta, v) = 2 * mu * ab *  ------------------
                                        rho²

   Where rho² = a² cos²(theta-alpha) + b² sin²(theta-alpha)
   
   and v0 = sqrt(x0²+y0²) * cos(tan⁻¹(y0/x0)-theta)

   To consider a thick beam, the Radon profile is integrated on the length of
   each pixel detector, using the function ROMBERG of YORICK with precision
   EPSILON.
   
   The function returns the sinogram of the Radon transform.

   SEE ALSO: trtt_single_ellipse_fan_beam_analytic_Radon_transform.
*/
{
    extern max_doublings;
    
    // Extract the geometrical parameters of the Ellipse
    x0 = Ellipse.x0 + xoff;
    y0 = Ellipse.y0 + yoff;
    a = Ellipse.a;
    b = Ellipse.b;
    alpha = Ellipse.alpha;
    mu = Ellipse.mu;

    E_sino = array(double, nv);
    
    for (q=1; q<=nv; ++q) {   
        E_sino_eval = trtt_ellipse_parallel_beam_analytic_Radon_profile_evaluator(theta,a,b,x0,y0,alpha);   
        E_sino(q) = mu * TRTT_WATER_ABSORPTION * romberg(E_sino_eval,v(q)-0.5*v_scl,v(q)+0.5*v_scl,epsilon);
    }    
    return E_sino;  
}
R_ellipse_single_PB = trtt_single_ellipse_parallel_beam_analytic_Radon_transform; // Define an easier alias

func trtt_single_ellipse_fan_beam_analytic_Radon_transform(Ellipse, v, nv, v_scl, theta, Rsc, Rsd, epsilon, xoff, yoff)
/* DOCUMENT trtt_single_ellipse_fan_beam_analytic_Radon_transform, Ellipse, v, nv, v_scl, theta, Rsc, Rsd, xoff, yoff

   This function calculates analytically the Radon transform, in FAN BEAM
   geometry, of an ellipse given a set of parameters defined in an ELLIPSE hash
   table :

   Ellipse
   |----- x0      // x-coordinate of the center in meters
   |----- y0      // y-coordinate of the center in meters
   |----- a       // length of the semimajor axis in meters
   |----- b       // length of the semiminor axis in meters
   |----- alpha   // angle of rotation of the major axis from the x-axis
   // in radians
   |----- mu      // attenuation (constant on the surface of the ellipse)

   The Radon transform is calculated on NV detector pixels at positions V, for
   an angle of projection THETA.

   The beam passing through the object is considered as having a thickness the
   width of which is the same as a pixel detector size V_SCL. It gives data
   compatible with a distance or spline driven projector.

   The algorithm looks for the length of intersection between the ellipse and the
   beam passing through the object. The analytic expression of the Radon profile
   of an ellipse with parameters given above is :

                                 sqrt(rho²-(v-v0)²) 
   R (theta, v) = 2 * mu * ab *  ------------------
                                        rho²

   Where rho² = a² cos²(theta-alpha) + b² sin²(theta-alpha)
   
   and v0 = sqrt(x0²+y0²) * cos(tan⁻¹(y0/x0)-theta)

   This result is for a "parallel beam" projection mode. For a "fan beam"
   projection, a change of variable is necessary:

   R (beta, v) = R_parallel (beta + gamma_v, Rsc * sin(gamma_v)) ]

   BETA is the projection angle and V is the position on the detector
   corresponding to the fan angle gamma_v.

   To consider a thick beam, the Radon profile is integrated on the length of
   each pixel detector, using the function ROMBERG of YORICK with precision
   EPSILON.
   
   The function returns the sinogram of the Radon transform.

   SEE ALSO: trtt_single_ellipse_parallel_beam_analytic_Radon_transform.
*/
{
    extern max_doublings;
    
    // Extract the geometrical parameters of the Ellipse
    x0 = Ellipse.x0 + xoff;
    y0 = Ellipse.y0 + yoff;
    a = Ellipse.a;
    b = Ellipse.b;
    alpha = Ellipse.alpha;
    mu = Ellipse.mu;

    E_sino = array(double, nv);
    
    for (q=1; q<=nv; ++q) {
        // factor = double(Rsc)/Rsd * abs(1.0/cos(atan(v(q),Rsd)));                
        E_sino_eval = trtt_ellipse_fan_beam_analytic_Radon_profile_evaluator(theta,a,b,x0,y0,alpha,Rsc,Rsd); 
        E_sino(q) = mu * TRTT_WATER_ABSORPTION * romberg(E_sino_eval,v(q)-0.5*v_scl,v(q)+0.5*v_scl,epsilon);
        // E_sino(q) = trtt_get_ellipse_fan_beam_analytic_Radon_profile(E_sino_eval, v(q))
    }       
    return E_sino;  
}
R_ellipse_single_FB = trtt_single_ellipse_fan_beam_analytic_Radon_transform; // Define an easier alias

func trtt_ellipses_parallel_beam_analytic_Radon_transform(Ellipses, v, nv, v_scl, theta, xoff, yoff, epsilon=)
/* DOCUMENT trtt_ellipses_parallel_beam_analytic_Radon_transform, Ellipses, v, nv, v_scl, theta, xoff, yoff, epsilon=
   
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
        // Calculate the sinogram for the given ellipse       
        Sinogram += R_ellipse_single_PB(E, v, nv, v_scl, theta, epsilon, xoff, yoff);      
    }
        
    return Sinogram;
}
R_ellipses_PB = trtt_ellipses_parallel_beam_analytic_Radon_transform; // Define an easier alias

func trtt_ellipses_fan_beam_analytic_Radon_transform(Ellipses, v, nv, v_scl, theta, Rsc, Rsd, xoff, yoff, epsilon=)
/* DOCUMENT trtt_ellipses_fan_beam_analytic_Radon_transform, Ellipses, v, nv, v_scl, theta, Rsc, Rsd, xoff, yoff, epsilon=
   
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
        // Calculate the sinogram for the given ellipse       
        Sinogram += R_ellipse_single_FB(E, v, nv, v_scl, theta, Rsc, Rsd, epsilon, xoff, yoff);      
    }

    return Sinogram;
}
R_ellipses_FB = trtt_ellipses_fan_beam_analytic_Radon_transform; // Define an easier alias

func shepp_sparrow_semi_dynamic_2D(Nnorm)
{
    DEG2RAD = pi/180;
    Ellipses = h_new(
        E1 = h_new( // bone
            x0 = 0.*Nnorm,  
            y0 = 0.*Nnorm, 
            a = 0.92*Nnorm, 
            b = 0.69*Nnorm, 
            alpha = pi/2,   
            mu = 1.0),
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
        E6 = h_new( // insert 2 a
            x0 = 0.21*Nnorm,
            y0 = 0.*Nnorm,
            a = 0.015*Nnorm,
            b = 0.015*Nnorm,
            alpha = 72.*DEG2RAD,
            mu = 0.1,
            tr_x=[-0.008*Nnorm, 1./5., -2], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            tr_y=[0.03*Nnorm, 1./5., -1]), // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
        E7 = h_new( // insert 2 b
            x0 = 0.23*Nnorm,
            y0 = 0.03*Nnorm,
            a = 0.015*Nnorm,
            b = 0.015*Nnorm,
            alpha = 72.*DEG2RAD,
            mu = 0.15,
            tr_x=[-0.008*Nnorm, 1./5., -2], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            tr_y=[0.03*Nnorm, 1./5., -1]), // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
        E8 = h_new( // blood/muscle
            x0 = 0.*Nnorm,
            y0 = 0.35*Nnorm,
            a = 0.25*Nnorm,
            b = 0.21*Nnorm,
            alpha = pi/2,
            mu = -0.01,
            tr_y=[0.016*Nnorm, 1./5., -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
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
            tr_x=[0.04*Nnorm, 1./5., -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            tr_y=[0.04*Nnorm, 1./5., -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]), // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            tr_alpha=[pi/12, 1./5., -1]),
        E11 = h_new( // insert 4
            x0 = -0.35*Nnorm,
            y0 = -0.65*Nnorm,
            a = 0.06*Nnorm,
            b = 0.02*Nnorm,
            alpha = 3.*pi/4,
            mu = 0.01,
            tr_x=[-0.02*Nnorm, 1./5., -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            tr_y=[0.02*Nnorm, 1./5., -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]), // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            tr_alpha=[-pi/12, 1./5., -1]),
        E12 = h_new( // insert 5
            x0 = 0.3*Nnorm,
            y0 = 0.6*Nnorm,
            a = 0.08*Nnorm,
            b = 0.04*Nnorm,
            alpha = -pi/4,
            mu = 0.025, // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            dis_a=[0.02*Nnorm, 1./5., -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            dis_b=[0.01*Nnorm, 1./5., -1]), // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
        E13 = h_new( // insert 6 a
            x0 = -0.55*Nnorm,
            y0 = 0.*Nnorm,
            a = 0.06*Nnorm,
            b = 0.02*Nnorm,
            alpha = pi/3,
            mu = 0.025, // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            tr_x=[0.02*Nnorm, 1./5., -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            tr_y=[0.04*Nnorm, 1./5., -1]), // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
        E14 = h_new( // insert 6 b
            x0 = -0.5*Nnorm,
            y0 = -0.05*Nnorm,
            a = 0.06*Nnorm,
            b = 0.02*Nnorm,
            alpha = pi/3,
            mu = -0.025, // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            tr_x=[-0.02*Nnorm, 1./5., -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            tr_y=[-0.04*Nnorm, 1./5., -1]), // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
        E15 = h_new( // insert 7
            x0 = 0.45*Nnorm,
            y0 = -0.25*Nnorm,
            a = 0.1*Nnorm,
            b = 0.06*Nnorm,
            alpha = pi/6,
            mu = -0.01, // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            dis_a=[-0.02*Nnorm, 1./5., -1], // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
            dis_b=[-0.01*Nnorm, 1./5., -1]), // [amplitude (in cm), frequency (in s^{-1}), degree periodic spline (-1 to have a sinus)]
        E16 = h_new( /* FIXME: "Support" ellipse : for absorption scaling */
            x0 = 0.*Nnorm,  
            y0 = 0.*Nnorm, 
            a = 0.92*Nnorm, 
            b = 0.69*Nnorm, 
            alpha = pi/2,   
            mu = 1.0)
            );
    Ellipses_list = h_keys(Ellipses);
    h_set, Ellipses, Ellipses_list = Ellipses_list;
    
    return Ellipses;
}

func shepp_sparrow_2D(Nnorm)
{
    DEG2RAD = pi/180;
    Ellipses = h_new(
        E1 = h_new( // bone
            x0 = 0.*Nnorm,  
            y0 = 0.*Nnorm, 
            a = 0.92*Nnorm, 
            b = 0.69*Nnorm, 
            alpha = pi/2,   
            mu = 1.0),
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
            mu = -0.16),              
        E4 = h_new( // blood/muscle
            x0 = -0.22*Nnorm,
            y0 = 0.*Nnorm,
            a = 0.41*Nnorm,
            b = 0.16*Nnorm,
            alpha = 108.*DEG2RAD,
            mu = -0.02),
        E5 = h_new( // blood/muscle
            x0 = 0.*Nnorm,
            y0 = 0.35*Nnorm,
            a = 0.25*Nnorm,
            b = 0.21*Nnorm,
            alpha = pi/2,
            mu = -0.01),
        E6 = h_new( // bone
            x0 = 0.*Nnorm,
            y0 = 0.1*Nnorm,
            a = 0.046*Nnorm,
            b = 0.046*Nnorm,
            alpha = 0.,
            mu = 0.2),
        E7 = h_new( // ~fat/muscle
            x0 = 0.*Nnorm,
            y0 = -0.1*Nnorm,
            a = 0.046*Nnorm,
            b = 0.046*Nnorm,
            alpha = 0.,
            mu = -0.06),
        E8 = h_new( // fat
            x0 = -0.08*Nnorm,
            y0 = -0.605*Nnorm,
            a = 0.046*Nnorm,
            b = 0.023*Nnorm,
            alpha = 0.,
            mu = -0.16),
        E9 = h_new( // bone
            x0 = 0.*Nnorm,
            y0 = -0.605*Nnorm,
            a = 0.023*Nnorm,
            b = 0.023*Nnorm,
            alpha = 0.,
            mu = 0.94),
        E10 = h_new( // blood/muscle
            x0 = 0.06*Nnorm,
            y0 = -0.605*Nnorm,
            a = 0.046*Nnorm,
            b = 0.023*Nnorm,
            alpha = pi/2,
            mu = -0.02),
        E11 = h_new(/*FIXME: "Support" ellipse : for absorption scaling */
            x0 = 0.*Nnorm,  
            y0 = 0.*Nnorm, 
            a = 0.92*Nnorm, 
            b = 0.69*Nnorm, 
            alpha = pi/2,   
            mu = 1.0)
        );
    Ellipses_list = h_keys(Ellipses);
    h_set, Ellipses, Ellipses_list = Ellipses_list;
    
    return Ellipses;
}


func shepp_logan_2D(Nnorm)
{
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
            mu = -0.98),
        E3 = h_new( // fat
            x0 = 0.22*Nnorm,
            y0 = 0.*Nnorm,
            a = 0.31*Nnorm,
            b = 0.11*Nnorm,
            alpha = 72.*DEG2RAD,
            mu = -0.02),              
        E4 = h_new( // blood/muscle
            x0 = -0.22*Nnorm,
            y0 = 0.*Nnorm,
            a = 0.41*Nnorm,
            b = 0.16*Nnorm,
            alpha = 108.*DEG2RAD,
            mu = -0.02),
        E5 = h_new( // blood/muscle
            x0 = 0.*Nnorm,
            y0 = 0.35*Nnorm,
            a = 0.25*Nnorm,
            b = 0.21*Nnorm,
            alpha = pi/2,
            mu = 0.01),
        E6 = h_new( // bone
            x0 = 0.*Nnorm,
            y0 = 0.1*Nnorm,
            a = 0.046*Nnorm,
            b = 0.046*Nnorm,
            alpha = 0.,
            mu = 0.01),
        E7 = h_new( // ~fat/muscle
            x0 = 0.*Nnorm,
            y0 = -0.1*Nnorm,
            a = 0.046*Nnorm,
            b = 0.046*Nnorm,
            alpha = 0.,
            mu = 0.01),
        E8 = h_new( // fat
            x0 = -0.08*Nnorm,
            y0 = -0.605*Nnorm,
            a = 0.046*Nnorm,
            b = 0.023*Nnorm,
            alpha = 0.,
            mu = 0.01),
        E9 = h_new( // bone
            x0 = 0.*Nnorm,
            y0 = -0.605*Nnorm,
            a = 0.023*Nnorm,
            b = 0.023*Nnorm,
            alpha = 0.,
            mu = 0.01),
        E10 = h_new( // blood/muscle
            x0 = 0.06*Nnorm,
            y0 = -0.605*Nnorm,
            a = 0.046*Nnorm,
            b = 0.023*Nnorm,
            alpha = pi/2,
            mu = 0.01)// ,
        // E11 = h_new(/*FIXME: "Support" ellipse : for absorption scaling */
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

func shepp_logan_2D_Jan(Nnorm)
{
    DEG2RAD = pi/180;
    Ellipses = h_new(
        E1 = h_new( // bone
            y0 = 0.*Nnorm,  
            x0 = 0.*Nnorm, 
            b = 0.92*Nnorm, 
            a = 0.69*Nnorm, 
            alpha = -pi/2,   
            mu = 2.0),
        E2 = h_new( // muscle
            y0 = 0.*Nnorm,
            x0 = -0.0184*Nnorm,
            b = 0.874*Nnorm,
            a = 0.6624*Nnorm,
            alpha = pi/2,
            mu = -0.98),
        E3 = h_new( // fat
            y0 = 0.22*Nnorm,
            x0 = 0.*Nnorm,
            b = 0.31*Nnorm,
            a = 0.11*Nnorm,
            alpha = -72.*DEG2RAD,
            mu = -0.02),              
        E4 = h_new( // blood/muscle
            y0 = -0.22*Nnorm,
            x0 = 0.*Nnorm,
            b = 0.41*Nnorm,
            a = 0.16*Nnorm,
            alpha = -108.*DEG2RAD,
            mu = -0.02),
        E5 = h_new( // blood/muscle
            y0 = 0.*Nnorm,
            x0 = 0.35*Nnorm,
            b = 0.25*Nnorm,
            a = 0.21*Nnorm,
            alpha = -pi/2,
            mu = 0.01),
        E6 = h_new( // bone
            y0 = 0.*Nnorm,
            x0 = 0.1*Nnorm,
            b = 0.046*Nnorm,
            a = 0.046*Nnorm,
            alpha = 0.,
            mu = 0.01),
        E7 = h_new( // ~fat/muscle
            y0 = 0.*Nnorm,
            x0 = -0.1*Nnorm,
            b = 0.046*Nnorm,
            a = 0.046*Nnorm,
            alpha = 0.,
            mu = 0.01),
        E8 = h_new( // fat
            y0 = -0.08*Nnorm,
            x0 = -0.605*Nnorm,
            b = 0.046*Nnorm,
            a = 0.023*Nnorm,
            alpha = 0.,
            mu = 0.01),
        E9 = h_new( // bone
            y0 = 0.*Nnorm,
            x0 = -0.605*Nnorm,
            b = 0.023*Nnorm,
            a = 0.023*Nnorm,
            alpha = 0.,
            mu = 0.01),
        E10 = h_new( // blood/muscle
            y0 = 0.06*Nnorm,
            x0 = -0.605*Nnorm,
            b = 0.046*Nnorm,
            a = 0.023*Nnorm,
            alpha = -pi/2,
            mu = 0.01)// ,
        // E11 = h_new(/*FIXME: "Support" ellipse : for absorption scaling */
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

func shepp_sparrow_semi_dynamic_2D_first_frame(Nnorm)
{
    DEG2RAD = pi/180;
    Ellipses = h_new(
        E1 = h_new( // bone
            x0 = 0.*Nnorm,  
            y0 = 0.*Nnorm, 
            a = 0.92*Nnorm, 
            b = 0.69*Nnorm, 
            alpha = pi/2,   
            mu = 1.0),
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
            mu = -0.16), 
        E4 = h_new( // fat
            x0 = -0.22*Nnorm,
            y0 = 0.*Nnorm,
            a = 0.41*Nnorm,
            b = 0.16*Nnorm,
            alpha = 108.*DEG2RAD,
            mu = -0.16),
        E5 = h_new( // insert
            x0 = -0.22*Nnorm,
            y0 = 0.*Nnorm,
            a = 0.035*Nnorm,
            b = 0.035*Nnorm,
            alpha = 72.*DEG2RAD,
            mu = 0.1),
        E6 = h_new( // blood/muscle
            x0 = 0.*Nnorm,
            y0 = 0.35*Nnorm,
            a = 0.25*Nnorm,
            b = 0.21*Nnorm,
            alpha = pi/2,
            mu = -0.01),
        E7 = h_new( // bone
            x0 = 0.*Nnorm,
            y0 = -0.605*Nnorm,
            a = 0.05*Nnorm,
            b = 0.05*Nnorm,
            alpha = 0.,
            mu = 0.94),
        E8 = h_new( /* FIXME: "Support" ellipse : for absorption scaling */
            x0 = 0.*Nnorm,  
            y0 = 0.*Nnorm, 
            a = 0.92*Nnorm, 
            b = 0.69*Nnorm, 
            alpha = pi/2,   
            mu = 1.0)
        );
    Ellipses_list = h_keys(Ellipses);
    h_set, Ellipses, Ellipses_list = Ellipses_list;
    
    return Ellipses;
}

func object_disk(Nnorm,nx,s_scl,d)
{
    DEG2RAD = pi/180;
    Ellipses = h_new(
        E8 = h_new( /* FIXME: "Support" ellipse : for absorption scaling */
            x0 = 0.*Nnorm,  
            y0 = 0.*Nnorm, 
            a = (d/(nx-1.))*Nnorm, 
            b = (d/(nx-1.))*Nnorm, 
            alpha = 0.0,   
            mu = 0.025)
        );
    Ellipses_list = h_keys(Ellipses);
    h_set, Ellipses, Ellipses_list = Ellipses_list;
    
    return Ellipses;
}

func object_tridisk(Nnorm)
{
    DEG2RAD = pi/180;
    Ellipses = h_new(
        E1 = h_new( /* FIXME: "Support" ellipse : for absorption scaling */
            x0 = 7.,  
            y0 = -7., 
            a = 12.0, 
            b = 4.0, 
            alpha = pi/4.0,   
            mu = 0.02),
        E2 = h_new( /* FIXME: "Support" ellipse : for absorption scaling */
            x0 = 6.,  
            y0 = 7., 
            a = 5.0, 
            b = 5.0, 
            alpha = 0.0,   
            mu = 0.025),
        E3 = h_new( /* FIXME: "Support" ellipse : for absorption scaling */
            x0 = -8.,  
            y0 = 0.0, 
            a = 7.0, 
            b = 7.0, 
            alpha = 0.0,   
            mu = 0.01)
        );
    Ellipses_list = h_keys(Ellipses);
    h_set, Ellipses, Ellipses_list = Ellipses_list;
    
    return Ellipses;
}
