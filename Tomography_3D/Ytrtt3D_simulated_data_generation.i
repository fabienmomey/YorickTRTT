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

func trtt_create_ellipses_phantom_3D(Ellipses, nx, ny, nz) 
/* DOCUMENT trtt_create_ellipses_phantom_3D, Ellipses, nx, ny, nz
   
   This function creates the 3D phantom image using the geometrical parameters
   of a 3D tomobj. A set of ellipses ELLIPSES is passed to the
   function. ELLIPSES is a hash table of the form:

   Ellipses           
   |
   |
   |----- E1
   |       |----- x0     // x-coordinate of the center in meters
   |       |----- y0     // y-coordinate of the center in meters
   |       |----- z0     // z-coordinate of the center in meters
   |       |----- a      // length of the x semi axis in meters
   |       |----- b      // length of the y semi axis in meters
   |       |----- c      // length of the z semi axis in meters 
   |       |----- phi    // phi Euler angle (in radians) (rotation about z-axis)
   |       |----- theta  // theta Euler angle (in radians) (rotation about x-axis)
   |       |----- psi    // psi Euler angle (in radians) (rotation about z-axis)
   |       |----- mu     // attenuation (constant on the surface of
                         // the ellipse)
   |
   |
   |----- E2
   |       |----- x0    
   |       |----- y0    
   |       |----- z0     
   |       |----- a     
   |       |----- b     
   |       |----- c    
   |       |----- phi   
   |       |----- theta  
   |       |----- psi    
   |       |----- mu     
   |
   |
   |----- E3
   :
   :  
   :
              
   This hash table requires a new entry per ellipse.

   This function is based on Matthias Schabel's MATLAB code for generating 3D
   Shepp-Logan phantoms.
   http://tomography.o-x-t.com/2008/04/13/3d-shepp-logan-phantom/
   http://www.mathworks.com/matlabcentral/fileexchange/9416

*/
{ 
    E_list = Ellipses.Ellipses_list;
    nE = numberof(E_list);

    /* Calculate coordinates */
    x = (indgen(nx)-((nx+1.0)/2.0))*(2.0/(nx-1.0));
    y = (indgen(ny)-((ny+1.0)/2.0))*(2.0/(ny-1.0));
    z = (indgen(nz)-((nz+1.0)/2.0))*(2.0/(nz-1.0));
    x = x(,-,-,-);
    y = y(-,,-,-);
    z = z(-,-,,-);

    shepp = array(double,nx,ny,nz,1);
    
    for (i=1; i<=nE; ++i)
    {
        // Get the current ellipse
        E = h_get(Ellipses, E_list(i));
                
        // Extract its geometrical parameters
        x0 = E.x0;
        y0 = E.y0;
        z0 = E.z0;
        asqr = (E.a)^2.0;
        bsqr = (E.b)^2.0;
        csqr = (E.c)^2.0;
        phi = E.phi;
        theta = E.theta;
        psi = E.psi;
        mu = E.mu;

        cphi = cos(phi);
        sphi = sin(phi);
        ctheta = cos(theta);
        stheta = sin(theta);
        cpsi = cos(psi);
        spsi = sin(psi);

        // Euler rotation matrix coefficients
        e1 = cpsi*cphi-ctheta*sphi*spsi;    e2 = cpsi*sphi+ctheta*cphi*spsi;    e3 = spsi*stheta;
        e4 = -spsi*cphi-ctheta*sphi*cpsi;   e5 = -spsi*sphi+ctheta*cphi*cpsi;   e6 = cpsi*stheta;
        e7 = stheta*sphi;                   e8 = -stheta*cphi;                  e9 = ctheta;

        xp = x-x0;
        yp = y-y0;
        zp = z-z0;
        
        w = e1*xp + e2*yp + e3*zp;
        v = e4*xp + e5*yp + e6*zp;
        u = e7*xp + e8*yp + e9*zp;

        // eq = (((w-x0)^2.0/asqr) + ((v-y0)^2.0/bsqr) + ((u-z0)^2.0/csqr) <= 1.0);
        eq = ((w^2.0/asqr) + (v^2.0/bsqr) + (u^2.0/csqr) <= 1.0);
        
        shepp += mu * eq;
    }

    return shepp;
}
trtt_create_phantom_3D = trtt_create_ellipses_phantom_3D;

func shepp_logan_3D(void)
{
    DEG2RAD = pi/180;
    Ellipses = h_new(
        E1 = h_new(
            x0 = 0., 
            y0 = 0.,
            z0 = 0.,
            a = 0.69, 
            b = 0.92,
            c = 0.81,
            phi = 0.0*DEG2RAD,
            theta = 0.0*DEG2RAD,
            psi = 0.0*DEG2RAD,
            mu = 2.0),
        E2 = h_new(
            x0 = 0., 
            y0 = -0.0184,
            z0 = 0.,
            a = 0.6624, 
            b = 0.874,
            c = 0.780,
            phi = 0.0*DEG2RAD,
            theta = 0.0*DEG2RAD,
            psi = 0.0*DEG2RAD,
            mu = -0.98),
        E3 = h_new(
            x0 = 0.22, 
            y0 = 0.,
            z0 = 0.,
            a = 0.11, 
            b = 0.31,
            c = 0.22,
            phi = -18.0*DEG2RAD,
            theta = 0.0*DEG2RAD,
            psi = 0.0*DEG2RAD,
            mu = -0.02),              
        E4 = h_new(
            x0 = -0.22, 
            y0 = 0.,
            z0 = 0.,
            a = 0.16, 
            b = 0.41,
            c = 0.28,
            phi = 18.0*DEG2RAD,
            theta = 0.0*DEG2RAD,
            psi = 0.0*DEG2RAD,
            mu = -0.02),
        E5 = h_new(
            x0 = 0., 
            y0 = 0.35,
            z0 = -0.15,
            a = 0.21, 
            b = 0.25,
            c = 0.41,
            phi = 0.0*DEG2RAD,
            theta = 0.0*DEG2RAD,
            psi = 0.0*DEG2RAD,
            mu = 0.01),
        E6 = h_new(
            x0 = 0., 
            y0 = 0.1,
            z0 = 0.25,
            a = 0.046, 
            b = 0.046,
            c = 0.05,
            phi = 0.0*DEG2RAD,
            theta = 0.0*DEG2RAD,
            psi = 0.0*DEG2RAD,
            mu = 0.01),
        E7 = h_new(
            x0 = 0., 
            y0 = -0.1,
            z0 = 0.25,
            a = 0.046, 
            b = 0.046,
            c = 0.05,
            phi = 0.0*DEG2RAD,
            theta = 0.0*DEG2RAD,
            psi = 0.0*DEG2RAD,
            mu = 0.01),
        E8 = h_new(
            x0 = -0.08, 
            y0 = -0.605,
            z0 = 0.,
            a = 0.046, 
            b = 0.023,
            c = 0.05,
            phi = 0.0*DEG2RAD,
            theta = 0.0*DEG2RAD,
            psi = 0.0*DEG2RAD,
            mu = 0.01),
        E9 = h_new(
            x0 = 0., 
            y0 = -0.606,
            z0 = 0.,
            a = 0.023, 
            b = 0.023,
            c = 0.02,
            phi = 0.0*DEG2RAD,
            theta = 0.0*DEG2RAD,
            psi = 0.0*DEG2RAD,
            mu = 0.01),
        E10 = h_new(
            x0 = 0.06, 
            y0 = -0.605,
            z0 = 0.,
            a = 0.023, 
            b = 0.046,
            c = 0.02,
            phi = 0.0*DEG2RAD,
            theta = 0.0*DEG2RAD,
            psi = 0.0*DEG2RAD,
            mu = 0.01)
        );
    Ellipses_list = h_keys(Ellipses);
    h_set, Ellipses, Ellipses_list = Ellipses_list;

    return Ellipses;
}
