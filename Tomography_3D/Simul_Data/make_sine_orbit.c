#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cone_defs.h"
#include "cone_types.h"
#include "orbits.h"
#include "scan_types.h"


/* Create orbit file for one of the cylindrical orbits */
/* List of orbits being implemented:
   x  helix:  2 turns over specified height
   [X] sine:   3 periods of sinewave
   cross:  2 intersecting "circles"
   dhelix: single turn of double helix, with caps
   cir_3l: circle and three line
*/

/* Detector is parallel to the z-axis which is center line of cylinder,
   detector moves in opposite direction to vertex point. 
   Most of the detector parameters are found in the file cone_defs.h */

/*  This program is for the SINE orbit  */


main(argc,argv)

     int argc;
     char *argv[];

{
  char *orbit_file;
  ORBIT orbit, bbb;
  TPARSCAN Tscan;
  int   k;
  double z,phi,cosphi,sinphi;
  int intervals;


  if(argc < 2) {
    printf("make_sine_orbit orbit_filename\n");
    return -1;
  }
  orbit_file = argv[1];

  /* The following are constants found in the file cyl_orbits.h */

  printf("Cylindrical orbit with %d vertex points\n",sine_points);
  printf("                half height:    %.2f   \n",1.0*hsine);

  printf("\nPreparing sine orbit file \n\n");

  intervals = sine_points;

  for(k=0; k< sine_points; ++k) {    

    //phi = 2*PI*k/intervals + 3.0/intervals;
    phi = 2*PI*k/intervals + PI/intervals; //pour compatibilité avec CB
    //phi = 19/18*PI*k/intervals + 19/36*PI/intervals; //pour simulation mire Ottawa
    cosphi = cos(phi);
    sinphi = sin(phi);
    z = sin(3*phi)*hsine; 

    orbit[k].focus[0] = sine_radius*cosphi;
    orbit[k].focus[1] = sine_radius*sinphi;
    orbit[k].focus[2] = z;
    orbit[k].focal_length = sine_focal;
    orbit[k].normal[0] = cosphi;
    orbit[k].normal[1] = sinphi;
    orbit[k].normal[2] = 0.0;
    orbit[k].u_axis[0] = -sinphi;
    orbit[k].u_axis[1] = cosphi;
    orbit[k].u_axis[2] = 0.0;
    orbit[k].v_axis[0] = 0.0;
    orbit[k].v_axis[1] = 0.0;
    orbit[k].v_axis[2] = 1.0;
    orbit[k].tangent[0] = -sine_radius*sinphi*2*PI/intervals;
    orbit[k].tangent[1] =  sine_radius*cosphi*2*PI/intervals;
    orbit[k].tangent[2] =  3.0*hsine*cos(3.0*phi)*2*PI/intervals;
    orbit[k].u_off     = u_offsino;
    orbit[k].v_off     = v_offsino + (z*sine_focal/(sine_radius*v_size));
    orbit[k].center[0] = (sine_radius-sine_focal) * cosphi;
    orbit[k].center[1] = (sine_radius-sine_focal) * sinphi;
    orbit[k].center[2] = z; 
    orbit[k].mindex = k;
    orbit[k].smoothing = 1.0;
    orbit[k].discontinuity = 0;
    orbit[k].neighbour = k+1;
    if (k == 0) orbit[k].discontinuity = 1;
    if (k == (sine_points-1)) orbit[k].neighbour = 0;    
  }

  // changement de paramétrisation: orbit vers T scan
  //ParIntExt_ParScan(orbit,sine_points,Tscan);//perturbation chaotique  
  // changement de paramétrisation: Tscan vers orbiT 
  //ParScan_ParIntExt(orbit,sine_points,Tscan);

  if( (write_orbit(orbit_file,orbit,sine_points)) != 0) return -1;

  /*reread to check*/

  k = read_orbit(orbit_file,bbb);

  return 0;
}
