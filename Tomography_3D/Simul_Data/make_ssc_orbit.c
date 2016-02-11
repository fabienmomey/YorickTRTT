#include <stdio.h>
#include <math.h>
#include "cone_defs.h"
#include "cone_types.h"
#include "orbits.h"

/* Create orbit file for a single circular trajectory, and optionally
   a constant M function file (i.e. |cos(mu)|) for the cone-beam FBP 
   algorithm.

   Detector parallel to the z-axis.*/

main(argc,argv)

int argc;
char *argv[];

{
  char *orbit_file;
  ORBIT orbit, bbb;
  int   k,r,v;
  double theta,costheta,sintheta;
  float min,max;
  float gm,dtheta;
  float smoo();

  if(argc < 2) {
   printf("make_ssc_orbit orbit_filename\n");
   return -1;
  }
  orbit_file = argv[1];
  printf("\nPreparing orbit file for a circular short-scan\n\n");

  printf("%d points ; radius : %.2f ; focal length : %.2f\n",cb_points,1.0*cb_radius,1.0*cb_focal);

  gm=30.0*PI/180.0;
  dtheta=(PI+2.0*gm)/cb_points;
  for(k=0;k< cb_points;++k) {
    
    theta = -PI/2.0-gm+k*dtheta;
    costheta = cos(theta);
    sintheta = sin(theta);
    
    orbit[k].focus[0] = cb_radius*costheta;
    orbit[k].focus[1] = cb_radius*sintheta;
    orbit[k].focus[2] = 0.0;
    orbit[k].focal_length = cb_focal;
    orbit[k].normal[0] = costheta;
    orbit[k].normal[1] = sintheta;
    orbit[k].normal[2] = 0.0;
    orbit[k].u_axis[0] = -sintheta;
    orbit[k].u_axis[1] = costheta;
    orbit[k].u_axis[2] = 0.0;
    orbit[k].tangent[0] = -cb_radius*dtheta*sintheta;
    orbit[k].tangent[1] =  cb_radius*dtheta*costheta;
    orbit[k].tangent[2] = 0.0;
    orbit[k].v_axis[0] = 0.0;
    orbit[k].v_axis[1] = 0.0;
    orbit[k].v_axis[2] = 1.0;
    orbit[k].u_off     = u_offsino;
    orbit[k].v_off     = v_offsino;
    orbit[k].center[0] = (cb_radius - cb_focal)*costheta;
    orbit[k].center[1] = (cb_radius - cb_focal)*sintheta;
    orbit[k].center[2] = 0.0;
    orbit[k].mindex = k;  /*All orbit points use the same M weight */
    orbit[k].smoothing = smoo(theta,dtheta,gm);
    orbit[k].neighbour = k+1;
    orbit[k].discontinuity = 0;
  }
  orbit[cb_points-1].neighbour = 0;
  orbit[0].discontinuity = 1;
  
  /*Write the orbit file */
  
  if( (write_orbit(argv[1],orbit,cb_points)) != 0) return -1;
  
  /*reread to check*/
  
  k = read_orbit(argv[1],bbb);

  return 0;
}

float smoo(theta,dtheta,gm)

float theta,dtheta,gm;

{
  float x,y;
  float p = 4.0;
 
  x=sin(0.5*theta*PI/(PI/2.0+dtheta/2.0+gm));
  y=1.0-x*x;
  return(y);
}
