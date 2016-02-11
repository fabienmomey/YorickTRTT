#include <stdio.h>
#include <math.h>
#include "cone_defs.h"
#include "cone_types.h"
#include "orbits.h"

/* Create orbit file for a single circular trajectory, and optionally
   a constant M function file (i.e. |cos(mu)|) for the cone-beam FBP algorithm.

   Detector parallel to the z-axis.
*/

main(argc,argv)

int argc;
char *argv[];

{
char *orbit_file, *weight_file;
LINOGRAM lino;
ORBIT orbit, bbb; 
int   k,r,v;
double theta,costheta,sintheta;
float p,w,mu,min,max;
float z; //ajout cath

if(argc < 2) {
   printf("make_CB_orbit orbit_filename [weight_file]\n");
   return -1;
}
orbit_file = argv[1];
if(argc > 2) {
   weight_file = argv[2];
   printf("\nPreparing orbit file and M weight file for a single circular orbit\n\n");
}
else  printf("\nPreparing orbit file for a single circular orbit\n\n");

printf("%d point orbit with radius : %.2f  ; detector focal length : %.2f\n",
cb_points,1.0*cb_radius,1.0*cb_focal);

/*prepare the weight sinogram (identical for all orbit points)*/
p = -cb_radius/(4.0*PI*cb_focal*cb_points);

if (argc > 2) {
   for(v=0;v<dlinoviews;++v) {
     if (v < linoviews) mu = atan(2.0*(v-(linoviews-1.0)/2.0)/linoviews);
     else mu = PI/2.0 + atan(2.0*(v-linoviews-(linoviews-1.0)/2.0)/linoviews);
     w =  p*fabs(cos(mu));
      for(r=0;r<linopixels;++r) lino[v][r] = w;
   }
}

 z=0;

for(k=0;k< cb_points;++k) {
//FIXME: Compatible TRTT
theta = 2*k*PI/cb_points;
/* theta = 2*k*PI/cb_points + PI/cb_points ; */
costheta = cos(theta);
sintheta = sin(theta);

orbit[k].focus[0] = cb_radius*costheta;
orbit[k].focus[1] = cb_radius*sintheta;
orbit[k].focus[2] = z;
orbit[k].focal_length = cb_focal;
orbit[k].normal[0] = costheta;
orbit[k].normal[1] = sintheta;
orbit[k].normal[2] = 0.0;
orbit[k].u_axis[0] = -sintheta;
orbit[k].u_axis[1] = costheta;
orbit[k].u_axis[2] = 0.0;
orbit[k].tangent[0] = -cb_radius*2.0*PI*sintheta/cb_points;
orbit[k].tangent[1] =  cb_radius*2.0*PI*costheta/cb_points;
orbit[k].tangent[2] = 0.0;
orbit[k].v_axis[0] = 0.0;
orbit[k].v_axis[1] = 0.0;
orbit[k].v_axis[2] = 1.0;
orbit[k].u_off     = u_offsino;
orbit[k].v_off     = v_offsino;
orbit[k].center[0] = (cb_radius - cb_focal)*costheta;
orbit[k].center[1] = (cb_radius - cb_focal)*sintheta;
orbit[k].center[2] = z;
orbit[k].mindex = 0;  /*All orbit points use the same M weight */
orbit[k].neighbour = k+1;
orbit[k].discontinuity = 0;
}
orbit[cb_points-1].neighbour = 0;
orbit[0].discontinuity = 1;

/* Write the M weight file. Only one record since it is symmetrical*/

if(argc > 2) {
   min = 0;
   max = 0;
   if( (write_image(weight_file,lino,(linopixels*dlinoviews),0,
        &min,&max,1,linopixels)) != 0) {
      printf("Error writing %s\n",weight_file);
      return -1;
   }
}

/*Write the orbit file */

if( (write_orbit(argv[1],orbit,cb_points)) != 0) return -1;

/*reread to check*/

k = read_orbit(argv[1],bbb);

return 0;
}
