#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cone_defs.h"
#include "cone_types.h"
#include "orbits.h"


/* Create orbit file and weight for a circle and helix trajectory */

/* First orbit in the z=0 plane */

/* Detector always parallel to the z-axis */


main(argc,argv)

int argc;
char *argv[];

{
char *orbit_file;
ORBIT orbit, bbb;
LINOGRAM lino;
int v,r;
float min,max;
int   k, k1;
double theta,costheta,sintheta;
float z, delz, shift;
float radius_tuy;
float smoo();


if(argc < 2) {
   printf("make_cirhel_orbit orbit_filename \n");
   return -1;
}
orbit_file = argv[1];
printf("\nPreparing orbit file for a circle and helix orbit\n");

printf("Circle and helix orbit. Circle : %d points, radius : %.2f, focal : %.2f\n",cb_points, cb_radius, cb_focal);
printf("Helix : %d points, radius : %.2f, focal : %.2f\n",
        helix_points,helix_radius,helix_focal);
printf("        - %.2f < z < + %.2f\n",hhelix,hhelix);    
printf("Circle sampling : %.2f mm, helix sampling (on z): %.2f mm\n", (2*PI*cb_radius/cb_points),(2*hhelix/(helix_points-1.0)));

radius_tuy=100.0;
printf("Tuy's condition satisfied for a sphere of radius %.2f mm\n",radius_tuy);
       

/*Circle*/
for(k=0;k< cb_points;++k) {

theta = 2*k*PI/cb_points;
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
orbit[k].tangent[0] = -cb_radius*2.0*PI*sintheta/cb_points;
orbit[k].tangent[1] =  cb_radius*2.0*PI*costheta/cb_points;
orbit[k].tangent[2] = 0.0;
orbit[k].u_axis[2] = 0.0;
orbit[k].v_axis[0] = 0.0;
orbit[k].v_axis[1] = 0.0;
orbit[k].v_axis[2] = 1.0;
orbit[k].u_off     = u_offsino;
orbit[k].v_off     = v_offsino;
orbit[k].center[0] = (cb_radius-cb_focal)*costheta;
orbit[k].center[1] = (cb_radius-cb_focal)*sintheta;
orbit[k].center[2] = 0.0;
orbit[k].mindex = k;
orbit[k].smoothing = 1.0;
orbit[k].discontinuity = 0;
orbit[k].neighbour = k+1;
if (k == 0) orbit[k].discontinuity = 1;
if(k == (cb_points-1)) orbit[k].neighbour = 0;

}


/*Helix*/

z = -hhelix;
delz = 2*hhelix/(helix_points-1.0); /*sampling distance along helix axis */

for(k=0;k< helix_points;++k) {

k1 = k + cb_points; /*position of record in orbit file*/

theta = 2*k*PI/(helix_points-1);
costheta = cos(theta);
sintheta = sin(theta);
orbit[k1].focus[0] = helix_radius*costheta;
orbit[k1].focus[1] = helix_radius*sintheta;
orbit[k1].focus[2] = z;
orbit[k1].focal_length = helix_focal;
orbit[k1].normal[0] = costheta;
orbit[k1].normal[1] = sintheta;
orbit[k1].normal[2] = 0.0;
orbit[k1].u_axis[0] = -sintheta;
orbit[k1].u_axis[1] = costheta;
orbit[k1].u_axis[2] = 0.0;
orbit[k1].v_axis[0] = 0.0;
orbit[k1].v_axis[1] = 0.0;
orbit[k1].v_axis[2] = 1.0;
orbit[k1].tangent[0] = -helix_radius*2.0*PI*sintheta/(helix_points-1);
orbit[k1].tangent[1] =  helix_radius*2.0*PI*costheta/(helix_points-1);
orbit[k1].tangent[2] =  delz;
orbit[k1].u_off     = u_offsino;
/*next point : centre of detector shifted axially to keep FOV centred */
/*orbit[k1].v_off     = v_offsino  + (z*cb_focal/(cb_radius*v_size));*/
orbit[k1].v_off     = v_offsino  + (z*cb_focal/(cb_radius*v_size));
orbit[k1].center[0] = (helix_radius-helix_focal)*costheta;
orbit[k1].center[1] = (helix_radius-helix_focal)*sintheta;
orbit[k1].center[2] = z;
orbit[k1].mindex = k1;
orbit[k1].smoothing = smoo(z,delz);
orbit[k1].discontinuity = 0;
orbit[k1].neighbour = k1+1;
if(k == 0) orbit[k1].discontinuity = 1;
if(k == (helix_points-1)) orbit[k1].neighbour = -1;

z += delz;
}

if( (write_orbit(orbit_file,orbit,(cb_points+helix_points))) != 0) return -1;

/*reread to check*/

k = read_orbit(orbit_file,bbb);

return 0;
}


float smoo(z,dz)

/*weight function along the helix axis, ensures smoothness at its
extremities */

float z,dz;

{
float x,y;
float p = 4.0;

x = sin(z*PI/((hhelix+dz/2.0)*2.0));
if(fabs(x) < 0.0001) y = 1.0;
else y = 1 - exp(p*log(fabs(x)));
return(y);
}



