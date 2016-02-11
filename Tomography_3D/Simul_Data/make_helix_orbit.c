#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cone_defs.h"
#include "cone_types.h"
#include "orbits.h"


/* Create orbit file and weight for a two turns helix */

/* Detector always parallel to the z-axis */

main(argc,argv)

int argc;
char *argv[];

{
char *orbit_file;
ORBIT orbit, bbb;
int v,r;
float min,max;
int   k, k1;
double theta,costheta,sintheta;
float z, delz;
float smoo();

if(argc < 2) {
   printf("make_helix_orbit orbit_filename\n");
   return -1;
}
orbit_file = argv[1];
printf("\nPreparing orbit file for a 2 turns helix\n\n");

printf("%d points, radius : %.2f ; focal : %.2f\n",
        helix_points, helix_radius, helix_focal);
printf("- %.2f < z < + %.2f\n",hhelix,hhelix);    

delz = 2*hhelix/(helix_points-1);
z = - hhelix;
 z = -80;


for(k=0;k< helix_points;++k) {

theta = 4*k*PI/(helix_points+1);
costheta = cos(theta);
sintheta = sin(theta);

orbit[k].focus[0] = helix_radius*costheta;
orbit[k].focus[1] = helix_radius*sintheta;
orbit[k].focus[2] = z;
orbit[k].focal_length = helix_focal;
orbit[k].normal[0] = costheta;
orbit[k].normal[1] = sintheta;
orbit[k].normal[2] = 0.0;
orbit[k].u_axis[0] = -sintheta;
orbit[k].u_axis[1] = costheta;
orbit[k].u_axis[2] = 0.0;
orbit[k].v_axis[0] = 0.0;
orbit[k].v_axis[1] = 0.0;
orbit[k].v_axis[2] = 1.0;
orbit[k].tangent[0] = -helix_radius*4.0*PI*sintheta/(helix_points+1);
orbit[k].tangent[1] =  helix_radius*4.0*PI*costheta/(helix_points+1);
orbit[k].tangent[2] =  delz;
orbit[k].u_off     = u_offsino;
//changement de cath : cb_focal => helix_focal et cb_radius => helix_radius
orbit[k].v_off     = v_offsino + (z*cb_focal/(cb_radius*v_size));
orbit[k].center[0] = (helix_radius-helix_focal)*costheta;
orbit[k].center[1] = (helix_radius-helix_focal)*sintheta;
orbit[k].center[2] = z;
orbit[k].mindex = k;
orbit[k].smoothing = smoo(z);
orbit[k].discontinuity = 0;
orbit[k].neighbour = k+1;
if (k == 0) orbit[k].discontinuity = 1;
if(k == (helix_points-1)) orbit[k].neighbour = -1;
z += delz;
}

if( (write_orbit(orbit_file,orbit,helix_points)) != 0) return -1;

/*reread to check*/

k = read_orbit(orbit_file,bbb);

return 0;
}

float smoo(z)

/*weight function along the line, ensures smoothness at its
extremities */

float z;

{
float x,y;
float p = 4.0;

x = sin(z*PI/(hhelix*2.0));
if(fabs(x) < 0.0001) y = 1.0;
else y = 1 - exp(p*log(fabs(x)));
return(y);
}
