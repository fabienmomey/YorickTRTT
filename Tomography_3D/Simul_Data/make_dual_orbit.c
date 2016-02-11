#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cone_defs.h"
#include "cone_types.h"
#include "orbits.h"


/* Create orbit file and weight for a dual circular trajectory */

/* First orbit in the z=0 plane, second orbit in the x=0 plane */

/* Detector parallel to the z-axis(first orbit) or
   parallel to the x_axis (second orbit) */


main(argc,argv)

int argc;
char *argv[];

{
char *orbit_file, *weight_file;
ORBIT orbit, bbb;
LINOGRAM lino;
int v,r;
int dual_points2;
float min,max;
int   k, k1;
double theta,costheta,sintheta;
float rescale;


if(argc < 2) {
   printf("make_dual_orbit orbit_filename [weight_filename]\n");
   return -1;
}
orbit_file = argv[1];
if (argc > 2) weight_file = argv[2];

if (argc > 2)
printf("\nPreparing orbit file and M weight file for a dual circular orbit\n");
else printf("\nPreparing orbit file for a dual circular orbit\n\n");

if((dual_points % 2) != 0) {
   printf("Number of orbit points should be even\n");
   return -1;
}
dual_points2 = dual_points/2;
rescale = (2*PI*dual_radius)/(dual_points2*dual_focal);

printf("Two orthogonal %d point orbits, radius : %.2f  ; detector focal length : %.2f\n",
dual_points2,1.0*dual_radius,1.0*dual_focal);



/*First orbit*/
for(k=0;k< dual_points2;++k) {

theta = 2*k*PI/dual_points2 ;
costheta = cos(theta);
sintheta = sin(theta);

orbit[k].focus[0] = dual_radius*costheta;
orbit[k].focus[1] = dual_radius*sintheta;
orbit[k].focus[2] = 0.0;
orbit[k].focal_length = dual_focal;
orbit[k].normal[0] = costheta;
orbit[k].normal[1] = sintheta;
orbit[k].normal[2] = 0.0;
orbit[k].u_axis[0] = -sintheta;
orbit[k].u_axis[1] = costheta;
orbit[k].u_axis[2] = 0.0;
orbit[k].tangent[0] = -dual_radius*2.0*PI*sintheta/dual_points2;
orbit[k].tangent[1] =  dual_radius*2.0*PI*costheta/dual_points2;
orbit[k].tangent[2] = 0.0;
orbit[k].v_axis[0] = 0.0;
orbit[k].v_axis[1] = 0.0;
orbit[k].v_axis[2] = 1.0;
orbit[k].u_off     = u_offsino;
orbit[k].v_off     = v_offsino;
orbit[k].center[0] = (dual_radius-dual_focal)*costheta;
orbit[k].center[1] = (dual_radius-dual_focal)*sintheta;
orbit[k].center[2] = 0.0;
orbit[k].mindex = k;
orbit[k].smoothing = 1.0;
orbit[k].neighbour = k+1;
orbit[k].discontinuity = 0;
if (k == 0) orbit[k].discontinuity = 1;
if (k == (dual_points2-1)) orbit[k].neighbour = 0;
if (argc > 2) {
  make_M_weight(lino,costheta,sintheta,rescale,orbit[k].u_off, orbit[k].v_off);
  min = 0;
  max = 0;
  if( (write_image(weight_file,lino,(linopixels*dlinoviews),k,&min,&max,dual_points,linopixels)) != 0) {
    printf("Error writing %s\n",weight_file);
    return -1;
  }
}

}

/*Second orbit*/
for(k=0;k< dual_points2;++k) {

theta = 2*k*PI/dual_points2 ;
costheta = cos(theta);
sintheta = sin(theta);
k1 = k+ dual_points2;
orbit[k1].focus[0] = 0;
orbit[k1].focus[1] = dual_radius*sintheta;
orbit[k1].focus[2] = dual_radius*costheta;
orbit[k1].focal_length = dual_focal;
orbit[k1].normal[0] = 0.0;
orbit[k1].normal[1] = sintheta;
orbit[k1].normal[2] = costheta;
orbit[k1].u_axis[0] = 0.0;
orbit[k1].u_axis[1] = costheta;
orbit[k1].u_axis[2] = -sintheta;
orbit[k1].tangent[0] = 0.0;
orbit[k1].tangent[1] = dual_radius*2.0*PI*costheta/dual_points2;
orbit[k1].tangent[2] = -dual_radius*2.0*PI*sintheta/dual_points2;
orbit[k1].v_axis[0] = 1.0;
orbit[k1].v_axis[1] = 0.0;
orbit[k1].v_axis[2] = 0.0;
orbit[k1].u_off     = u_offsino;
orbit[k1].v_off     = v_offsino;
orbit[k1].center[0] = 0.0;
orbit[k1].center[1] = (dual_radius-dual_focal)*sintheta;
orbit[k1].center[2] = (dual_radius-dual_focal)*costheta;
orbit[k1].mindex = k1;
orbit[k1].smoothing = 1.0;
orbit[k1].discontinuity = 0;
orbit[k1].neighbour = k1+1;
if (k == 0) orbit[k1].discontinuity = 1;
if (k == (dual_points2-1)) orbit[k1].neighbour = dual_points2;


if (argc > 2) {
  make_M_weight(lino,costheta,sintheta,rescale,
		orbit[k1].u_off, orbit[k1].v_off);
  min = 0;
  max = 0;
  if( (write_image(weight_file,lino,(linopixels*dlinoviews),(k1),&min,&max,dual_points,linopixels)) != 0) {
    printf("Error writing %s\n",weight_file);
    return -1;
  }
}

}

if( (write_orbit(orbit_file,orbit,dual_points)) != 0) return -1;

/*reread to check*/

k = read_orbit(orbit_file,bbb);

return 0;
}


make_M_weight(lino,costheta,sintheta,rescale,u_off,v_off)

LINOGRAM lino;
double costheta,sintheta;
float rescale;
float u_off,v_off;

/*Produce weight for the dual orbit trajectory*/

{
int v,r;
float mu,cosmu, sinmu;
float cosmum,xx,s;
float del,p;
float u[linopixels],wlino;
float centreoff;
int m = 4; /*free parameter*/

p = -rescale/(8*PI*PI);
del = costheta/dual_focal;
for(r=0;r<linopixels;++r) u[r]=u_size*del*(r-linopixels/2.0); 

for(v=0;v<dlinoviews;++v) {

  if (v < linoviews) {
    wlino=2.0*(v-(linoviews-1.0)/2.0)/linoviews;
    mu=atan(wlino);
  }
  else {
    wlino=2.0*(v-linoviews-(linoviews-1.0)/2.0)/linoviews;
    mu=PI/2.0 + atan(wlino);
  }
  wlino=sqrt(1.0+wlino*wlino);
  cosmu = cos(mu);
  sinmu = sin(mu);
  cosmum = exp(m*log(fabs(cosmu)));

  centreoff = (u_off - u_offsino)*u_size*cosmu + (v_off - v_offsino)*v_size*sinmu;
  centreoff *= del;

  for(r=0;r<linopixels;++r) {
    s = u[r]/wlino - centreoff - cosmu*sintheta;
    xx = 1.0 - s*s;
    if(xx > 0) xx = exp((m/2)*log(xx));
    else xx = 0.0;
    lino[v][r] = p*fabs(cosmu)*cosmum/(cosmum + xx);
   }
  
}
}
