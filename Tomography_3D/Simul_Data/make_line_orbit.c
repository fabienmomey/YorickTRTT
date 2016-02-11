#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cone_defs.h"
#include "cone_types.h"
#include "orbits.h"


/* Create orbit file and weight for a single line trajectory */

/* The line is at x=rline, y=0 */

/* Detector always parallel to the z-axis */

/* Relative weight of the line orbit is cb_radius**m multiplied by the function
   smoo(z) where z is the position of the source along the line orbit, and m is
   the power (usually = 4) in the general expression of the M function */


main(argc,argv)

int argc;
char *argv[];

{
char *orbit_file, *weight_file;
ORBIT orbit, bbb;
LINOGRAM lino;
int v,r;
float min,max;
int   k, k1;
double theta,costheta,sintheta;
float z, delz;
float radius_tuy;
float smoo();

if(argc < 2) {
   printf("make_line_orbit orbit_filename [weight_filename]\n");
   return -1;
}
orbit_file = argv[1];
if (argc > 2) weight_file = argv[2];
if (argc > 2)
printf("\nPreparing orbit file and M weight file for a line orbit\n");
else printf("\nPreparing orbit file for a line orbit\n\n");

printf("Line : %d points at y=0, x = %.2f, - %.2f < z < + %.2f,  focal : %.2f\n",
        line_points,rline,hline1,hline1,line_focal);

/*Line*/
z = -hline1;
delz = 2*hline1/(line_points-1.0); /*sampling distance along line */
for(k=0;k< line_points;++k) {

k1 = k; /*position of record in orbit file*/

orbit[k1].focus[0] = rline;
orbit[k1].focus[1] = 0;
orbit[k1].focus[2] = z;
orbit[k1].focal_length = line_focal;
orbit[k1].normal[0] = 1;
orbit[k1].normal[1] = 0.0;
orbit[k1].normal[2] = 0.0;
orbit[k1].u_axis[0] = 0.0;
orbit[k1].u_axis[1] = 1.0;
orbit[k1].u_axis[2] = 0.0;
orbit[k1].tangent[0] = 0.0;
orbit[k1].tangent[1] = 0.0;
orbit[k1].tangent[2] = delz;
orbit[k1].v_axis[0] = 0.0;
orbit[k1].v_axis[1] = 0.0;
orbit[k1].v_axis[2] = 1.0;
orbit[k1].u_off     = u_offsino;
/*next line : centre of detector shifted axially to keep FOV centred */
orbit[k1].v_off     = v_offsino  + (z*line_focal/(rline*v_size));
orbit[k1].center[0] = rline - line_focal;
orbit[k1].center[1] = 0.0;
orbit[k1].center[2] = orbit[k1].focus[2];
orbit[k1].mindex = k1;
orbit[k1].smoothing = smoo(z);
orbit[k1].discontinuity = 0;
orbit[k1].neighbour = k1+1;
if(k == 0) orbit[k1].discontinuity = 1;
if(k == (line_points-1)) orbit[k1].neighbour = -1;

if (argc > 2) {
  make_M_line(lino,z, orbit[k1].u_off, orbit[k1].v_off);
  min = 0;
  max = 0;
  if( (write_image(weight_file,lino,(linopixels*dlinoviews),k1,&min,&max,line_points,linopixels)) != 0) {
    printf("Error writing %s\n",weight_file);
    return -1;
  }
}

z += delz;
}

if( (write_orbit(orbit_file,orbit,line_points)) != 0) return -1;

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

x = sin(z*PI/(2.01*hline1));
if(fabs(x) < 0.0001) y = 1.0;
else y = 1 - exp(p*log(fabs(x)));
return(y);
}

make_M_line(lino,eta,u_off,v_off)

LINOGRAM lino;
float eta;
float u_off,v_off;

/*Produce weight for the circle and line orbit, this routine for line*/
/*Takes into account the fact that the weight is applied to a sinogram
  with radial samples shifted by cb_sino_offset (defined in cone_defs.h)*/

{
float smoo();
int v,r;
float mu,cosmu, sinmu;
float cosmum,sinmum,yy,xx,s;
float p,rescale;
float u[linopixels],wlino;
float centreoff;
int m = 4; /*free parameter*/

rescale = (2*hline1/(line_points-1.0))/line_focal;
p = -rescale/(4*PI*PI);
for(r=0;r<linopixels;++r) u[r]=u_size*(r-linopixels/2.0); 

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
  if(fabs(sinmu) > 0.00001) sinmum = exp(m*log(fabs(sinmu)));
  else sinmum = 0.0;

  centreoff = (u_off - u_offsino)*u_size*cosmu + (v_off - v_offsino)*v_size*sinmu;
  
  for(r=0;r<linopixels;++r) {
    s=u[r]/wlino - centreoff;
    xx = (s*s)/(line_focal*line_focal) + cosmu*cosmu;
    yy = ((s*rline/line_focal) + (eta*sinmu))/cb_radius; 
    xx = xx - yy*yy;

    /*For the current use of line data in the mix algorithm*/
    if(xx > 0.0001) { /*this plane also intersects the circular orbit*/
      lino[v][r] = 0.0;
    }
    else lino[v][r] = p*fabs(sinmu);

    /*For the optimal use of line data*/
    lino[v][r]=p*fabs(sinmu);
  }
}
}

