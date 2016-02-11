#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cone_defs.h"
#include "cone_types.h"
#include "orbits.h"


/* Create orbit file and weight for a curl trajectory 
      
                          |
                       |  |  |
                      _|__|__|_  
                     / |  |  | \
                     \_|  |  |_/
                       |  |  |
                       |  |  |
                       |     | 
*/

/* Circular basis of the curl in the z=0 plane from phi=pi/2-arcsin(r/R) 
   to phi=3*pi/2+arcsin(r/R). 
   First  vertical line at phi=pi/2-arcsin(r/R),   -hcurl < z < hcurl.
   Second vertical line at phi=3*pi/2+arcsin(r/R), -hcurl < z < hcurl.
   Third  vertical line at phi=pi,                 -hcurl < z < hcurl.
   R is curl_radius. r is the radius of the cylindrical or spherical support 
   of the image */ 


/* Detector always parallel to the z-axis */


main(argc,argv)

int argc;
char *argv[];

{
  char *orbit_file;
  ORBIT orbit, bbb;
  float min,max;
  int  i, k, k1;
  double phim,theta[3],dtheta,costheta,sintheta;
  float z, delz;
  float smoo_line(),smoo_circle();
  
  if(argc < 2) {
    printf("make_curl_orbit orbit_filename\n");
    return -1;
  }
  orbit_file = argv[1];
  
  phim=asin(object_radius/curl_radius);
  
  printf("\nPreparing orbit file for a curl orbit\n\n");
  
  printf("Curl orbit. Circular basis : %d points, radius : %.2f ; focal : %.2f\n", curlcircle_points, curl_radius, curl_focal);
  printf("First line : %d points at phi = %.2f, - %.2f < z < + %.2f,  focal : %.2f\n", curlline_points,(PI/2-phim),hcurl,hcurl,curl_focal);
  printf("Second line : %d points at phi = %.2f, - %.2f < z < %.2f,  focal : %.2f\n", curlline_points,(3*PI/2+phim),hcurl,hcurl,curl_focal);
  printf("Third line : %d points at phi = %.2f, - %.2f < z < %.2f,  focal : %.2f\n", curlline_points,PI,hcurl,hcurl,curl_focal);
  printf("Circular basis sampling : %.2f mm, Line sampling : %.2f mm\n", (curl_radius*(PI+2*phim)/(curlcircle_points-1)), (2*hcurl/(curlline_points-1)));
  
  
  /*Create the three lines*/

  theta[0]=PI/2-phim;
  theta[1]=3*PI/2+phim;
  theta[2]=PI;
  for(i=0;i<3;++i) {
    theta[i]=PI/2-phim;
    costheta=cos(theta[i]);
    sintheta=sin(theta[i]);
    z = -hcurl;
    delz = 2*hcurl/(curlline_points-1); /*sampling distance along line */

    for(k=0;k<curlline_points;++k) {

      k1=k+i*curlline_points;
      orbit[k1].focus[0] = curl_radius*costheta;
      orbit[k1].focus[1] = curl_radius*sintheta;
      orbit[k1].focus[2] = z;
      orbit[k1].focal_length = curl_focal;
      orbit[k1].normal[0] = costheta;
      orbit[k1].normal[1] = sintheta;
      orbit[k1].normal[2] = 0.0;
      orbit[k1].u_axis[0] = -sintheta;
      orbit[k1].u_axis[1] = costheta;
      orbit[k1].u_axis[2] = 0.0;
      orbit[k1].v_axis[0] = 0.0;
      orbit[k1].v_axis[1] = 0.0;
      orbit[k1].v_axis[2] = 1.0;
      orbit[k1].tangent[0] = 0.0;
      orbit[k1].tangent[1] = 0.0;
      orbit[k1].tangent[2] = delz;
      orbit[k1].u_off = u_offsino;
      orbit[k1].v_off = v_offsino  + (z*curl_focal/(curl_radius*v_size));
      orbit[k1].center[0] = (curl_radius-curl_focal)*costheta;
      orbit[k1].center[1] = (curl_radius-curl_focal)*sintheta;
      orbit[k1].center[2] = z;
      orbit[k1].mindex = k1;
      orbit[k1].smoothing = smoo_line(z);
      orbit[k1].discontinuity = 0;
      orbit[k1].neighbour = k1+1;
      if(k == 0) orbit[k1].discontinuity = 1;
      if(k == curlline_points-1) orbit[k1].neighbour = -1;

      z += delz;
    }
  }


  /*Circular basis*/

  dtheta=(PI+2*phim)/(curlcircle_points-1);
  theta[0] = PI/2-phim;
  for(k=0;k< curlcircle_points;++k) {
    
    k1=k+3*curlline_points;
    costheta = cos(theta[0]);
    sintheta = sin(theta[0]);
    
    orbit[k1].focus[0] = curl_radius*costheta;
    orbit[k1].focus[1] = curl_radius*sintheta;
    orbit[k1].focus[2] = 0.0;
    orbit[k1].focal_length = curl_focal;
    orbit[k1].normal[0] = costheta;
    orbit[k1].normal[1] = sintheta;
    orbit[k1].normal[2] = 0.0;
    orbit[k1].u_axis[0] = -sintheta;
    orbit[k1].u_axis[1] = costheta;
    orbit[k1].u_axis[2] = 0.0;
    orbit[k1].v_axis[0] = 0.0;
    orbit[k1].v_axis[1] = 0.0;
    orbit[k1].v_axis[2] = 1.0;
    orbit[k1].tangent[0] = -curl_radius*dtheta*sintheta;
    orbit[k1].tangent[1] =  curl_radius*dtheta*costheta;
    orbit[k1].tangent[2] = 0.0;
    orbit[k1].u_off     = u_offsino;
    orbit[k1].v_off     = v_offsino;
    orbit[k1].center[0] = (curl_radius-curl_focal)*costheta;
    orbit[k1].center[1] = (curl_radius-curl_focal)*sintheta;
    orbit[k1].center[2] = 0.0;
    orbit[k1].mindex = k1;
    orbit[k1].smoothing = smoo_circle(theta,phim);
    orbit[k1].discontinuity = 0;
    orbit[k1].neighbour = k1+1;
    if(k == 0) orbit[k1].discontinuity = 1;
    if(k == curlcircle_points-1) orbit[k1].neighbour = -1;

    theta[0]+=dtheta;
  }

write_orbit(orbit_file,orbit,(curlcircle_points+3*curlline_points));
}

/********************************************************/
float smoo_line(z)

/*weight function along the line, ensures smoothness at its
extremities */

float z;

{
double x,y;
double p = 2.0;

x = sin(z*PI/(2.1*hcurl));
if(fabs(x) < 0.0001) y = 1.0;
else y = 1 - exp(p*log(fabs(x)));
return (float) (y);
}

/********************************************************/
float smoo_circle(theta,phim)

/*weight function along the circle, ensures smoothness at its
extremities */

double theta,phim;

{
double x,y;
double p = 2.0;

x = sin( PI*(1.0+cos(theta))/(2.3*(1.0+sin(phim))) );
if(fabs(x) < 0.0001) y = 1.0;
else y = 1 - exp(p*log(fabs(x)));
return (float)(y);
}

