#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cone_defs.h"
#include "cone_types.h"
#include "orbits.h"


/* Create orbit file and weight for a circle and line_number lines trajectory */

/* First orbit in the z=0 plane */
/*line_number lines at x=[rline*cos(phi),rline*sin(phi)] with 
                    phi=(k*2*pi/line_number), k=0 ... line_number-1 */

/* Detector always parallel to the z-axis */

main(argc,argv)

int argc;
char *argv[];

{
char *orbit_file, *M_mix_file;
ORBIT orbit, bbb;
int v,r;
float min,max;
int   k, k1, line_index;
double theta,costheta,sintheta;
double phi, dphi, hline;
float z, delz;
float radius_tuy;
float smoo();
LINOGRAM lino;

if(argc < 2) {
   printf("make_cl_orbit orbit_filename [M_mix_filename]\n");
   return -1;
}
orbit_file = argv[1];
if(argc > 2) M_mix_file = argv[2];
if(argc >2) 
  printf("\nPreparing orbit file and M_mix file for a circle and %d lines\n\n", line_number);
else
  printf("\nPreparing orbit file for a circle and %d lines\n\n", line_number);

dphi=2.0*PI/line_number;
if (line_number == 1) hline = hline1;
if (line_number == 2) hline = hline2;
if (line_number > 2 ) hline = hline3;

printf("Circle : %d points, radius : %.2f ; focal : %.2f\n", cb_points, cb_radius, cb_focal);
printf("Lines : %d points, - %.2f < z < + %.2f,  focal : %.2f\n", line_points,hline,hline,line_focal);
for(line_index = 0; line_index < line_number; ++line_index) 
  printf("Azimuthal position of the line %d is %f\n", line_index, (line_index*dphi));
printf("Circle sampling : %.2f mm, Line sampling : %.2f mm\n", (2*PI*cb_radius/cb_points), (2*hline/line_points));

/*radius_tuy = hline2*cb_radius/(sqrt(cb_radius*cb_radius+hline2*hline2));
printf("Tuy's condition satisfied for a sphere of radius %.2f mm\n",radius_tuy);*/
       

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
orbit[k].u_axis[2] = 0.0;
orbit[k].v_axis[0] = 0.0;
orbit[k].v_axis[1] = 0.0;
orbit[k].v_axis[2] = 1.0;
orbit[k].tangent[0] = -cb_radius*2.0*PI*sintheta/cb_points;
orbit[k].tangent[1] =  cb_radius*2.0*PI*costheta/cb_points;
orbit[k].tangent[2] = 0.0;
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

/*Lines*/
for(line_index=0;line_index<line_number;++line_index){
  z = -hline;
  delz = 2*hline/(line_points-1.0); /*sampling distance along line */
  for(k=0;k< line_points;++k) {
    k1 = k + cb_points + line_index*line_points;  /*position of record in orbit file*/
    orbit[k1].focus[0] = rline*cos(line_index*dphi);
    orbit[k1].focus[1] = rline*sin(line_index*dphi);
    orbit[k1].focus[2] = z;
    orbit[k1].focal_length = line_focal;
    orbit[k1].normal[0] = cos(line_index*dphi);
    orbit[k1].normal[1] = sin(line_index*dphi);
    orbit[k1].normal[2] = 0.0;
    orbit[k1].u_axis[0] = -1.0*sin(line_index*dphi);
    orbit[k1].u_axis[1] = cos(line_index*dphi);
    orbit[k1].u_axis[2] = 0.0;
    orbit[k1].v_axis[0] = 0.0;
    orbit[k1].v_axis[1] = 0.0;
    orbit[k1].v_axis[2] = 1.0;
    orbit[k1].tangent[0] = 0.0;
    orbit[k1].tangent[1] = 0.0;
    orbit[k1].tangent[2] = delz;
    orbit[k1].u_off     = u_offsino;
    /*next line : centre of detector shifted axially to keep FOV centred */
    orbit[k1].v_off     = v_offsino  + (z*line_focal/(rline*v_size));
    orbit[k1].center[0] = (rline-line_focal)*cos(line_index*dphi);
    orbit[k1].center[1] = (rline-line_focal)*sin(line_index*dphi);
    orbit[k1].center[2] = orbit[k1].focus[2];
    orbit[k1].mindex = k + cb_points;
    orbit[k1].smoothing = smoo(z,hline,delz);
    orbit[k1].discontinuity = 0;
    orbit[k1].neighbour = k1+1;
    if(k == 0) orbit[k1].discontinuity = 1;
    if(k == (line_points-1)) orbit[k1].neighbour = -1;

    if (argc > 2 && line_index == 0) {
      make_M_line(lino, orbit[k1], z, delz,hline, dphi);
      min = 0;
      max = 0;
      if( (write_image(M_mix_file,lino,(linopixels*dlinoviews),(orbit[k1].mindex-cb_points),&min,&max,(cb_points+line_points),linopixels)) != 0) {
	printf("Error writing %s\n",M_mix_file);
	return -1;
      }
    }
    
    z += delz;
  }
}

if( (write_orbit(orbit_file,orbit,(cb_points+line_number*line_points))) != 0) return -1;

/*reread to check*/

k = read_orbit(orbit_file,bbb);

return 0;
}


/*weight function along the line, ensures smoothness at its
extremities */

float smoo(z,hline,delz)

float z,hline,delz;

{
  float x,y;
  float p=4.0;

  x=sin(z*PI/(2.0*hline+delz));
  if(fabs(x) < 0.0001) y = 1.0;
  else y = 1 - exp(p*log(fabs(x)));
  return(y);
}


make_M_line(M,orb,z,delz,hline,dphi)

LINOGRAM M;
DETECTOR orb;
float z,delz,hline;
double dphi;

{

int i,k;
int iu,imu;
double w,wm,u[linopixels],centreoff;
double s,zline;
double mu[dlinoviews],co_mu[dlinoviews],si_mu[dlinoviews],wlino[dlinoviews];
double focus[3],tangent[3],theta[3],th_lbd[3];
double Mcircle,make_Mcircle();
float smoo();

/* some initializations */


for(k=0;k<linopixels;++k) u[k]=u_size*(k-linopixels/2.0);

for(k=0;k<dlinoviews;++k) {
  if (k < linoviews) {
    wlino[k]=2.0*(k-(linoviews-1.0)/2.0)/linoviews;
    mu[k]=atan(wlino[k]);
  }
  else {
    wlino[k]=2.0*(k-linoviews-(linoviews-1.0)/2.0)/linoviews;
    mu[k]=PI/2.0 + atan(wlino[k]);
  }
  wlino[k]=1.0/sqrt(1.0+wlino[k]*wlino[k]);
  co_mu[k]=cos(mu[k]);
  si_mu[k]=sin(mu[k]);
}


/* compute and store the M function associated with the line */

for(i=0;i<3;++i) focus[i]=orb.focus[i];  
for(i=0;i<3;++i) tangent[i]=orb.tangent[i];  

for(imu=0;imu<dlinoviews;++imu)
  for(iu=0;iu < linopixels;++iu){
    
    centreoff = (orb.u_off-u_offsino)*u_size*co_mu[imu]+
      (orb.v_off-v_offsino)*v_size*si_mu[imu];
    w=sqrt((u[iu]*wlino[imu]-centreoff)*(u[iu]*wlino[imu]-centreoff)
	   +orb.focal_length*orb.focal_length);
    th_lbd[0]=orb.focal_length*co_mu[imu]/w;
    th_lbd[1]=orb.focal_length*si_mu[imu]/w;
    th_lbd[2]=(u[iu]*wlino[imu]-centreoff)/w;
    for(i=0;i<3;++i) theta[i]=0.0;
    for(i=0;i<3;++i) theta[i]+=th_lbd[0]*orb.u_axis[i];
    for(i=0;i<3;++i) theta[i]+=th_lbd[1]*orb.v_axis[i];
    for(i=0;i<3;++i) theta[i]+=th_lbd[2]*orb.normal[i];
    for(i=0;i<3;++i) if (theta[i] > 1.0) theta[i]=1.0;
    for(i=0;i<3;++i) if (theta[i] < -1.0) theta[i]=-1.0;
    s=0.0;
    for(i=0;i<3;++i) s+=(theta[i]*focus[i]); 
    
    Mcircle=make_Mcircle(s,theta);
    
    if(line_number == 1) M[imu][iu]=1.0;
    else {
      M[imu][iu]=orb.smoothing;
      for(k=1;k<line_number;++k) {
	if(fabs(theta[2]) > 0.01)
	  zline = (s-theta[0]*rline*cos(k*dphi)-theta[1]*rline*sin(k*dphi))/theta[2];
	else zline = 2*hline;
        if (fabs(zline) <= hline) M[imu][iu] += smoo(zline,hline,delz);
      }
      M[imu][iu]=orb.smoothing/M[imu][iu];
    }
    wm=0.0;
    for(i=0;i<3;++i) wm+=(theta[i]*tangent[i]);
    wm=fabs(wm);
    wm=wm*w/(orb.focal_length*orb.focal_length);
    wm*=(-1.0/(4.0*PI*PI));
    M[imu][iu] *= (wm*(1.0-Mcircle));
  }
}

double make_Mcircle(s,theta)

double s,theta[3];

{ double Mcircle;
  double th,alpha,mu;

  if (theta[2]<0.0) {
    s*=-1;
    theta[2]*=-1.0;
  }
  
  th=acos(theta[2]);

  if (fabs(s)>= (cb_radius*sin(th))) Mcircle=0.0;
  else {
    mu=asin(theta[2]*cb_radius/sqrt(cb_radius*cb_radius-s*s));
    alpha=PI/2-mu;
/*    Mcircle=1.0-exp(-opening*opening*alpha*alpha);*/
  }

  return Mcircle;
}  




