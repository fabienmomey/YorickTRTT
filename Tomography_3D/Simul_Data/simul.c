#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cone_defs.h"
#include "cone_types.h"
#include "simul.h"

/*new version December 1996 includes gaussian objects*/

/*---------------------------------------------------------------------------
The text file description of a phantom should have the following format :

phantom

object 1

ellipsoid
center 2.34 4.56 6.78
half-axis 4.56 6.777 0.123
angles 80.0 0.0

0.123


object 2

cylinder
radius 0.123
length 10
center -2 3 -4.56
angles 80.0 0.0
0.123

object 3

gaussian
center 2.34 4.56 6.78
sigma  4.56 6.777 0.123
angles 80.0 0.0

0.123


end

----------------------------------------------------------------------------*/

float line_integral(a,b,ob)

/*integral of object ob from point a to point b*/

VECTOR a,b;
OBJECT ob;

{
VECTOR del,ap,ar,br,bp;
double dels,determ,c1,c2,c0;
float l2; 
double path1,path2,z1,z2;
int i,j,k;
   if (ob.geometry== ellipsoid)
   {
      /*shift coordinates*/
      for(i=0;i<3;++i) {
         ap[i] = (a[i] - (ob.e.center[i]));
         bp[i] = (b[i] - (ob.e.center[i]));
      }

      /*rotate only if needed to speed-up*/
      if(((ob.e.theta) == 0) && ((ob.e.phi) == 0)) {
         for(i=0;i<3;++i) {
            ar[i] = ap[i];
            br[i] = bp[i];
         }
      }
      else {  /*rotation needed*/
         ar[1] = -ob.e.sp*ap[0] + ob.e.cp*ap[1];
         ar[0] = -ob.e.st*ap[2] + ob.e.ct*(ob.e.cp*ap[0] + ob.e.sp*ap[1]);
         ar[2] =  ob.e.ct*ap[2] + ob.e.st*(ob.e.cp*ap[0] + ob.e.sp*ap[1]);
         br[1] = -ob.e.sp*bp[0] + ob.e.cp*bp[1];
         br[0] = -ob.e.st*bp[2] + ob.e.ct*(ob.e.cp*bp[0] + ob.e.sp*bp[1]);
         br[2] =  ob.e.ct*bp[2] + ob.e.st*(ob.e.cp*bp[0] + ob.e.sp*bp[1]);
      }

      dels = 0;
      for(i=0;i<3;++i) {
         del[i] = br[i]-ar[i];
         dels += del[i]*del[i];
      }
      dels = sqrt(dels);

      for(i=0;i<3;++i) {
         del[i] /= (ob.e.half_axis[i]);
         ar[i]  /= (ob.e.half_axis[i]);
      }
      /*get the coefficient of the equation for the intersection of the
       line ab with the ellipsoid*/
      c0 = -1;
      c1 = 0;
      c2 = 0;
      for(i=0;i<3;++i) {
         c0 += ar[i]*ar[i];
         c1 += ar[i]*del[i];
         c2 += del[i]*del[i];
      }
      determ = c1*c1 - c0*c2;
      if(determ <= 0) return 0;
      else return ((ob.value)*2*dels*sqrt(determ)/c2);

   } /*end of ellipsoid*/

   else if(ob.geometry == cylinder) 
   {

      for(i=0;i<3;++i) {
         ap[i] = a[i] - ob.c.center[i];
         bp[i] = b[i] - ob.c.center[i];
      }
      /*rotate only if needed to speed-up*/
      if((ob.c.theta) == 0) {
         for(i=0;i<3;++i) {
            ar[i] = ap[i];
            br[i] = bp[i];
         }
      }
      else {  /*rotation needed*/
         ar[1] = -ob.c.sp*ap[0] + ob.c.cp*ap[1];
         ar[0] = -ob.c.st*ap[2] + ob.c.ct*(ob.c.cp*ap[0] + ob.c.sp*ap[1]);
         ar[2] =  ob.c.ct*ap[2] + ob.c.st*(ob.c.cp*ap[0] + ob.c.sp*ap[1]);
         br[1] = -ob.c.sp*bp[0] + ob.c.cp*bp[1];
         br[0] = -ob.c.st*bp[2] + ob.c.ct*(ob.c.cp*bp[0] + ob.c.sp*bp[1]);
         br[2] =  ob.c.ct*bp[2] + ob.c.st*(ob.c.cp*bp[0] + ob.c.sp*bp[1]);
      }
      /*rotated direction coef*/
      for(i=0;i<3;++i) del[i] = br[i]-ar[i];
      /*search intersection with infinite cylinder*/
      c0 = ar[0]*ar[0] + ar[1]*ar[1] - (ob.c.radius)*(ob.c.radius);
      c1 = ar[0]*del[0] + ar[1]*del[1];
      c2 = del[0]*del[0] + del[1]*del[1];
      determ = c1*c1 - c0*c2;
      if(determ > 0) 
        {
          /*intersections with infinite cyl. are a + path(b-a) */
          path1 = (-c1 + sqrt(determ))/c2;
          path2 = (-c1 - sqrt(determ))/c2;
          l2 = (ob.c.length)/2;
          z1 = ar[2] + path1*del[2];
          z2 = ar[2] + path2*del[2];
          if((z1 >= l2) && (z2 >= l2)) return 0;
          if((z1 <= -l2) && (z2 <= -l2)) return 0;
          /*possible intersection with upper and lower basis*/
          if(z1 > l2) path1 = (l2 - ar[2])/del[2];
          if(z2 < -l2) path2 = (-l2 - ar[2])/del[2];
          if(z1 < -l2) path1 = (-l2 - ar[2])/del[2];
          if(z2 > l2) path2 = (l2 - ar[2])/del[2];
          return ((ob.value)*sqrt(c2+del[2]*del[2])*fabs(path2-path1));       
        }
      else if(determ < 0.0) 
        {
          return 0; /*line does not intersect infinite cylinder*/
        }
      else if(c0 < 0)  /*line parallel to axis, within cylinder if c0 < 0*/
        {
          return ((ob.value)*(ob.c.length));
        }
      else return 0;  /*line parallel to axis, outside cylinder*/

   } /*end of cylinder*/

   if (ob.geometry== gaussian)
   {
      /*shift coordinates*/
      for(i=0;i<3;++i) {
         ap[i] = (a[i] - (ob.g.center[i]));
         bp[i] = (b[i] - (ob.g.center[i]));
      }

      /*rotate only if needed to speed-up*/
      if(((ob.g.theta) == 0) && ((ob.g.phi) == 0)) {
         for(i=0;i<3;++i) {
            ar[i] = ap[i];
            br[i] = bp[i];
         }
      }
      else {  /*rotation needed*/
         ar[1] = -ob.g.sp*ap[0] + ob.g.cp*ap[1];
         ar[0] = -ob.g.st*ap[2] + ob.g.ct*(ob.g.cp*ap[0] + ob.g.sp*ap[1]);
         ar[2] =  ob.g.ct*ap[2] + ob.g.st*(ob.g.cp*ap[0] + ob.g.sp*ap[1]);
         br[1] = -ob.g.sp*bp[0] + ob.g.cp*bp[1];
         br[0] = -ob.g.st*bp[2] + ob.g.ct*(ob.g.cp*bp[0] + ob.g.sp*bp[1]);
         br[2] =  ob.g.ct*bp[2] + ob.g.st*(ob.g.cp*bp[0] + ob.g.sp*bp[1]);
      }

      /*unit vector along line*/
      dels = 0;
      for(i=0;i<3;++i) {
         del[i] = br[i]-ar[i];
         dels += del[i]*del[i];
      }
      dels = sqrt(dels);
      for(i=0;i<3;++i) { 
	del[i] /= (dels*ob.g.sigma[i]);
	br[i]  /= (ob.g.sigma[i]);
      }

      c0 = 0.0;
      c1 = 0.0;
      c2 = 0.0;
      for(i=0;i<3;++i) {
         c0 += del[i]*del[i];
         c1 += del[i]*br[i]; 
         c2 += br[i]*br[i];
      }
      return ((ob.value)*sqrt(2.0*PI/c0)*exp(((c1*c1/c0) - c2)/2.0));

   } /*end of gaussian*/


   else return 0;
}

/*---------------------------------------------------------------------------*/
float phantom_line_integral(a,b,phantom,number)



VECTOR a,b;
OBJECT phantom[];
int number;

{
float sum;
int l;

sum = 0;
for(l=0;l<number;++l) sum += line_integral(a,b,phantom[l]);
return sum;

}
/*---------------------------------------------------------------------------*/


int within(a,ob)

/*return 1 if a within ob, and 0 otherwise */

VECTOR a;
OBJECT ob;

{
VECTOR ar,ap;
float rr;
int i;

if(ob.geometry == ellipsoid)
  {
     for(i=0;i<3;++i) ap[i] = (a[i] - (ob.e.center[i]));
      /*rotate*/
     ar[1] = -ob.e.sp*ap[0] + ob.e.cp*ap[1];
     ar[0] = -ob.e.st*ap[2] + ob.e.ct*(ob.e.cp*ap[0] + ob.e.sp*ap[1]);
     ar[2] =  ob.e.ct*ap[2] + ob.e.st*(ob.e.cp*ap[0] + ob.e.sp*ap[1]);
     rr= 0 ;
     for(i=0;i<3;++i)
        rr += ar[i]*ar[i]/((ob.e.half_axis[i])*(ob.e.half_axis[i]));
     if(rr <= 1) return 1;
     else return 0;
  }
else if(ob.geometry == cylinder)
  {
      for(i=0;i<3;++i) {
         ap[i] = a[i] - ob.c.center[i];
      }
      /*rotate*/
      ar[1] = -ob.c.sp*ap[0] + ob.c.cp*ap[1];
      ar[0] = -ob.c.st*ap[2] + ob.c.ct*(ob.c.cp*ap[0] + ob.c.sp*ap[1]);
      ar[2] =  ob.c.ct*ap[2] + ob.c.st*(ob.c.cp*ap[0] + ob.c.sp*ap[1]);
      rr = ar[0]*ar[0] + ar[1]*ar[1];
      if((rr <= (ob.c.radius)*(ob.c.radius)) 
          && (fabs(ar[2]) <= ((ob.c.length)/2))) return 1;
      else return 0;

  }
else if(ob.geometry == gaussian)
  {
     printf("Warning : within function meaningless for gaussians\n");
     return 0;
  }

else return 0;

}

/*---------------------------------------------------------------------------*/

read_phantom(filename,phantom,num_objects)

char *filename;
OBJECT phantom[];
int *num_objects;


{
FILE *f;
double x,y,z;
int k,l;
char s[80];

float scl_x = 0.5*(x_pixels-1)*x_size; 
float scl_y = 0.5*(y_pixels-1)*y_size; 
float scl_z = 0.5*(z_pixels-1)*z_size; 

*num_objects = 0;
l=-1; 
f = fopen(filename,"r");
if(f == NULL) {
   printf("File not opened\n");
   return -1;
}
else printf("File %s opened\n\n",filename);

fscanf(f,"%s",s);
while((strcmp(s,"phantom")) != 0) fscanf(f,"%s",s);

/* read one object*/
while((feof(f) == 0) && ((*num_objects) < MAX_OBJECTS)) {    
   while((strcmp(s,"object") != 0) && (strcmp(s,"end") != 0))   
      fscanf(f,"%s",s);
   if((strcmp(s,"end") == 0)) return 0;
   fscanf(f,"%s",s); 
   k = atoi(s);
/*   printf("Object nr %d ",k); just as a reference, not a sequential number*/
   fscanf(f,"%s",s);
   *num_objects += 1;
   l += 1;
   if((strcmp(s,"ellipsoid") == 0)) {
      phantom[l].geometry = ellipsoid;
      fscanf(f,"%s",s);
      if((strcmp(s,"center") != 0)) return -1;
      fscanf(f,"%s",s); phantom[l].e.center[0] = scl_x * (float) atof(s);
      fscanf(f,"%s",s); phantom[l].e.center[1] = scl_y * (float) atof(s);
      fscanf(f,"%s",s); phantom[l].e.center[2] = scl_z * (float) atof(s);
      fscanf(f,"%s",s);
      if((strcmp(s,"half-axis") != 0)) {
         printf("Error string : %s\n",s);
         return -1;
      }
      fscanf(f,"%s",s); phantom[l].e.half_axis[0] = scl_x * (float) atof(s);
      fscanf(f,"%s",s); phantom[l].e.half_axis[1] = scl_y * (float) atof(s);
      fscanf(f,"%s",s); phantom[l].e.half_axis[2] = scl_z * (float) atof(s);
      fscanf(f,"%s",s);
      if((strcmp(s,"angles") != 0)) return -1;
      fscanf(f,"%s",s); phantom[l].e.theta = (float) atof(s);
      fscanf(f,"%s",s); phantom[l].e.phi = (float) atof(s);
      fscanf(f,"%s",s); phantom[l].value = TRTT_WATER_ABSORPTION * (float) atof(s);
      /*evaluate trigonometric functions*/
      phantom[l].e.ct = cos(degree_to_rad*(phantom[l].e.theta));
      phantom[l].e.st = sin(degree_to_rad*(phantom[l].e.theta));
      phantom[l].e.cp = cos(degree_to_rad*(phantom[l].e.phi));
      phantom[l].e.sp = sin(degree_to_rad*(phantom[l].e.phi));
   }

   else if((strcmp(s,"cylinder") == 0)) {
      phantom[l].geometry = cylinder;
      fscanf(f,"%s",s);
      if((strcmp(s,"radius") != 0)) {
         printf("Error string (radius): %s\n",s);
         return -1;
      }
      fscanf(f,"%s",s); phantom[l].c.radius = (float) atof(s);
      fscanf(f,"%s",s);
      if((strcmp(s,"length") != 0)) {
         printf("Error string (length): %s\n",s);
         return -1;
      }
      fscanf(f,"%s",s); phantom[l].c.length = (float) atof(s);
      fscanf(f,"%s",s);
      if((strcmp(s,"center") != 0)) return -1;
      fscanf(f,"%s",s); phantom[l].c.center[0] = scl_x * (float) atof(s);
      fscanf(f,"%s",s); phantom[l].c.center[1] = scl_y * (float) atof(s);
      fscanf(f,"%s",s); phantom[l].c.center[2] = scl_z * (float) atof(s);
      fscanf(f,"%s",s);
      if((strcmp(s,"angles") != 0)) return -1;
      fscanf(f,"%s",s); phantom[l].c.theta = (float) atof(s);
      fscanf(f,"%s",s); phantom[l].c.phi = (float) atof(s);
      fscanf(f,"%s",s); phantom[l].value = TRTT_WATER_ABSORPTION * (float) atof(s);
      /*evaluate trigonometric functions*/
      phantom[l].c.ct = cos(degree_to_rad*(phantom[l].c.theta));
      phantom[l].c.st = sin(degree_to_rad*(phantom[l].c.theta));
      phantom[l].c.cp = cos(degree_to_rad*(phantom[l].c.phi));
      phantom[l].c.sp = sin(degree_to_rad*(phantom[l].c.phi));

   }
   else if((strcmp(s,"gaussian") == 0)) {
      phantom[l].geometry = gaussian;
      fscanf(f,"%s",s);
      if((strcmp(s,"center") != 0)) return -1;
      fscanf(f,"%s",s); phantom[l].g.center[0] = scl_x * (float) atof(s);
      fscanf(f,"%s",s); phantom[l].g.center[1] = scl_y * (float) atof(s);
      fscanf(f,"%s",s); phantom[l].g.center[2] = scl_z * (float) atof(s);
      fscanf(f,"%s",s);
      if((strcmp(s,"sigma") != 0)) {
         printf("Error string : %s\n",s);
         return -1;
      }
      fscanf(f,"%s",s); phantom[l].g.sigma[0] = (float) atof(s);
      fscanf(f,"%s",s); phantom[l].g.sigma[1] = (float) atof(s);
      fscanf(f,"%s",s); phantom[l].g.sigma[2] = (float) atof(s);
      fscanf(f,"%s",s);
      if((strcmp(s,"angles") != 0)) return -1;
      fscanf(f,"%s",s); phantom[l].g.theta = (float) atof(s);
      fscanf(f,"%s",s); phantom[l].g.phi = (float) atof(s);
      fscanf(f,"%s",s); phantom[l].value = TRTT_WATER_ABSORPTION * (float) atof(s);
      /*evaluate trigonometric functions*/
      phantom[l].g.ct = cos(degree_to_rad*(phantom[l].g.theta));
      phantom[l].g.st = sin(degree_to_rad*(phantom[l].g.theta));
      phantom[l].g.cp = cos(degree_to_rad*(phantom[l].g.phi));
      phantom[l].g.sp = sin(degree_to_rad*(phantom[l].g.phi));
   }
   else {
      printf("Unknown object type : %s\n",s);
      return -1;
   }


}

fclose(f);

}

/*---------------------------------------------------------------------------*/

display_phantom_data(phantom,num_objects)

OBJECT phantom[];
int num_objects;


{
int l;

printf("PHANTOM DEFINITION :\n\n\n");
for(l=0;l < num_objects;++l) {
   if(phantom[l].geometry == ellipsoid) {
      printf("  Object %d is an ellipsoid\n\n",(l+1));
      printf("     Value : %f\n",phantom[l].value);
      printf("     Center : %f , %f , %f\n",phantom[l].e.center[0],
                                            phantom[l].e.center[1],
                                            phantom[l].e.center[2]);
      printf("     Half axis : %f , %f , %f\n",phantom[l].e.half_axis[0],
                                               phantom[l].e.half_axis[1],
                                               phantom[l].e.half_axis[2]);  
      printf("     Polar angle of the third axis (degree) : %f\n",phantom[l].e.theta);
      printf("     Azimuth     of the third axis (degree) : %f\n\n\n",phantom[l].e.phi); 
   }
   else if(phantom[l].geometry == cylinder) {
      printf("  Object %d is a cylinder\n\n",(l+1));
      printf("     Value : %f\n",phantom[l].value);
      printf("     Center : %f , %f , %f\n",phantom[l].c.center[0],
                                            phantom[l].c.center[1],
                                            phantom[l].c.center[2]);
      printf("     Radius : %f\n",phantom[l].c.radius);
      printf("     Length : %f\n",phantom[l].c.length);
      printf("     Polar angle of the axis (degree) : %f\n",phantom[l].c.theta);
      printf("     Azimuth     of the axis (degree) : %f\n\n\n",phantom[l].c.phi);

   }
   else if(phantom[l].geometry == gaussian) {
      printf("  Object %d is a gaussian\n\n",(l+1));
      printf("     Value : %f\n",phantom[l].value);
      printf("     Center : %f , %f , %f\n",phantom[l].g.center[0],
                                            phantom[l].g.center[1],
                                            phantom[l].g.center[2]);
      printf("     standard dev. : %f , %f , %f\n",phantom[l].g.sigma[0],
                                                   phantom[l].g.sigma[1],
                                                   phantom[l].g.sigma[2]);  
      printf("     Polar angle of the third axis (degree) : %f\n",phantom[l].g.theta);
      printf("     Azimuth     of the third axis (degree) : %f\n\n\n",phantom[l].g.phi); 
   }
   else {
      printf("Unknown object in display_phantom_data");
   }
}

}


void  balls_center_projection(a,b,phantom,number,det,Tproj)
VECTOR a,b;
OBJECT phantom[];
int number;
DETECTOR det;
float *Tproj;
{
  int k;
  VECTOR ac,ap,p,bp;
  float l,l1,l2,u,v,eu,ev,n;
  int i;

  for (k=0;k<number;k++){

    l1=0;
    l2=0;
   
    for(i=0;i<3;i++){
      ac[i]=det.center[i]-a[i];
      ap[i]=phantom[k].e.center[i]-a[i];
      l1+=ac[i]*det.normal[i];
      l2+=ap[i]*det.normal[i];      
    }
    l=l1/l2;
    
    u=0;
    v=0;
    
    for(i=0;i<3;i++){
      p[i]=l*ap[i]+a[i];
      bp[i]=p[i]-b[i];
      u+=bp[i]*det.u_axis[i];
      v+=bp[i]*det.v_axis[i];      
    }
    
    //n=(float) rand();
    //eu=n/RAND_MAX-0.5;
    //n=(float) rand();
    //ev=n/RAND_MAX-0.5;
    //printf(" %.1f",eu);
    //Tproj[2*k]=u/u_size+eu;
    //Tproj[2*k+1]=v/v_size+ev;
    Tproj[2*k]=u/u_size;
    Tproj[2*k+1]=v/v_size;
    if ((k%2)==0) printf("\n");
    printf("E%d(%.9f,%.9f),  ",k,u,v);
  }
}

void  Err_BallCentProj_EllCentr(a,b,phantom,number,det,Terr,radius)
VECTOR a,b;
OBJECT phantom[];
int number;
DETECTOR det;
float *Terr;
float radius;
{
  int k;
  VECTOR ac,ap,p,bp;
  float l,l1,l2,u,v;
  int i;
  float magnification;

  for (k=0;k<number;k++){

    l1=0;
    l2=0;

    for(i=0;i<3;i++){
      ac[i]=det.center[i]-a[i];
      ap[i]=phantom[k].e.center[i]-a[i];
      l1+=ac[i]*det.normal[i];
      l2+=ap[i]*det.normal[i];
    }
    l=l1/l2;
    
    u=0;
    v=0;
    for(i=0;i<3;i++){
      p[i]=l*ap[i]+a[i];
      bp[i]=p[i]-b[i];
      u+=bp[i]*det.u_axis[i];
      v+=bp[i]*det.v_axis[i];
    }    
    u/=u_size;
    v/=v_size;

    magnification=(radius*radius)/(l2*l2-radius*radius);

    Terr[2*k]=(u-det.u_off)*magnification;
    Terr[2*k+1]=(v-det.v_off)*magnification;    
    //if ((k%2)==0) printf("\n");
    //printf("\t\tE%d(%.3e),  ",k,l2);
  }
}
