#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include "cone_defs.h"
#include "cone_types.h"
#include "orbits.h"
#include "scan_types.h"
#include <stdlib.h>

#define BUFFER  128
#define BUFFER_SNM 1024

/*routines for input output for the Kingston display.sun4 :
  write_image
  read_image
  routines to read Picker header-less files :
  read_SNM_image
  routines to read, write and print orbit description files :
  read_orbit
  write_orbit
  print_orbit
*/

int write_image(filename,data,number,index,pmin,pmax,iz,iy)

char 	*filename;	/*Output file name*/
float 	data[];		/*Array with the data to be written. 
			  Maybe be any type (PROJECTION, SINOGRAM, SLICE)*/
int	number; 	/*of floating point voxels to be written*/
int 	index;  	/*position in file, in records of "number" 
			  floating p.*/
float	*pmin,*pmax;	/*minimum and maximum*/
int     iz,iy;  	/*iz : number of records the file will contain*/
			/*iy : if one record is a 2D image, iy= im. size*/

{
   FILE 	*f;
   int  	k,first_to_write,left_to_write;
   long    	offset;
   float	a[BUFFER];  /*image written as floating numbers*/
   float	realdata =1.0;
   
   f = fopen(filename,"a+");
   if(f == NULL) return -1;
   /*position in file at index data */
   if(index == 0) offset = 0;
   else
                  offset = (BUFFER + index*number)*sizeof(float);
   /*we use BUFFER + because the header must be skipped */

   if( (fseek(f,offset,SEEK_SET)) != 0) {
      printf("Fseek problem (write_image)\n");
      return -1;
   }
   
   /*evaluate the minimum and maximum IF they have not been initialized*/
   if((*pmax) == 0) {
      *pmin = data[0];
      *pmax = data[0];
      for(k=0;k<number;++k) { 
   	if(*pmin > data[k]) *pmin = data[k];
	if(*pmax < data[k]) *pmax = data[k];
      }
   }
 
   /*write first the header if the first record is to be read*/
   if(index == 0) {
     for(k=0;k<BUFFER;++k)a[k] = 0.0;
      
      a[0] = *pmin;
      a[1] = (*pmax)*(1.01);
      a[2] = 1;
      a[3] = iz;
      a[4] = iy;
      a[5] = 513;
      a[10] = (float) 1.0*x_pixels;
      a[11] = (float) 1.0*y_pixels;
      a[12] = (float) 1.0*z_pixels;
      a[13] = x_size;
      a[14] = y_size;
      a[15] = z_size;
      a[16] = 1.0*x_off;
      a[17] = 1.0*y_off;
      a[18] = 1.0*z_off;
      a[19] = (float) 1.0*u_pixels;
      a[20] = (float) 1.0*v_pixels;
      a[31] = u_size;
      a[32] = v_size;
      a[33] = 1.0*u_offsino;
      a[34] = 1.0*v_offsino;
      a[35] = 1.0*stretch_factor;
      a[43] = dual_points;
      a[44] = dual_radius;
      a[45] = dual_focal;
      a[46] = cb_radius;
      a[47] = cb_focal;
      a[48] = cb_points;

      k = fwrite(a,sizeof(float),BUFFER,f);
      if(k != BUFFER) {
        printf("write_image : %d voxels written out of %d\n",k,BUFFER);
        return -1;
      }
   }   

   
   /* write buffers of data */
   left_to_write = number;
   first_to_write = 0;

   while(left_to_write >= BUFFER) {
      for(k=0;k<BUFFER;++k)
	  a[k] = data[first_to_write + k];
      k = fwrite(a,sizeof(float),BUFFER,f);
      if(k != BUFFER) {
	  printf("write_image : %d voxels written out of %d\n",k,BUFFER);
	   return -1;
      }
      first_to_write += BUFFER;
      left_to_write -= BUFFER;
   }
   /*last bit if necessary*/
   if(left_to_write > 0) {
      for(k=0;k<left_to_write;++k)
         a[k] = data[first_to_write + k];
      k = fwrite(a,sizeof(float),left_to_write,f);
      if(k != left_to_write) {
	  printf("write_image : %d voxels written out of %d\n",k,left_to_write);
	  return -1;
      }
   }
      
   fclose(f);
   return 0;

}

int read_image(filename,data,number,index,pmin,pmax,piz,piy)

char 	*filename;
float 	data[];
int	number; /*of floating point voxels to be read*/
int	index;  /*position in file, in number of images*/
float   *pmin,*pmax; /*not checked, read from header if index=0 */
int	*piz,*piy;  /*z and y dimensions, read from header if index=0*/

{
   FILE 	*f;
   int  	k,first_to_read,left_to_read;
   long		offset;
   float	a[BUFFER];

   f = fopen(filename,"r");
   if(f == NULL) return -1;

   if(index == 0) offset = 0;
   else
                  offset = (BUFFER + index*number)*sizeof(float);
   /*we use BUFFER + because the header must be skipped */

   if( (fseek(f,offset,SEEK_SET)) != 0) {
      printf("Fseek problem (read_image)\n");
      return -1;
   }

   /*read the header*/
   if(index == 0) {
      k = fread(a,sizeof(float),BUFFER,f);
      if(k != BUFFER) {
   	printf("read_image (header): %d voxels read instead of %d\n",k,BUFFER);
	   return -1;
      }
   *pmin = a[0];
   *pmax = a[1];
   if(a[2] != 1.0) printf("read_image : header flag is not real data\n");
   *piz = (int) a[3];
   *piy = (int) a[4];
   }

   /* read buffers of data */
   left_to_read = number;
   first_to_read = 0;

   while(left_to_read >= BUFFER) {
      k = fread(a,sizeof(float),BUFFER,f);
      if(k != BUFFER) {
   	printf("read_image : %d voxels read instead of %d\n",k,BUFFER);
	   return -1;
      }
      for(k=0;k<BUFFER;++k)
         data[first_to_read + k] = a[k];
      first_to_read += BUFFER;
      left_to_read -= BUFFER;
   }
   /*last bit if necessary*/
   if(left_to_read > 0) {
      k = fread(a,sizeof(float),left_to_read,f);
      if(k != left_to_read) {
   	printf("read_image : %d voxels read instead of %d\n",k,left_to_read);
	   return -1;
      }
      for(k=0;k<left_to_read;++k)
         data[first_to_read + k] = a[k];

   }
      
   fclose(f);
   return 0;

}

int display_header(filename,mode)

char 	*filename;
int 	mode;


{
   FILE 	*f;
   float	a[BUFFER];
   int k;

   f = fopen(filename,"r");
   if(f == NULL) return -1;

   /*read the header*/

   k = fread(a,sizeof(float),BUFFER,f);
   if(k != BUFFER) {
   	printf("read_image (header): %d voxels read instead of %d\n",k,BUFFER);
	return -1;
   }
   if(a[2] != 1.0) printf("read_image : header flag is not real data\n");


printf("Header signals : %d slices (%d pixels)\n",(int) a[3],(int) a[4]);
printf("Minimum : %f, Maximum : %f\n",a[0],a[1]);
if(mode == 1) {
printf("\nx y z pixels		:%d %d %d\n",(int) a[10],(int) a[11],(int) a[12]);
printf("pixel size		:%.2f %.2f %.2f\n",a[13],a[14],a[15]);
printf("origin at pixel		:%.2f %.2f %.2f\n",a[16],a[17],a[18]);
printf("projection pixels	:%d %d\n",(int) a[19],(int) a[20]);
printf("projection pixel size	:%.2f %.2f\n",a[31],a[32]);
printf("projection centre	:%.2f %.2f\n",a[33],a[34]);
printf("stretch factor		:%d\n",(int) a[35]);
printf("sinogram radial pixels	:%d\n",(int) a[36]);
printf("radial stretch factor	:%d\n",(int) a[37]);
printf("sinogram angular views	:%d\n",(int) a[38]);
printf("sinogram pixel size	:%.2f\n",a[39]);
printf("sinogram r centre	:%.2f\n",a[40]);
printf("sinogram interleave	:%.2f\n",a[41]);
printf("cb_sino_offset		:%.2f\n",a[42]);
printf("\nInclude file orbits.h :\n");
printf("Dual orbit : %d points, radius : %.2f, focal length : %.2f\n",(int) a[43],a[44],a[45]);
printf("CB   orbit : %d points, radius : %.2f, focal length : %.2f\n",(int) a[48],a[46],a[47]);
printf("PB   orbit : %d points, radius : %.2f, focal length : %.2f\n",(int) a[51],a[49],a[50]);
printf("PB-CB orbit: %d points, PB weight : %.2f, CB weight : %.2f\n",(int) a[54],a[53],a[52]);
}








}

int read_SNM_image(filename,data,number,index,pmin,pmax,piz,piy)

char 	*filename;
float 	data[];
int	number; /*of floating point voxels to be read*/
int	index;  /*position in file, in number of images*/
float   *pmin,*pmax; /*not checked, read from header*/
int	*piz,*piy;  /*z and y dimensions, only read from header and returned*/

{
   FILE 	*f;
   int  	k,first_to_read,left_to_read;
   long		offset;
   short int	a[BUFFER_SNM];

   f = fopen(filename,"r");
   if(f == NULL) return -1;

   if(index == 0) offset = 0;
   else
                  offset = (BUFFER_SNM + index*number)*sizeof(short int);
   /*we use BUFFER_SNM + because the 2K header must be skipped */



   if( (fseek(f,offset,SEEK_SET)) != 0) {
      printf("Fseek problem (read_image)\n");
      return -1;
   }

   /*read the header*/
   if(index == 0) {
      k = fread(a,sizeof(short int),BUFFER_SNM,f);
      if(k != BUFFER_SNM) {
   	printf("read_image (header): %d voxels read instead of %d\n",k,BUFFER);
	   return -1;
      }
   }

   /* read buffers of data */
   left_to_read = number;
   first_to_read = 0;

   while(left_to_read >= BUFFER_SNM) {
      k = fread(a,sizeof(short int),BUFFER_SNM,f);
      if(k != BUFFER_SNM) {
   	printf("read_image : %d voxels read instead of %d\n",k,BUFFER_SNM);
	   return -1;
      }
      for(k=0;k<BUFFER_SNM;++k)
         data[first_to_read + k] = a[k];
      first_to_read += BUFFER_SNM;
      left_to_read -= BUFFER_SNM;
   }
   /*last bit if necessary*/
   if(left_to_read > 0) {
      k = fread(a,sizeof(short int),left_to_read,f);
      if(k != left_to_read) {
   	printf("read_image : %d voxels read instead of %d\n",k,left_to_read);
	   return -1;
      }
      for(k=0;k<left_to_read;++k)
         data[first_to_read + k] = a[k];

   }
      
   fclose(f);
   return 0;

}




int write_orbit(filename,orbit,points)



char 	*filename;
ORBIT   orbit;
int 	points;

{
   FILE 	*f;
   int  	k;
 
   f = fopen(filename,"w");
   if(f == NULL) {
      printf("Could not open file %s (write_orbit)\n",filename);
      return -1;
   }
   if( (fwrite(orbit,sizeof(DETECTOR),points,f)) != points) {
      printf("Error in write_orbit\n");
      return -1;
   }
   fclose(f);
   return 0;

}
   

int read_orbit(filename,orbit)

/*returns the number of orbit points read*/

char 	*filename;
ORBIT   orbit;


{
   FILE 	*f;
   int  	k,i,points;
   DETECTOR     d;
   float        eps = 0.0001;
   float        nn,vv,uu,uv,nv,nu,tu,tv,r2;     
 
   f = fopen(filename,"r");
   if(f == NULL) {
      printf("Could not open file %s (read_orbit)\n",filename);
      return -1;
   }
   
   k=0;
   points = 0;
   i = 1;

   while((!feof(f)) && (i==1)) {
      i = fread(&orbit[points],sizeof(DETECTOR),1,f);
      if(i == 1)  ++ points;
   }
   printf("The orbit contains %d points\n",points);
   fclose(f);


   /*now check compatibility*/
   for(k=0;k<points;++k) {
      d = orbit[k];
      nn = -1;
      uu = -1;
      vv = -1;
      nv = 0;
      uv = 0;
      nu = 0;
      tu = 0;
      tv = 0;
      r2 = 0;
      for(i=0;i<3;++i) {
      nn += (d.normal[i])*(d.normal[i]);
      uu += (d.u_axis[i])*(d.u_axis[i]);
      vv += (d.v_axis[i])*(d.v_axis[i]);
      nv += (d.normal[i])*(d.v_axis[i]);
      nu += (d.normal[i])*(d.u_axis[i]);
      uv += (d.u_axis[i])*(d.v_axis[i]);
      tu += (d.u_axis[i])*(d.center[i] - d.focus[i]);
      tv += (d.v_axis[i])*(d.center[i] - d.focus[i]);
      r2 += (d.center[i] - d.focus[i])*(d.center[i] - d.focus[i]);
      }
      tu /= sqrt(r2);
      tv /= sqrt(r2);
      r2 = 1.0 - (r2/((d.focal_length)*(d.focal_length)));

      if(fabs(nn) > eps) printf("Incompatibility nn (read_orbit), k = %d\n",k);
      if(fabs(uu) > eps) printf("Incompatibility uu (read_orbit), k = %d\n",k);
      if(fabs(vv) > eps) printf("Incompatibility vv (read_orbit), k = %d\n",k);
      if(fabs(nv) > eps) printf("Incompatibility nv (read_orbit), k = %d\n",k);
      if(fabs(nu) > eps) printf("Incompatibility nu (read_orbit), k = %d\n",k);
      if(fabs(uv) > eps) printf("Incompatibility uv (read_orbit), k = %d\n",k);
      if(fabs(tu) > eps) printf("Incompatibility tu (read_orbit), k = %d\n",k);
      if(fabs(tv) > eps) printf("Incompatibility tv (read_orbit), k = %d\n",k);
      if(fabs(r2) > eps) printf("Incompatibility r2 (read_orbit), k = %d\n",k);
     
   }

   return points;

}
   

int print_orbit(orbit,k)

ORBIT orbit;
int   k;

{
DETECTOR d;

if( (k < 0) || (k >= MAX_ORB_POINTS)) {
   printf("Invalid k (print_orbit)\n");
   return -1;
}

d = orbit[k];

printf("\nDetector position number %d :\n",k);
printf("Focal length (mm) : %f\n",d.focal_length);
printf("Focal point        (mm):  %f  %f  %f\n",d.focus[0],d.focus[1],d.focus[2]);
printf("Normal to detector (mm):  %f  %f  %f\n",d.normal[0],d.normal[1],d.normal[2]);
printf("u axis             (mm):  %f  %f  %f\n",d.u_axis[0],d.u_axis[1],d.u_axis[2]);
printf("v axis             (mm):  %f  %f  %f\n",d.v_axis[0],d.v_axis[1],d.v_axis[2]);
printf("Center of detector (mm):  %f  %f  %f\n",d.center[0],d.center[1],d.center[2]);

}

Nnoise_v2(noise,phi,Tf,Ta,Tb,N)
     double *noise;
     double phi;
     double *Tf;
     double *Ta;
     double *Tb;
     int N;
{
  int i;

  *noise=0;
  for(i=0;i<N;i++)
    *noise +=Ta[i]*sin(Tf[i]*phi)+Tb[i]*cos(Tf[i]*phi);
}

Nnoise(noise)
     float *noise;
{
  *noise=(float)rand()*2/RAND_MAX -1;
}

ParIntExt_ParScan(orbit,Npts,Tscan)
ORBIT orbit;
int   Npts;
TPARSCAN Tscan;

{

  int k;
  VECTOR CH;
  float noise,a;
  int test_noise;
  double phi,eta_d,theta_d, phi_d,dphi_d;
  double Ceta,Seta,Ctheta,Stheta,Cphi,Sphi;
  
  
  for(k=0;k<Npts;k++){
    Tscan[k].focus[0]=orbit[k].focus[0];
    Tscan[k].focus[1]=orbit[k].focus[1];
    Tscan[k].focus[2]=orbit[k].focus[2];

    CH[0]=(orbit[k].u_off-u_pixels/2.)*orbit[k].u_axis[0] +
      (orbit[k].v_off-v_pixels/2.)*orbit[k].v_axis[0];
    CH[1]=(orbit[k].u_off-u_pixels/2.)*orbit[k].u_axis[1] +
      (orbit[k].v_off-v_pixels/2.)*orbit[k].v_axis[1];
    CH[2]=(orbit[k].u_off-u_pixels/2.)*orbit[k].u_axis[2] +
      (orbit[k].v_off-v_pixels/2.)*orbit[k].v_axis[2];

    Tscan[k].Dcenter[0]=orbit[k].center[0]-CH[0];
    Tscan[k].Dcenter[1]=orbit[k].center[1]-CH[1];
    Tscan[k].Dcenter[2]=orbit[k].center[2]-CH[2];

    Tscan[k].u_axis[0]=orbit[k].u_axis[0];
    Tscan[k].u_axis[1]=orbit[k].u_axis[1];
    Tscan[k].u_axis[2]=orbit[k].u_axis[2];

    Tscan[k].v_axis[0]=orbit[k].v_axis[0];
    Tscan[k].v_axis[1]=orbit[k].v_axis[1];
    Tscan[k].v_axis[2]=orbit[k].v_axis[2];

    Tscan[k].normal[0]=orbit[k].normal[0];
    Tscan[k].normal[1]=orbit[k].normal[1];
    Tscan[k].normal[2]=orbit[k].normal[2];

    
    //*********************************************************************
    //Perturbation de la trajectoire
    //*********************************************************************
    test_noise=1;   

    if(test_noise==1){

      //Noise on source position
      Nnoise(&noise); a=0.025*noise; a=0; /* 0.05 */     
      Tscan[k].focus[0] += a*Tscan[k].focus[0];
      Nnoise(&noise); a=0.025*noise; a=0; /* 0.05 */
      Tscan[k].focus[1] += a*Tscan[k].focus[1];
      Nnoise(&noise); a=0.3*noise; a=0;    /* 5 */
      Tscan[k].focus[2] += a;

      //Noise on detector orientation
      //phi=atan2(Tscan[k].normal[1],Tscan[k].normal[0]);
      phi = 2*PI*k/Npts + PI/Npts;
      Nnoise(&noise); noise=0; theta_d=2.5*PI/180.0*noise; /* 2.5 */
      Nnoise(&noise); if (k%2==0) noise=1; else noise=0;
      eta_d=22*PI/180.0*noise; /* 2.5 */
      Nnoise(&noise); noise=0; dphi_d=2.5*PI/180.0*noise; /* 2.5 */
      
      phi_d=phi+dphi_d;
      Ctheta=cos(theta_d);Stheta=sin(theta_d);
      Ceta=cos(eta_d);Seta=sin(eta_d);
      Cphi=cos(phi_d);Sphi=sin(phi_d);
      Tscan[k].normal[0] = Cphi*Ctheta;
      Tscan[k].normal[1] = Sphi*Ctheta;
      Tscan[k].normal[2] = Stheta;
      Tscan[k].u_axis[0] = -Ceta*Sphi-Seta*Cphi*Stheta;
      Tscan[k].u_axis[1] = Ceta*Cphi-Seta*Sphi*Stheta;
      Tscan[k].u_axis[2] = Seta*Ctheta;
      Tscan[k].v_axis[0] = Seta*Sphi-Ceta*Cphi*Stheta;
      Tscan[k].v_axis[1] = -Seta*Cphi -Ceta*Sphi*Stheta;
      Tscan[k].v_axis[2] = Ceta*Ctheta;

      //Noise on detector center
      Nnoise(&noise); a=0.05*noise; a=0;  /* 0.05 */
      Tscan[k].Dcenter[0] += a*Tscan[k].Dcenter[0];
      Nnoise(&noise); a=0.05*noise; a=0;  /* 0.05 */
      Tscan[k].Dcenter[1] += a*Tscan[k].Dcenter[1];
      Nnoise(&noise); a=5*noise; a=0; /* 0.05 */
      Tscan[k].Dcenter[2] += a;
    }    
  }
}


ParIntExt_ParScan_v2(orbit,Npts,Tscan,rayon)
ORBIT orbit;
int   Npts;
TPARSCAN Tscan;
double rayon;

{

  int k,i;
  VECTOR CH;
  double noise,a,Ta[10],Tb[10];
  int test_noise;
  double phi,eta_d,theta_d, phi_d,dphi_d;
  double Ceta,Seta,Ctheta,Stheta,Cphi,Sphi;
  double freq,rayon_p,Tfreq[10];
  
  
  for(k=0;k<Npts;k++){
    Tscan[k].focus[0]=orbit[k].focus[0];
    Tscan[k].focus[1]=orbit[k].focus[1];
    Tscan[k].focus[2]=orbit[k].focus[2];

    CH[0]=(orbit[k].u_off-u_pixels/2.)*orbit[k].u_axis[0] +
      (orbit[k].v_off-v_pixels/2.)*orbit[k].v_axis[0];
    CH[1]=(orbit[k].u_off-u_pixels/2.)*orbit[k].u_axis[1] +
      (orbit[k].v_off-v_pixels/2.)*orbit[k].v_axis[1];
    CH[2]=(orbit[k].u_off-u_pixels/2.)*orbit[k].u_axis[2] +
      (orbit[k].v_off-v_pixels/2.)*orbit[k].v_axis[2];

    Tscan[k].Dcenter[0]=orbit[k].center[0]-CH[0];
    Tscan[k].Dcenter[1]=orbit[k].center[1]-CH[1];
    Tscan[k].Dcenter[2]=orbit[k].center[2]-CH[2];

    Tscan[k].u_axis[0]=orbit[k].u_axis[0];
    Tscan[k].u_axis[1]=orbit[k].u_axis[1];
    Tscan[k].u_axis[2]=orbit[k].u_axis[2];

    Tscan[k].v_axis[0]=orbit[k].v_axis[0];
    Tscan[k].v_axis[1]=orbit[k].v_axis[1];
    Tscan[k].v_axis[2]=orbit[k].v_axis[2];

    Tscan[k].normal[0]=orbit[k].normal[0];
    Tscan[k].normal[1]=orbit[k].normal[1];
    Tscan[k].normal[2]=orbit[k].normal[2];

    
    //*********************************************************************
    //Perturbation REGULIERE de la trajectoire
    //*********************************************************************
    test_noise=1;   

    if(test_noise==1){
      //phi=atan2(Tscan[k].normal[1],Tscan[k].normal[0]);
      phi = 2*PI*k/Npts + PI/Npts;
      dphi_d=0;
      phi_d=phi+dphi_d;
      Cphi=cos(phi_d);Sphi=sin(phi_d);


      //Noise on source position

      Tfreq[0]=1; Tfreq[1]=2; Tfreq[2]=3;   
      Ta[0]=0.02; Ta[1]=0.005; Ta[2]=0.001;
      Tb[0]=0.01; Tb[1]=0.01; Tb[2]=0.005;
      Nnoise_v2(&noise,phi,Tfreq,Ta,Tb,3);     
      Tscan[k].focus[0] += rayon*Cphi*noise;

      Tfreq[0]=1; Tfreq[1]=2; Tfreq[2]=3;   
      Ta[0]=0.02; Ta[1]=0.001; Ta[2]=0.0;
      Tb[0]=0.005; Tb[1]=0.002; Tb[2]=0.0;
      Nnoise_v2(&noise,phi,Tfreq,Ta,Tb,2); 
      Tscan[k].focus[1] += rayon*Sphi*noise;
      
      
      Tfreq[0]=1; Tfreq[1]=2; 
      Ta[0]=0.03; Ta[1]=0.01;
      Tb[0]=0.01; Tb[1]=0.005;
      Nnoise_v2(&noise,phi,Tfreq,Ta,Tb,2);
      Tscan[k].focus[2] = rayon*noise;
      //printf("\nS(%d)=%f",k,noise);


      //Noise on detector orientation

      freq=2;
      noise=sin(freq*phi);
      theta_d=2.5*PI/180.0*noise; /* 5 */

      freq=1;
      noise=sin(freq*phi);
      eta_d=2.5*PI/180.0*noise; /* 2.5 */

      a=0.02;
      freq=2;
      noise=sin(freq*phi);
      dphi_d=a*2*PI/Npts*noise; /* a*dphi */
      
      Ctheta=cos(theta_d);Stheta=sin(theta_d);
      Ceta=cos(eta_d);Seta=sin(eta_d);      
      phi_d=phi+dphi_d;
      Cphi=cos(phi_d);Sphi=sin(phi_d);

      Tscan[k].normal[0] = Cphi*Ctheta;
      Tscan[k].normal[1] = Sphi*Ctheta;
      Tscan[k].normal[2] = Stheta;
      Tscan[k].u_axis[0] = -Ceta*Sphi-Seta*Cphi*Stheta;
      Tscan[k].u_axis[1] = Ceta*Cphi-Seta*Sphi*Stheta;
      Tscan[k].u_axis[2] = Seta*Ctheta;
      Tscan[k].v_axis[0] = Seta*Sphi-Ceta*Cphi*Stheta;
      Tscan[k].v_axis[1] = -Seta*Cphi -Ceta*Sphi*Stheta;
      Tscan[k].v_axis[2] = Ceta*Ctheta;

      //Noise on detector center
      freq=2;
      noise=sin(freq*phi);
      a=0.05*noise; /* 0.05 */
      Tscan[k].Dcenter[0] += a*Tscan[k].Dcenter[0];
      a=0.05*noise; /* 0.05 */
      Tscan[k].Dcenter[1] += a*Tscan[k].Dcenter[1];
      a=0.05*noise; /* 0.05 */
      Tscan[k].Dcenter[2] += a;
    }    
  }
}

ParScan_ParIntExt(orbit,Npts,Tscan)
ORBIT orbit;
int   Npts;
TPARSCAN Tscan;

{
  int k,i;
  float tmp1,tmp2;
  VECTOR CH;
  
  for(k=0;k<Npts;k++){
    orbit[k].focus[0]=Tscan[k].focus[0];
    orbit[k].focus[1]=Tscan[k].focus[1];
    orbit[k].focus[2]=Tscan[k].focus[2];

    orbit[k].normal[0]=Tscan[k].normal[0];
    orbit[k].normal[1]=Tscan[k].normal[1];
    orbit[k].normal[2]=Tscan[k].normal[2];

    orbit[k].u_axis[0]=Tscan[k].u_axis[0];
    orbit[k].u_axis[1]=Tscan[k].u_axis[1];
    orbit[k].u_axis[2]=Tscan[k].u_axis[2];

    orbit[k].v_axis[0]=Tscan[k].v_axis[0];
    orbit[k].v_axis[1]=Tscan[k].v_axis[1];
    orbit[k].v_axis[2]=Tscan[k].v_axis[2];

    tmp1=0;
    tmp2=0;
    for(i=0;i<3;i++){
      tmp1+=Tscan[k].focus[i]*Tscan[k].normal[i];
      tmp2+=Tscan[k].Dcenter[i]*Tscan[k].normal[i];
    }
    orbit[k].focal_length=tmp1-tmp2;

    orbit[k].center[0]=orbit[k].focus[0]-orbit[k].focal_length*orbit[k].normal[0];
    orbit[k].center[1]=orbit[k].focus[1]-orbit[k].focal_length*orbit[k].normal[1];
    orbit[k].center[2]=orbit[k].focus[2]-orbit[k].focal_length*orbit[k].normal[2];

    CH[0]= orbit[k].center[0]-Tscan[k].Dcenter[0];
    CH[1]= orbit[k].center[1]-Tscan[k].Dcenter[1];
    CH[2]= orbit[k].center[2]-Tscan[k].Dcenter[2];

    tmp1=0;
    tmp2=0;
    for(i=0;i<3;i++){
      tmp1+=CH[i]*Tscan[k].u_axis[i];
      tmp2+=CH[i]*Tscan[k].v_axis[i];
    }
    
    orbit[k].u_off=u_pixels/2. + tmp1*u_size;
    orbit[k].v_off=v_pixels/2. + tmp2*v_size;

  }
}
