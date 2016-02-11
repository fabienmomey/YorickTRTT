#include <stdio.h>
#include <math.h>
#include "cone_defs.h"
#include "cone_types.h"
#include "simul.h"
static IMAGE image;
main(argc,argv)

int argc;
char *argv[];

/*Simulates a 3D image from an object file */
{
int within();

char *object_file, *image_file;
OBJECT phantom[MAX_OBJECTS];
int number,n;
float x,y,z;
int i,j,k;
VECTOR a;
float max,min;


if(argc < 2) {
   printf("make_image object_file_name image_file_name\n");
   return -1;
}
object_file = argv[1];
image_file = argv[2];
printf("\nSIMULATE IMAGE %s FROM PHANTOM DESCRIBED IN %s\n\n",
       image_file,object_file);

if( (read_phantom(object_file,phantom,&number)) != 0) {
   printf("Error in read_phantom");
   return -1;
}
display_phantom_data(phantom,number);

min = 0;
max = 0;

for(k=0; k < z_pixels; ++k) {
    a[2] = -(z_off*z_size) + k*z_size;
    for(j=0; j < y_pixels; ++j) {
	a[1] = -(y_off*y_size) + j*y_size;
	for(i=0; i < x_pixels; ++i) { 
	    a[0] = -(x_off*x_size) + i*x_size;
	    image[k][j][i] = 0;
	    for(n=0; n< number;++n) 
		if(within(a,phantom[n]) == 1) image[k][j][i] += (phantom[n].value);
	    if(image[k][j][i] > max) max = image[k][j][i];
	    if(image[k][j][i] < min) min = image[k][j][i];  
	}
    }
}
/*
for(k=0;k<z_pixels;++k)
  for(j=0;j<y_pixels;++j)
    for(i=0;i<x_pixels;++i){
      if(image[k][j][i]>1.04) image[k][j][i]=1.04;
      if(image[k][j][i]<1.005) image[k][j][i]=1.005;
      printf("%f\n",(image[k][j][i]-1.005)/(1.04-1.005));}
*/
for(k=0;k<z_pixels;++k)
   if( (write_image(image_file,image[k],(x_pixels*y_pixels),k,&min,&max,z_pixels,y_pixels)) != 0) {
     printf("Error writing");
     return 0;
   }

}

