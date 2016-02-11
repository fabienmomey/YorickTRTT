#include <stdio.h>
#include <math.h>
#include "cone_defs.h"
#include "cone_types.h"
#include "simul.h"

main(argc,argv)

     int argc;
     char *argv[];

/*Simulates cone-beam data with motion*/
{
  float phantom_line_integral();	/*in file simul.c*/

  char *object_file; 		/*phantom geometric data*/
  OBJECT phantom[MAX_OBJECTS];
  int number;			/*number of objects in phantom*/

  char *data_file;		/*output file with all cone-beam proj*/
  int u_pixels_over = (int)(0.1*u_pixels);
  int v_pixels_over = (int)(0.1*v_pixels);
  int uover, vover;
  float *proj;

  char *orbit_file;		/*orbit and detector geometric data*/
  ORBIT orbit;
  int points;

  char *shift_name;		/*file with x,y,z shift (in mm) for each projection*/
  FILE *shift_file;
  char s[80];
  VECTOR shift;
  int shift_mode = 0;		/*flag = 1 if shift present*/

  VECTOR a,b;			/*line integrals from point a to b*/
  VECTOR b0, 			
    delu, 			/*move one detector pixel along u axis*/
    delv;		

  int k,i,u,v;
  float max,min,scale, r2;
  int truncate_flag = 0;

  if(argc < 4) {
    printf("sim_cone object_file_name orbit_file_name cone_data_name [shift_name]\n");
    return -1;
  }
  object_file = argv[1];
  orbit_file  = argv[2];
  data_file   = argv[3];

  printf("\nSIMULATE CB DATA %s FOR ORBIT %s FROM PHANTOM %s\n\n",
	 data_file,orbit_file,object_file);

  if(argc > 4) {
    shift_name = argv[4];
    shift_mode = 1;
    printf("\nObject shifts taken from file %s\n\n",shift_name);
    shift_file = fopen(shift_name,"r");
    if(shift_file == NULL) {
      printf("File %s not opened\n",shift_name);
      return -1;
    }
  }

  /* Allocate projections */
  proj = (float*)malloc(u_pixels_over*v_pixels_over*sizeof(float));
  if (proj == NULL) {
    printf("Error in allocating.");
    return -1;
  }
  for (i=0; i<(u_pixels_over*v_pixels_over); ++i)
      *(proj+i) = 0.0;

  /*get from file the geometric data of the phantom to be simulated*/
  if( (read_phantom(object_file,phantom,&number)) != 0) {
    printf("Error in read_phantom");
    return -1;
  }
  display_phantom_data(phantom,number);

  /*get from file the geometric characteristics of the orbit and detector*/
  points = read_orbit(orbit_file,orbit);
  if( points == 0) {
    printf("Problems reading the orbit file %s\n",orbit_file);
    return -1;
  }


  /*simulate successively the points cone beam projections*/
  for(k=0;k < points;++k) { 
    
    printf("\nProcessing projection %d",k);

    /* Re-initialize projection array */
    for (i=0; i<(u_pixels_over*v_pixels_over); ++i)
      *(proj+i) = 0.0;

    for(i=0;i<3;++i) {
      a[i] = orbit[k].focus[i];
      delu[i] = u_size*(orbit[k].u_axis[i]);
      delv[i] = v_size*(orbit[k].v_axis[i]);
      /* b0 is the vector of the left lower corner on the cone beam detector*/
      b0[i] = orbit[k].center[i] - orbit[k].u_off*delu[i] 
	- orbit[k].v_off*delv[i];
    }

    /*translate the focus and detector point by - shift */
    if(shift_mode == 1) {
      if(feof(shift_file) != 0) {
	printf("File %s contains less than %d vectors\n",shift_name,(k+1));
	return -1;
      }
      fscanf(shift_file,"%s",s); shift[0] = (float) atof(s);
      fscanf(shift_file,"%s",s); shift[1] = (float) atof(s);
      fscanf(shift_file,"%s",s); shift[2] = (float) atof(s);
      printf("%d %.2f %.2f %.2f\n",k,shift[0],shift[1],shift[2]);
      for(i=0;i<3;++i) {
	a[i]  -= shift[i];
	b0[i] -= shift[i];
      }
    }

    /* calculate one line integral from focus a to point on detector, b*/
    for(u=0;u < u_pixels;++u) {
	uover = u/10;
	for(i=0;i<3;++i) b[i] = b0[i];
	for(v=0;v < v_pixels; ++v) {
	    vover = v/10;
	    *(proj+u_pixels_over*vover+uover) += u_size*v_size*phantom_line_integral(a,b,phantom,number);    
	    /*move one detector pixel along v_axis*/
	    for(i=0;i<3;++i) b[i] += delv[i];  
	}
	/*before starting next line on the detector, increment b0 by one pixel along u*/
	for(i=0;i<3;++i) b0[i] += delu[i];
    }

    

    /* if(truncate_flag == 0) { */
    /*   truncate_flag = test_truncation(proj); */
    /*   if(truncate_flag == 1) { */
    /* 	printf("WARNING : the projection %d is probably truncated\n",k); */
    /* 	printf("No further warning will be issued\n\n"); */
    /*   } */
    /* } */

    /* Ottawa project: set 0 on pixel out of a disk of radius=9 inches */
    /* for(u=0;u < u_pixels;++u) { */
    /*  for(v=0;v < v_pixels; ++v) { */
    
    /* 	r2 = (u-u_offsino)*(u-u_offsino)*u_size*u_size + (v-v_offsino)*(v-v_offsino)*v_size*v_size; */
    /* 	if (r2>107.5*107.5) */
    /* 	  proj[u][v] = 0;	 */
    /*  }       */
    /*  }  */
    
    /*write one cone-beam projection */
    min = 0;
    max = 0;
    if( (write_image(data_file,proj,(u_pixels_over*v_pixels_over),k,&min,&max,points,v_pixels_over)) != 0) {
      printf("Error writing cone beam projection %d\n",k);
      return -1;
    }

  } /*next cone beam projection */
  
  /* Free memory */
  free(proj);

  return 0;

}

int test_truncation(proj)

     PROJECTION proj;

{
  int u,v;

  /*test whether the edge lines of the projection are zero*/

  for(u=0;u<u_pixels;++u) {
    if(proj[u][0] != 0) return 1;
    if(proj[u][v_pixels-1] != 0) return 1;
  }
  for(v=0;v<v_pixels;++v) {
    if(proj[0][v] != 0) return 1;
    if(proj[u_pixels-1][0] != 0) return 1;
  }
  return 0;
}
