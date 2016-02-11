/* Definition of data types for cone-beam reconstruction */
typedef float VECTOR[4];
typedef float FILTER[dlinopixels];
typedef float FILTER_CB[u_pixels2];
typedef float SLICE[y_pixels][x_pixels];
typedef SLICE IMAGE[z_pixels];
typedef float PROJECTION[u_pixels][v_pixels];
typedef float LARGE_PROJ[u_l][v_l];
typedef float LINOGRAM[dlinoviews][linopixels];
typedef float PADDED_HALF_LINO[linoviews][dlinopixels];
typedef struct
	{
         int discontinuity;
         int neighbour; 
	 float smoothing;
	 VECTOR focus;
         VECTOR tangent;
	 float  focal_length;
	 VECTOR normal; /*normal unit vector to detector*/
         VECTOR u_axis;
         VECTOR v_axis;
	 /* if the projection matrix was not bounded, the pixel
            u_off, v_off would be at the origin of the detector
            coordinate system*/
         float  u_off;  
         float  v_off;
	 VECTOR center; /*centre of the detector coordinate system*/
         int    mindex; /*index of the M file record to be used*/
	} DETECTOR;
typedef DETECTOR ORBIT[MAX_ORB_POINTS];

