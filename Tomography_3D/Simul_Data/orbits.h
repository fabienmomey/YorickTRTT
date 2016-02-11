/*Definition of constants for the orbit definitions*/

/*dual orbit*/
#define		dual_points	200
#define 	dual_radius     350.0
#define 	dual_focal	700.0

/*CB single orbit */

#define		cb_radius	100.0
#define		cb_focal        153.6
#define		cb_points       720


/*Lines orbit. Along x=[rline*cos(phi),rline*sin(phi)] from z=-hline to +hline
where phi=(k*2*pi/line_number), k=0 ... line_number-1 */ 

#define         line_number     3
#define		rline		300.0
#define 	hline1		208.0
#define 	hline2		120.0
#define 	hline3		100.0
#define		line_points	80
#define		line_focal	300.0 

/*helix */

#define         helix_points    18
#define         helix_radius    450.0
#define         helix_focal     600.0
#define         hhelix          50.0

/*sine */

#define         sine_points    104
#define         sine_radius    600.0
#define         sine_focal     1045.0
#define         hsine          0.0

/*curl*/

#define         curlline_points     20
#define         curlcircle_points   70
#define         curl_radius        300.0
#define         curl_focal         300.0
#define         hcurl              100.0
#define         object_radius      100.0

/*cyl orbit*/

#define         N_circles     5
#define         cyl_points    50
#define         cyl_radius    400.0
#define         cyl_focal     590.0
#define         cyl_h         100.0
