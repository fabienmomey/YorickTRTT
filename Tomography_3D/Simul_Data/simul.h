
/*version December 1996 includes gaussian*/

enum GEOMETRY{cylinder, ellipsoid, gaussian};
typedef struct {
   VECTOR 	center;
   VECTOR 	half_axis;
   float 	theta,phi;   /*axis orientation*/
   float 	ct,st,cp,sp; /*cos and sin of theta and phi*/
} ELLIPSOID;
typedef struct {
   float 	radius;
   float 	length;
   VECTOR 	center;
   float 	theta,phi;   /*axis orientation*/
   float 	ct,st,cp,sp; /*cos and sin of theta and phi*/
} CYLINDER;
typedef struct {
   VECTOR 	center;
   VECTOR 	sigma;
   float 	theta,phi;   /*axis orientation*/
   float 	ct,st,cp,sp; /*cos and sin of theta and phi*/
} GAUSSIAN;
typedef struct  {
   float 	value;
   enum GEOMETRY geometry;
   ELLIPSOID 	e;
   CYLINDER 	c; 
   GAUSSIAN     g; /*only one of the three last is used*/
} OBJECT;

