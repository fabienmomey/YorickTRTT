include, "Ytrtt2D.i", 1;

// data_dir = "/home/momey/Recherche_Tomographie/Data/data_CLB_09-11-2011/img_1.3.46.423632.135428.1320853170.4/";
// data_dir = "/home/momey/Recherche_Tomographie/Data/data_CLB_09-11-2011/img_1.3.46.423632.135428.1320855641.20/";
// data_dir = "/home/momey/Recherche_Tomographie/Data/data_CLB_09-11-2011/img_1.3.46.423632.135428.1320854043.7/";
// data_dir = "/home/momey/Recherche_Tomographie/Data/data_CLB_09-11-2011/img_1.3.46.423632.135428.1320854260.10/";
// data_dir = "/home/momey/Recherche_Tomographie/Data/data_CLB_patient_02-07-2012/img_1.3.46.423632.141000.1169042526.68/";
// data_dir= "/home/momey/Recherche_Tomographie/Data/data_CLB_09-11-2011/img_1.3.46.423632.135428.1320854260.10/";
// data_dir= "/home/momey/Recherche_Tomographie/Data/data_CLB_09-11-2011/img_1.3.46.423632.135428.1320854585.13/";
// data_dir= "/home/momey/Recherche_Tomographie/Data/data_CLB_09-11-2011/img_1.3.46.423632.135428.1320855182.16/";
// data_dir= "/home/momey/Recherche_Tomographie/Data/data_CLB_09-11-2011/img_1.3.46.423632.135428.1320855641.20/";

// data_dir= "/home/momey/Recherche_Tomographie/Data/tete/IMAGES/img_1.3.46.423632.141000.1263996445.20/";
data_dir= "/home/momey/Recherche_Tomographie/Data/data_CLB_patient_02-07-2012/img_1.3.46.423632.141000.1169042526.68/";


cmd1 = "dir -1 "+data_dir+"*.his > '/tmp/projections_list.txt'";
system, cmd1;

proj_filenames = open("/tmp/projections_list.txt");
projections_bin = rdfile(proj_filenames);
close, proj_filenames;

/* Get the number of projections */
ndata = numberof(projections_bin);
nu=nv=512;

/* Extract projections */
data = array(double,nu,nv,ndata);
for (k=1; k<=ndata; ++k) {
    data(..,k) = transpose(double(raw_read(projections_bin(k),short,512,512,offset=100) + (2.^16-1.0)));
}

/* Lookup table conversion */
/* FIXME: It makes the -log(I/I0) operation

   Raw data is encoded in UNSIGNED SHORT type : 0 -> 65535 (2^16 - 1)

   An increasing value means an increasing attenuation. The I0 intensity (no
   attenuation) is considered to be 2^16 = 65536. Thus the conversion of raw
   pixel values "x" into sinogram data values "q" is made as follows:

   q = log(I0) - log(I0-x)

   Values of 0 and 65535 mean that the pixel is saturated in intensity => these
   values are set at -1, to identify them as invalid data.

 */

I0 = 2.^16;
logI0 = log(I0);

for (k=1; k<=ndata; ++k) {
    local data_k; eq_nocopy, data_k, data(..,k);
    id=where(data_k<=0.0 | data_k>=(I0-1.0));
    id2=where(data_k>0.0 & data_k<(I0-1.0));
    if (is_array(id))
        data_k(id) = -1;
    if (is_array(id2))
        data_k(id2) = logI0 - log(I0-data_k(id2));
    mywindow, 0; fma; img_plot, data_k, cbar=1, vert=1, format="%.3f";
    pause, 1;
    data(..,k) = data_k;
}

/* Extract calibration parameters from geometry.rtk */
geom_file = "geometry.rtk";
cmd2 = "cat " + data_dir+geom_file +" | awk -F\"[><]\" '/SourceToIsocenterDistance/{printf \"%s\\n\", $3}; /SourceToDetectorDistance/{printf \"%s\\n\", $3}; /Angle|ProjectionOffsetX/{printf \"%s \", $3}; /ProjectionOffsetY/ {printf \"%s\\n\", $3}' > /tmp/calibration_param.txt";
system, cmd2;

calibraton_file = open("/tmp/calibration_param.txt");
Rsc = array(double,1);
Rsd = array(double,1);
n = read(calibraton_file, format="%f", Rsc);
n = read(calibraton_file, format="%f", Rsd);
param = array(double,3,ndata);
for (k=1; k<=ndata; ++k) {
    n = read(calibraton_file, format="%f %f %f", param(1,k), param(2,k), param(3,k));
}
close, calibraton_file;

DEG2RAD = pi/180.;
DATA = h_new(data = data,
             ndata = ndata,
             theta = DEG2RAD*param(1,),
             uoff = -1.0*param(3,),
             voff = -1.0*param(2,),
             Rsc = Rsc(1),
             Rsd = Rsd(1),
             nu = 512,
             nv = 512,
             u_scl = 409.6/512., // Dimensions détecteur : 40.96x40.96cm²
             v_scl = 409.6/512.
    );

yhd_save, data_dir+"Data4TRTT", DATA, overwrite=1n;

