data_dir = "/home/momey/Recherche_Tomographie/Data/data_CLB_09-11-2011/img_1.3.46.423632.135428.1320854585.13/";
proj_file = open(data_dir+"proj_list.txt");
projections=rdfile(proj_file);
close, file;

ndata = numberof(projections);
data = array(double,512,512,ndata);
for (k=1; k<=ndata; ++k) {
    write, format="projection n°%i\n", k;
    data(..,k) = double(raw_read(data_dir+projections(k),short,512,512,offset=100) + (2.^15));
}

for (k=1; k<=ndata; ++k) {
    mywindow, 0; fma; img_plot, data(..,k), cbar=1, vert=1, format="%.3f";
    pause, 10;
}
