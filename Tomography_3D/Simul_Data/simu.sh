rm -fr *.orbit *.data

make clean
make make_CB_orbit
make make_image
make sim_cone

./make_CB_orbit boule.orbit
./make_image boule.obj boule.data
./sim_cone boule.obj boule.orbit proj_boule.data
