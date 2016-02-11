CC=gcc -w -g
LIB = -static -L/usr/ucblib -lm -O4 
ALL : check_defs \
      sim_cone \
      sim_cone_centre \
      make_image \
      make_CB_orbit \
      make_ssc_orbit \
      make_circleline_orbit \
      make_line_orbit \
      make_c2l_orbit \
      make_sine_orbit \
      make_cirhel_orbit \
      make_cl_orbit \
      make_dual_orbit \
      make_helix_orbit \
      make_M_function \
      make_2d_filter \
      feldfft \
      cbbpj
check_defs : cone_defs.h cone_types.h check_defs.c
	$(CC) check_defs.c $(LIB) -o check_defs
	check_defs
sim_cone : cone_defs.h cone_types.h io.o simul.o sim_cone.c 
	$(CC) sim_cone.c io.o simul.o $(LIB) -o sim_cone
sim_cone_centre : cone_defs.h cone_types.h io.o simul.o sim_cone_centre.c 
	$(CC) sim_cone_centre.c io.o simul.o $(LIB) -o sim_cone_centre
make_image : cone_defs.h cone_types.h io.o simul.o make_image.c 
	$(CC) make_image.c io.o simul.o $(LIB) -o make_image
make_CB_orbit : cone_defs.h cone_types.h orbits.h io.o make_CB_orbit.c 
	$(CC) make_CB_orbit.c io.o $(LIB) -o make_CB_orbit
make_ssc_orbit : cone_defs.h cone_types.h orbits.h io.o make_ssc_orbit.c 
	$(CC) make_ssc_orbit.c io.o $(LIB) -o make_ssc_orbit
make_PB_orbit : cone_defs.h cone_types.h orbits.h io.o make_PB_orbit.c 
	$(CC) make_PB_orbit.c io.o $(LIB) -o make_PB_orbit
make_pbcb_orbit : cone_defs.h cone_types.h orbits.h io.o make_pbcb_orbit.c 
	$(CC) make_pbcb_orbit.c io.o $(LIB) -o make_pbcb_orbit
make_circleline_orbit : cone_defs.h cone_types.h orbits.h io.o make_circleline_orbit.c 
	$(CC) make_circleline_orbit.c io.o $(LIB) -o make_circleline_orbit
make_line_orbit : cone_defs.h cone_types.h orbits.h io.o make_line_orbit.c 
	$(CC) make_line_orbit.c io.o $(LIB) -o make_line_orbit
make_c2l_orbit : cone_defs.h cone_types.h orbits.h io.o make_c2l_orbit.c 
	$(CC) make_c2l_orbit.c io.o $(LIB) -o make_c2l_orbit
make_cirhel_orbit : cone_defs.h cone_types.h orbits.h io.o make_cirhel_orbit.c 
	$(CC) make_cirhel_orbit.c io.o $(LIB) -o make_cirhel_orbit
make_sine_orbit : cone_defs.h cone_types.h orbits.h scan_types.h io.o make_sine_orbit.c  
	$(CC) make_sine_orbit.c io.o $(LIB) -o make_sine_orbit
make_cl_orbit : cone_defs.h cone_types.h orbits.h io.o make_cl_orbit.c 
	$(CC) make_cl_orbit.c io.o $(LIB) -o make_cl_orbit
make_dual_orbit : cone_defs.h cone_types.h orbits.h io.o make_dual_orbit.c 
	$(CC) make_dual_orbit.c io.o $(LIB) -o make_dual_orbit
make_helix_orbit : cone_defs.h cone_types.h orbits.h io.o make_helix_orbit.c 
	$(CC) make_helix_orbit.c io.o $(LIB) -o make_helix_orbit
make_M_function : cone_defs.h cone_types.h io.o make_M_function.c 
	$(CC) make_M_function.c io.o $(LIB) -o make_M_function
make_2d_filter : cone_defs.h cone_types.h io.o fft.o make_2d_filter.c 
	$(CC) make_2d_filter.c io.o fft.o $(LIB) -o make_2d_filter
cbbpj : cone_defs.h cone_types.h lino.o io.o fft.o basic.o cone_basic.o chirpz.o cbbpj.c 
	$(CC) cbbpj.c lino.o io.o fft.o basic.o cone_basic.o chirpz.o $(LIB) -o cbbpj
feldfft : cone_defs.h cone_types.h io.o basic.o cone_basic.o fft.o feldfft.c
	$(CC) feldfft.c io.o basic.o cone_basic.o fft.o $(LIB) -o feldfft
lino.o :  cone_defs.h cone_types.h lino.c
	$(CC) -O4 -c lino.c
chirpz.o : chirpz.c
	$(CC) -O4 -c chirpz.c
io.o : cone_defs.h cone_types.h io.c
	$(CC) -c io.c
fft.o : fft.c
	$(CC) -O4 -c fft.c
simul.o : cone_defs.h cone_types.h simul.c
	$(CC) -O4 -c simul.c
basic.o : cone_defs.h cone_types.h basic.c
	$(CC) -O4 -c basic.c
cone_basic.o : cone_defs.h cone_types.h cone_basic.c
	$(CC) -O4 -c cone_basic.c
clean:
	rm -f core *~ *.o
