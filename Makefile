#******************************** G77/Linux Fortran ************************
FC     =       gfortran
EXTRA_OPT =     
FFLAGS  =       -O $(EXTRA_OPT) 
LDFLAGS =
time_it         = get_cpu_f95

#******************************** G77/Linux Fortran ************************
#FC     =       g77
#EXTRA_OPT =     -mpentium -malign-double -fforce-mem -fforce-addr \
#                -ffast-math -funroll-all-loops
## May want to experiment by adding the extra optimization flags to get
## better runtime. But then again, maybe not.
#FFLAGS  =       -O2 $(EXTRA_OPT) -ffloat-store
#LDFLAGS = 
#time_it         = get_cpu_sun

#************************  Intel  Intel Fortran ************************
#FC     =       ifort
#EXTRA_OPT =
##if debugging, use this instead
##EXTRA_OPT = -check bounds
#FFLAGS  =       -O2 -r8 
#LDFLAGS = 
#time_it         = get_cpu_sun

#******************************** PGI Fortran ************************
#FC      =       pgf77
#FFLAGS  =      -fast
#LDFLAGS =	-fast 
#time_it         = get_cpu_sun

#******************************** Sun Fortran ************************
#FC     =       f77
#FFLAGS  =      -fast -O
#LDFLAGS =	-fast -O
#time_it         = get_cpu_sun

#******************************** Lahey-Fujitsu lf95 ************************
#
#FC      =       lf95
#FFLAGS  =       --tpp --nsav -O --nwarn -c 
##FFLAGS  =       --tpp --nsav -O --nwarn -c --chk
#LDFLAGS =
#time_it         = get_cpu_sun

#****************************************************************************


OBJSB	=	bset.o \
		check.o \
		densdiskg.o \
		densenv.o \
		diskflux.o \
		etime_gfortran.o \
		envset.o \
		errmsg.o \
		findangle4.o \
		findangle5.o \
		find_wall.o \
		find_wall_1D.o \
		find_wall_2D.o \
		findangle3.o \
		getset.o \
		gridset.o \
		initarr.o \
		initp.o \
                initpout.o \
		initpspot.o \
		linterp.o \
                locate.o \
		namer.o \
		newdisktherm.o \
		newtts.o \
		output.o \
		peeloff_3d_B.o \
		phiface.o \
		plancknu.o \
		propagate_B.o \
		radface.o \
		ran2.o \
		reapar.o \
		read_16scat.o \
		rotate.o \
		rotate_back.o \
		samptabl.o \
		setup.o \
		setup_wave.o \
		splines.o \
		spotset.o \
		stokes16_B.o \
		stokespeel16_B.o \
		subs_g77.o \
                Rdist.o \
		tauint_3d_B.o \
		testgrid.o \
		thetaface.o \
		trilinfmat.o \
		vger_3d_B.o \
		waveset.o \
		wrimsg.o \
		zerod.o \

mcvger3d_B:		$(OBJSB)
		$(FC) $(OBJSB) $(LDFLAGS) -o mcvger3d_B

tarfile:;	tar cvf mcvger3d_B.tar *.f *.txt *.in NOTES.* Makefile oldsrc

clean:;		/bin/rm -f *.o

