
CC = cc
FC90 = gfortran
MPIF90 = 
FC = g77

## optimisation level
OPT = -O3
#OPT = -fast
#OPT = -cm -w95
#OPT = -cm -O3  -xW -w95 -ipo -tpp7


#LDFLAGS = -ipo  ## if ipo was used during compilation

#db = -g

# for lc at apac
#LAMHOME=/opt/lam-6.5.8/Intel
#INCLUDES = -I$(LAMHOME)/include

# for mpich at brecca.vpac
#MPIHOME = /usr/local/mpich-1.2.5-i-gm
#INCLUDES = -I$(MPIHOME)/include

## Compiler flags
FFLAGS = $(db) $(INCLUDES) $(OPT)\
# -warn unused -warn nousage -hpfwarn hpf \
# -warn declarations -warn argument_checking


FC90FLAGS = $(FFLAGS) 

# profiling while linking
#LDFLAGS = -p


MODS = globals.o interfaces.o

MODS = modules.o
LIBS = modules.o pcsi.o utils.o properties.o 

#gaussq.o dgamma.o d1mach.o
#SYSLIBS = -ldxml  -lmpi 
SYSLIBS =  -llapack -lm 

# lam in lc
#SYSLIBS =  -Vaxlib /opt/intel/mkl/lib/32/libmkl_lapack.a \
#          /opt/intel/mkl/lib/32/libmkl_p4.a \
#	  -L$(LAMHOME)/lib -llamf77mpi -lmpi -llam 

# for ifc in lc
#SYSLIBS =  -Vaxlib /usr/local/intel/mkl/lib/32/libmkl_p4.a \
#           /usr/local/intel/mkl/lib/32/libmkl_lapack.a


#  intel8 compiler in lc for using vml (mkl) library
#ifort -w -w95 -132 -D_LINUX -I../../include source/vserf.f  -L/opt/intel8-20050131/mkl72/lib/32 -lmkl_ia32 -lguide -lpthread -lm 

#SYSLIBS = -L/opt/intel8-20050131/mkl72/lib/32 \ -lmkl_solver -lmkl_lapack -lmkl_ia32 \ -lguide -lpthread -lm 


seqenscv-sp: dmodules.o gsipceq.o utils.o properties.o eq2ens-sp.o 
	$(FC90) -o $@ $(LDFLAGS) $+ $(SYSLIBS)
	strip $@


seqenscv: dmodules.f90 gsipceq.f90 utils.f90 properties.f90 eq2ens.f90 
	$(MPIF90) -o $@ $(FC90FLAGS) $(LDFLAGS) $+ $(SYSLIBS)
	strip $@


seqens-sp: dmodules.f90 gsipceq.f90 utils.f90 properties.f90 eqens-sp.f90 
	$(FC90) -o $@ $(FC90FLAGS) $(LDFLAGS) $+ $(SYSLIBS)
	strip $@



nseeded: dmodules.f90 gsipceq.f90 utils.f90 properties.f90 nseeded.f90 
	$(MPIF90) -o $@ $(FC90FLAGS) $(LDFLAGS) $+ $(SYSLIBS)
	strip $@

dltrajscv: dmodules.f90 gsipceq.f90 utils.f90 properties.f90 eq2ens.f90 
	$(MPIF90) -o $@ $(FC90FLAGS) $(LDFLAGS) $+ $(SYSLIBS)
	strip $@

ltrajs: modules.f90 gsipceq.f90 utils.f90 properties.f90 eqens.f90 
	$(MPIF90) -o $@ $(FC90FLAGS) $(LDFLAGS) $+ $(SYSLIBS)
	strip $@

a.out : fntest.f90
	$(FC) -o $@ $(FC90FLAGS) $(LDFLAGS) $+ $(SYSLIBS)


textra: textra.o hcoutils.o
	$(FC90) $(LDFLAGS) -o textra textra.o hcoutils.o 
	strip $@

fbspring: fb.f90 utils.f90
	$(FC90) -o $@ $(FC90FLAGS) $(LDFLAGS) $+ $(SYSLIBS)
	strip $@


.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .F .c

clean:
	rm *.o *.mod 

# Fixed source, if any:
.f.o:
	$(FC) -c $(FFLAGS) $<

.F.o:
	$(FC) -c $(CCOPTS) $<

.f.mod:
	$(FC) -c $(FC90FLAGS) $<

# Free form source:
.f90.o:
	$(FC90) -c $(FC90FLAGS) $<

.f90.mod:
	$(FC90) -c $(FC90FLAGS) $<

.c.o:
	$(CC) -c $(CFLAGS) $<
