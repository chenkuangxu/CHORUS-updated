FCOMP    = mpif90
OPTS     = -c -O3 
#OPTS     = -c -O0 -g -traceback -fpe:0 -check all -fpstkchk

SWP	= swplist
LINKOPTS = -O3 -o
OBJS = setup3d.o main.o readinput.o readGRIDinfo.o \
	connectivity.o connect_bdry.o map_proc_match.o \
	getIVCELL_proc.o map_interior_interface.o \
	getedges.o init_setup.o map_der_basefunc.o \
	calcjacob.o setinitialcond.o tecplotter3d.o \
	tecplotter3dsetup.o tecplot1.o getrusanovflux.o \
	getfluxvectors.o compflux.o interfaceflux_all.o \
	procintflux_all.o bcflux_e.o bcvisflux_e.o nablaQsp.o \
	getVfluxvectors.o getVfluxvectors_bot.o getVfluxvectors_top.o \
	compresid.o iterations.o write_data.o read_data_restart.o \
	transfinite_mapping.o calcjacob_transfinite.o \
	gauss_integration.o
	

sd3dhexa:$(OBJS)
	$(FCOMP) $(LINKOPTS) ./rundir/sd3d_sun $(OBJS)

clean:
	 rm  *.o *.mod
#	$(RM) $(EXEC) $(OBJS)

.SUFFIXES:
.SUFFIXES: .o .F .c .f .swp .f90

.F.o:
	$(FCOMP) $(OPTS) $(DEFINES) $<

.f.o:
	$(FCOMP)  $(OPTS) $(DEFINES) $<

.f90.o:
	$(FCOMP)  $(OPTS) $(DEFINES) $<

.F.swp:
	$(SWP) -c $(FFLAGS) $(QV_OPT) $(DEFINES) -WK,-cmp=$*.m $<

.f.swp:
	$(SWP) -c $(FFLAGS) $(QV_OPT) $(DEFINES) -WK,-cmp=$*.m $<

.c.swp:
	$(SWP) -c $(CFLAGS) $(QV_OPT) $(DEFINES) $<
