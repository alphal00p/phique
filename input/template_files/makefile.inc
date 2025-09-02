FC=
objdir=$(PWD)/obj/
moddir=$(PWD)/mod/
srcdir=$(PWD)/src/
libdir=$(PWD)/lib/
plotdir=$(PWD)/analysis/user/
gammaupcdir=$(PWD)/vendor/gammaUPC/
# uncomment the following line to use gdb to debug
# FFLAGS=-ggdb -g -fbounds-check -ffpe-trap=invalid,zero,overflow,underflow,denormal -Wall
objs=$(objdir)global_constants.o $(objdir)global.o $(objdir)qcd_constants.o $(objdir)qcd_setup.o $(objdir)qcd_coupling.o $(objdir)assertions.o $(objdir)runge_kutta.o $(objdir)newton_method.o $(objdir)coupling_global.o $(objdir)evaluate_couplings.o $(objdir)plot_NNLO.o $(objdir)MC_VEGAS.o $(objdir)photon_PDFs.o $(objdir)Constants.o $(objdir)Euler_gamma.o $(objdir)HarmonicSums.o $(objdir)potentialfunction.o $(objdir)phase_space_integrator.o $(objdir)driver.o
mod = $(moddir)/*.mod
exe = $(PWD)/bin/RunNNLO
linkexe = $(PWD)/RunNNLO
EXLIBS =
INCLUDE =
FFLAGS+=$(INCLUDE) -I$(gammaupcdir) -fPIC -w -fno-automatic -ffixed-line-length-132

all:    clean $(exe)
	ln -s $(exe) $(linkexe)

$(exe): makefile_hbook makefile_gammaUPC makefile_polygamma $(objs)
	@cp *.o $(objdir)
	@cp *.mod $(moddir)
	@rm *.o *.mod
	$(FC) -o $(exe) $(objs) -L$(libdir) -lhbook -lgammaUPC -lpolygamma $(EXLIBS)

$(objdir)phase_space_integrator.o $(moddir)phase_space_integrator.mod: $(srcdir)phase_space_integrator.f90 $(moddir)global.mod $(moddir)plot_NNLO.mod $(moddir)Constants.mod $(moddir)MC_VEGAS.mod $(moddir)photon_PDFs.mod $(moddir)global_constants.mod $(moddir)qcd_constants.mod $(moddir)newton_method.mod $(moddir)evaluate_couplings.mod $(moddir)potentialfunction.mod $(gammaupcdir)interpolation.mod $(gammaupcdir)nintlib.mod $(gammaupcdir)opticalglauber_geometry.mod $(gammaupcdir)elasticphotonphotonflux.mod
	$(FC) $(FFLAGS) -c $<

$(objdir)global.o $(moddir)global.mod: $(srcdir)global.f90
	$(FC) $(FFLAGS) -c $<

$(objdir)plot_NNLO.o $(moddir)plot_NNLO.mod: $(plotdir)plot_NNLO.f90 $(moddir)global.mod $(moddir)MC_VEGAS.mod
	$(FC) $(FFLAGS) -c $<

$(objdir)photon_PDFs.o $(moddir)photon_PDFs.mod: $(srcdir)photon_PDFs.f90 $(moddir)global.mod
	$(FC) $(FFLAGS) -c $<

$(objdir)coupling_global.o $(moddir)coupling_global.mod: $(srcdir)coupling_global.f90 $(moddir)qcd_coupling.mod
	$(FC) $(FFLAGS) -c $<

$(objdir)evaluate_couplings.o $(moddir)evaluate_couplings.mod: $(srcdir)evaluate_couplings.f90 $(moddir)global.mod $(moddir)global_constants.mod $(moddir)qcd_constants.mod $(moddir)coupling_global.mod
	$(FC) $(FFLAGS) -c $<

$(objdir)global_constants.o $(moddir)global_constants.mod: $(srcdir)global_constants.f90
	$(FC) $(FFLAGS) -c $<

$(objdir)Constants.o $(moddir)Constants.mod: $(srcdir)Constants.f90 $(moddir)global.mod
	$(FC) $(FFLAGS) -c $<

$(objdir)Euler_gamma.o $(moddir)Euler_gamma.mod: $(srcdir)Euler_gamma.f90
	$(FC) $(FFLAGS) -c $<

$(objdir)HarmonicSums.o $(moddir)HarmonicSums.mod: $(srcdir)HarmonicSums.f90 $(moddir)Euler_gamma.mod
	$(FC) $(FFLAGS) -c $<

$(objdir)potentialfunction.o $(moddir)potentialfunction.mod: $(srcdir)potentialfunction.f90 $(moddir)HarmonicSums.mod
	$(FC) $(FFLAGS) -c $<

$(objdir)qcd_constants.o $(moddir)qcd_constants.mod: $(srcdir)qcd_constants.f90
	$(FC) $(FFLAGS) -c $<

$(objdir)qcd_setup.o $(moddir)qcd_setup.mod: $(srcdir)qcd_setup.f90 $(moddir)qcd_constants.mod $(moddir)global_constants.mod
	$(FC) $(FFLAGS) -c $<

$(objdir)qcd_coupling.o $(moddir)qcd_coupling.mod: $(srcdir)qcd_coupling.f90 $(moddir)assertions.mod $(moddir)global_constants.mod $(moddir)qcd_constants.mod $(moddir)qcd_setup.mod $(moddir)runge_kutta.mod $(moddir)newton_method.mod
	$(FC) $(FFLAGS) -c $<

$(objdir)assertions.o $(moddir)assertions.mod: $(srcdir)assertions.f90 $(moddir)global_constants.mod
	$(FC) $(FFLAGS) -c $<

$(objdir)runge_kutta.o $(moddir)runge_kutta.mod: $(srcdir)runge_kutta.f90 $(moddir)global_constants.mod
	$(FC) $(FFLAGS) -c $<

$(objdir)newton_method.o $(moddir)newton_method.mod: $(srcdir)newton_method.f90 $(moddir)global_constants.mod
	$(FC) $(FFLAGS) -c $<

$(objdir)MC_VEGAS.o $(moddir)MC_VEGAS.mod: $(srcdir)MC_VEGAS.f90
	$(FC) $(FFLAGS) -c $<

$(objdir)driver.o: $(srcdir)driver.f90 $(moddir)phase_space_integrator.mod 
	$(FC) $(FFLAGS) -c $<

makefile_hbook :
	(cd ./analysis/hbook; make clean -f makefile_hbook ; make -f makefile_hbook; rm *.o; cd ../../)
makefile_gammaUPC :
	(cd ./vendor/gammaUPC; make clean -f makefile_gammaUPC ; make -f makefile_gammaUPC; rm *.o; cd ../../)
makefile_polygamma :
	(cd ./vendor/polygamma; make clean -f makefile_polygamma; make -f makefile_polygamma; rm *.o; cd ../../)

clean:
	@(cd ./analysis/hbook; make clean -f makefile_hbook; cd ../../)
	@(cd ./vendor/gammaUPC; make clean -f makefile_gammaUPC; cd ../../)
	@(cd ./vendor/polygamma; make clean -f makefile_polygamma; cd ../../)
	rm -rf $(objs) $(mod) $(exe) $(linkexe) $(PWD)/tmp/* *~ $(PWD)/input/*~ $(PWD)/src/*~ $(PWD)/output/* fort.* ./bin/*
end
