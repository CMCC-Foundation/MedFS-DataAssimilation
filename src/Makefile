SHELL = /bin/sh

##########################################################
#
#        Makefile for GLoceanvar
#
##########################################################

# ARCH = ifort_omp
# ARCH = sx9_omp_notrace
# ARCH = moon
# ARCH = athena_dbg
# ARCH = athena
ARCH = med
#ARCH = athena_mpi_dbg

include ./include.$(ARCH)

EXEC = ../bin/var_3d_$(ARCH).x
EXTA = ../bin/timeadd.x
EXTD = ../bin/timediff.x
EXTR = ../bin/prog_rho.x
LIBRARY = liboceanvar_$(ARCH).a
LIBLINK = -L. -loceanvar_$(ARCH)

ifeq ($(DEBUG),1) 
	F_FF = $(F_D)
else
	F_FF = $(F_O)
endif

OBJSTR  =  \
	set_knd.o\
	mod_obsdef.o\
	mod_adtest.o\
	drv_str.o\
	obs_str.o\
	grd_str.o\
	bmd_str.o\
	eof_str.o\
	ctl_str.o\
	mod_myfrtprof.o\
	mod_mynetcdf.o\
	mod_ceres.o\
        sdl_module.o\
        biascorr.o\
	sbcblk_core.o\
	tlad.o\
	tlad_test.o\
	tlad_vars.o\
	traadv_cen2h.o\
	traadv_cen2v.o\
	traldf_lap.o\
	trazdf_exp.o\
	trazdf_imp.o\
 	mod_iounits.o\
 	mod_run.o\
 	oceantools.o\
 	check_time.o\
 	gcdist.o\
        mod_mpirel.o\
        defprec.o\
        errioipsl.o\
        stringop.o\
        mathelp.o\
        calendar.o\
        mod_obshandling.o\
        mod_recfilter.o\
        read_driftsst.o\
        read_avhrr_sst_nc.o\
        readdrifterfile.o\
        read_amsr2_sst_nc.o\
	read_noaa_sst_nc.o\
        random_numbers.o\
        mod_lbclnk.o\
        mod_rdfield.o\
        mod_readfg.o\
	mod_writecorr.o\
	procsplit.o\
        filtrw.o\
        ecsort.o\
        mod_extgrid.o\
        mod_allocobs.o\
        mod_readobas.o\
        mod_readobas103.o\
        mod_debugtmp.o\
        mod_drwn.o\
        mod_ice.o\
        mod_bal.o\
        mod_bmd.o\
        mod_phinterp.o\
        mod_weakly.o\
        mod_hybrid.o\
        mod_bal_airsea.o\
	mod_eofs3d.o\
	mod_track_correl.o\
	mod_prof_correl.o\
	mod_sort.o\
        mod_stats.o
        
OBJS    =  \
	srfbal_ad.o\
	srfbal_tl.o\
	readair.o\
	readairfile.o\
	init0.o\
	ini_bal.o\
	bgpert.o\
	slamdtundef.o\
        rmvshallowins.o\
	prepsss.o\
	balance.o\
	balance_ad.o\
	dyn_height.o\
	dyn_height_ad.o\
	reg_rej.o\
	regrndrej.o\
	init_random_seed.o\
	dens_tl.o\
	dens_ad.o\
	matinv.o\
	matinv2.o\
        findfgtimes.o \
        findfgtimes_tl.o \
        readfg4df.o \
        readfgetab.o \
        readfgfields.o \
        readfgfields_a.o \
        read_firstguess.o \
        deallbg.o \
        rf_init.o \
        init_huber.o \
        coastrej.o \
        sla_bcorr.o \
        cg1.o \
        setfgfiles.o \
        obsstat.o \
        print_eofs.o \
        setfgflayout.o \
        setmppneighb.o \
    	io_obs.o \
    	writeobs102.o \
    	readobs103.o \
    	readobs103dd.o \
    	readobs103bin.o \
    	bgsde_obs.o \
	unbias_sla.o\
	unbias_sla_at.o\
        check_slabg.o\
	termin0.o\
	termin1.o\
	read_bineof.o\
	mdeofs.o\
	eigenv.o\
	readsst.o\
	rddns.o\
	rddns2.o\
	read_amsre_day.o\
	int_par_sst.o\
	obs_sst.o\
	obs_sst_tl.o\
	obs_sst_ad.o\
	sulogf.o\
        sumpi.o\
	abor1.o\
	abor1fl.o\
        posnam.o\
	routines.o\
	lbfgsb.o\
	cgfam.o\
	cgsearch.o\
	linpack.o\
	timer.o\
	def_nml.o\
	sav_itr.o\
	ini_itr.o\
	ini_bmd.o\
	rdgrds.o\
        handle_err.o\
	set_stdgrd.o\
	get_obs.o\
	get_obs_sla.o\
	get_obs_arg.o\
	get_obs_xbt.o\
	get_obs_gld.o\
	get_obs_vel.o\
	get_obs_trj.o\
	int_par_sla.o\
	int_par.o\
	obs_vec.o\
	sucov.o\
        reo4min.o\
	sueof.o\
	surf.o\
	surf_lv.o\
	surf_lvf.o\
	surfsr1z.o\
  	read_rfceg.o\
  	read_gtspp.o\
  	read_argo.o\
  	read_arglevs.o\
  	test_machp.o\
  	redcheck.o\
	ini_cfn.o\
	min_cfn.o\
	costf.o\
	cnv_ctv.o\
	ver_hor.o\
	ver_horf.o\
	ver_hor202.o\
	get_byg.o\
	veof.o\
	veoff.o\
	bar_mod.o\
	get_vel.o\
	div_dmp.o\
	bar_mod_ad.o\
	invrt.o\
	invrt_ad.o\
	shapiro_1d.o\
	obsop.o\
	obs_sla.o\
	obs_arg.o\
	obs_xbt.o\
	obs_gld.o\
	obs_vel.o\
	obs_trj.o\
	resid.o\
	res_inc.o\
	obsop_ad.o\
	obs_sla_ad.o\
	obs_arg_ad.o\
	obs_xbt_ad.o\
	obs_gld_ad.o\
	obs_vel_ad.o\
	obs_trj_ad.o\
	veof_ad.o\
	veoff_ad.o\
	ver_hor_ad.o\
	ver_horf_ad.o\
	ver_hor_ad202.o\
	get_vel_ad.o\
	prepprof.o\
	prep_extgrid.o\
	prep_extgridf.o\
	int_bacval_sla.o\
	div_dmp_ad.o\
	get_byg_ad.o\
	mod_trj_tl.o\
	mod_trj_ad.o\
	cnv_ctv_ad.o\
	cnv_inn.o\
	hrinterp.o\
	wrt_out.o\
	wrt_obs.o\
	wrt_obs_screen.o\
        varjob.o\
        setup0.o\
        m1qn3.o\
        readNetCDF_SLA_ATP.o\
        readnetcdf_sla_new.o\
        readnetcdf_sla_tap.o\
        read_remss_mwoisst_nc.o\
        read_l3s_sst_nc.o\
        read_l3_sst_nc.o\
        readsla.o\
        obs_sla2e.o\
        obs_sla2e_ad.o\
        obs_sla2.o\
        obs_sla2_ad.o\
        prep_lha.o\
        writecorr2d.o\
        wranal.o\
        wrmdt.o\
        wraninc.o\
        wraninchr.o\
        wraninctl.o\
        wranincad.o\
        check_err.o\
        gentrbk.o\
	readen3.o\
	readcrt.o\
	readcrtfile.o\
	readcora.o\
	readnrt.o\
	readen3file.o\
	readcorfile.o\
	readicoads.o\
	readicoadsfile.o\
	readicoads_nc.o\
	readicoadsfilenc.o\
	readwod.o\
	readwodfile.o\
	readinslevs.o\
	assign_slaerr.o\
	assign_ssterr.o\
        readobs.o\
        readobsa.o\
        readobsa103.o\
	int_par_ins.o\
	int_par_ins_2.o\
	int_par_sla_2.o\
	int_par_sst_2.o\
	int_par_sss_2.o\
	obs_ins.o\
	obs_ins_ad.o\
	obs_ins_tl.o\
	obs_vec2.o\
	int_par2.o\
	obserrors.o\
	assign_inserr.o\
	assign_icoerr.o\
	thinout_sla.o\
	thinout_sst.o\
	climcheck.o\
	read_clim.o\
	bgqcheck.o\
	vertcheck.o\
	bgqcheckp.o\
	read_bgerrs.o\
	usageov.o\
	calcdeta.o\
	thinout_ins_a.o\
	thinout_ins_b.o\
	thinning.o\
	repeve.o\
	blacklist.o\
        writencobs.o\
        writencobss.o\
        sfiltatslal.o\
        read_remss_swath_amsre_nc.o\
	assign_ssserr.o \
	int_par_sss.o \
	obs_sss.o \
	obs_sss_ad.o \
	obs_sss_tl.o \
	readsss.o \
	thinout_sss.o \
	read_aquarius.o \
	read_sacd.o \
	slaatcorr.o \
	getsstmb.o \
      	compalpha.o \
      	compalpha_lv.o \
      	compalpha_lvm.o \
      	compsc.o \
      	compsc_lv.o \
      	defextgrd.o \
      	defextgrd_lv.o \
      	defextgrd_lvm.o \
      	defextgrdsrz.o \
      	read_cr.o \
      	read_crz.o \
      	read_rfc.o \
      	read_rfcz.o \
      	def_coef.o \
      	rcfl_2.o \
      	rcfl_2_ad.o \
      	rcfl_x_2d.o \
      	rcfl_x_ad_2d.o \
      	rcfl_x_ad.o \
      	rcfl_x.o \
      	rcfl_xz_ad.o \
      	rcfl_xz.o \
      	rcfl_y_2d.o \
      	rcfl_y_ad_2d.o \
      	rcfl_y_ad.o \
      	rcfl_yz_ad.o \
      	rcfl_y_ad_old.o \
      	rcfl_y.o \
      	rcfl_yz.o \
      	rcfl_y_old.o \
      	rcflf_x.o \
      	rcflf_y.o \
      	rcflf_x_ad.o \
      	rcflf_y_ad.o \
      	rdfdf.o \
      	read_crf1.o \
      	rearralpha.o \
      	rearralpha_lv.o \
      	rearralpha_lvm.o \
      	rearralphasrz.o \
      	xind1.o \
      	check.o \
      	setuprfcr.o \
      	index_sort.o \
      	read_nceof.o \
      	write_nceof.o \
      	cofdis_2d.o \
      	haversine.o \
      	prep_c202.o \
      	read_slamis.o \
        master.o

AUXEC   =  \
	ecqsort.o\
	rsort32.o\
	rsort64.o

TIMEFL	= \
	mod_timefinc.o\
	adddate.o\
	addtime.o\
	addtimer.o\
	diffdate.o\
	difftime.o\
	difftimer.o\
	toupperdt.o

AUXVAR  = \
	CDFP_PublicReadAlongTrackProduct.o

.SUFFIXES: .o .F90 .F .c .F90r .F90s


all:  $(EXEC)
	@echo 
	@echo "****"
	@echo "All the executables have been compiled!"
	@echo "****"
	@echo 

logf :
	$(RM) $(LOGFILE)
	
$(LIBRARY) : $(OBJSTR) $(OBJS) $(AUXEC) $(AUXVAR) $(TIMEFL)

$(EXEC) : $(OBJSTR) $(OBJS) $(AUXEC) $(AUXVAR) $(TIMEFL)
	$(F_C) $(F_F) $(F_FF) $(P_P) master.F90 $(EXTINC) -o master.o
	$(F_L) $(L_O) -o $(EXEC) master.o $(LIBLINK) $(EXTLIB)

$(EXTA)	: $(TIMEFL)
	$(F_C) $(F_F) $(F_FF) $(P_P) timeadd.F90 -o timeadd.o
	$(F_L) $(F_F) $(F_FF) $(L_O) -o $(EXTA) $(TIMEFL) timeadd.o

$(EXTD)	: $(TIMEFL)
	$(F_C) $(F_F) $(F_FF) $(P_P) timediff.F90 -o timediff.o
	$(F_L) $(F_F) $(F_FF) $(L_O) -o $(EXTD) $(TIMEFL) timediff.o

$(EXTR)	: $(TIMEFL)
	$(F_C) $(F_F) $(F_FF) $(P_P) prog_rho.F90 -o prog_rho.o
	$(F_L) $(F_F) $(F_FF) $(L_O) -o $(EXTR) prog_rho.o \
	-L. -loceanvar_$(ARCH) $(EXTLIB)

.DEFAULTS:
%.o : %.mod
.F90.o :  $(OBJSTR)
	$(F_C) $(F_F) $(F_FF) $(P_P) $*.F90 $(EXTINC) -o $*.o
	$(AR_I) $(LIBRARY) $*.o	

.F90s.o :
	mv $*.F90s $*.F90
	$(F_CS) $(F_F) $(F_F77) -g $(NOOPT) $(P_P) $*.F90 $(EXTINC) -o $*.o
	$(AR_I) $(LIBRARY) $*.o	
	mv $*.F90 $*.F90s

.F90r.o :
	mv $*.F90r $*.F90
	$(F_C) $(F_F) $(F_F77) -g $(NOOPT) $(P_P) $*.F90 $(EXTINC) -o $*.o
	$(AR_I) $(LIBRARY) $*.o	
	mv $*.F90 $*.F90r

.F.o :           
	$(F_77) $(F_F) $(F_F77) $(F_FF) $(P_P) $*.F $(EXTINC) -o $*.o
	$(AR_I) $(LIBRARY) $*.o	

.c.o :           
	$(C_C) $(C_F) $*.c $(EXTINC) -o $*.o
	$(AR_I) $(LIBRARY) $*.o	

clean:
	$(RM) *.o *.mod *.lst ../build/*

depend:
	makedepf90 *.F90 -u netcdf > dependencies.inc

clean_all:	clean
	$(RM) $(EXEC) $(LIBRARY)

# include dependencies
include dependencies.inc
# DO NOT DELETE THIS LINE - used by make depend
