# main_freq = main_lag_freq.f90
main_freq = main_lag_freq_multiobs2.f90
main_ener = main_lag_ener_multiobs2.f90
main_ener_EK = main_lag_ener_multiobs2_EK.f90
# main_ener = main_lag_ener.f90
lc2pds    = main_lc_pds.f90
lc2pds_nu = main_lc_pds_nustar.f90

libs = -L ~/Software/cfitsio/ \
-L/Users/gullo/Software/heasoft-6.28/x86_64-apple-darwin18.7.0/lib\
-lcfitsio -lXSFunctions -lXSModel -lXSUtil -lXS -lCCfits_2.5

line_limit    = -ffree-line-length-none
# other_options = -fno-second-underscore -fno-automatic -fPIC

incs          = $(line_limit) $(other_options)

optimization  = -O3
debug         = -pedantic -Wall
comp          = gfortran 


# set compiler
comp          = gfortran $(optimization) $(debug) 

modulae = mod_dyn_lc.f90 mod_rebin.f90 mod_random.f90 

modul:
	$(comp) $(incs) $(libs) -c $(modulae)
objs_freq:
	$(comp) $(incs) $(libs) -c $(main_freq)
objs_ener:
	$(comp) $(incs) $(libs) -c $(main_ener)
objs_ener_EK:
	$(comp) $(incs) $(libs) -c $(main_ener_EK)
objs_lc:
	$(comp) $(incs) $(libs) -c $(lc2pds)
objs_lc_nu:
	$(comp) $(incs) $(libs) -c $(lc2pds_nu)

freq: clean modul objs_freq
	$(comp)  $(incs) $(libs) *.o -o freq

ener: clean modul objs_ener
	$(comp)  $(incs) $(libs) *.o -o ener

ener_EK: clean modul objs_ener_EK
	$(comp)  $(incs) $(libs) *.o -o ek_ener

lc2pds: clean modul objs_lc 
	$(comp)  $(incs) $(libs) *.o -o lc2pds

lc2pds_nu: clean modul objs_lc_nu
	$(comp)  $(incs) $(libs) *.o -o lc2pds_nu

clean:
	rm -f *.o *~ subroutines/*~ ener ek_ener freq lc2pds lc2pds_nu

clean_products:
	rm -f lag_freq* lag_ener* freq_cross_* ener_* fort.* products/* products/printed_lcs/*
