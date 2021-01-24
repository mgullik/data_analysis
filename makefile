main_freq = main_lag_freq.f90
main_ener = main_lag_ener.f90
main_old  = re_im_freq.f90

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

freq: clean_products modul objs_freq
	$(comp)  $(incs) $(libs) *.o -o freq

ener: clean_products modul objs_ener
	$(comp)  $(incs) $(libs) *.o -o ener

old:
	$(comp) $(main_old) $(incs) $(libs) -o old_lag_freq

clean:
	rm -f *.o *~ subroutines/*~ 

clean_products:
	rm -f *.o *~ subroutines/*~ lag_freq_err_prop.dat lag_freq_err_cohe.dat cross_spec_imaginary_vs_freq.dat cross_spec_real_vs_freq.dat freq ener old_freq
