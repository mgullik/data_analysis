subroutine print_lag_ener(rc_avefq_en, ic_avefq_en, err_std_rc, &
     err_std_ic, err_Aform_rc_ic, lag_en, err_Aform_lag, &
     err_cohe_lag, err_prop_lag)
  use dyn_lc 
  implicit none

  double precision, intent(IN) :: rc_avefq_en(freq_num, en_num), &
                         ic_avefq_en(freq_num, en_num), &
                         err_std_rc(freq_num, en_num), &
                         err_std_ic(freq_num, en_num), &
                         err_Aform_rc_ic(freq_num, en_num), &
                         lag_en (freq_num, en_num), &
                         err_Aform_lag(freq_num, en_num), &
                         err_cohe_lag(freq_num, en_num), &
                         err_prop_lag(freq_num, en_num)
  integer             :: jj, k
  character (len=200) :: filename, name_base, name_base_2
  
   filename = 'products/ener_rc.dat'
   open(84, file = filename)
   filename = 'products/ener_rc_err_Aform.dat'
   open(840, file = filename)
   filename = 'products/ener_rc_err_std.dat'
   open(841, file = filename)
   filename = 'products/ener_ic.dat'
   open(85, file = filename)
   filename = 'products/ener_ic_err_Aform.dat'
   open(850, file = filename)
   filename = 'products/ener_ic_err_std.dat'
   open(851, file = filename)
   filename = 'products/lag_ener.dat'
   open(86, file = filename)
   filename = 'products/lag_ener_err_Aform.dat'
   open(860, file = filename)
   filename = 'products/lag_ener_err_cohe.dat'
   open(861, file = filename)
   filename = 'products/lag_ener_err_prop.dat'
   open(862, file = filename)
  
!printing different stuff
   write(84,  *) 'skip on'
   write(840, *) 'skip on'
   write(841, *) 'skip on'
   write(85,  *) 'skip on'
   write(850, *) 'skip on'
   write(851, *) 'skip on'
   write(86,  *) 'skip on'
   write(860, *) 'skip on'
   write(861, *) 'skip on'
   write(862, *) 'skip on'
   write(84,  *) 'read serr  1, 2, 3'
   write(840, *) 'read serr  1, 2'
   write(841, *) 'read serr  1, 2'
   write(85,  *) 'read serr  1, 2, 3'
   write(850, *) 'read serr  1, 2'
   write(851, *) 'read serr  1, 2'
   write(86,  *) 'read serr  1, 2, 3, 4'
   write(860, *) 'read serr  1, 2, 3, 4'
   write(861, *) 'read serr  1, 2'
   write(862, *) 'read serr  1, 2'

   do jj = 1, freq_num 
      do k = 1, en_num 
         write(84, *)  (r_en_bin(k) + l_en_bin(k)) * 0.5, (r_en_bin(k) - l_en_bin(k)) * 0.5, rc_avefq_en(jj, k), err_Aform_rc_ic(jj, k), rc_avefq_en(jj, k), err_std_rc(jj, k)
         write(840, *)  (r_en_bin(k) + l_en_bin(k)) * 0.5, (r_en_bin(k) - l_en_bin(k)) * 0.5, rc_avefq_en(jj, k), err_Aform_rc_ic(jj, k)
         write(841, *)  (r_en_bin(k) + l_en_bin(k)) * 0.5, (r_en_bin(k) - l_en_bin(k)) * 0.5, rc_avefq_en(jj, k), err_std_rc(jj, k)
         
         write(85, *)  (r_en_bin(k) + l_en_bin(k)) * 0.5, (r_en_bin(k) - l_en_bin(k)) * 0.5, ic_avefq_en(jj, k), err_Aform_rc_ic(jj, k), ic_avefq_en(jj, k), err_std_ic(jj, k)
         write(850, *)  (r_en_bin(k) + l_en_bin(k)) * 0.5, (r_en_bin(k) - l_en_bin(k)) * 0.5, ic_avefq_en(jj, k), err_Aform_rc_ic(jj, k)
         write(851, *)  (r_en_bin(k) + l_en_bin(k)) * 0.5, (r_en_bin(k) - l_en_bin(k)) * 0.5, ic_avefq_en(jj, k), err_std_ic(jj, k)

         write(86, *)  (r_en_bin(k) + l_en_bin(k)) * 0.5, (r_en_bin(k) - l_en_bin(k)) * 0.5, lag_en(jj, k), err_Aform_lag(jj, k), lag_en(jj, k), err_cohe_lag(jj, k), lag_en(jj, k), err_prop_lag(jj, k)
         write(860, *)  (r_en_bin(k) + l_en_bin(k)) * 0.5, (r_en_bin(k) - l_en_bin(k)) * 0.5, lag_en(jj, k), err_Aform_lag(jj, k)
         write(861, *)  (r_en_bin(k) + l_en_bin(k)) * 0.5, (r_en_bin(k) - l_en_bin(k)) * 0.5, lag_en(jj, k), err_cohe_lag(jj, k)
         write(862, *)  (r_en_bin(k) + l_en_bin(k)) * 0.5, (r_en_bin(k) - l_en_bin(k)) * 0.5, lag_en(jj, k), err_prop_lag(jj, k) 

      enddo
      write(84,  *) 'no no no no'
      write(840, *) 'no no no no'
      write(841, *) 'no no no no'
      write(85,  *) 'no no no no'
      write(850, *) 'no no no no'
      write(851, *) 'no no no no'
      write(86,  *) 'no no no no'
      write(860, *) 'no no no no'
      write(861, *) 'no no no no'
      Write(862, *) 'no no no no'
   enddo


   write(84 , *) 'scr white'
   write(84 , *) 'log x on'
   write(84 , *) 'la y Real Cross'
   write(84 , *) 'la x Energy (KeV)'
   write(84 , *) 'lw 5'
   write(84 , *) 'tim off'
   write(840, *) 'scr white'
   write(840, *) 'log x on'
   write(840, *) 'la y Real Cross'
   write(840, *) 'la x Energy (KeV)'
   write(840, *) 'lw 5'
   write(840, *) 'tim off'
   write(841, *) 'scr white'
   write(841, *) 'log x on'
   write(841, *) 'la y Real Cross'
   write(841, *) 'la x Energy (KeV)'
   write(841, *) 'lw 5'
   write(841, *) 'tim off'
   write(85 , *) 'scr white'
   write(85 , *) 'log x on'
   write(85 , *) 'la y Imaginary Cross'
   write(85 , *) 'la x Energy (KeV)'
   write(85 , *) 'lw 5'
   write(85 , *) 'tim off'
   write(850, *) 'scr white'
   write(850, *) 'log x on'
   write(850, *) 'la y Imaginary Cross'
   write(850, *) 'la x Energy (KeV)'
   write(850, *) 'lw 5'
   write(850, *) 'tim off'
   write(851, *) 'scr white'
   write(851, *) 'log x on'
   write(851, *) 'la y Imaginary Cross'
   write(851, *) 'la x Energy (KeV)'
   write(851, *) 'lw 5'
   write(851, *) 'tim off'

   write(86 , *) 'scr white'
   write(86 , *) 'log x on'
   write(86 , *) 'la y Lag (s)'
   write(86 , *) 'la x Energy (KeV)'
   write(86 , *) 'lw 5'
   write(86 , *) 'tim off'
   write(860, *) 'scr white'
   write(860, *) 'log x on'
   write(860, *) 'la y Lag (s)'
   write(860, *) 'la x Energy (KeV)'
   write(860, *) 'lw 5'
   write(860, *) 'tim off'
   write(861, *) 'scr white'
   write(861, *) 'log x on'
   write(861, *) 'la y Lag (s)'
   write(861, *) 'la x Energy (KeV)'
   write(861, *) 'lw 5'
   write(861, *) 'tim off'

   close(84 ) 
   close(840) 
   close(841) 
   close(85 ) 
   close(850) 
   close(851) 
   close(86 ) 
   close(860) 
   close(861) 
   close(862) 

  write(*,*)"----------------------------------------------------------------"
  write(*,*)"Plots:"
  write(*,*)"ener_rc.dat ..... Re part of the cross spectrum (1:Adam error; 2:std error)"
  write(*,*)"ener_ic.dat ..... Im part of the cross spectrum (1:Adam error; 2:std error)"
  write(*,*)"lag_ener.dat ..... lag spectrum (1:Adam error; 2: coherence error; 3:prop error)"
  write(*,*)"----------------------------------------------------------------"


!XSPEC stuff  (use flx2xsp to create the pha file)
   do jj = 1, freq_num 
      write(name_base  , '(A,I1,A)') 'products/rc',jj - 1,'.dat'  
      write(name_base_2, '(A,I1,A)') 'products/ic',jj - 1,'.dat'  
      open(11, file = trim(name_base))
      open(12, file = trim(name_base_2))
       write(11, *) '0.1 0.3 0.01 0.01 '
      write(12, *) '0.1 0.3 0.01 0.01 '
      do k = 1, en_num 
         write(11, *)  l_en_bin(k), r_en_bin(k),  (r_en_bin(k) - l_en_bin(k)) * rc_avefq_en(jj, k), (r_en_bin(k) - l_en_bin(k)) * err_std_rc(jj, k)
         write(12, *)  l_en_bin(k), r_en_bin(k),  (r_en_bin(k) - l_en_bin(k)) * ic_avefq_en(jj, k), (r_en_bin(k) - l_en_bin(k)) * err_std_ic(jj, k)
      enddo
      write(11, *) '10.0 20. 0.01 0.01 '
      write(12, *) '10.0 20. 0.01 0.01'
      close(11)
      close(12)
   enddo
   write(*,*) 'Use flx2xsp to create the pha files from rc and ic'

! !XSPEC lags 
   do jj = 1, freq_num 
      write(name_base  , '(A,I1,A)') 'products/lag',jj - 1,'.dat'  
      open(11, file = trim(name_base))
      write(11, *) '0.1 0.3 0.01 0.01 '
      do k = 1, en_num 
         write(11, *)  l_en_bin(k), r_en_bin(k),  (r_en_bin(k) - l_en_bin(k)) * lag_en(jj, k), (r_en_bin(k) - l_en_bin(k)) *  err_cohe_lag(jj, k)
      enddo
      write(11, *) '10.0 20. 0.01 0.01 '
      close(11)
   enddo
   write(*,*) 'Use flx2xsp to create the pha files from the lag file'




end subroutine print_lag_ener
