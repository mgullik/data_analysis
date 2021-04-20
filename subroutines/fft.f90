!------------------------------------------------------------------------
      subroutine ncperiodogram_no_norm(ht, st, rc, ic, dim)
! Calculates a cross spectrum between the time series ht(int_len_dim) and st(int_len_dim)
! Phase is such that +ve lag corresponds to ht lagging st
! ***MODIFIED*** from Press et al (1992) DFT definition, this is S^*(\nu)H(\nu)
        ! use dyn_lc
        implicit none
        integer, intent(IN)  :: dim
        double precision , intent(IN)  :: ht(dim), st(dim) 
        double precision , intent(OUT) :: rc(dim / 2), ic(dim / 2) 
        integer              :: j
        double precision , allocatable :: datah(:), datas(:)
        if (.not. allocated(datah)) allocate(datah(2 * dim))
        if (.not. allocated(datas)) allocate(datas(2 * dim))

        do j = 1, dim         
           datah(2 * j - 1) = ht(j)
           datah(2 * j)     = 0.d0
           datas(2 * j - 1) = st(j)
           datas(2 * j)     = 0.d0
        end do
        call four1(datah, dim, 1)
        call four1(datas, dim, 1)
        
        do j = 1, dim / 2
           rc(j) = datah(2 * j + 1) * datas(2 * j + 1) + datah(2 * j + 2) * datas(2 * j + 2)
           ic(j) = datah(2 * j + 2) * datas(2 * j + 1) - datah(2 * j + 1) * datas(2 * j + 2)

        end do
        if (allocated(datah)) deallocate(datah)
        if (allocated(datas)) deallocate(datas)
        return
      end subroutine ncperiodogram_no_norm
!------------------------------------------------------------------------

!------------------------------------------------------------------------
      subroutine periodogram_no_norm(time_series, pw, dim)
! Calculates power spectrum of the time series lc(dim)
! ***MODIFIED*** from Press et al (1992) DFT definition
        ! use dyn_lc
        implicit none
        integer, intent(IN)  :: dim
        double precision , intent(IN)  :: time_series(dim)
        double precision , intent(OUT) :: pw(dim / 2)

        integer              :: j
        double precision ,allocatable  :: datah(:)

        if (.not. allocated(datah)) allocate(datah(2 * dim))
        ! sum = 0.d0
        

        do j = 1, dim
           datah(2 * j - 1) = time_series(j)
           datah(2 * j)   = 0.d0
        end do

        call four1(datah, dim, 1)

        do j = 1, dim / 2
          pw(j) = (datah(2 * j + 1)**2 + datah(2 * j + 2 )**2)
        end do
        
        if (allocated(datah)) deallocate(datah)
        return
      end subroutine periodogram_no_norm
!------------------------------------------------------------------------

!------------------------------------------------------------------------
      subroutine ncperiodogram_frac_rms(ht, st, rc, ic, dt, dim)
! Calculates a cross spectrum between the time series ht(int_len_dim) and st(int_len_dim)
! In fractional  rms normalisation
! Phase is such that +ve lag corresponds to ht lagging st
! ***MODIFIED*** from Press et al (1992) DFT definition, this is S^*(\nu)H(\nu)
        ! use dyn_lc
        implicit none
        integer, intent(IN)  :: dim
        double precision , intent(IN)  :: dt, ht(dim), st(dim) 
        double precision , intent(OUT) :: rc(dim / 2), ic(dim / 2) 
        integer              :: j
        double precision               :: meanh, means
        double precision , allocatable :: datah(:), datas(:)
        if (.not. allocated(datah)) allocate(datah(2 * dim))
        if (.not. allocated(datas)) allocate(datas(2 * dim))

        do j = 1, dim         
           datah(2 * j - 1) = ht(j)
           datah(2 * j)     = 0.d0
           datas(2 * j - 1) = st(j)
           datas(2 * j)     = 0.d0
        end do
        call four1(datah, dim, 1)
        call four1(datas, dim, 1)

        meanh = datah(1) / dim
        means = datas(1) / dim
        
        do j = 1, dim / 2
           rc(j) = datah(2 * j + 1) * datas(2 * j + 1) + datah(2 * j + 2) * datas(2 * j + 2)
           rc(j) = rc(j) * 2.d0 * dt / (dble(dim) * meanh * means) 
           ic(j) = datah(2 * j + 2) * datas(2 * j + 1) - datah(2 * j + 1) * datas(2 * j + 2)
           ic(j) = ic(j) *  2.d0 * dt / (dble(dim) * meanh * means) 

        end do
        if (allocated(datah)) deallocate(datah)
        if (allocated(datas)) deallocate(datas)
        return
      end subroutine ncperiodogram_frac_rms
!------------------------------------------------------------------------

!------------------------------------------------------------------------
      subroutine periodogram_frac_rms(time_series, pw, dt, dim)
! Calculates power spectrum of the time series lc(dim)
! In fractional  rms normalisation
! ***MODIFIED*** from Press et al (1992) DFT definition
        ! use dyn_lc
        implicit none
        integer, intent(IN)  :: dim
        double precision , intent(IN)  :: dt, time_series(dim)
        double precision , intent(OUT) :: pw(dim / 2)

        integer              :: j
        double precision               :: mean !sum, var
        double precision ,allocatable  :: datah(:)

        if (.not. allocated(datah)) allocate(datah(2 * dim))
        ! sum = 0.0
        

        do j = 1, dim
           datah(2 * j - 1) = time_series(j)
           datah(2 * j)   = 0.d0
        end do

        call four1(datah, dim, 1)
        mean = datah(1) / dim

        do j = 1, dim / 2
          pw(j) = (datah(2 * j + 1)**2 + datah(2 * j + 2 )**2) * 2. * dt / (float(dim) * mean**2) 
        end do
        
        if (allocated(datah)) deallocate(datah)
        return
      end subroutine periodogram_frac_rms
!------------------------------------------------------------------------

!------------------------------------------------------------------------
      subroutine ncperiodogram(ht, st, rc, ic, dt, dim)
! Calculates a cross spectrum between the time series ht(int_len_dim) and st(int_len_dim)
! In absolute  rms normalisation
! Phase is such that +ve lag corresponds to ht lagging st
! ***MODIFIED*** from Press et al (1992) DFT definition, this is S^*(\nu)H(\nu)
        ! use dyn_lc
        implicit none
        integer, intent(IN)  :: dim
        double precision , intent(IN)  :: dt, ht(dim), st(dim) 
        double precision , intent(OUT) :: rc(dim / 2), ic(dim / 2) 
        integer              :: j
        double precision               :: meanh, means
        double precision , allocatable :: datah(:), datas(:)
        if (.not. allocated(datah)) allocate(datah(2 * dim))
        if (.not. allocated(datas)) allocate(datas(2 * dim))

        do j = 1, dim         
           datah(2 * j - 1) = ht(j)
           datah(2 * j)     = 0.d0
           datas(2 * j - 1) = st(j)
           datas(2 * j)     = 0.d0
        end do
        call four1(datah, dim, 1)
        call four1(datas, dim, 1)

        meanh = datah(1) / dim
        means = datas(1) / dim
        
        do j = 1, dim / 2
           rc(j) = datah(2 * j + 1) * datas(2 * j + 1) + datah(2 * j + 2) * datas(2 * j + 2)
           rc(j) = rc(j) * 2.d0 * dt / (dble(dim) ) 
           ic(j) = datah(2 * j + 2) * datas(2 * j + 1) - datah(2 * j + 1) * datas(2 * j + 2)
           ic(j) = ic(j) *  2.d0 * dt / (dble(dim)) 

        end do
        if (allocated(datah)) deallocate(datah)
        if (allocated(datas)) deallocate(datas)
        return
    end subroutine ncperiodogram
!------------------------------------------------------------------------

!------------------------------------------------------------------------
      subroutine periodogram(time_series, pw, dt, dim)
! Calculates power spectrum between the time series lc(dim)
! In absolute  rms normalisation
! Phase is such that +ve lag corresponds to ht lagging st
! ***MODIFIED*** from Press et al (1992) DFT definition
        ! use dyn_lc
        implicit none
        integer, intent(IN)  :: dim
        double precision , intent(IN)  :: dt, time_series(dim)
        double precision , intent(OUT) :: pw(dim / 2)

        integer              :: j
        ! double precision               :: mean, sum, var
        double precision ,allocatable  :: datah(:)

        if (.not. allocated(datah)) allocate(datah(2 * dim))
        ! sum = 0.0
        do j = 1, dim
           ! sum = sum + time_series(j)**2
           datah(2 * j - 1) = time_series(j)
           datah(2 * j)   = 0.d0
        end do

        call four1(datah, dim, 1)
        ! mean = datah(1) / dim
        ! write(*,*) 'mean ', mean
        
        do j = 1, dim / 2
           pw(j) = (datah(2 * j + 1)**2 + datah(2 * j + 2 )**2) * 2.d0 * dt / (float(dim)) 
        end do

!Check the Parseval theorem
        ! sum = 0.d0
        ! do j = 1, dim
        !    sum = sum + time_series(j)**2
        ! enddo        
        ! var = 0.d0
        ! do j = 0, dim - 1
        !    !           write(*,*) 2 * j + 1, 2 * j + 2           
        !    var = var + datah(2 * j + 1)**2 + datah(2 * j + 2 )**2
        !    ! write(111,*) j, datah(2 * j + 1)**2 + datah(2 * j + 2 )**2
        ! enddo
        ! var = var / dim

        ! write(*,*)  'Parseval theorem'
        ! write(*,*)  sum, var

        ! ! var = 0.0
        ! ! do j = 1, dim / 2 
        ! !    var = var + pw(j)
        ! ! enddo
        ! ! var = (2 * var + datah(1)**2)  / dim
        ! ! write(*,*)  'Parcival theorem'
        ! ! write(*,*)  sum, var
        
        if (allocated(datah)) deallocate(datah)
        return
      end subroutine periodogram
!------------------------------------------------------------------------


!------------------------------------------------------------------------
      subroutine periodogram_leahy(ht, pw, dt, dim)
! Calculates power spectrum between the time series ht(int_len_dim)
! In leahy normalisation 
! Phase is such that +ve lag corresponds to ht lagging st
! ***MODIFIED*** from Press et al (1992) DFT definition
        ! use dyn_lc
        implicit none
        integer          , intent(IN)  :: dim
        double precision , intent(IN)  :: ht(dim), dt
        double precision , intent(OUT) :: pw(dim / 2)

        integer              :: j
        ! double precision               :: sum!, var sum2, par
        double precision ,allocatable  :: datah(:)

        if (.not. allocated(datah)) allocate(datah(2 * dim))
        ! sum  = 0.0
        ! sum2 = 0.0
        do j = 1, dim
           ! sum  = sum  + ht(j)
           ! sum2 = sum2 + ht(j)**2
           datah(2 * j - 1) = ht(j)
           datah(2 * j)   = 0.d0
        end do
        call four1(datah, dim, 1)
        
        ! write(*,*) ' number of photons: ', sum, datah(1)

        do j = 1, dim / 2
           pw(j) = (datah(2 * j + 1)**2 + datah(2 * j + 2 )**2) * 2.d0 * dt / datah(1)
           ! write(11, *) j, pw(j)
        end do
        ! write(11, *) 'no no'
        
!Check the Parseval theorem
        ! par = 0.d0
        ! do j = 0, dim - 1
        !    !           write(*,*) 2 * j + 1, 2 * j + 2           
        !    par = par + datah(2 * j + 1)**2 + datah(2 * j + 2 )**2
        !    ! write(111,*) j, datah(2 * j + 1)**2 + datah(2 * j + 2 )**2
        ! enddo
        ! par = par / dble(dim)

        ! write(*,*)  'Parseval theorem'
        ! write(*,*)  sum2, par

! !CHECK THE VARIANCE CALCULATION        

! !First the variance from the light curve
!         var = (sum2 / dble(dim)) - (sum / dble(dim))**2 
!         write(*,*) 'variance from lc', var
        
! !Then variance from the fourier transform 
!         var = 0.d0
!         do j = 1, dim - 1
!            var = var + datah(2 * j + 1)**2 + datah(2 * j + 2 )**2
!         enddo
!         var = var / dble(dim)
!         write(*,*) 'variance from ft', var

!         ! do j = 0, dim - 1
!         !    write(10,*) j, datah(2 * j + 1)
!         ! enddo
!         ! write(10, *) 'no no'
        
        if (allocated(datah)) deallocate(datah)
        return
      end subroutine periodogram_leahy
!------------------------------------------------------------------------


!-----------------------------------------------------------------------
      SUBROUTINE four1(data,nn,isign)
      INTEGER isign,nn
      DOUBLE PRECISION data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      DOUBLE PRECISION tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      ! write(*,*) 'ciao ciao '

      ! do i = 1, nn
      !    write(*,*) i, data(2 * i - 1), data(2 * i)
      ! enddo
      ! write(*,*) 'ciao ciao '
        
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      end
!---------------------------------------------------------------

!----------------------------------------------------------------
      subroutine ncperiodogram_no_norm_real(ht, st, rc, ic, dim)
! Calculates a cross spectrum between the time series ht(int_len_dim) and st(int_len_dim)
! Phase is such that +ve lag corresponds to ht lagging st
! ***MODIFIED*** from Press et al (1992) DFT definition, this is S^*(\nu)H(\nu)
        use dyn_lc
        implicit none
        integer, intent(IN)  :: dim
        real   , intent(IN)  :: ht(dim), st(dim) 
        real   , intent(OUT) :: rc(dim / 2), ic(dim / 2) 
        integer              :: j
        real   , allocatable :: datah(:), datas(:)
        if (.not. allocated(datah)) allocate(datah(2 * dim))
        if (.not. allocated(datas)) allocate(datas(2 * dim))

        do j = 1, dim         
           datah(2 * j - 1) = ht(j)
           datah(2 * j)     = 0.0
           datas(2 * j - 1) = st(j)
           datas(2 * j)     = 0.0
        end do
        call four1_real(datah, dim, 1)
        call four1_real(datas, dim, 1)
        
        do j = 1, dim / 2
           rc(j) = datah(2 * j + 1) * datas(2 * j + 1) + datah(2 * j + 2) * datas(2 * j + 2)
           ic(j) = datah(2 * j + 2) * datas(2 * j + 1) - datah(2 * j + 1) * datas(2 * j + 2)

        end do
        if (allocated(datah)) deallocate(datah)
        if (allocated(datas)) deallocate(datas)
        return
      end subroutine ncperiodogram_no_norm_real
!-----------------------------------------------------------------

!-----------------------------------------------------------------
      subroutine periodogram_no_norm_real(time_series, pw, dim)
! Calculates power spectrum of the time series lc(dim)
! ***MODIFIED*** from Press et al (1992) DFT definition
        use dyn_lc
        implicit none
        integer, intent(IN)  :: dim
        real   , intent(IN)  :: time_series(dim)
        real   , intent(OUT) :: pw(dim / 2)

        integer              :: j
        real   ,allocatable  :: datah(:)

        if (.not. allocated(datah)) allocate(datah(2 * dim))
        ! sum = 0.0
        

        do j = 1, dim
           datah(2 * j - 1) = time_series(j)
           datah(2 * j)   = 0.0
        end do

        call four1_real(datah, dim, 1)

        do j = 1, dim / 2
          pw(j) = (datah(2 * j + 1)**2 + datah(2 * j + 2 )**2)
        end do
        
        if (allocated(datah)) deallocate(datah)
        return
      end subroutine periodogram_no_norm_real
!-------------------------------------------------------------------

!-------------------------------------------------------------------
      subroutine ncperiodogram_frac_rms_real(ht, st, rc, ic, dim)
! Calculates a cross spectrum between the time series ht(int_len_dim) and st(int_len_dim)
! In fractional  rms normalisation
! Phase is such that +ve lag corresponds to ht lagging st
! ***MODIFIED*** from Press et al (1992) DFT definition, this is S^*(\nu)H(\nu)
        use dyn_lc
        implicit none
        integer, intent(IN)  :: dim
        real   , intent(IN)  :: ht(dim), st(dim) 
        real   , intent(OUT) :: rc(dim / 2), ic(dim / 2) 
        integer              :: j
        real                 :: meanh, means
        real   , allocatable :: datah(:), datas(:)
        if (.not. allocated(datah)) allocate(datah(2 * dim))
        if (.not. allocated(datas)) allocate(datas(2 * dim))

        do j = 1, dim         
           datah(2 * j - 1) = ht(j)
           datah(2 * j)     = 0.0
           datas(2 * j - 1) = st(j)
           datas(2 * j)     = 0.0
        end do
        call four1_real(datah, dim, 1)
        call four1_real(datas, dim, 1)

        meanh = datah(1) / dim
        means = datas(1) / dim
        
        do j = 1, dim / 2
           rc(j) = datah(2 * j + 1) * datas(2 * j + 1) + datah(2 * j + 2) * datas(2 * j + 2)
           rc(j) = rc(j) * 2. * real(dt) / (real(dim) * meanh * means) 
           ic(j) = datah(2 * j + 2) * datas(2 * j + 1) - datah(2 * j + 1) * datas(2 * j + 2)
           ic(j) = ic(j) *  2. * real(dt) / (real(dim) * meanh * means) 

        end do
        if (allocated(datah)) deallocate(datah)
        if (allocated(datas)) deallocate(datas)
        return
      end subroutine ncperiodogram_frac_rms_real
!-------------------------------------------------------------------

!-------------------------------------------------------------------
      subroutine periodogram_frac_rms_real(time_series, pw, dim)
! Calculates power spectrum of the time series lc(dim)
! In fractional  rms normalisation
! ***MODIFIED*** from Press et al (1992) DFT definition
        use dyn_lc
        implicit none
        integer, intent(IN)  :: dim
        real   , intent(IN)  :: time_series(dim)
        real   , intent(OUT) :: pw(dim / 2)

        integer              :: j
        real                 :: mean !sum, var
        real   ,allocatable  :: datah(:)

        if (.not. allocated(datah)) allocate(datah(2 * dim))
        ! sum = 0.0
        

        do j = 1, dim
           datah(2 * j - 1) = time_series(j)
           datah(2 * j)   = 0.0
        end do

        call four1_real(datah, dim, 1)
        mean = datah(1) / dim

        do j = 1, dim / 2
          pw(j) = (datah(2 * j + 1)**2 + datah(2 * j + 2 )**2) * 2. * real(dt) / (float(dim) * mean**2) 
        end do
        
        if (allocated(datah)) deallocate(datah)
        return
      end subroutine periodogram_frac_rms_real
!------------------------------------------------------------------

!------------------------------------------------------------------
      subroutine ncperiodogram_real(ht, st, rc, ic, dim)
! Calculates a cross spectrum between the time series ht(int_len_dim) and st(int_len_dim)
! In absolute  rms normalisation
! Phase is such that +ve lag corresponds to ht lagging st
! ***MODIFIED*** from Press et al (1992) DFT definition, this is S^*(\nu)H(\nu)
        use dyn_lc
        implicit none
        integer, intent(IN)  :: dim
        real   , intent(IN)  :: ht(dim), st(dim) 
        real   , intent(OUT) :: rc(dim / 2), ic(dim / 2) 
        integer              :: j
        real                 :: meanh, means
        real   , allocatable :: datah(:), datas(:)
        if (.not. allocated(datah)) allocate(datah(2 * dim))
        if (.not. allocated(datas)) allocate(datas(2 * dim))

        do j = 1, dim         
           datah(2 * j - 1) = ht(j)
           datah(2 * j)     = 0.0
           datas(2 * j - 1) = st(j)
           datas(2 * j)     = 0.0
        end do
        call four1_real(datah, dim, 1)
        call four1_real(datas, dim, 1)

        meanh = datah(1) / dim
        means = datas(1) / dim
        
        do j = 1, dim / 2
           rc(j) = datah(2 * j + 1) * datas(2 * j + 1) + datah(2 * j + 2) * datas(2 * j + 2)
           rc(j) = rc(j) * 2.0 * real(dt) / (real(dim) ) 
           ic(j) = datah(2 * j + 2) * datas(2 * j + 1) - datah(2 * j + 1) * datas(2 * j + 2)
           ic(j) = ic(j) *  2.0 * real(dt) / (real(dim)) 

        end do
        if (allocated(datah)) deallocate(datah)
        if (allocated(datas)) deallocate(datas)
        return
      end subroutine ncperiodogram_real
!-------------------------------------------------------------------

!-------------------------------------------------------------------
      subroutine periodogram_real(time_series, pw, dim)
! Calculates power spectrum between the time series lc(dim)
! In absolute  rms normalisation
! Phase is such that +ve lag corresponds to ht lagging st
! ***MODIFIED*** from Press et al (1992) DFT definition
        use dyn_lc
        implicit none
        integer, intent(IN)  :: dim
        real   , intent(IN)  :: time_series(dim)
        real   , intent(OUT) :: pw(dim / 2)

        integer              :: j
        real                 :: mean !sum, var
        real   ,allocatable  :: datah(:)

        if (.not. allocated(datah)) allocate(datah(2 * dim))
        ! sum = 0.0
        

        do j = 1, dim
           ! sum = sum + ht(j)**2
           datah(2 * j - 1) = time_series(j)
           datah(2 * j)   = 0.0
        end do

        call four1_real(datah, dim, 1)
        mean = datah(1) / dim

        do j = 1, dim / 2
           pw(j) = (datah(2 * j + 1)**2 + datah(2 * j + 2 )**2) * 2. * real(dt) / (float(dim)) 
        end do

!Check the Parcival theorem
!         var = 0.0
!         do j = 0, dim - 1
! !           write(*,*) 2 * j + 1, 2 * j + 2
!            var = var + datah(2 * j + 1)**2 + datah(2 * j + 2 )**2
!            write(111,*) j, datah(2 * j + 1)**2 + datah(2 * j + 2 )**2
!         enddo
!         var = var / dim
!         write(*,*)  'Parcival theorem'
!         write(*,*)  sum, var

!         var = 0.0
!         do j = 1, dim / 2 
!            var = var + pw(j)
!         enddo
!         var = (2 * var + datah(1)**2)  / dim
!         write(*,*)  'Parcival theorem'
!         write(*,*)  sum, var
        
        if (allocated(datah)) deallocate(datah)
        return
      end subroutine periodogram_real
!-------------------------------------------------------------------


!-------------------------------------------------------------------
      subroutine periodogram_leahy_real(ht, pw)
! Calculates power spectrum between the time series ht(int_len_dim)
! In leahy normalisation 
! Phase is such that +ve lag corresponds to ht lagging st
! ***MODIFIED*** from Press et al (1992) DFT definition
        use dyn_lc
        implicit none
        real   , intent(IN)  :: ht(int_len_dim)
        real   , intent(OUT) :: pw(int_len_dim / 2)

        integer              :: j
        real                 :: mean
        real   ,allocatable  :: datah(:)

        if (.not. allocated(datah)) allocate(datah(2 * int_len_dim))
        ! sum = 0.0
        do j = 1, int_len_dim
           ! sum = sum + ht(j)**2
           datah(2 * j - 1) = ht(j)
           datah(2 * j)   = 0.0
        end do
        call four1_real(datah, int_len_dim, 1)
        mean = datah(1) / int_len_dim
        do j = 1, int_len_dim / 2
           pw(j) = (datah(2 * j + 1)**2 + datah(2 * j + 2 )**2) * 2.0 / datah(1)  
        end do
        
        if (allocated(datah)) deallocate(datah)
        return
      end subroutine periodogram_leahy_real
!-----------------------------------------------------------------

!-----------------------------------------------------------------
      SUBROUTINE four1_real(data,nn,isign)
      INTEGER isign,nn
      REAL data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      ! write(*,*) 'ciao ciao '

      ! do i = 1, nn
      !    write(*,*) i, data(2 * i - 1), data(2 * i)
      ! enddo
      ! write(*,*) 'ciao ciao '
        
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      end
!-----------------------------------------------------------------------




!---------------------------------------------------------------------!
subroutine cross_FT(lc1, lc2, rc, ic, NN, dt)
! This subroutine calculates the cross-spectrum between lc1 and lc2
! lc2 is the complex cojugate light curve (reference light curve)

! re1, im1, re2, im2 are the real and imaginary part to calculare the 
! Fourier transform 
  implicit none 
  integer         , intent(IN)    :: NN
  double precision, intent(IN)    :: lc1(NN), lc2(NN), dt
  double precision, intent(OUT)   :: rc(NN / 2), ic(NN / 2)
  
  integer                         :: i
  double precision, allocatable   :: re1(:), im1(:), re2(:), im2(:)

  allocate (re1(NN / 2))
  allocate (im1(NN / 2))
  allocate (re2(NN / 2))
  allocate (im2(NN / 2))
 
  call FT_not_fast(lc1, re1, im1, NN)
  call FT_not_fast(lc2, re2, im2, NN)

  do i = 1, NN / 2    
     rc(i) = (re1(i) * re2(i) + im1(i) * im2(i)) *  2 * dt / (real(NN))
     ic(i) = (im1(i) * re2(i) - re1(i) * im2(i)) *  2 * dt / (real(NN))
  enddo

   end subroutine Cross_FT
!---------------------------------------------------------------------!


!---------------------------------------------------------------------!
subroutine power_FT(lc1, pw, NN, dt)
! This subroutine calculates the power-spectrum of lc1
! re1, im1 are the real and imaginary part to calculare the Fourier transform 
  implicit none 
  integer         , intent(IN)    :: NN
  double precision, intent(IN)    :: lc1(NN), dt 
  double precision, intent(OUT)   :: pw(NN / 2)
  
  integer                         :: i
  double precision, allocatable   :: re1(:), im1(:)

  allocate (re1(NN / 2))
  allocate (im1(NN / 2))
 
  call FT_not_fast(lc1, re1, im1, NN)

  do i = 1, NN / 2    
     pw(i) = (re1(i) * re1(i) + im1(i) * im1(i) ) *  2 * dt / (real(NN))
  enddo

end subroutine Power_FT
!---------------------------------------------------------------------!


!---------------------------------------------------------------------!
subroutine FT_not_fast(time_series, re, im, NN)
! This subroutine calculates the Fourier transfor of a series time_series with NN data points 
! re and im are respectively the real and the imaginary part of the FT
! re and im don't store the average of the light curve, so they have NN/2
! The first frequency is nu = 1/(NN*dt) the last is the Nyquist frequency nuNy = 1/(2*dt)
  implicit none 
  integer         , intent(IN)    :: NN
  double precision, intent(IN)    :: time_series(NN)
  double precision, intent(INOUT) :: re(NN / 2), im(NN / 2)
  
  integer                         :: n, k
  double precision                :: arg
  double precision, parameter     :: pi = acos(-1.0)
  

  do n = 1, NN / 2  
     re(n) = 0.0 
     im(n) = 0.0 

     do k = 1, NN 
        arg = 2.0 * pi * real(k * n) / real(NN)
        re(n) = re(n) + time_series(k) * cos(arg)
        im(n) = im(n) + time_series(k) * sin(arg)
     enddo
     re(n) = re(n) / real(NN)
     im(n) = im(n) / real(NN)
  enddo

end subroutine FT_not_fast
!---------------------------------------------------------------------!


!-----------------------------------------------------------------------
subroutine SFT(at,nnmax,nn,Af)
  implicit none
  integer nnmax,nn,j,k
  real at(nnmax)
  complex Af(nnmax/2)
  real pi,arg
  pi = acos(-1.0)
  do j = 1,nn/2
     Af(j) = 0.0
     do k = 1,nn
        arg   = 2.0 * pi * real( j * k ) / real(nn)
        Af(j) = Af(j) + at(k) * complex( cos(arg) , sin(arg) )
     end do
     Af(j) = Af(j) / real(nn)
  end do
  return
end subroutine SFT
!-----------------------------------------------------------------------
