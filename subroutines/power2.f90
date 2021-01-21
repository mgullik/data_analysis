!-----------------------------------------------------------------------
      function power2(number)
! Returns TRUE if the number is a power of 2, if not it returns FALSE 
        implicit none
        integer :: number,x
        logical :: power2
        x = number
           ! write(*,*) x
        if (x .eq. 0)then 
           power2 = .false.
           return
        endif
        
        do while(MOD(x,2) .eq. 0)
           x = x / 2
        end do
        if (x .gt. 1) then
           power2 =.false.
           return
        else
           power2 = .true.
        endif
        
      end function power2
!-----------------------------------------------------------------------
