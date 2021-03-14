subroutine load_single_lc(filename)
  use dyn_lc
  implicit none 

  character (len=500)   :: filename

    call extract_lc(filename)
    
    if (allocated(lc)) then 
       call split_lc()
    else
       write(*,*) '    !!! NO LIGHT CURVE !!! '
       write(*,*) '        STOP HERE   '
       stop
    endif


    
  end subroutine load_single_lc
