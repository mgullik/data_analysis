  function funUniformSingle() result(randUni)
    implicit none
    real randUni
    call random_seed
    call random_number(randUni) 
  end function funUniformSingle

!Poisson function -- returns a single Poisson random variable
  function funPoissonSingle(lambda) result(randPoisson)
    implicit none
    real, intent(in) :: lambda !input
    real :: exp_lambda !constant for terminating loop
    real :: randUni !uniform variable 
    real :: prodUni !product of uniform variables
    integer :: randPoisson !Poisson variable

    !declare functions
    real funUniformSingle
    
    exp_lambda= exp(-lambda) !calculate constant

    !initialize variables
    randPoisson = -1
    prodUni = 1
    do while (prodUni > exp_lambda)       
       randUni = funUniformSingle() !generate uniform variable
       prodUni = prodUni * randUni !update product
       randPoisson = randPoisson + 1 !increase Poisson variable
    end do
  end function funPoissonSingle
