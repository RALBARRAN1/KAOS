module PoissonRNG

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	3.3- POISSON RANDOM NUMBER GENERATOR:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! COMPUTE POISSON RANDOM VARIATES FROM UNIFORM VARIATES:

! ----------------------------------------------------

	subroutine PoissonRNGSub

		! ----------------------------------------------------

		ExpPoisson(1)= exp(-MeanPoisson(1))
    RandPoisson(1)= -1d0
    UniformProd(1)= 1d0
    do while (UniformProd(1) > ExpPoisson(1))
      call random_number(RandUni)
      UniformProd(1)= UniformProd(1)*RandUni(1)
      RandPoisson(1)= RandPoisson(1)+ 1d0
    end do

		! ----------------------------------------------------

	end subroutine PoissonRNGSub

end module PoissonRNG
