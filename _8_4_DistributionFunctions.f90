module DistributionFunctions

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.4 DISTRIBUTION FUNCTIONS:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use IonDistributionFunctions
use ENADistributionFunctions

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! COMPUTE DISTRIBUTION FUNCTIONS [s^3/m^6] FOR ALL GRIDS:

	subroutine DistributionFunctionsSub

		! ----------------------------------------------------
		! ----------------------------------------------------

		call IonDistributionFunctionsSub

		!if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
		!	call ENADistributionFunctionsSub
		!end if

		! ----------------------------------------------------
		! ----------------------------------------------------

		if (rank == 0) then
			if (n == 1) then
				call cpu_time(S84End)
				write(S84string, '(F10.4)')  S84End
				write(*, *) trim('%% 8.2- REAL-TIME= ' // adjustl(S84string)) // &
					trim(' s. INITIAL DISTRIBUTION FUNCTIONS COMPLETE %%')
			end if
		end if

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine DistributionFunctionsSub

end module DistributionFunctions
