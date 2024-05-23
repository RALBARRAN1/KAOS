module DensityProfileA

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	2- DENSITY PROFILE A:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use DensityProfileA1

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! SET FUNCTIONAL FORM OF DENSITY PROFILE ONTO CONFIG. SPACE CELL-CENTERED .m files:

	subroutine DensityProfileASub

		! ----------------------------------------------------

		call DensityProfileA1Sub

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR CONSISTENT TOTAL PARTICLE NUMBER:

		! ----------------------------------------------------

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				if (SpecieT(s)%FluxTubeT(f)%NsT(1) /= sum(SpecieT(s)%FluxTubeT(f)%NsFApT(:))) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT TOTAL PARTICLE NUMBER= ', SpecieT(s)%FluxTubeT(f)%NsT(1), &
						' AND FLUX-TUBE SUMMATION= ', sum(SpecieT(s)%FluxTubeT(f)%NsFARpT(:)), &
						' FOR SPECIE= ', s, ', AND FLUX TUBE= ', f, &
						' IN DENSITY PROFILE A SUBROUTINE' // achar(27) // '[0m.'
				end if

				if (rank == 0) then
					if (sum(SpecieT(s)%FluxTubeT(f)%NsFARpT(:)) /= SpecieT(s)%FluxTubeT(f)%NsRRT(1)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' INCONSISTENT CONFIG-SPACE PARTICLE NUMBER SUMMATION= ', sum(SpecieT(s)%FluxTubeT(f)%NsFARpT(:)), &
							' AND TOTAL PARTICLE NUMBER= ', SpecieT(s)%FluxTubeT(f)%NsRRT(1), &
							' FOR SPECIE= ', s, ', AND FLUX TUBE= ', f, &
							' IN DENSITY PROFILE A SUBROUTINE' // achar(27) // '[0m.'
					end if
				end if
			end do
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

		! PRINT TERMINAL EXPORT PARAMETERS:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				if (rank == 0) then
					write(*, *) 'Cluster Initial Ion Macro-particle Number per Config-Space Grid Cell= ', &
						SpecieT(s)%FluxTubeT(f)%NsFARpT(:)
					write(*, *)
					write(paramstring, '(i10)') sum(SpecieT(s)%FluxTubeT(f)%NsFARpT(:))
					write(*, *) trim('Cluster Total Initial Ion Macro-particle Number= ' // adjustl(paramstring))
					write(paramstring, '(i10)') nint(SpecieT(s)%FluxTubeT(f)%nsnormCLBGT(1))
					write(*, *) trim('Cluster LB Nominal Injection Macro-particle Number= ' // adjustl(paramstring))
					write(paramstring, '(i10)') nint(SpecieT(s)%FluxTubeT(f)%nsnormCUBGT(1))
					write(*, *) trim('Cluster UB Nominal Injection Macro-particle Number= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%NsT(1)
					write(*, *) trim('Root Rank Total Initial Ion Macro-particle Number= ' // adjustl(paramstring))

					write(paramstring, '(D10.4)') SpecieT(s)%FluxTubeT(f)%nsnormfacT(1)
					write(*, *) trim('Macro-particle Normalization Constant= ' // adjustl(paramstring))
				end if
			end do
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				if (rank == 0) then
					call cpu_time(S2End)
					write(S2string, '(F10.4)')  S2End
					write(*, *) trim('%% 2- REAL-TIME= ' // adjustl(S2string)) // &
						trim(' s. DENSITY PROFILE A COMPLETE %%')
				end if
			end do
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine DensityProfileASub

end module DensityProfileA
