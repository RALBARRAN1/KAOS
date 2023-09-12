module DensityProfileB

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	4- DENSITY PROFILE B:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! SET DENSITY DISTRIBUTION ONTO CONFIG. SPACE CELL-CENTERED .m files:

	subroutine DensityProfileBSub

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				do Qind= NqLB(1), NqUB(1), 1

					! ----------------------------------------------------

					allocate(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%pniICT &
						(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1)), &
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniICT &
						(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1)))

					! ----------------------------------------------------

					do FAindIC= 1, SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1), 1

						! Set initial dipole coordinates
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%pniICT(FAindIC)= &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%pGCT(1)
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniICT(FAindIC)= &
							SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1)

						! ----------------------------------------------------

						! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES WITHIN
						! CONFIGURATION-SPACE GRID BOUNDARIES:

						if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%pniICT(FAindIC))) &
							.eqv. .true.) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%pniICT(:)) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' pniICT HAS', &
								' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
								f, ', Qind= ', Qind, ', AND FAindIC= ', FAindIC, &
								' IN DENSITY PROFILE B SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniICT(FAindIC))) &
							.eqv. .true.) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniICT(:)) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qniICT HAS', &
								' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
								f, ', Qind= ', Qind, ', AND FAindIC= ', FAindIC, &
								' IN DENSITY PROFILE B SUBROUTINE' // achar(27) // '[0m.'
						end if

						if (SpecieT(s)%FluxTubeT(f)%QCellT(1)%qGLT(1) <= 0) then
							if (((Qind == NqLB(1)) .and. (SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniICT(FAindIC) &
								< SpecieT(s)%FluxTubeT(f)%QCellT(SpecieT(s)%FluxTubeT(f)%NqICAT(1)+ &
								Qind- 1)%qGLT(1))) .or. ((Qind == NqLB(1)) .and. &
								(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniICT(FAindIC) > &
								SpecieT(s)%FluxTubeT(f)%QCellT(SpecieT(s)%FluxTubeT(f)%NqICAT(1)+ &
								Qind- 1)%qGHT(1)))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qniICT= ', &
									SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniICT(FAindIC), &
									' VALUE OUT OF CONFIG-SPACE GRID CELL FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND FAindIC= ', &
									FAindIC, ' IN DENSITY PROFILE B SUBROUTINE' &
									// achar(27) // '[0m.'
							end if
							if (((Qind /= NqLB(1)) .and. (SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniICT(FAindIC) &
								<= SpecieT(s)%FluxTubeT(f)%QCellT(SpecieT(s)%FluxTubeT(f)%NqICAT(1)+ &
								Qind- 1)%qGLT(1))) .or. ((Qind /= NqLB(1)) .and. &
								(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniICT(FAindIC) > &
								SpecieT(s)%FluxTubeT(f)%QCellT(SpecieT(s)%FluxTubeT(f)%NqICAT(1)+ &
								Qind- 1)%qGHT(1)))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qniICT= ', &
									SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniICT(FAindIC), &
									' VALUE OUT OF CONFIG-SPACE GRID CELL FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND FAindIC= ', &
									FAindIC, ' IN DENSITY PROFILE B SUBROUTINE' &
									// achar(27) // '[0m.'
							end if
						end if

						if (SpecieT(s)%FluxTubeT(f)%QCellT(1)%qGLT(1) > 0) then
							if (((Qind == NqLB(1)) .and. (SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniICT(FAindIC) &
								> SpecieT(s)%FluxTubeT(f)%QCellT(SpecieT(s)%FluxTubeT(f)%NqICAT(1)+ &
								Qind- 1)%qGLT(1))) .or. ((Qind == NqLB(1)) .and. &
								(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniICT(FAindIC) < &
								SpecieT(s)%FluxTubeT(f)%QCellT(SpecieT(s)%FluxTubeT(f)%NqICAT(1)+ &
								Qind- 1)%qGHT(1)))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qniICT= ', &
									SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniICT(FAindIC), &
									' VALUE OUT OF CONFIG-SPACE GRID CELL FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND FAindIC= ', &
									FAindIC, ' IN DENSITY PROFILE B SUBROUTINE' &
									// achar(27) // '[0m.'
							end if
							if (((Qind /= NqLB(1)) .and. (SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniICT(FAindIC) &
								>= SpecieT(s)%FluxTubeT(f)%QCellT(SpecieT(s)%FluxTubeT(f)%NqICAT(1)+ &
								Qind- 1)%qGLT(1))) .or. ((Qind /= NqLB(1)) .and. &
								(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniICT(FAindIC) < &
								SpecieT(s)%FluxTubeT(f)%QCellT(SpecieT(s)%FluxTubeT(f)%NqICAT(1)+ &
								Qind- 1)%qGHT(1)))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qniICT= ', &
									SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniICT(FAindIC), &
									' VALUE OUT OF CONFIG-SPACE GRID CELL FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND FAindIC= ', &
									FAindIC, ' IN DENSITY PROFILE B SUBROUTINE' &
									// achar(27) // '[0m.'
							end if
						end if

						! ----------------------------------------------------

					end do

					! ----------------------------------------------------

				end do
			end do
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

		! PRINT TERMINAL EXPORT PARAMETERS:

		! do s= 1, Stot, 1
! 			do f= 1, SpecieT(s)%NfT(1), 1
! 				do Qind= NqLB(1), NqUB(1), 1
! 					if (rank == 0) then
! 						write(*, *) 'pniICT= ', s, f, Qind, SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%pniICT(:)
! 						write(*, *) 'qniICT= ', s, f, Qind, SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniICT(:)
! 						write(*, *) 'qniAT= ', s, f, Qind, SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniAT(1)
! 						write(*, *) 'qniBT= ', s, f, Qind, SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniBT(1)
! 					end if
! 				end do
! 			end do
! 		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				if (rank == 0) then
					call cpu_time(S4End)
					write(S4string, '(i10)')  nint(S4End)
					write(*, *) trim('%% 4- RANK= ' // adjustl(rankstring)) // &
						trim(', REAL-TIME= ' // adjustl(S4string)) // &
						trim(' s. DENSITY PROFILE B COMPLETE %%')
				end if
			end do
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine DensityProfileBSub

end module DensityProfileB
