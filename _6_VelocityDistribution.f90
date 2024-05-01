module VelocityDistribution

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	6- VELOCITY DISTRIBUTION:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use MBVelocityDistribution

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! ALLOCATE ALL ION PARTICLE INITIAL VELOCITIES BY CARTESIAN MB DISTRIBUTION:
! Note: Make initial velocity distribution flux tube and config-space grid cell dependent
! to accomodate a position dependent initial particle species temperature.

	subroutine VelocityDistributionSub

		! ----------------------------------------------------

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				ICflag(1)= 1d0

				do Qind= NqLB(1), NqUB(1), 1

					! ----------------------------------------------------

					allocate(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%Vperp1ICT &
						(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1)), &
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%Vperp2ICT &
						(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1)), &
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VperpICT &
						(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1)), &
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VparICT &
						(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1)), &
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VxICT &
						(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1)), &
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VyICT &
						(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1)), &
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VzICT &
						(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1)))

					do FAindIC= 1, SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1), 1

						! ----------------------------------------------------

						thetaMB(1)= SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)% &
							thetafinalICT(FAindIC)
						phiMB(1)= SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)% &
							phifinalICT(FAindIC)
						ellMB(1)= SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)% &
							ellfinalICT(FAindIC)

						if (SpecieT(1)%FluxTubeT(1)%SPINUPflagT(1) == 1) then
							TsMB(1)= SpecieT(s)%FluxTubeT(f)%TsT(1, SpecieT(s)%FluxTubeT(f)%NqICAT(1)+ Qind- 1)
						end if
						if (SpecieT(1)%FluxTubeT(1)%SPINUPflagT(1) == 0) then
							TsMB(1)= SpecieT(s)%FluxTubeT(f)%QCellT(SpecieT(s)% &
				  			FluxTubeT(f)%NqICAT(1)+ Qind- 1)%TemperatureInputT(1)
						end if

						! ----------------------------------------------------

						call MBVelocityDistributionSub

						! ----------------------------------------------------

						! Create nested derived data types
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%Vperp1ICT(FAindIC)= Vperp1MB(1)
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%Vperp2ICT(FAindIC)= Vperp2MB(1)
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VperpICT(FAindIC)= VperpMB(1)
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VparICT(FAindIC)= VparMB(1)
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VxICT(FAindIC)= VxMB(1)
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VyICT(FAindIC)= VyMB(1)
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VzICT(FAindIC)= VzMB(1)

						! ----------------------------------------------------

            ! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSION, SIZES AND
						! FINITE VALUES:

						if ((Vperp1MB(1) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%Vperp1ICT(FAindIC)) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%Vperp1ICT(FAindIC))) &
							.eqv. .true.) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%Vperp1ICT(:)) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp1MBT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND FAindIC= ', FAindIC, &
								' IN VELOCITY DISTRIBUTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((Vperp2MB(1) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%Vperp2ICT(FAindIC)) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%Vperp2ICT(FAindIC))) &
							.eqv. .true.) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%Vperp2ICT(:)) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp2MBT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND FAindIC= ', FAindIC, &
								' IN VELOCITY DISTRIBUTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((VperpMB(1) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VperpICT(FAindIC)) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VperpICT(FAindIC))) &
							.eqv. .true.) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VperpICT(:)) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VperpMBT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND FAindIC= ', FAindIC, &
								' IN VELOCITY DISTRIBUTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((VparMB(1) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VparICT(FAindIC)) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VparICT(FAindIC))) &
							.eqv. .true.) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VparICT(:)) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VparMBT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND FAindIC= ', FAindIC, &
								' IN VELOCITY DISTRIBUTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((VxMB(1) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VxICT(FAindIC)) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VxICT(FAindIC))) &
							.eqv. .true.) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VxICT(:)) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VxMBT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND FAindIC= ', FAindIC, &
								' IN VELOCITY DISTRIBUTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((VyMB(1) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VyICT(FAindIC)) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VyICT(FAindIC))) &
							.eqv. .true.) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VyICT(:)) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VyMBT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND FAindIC= ', FAindIC, &
								' IN VELOCITY DISTRIBUTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((VzMB(1) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VzICT(FAindIC)) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VzICT(FAindIC))) &
							.eqv. .true.) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VzICT(:)) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VzMBT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND FAindIC= ', FAindIC, &
								' IN VELOCITY DISTRIBUTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						! ----------------------------------------------------

					end do
				end do

				ICflag(1)= 0d0

			end do
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

		! PRINT TERMINAL EXPORT PARAMETERS:

		! do s= 1, Stot, 1
! 			do f= 1, SpecieT(s)%NfT(1), 1
! 				do Qind= NqLB(1), NqUB(1), 1
! 					if (rank == 0) then
! 						! write(*, *) 'Vperp1ICT= ', &
! ! 						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%Vperp1ICT(:)
! ! 						write(*, *) 'Vperp2ICT= ', &
! ! 						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%Vperp2ICT(:)
! ! 						write(*, *) 'VperpICT= ', &
! ! 						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VperpICT(:)
! 						write(*, *) 'VparICT= ', &
! 						sum(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VparICT(:))/SpecieT(s)%FluxTubeT(f)%NsT(1)
! 						! write(*, *) 'VxICT= ', &
! ! 						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VxICT(:)
! ! 						write(*, *) 'VyICT= ', &
! ! 						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VyICT(:)
! ! 						write(*, *) 'VzICT= ', &
! ! 						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VzICT(:)
! 					end if
! 				end do
! 			end do
! 		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				if (rank == 0) then
					call cpu_time(S6End)
					write(S6string, '(i10)')  nint(S6End)
					write(*, *) trim('%% 6- RANK= ' // adjustl(rankstring)) // &
						trim(', REAL-TIME= ' // adjustl(S6string)) // &
						trim(' s. VELOCITY DISTRIBUTION COMPLETE %%')
				end if
			end do
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine VelocityDistributionSub

end module VelocityDistribution
