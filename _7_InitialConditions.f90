module InitialConditions

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	7- INITIAL CONDITIONS:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use RK4DipolePolynomialSolver

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! ALLOCATE INITIAL CONDITIONS TO ALL PARTICLES:

	subroutine InitialConditionsSub

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				! ----------------------------------------------------

          allocate(SpecieT(s)%FluxTubeT(f)%x0T(SpecieT(s)%FluxTubeT(f)%NsT(1)), &
						SpecieT(s)%FluxTubeT(f)%y0T(SpecieT(s)%FluxTubeT(f)%NsT(1)), &
						SpecieT(s)%FluxTubeT(f)%z0T(SpecieT(s)%FluxTubeT(f)%NsT(1)), &
						SpecieT(s)%FluxTubeT(f)%r0T(SpecieT(s)%FluxTubeT(f)%NsT(1)), &
						SpecieT(s)%FluxTubeT(f)%theta0T(SpecieT(s)%FluxTubeT(f)%NsT(1)), &
						SpecieT(s)%FluxTubeT(f)%phi0T(SpecieT(s)%FluxTubeT(f)%NsT(1)), &
						SpecieT(s)%FluxTubeT(f)%q0T(SpecieT(s)%FluxTubeT(f)%NsT(1)), &
						SpecieT(s)%FluxTubeT(f)%p0T(SpecieT(s)%FluxTubeT(f)%NsT(1)), &
						SpecieT(s)%FluxTubeT(f)%Vperp10T(SpecieT(s)%FluxTubeT(f)%NsT(1)), &
						SpecieT(s)%FluxTubeT(f)%Vperp20T(SpecieT(s)%FluxTubeT(f)%NsT(1)), &
						SpecieT(s)%FluxTubeT(f)%Vperp0T(SpecieT(s)%FluxTubeT(f)%NsT(1)), &
						SpecieT(s)%FluxTubeT(f)%Vpar0T(SpecieT(s)%FluxTubeT(f)%NsT(1)), &
						SpecieT(s)%FluxTubeT(f)%Vx0T(SpecieT(s)%FluxTubeT(f)%NsT(1)), &
						SpecieT(s)%FluxTubeT(f)%Vy0T(SpecieT(s)%FluxTubeT(f)%NsT(1)), &
						SpecieT(s)%FluxTubeT(f)%Vz0T(SpecieT(s)%FluxTubeT(f)%NsT(1)))

				! ----------------------------------------------------

        do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
					do FAindIC= 1, SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1), 1

						SpecieT(s)%FluxTubeT(f)%x0T(sum(SpecieT(s)%FluxTubeT(f)% &
							NsFApT(1:Qind- 1))+ FAindIC)= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%xfinalICT(FAindIC)
						SpecieT(s)%FluxTubeT(f)%y0T(sum(SpecieT(s)%FluxTubeT(f)% &
							NsFApT(1:Qind- 1))+ FAindIC)= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%yfinalICT(FAindIC)
						SpecieT(s)%FluxTubeT(f)%z0T(sum(SpecieT(s)%FluxTubeT(f)% &
							NsFApT(1:Qind- 1))+ FAindIC)= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%zfinalICT(FAindIC)
						SpecieT(s)%FluxTubeT(f)%r0T(sum(SpecieT(s)%FluxTubeT(f)% &
							NsFApT(1:Qind- 1))+ FAindIC)= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%rfinalICT(FAindIC)
						SpecieT(s)%FluxTubeT(f)%theta0T(sum(SpecieT(s)%FluxTubeT(f)% &
							NsFApT(1:Qind- 1))+ FAindIC)= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%thetafinalICT(FAindIC)
						SpecieT(s)%FluxTubeT(f)%phi0T(sum(SpecieT(s)%FluxTubeT(f)% &
							NsFApT(1:Qind- 1))+ FAindIC)= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%phifinalICT(FAindIC)
						SpecieT(s)%FluxTubeT(f)%q0T(sum(SpecieT(s)%FluxTubeT(f)% &
							NsFApT(1:Qind- 1))+ FAindIC)= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qfinalICT(FAindIC)
						SpecieT(s)%FluxTubeT(f)%p0T(sum(SpecieT(s)%FluxTubeT(f)% &
							NsFApT(1:Qind- 1))+ FAindIC)= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%pfinalICT(FAindIC)
						SpecieT(s)%FluxTubeT(f)%Vperp10T(sum(SpecieT(s)%FluxTubeT(f)% &
							NsFApT(1:Qind- 1))+ FAindIC)= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%Vperp1ICT(FAindIC)
						SpecieT(s)%FluxTubeT(f)%Vperp20T(sum(SpecieT(s)%FluxTubeT(f)% &
							NsFApT(1:Qind- 1))+ FAindIC)= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%Vperp2ICT(FAindIC)
						SpecieT(s)%FluxTubeT(f)%Vperp0T(sum(SpecieT(s)%FluxTubeT(f)% &
							NsFApT(1:Qind- 1))+ FAindIC)= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VperpICT(FAindIC)
						SpecieT(s)%FluxTubeT(f)%Vpar0T(sum(SpecieT(s)%FluxTubeT(f)% &
							NsFApT(1:Qind- 1))+ FAindIC)= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VparICT(FAindIC)
						SpecieT(s)%FluxTubeT(f)%Vx0T(sum(SpecieT(s)%FluxTubeT(f)% &
							NsFApT(1:Qind- 1))+ FAindIC)= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VxICT(FAindIC)
						SpecieT(s)%FluxTubeT(f)%Vy0T(sum(SpecieT(s)%FluxTubeT(f)% &
							NsFApT(1:Qind- 1))+ FAindIC)= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VyICT(FAindIC)
						SpecieT(s)%FluxTubeT(f)%Vz0T(sum(SpecieT(s)%FluxTubeT(f)% &
							NsFApT(1:Qind- 1))+ FAindIC)= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%VzICT(FAindIC)

					end do
				end do

				! ----------------------------------------------------

			end do
		end do

		! ----------------------------------------------------

    ! SET ALL INITIAL CONDITION DIAGNOSTICS:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1

					! ----------------------------------------------------

					! DIAGNOSTIC FLAG FOR ALL INITIAL VALUES WITHIN PHASE-SPACE GRID:

! 					get max limits to compare here

					! if ((SpecieT(s)%FluxTubeT(f)%q0T(j) <= &
! 						SpecieT(s)%FluxTubeT(f)%qGLT(1, 1)) .or. &
! 						(SpecieT(s)%FluxTubeT(f)%q0T(j) > SpecieT(s)%FluxTubeT(f)%QCellT( &
! 						((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))%qGHT(1))) then
! 						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
! 							' q0T= ', SpecieT(s)%FluxTubeT(f)%q0T(j), ' VALUE OUT OF CONFIG-SPACE', &
! 							' GRID FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND PARTICLE= ', j, &
! 							' IN INITIAL CONDITIONS SUBROUTINE' // achar(27) // '[0m.'
! 					end if
!
! 					if ((SpecieT(s)%FluxTubeT(f)%Vperp0T(j) <= &
! 						SpecieT(s)%FluxTubeT(f)%QCellT(1)%VCellT(1, 1)%VperpGLT(1)) .or. &
! 						(SpecieT(s)%FluxTubeT(f)%Vperp0T(j) > SpecieT(s)%FluxTubeT(f)% &
! 						QCellT(1)%VCellT(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
! 							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1))%VperpGHT(1))) then
! 						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
! 							' Vperp0T= ', SpecieT(s)%FluxTubeT(f)%Vperp0T(j), ' VALUE OUT OF', &
! 							' VEL-SPACE GRID FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
! 							', AND PARTICLE= ', j, ' IN INITIAL CONDITIONS SUBROUTINE' &
! 							// achar(27) // '[0m.'
! 					end if
!
! 					if ((SpecieT(s)%FluxTubeT(f)%Vpar0T(j) <= &
! 						SpecieT(s)%FluxTubeT(f)%QCellT(1)%VCellT(1, 1)%VparGLT(1)) .or. &
! 						(SpecieT(s)%FluxTubeT(f)%Vperp0T(j) > SpecieT(s)%FluxTubeT(f)% &
! 						QCellT(1)%VCellT(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
! 							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1))%VparGHT(1))) then
! 						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
! 							' Vpar0T= ', SpecieT(s)%FluxTubeT(f)%Vpar0T(j), ' VALUE OUT OF', &
! 							' VEL-SPACE GRID FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
! 							', AND PARTICLE= ', j, ' IN INITIAL CONDITIONS SUBROUTINE' &
! 							// achar(27) // '[0m.'
! 					end if

					! ----------------------------------------------------

					! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

          if ((isnan(real(SpecieT(s)%FluxTubeT(f)%x0T(j))) .eqv. .true.) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%x0T(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' x0T HAS', &
							' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', AND PARTICLE= ', j, ' IN INITIAL CONDITIONS SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if ((isnan(real(SpecieT(s)%FluxTubeT(f)%y0T(j))) .eqv. .true.) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%y0T(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' y0T HAS', &
							' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', AND PARTICLE= ', j, ' IN INITIAL CONDITIONS SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if ((isnan(real(SpecieT(s)%FluxTubeT(f)%z0T(j))) .eqv. .true.) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%z0T(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' z0T HAS', &
							' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', AND PARTICLE= ', j, ' IN INITIAL CONDITIONS SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if ((isnan(real(SpecieT(s)%FluxTubeT(f)%r0T(j))) .eqv. .true.) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%r0T(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' r0T HAS', &
							' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', AND PARTICLE= ', j, ' IN INITIAL CONDITIONS SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if ((isnan(real(SpecieT(s)%FluxTubeT(f)%theta0T(j))) .eqv. .true.) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%theta0T(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' theta0T HAS', &
							' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', AND PARTICLE= ', j, ' IN INITIAL CONDITIONS SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if ((isnan(real(SpecieT(s)%FluxTubeT(f)%phi0T(j))) .eqv. .true.) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%phi0T(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' phi0T HAS', &
							' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', AND PARTICLE= ', j, ' IN INITIAL CONDITIONS SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if ((isnan(real(SpecieT(s)%FluxTubeT(f)%q0T(j))) .eqv. .true.) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%q0T(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' q0T HAS', &
							' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', AND PARTICLE= ', j, ' IN INITIAL CONDITIONS SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if ((isnan(real(SpecieT(s)%FluxTubeT(f)%p0T(j))) .eqv. .true.) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%p0T(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' p0T HAS', &
							' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', AND PARTICLE= ', j, ' IN INITIAL CONDITIONS SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if ((isnan(real(SpecieT(s)%FluxTubeT(f)%Vperp10T(j))) .eqv. .true.) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%Vperp10T(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp10T HAS', &
							' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', AND PARTICLE= ', j, ' IN INITIAL CONDITIONS SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if ((isnan(real(SpecieT(s)%FluxTubeT(f)%Vperp20T(j))) .eqv. .true.) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%Vperp20T(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp20T HAS', &
							' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', AND PARTICLE= ', j, ' IN INITIAL CONDITIONS SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if ((isnan(real(SpecieT(s)%FluxTubeT(f)%Vperp0T(j))) .eqv. .true.) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%Vperp0T(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp0T HAS', &
							' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', AND PARTICLE= ', j, ' IN INITIAL CONDITIONS SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if ((isnan(real(SpecieT(s)%FluxTubeT(f)%Vpar0T(j))) .eqv. .true.) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%Vpar0T(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vpar0T HAS', &
							' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', AND PARTICLE= ', j, ' IN INITIAL CONDITIONS SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if ((isnan(real(SpecieT(s)%FluxTubeT(f)%Vx0T(j))) .eqv. .true.) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%Vx0T(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vx0T HAS', &
							' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', AND PARTICLE= ', j, ' IN INITIAL CONDITIONS SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if ((isnan(real(SpecieT(s)%FluxTubeT(f)%Vy0T(j))) .eqv. .true.) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%Vy0T(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vy0T HAS', &
							' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', AND PARTICLE= ', j, ' IN INITIAL CONDITIONS SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if ((isnan(real(SpecieT(s)%FluxTubeT(f)%Vz0T(j))) .eqv. .true.) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%Vz0T(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vz0T HAS', &
							' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', AND PARTICLE= ', j, ' IN INITIAL CONDITIONS SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

        end do

			! ----------------------------------------------------

			end do
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

    ! PRINT TERMINAL EXPORT PARAMETERS:

		! do s= 1, Stot, 1
! 			do f= 1, SpecieT(s)%NfT(1), 1
! 				if (rank == 0) then
! 					write(*, *) 'x0T sf= ', s, f, SpecieT(s)%FluxTubeT(f)%x0T(:)
! 					write(*, *) 'y0T sf= ', s, f, SpecieT(s)%FluxTubeT(f)%y0T(:)
! 					write(*, *) 'z0T sf= ', s, f, SpecieT(s)%FluxTubeT(f)%z0T(:)
! 					write(*, *) 'r0T sf= ', s, f, SpecieT(s)%FluxTubeT(f)%r0T(:)
! 					write(*, *) 'theta0T sf= ', s, f, SpecieT(s)%FluxTubeT(f)%theta0T(:)
! 					write(*, *) 'phi0T sf= ', s, f, SpecieT(s)%FluxTubeT(f)%phi0T(:)
! 					write(*, *) 'q0T sf= ', s, f, &
! 						SpecieT(s)%FluxTubeT(f)%q0T(SpecieT(s)%FluxTubeT(f)%NsT(1))
! 					write(*, *) 'p0T sf= ', s, f, SpecieT(s)%FluxTubeT(f)%p0T(:)
! 					write(*, *) 'phid0T sf= ', s, f, SpecieT(s)%FluxTubeT(f)%phid0T(:)
! 					write(*, *) 'Vperp10T sf= ', s, f, SpecieT(s)%FluxTubeT(f)%Vperp10T(:)
! 					write(*, *) 'Vperp20T sf= ', s, f, SpecieT(s)%FluxTubeT(f)%Vperp20T(:)
! 					write(*, *) 'Vperp0T sf= ', s, f, &
! 						SpecieT(s)%FluxTubeT(f)%Vperp0T(:)
! 					write(*, *) 'Vpar0T sf= ', s, f, &
! 						SpecieT(s)%FluxTubeT(f)%Vpar0T(:)
! 					write(*, *) 'Vx0T sf= ', s, f, SpecieT(s)%FluxTubeT(f)%Vx0T(:)
! 					write(*, *) 'Vy0T sf= ', s, f, SpecieT(s)%FluxTubeT(f)%Vy0T(:)
! 					write(*, *) 'Vz0T sf= ', s, f, SpecieT(s)%FluxTubeT(f)%Vz0T(:)
! 				end if
! 			end do
! 		end do

		! do s= 1, Stot, 1
! 			do f= 1, SpecieT(s)%NfT(1), 1
! 				if (SpecieT(s)%FluxTubeT(f)%QCellT(1)%PHASEDISTRIBflagT(1) == 1) then
! 					if (rank == 0) then
! 						write(*, *) 'NqqReNormE1T= ', SpecieT(s)%FluxTubeT(f)%NqqReNormE1T(1)
! 						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
! 							write(*, *) 'NqReNorm1pT= ', SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqReNorm1pT(1)
! 							write(*, *) 'NqReNorm1T= ', SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqReNorm1T(1)
! 							do Vperpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), 1
! 								do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1
! 									write(*, *) 'NphVperpReNormE1T= ', &
! 										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT &
! 											(Vperpind, Vparind)%NphVperpReNormE1T(1)
! 									write(*, *) 'NphVparReNormE1T= ', &
! 										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT &
! 											(Vperpind, Vparind)%NphVparReNormE1T(1)
! 									write(*, *) 'NphVperpVparReNormE1T= ', &
! 										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT &
! 											(Vperpind, Vparind)%NphVperpVparReNormE1T(1)
! 									write(*, *) 'NphReNorm1pT= ', &
! 										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT &
! 											(Vperpind, Vparind)%NphReNorm1pT(1)
! 									write(*, *) 'NphReNorm1T= ', &
! 										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT &
! 											(Vperpind, Vparind)%NphReNorm1T(1)
! 								end do
! 							end do
! 						end do
! 					end if
! 				end if
! 			end do
! 		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				if (rank == 0) then
					call cpu_time(S7End)
					write(S7string, '(i10)')  nint(S7End)
					write(*, *) trim('%% 7- RANK= ' // adjustl(rankstring)) // &
						trim(', REAL-TIME= ' // adjustl(S7string)) // &
						trim(' s. INITIAL CONDITIONS COMPLETE %%')
				end if
			end do
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine InitialConditionsSub

end module InitialConditions
