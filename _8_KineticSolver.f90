module KineticSolver

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8- 3D KINETIC SOLVER:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use DataTypeAllocB
use IonNeutralCollisions
use BoundaryConditions
use Convection
use KineticSolverA
use KineticSolverB

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! USE ALL PARTICLE INITIAL PHASE-SPACE CONDITIONS TO GET POSITION-DEPENDENT ACCELERATION VALUES
! AND INTEGRATE THEM BY RK4 METHOD FOR 3D CARTESIAN PHASE-SPACE VALUES:

	subroutine KineticSolverSub

		! ----------------------------------------------------

		! BEGIN KINETIC SOLVER:

		! ----------------------------------------------------

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				! ----------------------------------------------------

				! ALLOCATE PARTICLE NUMBER-INDEPENDENT DERIVED DATA TYPES:

				call DataTypeAllocBSub

				! ----------------------------------------------------

				! BEGIN KINETIC SOLVER ROUTINE:

				do n= 1, SpecieT(s)%FluxTubeT(f)%NtT(1), dt

					! ----------------------------------------------------

					! UPDATE ELECTRON TEMPERATURE IN TIME:

					do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
						if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
							(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn)))) then
							do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
								if ((n == 1) .and. (nn == 1)) then
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn)= SpecieT(s)%FluxTubeT(f)%Te0T(Qind)
								else
									if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn) < dNTeEND) then
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn- 1)+ dNTe
									end if
									if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn) > dNTeEND) then
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn- 1)
									end if
								end if
							end do
						end if
					end do

					! ----------------------------------------------------

					! PERFORM ENA PRODUCTION AND PARTICLE DE-MAGNETIZATION ON APPROPRIATE COLLISION FREQUENCY:

					if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
						if (n /= 1d0) then
							do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1), 1
								if ((nn /= 1d0) .and. (n >= (nn- 2)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn)+ 1d0) .and. &
									(n <= (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn))) then

									! ----------------------------------------------------

									do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

										! ----------------------------------------------------

										if (SpecieT(s)%FluxTubeT(f)%nuIonNeutPoiT(nn- 1, Qind) > 0d0) then

											! ----------------------------------------------------

											ENAcount(1)= 0d0
											ENAloopj: do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
												if ((ENAflag(j) .eqv. .false.) .and. (Qindk1(j) == Qind)) then
													do while (ENAcount(1) <= nint(SpecieT(s)%FluxTubeT(f)%nuIonNeutPoiT(nn- 1, Qind)))

														! ----------------------------------------------------

														call IonNeutralCollisionsSub

														! ----------------------------------------------------

														if (ENAcount(1) < nint(SpecieT(s)%FluxTubeT(f)%nuIonNeutPoiT(nn- 1, Qind))) then
															ENAcount(1)= ENAcount(1)+ 1d0
														end if
														if (ENAcount(1) == nint(SpecieT(s)%FluxTubeT(f)%nuIonNeutPoiT(nn- 1, Qind))) then
															exit ENAloopj
														end if

														! ----------------------------------------------------

														! DIAGNOSTIC FLAGS FOR ENA PRODUCTION NUMBER:

														if (ENAcount(1) > nint(SpecieT(s)%FluxTubeT(f)%nuIonNeutPoiT(nn- 1, Qind))) then
															write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
																' INCONSISTENT ENA PRODUCTION NUMBER ', &
																'	ENAcount= ', ENAcount(1), &
																' AND nuIonNeutPoiT= ', nint(SpecieT(s)%FluxTubeT(f)%nuIonNeutPoiT(nn- 1, Qind)), &
																' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
																' PARTICLE NUMBER= ', j, ', AND TIME-STEP= ', n, ' IN KINETIC SOLVER ', &
																' SUBROUTINE' // achar(27) // '[0m.'
														end if

														! ----------------------------------------------------

														cycle ENAloopj

													end do
												else
													cycle ENAloopj
												end if
											end do ENAloopj

											! ----------------------------------------------------

										end if

										! ----------------------------------------------------

									end do

									! ----------------------------------------------------

								end if
							end do
						end if
					end if

					! ----------------------------------------------------

					! COMPUTE PARTICLE LOSS AND LOWER/UPPER BOUNDARY INJECTION:

					call BoundaryConditionsSub

					! ----------------------------------------------------

					! UPDATE TIME PARAMETER AT INITIAL AND ALL OTHER TIME:

					if (n == 1) then
						Time(1)= (n- 1d0)*SpecieT(s)%FluxTubeT(f)%hT(1)
					else if (n /= 1) then
						Time(1)= TimeN(1)
					end if

					TimeN(1)= (n+ dt)*SpecieT(s)%FluxTubeT(f)%hT(1)

					! ----------------------------------------------------

					! DIAGNOSTICS FOR CONSISTENT INITIAL AND FINAL TIME PARAMETER:

					if ((n == 1) .and. (Time(1) /= 0d0)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INITIAL TIME DOES', &
							' NOT EQUAL A+ 1 FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND TIME-STEP= ', &
							n, ' IN KINETIC SOLVER SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((n == SpecieT(s)%FluxTubeT(f)%NtT(1)) .and. &
						(abs(Time(1)- SpecieT(s)%FluxTubeT(f)%BT(1)) > 1d-6)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' FINAL TIME DOES', &
							' NOT EQUAL B FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND TIME-STEP= ', &
							n, ' IN KINETIC SOLVER SUBROUTINE' // achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

					! CONVECT FLUX-TUBE GRID AND TRANSLATE PARTICLES IN EXB:

					call ConvectionSub

					! ----------------------------------------------------

					! RUN KINETIC CODE (3D CARTESIAN RK4) AT INITIAL TIME AND ALL OTHER TIME:

					call KineticSolverASub
					call KineticSolverBSub

					! ----------------------------------------------------

				end do

				! ----------------------------------------------------

			end do
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				if (rank == 0) then
					call cpu_time(S8End)
					write(S8string, '(i10)')  nint(S8End)
					write(*, *) trim('%% 8- RANK= ' // adjustl(rankstring)) // &
						trim(', REAL-TIME= ' // adjustl(S8string)) // &
						trim(' s. KINETIC SOLVER COMPLETE %%')
				end if
			end do
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine KineticSolverSub

end module KineticSolver
