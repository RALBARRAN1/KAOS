module KineticSolver

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8 3D KINETIC SOLVER:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use DataTypeAllocB
use IonNeutralCollisions
use BoundaryConditions
!use ConvecConfigGridGenerator
!use ConvecVelGridGenerator
use ConfigGridGenerator
use VelGridGenerator
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
							(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1)))) then
							do Qind= NqLB(1), NqUB(1), 1
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
								if ((nn /= 1d0) .and. (n >= (nn- 2)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1)+ 1d0) .and. &
									(n <= (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1))) then

									! ----------------------------------------------------

									do Qind= NqLB(1), NqUB(1), 1

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

					! UPDATE PHASE-SPACE GRID WITHOUT CONVECTING FLUX-TUBES:

					if (SpecieT(s)%FluxTubeT(f)%CONVECTIONflagT(1) == 0) then
						do nn= 2, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1

							SpecieT(s)%FluxTubeT(f)%qGCT(nn, :)= SpecieT(s)%FluxTubeT(f)%qGCT(1, :)
							SpecieT(s)%FluxTubeT(f)%hqCT(nn, :)= SpecieT(s)%FluxTubeT(f)%hqCT(1, :)
							SpecieT(s)%FluxTubeT(f)%dpCT(nn, :)= SpecieT(s)%FluxTubeT(f)%dpCT(1, :)
							SpecieT(s)%FluxTubeT(f)%dqCT(nn, :)= SpecieT(s)%FluxTubeT(f)%dqCT(1, :)
							SpecieT(s)%FluxTubeT(f)%dphiCT(nn, :)= SpecieT(s)%FluxTubeT(f)%dphiCT(1, :)
							SpecieT(s)%FluxTubeT(f)%rGCT(nn, :)= SpecieT(s)%FluxTubeT(f)%rGCT(1, :)
							SpecieT(s)%FluxTubeT(f)%phiGCT(nn, :)= SpecieT(s)%FluxTubeT(f)%phiGCT(1, :)
							SpecieT(s)%FluxTubeT(f)%thetaGCT(nn, :)= SpecieT(s)%FluxTubeT(f)%thetaGCT(1, :)
							SpecieT(s)%FluxTubeT(f)%ellGCT(nn, :)= SpecieT(s)%FluxTubeT(f)%ellGCT(1, :)
							SpecieT(s)%FluxTubeT(f)%qGLT(nn, :)= SpecieT(s)%FluxTubeT(f)%qGLT(1, :)
							SpecieT(s)%FluxTubeT(f)%qGHT(nn, :)= SpecieT(s)%FluxTubeT(f)%qGHT(1, :)
							SpecieT(s)%FluxTubeT(f)%pGCT(nn, :)= SpecieT(s)%FluxTubeT(f)%pGCT(1, :)
							SpecieT(s)%FluxTubeT(f)%d3xCT(nn, :)= SpecieT(s)%FluxTubeT(f)%d3xCT(1, :)
							SpecieT(s)%FluxTubeT(f)%TsPerpT(nn, :)= SpecieT(s)%FluxTubeT(f)%TsPerpT(1, :)
							SpecieT(s)%FluxTubeT(f)%TsParT(nn, :)= SpecieT(s)%FluxTubeT(f)%TsParT(1, :)
							SpecieT(s)%FluxTubeT(f)%TsT(nn, :)= SpecieT(s)%FluxTubeT(f)%TsT(1, :)

							if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
								SpecieT(s)%FluxTubeT(f)%nsnormCNeutT(nn, :)= SpecieT(s)%FluxTubeT(f)%nsnormCNeutT(1, :)
							end if

						end do
					end if

					! ----------------------------------------------------

					! UPDATE PHASE-SPACE GRID WITH CONVECTING FLUX-TUBES:

					if (SpecieT(s)%FluxTubeT(f)%CONVECTIONflagT(1) == 1) then
						do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
							if ((n /= 1) .and. (nn /= 1) .and. &
								(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1))) then

									! Translate Q cells in L-shell (p, phi)

									!call ConvecConfigGridGeneratorSub
									!call ConvecVelGridGeneratorSub

									!call ConfigGridGeneratorSub
									!call VelGridGeneratorSub

									! ----------------------------------------------------

									! SET BOUNDARY GHOST CELL DENSITIES:

									! ----------------------------------------------------

									if (SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 0) then
										SpecieT(s)%FluxTubeT(f)%d3xCLBT(1)= &
											SpecieT(s)%FluxTubeT(f)%QCell0T(NqLB(1)+ 1)%d3xC0T(1)
										SpecieT(s)%FluxTubeT(f)%sigmaLBT(1)= &
											SpecieT(s)%FluxTubeT(f)%QCell0T(NqLB(1)+ 1)%d3xC0T(1)/ &
											(SpecieT(s)%FluxTubeT(f)%QCell0T(NqLB(1)+ 1)%hqC0T(1)* &
											SpecieT(s)%FluxTubeT(f)%QCell0T(NqLB(1)+ 1)%dqC0T(1))
									end if

									if (SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 1) then
										SpecieT(s)%FluxTubeT(f)%d3xCLBT(1)= &
											SpecieT(s)%FluxTubeT(f)%QCell0T(NqLB(1)+ 1)%d3xC0T(1)
										SpecieT(s)%FluxTubeT(f)%sigmaLBT(1)= &
											SpecieT(s)%FluxTubeT(f)%QCell0T(NqLB(1)+ 1)%d3xC0T(1)/ &
											(SpecieT(s)%FluxTubeT(f)%QCell0T(NqLB(1)+ 1)%hqC0T(1)* &
											SpecieT(s)%FluxTubeT(f)%QCell0T(NqLB(1)+ 1)%dqC0T(1))
									end if

									if (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 0) then
										SpecieT(s)%FluxTubeT(f)%d3xCUBT(1)= &
											SpecieT(s)%FluxTubeT(f)%QCell0T(NqUB(1)+ 1)%d3xC0T(1)
										SpecieT(s)%FluxTubeT(f)%sigmaUBT(1)= &
											SpecieT(s)%FluxTubeT(f)%QCell0T(NqUB(1)+ 1)%d3xC0T(1)/ &
											(SpecieT(s)%FluxTubeT(f)%QCell0T(NqUB(1)+ 1)%hqC0T(1)* &
											SpecieT(s)%FluxTubeT(f)%QCell0T(NqUB(1)+ 1)%dqC0T(1))
									end if

									if (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 1) then
										SpecieT(s)%FluxTubeT(f)%d3xCUBT(1)= &
											SpecieT(s)%FluxTubeT(f)%QCell0T(NqUB(1)+ 1)%d3xC0T(1)
										SpecieT(s)%FluxTubeT(f)%sigmaUBT(1)= &
											SpecieT(s)%FluxTubeT(f)%QCell0T(NqUB(1)+ 1)%d3xC0T(1)/ &
											(SpecieT(s)%FluxTubeT(f)%QCell0T(NqUB(1)+ 1)%hqC0T(1)* &
											SpecieT(s)%FluxTubeT(f)%QCell0T(NqUB(1)+ 1)%dqC0T(1))
									end if

									SpecieT(s)%FluxTubeT(f)%LBNominalDensityT(1)= nsnormCLB(1)
									SpecieT(s)%FluxTubeT(f)%UBNominalDensityT(1)= nsnormCUB(1)

									! ----------------------------------------------------

									! DIAGNOSTIC FLAGS CONSISTENT LOWER AND UPPER BOUNDARY GHOST CELL DENSITY:

									if (SpecieT(1)%FluxTubeT(1)%LBCONDITIONflagT(1) == 1) then
										if (nsnormCLB(1) == 0d0) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' ZERO LOWER BOUNDARY DENSITY= ', nsnormCLB(1), &
												' FOR SPECIE= ', s, ', AND FLUX TUBE= ', f, &
												' IN DENSITY PROFILE A SUBROUTINE' // achar(27) // '[0m.'
										end if

										if (SpecieT(s)%FluxTubeT(f)%LBNominalDensityT(1) == 0d0) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NOMINAL LOWER BOUNDARY DENSITY= ', SpecieT(s)%FluxTubeT(f)%LBNominalDensityT(1), &
												' FOR SPECIE= ', s, ', AND FLUX TUBE= ', f, &
												' IN DENSITY PROFILE A SUBROUTINE' // achar(27) // '[0m.'
										end if
									end if

									if (SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1) == 1) then
										if (rank == 0) then
											do Qind= NqLB(1), NqUB(1)- 1, 1
												if (SpecieT(s)%FluxTubeT(f)%NsFARpT(Qind) <= SpecieT(s)%FluxTubeT(f)%NsFARpT(Qind+ 1)) then
													write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
														' NON-DECREASING PARTICLE NUMBERS WITH ALTITUDE', &
														', NsFARpT(Qind)= ', SpecieT(s)%FluxTubeT(f)%NsFARpT(Qind), &
														', NsFARpT(Qind+ 1)= ', SpecieT(s)%FluxTubeT(f)%NsFARpT(Qind+ 1), &
														' FOR SPECIE= ', s, ', AND FLUX TUBE= ', f, &
														' IN DENSITY PROFILE A SUBROUTINE' // achar(27) // '[0m.'
												end if
											end do
										end if
									end if

									! ----------------------------------------------------

									! RE-INDEX CONFIGURATION SPACE GRID FOR A NON-COMPUTATIONAL LOWER BOUNDARY GHOST CELL:

									! ----------------------------------------------------

									do Qind= NqLB(1), NqUB(1), 1

										SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%qGC0T(1)
										SpecieT(s)%FluxTubeT(f)%hqCT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%hqC0T(1)
										SpecieT(s)%FluxTubeT(f)%dpCT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%dpC0T(1)
										SpecieT(s)%FluxTubeT(f)%dqCT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%dqC0T(1)
										SpecieT(s)%FluxTubeT(f)%dphiCT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%dphiC0T(1)
										SpecieT(s)%FluxTubeT(f)%rGCT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%rGC0T(1)
										SpecieT(s)%FluxTubeT(f)%phiGCT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%phiGC0T(1)
										SpecieT(s)%FluxTubeT(f)%thetaGCT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%thetaGC0T(1)
										SpecieT(s)%FluxTubeT(f)%ellGCT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%ellGC0T(1)
										SpecieT(s)%FluxTubeT(f)%qGLT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%qGL0T(1)
										SpecieT(s)%FluxTubeT(f)%qGHT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%qGH0T(1)
										SpecieT(s)%FluxTubeT(f)%pGCT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%pGC0T(1)
										SpecieT(s)%FluxTubeT(f)%d3xCT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%d3xC0T(1)
										SpecieT(s)%FluxTubeT(f)%TsPerpT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%TsPerp0T(1)
										SpecieT(s)%FluxTubeT(f)%TsParT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%TsPar0T(1)
										SpecieT(s)%FluxTubeT(f)%TsT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%Ts0T(1)

										if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
											SpecieT(s)%FluxTubeT(f)%nsnormCNeutT(nn, Qind)= nsnormCNeut0(Qind+ 1)
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR CONSISTENT PHI VALUE:

										if (SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(1) /= SpecieT(s)%FluxTubeT(f)%QCell0T(1)%phiGC0T(1)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT INITIAL phiGC0T VALUE= ', &
												SpecieT(s)%FluxTubeT(f)%QCell0T(1)%phiGC0T(1), ' AND phiGC0T VALUE= ', &
												SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(1), &
												' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND Q CELL= ', Qind, &
												' IN KINETIC SOLVER SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end do

									! ----------------------------------------------------

							end if
						end do
					end if

					!do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
					!	if ((n /= 1) .and. (nn /= 1) .and. (n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1))) then
					! update grid with conservation of f and particle number,
					! update LBNominalDensityT and UBNominalDensityT
					! update all current particle positions and velocities with convection and with betatron acceleration
					!do j= 1, (NsTK(1)- dNsTK2(1)- dNsTK3(1)), 1
					! xN(j), yN(j), zN(j), Vperp1N(j), Vperp2N(j), VperpN(j), VxN(j), VyN(j), VzN(j) for ions and ENAs.
					! Note: all new LB and UB injected particles are reset on correct Lshell in KineticSolverB given
					! SpecieT(s)%FluxTubeT(f)%pGCT(nn, Qind) and SpecieT(s)%FluxTubeT(f)%phiGCT(nn, Qind)
					!1- Conserve f and Update GRID
					!2- Update LB and UB densities and ensure robust statistics
					!3- Translate all particle positions in (p, q, phi)
					!4- Update Vperp values with betatron ACCELERATION
					!5- On new field line, update AEAmagN(j), AGmagN(j), AEPmagN(j)
					!6- In KineticUpdateA, KineticUpdateB, KineticUpdateC add parallel centrifual acceleration

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
