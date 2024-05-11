module SecondPerpIonMoment

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.6.5 SECOND PERPENDICULAR ION MOMENT:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use NonUniform2DIntegrator
use NonUniform3DIonIntegrator

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! COMPUTE SECOND PERPENDICULAR ION MOMENT (PERPENDICULAR ENERGY) OF PHASE-SPACE DISTRIBUTION FUNCTIONS:

	subroutine SecondPerpIonMomentSub

		! ----------------------------------------------------

		if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

			do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
				if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
					(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacT(1:nn- 1))))) then

					! ----------------------------------------------------

					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

						! ----------------------------------------------------

						! COMPUTE ION PERPENDICULAR 1 PLASMA ENERGY (PERPENDICULAR 1 SECOND MOMENT):

						! ----------------------------------------------------

						if ((SpecieT(s)%FluxTubeT(f)% &
							PHASEENERGYPERPIONMOMENTflagT(1) == 1) .and. &
							(rank == 0)) then

							! ----------------------------------------------------

							allocate(ggg2Perp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), &
								SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
								SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), &
								SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1), &
								Sum12Perp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), &
								SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
								SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), &
								SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))

							allocate(II(1), MM(1))

							! ----------------------------------------------------

							do Vperp1ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), 1
								do Vperp2ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), 1
									do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1

										! ----------------------------------------------------

										if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn) == 0d0) then

											ggg2Perp(Vperp1ind, Vperp2ind, Vparind, nn)= 0d0

										else if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn) /= 0d0) then

											if (SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTCENTERflagT(1) == 0) then
												ggg2Perp(Vperp1ind, Vperp2ind, Vparind, nn)= &
													((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
													V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GCT(1))**2d0)* &
													SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
													V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn)
											end if

											if (SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTCENTERflagT(1) == 1) then
												ggg2Perp(Vperp1ind, Vperp2ind, Vparind, nn)= &
													((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
													V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GCT(1)- &
													SpecieT(s)%FluxTubeT(f)%M1Perp1phRT(nn, Qind))**2d0)* &
													SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
													V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn)

												if (((isnan(real(ggg2Perp(Vperp1ind, Vperp2ind, Vparind, nn))) &
													.eqv. .true.)) .or. (ggg2Perp(Vperp1ind, Vperp2ind, Vparind, nn) == 0d0)) then
													ggg2Perp(Vperp1ind, Vperp2ind, Vparind, nn)= &
														((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
														V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GCT(1))**2d0)* &
														SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
														V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn)
												end if

											end if
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

										if ((isnan(real(ggg2Perp(Vperp1ind, Vperp1ind, Vparind, nn))) &
											.eqv. .true.) .or. (size(ggg2Perp(:, :, :, :)) /= &
											(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1))* &
											(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1))* &
											(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1))* &
											(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' INTEGRAND HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
												s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', &
												Vperp1ind, ', Vperp2ind= ', &
												Vperp2ind, ', Vparind= ', Vparind, ', AND', &
												' STATISTICAL TIME-STEP= ', nn, ' IN SECOND', &
												' PERPENDICULAR 1 ION MOMENT SUBROUTINE' &
												// achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR CONSISTENT INTEGRANDS:

										if (((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn) /= 0d0) .and. &
											(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GCT(1) /= 0d0) .and. &
											(ggg2Perp(Vperp1ind, Vperp2ind, Vparind, nn) == 0d0)) .or. &
											((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn) == 0d0) .and. &
											(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GCT(1) == 0d0) .and. &
											(ggg2Perp(Vperp1ind, Vperp2ind, Vparind, nn) /= 0d0))) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' INCONSISTENT INTEGRAND TEST FOR SPECIE= ', &
												s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', &
												Vperp1ind, ', Vperp2ind= ', &
												Vperp2ind, ', Vparind= ', Vparind, ', AND', &
												' STATISTICAL TIME-STEP= ', nn, ' IN SECOND', &
												' PERPENDICULAR 1 ION MOMENT SUBROUTINE' &
												// achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end do
								end do
							end do

							! ----------------------------------------------------

							! Compute 3D velocity space integration
							call NonUniform3DIonIntegratorSub

							! Second perpendicular 1 ion distribution function moment [J]

							if (SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) == 0d0) then
								SpecieT(s)%FluxTubeT(f)%M2Perp1phRT(nn, Qind)= 0d0
							else if (SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) /= 0d0) then
								SpecieT(s)%FluxTubeT(f)%M2Perp1phRT(nn, Qind)= &
									(SpecieT(s)%msT(1)/2d0)* &
									(MM(1)/SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind))
							end if

							! ----------------------------------------------------

							! DIAGNOSTIC FLAGS FOR MOMENT CONSISTENCY:

							if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn) == 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M2Perp1phRT(nn, Qind) /= 0d0)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' INCONSISTENT M2Perp1phRT and NqRT VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN SECOND PERPENDICULAR 1 ION MOMENT', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((sum(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%F2PerpphRTp(:, :, :, nn)) == 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M2Perp1phRT(nn, Qind) /= 0d0)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' INCONSISTENT M2Perp1phRT and F2PerpphRTp SUMMATION FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN SECOND PERPENDICULAR 1 ION MOMENT', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if (((sum(ggg2Perp(:, :, :, nn)) /= 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M2Perp1phRT(nn, Qind) == 0d0)) .or. &
								((sum(ggg2Perp(:, :, :, nn)) == 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M2Perp1phRT(nn, Qind) /= 0d0))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' INCONSISTENT M2Perp1phRT and INTEGRAND SUMMATION FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN SECOND PERPENDICULAR 1 ION MOMENT', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							! ----------------------------------------------------

							deallocate(ggg2Perp, Sum12Perp)

							deallocate(II, MM)

							! ----------------------------------------------------

						end if

						! ----------------------------------------------------

					end do
				end if
			end do

			do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
				if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
					(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacT(1:nn- 1))))) then

					! ----------------------------------------------------

					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

						! ----------------------------------------------------

						! COMPUTE ION PERPENDICULAR 2 PLASMA ENERGY (PERPENDICULAR 1 SECOND MOMENT):

						! ----------------------------------------------------

						if ((SpecieT(s)%FluxTubeT(f)% &
							PHASEENERGYPERPIONMOMENTflagT(1) == 1) .and. &
							(rank == 0)) then

							! ----------------------------------------------------

							allocate(ggg2Perp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), &
								SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
								SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), &
								SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1), &
								Sum12Perp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), &
								SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
								SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), &
								SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))

							allocate(II(1), MM(1))

							! ----------------------------------------------------

							do Vperp1ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), 1
								do Vperp2ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), 1
									do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1

										! ----------------------------------------------------

										if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn) == 0d0) then

											ggg2Perp(Vperp1ind, Vperp2ind, Vparind, nn)= 0d0

										else if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn) /= 0d0) then

											if (SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTCENTERflagT(1) == 0) then
												ggg2Perp(Vperp1ind, Vperp2ind, Vparind, nn)= &
													((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
													V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GCT(1))**2d0)* &
													SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
													V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn)
											end if

											if (SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTCENTERflagT(1) == 1) then
												ggg2Perp(Vperp1ind, Vperp2ind, Vparind, nn)= &
													((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
													V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GCT(1)- &
													SpecieT(s)%FluxTubeT(f)%M1Perp2phRT(nn, Qind))**2d0)* &
													SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
													V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn)

												if (((isnan(real(ggg2Perp(Vperp1ind, Vperp2ind, Vparind, nn))) &
													.eqv. .true.)) .or. (ggg2Perp(Vperp1ind, Vperp2ind, Vparind, nn) == 0d0)) then
													ggg2Perp(Vperp1ind, Vperp2ind, Vparind, nn)= &
														((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
														V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GCT(1))**2d0)* &
														SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
														V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn)
												end if

											end if
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

										if ((isnan(real(ggg2Perp(Vperp1ind, Vperp1ind, Vparind, nn))) &
											.eqv. .true.) .or. (size(ggg2Perp(:, :, :, :)) /= &
											(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1))* &
											(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1))* &
											(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1))* &
											(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' INTEGRAND HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
												s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', &
												Vperp1ind, ', Vperp2ind= ', &
												Vperp2ind, ', Vparind= ', Vparind, ', AND', &
												' STATISTICAL TIME-STEP= ', nn, ' IN SECOND', &
												' PERPENDICULAR 2 ION MOMENT SUBROUTINE' &
												// achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR CONSISTENT INTEGRANDS:

										if (((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn) /= 0d0) .and. &
											(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GCT(1) /= 0d0) .and. &
											(ggg2Perp(Vperp1ind, Vperp2ind, Vparind, nn) == 0d0)) .or. &
											((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn) == 0d0) .and. &
											(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GCT(1) == 0d0) .and. &
											(ggg2Perp(Vperp1ind, Vperp2ind, Vparind, nn) /= 0d0))) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' INCONSISTENT INTEGRAND TEST FOR SPECIE= ', &
												s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', &
												Vperp1ind, ', Vperp2ind= ', &
												Vperp2ind, ', Vparind= ', Vparind, ', AND', &
												' STATISTICAL TIME-STEP= ', nn, ' IN SECOND', &
												' PERPENDICULAR 1 ION MOMENT SUBROUTINE' &
												// achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end do
								end do
							end do

							! ----------------------------------------------------

							! Compute 3D velocity space integration
							call NonUniform3DIonIntegratorSub

							! Second perpendicular 2 ion distribution function moment [J]

							if (SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) == 0d0) then
								SpecieT(s)%FluxTubeT(f)%M2Perp2phRT(nn, Qind)= 0d0
							else if (SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) /= 0d0) then
								SpecieT(s)%FluxTubeT(f)%M2Perp2phRT(nn, Qind)= &
									(SpecieT(s)%msT(1)/2d0)* &
									(MM(1)/SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind))
							end if

							! ----------------------------------------------------

							! DIAGNOSTIC FLAGS FOR MOMENT CONSISTENCY:

							if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn) == 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M2Perp2phRT(nn, Qind) /= 0d0)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' INCONSISTENT M2Perp2phRT and NqRT VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN SECOND PERPENDICULAR 2 ION MOMENT', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((sum(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%F2PerpphRTp(:, :, :, nn)) == 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M2Perp2phRT(nn, Qind) /= 0d0)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' INCONSISTENT M2Perp2phRT and F2PerpphRTp SUMMATION FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN SECOND PERPENDICULAR 2 ION MOMENT', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if (((sum(ggg2Perp(:, :, :, nn)) /= 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M2Perp2phRT(nn, Qind) == 0d0)) .or. &
								((sum(ggg2Perp(:, :, :, nn)) == 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M2Perp2phRT(nn, Qind) /= 0d0))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' INCONSISTENT M2Perp2phRT and INTEGRAND SUMMATION FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN SECOND PERPENDICULAR 2 ION MOMENT', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							! ----------------------------------------------------

							deallocate(ggg2Perp, Sum12Perp)

							deallocate(II, MM)

							! ----------------------------------------------------

						end if

						! ----------------------------------------------------

					end do
				end if
			end do

		end if

		if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 0) then

			do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
				if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
					(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacT(1:nn- 1))))) then

					! ----------------------------------------------------

					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

						! ----------------------------------------------------

						! COMPUTE ION PERPENDICULAR PLASMA ENERGY (PERPENDICULAR SECOND MOMENT):

						! ----------------------------------------------------

						if ((SpecieT(s)%FluxTubeT(f)% &
							PHASEENERGYPERPIONMOMENTflagT(1) == 1) .and. &
							(rank == 0)) then

							! ----------------------------------------------------

							allocate(ggg(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
								SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), &
								SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1), &
								Sum1(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
								SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), &
								SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))

							allocate(II(1), MM(1))

							! ----------------------------------------------------

							do Vperpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), 1
								do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1

									! ----------------------------------------------------

									if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%FphRT(nn) == 0d0) then

										ggg(Vperpind, Vparind, nn)= 0d0

									else if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%FphRT(nn) /= 0d0) then

										if (SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTCENTERflagT(1) == 0) then
											ggg(Vperpind, Vparind, nn)= &
												(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
												VCellT(Vperpind, Vparind)%VperpGCT(1)* &
												((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
												VCellT(Vperpind, Vparind)%VperpGCT(1))**2d0)* &
												SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												VCellT(Vperpind, Vparind)%FphRT(nn))
										end if
										if (SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTCENTERflagT(1) == 1) then
											ggg(Vperpind, Vparind, nn)= &
												(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
												VCellT(Vperpind, Vparind)%VperpGCT(1)* &
												((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
												VCellT(Vperpind, Vparind)%VperpGCT(1)- &
												SpecieT(s)%FluxTubeT(f)%M1PerpphRT(nn, Qind))**2d0)* &
												SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												VCellT(Vperpind, Vparind)%FphRT(nn))

											if (((isnan(real(ggg(Vperpind, Vparind, nn))) &
												.eqv. .true.)) .or. (ggg(Vperpind, Vparind, nn) == 0d0)) then
												ggg(Vperpind, Vparind, nn)= &
													(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
													VCellT(Vperpind, Vparind)%VperpGCT(1)* &
													((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
													VCellT(Vperpind, Vparind)%VperpGCT(1))**2d0)* &
													SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
													VCellT(Vperpind, Vparind)%FphRT(nn))
											end if

										end if
									end if

									! ----------------------------------------------------

									! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

									if ((isnan(real(ggg(Vperpind, Vparind, nn))) &
										.eqv. .true.) .or. (size(ggg(:, :, :)) /= &
										(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1))* &
										(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1))* &
										(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' INTEGRAND HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', &
											Vperpind, ', Vparind= ', Vparind, ', AND', &
											' STATISTICAL TIME-STEP= ', nn, ' IN SECOND', &
											' PERPENDICULAR ION MOMENT SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									! ----------------------------------------------------

									! DIAGNOSTIC FLAGS FOR CONSISTENT INTEGRANDS:

									if (((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%FphRT(nn) /= 0d0) .and. &
										(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										VCellT(Vperpind, Vparind)%VperpGCT(1) /= 0d0) .and. &
										((SpecieT(s)%FluxTubeT(f)%M1PerpphRT(nn, Qind)- &
										SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										VCellT(Vperpind, Vparind)%VperpGCT(1)) /= 0d0) .and. &
										((SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn, Qind)- &
										SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										VCellT(Vperpind, Vparind)%VparGCT(1)) /= 0d0) .and. &
										(ggg(Vperpind, Vparind, nn) == 0d0)) .or. &
										((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%FphRT(nn) == 0d0) .and. &
										(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										VCellT(Vperpind, Vparind)%VperpGCT(1) == 0d0) .and. &
										((SpecieT(s)%FluxTubeT(f)%M1PerpphRT(nn, Qind)- &
										SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										VCellT(Vperpind, Vparind)%VperpGCT(1)) == 0d0) .and. &
										((SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn, Qind)- &
										SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										VCellT(Vperpind, Vparind)%VparGCT(1)) == 0d0) .and. &
										(ggg(Vperpind, Vparind, nn) /= 0d0))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' INCONSISTENT INTEGRAND TEST FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', &
											Vperpind, ', Vparind= ', Vparind, ', AND', &
											' STATISTICAL TIME-STEP= ', nn, ' IN SECOND', &
											' PERPENDICULAR ION MOMENT SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									! ----------------------------------------------------

								end do
							end do

							! ----------------------------------------------------

							! Compute 2D velocity space integration
							call NonUniform2DIntegratorSub

							! Second perpendicular ion distribution function moment [J]

							if (SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) == 0d0) then
								SpecieT(s)%FluxTubeT(f)%M2PerpphRT(nn, Qind)= 0d0
							else if (SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) /= 0d0) then
								SpecieT(s)%FluxTubeT(f)%M2PerpphRT(nn, Qind)= &
									(SpecieT(s)%msT(1)/2d0)* &
									(MM(1)/SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind))
							end if

							! ----------------------------------------------------

							! DIAGNOSTIC FLAGS FOR MOMENT CONSISTENCY:

							if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn) == 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M2PerpphRT(nn, Qind) /= 0d0)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' INCONSISTENT M2PerpphRT and NqRT VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN SECOND PERPENDICULAR ION MOMENT', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((sum(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%FphRTp(:, :, nn)) == 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M2PerpphRT(nn, Qind) /= 0d0)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' INCONSISTENT M2PerpphRT and FphRTp SUMMATION FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN SECOND PERPENDICULAR ION MOMENT', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if (((sum(ggg(:, :, nn)) /= 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M2PerpphRT(nn, Qind) == 0d0)) .or. &
								((sum(ggg(:, :, nn)) == 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M2PerpphRT(nn, Qind) /= 0d0))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' INCONSISTENT M2PerpphRT and INTEGRAND SUMMATION FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN SECOND PERPENDICULAR ION MOMENT', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							! ----------------------------------------------------

							deallocate(ggg, Sum1)

							deallocate(II, MM)

							! ----------------------------------------------------

						end if

						! ----------------------------------------------------

					end do
				end if
			end do

		end if

		! ----------------------------------------------------

		! do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
! 			if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
! 				(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacT(1:nn- 1))))) then
! 				do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
! 					if ((SpecieT(s)%FluxTubeT(f)% &
! 						PHASEENERGYPERPIONMOMENTflagT(1) == 1) .and. &
! 						(rank == 0)) then
! 						if (SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) /= 0) then
! 							write(*, *) 'M2PerpphRT s, f, nn, Qind= ', s, f, nn, Qind, &
! 								SpecieT(s)%FluxTubeT(f)%M2PerpphRT(nn, Qind)
! 							write(*, *) 'M2Perp REFR s, f, nn, Qind= ', s, f, nn, Qind, &
! 								SpecieT(s)%FluxTubeT(f)%EnergyPerpREFRT(nn, Qind)
! 						end if
! 					end if
! 				end do
! 			end if
! 		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

		if (rank == 0) then
			if (n == 1) then
				call cpu_time(S865End)
				write(S865string, '(i10)')  nint(S865End)
				write(*, *) trim('%% 8.7- RANK= ' // adjustl(rankstring)) // &
					trim(', REAL-TIME= ' // adjustl(S865string)) // &
					trim(' s. INITIAL SECOND PERPENDICULAR ION MOMENT COMPLETE %%')
			end if
		end if

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine SecondPerpIonMomentSub

end module SecondPerpIonMoment
