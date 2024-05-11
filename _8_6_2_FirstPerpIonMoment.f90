module FirstPerpIonMoment

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.6.2 FIRST PERPENDICULAR ION MOMENT:	%%%%%%

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

! COMPUTE FIRST PERPENDICULAR ION MOMENT (PERPENDICULAR VELOCITY) OF PHASE-SPACE DISTRIBUTION FUNCTIONS:

	subroutine FirstPerpIonMomentSub

		! ----------------------------------------------------

		if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

			do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
				if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
					(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1))))) then

					! ----------------------------------------------------

					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

						! ----------------------------------------------------

						! COMPUTE ION PLASMA PERPENDICULAR 1 VELOCITY (FIRST PERPENDICULAR 1 MOMENT):

						! ----------------------------------------------------

						if ((SpecieT(s)%FluxTubeT(f)% &
							PHASEVELPERPIONMOMENTflagT(1) == 1) .and. &
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

											ggg2Perp(Vperp1ind, Vperp2ind, Vparind, nn)= &
												SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
												V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GCT(1)* &
												SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn)

										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

										if ((isnan(real(ggg2Perp(Vperp1ind, Vperp2ind, Vparind, nn))) &
											.eqv. .true.) .or. (size(ggg2Perp(:, :, :, :)) /= &
											(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1))* &
											(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1))* &
											(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1))* &
											(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' INTEGRAND HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
												s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', &
												Vperp1ind, ', Vperp2ind= ', &
												Vperp2ind,', Vparind= ', Vparind, ', AND', &
												' STATISTICAL TIME-STEP= ', nn, ' IN FIRST', &
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
												' INCONSISTENT INTEGRAND VALUE FOR SPECIE= ', &
												s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', &
												Vperp1ind, ', Vperp2ind= ', &
												Vperp2ind, ', Vparind= ', Vparind, ', AND', &
												' STATISTICAL TIME-STEP= ', nn, ' IN FIRST', &
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

							! First perpendicular 1 ion distribution function moment [m/s]

							if (SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) == 0d0) then
								SpecieT(s)%FluxTubeT(f)%M1Perp1phRT(nn, Qind)= 0d0
							else if (SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) /= 0d0) then
								SpecieT(s)%FluxTubeT(f)%M1Perp1phRT(nn, Qind)= &
									MM(1)/SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind)
							end if

							! ----------------------------------------------------

							! DIAGNOSTIC FLAGS FOR MOMENT CONSISTENCY:

							if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn) == 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M1Perp1phRT(nn, Qind) /= 0d0)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' INCONSISTENT M1Perp1phRT and NqRT VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN FIRST PERPENDICULAR 1 ION MOMENT', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((sum(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%F2PerpphRTp(:, :, :, nn)) == 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M1Perp1phRT(nn, Qind) /= 0d0)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' INCONSISTENT M1Perp1phRT and F2PerpphRTp SUMMATION FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN FIRST PERPENDICULAR 1 ION MOMENT', &
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
					(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1))))) then

					! ----------------------------------------------------

					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

						! ----------------------------------------------------

						! COMPUTE ION PLASMA PERPENDICULAR 2 VELOCITY (FIRST PERPENDICULAR 2 MOMENT):

						! ----------------------------------------------------

						if ((SpecieT(s)%FluxTubeT(f)% &
							PHASEVELPERPIONMOMENTflagT(1) == 1) .and. &
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

											ggg2Perp(Vperp1ind, Vperp2ind, Vparind, nn)= &
												SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
												V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GCT(1)* &
												SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn)

										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

										if ((isnan(real(ggg2Perp(Vperp1ind, Vperp2ind, Vparind, nn))) &
											.eqv. .true.) .or. (size(ggg2Perp(:, :, :, :)) /= &
											(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1))* &
											(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1))* &
											(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1))* &
											(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' INTEGRAND HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
												s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', &
												Vperp1ind, ', Vperp2ind= ', &
												Vperp2ind,', Vparind= ', Vparind, ', AND', &
												' STATISTICAL TIME-STEP= ', nn, ' IN FIRST', &
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
												' INCONSISTENT INTEGRAND VALUE FOR SPECIE= ', &
												s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', &
												Vperp1ind, ', Vperp2ind= ', &
												Vperp2ind, ', Vparind= ', Vparind, ', AND', &
												' STATISTICAL TIME-STEP= ', nn, ' IN FIRST', &
												' PERPENDICULAR 2 ION MOMENT SUBROUTINE' &
												// achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end do
								end do
							end do

							! ----------------------------------------------------

							! Compute 3D velocity space integration
							call NonUniform3DIonIntegratorSub

							! First perpendicular 2 ion distribution function moment [m/s]

							if (SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) == 0d0) then
								SpecieT(s)%FluxTubeT(f)%M1Perp2phRT(nn, Qind)= 0d0
							else if (SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) /= 0d0) then
								SpecieT(s)%FluxTubeT(f)%M1Perp2phRT(nn, Qind)= &
									MM(1)/SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind)
							end if

							! ----------------------------------------------------

							! DIAGNOSTIC FLAGS FOR MOMENT CONSISTENCY:

							if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn) == 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M1Perp2phRT(nn, Qind) /= 0d0)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' INCONSISTENT M1Perp2phRT and NqRT VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN FIRST PERPENDICULAR 2 ION MOMENT', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((sum(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%F2PerpphRTp(:, :, :, nn)) == 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M1Perp2phRT(nn, Qind) /= 0d0)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' INCONSISTENT M1Perp2phRT and F2PerpphRTp SUMMATION FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN FIRST PERPENDICULAR 2 ION MOMENT', &
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
					(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1))))) then

					! ----------------------------------------------------

					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

						! ----------------------------------------------------

						! COMPUTE ION PLASMA PERPENDICULAR VELOCITY (FIRST PERPENDICULAR MOMENT):

						! ----------------------------------------------------

						if ((SpecieT(s)%FluxTubeT(f)% &
							PHASEVELPERPIONMOMENTflagT(1) == 1) .and. &
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

										ggg(Vperpind, Vparind, nn)= &
											(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
											VCellT(Vperpind, Vparind)%VperpGCGT(1)**2d0)* &
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											VCellT(Vperpind, Vparind)%FphRT(nn)

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
											' STATISTICAL TIME-STEP= ', nn, ' IN FIRST', &
											' PERPENDICULAR ION MOMENT SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									! ----------------------------------------------------

									! DIAGNOSTIC FLAGS FOR CONSISTENT INTEGRANDS:

									if (((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%FphRT(nn) /= 0d0) .and. &
										(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										VCellT(Vperpind, Vparind)%VperpGCGT(1) /= 0d0) .and. &
										(ggg(Vperpind, Vparind, nn) == 0d0)) .or. &
										((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%FphRT(nn) == 0d0) .and. &
										(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										VCellT(Vperpind, Vparind)%VperpGCGT(1) == 0d0) .and. &
										(ggg(Vperpind, Vparind, nn) /= 0d0))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' INCONSISTENT INTEGRAND VALUE FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', &
											Vperpind, ', Vparind= ', Vparind, ', AND', &
											' STATISTICAL TIME-STEP= ', nn, ' IN FIRST', &
											' PERPENDICULAR ION MOMENT SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									if (((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										VCellT(Vperpind, Vparind)%VperpGCGT(1) /= 0d0) .and. &
										(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										VCellT(Vperpind, Vparind)%VparGCGT(1) /= 0d0)) .and. &
										(((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%FphRT(nn) == 0d0) .and. &
										(ggg(Vperpind, Vparind, nn) /= 0d0)) .or. &
										((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%FphRT(nn) /= 0d0) .and. &
										(ggg(Vperpind, Vparind, nn) == 0d0)))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' INCONSISTENT INTEGRAND AND FphRT VALUE FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', &
											Vperpind, ', Vparind= ', Vparind, ', AND', &
											' STATISTICAL TIME-STEP= ', nn, ' IN FIRST', &
											' PERPENDICULAR ION MOMENT SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									if (((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										VCellT(Vperpind, Vparind)%VperpGCGT(1) /= 0d0) .and. &
										(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										VCellT(Vperpind, Vparind)%VparGCGT(1) /= 0d0) .and. &
										(ggg(Vperpind, Vparind, nn) /= 0d0) .and. &
										(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%g0phRT(nn) == 0d0)) .or. &
										((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										VCellT(Vperpind, Vparind)%VperpGCGT(1) /= 0d0) .and. &
										(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										VCellT(Vperpind, Vparind)%VparGCGT(1) /= 0d0) .and. &
										(ggg(Vperpind, Vparind, nn) == 0d0) .and. &
										(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%g0phRT(nn) /= 0d0))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' INCONSISTENT INTEGRAND AND g0phRT VALUE FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', &
											Vperpind, ', Vparind= ', Vparind, ', AND', &
											' STATISTICAL TIME-STEP= ', nn, ' IN FIRST', &
											' PERPENDICULAR ION MOMENT SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									! ----------------------------------------------------

								end do
							end do

							! ----------------------------------------------------

							! Compute 2D velocity space integration
							call NonUniform2DIntegratorSub

							! First perpendicular ion distribution function moment [m/s]

							if (SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) == 0d0) then
								SpecieT(s)%FluxTubeT(f)%M1PerpphRT(nn, Qind)= 0d0
							else if (SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) /= 0d0) then
								SpecieT(s)%FluxTubeT(f)%M1PerpphRT(nn, Qind)= &
									MM(1)/SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind)
							end if

							! ----------------------------------------------------

							! DIAGNOSTIC FLAGS FOR MOMENT CONSISTENCY:

							if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn) == 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M1PerpphRT(nn, Qind) /= 0d0)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' INCONSISTENT M1PerpphRT and NqRT VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN FIRST PERPENDICULAR ION MOMENT', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((sum(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%FphRTp(:, :, nn)) == 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M1PerpphRT(nn, Qind) /= 0d0)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' INCONSISTENT M1PerpphRT and FphRTp SUMMATION FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN FIRST PERPENDICULAR ION MOMENT', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if (((sum(ggg(:, :, nn)) /= 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M1PerpphRT(nn, Qind) == 0d0)) .or. &
								((sum(ggg(:, :, nn)) == 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M1PerpphRT(nn, Qind) /= 0d0))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' INCONSISTENT M1PerpphRT and INTEGRAND SUMMATION FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN FIRST PERPENDICULAR ION MOMENT', &
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
! 				(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1))))) then
! 				do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
! 					if ((SpecieT(s)%FluxTubeT(f)% &
! 						PHASEVELPERPIONMOMENTflagT(1) == 1) .and. &
! 						(rank == 0)) then
! 						if (SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) /= 0) then
! 							write(*, *) 'M1PerpphRT s, f, nn, Qind= ', s, f, nn, Qind, &
! 								SpecieT(s)%FluxTubeT(f)%M1PerpphRT(nn, Qind)
! 							write(*, *) 'M1Perp REFR s, f, nn, Qind= ', s, f, nn, Qind, &
! 								SpecieT(s)%FluxTubeT(f)%VperpREFRT(nn, Qind)
! 						end if
! 					end if
! 				end do
! 			end if
! 		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

		if (rank == 0) then
			if (n == 1) then
				call cpu_time(S862End)
				write(S862string, '(i10)')  nint(S862End)
				write(*, *) trim('%% 8.4- RANK= ' // adjustl(rankstring)) // &
					trim(', REAL-TIME= ' // adjustl(S862string)) // &
					trim(' s. INITIAL FIRST PERPENDICULAR ION MOMENT COMPLETE %%')
			end if
		end if

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine FirstPerpIonMomentSub

end module FirstPerpIonMoment
