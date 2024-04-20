module ZerothIonMoment

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.6.1 ZEROTH ION MOMENT:	%%%%%%

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

! COMPUTE ZEROTH MOMENT (PLASMA DENSITY) OF ION DISTRIBUTION FUNCTIONS:

	subroutine ZerothIonMomentSub

		! ----------------------------------------------------

		if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

			do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
				if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
					(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1)))) then

					! ----------------------------------------------------

					do Qind= NqLB(1), NqUB(1), 1

						! ----------------------------------------------------

						! COMPUTE ION PLASMA DENSITY (ZEROTH MOMENT):

						! ----------------------------------------------------

						if ((SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT(1) == 1) .and. &
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

											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%g02PerpphRT(nn)= 0d0

										else if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn) /= 0d0) then

											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%g02PerpphRT(nn)= &
												SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn)

										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

										if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%g02PerpphRT(nn))) &
											.eqv. .true.) .or. &
											(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%g02PerpphRT(:)) /= &
											(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' g02PerpphRT HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
												s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', &
												Vperp1ind, ', Vperp2ind= ', &
												Vperp2ind, ', Vparind= ', Vparind, ', AND', &
												' STATISTICAL TIME-STEP= ', nn, ' IN ZEROTH', &
												' ION MOMENT SUBROUTINE' &
												// achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR CONSISTENT INTEGRANDS:

										if (((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn) /= 0d0) .and. &
											(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%g02PerpphRT(nn) == 0d0)) .or. &
											((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn) == 0d0) .and. &
											(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%g02PerpphRT(nn) /= 0d0))) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' INCONSISTENT g02PerpphRT VALUE FOR SPECIE= ', &
												s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', &
												Vperp1ind, ', Vperp2ind= ', &
												Vperp2ind, ', Vparind= ', Vparind, ', AND', &
												' STATISTICAL TIME-STEP= ', nn, ' IN ZEROTH', &
												' ION MOMENT SUBROUTINE' &
												// achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

										! Select velocity space integrand
										ggg2Perp(Vperp1ind, Vperp2ind, Vparind, nn)= &
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%g02PerpphRT(nn)

										! ----------------------------------------------------

										! DIAGNOSTIC FLAG FOR PROPER INTEGRAND CONVERSION:

										if (abs(ggg2Perp(Vperp1ind, Vperp2ind, Vparind, nn)- &
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%g02PerpphRT(nn)) /= 0d0) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' ggg2Perp VALUE DOES NOT PROPERLY CONVERT FOR SPECIE= ', &
												s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
												', AND STATISTICAL TIME-STEP= ', nn, ' IN ZEROTH', &
												' ION MOMENT SUBROUTINE' &
												// achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end do
								end do
							end do

							! ----------------------------------------------------

							! Compute 3D velocity space integration
							call NonUniform3DIonIntegratorSub

							! Zeroth ion distribution function moment [m^-3]

							if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn) == 0d0) then
								SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind)= 0d0
							else if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn) /= 0d0) then
								SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind)= MM(1)
							end if

							! ----------------------------------------------------

							! DIAGNOSTIC FLAGS FOR MOMENT CONSISTENCY:

							if ((Qind == NqLB(1)) .and. (SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) == 0d0)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' ZERO DENSITY AT LOWER BOUNDARY FOR SPECIE= ', &
									SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind), &
									s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
									', AND STATISTICAL TIME-STEP= ', nn, ' IN ZEROTH', &
									' ION MOMENT SUBROUTINE' &
									// achar(27) // '[0m.'
							end if

							if (((sum(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%F2PerpphRTp(:, :, :, nn)) /= 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) == 0d0)) .or. &
								((sum(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%F2PerpphRTp(:, :, :, nn)) == 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) /= 0d0))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' INCONSISTENT M0phRT and F2PerpphRTp SUMMATION FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN ZEROTH ION MOMENT', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if (((sum(ggg2Perp(:, :, :, nn)) /= 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) == 0d0)) .or. &
								((sum(ggg2Perp(:, :, :, nn)) == 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) /= 0d0))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' INCONSISTENT M0phRT and INTEGRAND SUMMATION FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN ZEROTH ION MOMENT', &
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
					(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1)))) then

					! ----------------------------------------------------

					do Qind= NqLB(1), NqUB(1), 1

						! ----------------------------------------------------

						! COMPUTE ION PLASMA DENSITY (ZEROTH MOMENT):

						! ----------------------------------------------------

						if ((SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT(1) == 1) .and. &
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

										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%g0phRT(nn)= 0d0

									else if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%FphRT(nn) /= 0d0) then

										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											VCellT(Vperpind, Vparind)%g0phRT(nn)= &
											SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
											VCellT(Vperpind, Vparind)%VperpGCT(1)* &
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											VCellT(Vperpind, Vparind)%FphRT(nn)

									end if

									! ----------------------------------------------------

									! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

									if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%g0phRT(nn))) &
										.eqv. .true.) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%g0phRT(:)) /= &
										(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' g0phRT HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', &
											Vperpind, ', Vparind= ', Vparind, ', AND', &
											' STATISTICAL TIME-STEP= ', nn, ' IN ZEROTH', &
											' ION MOMENT SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									! ----------------------------------------------------

									! DIAGNOSTIC FLAGS FOR CONSISTENT INTEGRANDS:

									if (((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%FphRT(nn) /= 0d0) .and. &
										(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										VCellT(Vperpind, Vparind)%VperpGCT(1) /= 0d0) .and. &
										(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%g0phRT(nn) == 0d0)) .or. &
										((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%FphRT(nn) == 0d0) .and. &
										(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										VCellT(Vperpind, Vparind)%VperpGCT(1) == 0d0) .and. &
										(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%g0phRT(nn) /= 0d0))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' INCONSISTENT g0phRT VALUE FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', &
											Vperpind, ', Vparind= ', Vparind, ', AND', &
											' STATISTICAL TIME-STEP= ', nn, ' IN ZEROTH', &
											' ION MOMENT SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									if (((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										VCellT(Vperpind, Vparind)%VperpGCT(1) /= 0d0) .and. &
										(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										VCellT(Vperpind, Vparind)%VparGCT(1) /= 0d0)) .and. &
										(((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%FphRT(nn) == 0d0) .and. &
										(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%g0phRT(nn) /= 0d0)) .or. &
										((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%FphRT(nn) /= 0d0) .and. &
										(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%g0phRT(nn) == 0d0)))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' INCONSISTENT g0phRT AND FphRT VALUE FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', &
											Vperpind, ', Vparind= ', Vparind, ', AND', &
											' STATISTICAL TIME-STEP= ', nn, ' IN ZEROTH', &
											' ION MOMENT SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									! ----------------------------------------------------

									! Select velocity space integrand
									ggg(Vperpind, Vparind, nn)= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%g0phRT(nn)

									! ----------------------------------------------------

									! DIAGNOSTIC FLAG FOR PROPER INTEGRAND CONVERSION:

									if (abs(ggg(Vperpind, Vparind, nn)- &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%g0phRT(nn)) /= 0d0) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' ggg VALUE DOES NOT PROPERLY CONVERT FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
											', AND STATISTICAL TIME-STEP= ', nn, ' IN ZEROTH', &
											' ION MOMENT SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									! ----------------------------------------------------

								end do
							end do

							! ----------------------------------------------------

							! Compute 2D velocity space integration
							call NonUniform2DIntegratorSub

							! Zeroth ion distribution function moment [m^-3]

							if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn) == 0d0) then
								SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind)= 0d0

								! ----------------------------------------------------

								! DIAGNOSTIC FLAG FOR ZERO LOWER BOUNDARY DENSITY:

								if ((Qind == NqLB(1)) .and. (SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) == 0d0)) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' ZERO DENSITY AT LOWER BOUNDARY FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
										', AND STATISTICAL TIME-STEP= ', nn, ' IN ZEROTH', &
										' ION MOMENT SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								! ----------------------------------------------------

							else if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn) /= 0d0) then
								SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind)= MM(1)
							end if

							! ----------------------------------------------------

							! DIAGNOSTIC FLAGS FOR MOMENT CONSISTENCY:

							if (SpecieT(s)%FluxTubeT(f)%IONNOISEflagT(1) == 0) then
								if (((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn) == 0d0) .and. &
									(SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) /= 0d0)) .or. &
									((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn) /= 0d0) .and. &
									(SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) == 0d0))) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' INCONSISTENT M0phRT= ', SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind), &
										', and NqRT= ', SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn), &
										' VALUE FOR SPECIE= ', s, &
										', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
										' TIME-STEP= ', nn, ' IN ZEROTH ION MOMENT SUBROUTINE' &
										// achar(27) // '[0m.'
								end if
							end if

							if (((sum(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%FphRTp(:, :, nn)) /= 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) == 0d0)) .or. &
								((sum(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%FphRTp(:, :, nn)) == 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) /= 0d0))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' INCONSISTENT M0phRT and FphRTp SUMMATION FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN ZEROTH ION MOMENT', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if (((sum(ggg(:, :, nn)) /= 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) == 0d0)) .or. &
								((sum(ggg(:, :, nn)) == 0d0) .and. &
								(SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) /= 0d0))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' INCONSISTENT M0phRT and INTEGRAND SUMMATION FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN ZEROTH ION MOMENT', &
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
		! ----------------------------------------------------

		if (rank == 0) then
			if (n == 1) then
				call cpu_time(S861End)
				write(S861string, '(i10)')  nint(S861End)
				write(*, *) trim('%% 8.3- RANK= ' // adjustl(rankstring)) // &
					trim(', REAL-TIME= ' // adjustl(S861string)) // &
					trim(' s. INITIAL ZEROTH ION MOMENT COMPLETE %%')
			end if
		end if

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine ZerothIonMomentSub

end module ZerothIonMoment
