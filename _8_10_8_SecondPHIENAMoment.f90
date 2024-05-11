module SecondPHIENAMoment

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.10.8 SECOND PHI-COMPONENT ENA MOMENT:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use NonUniform3DIntegrator

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! COMPUTE SECOND ENA MOMENT (PHI-ENERGY) OF PHASE-SPACE DISTRIBUTION FUNCTIONS:

	subroutine SecondPHIENAMomentSub

		! ----------------------------------------------------

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
				(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacT(1:nn- 1))))) then

				! ----------------------------------------------------

				do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

					! ----------------------------------------------------

					! COMPUTE ENA PHI-ENERGY (SECOND MOMENT):

					! ----------------------------------------------------

					if (rank == 0) then

						! ----------------------------------------------------

						allocate(gggENA(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1), &
							SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1), &
							Sum1ENA(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1), &
							SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))

						allocate(IIENA(1), MMENA(1))

						! ----------------------------------------------------

						do Vpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), 1
							do Vqind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), 1
								do Vphiind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1), 1

									! ----------------------------------------------------

									if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%FphENART(nn) == 0) then

										gggENA(Vpind, Vqind, Vphiind, nn)= 0d0

									else if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%FphENART(nn) /= 0) then

										gggENA(Vpind, Vqind, Vphiind, nn)= &
											((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
											V3CellT(Vpind, Vqind, Vphiind)%hVpCT(1))* &
											(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
											V3CellT(Vpind, Vqind, Vphiind)%hVqCT(1))* &
											(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
											V3CellT(Vpind, Vqind, Vphiind)%hVphiCT(1)))* &
											(abs((SpecieT(s)%FluxTubeT(f)%M1PHIphENART(nn, Qind)- &
											SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
											V3CellT(Vpind, Vqind, Vphiind)%VphiGCT(1))**2d0))* &
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V3CellT(Vpind, Vqind, Vphiind)%FphENART(nn)

									end if

									! ----------------------------------------------------

									! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

									if ((isnan(real(gggENA(Vpind, Vqind, Vphiind, nn))) &
										.eqv. .true.) .or. (size(gggENA(:, :, :, :)) /= &
										(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1))* &
										(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1))* &
										(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1))* &
										(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' INTEGRAND HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vpind= ', &
											Vpind, ', Vqind= ', Vqind, ', Vphiind= ', Vphiind, &
											', AND STATISTICAL TIME-STEP= ', nn, ' IN', &
											' SECOND PHI-COMPONENT ENA MOMENT SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									! ----------------------------------------------------

									! DIAGNOSTIC FLAGS FOR CONSISTENT INTEGRANDS:

									if (((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%FphENART(nn) /= 0) .and. &
										((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										V3CellT(Vpind, Vqind, Vphiind)%hVpCT(1))* &
										(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										V3CellT(Vpind, Vqind, Vphiind)%hVqCT(1))* &
										(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										V3CellT(Vpind, Vqind, Vphiind)%hVphiCT(1)) /= 0) .and. &
										((abs((SpecieT(s)%FluxTubeT(f)%M1PHIphENART(nn, Qind)- &
										SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										V3CellT(Vpind, Vqind, Vphiind)%VphiGCT(1))**2d0)) /= 0) .and. &
										(gggENA(Vpind, Vqind, Vphiind, nn) == 0)) .or. &
										((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%FphENART(nn) == 0) .and. &
										((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										V3CellT(Vpind, Vqind, Vphiind)%hVpCT(1))* &
										(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										V3CellT(Vpind, Vqind, Vphiind)%hVqCT(1))* &
										(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										V3CellT(Vpind, Vqind, Vphiind)%hVphiCT(1)) == 0) .and. &
										((abs((SpecieT(s)%FluxTubeT(f)%M1PHIphENART(nn, Qind)- &
										SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										V3CellT(Vpind, Vqind, Vphiind)%VphiGCT(1))**2d0)) == 0) .and. &
										(gggENA(Vpind, Vqind, Vphiind, nn) /= 0))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' INCONSISTENT INTEGRAND VALUE FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vpind= ', &
											Vpind, ', Vqind= ', Vqind, ', Vphiind= ', Vphiind, &
											', AND STATISTICAL TIME-STEP= ', nn, ' IN', &
											' SECOND PHI-COMPONENT ENA MOMENT SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									if ((((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										V3CellT(Vpind, Vqind, Vphiind)%hVpCT(1))* &
										(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										V3CellT(Vpind, Vqind, Vphiind)%hVqCT(1))* &
										(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										V3CellT(Vpind, Vqind, Vphiind)%hVphiCT(1)) /= 0) .and. &
										((abs((SpecieT(s)%FluxTubeT(f)%M1PHIphENART(nn, Qind)- &
										SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										V3CellT(Vpind, Vqind, Vphiind)%VphiGCT(1))**2d0)) /= 0) .and. &
										(gggENA(Vpind, Vqind, Vphiind, nn) /= 0) .and. &
										(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%g0phENART(nn) == 0)) .or. &
										(((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										V3CellT(Vpind, Vqind, Vphiind)%hVpCT(1))* &
										(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										V3CellT(Vpind, Vqind, Vphiind)%hVqCT(1))* &
										(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										V3CellT(Vpind, Vqind, Vphiind)%hVphiCT(1)) /= 0) .and. &
										((abs((SpecieT(s)%FluxTubeT(f)%M1PHIphENART(nn, Qind)- &
										SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
										V3CellT(Vpind, Vqind, Vphiind)%VphiGCT(1))**2d0)) /= 0) .and. &
										(gggENA(Vpind, Vqind, Vphiind, nn) == 0) .and. &
										(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%g0phENART(nn) /= 0))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' INCONSISTENT INTEGRAND AND g0phENART VALUE FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vpind= ', &
											Vpind, ', Vqind= ', Vqind, ', Vphiind= ', Vphiind, &
											', AND STATISTICAL TIME-STEP= ', nn, ' IN', &
											' SECOND PHI-COMPONENT ENA MOMENT SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									! ----------------------------------------------------

								end do
							end do
						end do

						! ----------------------------------------------------

						! Compute 3D velocity space integration
						call NonUniform3DIntegratorSub

						! Select distribution function moment
						if (SpecieT(s)%FluxTubeT(f)%M0phENART(nn, Qind) == 0) then
							SpecieT(s)%FluxTubeT(f)%M2PHIphENART(nn, Qind)= 0d0
						else if (SpecieT(s)%FluxTubeT(f)%M0phENART(nn, Qind) /= 0) then
							SpecieT(s)%FluxTubeT(f)%M2PHIphENART(nn, Qind)= &
								(6.242d18)*(SpecieT(s)%msT(1)/2d0)* &
								(MMENA(1)/SpecieT(s)%FluxTubeT(f)%M0phENART(nn, Qind))
						end if

						! ----------------------------------------------------

						! DIAGNOSTIC FLAGS FOR MOMENT CONSISTENCY:

						if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqENART(nn) == 0) .and. &
							(SpecieT(s)%FluxTubeT(f)%M2PHIphENART(nn, Qind) /= 0)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
								' INCONSISTENT M2PHIphENART and NqENART VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
								' TIME-STEP= ', nn, ' IN SECOND PHI-COMPONENT ENA MOMENT', &
								' SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((sum(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%FphENARTp(:, :, :, nn)) == 0) .and. &
							(SpecieT(s)%FluxTubeT(f)%M2PHIphENART(nn, Qind) /= 0)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
								' INCONSISTENT M2PHIphENART and FphENARTp SUMMATION FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
								' TIME-STEP= ', nn, ' IN SECOND PHI-COMPONENT ENA MOMENT', &
								' SUBROUTINE' // achar(27) // '[0m.'
						end if

						if (((sum(gggENA(:, :, :, nn)) /= 0) .and. &
							(SpecieT(s)%FluxTubeT(f)%M2PHIphENART(nn, Qind) == 0)) .or. &
							((sum(gggENA(:, :, :, nn)) == 0) .and. &
							(SpecieT(s)%FluxTubeT(f)%M2PHIphENART(nn, Qind) /= 0))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
								' INCONSISTENT M2PHIphENART and INTEGRAND SUMMATION FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
								' TIME-STEP= ', nn, ' IN SECOND PHI-COMPONENT ENA MOMENT', &
								' SUBROUTINE' // achar(27) // '[0m.'
						end if

						! ----------------------------------------------------

						deallocate(gggENA, Sum1ENA)

						deallocate(IIENA, MMENA)

						! ----------------------------------------------------

					end if

					! ----------------------------------------------------

				end do
			end if
		end do

		! ----------------------------------------------------

		! do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
! 			if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
! 				(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacT(1:nn- 1))))) then
! 				do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
! 					if (rank == 0) then
! 						if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqENART(nn) /= 0) then
! 							write(*, *) 'M2PHIphENART s, f, nn, Qind= ', s, f, nn, Qind, &
! 								SpecieT(s)%FluxTubeT(f)%M2PHIphENART(nn, Qind)
! 							write(*, *) 'M2PHIENA REFR s, f, nn, Qind= ', s, f, nn, Qind, &
! 								SpecieT(s)%FluxTubeT(f)%EnergyPHIENAREFRT(nn, Qind)
! 						end if
! 					end if
! 				end do
! 			end if
! 		end do

		! ----------------------------------------------------

	end subroutine SecondPHIENAMomentSub

end module SecondPHIENAMoment
