module IonDistributionFunctions

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.4.1 ION DISTRIBUTION FUNCTIONS:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! COMPUTE ION DISTRIBUTION FUNCTIONS [s^3/m^6] FOR ALL GRIDS:

	subroutine IonDistributionFunctionsSub

		! ----------------------------------------------------
		! ----------------------------------------------------

		! COMPUTE ALL ION DISTRIBUTION FUNCTIONS:

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
				(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1)))) then

				! ----------------------------------------------------

				if (rank == 0) then
					do Qind= NqLB(1), NqUB(1), 1
						if (SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1) == 1) then

							if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

								do Vperp1ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), 1
									do Vperp2ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), 1
										do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1

											! ----------------------------------------------------

											if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn) == 0d0) then
												SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
													V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn)= 0d0
											else if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn) /= 0d0) then
												SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
													V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn)= &
													SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
													V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%N2PerpphRT(nn)/ &
													SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
													V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%d3vCT(1)
											end if

											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												F2PerpphRTp(Vperp1ind, Vperp2ind, Vparind, nn)= &
												SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn)

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

											if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn))) &
												.eqv. .true.) .or. &
												(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(:)) /= &
												(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' F2PerpphRT HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
													s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', &
													Vperp1ind, ', Vperp2ind= ', Vperp2ind, ', Vparind= ', Vparind, ', AND', &
													' STATISTICAL TIME-STEP= ', nn, ' IN ION DISTRIBUTION', &
													' FUNCTIONS SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO AND ZERO N2PerpphRT AND F2PerpphRT ELEMENTS:

											if (((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn) /= 0) .and. &
												(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%N2PerpphRT(nn) == 0)) &
												.or. ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn) == 0) .and. &
												(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%N2PerpphRT(nn) /= 0))) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' ZERO AND NON-ZERO N2PerpphRT AND F2PerpphRT ELEMENTS FOR SPECIE= ', &
													s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', &
													Vperp1ind, ', Vperp2ind= ', Vperp2ind, ', Vparind= ', Vparind, ', AND', &
													' STATISTICAL TIME-STEP= ', nn, ' IN ION DISTRIBUTION', &
													' FUNCTIONS SUBROUTINE' // achar(27) // '[0m.'
											end if

											if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT(nn) < 0) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' NEGATIVE F2PerpphRT ELEMENT FOR SPECIE= ', s, &
													', FLUX TUBE= ', f, ', Qind= ', Qind, &
													', Vperp1ind= ', Vperp1ind, ', Vperp2ind= ', Vperp2ind, ', Vparind= ', &
													Vparind, ', AND STATISTICAL TIME-STEP= ', nn, &
													' IN ION DISTRIBUTION FUNCTIONS SUBROUTINE' &
													// achar(27) // '[0m.'
											end if

											if (abs(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%N2PerpphRT(nn)- &
												SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												F2PerpphRTp(Vperp1ind, Vperp2ind, Vparind, nn)*SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
												V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%d3vCT(1)) >= 1d0) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' F2PerpphRT= ', &
													SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
													F2PerpphRTp(Vperp1ind, Vperp2ind, Vparind, nn)* &
													SpecieT(s)%FluxTubeT(f)%QCellT(1)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%d3vCT(1), &
													' HAS INCORRECT INVERSE WITH N2PerpphRT= ', &
													SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
													V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%N2PerpphRT(nn), ' FOR SPECIE= ', &
													s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', &
													Vperp1ind, ', Vperp2ind= ', Vperp2ind, ', Vparind= ', Vparind, ', AND', &
													' STATISTICAL TIME-STEP= ', nn, ' IN ION DISTRIBUTION', &
													' FUNCTIONS SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end do
									end do
								end do

							else

								do Vperpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), 1
									do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1

										! ----------------------------------------------------

										if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn) == 0d0) then
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												VCellT(Vperpind, Vparind)%FphRT(nn)= 0d0
										else if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn) /= 0d0) then
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												VCellT(Vperpind, Vparind)%FphRT(nn)= &
												SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												VCellT(Vperpind, Vparind)%NphRT(nn)/ &
												SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
												VCellT(Vperpind, Vparind)%d3vCT(1)
										end if

										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											FphRTp(Vperpind, Vparind, nn)= &
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											VCellT(Vperpind, Vparind)%FphRT(nn)

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

										if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											VCellT(Vperpind, Vparind)%FphRT(nn))) &
											.eqv. .true.) .or. &
											(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											VCellT(Vperpind, Vparind)%FphRT(:)) /= &
											(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' FphRT HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
												s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', &
												Vperpind, ', Vparind= ', Vparind, ', AND', &
												' STATISTICAL TIME-STEP= ', nn, ' IN ION DISTRIBUTION', &
												' FUNCTIONS SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO AND ZERO NphRT AND FphRT ELEMENTS:

										if (((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											VCellT(Vperpind, Vparind)%FphRT(nn) /= 0) .and. &
											(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											VCellT(Vperpind, Vparind)%NphRT(nn) == 0)) &
											.or. ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											VCellT(Vperpind, Vparind)%FphRT(nn) == 0) .and. &
											(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											VCellT(Vperpind, Vparind)%NphRT(nn) /= 0))) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' ZERO AND NON-ZERO NphRT AND FphRT ELEMENTS FOR SPECIE= ', &
												s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', &
												Vperpind, ', Vparind= ', Vparind, ', AND', &
												' STATISTICAL TIME-STEP= ', nn, ' IN ION DISTRIBUTION', &
												' FUNCTIONS SUBROUTINE' // achar(27) // '[0m.'
										end if

										if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											VCellT(Vperpind, Vparind)%FphRT(nn) < 0) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NEGATIVE FphRT ELEMENT FOR SPECIE= ', s, &
												', FLUX TUBE= ', f, ', Qind= ', Qind, &
												', Vperpind= ', Vperpind, ', Vparind= ', &
												Vparind, ', AND STATISTICAL TIME-STEP= ', nn, &
												' IN ION DISTRIBUTION FUNCTIONS SUBROUTINE' &
												// achar(27) // '[0m.'
										end if

										if (abs(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											VCellT(Vperpind, Vparind)%NphRT(nn)- &
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											FphRTp(Vperpind, Vparind, nn)*SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
											VCellT(Vperpind, Vparind)%d3vCT(1)) >= 1d0) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' FphRT= ', &
												SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												FphRTp(Vperpind, Vparind, nn)* &
												SpecieT(s)%FluxTubeT(f)%QCellT(1)%VCellT(Vperpind, Vparind)%d3vCT(1), &
												' HAS INCORRECT INVERSE WITH NphRT= ', &
												SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
													VCellT(Vperpind, Vparind)%NphRT(nn), ' FOR SPECIE= ', &
												s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', &
												Vperpind, ', Vparind= ', Vparind, ', AND', &
												' STATISTICAL TIME-STEP= ', nn, ' IN ION DISTRIBUTION', &
												' FUNCTIONS SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end do
								end do

							end if

						end if
					end do
				end if

				! ----------------------------------------------------

			end if
		end do

		! ----------------------------------------------------

	end subroutine IonDistributionFunctionsSub

end module IonDistributionFunctions
