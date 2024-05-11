module ENADistributionFunctions

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.4.2 ENA DISTRIBUTION FUNCTIONS:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! COMPUTE ENA DISTRIBUTION FUNCTIONS [s^3/m^6] FOR ALL GRIDS:

	subroutine ENADistributionFunctionsSub

		! ----------------------------------------------------
		! ----------------------------------------------------

		! COMPUTE ALL ENA DISTRIBUTION FUNCTIONS:

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
				(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1))))) then

				! ----------------------------------------------------

				if (rank == 0) then
					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
						do Vpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), 1
							do Vqind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), 1
								do Vphiind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1), 1

									if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqENART(nn) == 0d0) then
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V3CellT(Vpind, Vqind, Vphiind)%FphENART(nn)= 0d0
									else if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqENART(nn) &
										/= 0d0) then
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V3CellT(Vpind, Vqind, Vphiind)%FphENART(nn)= &
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V3CellT(Vpind, Vqind, Vphiind)%NphENART(nn)/ &
											SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
											V3CellT(Vpind, Vqind, Vphiind)%d33vCT(1)
									end if

									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										FphENARTp(Vpind, Vqind, Vphiind, nn)= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%FphENART(nn)

									! ----------------------------------------------------

									! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

									if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%FphENART(nn))) &
										.eqv. .true.) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%FphENART(:)) /= &
										(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' FphENART HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vpind= ', &
											Vpind, ', Vqind= ', Vqind, ', Vphiind= ', Vphiind, &
											', AND STATISTICAL TIME-STEP= ', nn, &
											' IN ENA DISTRIBUTION', &
											' FUNCTIONS SUBROUTINE' // achar(27) // '[0m.'
									end if

									! ----------------------------------------------------

									! DIAGNOSTIC FLAGS FOR NON-ZERO AND ZERO NphENART AND
									! FphENART ELEMENTS:

									if (((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%FphENART(nn) /= 0) .and. &
										(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%NphENART(nn) == 0)) &
										.or. ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%FphENART(nn) == 0) .and. &
										(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%NphENART(nn) /= 0))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' ZERO AND NON-ZERO NphENART AND FphENART ELEMENTS', &
											' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
											', Qind= ', Qind, ', Vpind= ', &
											Vpind, ', Vqind= ', Vqind, ', Vphiind= ', Vphiind, &
											', AND STATISTICAL TIME-STEP= ', nn, &
											' IN ENA DISTRIBUTION', &
											' FUNCTIONS SUBROUTINE' // achar(27) // '[0m.'
									end if

									if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%FphENART(nn) < 0) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' NEGATIVE FphENART ELEMENT FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vpind= ', &
											Vpind, ', Vqind= ', Vqind, ', Vphiind= ', Vphiind, &
											', AND STATISTICAL TIME-STEP= ', nn, &
											' IN ENA DISTRIBUTION', &
											' FUNCTIONS SUBROUTINE' // achar(27) // '[0m.'
									end if

									! ----------------------------------------------------

								end do
							end do
						end do
					end do
				end if

				! ----------------------------------------------------

			end if
		end do

		! ----------------------------------------------------

	end subroutine ENADistributionFunctionsSub

end module ENADistributionFunctions
