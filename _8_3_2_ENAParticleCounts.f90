module ENAParticleCounts

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.3.2 ENA PARTICLE COUNTS:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use MPIReduceSum

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! COMPUTE STATISTICAL GRID COUNTS FOR ALL GRIDS:

	subroutine ENAParticleCountsSub

		! ----------------------------------------------------
		! ----------------------------------------------------

		! COMPUTE STATISTICAL GRID COUNTS FOR ENA CONFIG-SPACE GRIDS:

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
				(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn)))) then

				! ----------------------------------------------------

				! COMPUTE LOGICAL FLAGS FOR CONFIG-SPACE STATISTICAL BINNING:

				if (SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1) == 1) then

					! ----------------------------------------------------

					! COMPUTE RAW ENA CONFIG-SPACE COUNTS:

					if ((n == 1) .and. (nn == 1)) then
						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
							SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind)= 0d0
						end do
					end if

					if (((n /= 1) .and. (nn /= 1) .and. &
						(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn)))) then
						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
							jcount= 0d0
							SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind)= jcount
							jloopec: do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
								if (ENAflag(j) .eqv. .true.) then
									if (Qind == Qindk1(j)) then

										jcount= jcount+ 1d0
										SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind)= jcount

										cycle jloopec

									end if
								end if
							end do jloopec
						end do
					end if

					! ----------------------------------------------------

					! DIAGNOSTIC FLAGS FOR TOTAL PARTICLE NUMBER CONSERVATION IN CONFIG-SPACE:

					! ----------------------------------------------------

					if ((nn /= 1) .and. (n /= 1)) then

						! ----------------------------------------------------

						!if ((sum(SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, :))+ &
						!	sum(SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, :)) /= &
						!	SpecieT(s)%FluxTubeT(f)%NsnT(nn)) .and. (dNstK1(1) == 0)) then
						!	write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' TOTAL PARTICLE', &
						!		' NUMBER= ', sum(SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, :))+ &
						!		sum(SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, :)), &
						!		', NsT= ', SpecieT(s)%FluxTubeT(f)%NsnT(nn), &
						!		' NOT CONSERVED IN CONFIG-SPACE GRID FOR SPECIE= ', s, &
						!		', FLUX TUBE= ', f, ', AND STATISTICAL TIME-STEP= ', nn, &
						!		' IN CONFIG PARTICLE COUNTS SUBROUTINE' // achar(27) // '[0m.'
						!end if

						! ----------------------------------------------------

					end if

					! ----------------------------------------------------

					! REDUCE ALL STATISTICAL PARTICLE COUNTS IN CONFIG-SPACE TO MPI ROOT RANK (0):

					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

						! ----------------------------------------------------

						if ((nn /= 1) .and. (n /= 1)) then

							! ----------------------------------------------------

							SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqENAT(nn)= &
								SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind)* &
								(SpecieT(s)%FluxTubeT(f)%nsnormfacT(1)/ &
								SpecieT(s)%FluxTubeT(f)%d3xCT(nn, Qind))
							SpecieT(s)%FluxTubeT(f)%NqENATp(nn, Qind)= &
								SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqENAT(nn)

							call mpi_barrier(MPI_COMM_WORLD, ierr)
							NqENA(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqENAT(nn)
							MPIRedSumIn(1)= NqENA(1)
							call MPIReduceSumSub
							NqENAR(1)= MPIRedSumOut(1)

							!call mpi_reduce(NqENA(1), NqENAR(1), 1, &
							!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

							if (rank == 0) then
								SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqENART(nn)= NqENAR(1)

								SpecieT(s)%FluxTubeT(f)%NqENARTp(nn, Qind)= &
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqENART(nn)
							end if

							! ----------------------------------------------------

						end if

						! ----------------------------------------------------

					end do

					! ----------------------------------------------------

					! DIAGNOSTIC FLAG FOR NaN VALUES OF NqRT AND CONSISTENCY:

          ! ----------------------------------------------------

					if (rank == 0) then
						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

							! ----------------------------------------------------

							if ((nn /= 1) .and. (n /= 1)) then

								! ----------------------------------------------------

								if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqENART(nn))) &
									.eqv. .true.) .or. &
									(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqENART(:)) /= &
									(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqENART HAS', &
										' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
										', Qind= ', Qind, ', AND STATISTICAL TIME-STEP= ', nn, &
										' IN CONFIG PARTICLE COUNTS SUBROUTINE' // achar(27) // '[0m.'
								end if

								if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqENAT(nn) /= 0 &
									.and. SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqENART(nn) == 0)) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' ZERO AND NON-ZERO NqENAT AND NqENART ELEMENTS DO NOT MATCH', &
										' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
										', AND STATISTICAL TIME-STEP= ', nn, ' IN CONFIG PARTICLE COUNTS', &
										' SUBROUTINE' // achar(27) // '[0m.'
								end if

								if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqENART(nn) < 0) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' NEGATIVE NqENART ELEMENTS FOR SPECIE= ', s, &
										', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
										' TIME-STEP= ', nn, ' IN CONFIG PARTICLE COUNTS SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								! ----------------------------------------------------

							end if

							! ----------------------------------------------------

						end do
					end if

					! ----------------------------------------------------

				end if

				! ----------------------------------------------------

			end if
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

		! COMPUTE STATISTICAL GRID COUNTS FOR ENA VEL-SPACE GRIDS:

		! ----------------------------------------------------
		! ----------------------------------------------------

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
				(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn)))) then

				! ----------------------------------------------------

				! COMPUTE LOGICAL FLAGS FOR ENA PHASE-SPACE STATISTICAL BINNING:

				do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

					! ----------------------------------------------------

					if (SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind) == 0d0) then
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
							NphReNormENAT(nn, :, :, :)= 0d0
					else if (SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind) /= 0d0) then

						! ----------------------------------------------------

						! COMPUTE RAW ENA PHASE-SPACE COUNTS:

						if ((n /= 1) .and. (nn /= 1) .and. &
							(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn))) then
							do Vpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), 1
								do Vqind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), 1
									do Vphiind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1), 1
										jcount= 0d0
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												NphReNormENAT(nn, Vpind, Vqind, Vphiind)= jcount
										jloopev: do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
											if (ENAflag(j) .eqv. .true.) then
												if ((Vpind == Vpindk1(j)) .and. (Vqind == Vqindk1(j)) &
													.and. (Vphiind == Vphiindk1(j)) .and. (Qind == Qindk1(j))) then

														jcount= jcount+ 1d0
														SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
																NphReNormENAT(nn, Vpind, Vqind, Vphiind)= jcount

														cycle jloopev

												end if
											end if
										end do jloopev
									end do
								end do
							end do
						end if

						! ----------------------------------------------------

					end if

					! ----------------------------------------------------

				end do

				! ----------------------------------------------------

				do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

					! ----------------------------------------------------

					! DIAGNOSTIC FLAGS FOR CONSISTENCY:

					!if ((sum(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
					!	NphReNormENAT(nn, :, :, :))- SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind)) > 1d0) then
					!	write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					!		' UNEQUAL NqReNormENAT= ', SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind), &
					!		' AND NphReNormENAT SUM= ', sum(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
					!			NphReNormENAT(nn, :, :, :)), &
					!		' FOR SPECIE= ', s, &
					!		', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
					!		' TIME-STEP= ', nn, ' IN ENA VEL PARTICLE COUNTS SUBROUTINE' &
					!		// achar(27) // '[0m.'
					!end if

					! ----------------------------------------------------

					! REDUCE ALL RAW ENA COUNTS IN PHASE-SPACE TO MPI ROOT RANK (0):

					do Vpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), 1
						do Vqind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), 1
							do Vphiind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1), 1

								! ----------------------------------------------------

								call mpi_barrier(MPI_COMM_WORLD, ierr)
								NphENA(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									NphReNormENAT(nn, Vpind, Vqind, Vphiind)
								MPIRedSumIn(1)= NphENA(1)
								call MPIReduceSumSub
								NphENAR(1)= MPIRedSumOut(1)

								!call mpi_reduce(NphENA(1), NphENAR(1), 1, &
								!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

								if (rank == 0) then
									if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqENART(nn) /= 0) then

									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										NphReNormENART(nn, Vpind, Vqind, Vphiind)= &
											NphENAR(1)

									else if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqENART(nn) == 0) then

										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											NphReNormENART(nn, Vpind, Vqind, Vphiind)= 0d0

									end if

									! ----------------------------------------------------

									! FILTER STATISTICAL NOISE FOR UN-NORMALIZED ENA COUNTS:

									if (SpecieT(1)%FluxTubeT(1)%ENANOISEflagT(1) == 1) then
										if (((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											NphReNormENART(nn, Vpind, Vqind, Vphiind))**(-1d0/2d0)) > &
											SpecieT(s)%FluxTubeT(f)%ENANoiseLimitT(1)) then
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												NphReNormENART(nn, Vpind, Vqind, Vphiind)= 0d0
										end if
									end if

									! ----------------------------------------------------

									if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										NphReNormENART(nn, Vpind, Vqind, Vphiind) == 0d0) then
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V3CellT(Vpind, Vqind, Vphiind)%NphENART(nn)=	0d0
									else if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										NphReNormENART(nn, Vpind, Vqind, Vphiind) /= 0d0) then
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V3CellT(Vpind, Vqind, Vphiind)%NphENART(nn)= &
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											NphReNormENART(nn, Vpind, Vqind, Vphiind)* &
											(SpecieT(s)%FluxTubeT(f)%nsnormfacT(1)/ &
											SpecieT(s)%FluxTubeT(f)%d3xCT(nn, Qind))
									end if

									SpecieT(s)%FluxTubeT(f)% &
										QCellT(Qind)%NphENARTp(Vpind, Vqind, Vphiind, nn)= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%NphENART(nn)

									! ----------------------------------------------------

									! DIAGNOSTIC FLAGS FOR CONSISTENCY:

									if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										NphReNormENART(nn, Vpind, Vqind, Vphiind) /= 0 .and. &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%NphENART(nn) == 0) .or. &
										(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										NphReNormENART(nn, Vpind, Vqind, Vphiind) == 0 .and. &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%NphENART(nn) /= 0)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' ZERO AND NON-ZERO NphENART AND NphReNormENART', &
											' ELEMENTS DO NOT MATCH FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vpind= ', &
											Vpind, ', Vqind= ', Vqind, ', AND Vphiind= ', &
											Vphiind, ', AND STATISTICAL TIME-STEP= ', nn, &
											' IN ENA VEL PARTICLE COUNTS SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%NphENART(nn) < 0) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' NEGATIVE NphENART ELEMENTS FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vpind= ', &
											Vpind, ', Vqind= ', Vqind, ', AND Vphiind= ', &
											Vphiind, ', AND STATISTICAL TIME-STEP= ', nn, &
											' IN ENA VEL PARTICLE COUNTS SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										NphReNormENART(nn, Vpind, Vqind, Vphiind))) .eqv. .true.) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										NphReNormENART(:, :, :, :)) /= &
										(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)* &
										SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1)* &
										SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1)* &
										SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' NphReNormENART HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vpind= ', &
											Vpind, ', Vqind= ', Vqind, ', AND Vphiind= ', &
											Vphiind, ', AND STATISTICAL TIME-STEP= ', nn, &
											' IN ENA VEL PARTICLE COUNTS SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%NphENART(nn))) .eqv. .true.) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%NphENART(:)) /= &
										(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NphENART', &
											' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vpind= ', &
											Vpind, ', Vqind= ', Vqind, ', AND Vphiind= ', &
											Vphiind, ', AND STATISTICAL TIME-STEP= ', nn, &
											' IN ENA VEL PARTICLE COUNTS SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									! ----------------------------------------------------

								end if

								! ----------------------------------------------------

							end do
						end do
					end do

					! ----------------------------------------------------

				end do

				! ----------------------------------------------------

			end if
		end do

		if (rank == 0) then
			do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
				if (((n /= 1) .and. (nn /= 1) .and. &
					(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn)))) then
					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

						if (sum(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
							NphReNormENART(nn, :, :, :)) /= 0d0) then

							write(*, *) 'NphReNormENART= ', nn, Qind, sum(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								NphReNormENART(nn, :, :, :))

						end if

					end do
				end if
			end do
		end if

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine ENAParticleCountsSub

end module ENAParticleCounts
