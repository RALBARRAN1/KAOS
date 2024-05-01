module IonParticleCounts

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.3.1 ION PARTICLE COUNTS:	%%%%%%

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

	subroutine IonParticleCountsSub

		! ----------------------------------------------------
		! ----------------------------------------------------

		! COMPUTE STATISTICAL GRID COUNTS FOR ION CONFIG-SPACE GRIDS:

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
				(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1)))) then

				! ----------------------------------------------------

				! COMPUTE LOGICAL FLAGS FOR CONFIG-SPACE STATISTICAL BINNING:

				! ----------------------------------------------------

				! COMPUTE RAW ION CONFIG-SPACE COUNTS:

				if ((n == 1) .and. (nn == 1)) then
					do Qind= NqLB(1), NqUB(1), 1
						jcount= 0d0
						SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind)= jcount
						jloopic1: do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
							if (ENAflag(j) .eqv. .false.) then
								if (Qind == Qindk0(j)) then

									jcount= jcount+ 1d0
									SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind)= jcount

									cycle jloopic1

								end if
							end if
						end do jloopic1
					end do
				end if
				if ((n /= 1) .and. (nn /= 1) .and. &
					(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1))) then
					do Qind= NqLB(1), NqUB(1), 1
						jcount= 0d0
						SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind)= jcount
						jloopic2: do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
							if (ENAflag(j) .eqv. .false.) then
								if (Qind == Qindk1(j)) then

									jcount= jcount+ 1d0
									SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind)= jcount

									cycle jloopic2

								end if
							end if
						end do jloopic2
					end do
				end if

				! ----------------------------------------------------

				! DIAGNOSTIC FLAGS FOR TOTAL PARTICLE NUMBER CONSERVATION IN CONFIG-SPACE:

				! ----------------------------------------------------

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 0) then

					! ----------------------------------------------------

					!if ((sum(SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, :)) /= &
					!	SpecieT(s)%FluxTubeT(f)%NsnT(nn)) .and. (dNsTK1(1) == 0)) then
					!	write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' TOTAL PARTICLE', &
					!		' NUMBER= ', sum(SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, :)), &
					!		', NsnT= ', SpecieT(s)%FluxTubeT(f)%NsnT(nn), &
					!		', NsT= ', SpecieT(s)%FluxTubeT(f)%NsT(1), &
					!		', dNsTK1= ', dNsTK1(1), &
					!		', dNsTK2= ', dNsTK2(1), ', dNsTK3= ', dNsTK3(1), &
					!		' NOT CONSERVED IN CONFIG-SPACE GRID FOR SPECIE= ', s, &
					!		', FLUX TUBE= ', f, ', AND STATISTICAL TIME-STEP= ', nn, &
					!		' IN CONFIG PARTICLE COUNTS SUBROUTINE' // achar(27) // '[0m.'
					!end if

					! ----------------------------------------------------

				end if

				! ----------------------------------------------------

				! REDUCE ALL STATISTICAL PARTICLE COUNTS IN CONFIG-SPACE TO MPI ROOT RANK (0):

				do Qind= NqLB(1), NqUB(1), 1

					! ----------------------------------------------------

					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqT(nn)= &
						SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind)* &
						(SpecieT(s)%FluxTubeT(f)%nsnormfacT(1)/ &
						SpecieT(s)%FluxTubeT(f)%d3xCT(nn, Qind))
					SpecieT(s)%FluxTubeT(f)%NqTp(nn, Qind)= &
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqT(nn)

					call mpi_barrier(MPI_COMM_WORLD, ierr)
					Nq(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqT(nn)
					MPIRedSumIn(1)= Nq(1)
					call MPIReduceSumSub
					NqR(1)= MPIRedSumOut(1)

					!call mpi_reduce(Nq(1), NqR(1), 1, &
					!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

					if (rank == 0) then
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn)= NqR(1)
						SpecieT(s)%FluxTubeT(f)%NqRTp(nn, Qind)= &
							SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn)
					end if

					! ----------------------------------------------------

				end do

				! ----------------------------------------------------

				! DIAGNOSTIC FLAG FOR NaN VALUES OF NqRT AND CONSISTENCY:

				if (rank == 0) then
					do Qind= NqLB(1), NqUB(1), 1
						if (Qind == NqLB(1)) then
							if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn) == 0d0) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' ZERO LOWER BOUNDARY ION COUNT NqRT= ', &
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn), &
									' SPECIE= ', s, ', FLUX TUBE= ', f, &
									', Qind= ', Qind, ', AND STATISTICAL TIME-STEP= ', nn, &
									' IN ION PARTICLE COUNTS SUBROUTINE' // achar(27) // '[0m.'
							end if
						end if
						if (nn == 1) then
							!if (Qind == NqLB(1)) then
							!	if (abs(nint(SpecieT(s)%FluxTubeT(f)%NqRTp(nn, Qind)*(SpecieT(s)%FluxTubeT(f)%d3xCT(nn, Qind)/ &
							!		SpecieT(s)%FluxTubeT(f)%nsnormfacT(1)))- &
							!		nint(SpecieT(s)%FluxTubeT(f)%LBNominalDensityT(1)*(SpecieT(s)%FluxTubeT(f)%d3xCLBT(1)/ &
							!		SpecieT(s)%FluxTubeT(f)%nsnormfacT(1)))) > 1d0) then
							!		write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							!			' CONFIG-SPACE GRID COUNTS= ', &
							!			nint(SpecieT(s)%FluxTubeT(f)%NqRTp(nn, Qind)* &
							!			(SpecieT(s)%FluxTubeT(f)%d3xCT(nn, Qind)/ &
							!			SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))), &
							!			' ARE UNEQUAL TO LB NOMINAL DENSITY= ', &
							!			nint(SpecieT(s)%FluxTubeT(f)%LBNominalDensityT(1)* &
							!			(SpecieT(s)%FluxTubeT(f)%d3xCLBT(1)/ &
							!			SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))), &
							!			' FOR SPECIE= ', s, ', AND FLUX TUBE= ', f, &
							!			' IN ION PARTICLE COUNTS SUBROUTINE'	// achar(27) // '[0m.'
							!	end if
							!end if
						!	if (Qind == NqUB(1)) then
						!		if (abs(nint(SpecieT(s)%FluxTubeT(f)%NqRTp(nn, Qind)*(SpecieT(s)%FluxTubeT(f)%d3xCT(nn, Qind)/ &
						!		SpecieT(s)%FluxTubeT(f)%nsnormfacT(1)))- &
						!		nint(SpecieT(s)%FluxTubeT(f)%UBNominalDensityT(1)*(SpecieT(s)%FluxTubeT(f)%d3xCUBT(1)/ &
						!		SpecieT(s)%FluxTubeT(f)%nsnormfacT(1)))) > 1d0) then
						!			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						!				' CONFIG-SPACE GRID COUNTS= ', &
						!				nint(SpecieT(s)%FluxTubeT(f)%NqRTp(nn, Qind)* &
						!				(SpecieT(s)%FluxTubeT(f)%d3xCT(nn, Qind)/ &
						!				SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))), &
						!				' ARE UNEQUAL TO UB NOMINAL DENSITY= ', &
						!				nint(SpecieT(s)%FluxTubeT(f)%UBNominalDensityT(1)* &
						!				(SpecieT(s)%FluxTubeT(f)%d3xCUBT(1)/ &
						!				SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))), &
						!				' FOR SPECIE= ', s, ', AND FLUX TUBE= ', f, &
						!				' IN ION PARTICLE COUNTS SUBROUTINE'	// achar(27) // '[0m.'
						!		end if
						!	end if
							if (SpecieT(s)%FluxTubeT(f)%qGLT(nn, 1) <= 0) then
								if (SpecieT(s)%FluxTubeT(f)%NqICAT(1) == Qind) then
									if (nint(SpecieT(s)%FluxTubeT(f)%NqRTp(nn, Qind)* &
										(SpecieT(s)%FluxTubeT(f)%d3xCT(nn, Qind)/ &
										SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))) /= &
										SpecieT(s)%FluxTubeT(f)%QCellICT(Qind- &
										SpecieT(s)%FluxTubeT(f)%NqICAT(1)+ 1)%NsFARRT(1)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' TEST CONFIG-SPACE GRID COUNTS= ', &
											nint(SpecieT(s)%FluxTubeT(f)%NqRTp(nn, Qind)* &
											(SpecieT(s)%FluxTubeT(f)%d3xCT(nn, Qind)/ &
											SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))), &
											' AT INITIAL TIME ARE UNEQUAL', &
											' TO INITIAL Q CELL POPULATION= ', &
											SpecieT(s)%FluxTubeT(f)%QCellICT(Qind- &
											SpecieT(s)%FluxTubeT(f)%NqICAT(1)+ 1)%NsFARRT(1), &
											' FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
											' IN ION PARTICLE COUNTS SUBROUTINE'	// achar(27) // '[0m.'
									end if
								end if
							end if
							if (SpecieT(s)%FluxTubeT(f)%qGLT(nn, 1) > 0) then
								if (SpecieT(s)%FluxTubeT(f)%NqICAT(1) == Qind) then
									if (nint(SpecieT(s)%FluxTubeT(f)%NqRTp(nn, Qind)* &
										(SpecieT(s)%FluxTubeT(f)%d3xCT(nn, Qind)/ &
										SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))) /= &
										SpecieT(s)%FluxTubeT(f)%QCellICT(Qind- &
										SpecieT(s)%FluxTubeT(f)%NqICAT(1)+ 1)%NsFARRT(1)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' TEST CONFIG-SPACE GRID COUNTS= ', &
											nint(SpecieT(s)%FluxTubeT(f)%NqRTp(nn, Qind)* &
											(SpecieT(s)%FluxTubeT(f)%d3xCT(nn, Qind)/ &
											SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))), &
											' AT INITIAL TIME ARE UNEQUAL', &
											' TO INITIAL Q CELL POPULATION= ', &
											SpecieT(s)%FluxTubeT(f)%QCellICT(Qind- &
											SpecieT(s)%FluxTubeT(f)%NqICAT(1)+ 1)%NsFARRT(1), &
											' FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
											' IN ION PARTICLE COUNTS SUBROUTINE'	// achar(27) // '[0m.'
									end if
								end if
							end if
						end if

						if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn))) &
							.eqv. .true.) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(:)) /= &
							(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqRT HAS', &
								' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', Qind= ', Qind, ', AND STATISTICAL TIME-STEP= ', nn, &
								' IN ION PARTICLE COUNTS SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqT(nn) /= 0 &
							.and. SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn) == 0)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
								' ZERO AND NON-ZERO NqT AND NqRT ELEMENTS DO NOT MATCH', &
								' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
								', AND STATISTICAL TIME-STEP= ', nn, ' IN ION PARTICLE COUNTS', &
								' SUBROUTINE' // achar(27) // '[0m.'
						end if

						if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn) < 0) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
								' NEGATIVE NqRT ELEMENTS FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
								' TIME-STEP= ', nn, ' IN ION PARTICLE COUNTS SUBROUTINE' &
								// achar(27) // '[0m.'
						end if

						! ----------------------------------------------------

					end do
				end if

				! ----------------------------------------------------

			end if
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

		! COMPUTE STATISTICAL GRID COUNTS FOR 3D AND 2D ION VEL-SPACE GRIDS:

		! ----------------------------------------------------
		! ----------------------------------------------------

		if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

			do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
				if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
					(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1)))) then

					! ----------------------------------------------------

					! COMPUTE LOGICAL FLAGS FOR PHASE-SPACE STATISTICAL BINNING:

					do Qind= NqLB(1), NqUB(1), 1

						! ----------------------------------------------------

						if (SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind) == 0d0) then
							SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								N2PerpphReNormT(nn, :, :, :)= 0d0
						else if (SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind) /= 0d0) then

							! ----------------------------------------------------

							! COMPUTE RAW ION PHASE-SPACE COUNTS:

							if ((n == 1) .and. (nn == 1)) then
								do Vperp1ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), 1
									do Vperp2ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), 1
										do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1
											jcount= 0d0
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
													N2PerpphReNormT(nn, Vperp1ind, Vperp2ind, Vparind)= jcount
											jloopiv1: do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
												if (ENAflag(j) .eqv. .false.) then
													if ((Vperp1ind == Vperp1indk0(j)) .and. (Vperp2ind == Vperp2indk0(j)) &
													 	.and. (Vparind == Vparindk0(j)) .and. (Qind == Qindk0(j))) then

														jcount= jcount+ 1d0
														SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
																N2PerpphReNormT(nn, Vperp1ind, Vperp2ind, Vparind)= jcount

														cycle jloopiv1

													end if
												end if
											end do jloopiv1
										end do
									end do
								end do
							end if
							if ((n /= 1) .and. (nn /= 1) .and. &
								(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1))) then
								do Vperp1ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), 1
									do Vperp2ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), 1
										do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1
											jcount= 0d0
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
													N2PerpphReNormT(nn, Vperp1ind, Vperp2ind, Vparind)= jcount
											jloopiv2: do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
												if (ENAflag(j) .eqv. .false.) then
													if ((Vperp1ind == Vperp1indk1(j)) .and. (Vperp2ind == Vperp2indk1(j)) &
												 		.and. (Vparind == Vparindk1(j)) .and. (Qind == Qindk1(j))) then

														jcount= jcount+ 1d0
														SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
																N2PerpphReNormT(nn, Vperp1ind, Vperp2ind, Vparind)= jcount

														cycle jloopiv2

													end if
												end if
											end do jloopiv2
										end do
									end do
								end do
							end if

							! ----------------------------------------------------

						end if

						! ----------------------------------------------------

					end do

					! ----------------------------------------------------

					do Qind= NqLB(1), NqUB(1), 1

						! ----------------------------------------------------

						! DIAGNOSTIC FLAGS FOR CONSISTENCY:

						!if ((sum(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
						!	N2PerpphReNormT(nn, :, :, :))- &
						!	SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind)) > 1d0) then
						!	write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						!		' UNEQUAL NqReNormT= ', SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind), &
						!		' AND N2PerpphReNormT SUM= ', sum(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
						!			N2PerpphReNormT(nn, :, :, :)), &
						!		' FOR SPECIE= ', s, &
						!		', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
						!		' TIME-STEP= ', nn, ' IN ION PARTICLE COUNTS SUBROUTINE' &
						!		// achar(27) // '[0m.'
						!end if

						! ----------------------------------------------------

						! REDUCE ALL RAW PARTICLE COUNTS IN PHASE-SPACE TO MPI ROOT RANK (0):

						do Vperp1ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), 1
							do Vperp2ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), 1
								do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1

									! ----------------------------------------------------

									call mpi_barrier(MPI_COMM_WORLD, ierr)
									N2Perpph(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										N2PerpphReNormT(nn, Vperp1ind, Vperp2ind, Vparind)
									MPIRedSumIn(1)= N2Perpph(1)
									call MPIReduceSumSub
									N2PerpphR(1)= MPIRedSumOut(1)

									!call mpi_reduce(N2Perpph(1), N2PerpphR(1), 1, &
									!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

									if (rank == 0) then
										if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn) /= 0) then

											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%N2PerpphReNormRT(nn, Vperp1ind, Vperp2ind, Vparind)= &
													N2PerpphR(1)

										else if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn) == 0) then

											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												N2PerpphReNormRT(nn, Vperp1ind, Vperp2ind, Vparind)= 0d0

										end if

										! ----------------------------------------------------

										! FILTER STATISTICAL NOISE FOR UN-NORMALIZED ION COUNTS:

										if (SpecieT(s)%FluxTubeT(f)%IONNOISEflagT(1) == 1) then
											if (((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												N2PerpphReNormRT(nn, Vperp1ind, Vperp2ind, Vparind))**(-1d0/2d0)) > &
												SpecieT(s)%FluxTubeT(f)%IonNoiseLimitT(1)) then
												SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
													N2PerpphReNormRT(nn, Vperp1ind, Vperp2ind, Vparind)= 0d0
											end if
										end if

										! ----------------------------------------------------

										if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											N2PerpphReNormRT(nn, Vperp1ind, Vperp2ind, Vparind) == 0d0) then
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%N2PerpphRT(nn)= 0d0
										else if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											N2PerpphReNormRT(nn, Vperp1ind, Vperp2ind, Vparind) /= 0d0) then
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%N2PerpphRT(nn)= &
												SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												N2PerpphReNormRT(nn, Vperp1ind, Vperp2ind, Vparind)* &
												(SpecieT(s)%FluxTubeT(f)%nsnormfacT(1)/ &
												SpecieT(s)%FluxTubeT(f)%d3xCT(nn, Qind))
										end if

										SpecieT(s)%FluxTubeT(f)% &
											QCellT(Qind)%N2PerpphRTp(Vperp1ind, Vperp2ind, Vparind, nn)= &
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%N2PerpphRT(nn)

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR CONSISTENCY:

										if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											N2PerpphReNormRT(nn, Vperp1ind, Vperp2ind, Vparind) /= 0 .and. &
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%N2PerpphRT(nn) == 0) .or. &
											(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											N2PerpphReNormRT(nn, Vperp1ind, Vperp2ind, Vparind) == 0 .and. &
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%N2PerpphRT(nn) /= 0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' ZERO AND NON-ZERO N2PerpphRT AND N2PerpphReNormRT ELEMENTS DO NOT MATCH', &
												' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
												', Vperp1ind= ', Vperp1ind, ', Vperp2ind= ', Vperp2ind, ', Vparind= ', Vparind, &
												', AND STATISTICAL TIME-STEP= ', nn, &
												' IN ION PARTICLE COUNTS', &
												' SUBROUTINE' // achar(27) // '[0m.'
										end if

										if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											N2PerpphReNormRT(nn, Vperp1ind, Vperp2ind, Vparind))) .eqv. .true.) .or. &
											(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											N2PerpphReNormRT(:, :, :, :)) /= (SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)* &
											SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1)* &
											SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1)* &
											SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1))) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' N2PerpphReNormRT HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
												s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
												' TIME-STEP= ', nn, ' IN ION PARTICLE COUNTS SUBROUTINE' &
												// achar(27) // '[0m.'
										end if

										if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%N2PerpphRT(nn))) .eqv. .true.) .or. &
											(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%N2PerpphRT(:)) /= &
											(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' N2PerpphRT HAS', &
												' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
												', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, ', Vperp2ind= ', Vperp2ind, &
												', Vparind= ', Vparind, ', AND STATISTICAL TIME-STEP= ', nn, &
												' IN ION PARTICLE COUNTS SUBROUTINE' // achar(27) // '[0m.'
										end if

										if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%N2PerpphRT(nn) < 0) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NEGATIVE N2PerpphRT ELEMENTS FOR SPECIE= ', s, &
												', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', &
												Vperp1ind, ', Vperp2ind= ', Vperp2ind, ', Vparind= ', Vparind, &
												', AND STATISTICAL TIME-STEP= ', nn, ' IN ION PARTICLE COUNTS SUBROUTINE' &
												// achar(27) // '[0m.'
										end if

									end if

									! ----------------------------------------------------

								end do
							end do
						end do

						! ----------------------------------------------------

						if (rank == 0) then
							if (Qind == NqLB(1)) then
								!do Vperp1ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), 1
								!	do Vperp2ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), 1
								!		do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1
								!			if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%N2PerpphRTp(Vperp1ind, Vperp2ind, Vparind, nn) == 0d0) then
								!				write(*, *) SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%N2PerpphRTp(Vperp1ind, Vperp2ind, Vparind, nn)
								!			end if
								!		end do
								!	end do
								!end do
								if (sum(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%N2PerpphRTp(:, :, :, nn)) == 0d0) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' ZERO LOWER BOUNDARY ION COUNT SUM N2PerpphRTp= ', &
										sum(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%N2PerpphRTp(:, :, :, nn)), &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn), &
										' SPECIE= ', s, ', FLUX TUBE= ', f, &
										', Qind= ', Qind, ', AND STATISTICAL TIME-STEP= ', nn, &
										' IN ION PARTICLE COUNTS SUBROUTINE' // achar(27) // '[0m.'
								end if
							end if
						end if

					end do

					! ----------------------------------------------------

				end if
			end do

		end if

		! ----------------------------------------------------
		! ----------------------------------------------------

		! COMPUTE STATISTICAL GRID COUNTS FOR 2D ION VEL-SPACE GRIDS:

		! ----------------------------------------------------
		! ----------------------------------------------------

		if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 0) then

			do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
				if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
					(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1)))) then

					! ----------------------------------------------------

					! COMPUTE LOGICAL FLAGS FOR PHASE-SPACE STATISTICAL BINNING:

					do Qind= NqLB(1), NqUB(1), 1

						! ----------------------------------------------------

						if (SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind) == 0d0) then
							SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								NphReNormT(nn, :, :)= 0d0
						else if (SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind) /= 0d0) then

							! ----------------------------------------------------

							! COMPUTE RAW ION PHASE-SPACE COUNTS:

							if ((n == 1) .and. (nn == 1)) then
								do Vperpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), 1
									do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1
										jcount= 0d0
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												NphReNormT(nn, Vperpind, Vparind)= jcount
										jloopiv1p: do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
											if (ENAflag(j) .eqv. .false.) then
												if ((Vperpind == Vperpindk0(j)) .and. (Vparind == Vparindk0(j)) &
													.and. (Qind == Qindk0(j))) then

													jcount= jcount+ 1d0
													SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
															NphReNormT(nn, Vperpind, Vparind)= jcount

													cycle jloopiv1p

												end if
											end if
										end do jloopiv1p
									end do
								end do
							end if
							if ((n /= 1) .and. (nn /= 1) .and. &
								(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1))) then
								do Vperpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), 1
									do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1
										jcount= 0d0
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												NphReNormT(nn, Vperpind, Vparind)= jcount
										jloopiv2p: do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
											if (ENAflag(j) .eqv. .false.) then
											if ((Vperpind == Vperpindk1(j)) .and. (Vparind == Vparindk1(j)) &
												.and. (Qind == Qindk1(j))) then

													jcount= jcount+ 1d0
													SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
															NphReNormT(nn, Vperpind, Vparind)= jcount

													cycle jloopiv2p

												end if
											end if
										end do jloopiv2p
									end do
								end do
							end if

							! ----------------------------------------------------

						end if

						! ----------------------------------------------------

					end do

					! ----------------------------------------------------

					do Qind= NqLB(1), NqUB(1), 1

						! ----------------------------------------------------

						! DIAGNOSTIC FLAGS FOR CONSISTENCY:

						if ((sum(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
							NphReNormT(nn, :, :))- &
							SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind)) > 1d0) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
								' UNEQUAL NqReNormT= ', SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind), &
								' AND NphReNormT SUM= ', sum(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									NphReNormT(nn, :, :)), &
								' FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
								' TIME-STEP= ', nn, ' IN ION VEL PARTICLE COUNTS SUBROUTINE' &
								// achar(27) // '[0m.'
						end if

						! ----------------------------------------------------

						! REDUCE ALL RAW PARTICLE COUNTS IN PHASE-SPACE TO MPI ROOT RANK (0):

						do Vperpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), 1
							do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1

								! ----------------------------------------------------

								call mpi_barrier(MPI_COMM_WORLD, ierr)
								Nph(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									NphReNormT(nn, Vperpind, Vparind)
								MPIRedSumIn(1)= Nph(1)
								call MPIReduceSumSub
								NphR(1)= MPIRedSumOut(1)

								!call mpi_reduce(Nph(1), NphR(1), 1, &
								!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

								if (rank == 0) then
									if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn) /= 0) then

										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NphReNormRT(nn, Vperpind, Vparind)= &
												NphR(1)

									else if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn) == 0) then

										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											NphReNormRT(nn, Vperpind, Vparind)= 0d0

									end if

									! ----------------------------------------------------

									! FILTER STATISTICAL NOISE FOR UN-NORMALIZED ION COUNTS:

									if (SpecieT(s)%FluxTubeT(f)%IONNOISEflagT(1) == 1) then
										if (((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											NphReNormRT(nn, Vperpind, Vparind))**(-1d0/2d0)) > &
											SpecieT(s)%FluxTubeT(f)%IonNoiseLimitT(1)) then
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
												NphReNormRT(nn, Vperpind, Vparind)= 0d0
										end if
									end if

									! ----------------------------------------------------

									if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										NphReNormRT(nn, Vperpind, Vparind) == 0d0) then
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											VCellT(Vperpind, Vparind)%NphRT(nn)=	0d0
									else if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										NphReNormRT(nn, Vperpind, Vparind) /= 0d0) then
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											VCellT(Vperpind, Vparind)%NphRT(nn)= &
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
											NphReNormRT(nn, Vperpind, Vparind)* &
											(SpecieT(s)%FluxTubeT(f)%nsnormfacT(1)/ &
											SpecieT(s)%FluxTubeT(f)%d3xCT(nn, Qind))
									end if

									SpecieT(s)%FluxTubeT(f)% &
										QCellT(Qind)%NphRTp(Vperpind, Vparind, nn)= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%NphRT(nn)

									! ----------------------------------------------------

									! DIAGNOSTIC FLAGS FOR CONSISTENCY:

									if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										NphReNormRT(nn, Vperpind, Vparind) /= 0 .and. &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%NphRT(nn) == 0) .or. &
										(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										NphReNormRT(nn, Vperpind, Vparind) == 0 .and. &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%NphRT(nn) /= 0)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' ZERO AND NON-ZERO NphRT AND NphReNormRT ELEMENTS DO NOT MATCH', &
											' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
											', Vperpind= ', Vperpind, ', Vparind= ', Vparind, &
											', AND STATISTICAL TIME-STEP= ', nn, &
											' IN ION VEL PARTICLE COUNTS', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										NphReNormRT(nn, Vperpind, Vparind))) .eqv. .true.) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										NphReNormRT(:, :, :)) /= (SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)* &
										SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1)* &
										SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' NphReNormRT HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
											' TIME-STEP= ', nn, ' IN ION VEL PARTICLE COUNTS SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%NphRT(nn))) .eqv. .true.) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%NphRT(:)) /= &
										(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NphRT HAS', &
											' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
											', Qind= ', Qind, ', Vperpind= ', Vperpind, ', Vparind= ', &
											Vparind, ', AND STATISTICAL TIME-STEP= ', nn, &
											' IN ION VEL PARTICLE COUNTS SUBROUTINE' // achar(27) // '[0m.'
									end if

									if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										VCellT(Vperpind, Vparind)%NphRT(nn) < 0) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' NEGATIVE NphRT ELEMENTS FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', &
											Vperpind, ', Vparind= ', Vparind, ', AND STATISTICAL', &
											' TIME-STEP= ', nn, ' IN ION VEL PARTICLE COUNTS SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

								end if

								! ----------------------------------------------------

							end do
						end do

						! ----------------------------------------------------

					end do

					! ----------------------------------------------------

				end if
			end do

		end if

		!if (rank == 0) then
		!	do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
		!		if (((n /= 1) .and. (nn /= 1) .and. &
		!			(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1)))) then
		!			do Qind= NqLB(1), NqUB(1), 1

		!				if (sum(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
		!					N2PerpphReNormRT(nn, :, :, :)) /= 0d0) then

		!					write(*, *) 'N2PerpphReNormRT= ', nn, Qind, sum(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
		!						N2PerpphReNormRT(nn, :, :, :))

		!				end if

		!			end do
		!		end if
		!	end do
		!end if

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine IonParticleCountsSub

end module IonParticleCounts
