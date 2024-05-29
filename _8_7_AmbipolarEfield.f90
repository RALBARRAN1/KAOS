module AmbipolarEfield

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.7 AMBIPOLAR ELECTRIC FIELD:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use MomentFilter

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! COMPUTE FIELD-ALIGNED AMBIPOLAR ELECTRIC FIELD [V/m] FOR ALL PARTICLES:

	subroutine AmbipolarEfieldSub

		! ----------------------------------------------------

		! APPLY NOISE FILTER ON ION DENSITY AND PARALLEL VELOCITY MOMENTS FOR AMBIPOLAR E FIELD CALCULATION:

		! ----------------------------------------------------

		if ((SpecieT(1)%FluxTubeT(1)%MOMENTFILTERflagT(1) == 1) .and. (rank == 0)) then

			! ----------------------------------------------------

			! FILTER ION DENSITY MOMENT ACCORDING TO MOVING AVERAGE POINT:

			if (SpecieT(s)%FluxTubeT(f)%NqICT(1) >= nint((M0MAfilterPt- 1d0)/2d0)) then
				do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
					if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
						(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1))))) then

						MAfilterPt(1)= M0MAfilterPt
						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
							SpecieT(s)%FluxTubeT(f)%MomentFiltInT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%M0phRT(nn- 1, Qind)
						end do
						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
							call MomentFilterSub
						end do
						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
							SpecieT(s)%FluxTubeT(f)%M0FiltAvrgRT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%MomentFiltOutT(nn, Qind)
						end do

					end if
				end do
			end if

			! ----------------------------------------------------

			! FILTER ION PARALLEL VELOCITY MOMENT ACCORDING TO MOVING AVERAGE POINT:

			if (SpecieT(s)%FluxTubeT(f)%NqICT(1) >= nint((M1ParMAfilterPt- 1d0)/2d0)) then
				do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
					if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
						(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1))))) then

						MAfilterPt(1)= M1ParMAfilterPt
						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
							SpecieT(s)%FluxTubeT(f)%MomentFiltInT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn- 1, Qind)
						end do
						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
							call MomentFilterSub
						end do
						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
							SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%MomentFiltOutT(nn, Qind)
						end do

					end if
				end do
			end if

			! ----------------------------------------------------

		end if

		! ----------------------------------------------------

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
				(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1))))) then

				! ----------------------------------------------------

				if ((SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1) == 1) .or. &
					((SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1) == 0) .and. &
					((SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 1) .or. (dNTe /= 0d0)) .and. &
					((n /= 1) .and. (nn /= 1) .and. &
					(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1)))))) then

					! ----------------------------------------------------

					! COMPUTE AMBIPOLAR ELECTRIC FIELD MAGNITUDES:

					! ----------------------------------------------------

					if ((rank == 0) .and. &
						(SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 1)) then

						! ----------------------------------------------------

						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

							if (SpecieT(s)%FluxTubeT(f)%M0phRT(nn- 1, Qind) == 0d0) then
								SpecieT(s)%FluxTubeT(f)%LambdaDRT(nn, Qind)= 0d0
							else
								SpecieT(s)%FluxTubeT(f)%LambdaDRT(nn, Qind)= &
									sqrt((epsilon0*kB*SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn))/ &
									(SpecieT(s)%FluxTubeT(f)%M0phRT(nn- 1, Qind)*SpecieT(s)%qsT(1)**2d0))
							end if

							! ----------------------------------------------------

							! ENSURE FLUX-TUBE GRID SEPARATION IS LARGER THAN DEBYE LENGTH:

							if (SpecieT(s)%FluxTubeT(f)%M0phRT(nn- 1, Qind) /= 0d0) then
								if ((SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)* &
									SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind)) < &
									SpecieT(s)%FluxTubeT(f)%LambdaDRT(nn, Qind)) then
		 								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
		 									' Debye Length [m]= ', SpecieT(s)%FluxTubeT(f)%LambdaDRT(nn, Qind), &
		 									' IS GREATER THAN FLUX-TUBE GRID SEPARATION= ', &
											(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)* &
											SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind)), &
											' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
											', AND MASTER TIME-STEP= ', nn, &
											' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' &
		 									// achar(27) // '[0m.'
								end if
							end if

							! ----------------------------------------------------

						end do

						! ----------------------------------------------------

						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

							! ----------------------------------------------------

							if (SpecieT(s)%FluxTubeT(f)%EAINERTIALflagT(1) == 0) then
								EAInertialR(1)= 0d0

								SpecieT(s)%FluxTubeT(f)%EAInertialRT(nn, Qind)= abs(EAInertialR(1))
							end if
							if (SpecieT(s)%FluxTubeT(f)%EAPRESSUREflagT(1) == 0) then
								EAPressureR(1)= 0d0

								if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
									SpecieT(s)%FluxTubeT(f)%EAPressureRT(nn, Qind)= abs(EAPressureR(1))
								end if
								if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
									SpecieT(s)%FluxTubeT(f)%EAPressureRT(nn, Qind)= EAPressureR(1)
								end if
							end if
							if (SpecieT(s)%FluxTubeT(f)%M0phRT(nn- 1, Qind) == 0d0) then
								EAInertialR(1)= 0d0
								EAPressureR(1)= 0d0

								SpecieT(s)%FluxTubeT(f)%EAInertialRT(nn, Qind)= abs(EAInertialR(1))
								if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
									SpecieT(s)%FluxTubeT(f)%EAPressureRT(nn, Qind)= abs(EAPressureR(1))
								end if
								if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
									SpecieT(s)%FluxTubeT(f)%EAPressureRT(nn, Qind)= EAPressureR(1)
								end if
							end if

							! ----------------------------------------------------

							! COMPUTE INERTIAL TERM OF AMBIPOLAR ELECTRIC FIELD:

							! ----------------------------------------------------

							! Compute Five-Point Left Endpoint Derivative
							if ((SpecieT(s)%FluxTubeT(f)%M0phRT(nn- 1, Qind) /= 0d0) &
								.and. (Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) &
								.and. (SpecieT(s)%FluxTubeT(f)%EAINERTIALflagT(1) == 1)) then
								if ((n == 1) .and. (nn == 1)) then

									EAInertialR1(1)= (melec/SpecieT(s)%qsT(1))* &
										abs(SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn- 1, Qind))/ &
										(SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)*SpecieT(s)%hT)
									EAInertialR2(1)= (melec/SpecieT(s)%qsT(1))* &
										abs(SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn- 1, Qind))* &
										(1d0/(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
										((abs(-25d0*SpecieT(s)%FluxTubeT(f)% &
										M1ParphRT(nn- 1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
										M1ParphRT(nn- 1, Qind+ 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
										M1ParphRT(nn- 1, Qind+ 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
										M1ParphRT(nn- 1, Qind+ 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
										M1ParphRT(nn- 1, Qind+ 4)))/ &
										(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind)))

									EAInertialR(1)= EAInertialR1(1)+ EAInertialR2(1)

									SpecieT(s)%FluxTubeT(f)%EAInertialRT(nn, Qind)= abs(EAInertialR(1))

								end if
								if ((n /= 1) .and. (nn /= 1)) then
									if (SpecieT(1)%FluxTubeT(1)%MOMENTFILTERflagT(1) == 0) then
										if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 1) then

											EAInertialR1(1)= (melec/SpecieT(s)%qsT(1))*&
												abs(SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn- 1, Qind)- &
												SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn- 2, Qind))/ &
												(SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)*SpecieT(s)%hT)
											EAInertialR2(1)= (melec/SpecieT(s)%qsT(1))* &
												abs(SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn- 1, Qind))* &
												(1d0/(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
												((abs(-25d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(nn- 1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(nn- 1, Qind+ 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(nn- 1, Qind+ 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(nn- 1, Qind+ 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(nn- 1, Qind+ 4)))/ &
												(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind)))

										end if
										if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 0) then

											EAInertialR1(1)= (melec/SpecieT(s)%qsT(1))*&
												abs(SpecieT(s)%FluxTubeT(f)%M1ParphRT(2, Qind)- &
												SpecieT(s)%FluxTubeT(f)%M1ParphRT(1, Qind))/ &
												(SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)*SpecieT(s)%hT)
											EAInertialR2(1)= (melec/SpecieT(s)%qsT(1))* &
												abs(SpecieT(s)%FluxTubeT(f)%M1ParphRT(1, Qind))* &
												(1d0/(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
												((abs(-25d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(1, Qind+ 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(1, Qind+ 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(1, Qind+ 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(1, Qind+ 4)))/ &
												(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind)))

										end if
									end if
									if (SpecieT(1)%FluxTubeT(1)%MOMENTFILTERflagT(1) == 1) then
										if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 1) then

											EAInertialR1(1)= (melec/SpecieT(s)%qsT(1))*&
												abs(SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(nn, Qind)- &
												SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(nn- 1, Qind))/ &
												(SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)*SpecieT(s)%hT)
											EAInertialR2(1)= (melec/SpecieT(s)%qsT(1))* &
												abs(SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(nn, Qind))* &
												(1d0/(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
												((abs(-25d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(nn, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(nn, Qind+ 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(nn, Qind+ 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(nn, Qind+ 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(nn, Qind+ 4)))/ &
												(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind)))

										end if
										if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 0) then

											EAInertialR1(1)= (melec/SpecieT(s)%qsT(1))*&
												abs(SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(2, Qind)- &
												SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(1, Qind))/ &
												(SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)*SpecieT(s)%hT)
											EAInertialR2(1)= (melec/SpecieT(s)%qsT(1))* &
												abs(SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(1, Qind))* &
												(1d0/(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
												((abs(-25d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(1, Qind+ 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(1, Qind+ 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(1, Qind+ 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(1, Qind+ 4)))/ &
												(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind)))

										end if
									end if

									EAInertialR(1)= EAInertialR1(1)+ EAInertialR2(1)

									SpecieT(s)%FluxTubeT(f)%EAInertialRT(nn, Qind)= abs(EAInertialR(1))

								end if
							end if

							! Compute Five-Point Right Endpoint Derivative
							if ((SpecieT(s)%FluxTubeT(f)%M0phRT(nn- 1, Qind) /= 0d0) &
								.and. (Qind == SpecieT(s)%FluxTubeT(f)%NqUBT(1)) &
								.and. (SpecieT(s)%FluxTubeT(f)%EAINERTIALflagT(1) == 1)) then
								if ((n == 1) .and. (nn == 1)) then

									!EAInertialR1(1)= (melec/SpecieT(s)%qsT(1))* &
									!	abs(SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn- 1, Qind))/ &
									!	(SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)*SpecieT(s)%hT)
									!EAInertialR2(1)= (melec/SpecieT(s)%qsT(1))* &
									!	abs(SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn- 1, Qind))* &
									!	(1d0/(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
									!	((abs(-25d0*SpecieT(s)%FluxTubeT(f)% &
									!	M1ParphRT(nn- 1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
									!	M1ParphRT(nn- 1, Qind- 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
									!	M1ParphRT(nn- 1, Qind- 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
									!	M1ParphRT(nn- 1, Qind- 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
									!	M1ParphRT(nn- 1, Qind- 4)))/ &
									!	(-12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind)))

									!EAInertialR(1)= EAInertialR1(1)+ EAInertialR2(1)

									EAInertialR(1)= 0d0

									SpecieT(s)%FluxTubeT(f)%EAInertialRT(nn, Qind)= abs(EAInertialR(1))

								end if
								if ((n /= 1) .and. (nn /= 1)) then
									if (SpecieT(1)%FluxTubeT(1)%MOMENTFILTERflagT(1) == 0) then
										if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 1) then

											EAInertialR1(1)= (melec/SpecieT(s)%qsT(1))*&
												abs(SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn- 1, Qind)- &
												SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn- 2, Qind))/ &
												(SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)*SpecieT(s)%hT)
											EAInertialR2(1)= (melec/SpecieT(s)%qsT(1))* &
												abs(SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn- 1, Qind))* &
												(1d0/(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
												((abs(-25d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(nn- 1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(nn- 1, Qind- 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(nn- 1, Qind- 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(nn- 1, Qind- 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(nn- 1, Qind- 4)))/ &
												(-12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind)))

										end if
										if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 0) then

											EAInertialR1(1)= (melec/SpecieT(s)%qsT(1))*&
												abs(SpecieT(s)%FluxTubeT(f)%M1ParphRT(2, Qind)- &
												SpecieT(s)%FluxTubeT(f)%M1ParphRT(1, Qind))/ &
												(SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)*SpecieT(s)%hT)
											EAInertialR2(1)= (melec/SpecieT(s)%qsT(1))* &
												abs(SpecieT(s)%FluxTubeT(f)%M1ParphRT(1, Qind))* &
												(1d0/(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
												((abs(-25d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(1, Qind- 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(1, Qind- 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(1, Qind- 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(1, Qind- 4)))/ &
												(-12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind)))

										end if
									end if
									if (SpecieT(1)%FluxTubeT(1)%MOMENTFILTERflagT(1) == 1) then
										if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 1) then

											EAInertialR1(1)= (melec/SpecieT(s)%qsT(1))*&
												abs(SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(nn, Qind)- &
												SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(nn- 1, Qind))/ &
												(SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)*SpecieT(s)%hT)
											EAInertialR2(1)= (melec/SpecieT(s)%qsT(1))* &
												abs(SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(nn, Qind))* &
												(1d0/(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
												((abs(-25d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(nn, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(nn, Qind- 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(nn, Qind- 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(nn, Qind- 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(nn, Qind- 4)))/ &
												(-12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind)))

										end if
										if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 0) then

											EAInertialR1(1)= (melec/SpecieT(s)%qsT(1))*&
												abs(SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(2, Qind)- &
												SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(1, Qind))/ &
												(SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)*SpecieT(s)%hT)
											EAInertialR2(1)= (melec/SpecieT(s)%qsT(1))* &
												abs(SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(1, Qind))* &
												(1d0/(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
												((abs(-25d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(1, Qind- 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(1, Qind- 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(1, Qind- 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(1, Qind- 4)))/ &
												(-12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind)))

										end if
									end if

									!EAInertialR(1)= EAInertialR1(1)+ EAInertialR2(1)

									EAInertialR(1)= 0d0

									SpecieT(s)%FluxTubeT(f)%EAInertialRT(nn, Qind)= abs(EAInertialR(1))

								end if
							end if

							! Compute Three-Point Midpoint Derivative
							if ((SpecieT(s)%FluxTubeT(f)%M0phRT(nn- 1, Qind) /= 0d0) &
								.and. ((Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1) .or. (Qind == SpecieT(s)%FluxTubeT(f)%NqUBT(1)- 1)) &
								.and. (SpecieT(s)%FluxTubeT(f)%EAINERTIALflagT(1) == 1)) then
								if ((n == 1) .and. (nn == 1)) then

									EAInertialR1(1)= (melec/SpecieT(s)%qsT(1))* &
										abs(SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn- 1, Qind))/ &
										(SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)*SpecieT(s)%hT)
									EAInertialR2(1)= (melec/SpecieT(s)%qsT(1))* &
										abs(SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn- 1, Qind))* &
										(1d0/(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
										((abs(SpecieT(s)%FluxTubeT(f)% &
										M1ParphRT(nn- 1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
										M1ParphRT(nn- 1, Qind- 1)))/ &
										(2d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind)))

									EAInertialR(1)= EAInertialR1(1)+ EAInertialR2(1)

									SpecieT(s)%FluxTubeT(f)%EAInertialRT(nn, Qind)= abs(EAInertialR(1))

								end if
								if ((n /= 1) .and. (nn /= 1)) then
									if (SpecieT(1)%FluxTubeT(1)%MOMENTFILTERflagT(1) == 0) then
										if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 1) then

											EAInertialR1(1)= (melec/SpecieT(s)%qsT(1))*&
												abs(SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn- 1, Qind)- &
												SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn- 2, Qind))/ &
												(SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)*SpecieT(s)%hT)
											EAInertialR2(1)= (melec/SpecieT(s)%qsT(1))* &
												abs(SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn- 1, Qind))* &
												(1d0/(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
												((abs(SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(nn- 1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(nn- 1, Qind- 1)))/ &
												(2d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind)))

										end if
										if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 0) then

											EAInertialR1(1)= (melec/SpecieT(s)%qsT(1))*&
												abs(SpecieT(s)%FluxTubeT(f)%M1ParphRT(2, Qind)- &
												SpecieT(s)%FluxTubeT(f)%M1ParphRT(1, Qind))/ &
												(SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)*SpecieT(s)%hT)
											EAInertialR2(1)= (melec/SpecieT(s)%qsT(1))* &
												abs(SpecieT(s)%FluxTubeT(f)%M1ParphRT(1, Qind))* &
												(1d0/(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
												((abs(SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(1, Qind- 1)))/ &
												(2d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind)))

										end if
									end if
									if (SpecieT(1)%FluxTubeT(1)%MOMENTFILTERflagT(1) == 1) then
										if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 1) then

											EAInertialR1(1)= (melec/SpecieT(s)%qsT(1))*&
												abs(SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(nn, Qind)- &
												SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(nn- 1, Qind))/ &
												(SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)*SpecieT(s)%hT)
											EAInertialR2(1)= (melec/SpecieT(s)%qsT(1))* &
												abs(SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(nn, Qind))* &
												(1d0/(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
												((abs(SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(nn, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(nn, Qind- 1)))/ &
												(2d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind)))

										end if
										if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 0) then

											EAInertialR1(1)= (melec/SpecieT(s)%qsT(1))*&
												abs(SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(2, Qind)- &
												SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(1, Qind))/ &
												(SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)*SpecieT(s)%hT)
											EAInertialR2(1)= (melec/SpecieT(s)%qsT(1))* &
												abs(SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(1, Qind))* &
												(1d0/(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
												((abs(SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(1, Qind- 1)))/ &
												(2d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind)))

										end if
									end if

									EAInertialR(1)= EAInertialR1(1)+ EAInertialR2(1)

									SpecieT(s)%FluxTubeT(f)%EAInertialRT(nn, Qind)= abs(EAInertialR(1))

								end if
							end if

							! Compute Five-Point Midpoint Derivative
							if ((SpecieT(s)%FluxTubeT(f)%M0phRT(nn- 1, Qind) /= 0d0) &
								.and. (Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. (Qind /= SpecieT(s)%FluxTubeT(f)%NqUBT(1)) &
								.and. (Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1) .and. (Qind /= SpecieT(s)%FluxTubeT(f)%NqUBT(1)- 1) &
								.and. (SpecieT(s)%FluxTubeT(f)%EAINERTIALflagT(1) == 1)) then
								if ((n == 1) .and. (nn == 1)) then

									EAInertialR1(1)= (melec/SpecieT(s)%qsT(1))* &
										abs(SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn- 1, Qind))/ &
										(SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)*SpecieT(s)%hT)
									EAInertialR2(1)= (melec/SpecieT(s)%qsT(1))* &
										abs(SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn- 1, Qind))* &
										(1d0/(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
										(abs(SpecieT(s)%FluxTubeT(f)% &
										M1ParphRT(nn- 1, Qind- 2)- 8d0*SpecieT(s)%FluxTubeT(f)% &
										M1ParphRT(nn- 1, Qind- 1)+ 8d0*SpecieT(s)%FluxTubeT(f)% &
										M1ParphRT(nn- 1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
										M1ParphRT(nn- 1, Qind+ 2)))/ &
										(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))

									EAInertialR(1)= EAInertialR1(1)+ EAInertialR2(1)

									SpecieT(s)%FluxTubeT(f)%EAInertialRT(nn, Qind)= abs(EAInertialR(1))

								end if
								if ((n /= 1) .and. (nn /= 1)) then
									if (SpecieT(1)%FluxTubeT(1)%MOMENTFILTERflagT(1) == 0) then
										if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 1) then

											EAInertialR1(1)= (melec/SpecieT(s)%qsT(1))*&
												abs(SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn- 1, Qind)- &
												SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn- 2, Qind))/ &
												(SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)*SpecieT(s)%hT)
											EAInertialR2(1)= (melec/SpecieT(s)%qsT(1))* &
												abs(SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn- 1, Qind))* &
												(1d0/(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
												(abs(SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(nn- 1, Qind- 2)- 8d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(nn- 1, Qind- 1)+ 8d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(nn- 1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(nn- 1, Qind+ 2)))/ &
												(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))

										end if
										if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 0) then

											EAInertialR1(1)= (melec/SpecieT(s)%qsT(1))*&
												abs(SpecieT(s)%FluxTubeT(f)%M1ParphRT(2, Qind)- &
												SpecieT(s)%FluxTubeT(f)%M1ParphRT(1, Qind))/ &
												(SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)*SpecieT(s)%hT)
											EAInertialR2(1)= (melec/SpecieT(s)%qsT(1))* &
												abs(SpecieT(s)%FluxTubeT(f)%M1ParphRT(1, Qind))* &
												(1d0/(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
												(abs(SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(1, Qind- 2)- 8d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(1, Qind- 1)+ 8d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
												M1ParphRT(1, Qind+ 2)))/ &
												(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))

										end if
									end if
									if (SpecieT(1)%FluxTubeT(1)%MOMENTFILTERflagT(1) == 1) then
										if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 1) then

											EAInertialR1(1)= (melec/SpecieT(s)%qsT(1))*&
												abs(SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(nn, Qind)- &
												SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(nn- 1, Qind))/ &
												(SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)*SpecieT(s)%hT)
											EAInertialR2(1)= (melec/SpecieT(s)%qsT(1))* &
												abs(SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(nn, Qind))* &
												(1d0/(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
												(abs(SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(nn, Qind- 2)- 8d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(nn, Qind- 1)+ 8d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(nn, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(nn, Qind+ 2)))/ &
												(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))

										end if
										if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 0) then

											EAInertialR1(1)= (melec/SpecieT(s)%qsT(1))*&
												abs(SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(2, Qind)- &
												SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(1, Qind))/ &
												(SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)*SpecieT(s)%hT)
											EAInertialR2(1)= (melec/SpecieT(s)%qsT(1))* &
												abs(SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(1, Qind))* &
												(1d0/(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
												(abs(SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(1, Qind- 2)- 8d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(1, Qind- 1)+ 8d0*SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
												M1ParFiltAvrgRT(1, Qind+ 2)))/ &
												(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))

										end if
									end if

									EAInertialR(1)= EAInertialR1(1)+ EAInertialR2(1)

									SpecieT(s)%FluxTubeT(f)%EAInertialRT(nn, Qind)= abs(EAInertialR(1))

								end if
							end if

							! ----------------------------------------------------

						end do

						! ----------------------------------------------------

						! COMPUTE PRESSURE TERM OF AMBIPOLAR ELECTRIC FIELD:

						! ----------------------------------------------------

						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

							! ----------------------------------------------------

							! Compute Five-Point Left Endpoint Derivative
							if ((SpecieT(s)%FluxTubeT(f)%M0phRT(nn- 1, Qind) /= 0d0) &
								.and. (Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) &
								.and. (SpecieT(s)%FluxTubeT(f)%EAPRESSUREflagT(1) == 1)) then

								if (SpecieT(1)%FluxTubeT(1)%MOMENTFILTERflagT(1) == 0) then
									if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 1) then

										EAPressureR1(1)= (1d0/abs(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
											(kB*SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn)/ &
											(SpecieT(s)%qsT(1)*SpecieT(s)%FluxTubeT(f)%M0phRT(nn- 1, Qind)))

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
											EAPressureR2(1)= (abs(-25d0*SpecieT(s)%FluxTubeT(f)% &
												M0phRT(nn- 1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
												M0phRT(nn- 1, Qind+ 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
												M0phRT(nn- 1, Qind+ 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
												M0phRT(nn- 1, Qind+ 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
												M0phRT(nn- 1, Qind+ 4)))/ &
												(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
										end if

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then ! S. Magnetic Hemisphere
												EAPressureR2(1)= -((-25d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind+ 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind+ 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind+ 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind+ 4)))/ &
													(-12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then ! N. Magnetic Hemisphere
												EAPressureR2(1)= -((-25d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind+ 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind+ 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind+ 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind+ 4)))/ &
													(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
										end if

									end if
									if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 0) then

										EAPressureR1(1)= (1d0/abs(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
											(kB*SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn)/ &
											(SpecieT(s)%qsT(1)*SpecieT(s)%FluxTubeT(f)%M0phRT(1, Qind)))

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
											EAPressureR2(1)= (abs(-25d0*SpecieT(s)%FluxTubeT(f)% &
												M0phRT(1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
												M0phRT(1, Qind+ 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
												M0phRT(1, Qind+ 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
												M0phRT(1, Qind+ 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
												M0phRT(1, Qind+ 4)))/ &
												(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
										end if

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then ! S. Magnetic Hemisphere
												EAPressureR2(1)= -((-25d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind+ 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind+ 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind+ 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind+ 4)))/ &
													(-12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then ! N. Magnetic Hemisphere
												EAPressureR2(1)= -((-25d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind+ 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind+ 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind+ 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind+ 4)))/ &
													(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
										end if

									end if
								end if
								if (SpecieT(1)%FluxTubeT(1)%MOMENTFILTERflagT(1) == 1) then
									if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 1) then

										EAPressureR1(1)= (1d0/abs(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
											(kB*SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn)/ &
											(SpecieT(s)%qsT(1)*SpecieT(s)%FluxTubeT(f)%M0FiltAvrgRT(nn, Qind)))

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
											EAPressureR2(1)= (abs(-25d0*SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(nn, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(nn, Qind+ 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(nn, Qind+ 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(nn, Qind+ 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(nn, Qind+ 4)))/ &
												(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
										end if

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then ! S. Magnetic Hemisphere
												EAPressureR2(1)= -((-25d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind+ 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind+ 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind+ 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind+ 4)))/ &
													(-12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then ! N. Magnetic Hemisphere
												EAPressureR2(1)= -((-25d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind+ 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind+ 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind+ 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind+ 4)))/ &
													(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
										end if

									end if
									if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 0) then

										EAPressureR1(1)= (1d0/abs(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
											(kB*SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn)/ &
											(SpecieT(s)%qsT(1)*SpecieT(s)%FluxTubeT(f)%M0FiltAvrgRT(1, Qind)))

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
											EAPressureR2(1)= (abs(-25d0*SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(1, Qind+ 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(1, Qind+ 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(1, Qind+ 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(1, Qind+ 4)))/ &
												(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
										end if

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then ! S. Magnetic Hemisphere
												EAPressureR2(1)= -((-25d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind+ 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind+ 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind+ 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind+ 4)))/ &
													(-12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then ! N. Magnetic Hemisphere
												EAPressureR2(1)= -((-25d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind+ 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind+ 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind+ 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind+ 4)))/ &
													(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
										end if

									end if
								end if

								EAPressureR(1)= EAPressureR1(1)*EAPressureR2(1)

								if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
									SpecieT(s)%FluxTubeT(f)%EAPressureRT(nn, Qind)= abs(EAPressureR(1))
								end if
								if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
									SpecieT(s)%FluxTubeT(f)%EAPressureRT(nn, Qind)= EAPressureR(1)
								end if

							end if

							! Compute Five-Point Right Endpoint Derivative
							if ((SpecieT(s)%FluxTubeT(f)%M0phRT(nn- 1, Qind) /= 0d0) &
								.and. (Qind == SpecieT(s)%FluxTubeT(f)%NqUBT(1)) &
								.and. (SpecieT(s)%FluxTubeT(f)%EAPRESSUREflagT(1) == 1)) then

								if (SpecieT(1)%FluxTubeT(1)%MOMENTFILTERflagT(1) == 0) then
									if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 1) then

										EAPressureR1(1)= (1d0/abs(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
											(kB*SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn)/ &
											(SpecieT(s)%qsT(1)*SpecieT(s)%FluxTubeT(f)%M0phRT(nn- 1, Qind)))

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
											EAPressureR2(1)= (abs(-25d0*SpecieT(s)%FluxTubeT(f)% &
												M0phRT(nn- 1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
												M0phRT(nn- 1, Qind- 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
												M0phRT(nn- 1, Qind- 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
												M0phRT(nn- 1, Qind- 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
												M0phRT(nn- 1, Qind- 4)))/ &
												(-12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
										end if

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then ! S. Magnetic Hemisphere
												EAPressureR2(1)= -((-25d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind- 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind- 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind- 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind- 4)))/ &
													(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then ! N. Magnetic Hemisphere
												EAPressureR2(1)= -((-25d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind- 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind- 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind- 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind- 4)))/ &
													(-12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
										end if

									end if
									if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 0) then

										EAPressureR1(1)= (1d0/abs(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
											(kB*SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn)/ &
											(SpecieT(s)%qsT(1)*SpecieT(s)%FluxTubeT(f)%M0phRT(1, Qind)))

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
											EAPressureR2(1)= (abs(-25d0*SpecieT(s)%FluxTubeT(f)% &
												M0phRT(1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
												M0phRT(1, Qind- 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
												M0phRT(1, Qind- 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
												M0phRT(1, Qind- 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
												M0phRT(1, Qind- 4)))/ &
												(-12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
										end if

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then ! S. Magnetic Hemisphere
												EAPressureR2(1)= -((-25d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind- 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind- 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind- 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind- 4)))/ &
													(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then ! N. Magnetic Hemisphere
												EAPressureR2(1)= -((-25d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind- 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind- 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind- 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind- 4)))/ &
													(-12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
										end if

									end if
								end if
								if (SpecieT(1)%FluxTubeT(1)%MOMENTFILTERflagT(1) == 1) then
									if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 1) then

										EAPressureR1(1)= (1d0/abs(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
											(kB*SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn)/ &
											(SpecieT(s)%qsT(1)*SpecieT(s)%FluxTubeT(f)%M0FiltAvrgRT(nn, Qind)))

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
											EAPressureR2(1)= (abs(-25d0*SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(nn, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(nn, Qind- 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(nn, Qind- 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(nn, Qind- 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(nn, Qind- 4)))/ &
												(-12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
										end if

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then ! S. Magnetic Hemisphere
												EAPressureR2(1)= -((-25d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind- 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind- 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind- 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind- 4)))/ &
													(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then ! N. Magnetic Hemisphere
												EAPressureR2(1)= -((-25d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind- 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind- 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind- 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind- 4)))/ &
													(-12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
										end if

									end if
									if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 0) then

										EAPressureR1(1)= (1d0/abs(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
											(kB*SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn)/ &
											(SpecieT(s)%qsT(1)*SpecieT(s)%FluxTubeT(f)%M0FiltAvrgRT(1, Qind)))

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
											EAPressureR2(1)= (abs(-25d0*SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(1, Qind- 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(1, Qind- 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(1, Qind- 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(1, Qind- 4)))/ &
												(-12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
										end if

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then ! S. Magnetic Hemisphere
												EAPressureR2(1)= -((-25d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind- 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind- 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind- 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind- 4)))/ &
													(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then ! N. Magnetic Hemisphere
												EAPressureR2(1)= -((-25d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind)+ 48d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind- 1)- 36d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind- 2)+ 16d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind- 3)- 3d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind- 4)))/ &
													(-12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
										end if

									end if
								end if

								EAPressureR(1)= EAPressureR1(1)*EAPressureR2(1)

								if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
									SpecieT(s)%FluxTubeT(f)%EAPressureRT(nn, Qind)= abs(EAPressureR(1))
								end if
								if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
									SpecieT(s)%FluxTubeT(f)%EAPressureRT(nn, Qind)= EAPressureR(1)
								end if

							end if

							! Compute Three-Point Midpoint Derivative
							if ((SpecieT(s)%FluxTubeT(f)%M0phRT(nn- 1, Qind) /= 0d0) &
								.and. ((Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1) .or. (Qind == SpecieT(s)%FluxTubeT(f)%NqUBT(1)- 1)) &
								.and. (SpecieT(s)%FluxTubeT(f)%EAPRESSUREflagT(1) == 1)) then

								if (SpecieT(1)%FluxTubeT(1)%MOMENTFILTERflagT(1) == 0) then
									if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 1) then

										EAPressureR1(1)= (1d0/abs(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
											(kB*SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn)/ &
											(SpecieT(s)%qsT(1)*SpecieT(s)%FluxTubeT(f)%M0phRT(nn- 1, Qind)))

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
											EAPressureR2(1)= (abs(SpecieT(s)%FluxTubeT(f)% &
												M0phRT(nn- 1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
												M0phRT(nn- 1, Qind- 1)))/ &
												(2d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
										end if

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then ! S. Magnetic Hemisphere
												EAPressureR2(1)= -((SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind- 1)))/ &
													(-2d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then ! N. Magnetic Hemisphere
												EAPressureR2(1)= -((SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind- 1)))/ &
													(2d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
										end if

									end if
									if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 0) then

										EAPressureR1(1)= (1d0/abs(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
											(kB*SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn)/ &
											(SpecieT(s)%qsT(1)*SpecieT(s)%FluxTubeT(f)%M0phRT(1, Qind)))

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
											EAPressureR2(1)= (abs(SpecieT(s)%FluxTubeT(f)% &
												M0phRT(1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
												M0phRT(1, Qind- 1)))/ &
												(2d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
										end if

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then ! S. Magnetic Hemisphere
												EAPressureR2(1)= -((SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind- 1)))/ &
													(-2d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then ! N. Magnetic Hemisphere
												EAPressureR2(1)= -((SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind- 1)))/ &
													(2d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
										end if

									end if
								end if
								if (SpecieT(1)%FluxTubeT(1)%MOMENTFILTERflagT(1) == 1) then
									if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 1) then

										EAPressureR1(1)= (1d0/abs(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
											(kB*SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn)/ &
											(SpecieT(s)%qsT(1)*SpecieT(s)%FluxTubeT(f)%M0FiltAvrgRT(nn, Qind)))

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
											EAPressureR2(1)= (abs(SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(nn, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(nn, Qind- 1)))/ &
												(2d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
										end if

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then ! S. Magnetic Hemisphere
												EAPressureR2(1)= -((SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind- 1)))/ &
													(-2d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then ! N. Magnetic Hemisphere
												EAPressureR2(1)= -((SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind- 1)))/ &
													(2d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
										end if

									end if
									if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 0) then

										EAPressureR1(1)= (1d0/abs(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
											(kB*SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn)/ &
											(SpecieT(s)%qsT(1)*SpecieT(s)%FluxTubeT(f)%M0FiltAvrgRT(1, Qind)))

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
											EAPressureR2(1)= (abs(SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(1, Qind- 1)))/ &
												(2d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
										end if

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then ! S. Magnetic Hemisphere
												EAPressureR2(1)= -((SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind- 1)))/ &
													(-2d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then ! N. Magnetic Hemisphere
												EAPressureR2(1)= -((SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind- 1)))/ &
													(2d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
										end if

									end if
								end if

								EAPressureR(1)= EAPressureR1(1)*EAPressureR2(1)

								if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
									SpecieT(s)%FluxTubeT(f)%EAPressureRT(nn, Qind)= abs(EAPressureR(1))
								end if
								if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
									SpecieT(s)%FluxTubeT(f)%EAPressureRT(nn, Qind)= EAPressureR(1)
								end if

							end if

							! Compute Five-Point Midpoint Derivative
							if ((SpecieT(s)%FluxTubeT(f)%M0phRT(nn- 1, Qind) /= 0d0) &
								.and. (Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. (Qind /= SpecieT(s)%FluxTubeT(f)%NqUBT(1)) &
								.and. (Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1) .and. (Qind /= SpecieT(s)%FluxTubeT(f)%NqUBT(1)- 1) &
								.and. (SpecieT(s)%FluxTubeT(f)%EAPRESSUREflagT(1) == 1)) then

								if (SpecieT(1)%FluxTubeT(1)%MOMENTFILTERflagT(1) == 0) then
									if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 1) then

										EAPressureR1(1)= (1d0/abs(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
											(kB*SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn)/ &
											(SpecieT(s)%qsT(1)*SpecieT(s)%FluxTubeT(f)%M0phRT(nn- 1, Qind)))

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
											EAPressureR2(1)= (abs(SpecieT(s)%FluxTubeT(f)% &
												M0phRT(nn- 1, Qind- 2)- 8d0*SpecieT(s)%FluxTubeT(f)% &
												M0phRT(nn- 1, Qind- 1)+ 8d0*SpecieT(s)%FluxTubeT(f)% &
												M0phRT(nn- 1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
												M0phRT(nn- 1, Qind+ 2)))/ &
												(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
										end if

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then ! S. Magnetic Hemisphere
												EAPressureR2(1)= -((SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind- 2)- 8d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind- 1)+ 8d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind+ 2)))/ &
													(-12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then ! N. Magnetic Hemisphere
												EAPressureR2(1)= -((SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind- 2)- 8d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind- 1)+ 8d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
													M0phRT(nn- 1, Qind+ 2)))/ &
													(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
										end if

									end if
									if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 0) then

										EAPressureR1(1)= (1d0/abs(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
											(kB*SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn)/ &
											(SpecieT(s)%qsT(1)*SpecieT(s)%FluxTubeT(f)%M0phRT(1, Qind)))

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
											EAPressureR2(1)= (abs(SpecieT(s)%FluxTubeT(f)% &
												M0phRT(1, Qind- 2)- 8d0*SpecieT(s)%FluxTubeT(f)% &
												M0phRT(1, Qind- 1)+ 8d0*SpecieT(s)%FluxTubeT(f)% &
												M0phRT(1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
												M0phRT(1, Qind+ 2)))/ &
												(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
										end if

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then ! S. Magnetic Hemisphere
												EAPressureR2(1)= -((SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind- 2)- 8d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind- 1)+ 8d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind+ 2)))/ &
													(-12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then ! N. Magnetic Hemisphere
												EAPressureR2(1)= -((SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind- 2)- 8d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind- 1)+ 8d0*SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
													M0phRT(1, Qind+ 2)))/ &
													(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
										end if

									end if
								end if
								if (SpecieT(1)%FluxTubeT(1)%MOMENTFILTERflagT(1) == 1) then
									if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 1) then

										EAPressureR1(1)= (1d0/abs(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
											(kB*SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn)/ &
											(SpecieT(s)%qsT(1)*SpecieT(s)%FluxTubeT(f)%M0FiltAvrgRT(nn, Qind)))

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
											EAPressureR2(1)= (abs(SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(nn, Qind- 2)- 8d0*SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(nn, Qind- 1)+ 8d0*SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(nn, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(nn, Qind+ 2)))/ &
												(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
										end if

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then ! S. Magnetic Hemisphere
												EAPressureR2(1)= -((SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind- 2)- 8d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind- 1)+ 8d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind+ 2)))/ &
													(-12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then ! N. Magnetic Hemisphere
												EAPressureR2(1)= -((SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind- 2)- 8d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind- 1)+ 8d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(nn, Qind+ 2)))/ &
													(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
										end if

									end if
									if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 0) then

										EAPressureR1(1)= (1d0/abs(SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)))* &
											(kB*SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn)/ &
											(SpecieT(s)%qsT(1)*SpecieT(s)%FluxTubeT(f)%M0FiltAvrgRT(1, Qind)))

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
											EAPressureR2(1)= (abs(SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(1, Qind- 2)- 8d0*SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(1, Qind- 1)+ 8d0*SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
												M0FiltAvrgRT(1, Qind+ 2)))/ &
												(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
										end if

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then ! S. Magnetic Hemisphere
												EAPressureR2(1)= -((SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind- 2)- 8d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind- 1)+ 8d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind+ 2)))/ &
													(-12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
											if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then ! N. Magnetic Hemisphere
												EAPressureR2(1)= -((SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind- 2)- 8d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind- 1)+ 8d0*SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind+ 1)- SpecieT(s)%FluxTubeT(f)% &
													M0FiltAvrgRT(1, Qind+ 2)))/ &
													(12d0*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind))
											end if
										end if

									end if
								end if

								EAPressureR(1)= EAPressureR1(1)*EAPressureR2(1)

								if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
									SpecieT(s)%FluxTubeT(f)%EAPressureRT(nn, Qind)= abs(EAPressureR(1))
								end if
								if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
									SpecieT(s)%FluxTubeT(f)%EAPressureRT(nn, Qind)= EAPressureR(1)
								end if

							end if

						end do

						! ----------------------------------------------------

						! Cutoff boundary self-consistent ambipolar field values:
						if (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 1) then
							do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
								if ((SpecieT(s)%FluxTubeT(f)%M0phRT(nn- 1, Qind) /= 0d0) &
									.and. ((Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .or. &
									(Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1) .or. (Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 2)) &
									.and. (SpecieT(s)%FluxTubeT(f)%EAPRESSUREflagT(1) == 1)) then

									SpecieT(s)%FluxTubeT(f)%EAPressureRT(nn, Qind)= &
										SpecieT(s)%FluxTubeT(f)%EAPressureRT(nn, SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3)

								end if
							end do
						end if

						! ----------------------------------------------------

						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

							! ----------------------------------------------------

							! Save electron inertial and pressure terms as derived data types:
							SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind)= &
							SpecieT(s)%FluxTubeT(f)%EAInertialRT(nn, Qind)+ &
								SpecieT(s)%FluxTubeT(f)%EAPressureRT(nn, Qind)

							! ----------------------------------------------------

							! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

							if ((isnan(real(SpecieT(s)%FluxTubeT(f)%EAInertialRT(nn, Qind))) &
								.eqv. .true.) .or. (size(SpecieT(s)%FluxTubeT(f)%EAInertialRT(:, :)) &
								/= (SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)* &
								((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' EAInertialRT= ', SpecieT(s)%FluxTubeT(f)%EAInertialRT(nn, Qind), &
									' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND MASTER', &
									' TIME-STEP= ', nn, ' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' &
									// achar(27) // '[0m.'
							end if

							if ((isnan(real(SpecieT(s)%FluxTubeT(f)%EAPressureRT(nn, Qind))) &
								.eqv. .true.) .or. (size(SpecieT(s)%FluxTubeT(f)%EAPressureRT(:, :)) &
								/= (SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)* &
								((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' EAPressureRT= ', SpecieT(s)%FluxTubeT(f)%EAPressureRT(nn, Qind), &
									' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND MASTER', &
									' TIME-STEP= ', nn, ' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' &
									// achar(27) // '[0m.'
							end if

							if ((isnan(real(SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind))) &
								.eqv. .true.) .or. (size(SpecieT(s)%FluxTubeT(f)%EAmagRT(:, :)) &
								/= (SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)* &
								((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' EAmagRT= ', SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind), &
									' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND MASTER', &
									' TIME-STEP= ', nn, ' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' &
									// achar(27) // '[0m.'
							end if

							if ((nn > 1) .and. (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 0) .and. (dNTe == 0d0)) then
								if (SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind) /= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn- 1, Qind)) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' EAmagRT= ', SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind), &
										' IS NOT CONSTANT IN TIME FOR SPECIE= ', s, &
										', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND MASTER', &
										' TIME-STEP= ', nn, ' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' &
										// achar(27) // '[0m.'
								end if
							end if

							! ----------------------------------------------------

						end do

						! ----------------------------------------------------

					end if

					! ----------------------------------------------------

				end if

				! ----------------------------------------------------

			end if
		end do

		! ----------------------------------------------------

		! SET INPUT AMBIPOLAR FIELD FROM PREVIOUS SIMULATION:

		! ----------------------------------------------------

		if ((SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1) == 0) .and. (rank == 0) &
			.and. ((SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 1) .or. (dNTe /= 0d0))) then
			do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
				if ((n == 1) .and. (nn == 1)) then
					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
						SpecieT(s)%FluxTubeT(f)%EAInertialRT(nn, Qind)= &
							SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%EAInertialInputT(1)
						SpecieT(s)%FluxTubeT(f)%EAPressureRT(nn, Qind)= &
							SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%EAPressureInputT(1)
						SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind)= &
							SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%EAmagInputT(1)
					end do
				end if
			end do
		end if

		if ((SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1) == 0) .and. (rank == 0) &
			.and. (SpecieT(1)%FluxTubeT(1)%EAMBSELFCONSISTflagT(1) == 0) .and. (dNTe == 0d0)) then
			do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
				if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
					(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1))))) then
					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
						SpecieT(s)%FluxTubeT(f)%EAInertialRT(nn, Qind)= &
							SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%EAInertialInputT(1)
						SpecieT(s)%FluxTubeT(f)%EAPressureRT(nn, Qind)= &
							SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%EAPressureInputT(1)
						SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind)= &
							SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%EAmagInputT(1)
					end do
				end if
			end do
		end if

		! ----------------------------------------------------

		! SET OUTPUT AMBIPOLAR FIELD FOR NEXT SIMULATION:

		! ----------------------------------------------------

		if ((SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1) == 1) .and. (rank == 0)) then
			do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
				if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
					(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1))))) then
					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
						SpecieT(s)%FluxTubeT(f)%EAInertialOutputRT(Qind)= SpecieT(s)%FluxTubeT(f)%EAInertialRT(nn, Qind)
						SpecieT(s)%FluxTubeT(f)%EAPressureOutputRT(Qind)= SpecieT(s)%FluxTubeT(f)%EAPressureRT(nn, Qind)
						SpecieT(s)%FluxTubeT(f)%EAmagOutputRT(Qind)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind)
					end do
				end if
			end do
		end if

		! ----------------------------------------------------

		! BROADCAST AMBIPOLAR ELECTRIC FIELD MAGNITUDES ACROSS ALL MPI RANKS:

		! ----------------------------------------------------

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
				(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1))))) then
				if (SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 1) then

					call mpi_barrier(MPI_COMM_WORLD, ierr)
					call mpi_bcast(SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, :), &
						((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1), &
						MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

				end if
			end if
		end do

		! ----------------------------------------------------

		! COMPUTE ALL PARTICLE AMBIPOLAR ELECTRIC FORCE ACCELERATION MAGNITUDES ON MASTER TIME-STEPS:

		! ----------------------------------------------------

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
				(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1))))) then

				! ----------------------------------------------------

				if (SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 1) then

					! ----------------------------------------------------

					! Compute Ambipolar E field for Southern Magnetic Hemisphere
					if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then
						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
							do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1

								! ----------------------------------------------------

								if (nn == 1) then
									if (Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) then
										if ((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) <= &
		                  SpecieT(s)%FluxTubeT(f)%q0T(j)) .and. &
		                  (SpecieT(s)%FluxTubeT(f)%q0T(j) <= &
		                  SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind))) then

											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
												EAmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind))
											end if
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
												EAmagInterp(1)= (SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind))
											end if

											if (EAmagInterp(1) == 0d0) then
												AEAmagN(j)= 0d0
											else
			                  if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
													AEAmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
												if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
													AEAmagN(j)= ((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEAmagN(j) == 0) .and. (EAmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' MASTER TIME-STEP HAS ZERO AEAmagN= ', AEAmagN(j), &
													', EAmagInterp= ', EAmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', MASTER TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								if (nn == 1) then
									if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. (Qind /= SpecieT(s)%FluxTubeT(f)%NqUBT(1))) then
										if ((SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind- 1) < &
		                  SpecieT(s)%FluxTubeT(f)%q0T(j)) .and. &
		                  (SpecieT(s)%FluxTubeT(f)%q0T(j) <= &
		                  SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind))) then

											xLinInterp(1)= SpecieT(s)%FluxTubeT(f)%q0T(j)
											xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind- 1)
											xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind)
											yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind- 1)
											yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind)

											call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
												yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
												EAmagInterp(1)= abs(yLinInterp(1))
											end if
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
												EAmagInterp(1)= (yLinInterp(1))
											end if

											if (EAmagInterp(1) == 0d0) then
												AEAmagN(j)= 0d0
											else
			                  if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
													AEAmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
												if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
													AEAmagN(j)= ((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEAmagN(j) == 0) .and. (EAmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' MASTER TIME-STEP HAS ZERO AEAmagN= ', AEAmagN(j), &
													', EAmagInterp= ', EAmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', MASTER TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								if (nn == 1) then
									if (Qind == SpecieT(s)%FluxTubeT(f)%NqUBT(1)) then
										if ((SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind- 1) < &
											SpecieT(s)%FluxTubeT(f)%q0T(j)) .and. &
											(SpecieT(s)%FluxTubeT(f)%q0T(j) < &
											SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind))) then

											xLinInterp(1)= SpecieT(s)%FluxTubeT(f)%q0T(j)
											xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind- 1)
											xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind)
											yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind- 1)
											yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind)

											call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
												yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
												EAmagInterp(1)= abs(yLinInterp(1))
											end if
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
												EAmagInterp(1)= (yLinInterp(1))
											end if

											if (EAmagInterp(1) == 0d0) then
												AEAmagN(j)= 0d0
											else
												if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
													AEAmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
												if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
													AEAmagN(j)= ((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEAmagN(j) == 0) .and. (EAmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' MASTER TIME-STEP HAS ZERO AEAmagN= ', AEAmagN(j), &
													', EAmagInterp= ', EAmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', MASTER TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if

										if ((SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind) <= &
		                  SpecieT(s)%FluxTubeT(f)%q0T(j)) .and. &
		                  (SpecieT(s)%FluxTubeT(f)%q0T(j) <= &
		                  SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind))) then

											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
												EAmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind))
											end if
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
												EAmagInterp(1)= (SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind))
											end if

											if (EAmagInterp(1) == 0d0) then
												AEAmagN(j)= 0d0
											else
			                  if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
													AEAmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
												if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
													AEAmagN(j)= ((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEAmagN(j) == 0) .and. (EAmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' MASTER TIME-STEP HAS ZERO AEAmagN= ', AEAmagN(j), &
													', EAmagInterp= ', EAmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', MASTER TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								! ----------------------------------------------------

								if (nn /= 1) then
									if (Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) then
										if ((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) <= &
		                  qk4(j)) .and. (qk4(j) <= &
		                  SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind))) then

											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
												EAmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind))
											end if
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
												EAmagInterp(1)= (SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind))
											end if

											if (EAmagInterp(1) == 0d0) then
												AEAmagN(j)= 0d0
											else
			                  if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
													AEAmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
												if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
													AEAmagN(j)= ((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEAmagN(j) == 0) .and. (EAmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' MASTER TIME-STEP HAS ZERO AEAmagN= ', AEAmagN(j), &
													', EAmagInterp= ', EAmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', MASTER TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								if (nn /= 1) then
									if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. (Qind /= SpecieT(s)%FluxTubeT(f)%NqUBT(1))) then
										if ((SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind- 1) < &
		                  qk4(j)) .and. (qk4(j) <= &
		                  SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind))) then

											xLinInterp(1)= qk4(j)
											xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind- 1)
											xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind)
											yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind- 1)
											yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind)

											call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
												yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
												EAmagInterp(1)= abs(yLinInterp(1))
											end if
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
												EAmagInterp(1)= (yLinInterp(1))
											end if

											if (EAmagInterp(1) == 0d0) then
												AEAmagN(j)= 0d0
											else
			                  if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
													AEAmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
												if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
													AEAmagN(j)= ((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEAmagN(j) == 0) .and. (EAmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' MASTER TIME-STEP HAS ZERO AEAmagN= ', AEAmagN(j), &
													', EAmagInterp= ', EAmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', MASTER TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								if (nn /= 1) then
									if (Qind == SpecieT(s)%FluxTubeT(f)%NqUBT(1)) then
										if ((SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind- 1) < &
											qk4(j)) .and. (qk4(j) < &
											SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind))) then

											xLinInterp(1)= qk4(j)
											xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind- 1)
											xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind)
											yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind- 1)
											yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind)

											call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
												yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
												EAmagInterp(1)= abs(yLinInterp(1))
											end if
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
												EAmagInterp(1)= (yLinInterp(1))
											end if

											if (EAmagInterp(1) == 0d0) then
												AEAmagN(j)= 0d0
											else
												if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
													AEAmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
												if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
													AEAmagN(j)= ((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEAmagN(j) == 0) .and. (EAmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' MASTER TIME-STEP HAS ZERO AEAmagN= ', AEAmagN(j), &
													', EAmagInterp= ', EAmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', MASTER TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if

										if ((SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind) <= &
		                  qk4(j)) .and. (qk4(j) <= &
		                  SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind))) then

											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
												EAmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind))
											end if
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
												EAmagInterp(1)= (SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind))
											end if

											if (EAmagInterp(1) == 0d0) then
												AEAmagN(j)= 0d0
											else
			                  if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
													AEAmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
												if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
													AEAmagN(j)= ((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEAmagN(j) == 0) .and. (EAmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' MASTER TIME-STEP HAS ZERO AEAmagN= ', AEAmagN(j), &
													', EAmagInterp= ', EAmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', MASTER TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								! ----------------------------------------------------

							end do
						end do
					end if

					! ----------------------------------------------------

					! Compute Ambipolar E field for Northern Magnetic Hemisphere
					if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then
						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
							do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1

								! ----------------------------------------------------

								if (nn == 1) then
									if (Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) then
										if ((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) >= &
		                  SpecieT(s)%FluxTubeT(f)%q0T(j)) .and. &
		                  (SpecieT(s)%FluxTubeT(f)%q0T(j) >= &
		                  SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind))) then

											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
												EAmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind))
											end if
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
												EAmagInterp(1)= (SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind))
											end if

											if (EAmagInterp(1) == 0d0) then
												AEAmagN(j)= 0d0
											else
			                  if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
													AEAmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
												if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
													AEAmagN(j)= ((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEAmagN(j) == 0) .and. (EAmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' MASTER TIME-STEP HAS ZERO AEAmagN= ', AEAmagN(j), &
													', EAmagInterp= ', EAmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', MASTER TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								if (nn == 1) then
									if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. (Qind /= SpecieT(s)%FluxTubeT(f)%NqUBT(1))) then
										if ((SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind- 1) > &
		                  SpecieT(s)%FluxTubeT(f)%q0T(j)) .and. &
		                  (SpecieT(s)%FluxTubeT(f)%q0T(j) >= &
		                  SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind))) then

											xLinInterp(1)= SpecieT(s)%FluxTubeT(f)%q0T(j)
											xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind- 1)
											xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind)
											yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind- 1)
											yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind)

											call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
												yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
												EAmagInterp(1)= abs(yLinInterp(1))
											end if
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
												EAmagInterp(1)= (yLinInterp(1))
											end if

											if (EAmagInterp(1) == 0d0) then
												AEAmagN(j)= 0d0
											else
			                  if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
													AEAmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
												if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
													AEAmagN(j)= ((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEAmagN(j) == 0) .and. (EAmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' MASTER TIME-STEP HAS ZERO AEAmagN= ', AEAmagN(j), &
													', EAmagInterp= ', EAmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', MASTER TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								if (nn == 1) then
									if (Qind == SpecieT(s)%FluxTubeT(f)%NqUBT(1)) then
										if ((SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind- 1) > &
											SpecieT(s)%FluxTubeT(f)%q0T(j)) .and. &
											(SpecieT(s)%FluxTubeT(f)%q0T(j) > &
											SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind))) then

											xLinInterp(1)= SpecieT(s)%FluxTubeT(f)%q0T(j)
											xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind- 1)
											xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind)
											yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind- 1)
											yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind)

											call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
												yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
												EAmagInterp(1)= abs(yLinInterp(1))
											end if
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
												EAmagInterp(1)= (yLinInterp(1))
											end if

											if (EAmagInterp(1) == 0d0) then
												AEAmagN(j)= 0d0
											else
												if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
													AEAmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
												if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
													AEAmagN(j)= ((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEAmagN(j) == 0) .and. (EAmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' MASTER TIME-STEP HAS ZERO AEAmagN= ', AEAmagN(j), &
													', EAmagInterp= ', EAmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', MASTER TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if

										if ((SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind) >= &
		                  SpecieT(s)%FluxTubeT(f)%q0T(j)) .and. &
		                  (SpecieT(s)%FluxTubeT(f)%q0T(j) >= &
		                  SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind))) then

											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
												EAmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind))
											end if
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
												EAmagInterp(1)= (SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind))
											end if

											if (EAmagInterp(1) == 0d0) then
												AEAmagN(j)= 0d0
											else
			                  if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
													AEAmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
												if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
													AEAmagN(j)= ((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEAmagN(j) == 0) .and. (EAmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' MASTER TIME-STEP HAS ZERO AEAmagN= ', AEAmagN(j), &
													', EAmagInterp= ', EAmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', MASTER TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								! ----------------------------------------------------

								if (nn /= 1) then
									if (Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) then
										if ((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) >= &
		                  qk4(j)) .and. (qk4(j) >= &
		                  SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind))) then

											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
												EAmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind))
											end if
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
												EAmagInterp(1)= (SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind))
											end if

											if (EAmagInterp(1) == 0d0) then
												AEAmagN(j)= 0d0
											else
			                  if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
													AEAmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
												if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
													AEAmagN(j)= ((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEAmagN(j) == 0) .and. (EAmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' MASTER TIME-STEP HAS ZERO AEAmagN= ', AEAmagN(j), &
													', EAmagInterp= ', EAmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', MASTER TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								if (nn /= 1) then
									if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. (Qind /= SpecieT(s)%FluxTubeT(f)%NqUBT(1))) then
										if ((SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind- 1) > &
		                  qk4(j)) .and. (qk4(j) >= &
		                  SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind))) then

											xLinInterp(1)= qk4(j)
											xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind- 1)
											xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind)
											yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind- 1)
											yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind)

											call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
												yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
												EAmagInterp(1)= abs(yLinInterp(1))
											end if
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
												EAmagInterp(1)= (yLinInterp(1))
											end if

											if (EAmagInterp(1) == 0d0) then
												AEAmagN(j)= 0d0
											else
			                  if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
													AEAmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
												if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
													AEAmagN(j)= ((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEAmagN(j) == 0) .and. (EAmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' MASTER TIME-STEP HAS ZERO AEAmagN= ', AEAmagN(j), &
													', EAmagInterp= ', EAmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', MASTER TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								if (nn /= 1) then
									if (Qind == SpecieT(s)%FluxTubeT(f)%NqUBT(1)) then
										if ((SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind- 1) > &
											qk4(j)) .and. (qk4(j) > &
											SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind))) then

											xLinInterp(1)= qk4(j)
											xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind- 1)
											xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind)
											yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind- 1)
											yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind)

											call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
												yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
												EAmagInterp(1)= abs(yLinInterp(1))
											end if
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
												EAmagInterp(1)= (yLinInterp(1))
											end if

											if (EAmagInterp(1) == 0d0) then
												AEAmagN(j)= 0d0
											else
												if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
													AEAmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
												if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
													AEAmagN(j)= ((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEAmagN(j) == 0) .and. (EAmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' MASTER TIME-STEP HAS ZERO AEAmagN= ', AEAmagN(j), &
													', EAmagInterp= ', EAmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', MASTER TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if

										if ((SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind) >= &
		                  qk4(j)) .and. (qk4(j) >= &
		                  SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind))) then

											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
												EAmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind))
											end if
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
												EAmagInterp(1)= (SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind))
											end if

											if (EAmagInterp(1) == 0d0) then
												AEAmagN(j)= 0d0
											else
												if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
				                  AEAmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
												if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
													AEAmagN(j)= ((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
												end if
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEAmagN(j) == 0) .and. (EAmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' MASTER TIME-STEP HAS ZERO AEAmagN= ', AEAmagN(j), &
													', EAmagInterp= ', EAmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', MASTER TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								! ----------------------------------------------------

							end do
						end do
					end if

					! ----------------------------------------------------

					! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

					do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
						if ((isnan(real(AEAmagN(j))) .eqv. .true.) .or. &
							(size(AEAmagN(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAmagN= ', &
								AEAmagN(j), ' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
								' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
						end if
					end do

					! ----------------------------------------------------

				else if (SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 0) then

					AEAmagN(:)= 0d0

				end if

				! ----------------------------------------------------

			end if
		end do

		! ----------------------------------------------------

		! COMPUTE ALL PARTICLE AMBIPOLAR ELECTRIC FORCE ACCELERATION MAGNITUDES ON NON-MASTER TIME-STEPS:

		! ----------------------------------------------------

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if ((nn /= 1d0) .and. (((nn == 2d0) .and. ((n > sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 2))+ 1d0) .and. &
				(n < sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1))))) .or. &
				((nn > 2d0) .and. ((n >= sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 2))+ 1d0) .and. &
				(n < sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1))))))) then

				! ----------------------------------------------------

				if (SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 1) then

					! ----------------------------------------------------

					! Compute Ambipolar E field for Southern Magnetic Hemisphere
					if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn- 1, 1) <= 0) then
						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
							do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1

								! ----------------------------------------------------

								if (Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) then
									if ((SpecieT(s)%FluxTubeT(f)%qGLGT(nn- 1, Qind) <= &
	                  qk4(j)) .and. (qk4(j) <= &
	                  SpecieT(s)%FluxTubeT(f)%qGCGT(nn- 1, Qind))) then

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
											EAmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EAmagRT(nn- 1, Qind))
										end if
										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
											EAmagInterp(1)= (SpecieT(s)%FluxTubeT(f)%EAmagRT(nn- 1, Qind))
										end if

										if (EAmagInterp(1) == 0d0) then
											AEAmagN(j)= 0d0
										else
		                  if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
												AEAmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
											end if
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
												AEAmagN(j)= ((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
											end if
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

										if ((AEAmagN(j) == 0) .and. (EAmagInterp(1) /= 0d0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NON-MASTER TIME-STEP HAS ZERO AEAmagN= ', AEAmagN(j), &
												', EAmagInterp= ', EAmagInterp(1), &
												' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
												f, ', Q GRID CELL= ', Qind, ', TIME-STEP= ', n, &
												', AND PARTICLE= ', j, &
												' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end if
								end if

								if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. (Qind /= SpecieT(s)%FluxTubeT(f)%NqUBT(1))) then
									if ((SpecieT(s)%FluxTubeT(f)%qGCGT(nn- 1, Qind- 1) < &
	                  qk4(j)) .and. (qk4(j) <= &
	                  SpecieT(s)%FluxTubeT(f)%qGCGT(nn- 1, Qind))) then

										xLinInterp(1)= qk4(j)
										xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%qGCGT(nn- 1, Qind- 1)
										xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%qGCGT(nn- 1, Qind)
										yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn- 1, Qind- 1)
										yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn- 1, Qind)

										call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
											yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
											EAmagInterp(1)= abs(yLinInterp(1))
										end if
										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
											EAmagInterp(1)= (yLinInterp(1))
										end if

										if (EAmagInterp(1) == 0d0) then
											AEAmagN(j)= 0d0
										else
		                  if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
												AEAmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
											end if
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
												AEAmagN(j)= ((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
											end if
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

										if ((AEAmagN(j) == 0) .and. (EAmagInterp(1) /= 0d0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NON-MASTER TIME-STEP HAS ZERO AEAmagN= ', AEAmagN(j), &
												', EAmagInterp= ', EAmagInterp(1), &
												' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
												f, ', Q GRID CELL= ', Qind, ', TIME-STEP= ', n, &
												', AND PARTICLE= ', j, &
												' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end if
								end if

								if (Qind == SpecieT(s)%FluxTubeT(f)%NqUBT(1)) then
									if ((SpecieT(s)%FluxTubeT(f)%qGCGT(nn- 1, Qind- 1) < &
										qk4(j)) .and. (qk4(j) < &
										SpecieT(s)%FluxTubeT(f)%qGCGT(nn- 1, Qind))) then

										xLinInterp(1)= qk4(j)
										xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%qGCGT(nn- 1, Qind- 1)
										xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%qGCGT(nn- 1, Qind)
										yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn- 1, Qind- 1)
										yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn- 1, Qind)

										call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
											yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
											EAmagInterp(1)= abs(yLinInterp(1))
										end if
										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
											EAmagInterp(1)= (yLinInterp(1))
										end if

										if (EAmagInterp(1) == 0d0) then
											AEAmagN(j)= 0d0
										else
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
												AEAmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
											end if
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
												AEAmagN(j)= ((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
											end if
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

										if ((AEAmagN(j) == 0) .and. (EAmagInterp(1) /= 0d0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NON-MASTER TIME-STEP HAS ZERO AEAmagN= ', AEAmagN(j), &
												', EAmagInterp= ', EAmagInterp(1), &
												' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
												f, ', Q GRID CELL= ', Qind, ', TIME-STEP= ', n, &
												', AND PARTICLE= ', j, &
												' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end if

									if ((SpecieT(s)%FluxTubeT(f)%qGCGT(nn- 1, Qind) <= &
	                  qk4(j)) .and. (qk4(j) <= &
	                  SpecieT(s)%FluxTubeT(f)%qGHGT(nn- 1, Qind))) then

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
											EAmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EAmagRT(nn- 1, Qind))
										end if
										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
											EAmagInterp(1)= (SpecieT(s)%FluxTubeT(f)%EAmagRT(nn- 1, Qind))
										end if

										if (EAmagInterp(1) == 0d0) then
											AEAmagN(j)= 0d0
										else
		                  if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
												AEAmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
											end if
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
												AEAmagN(j)= ((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
											end if
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

										if ((AEAmagN(j) == 0) .and. (EAmagInterp(1) /= 0d0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NON-MASTER TIME-STEP HAS ZERO AEAmagN= ', AEAmagN(j), &
												', EAmagInterp= ', EAmagInterp(1), &
												' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
												f, ', Q GRID CELL= ', Qind, ', TIME-STEP= ', n, &
												', AND PARTICLE= ', j, &
												' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end if
								end if

								! ----------------------------------------------------

							end do
						end do
					end if

					! ----------------------------------------------------

					! Compute Ambipolar E field for Northern Magnetic Hemisphere
					if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn- 1, 1) > 0) then
						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
							do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1

								! ----------------------------------------------------

								if (Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) then
									if ((SpecieT(s)%FluxTubeT(f)%qGLGT(nn- 1, Qind) >= &
	                  qk4(j)) .and. (qk4(j) >= &
	                  SpecieT(s)%FluxTubeT(f)%qGCGT(nn- 1, Qind))) then

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
											EAmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EAmagRT(nn- 1, Qind))
										end if
										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
											EAmagInterp(1)= (SpecieT(s)%FluxTubeT(f)%EAmagRT(nn- 1, Qind))
										end if

										if (EAmagInterp(1) == 0d0) then
											AEAmagN(j)= 0d0
										else
		                  if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
												AEAmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
											end if
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
												AEAmagN(j)= ((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
											end if
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

										if ((AEAmagN(j) == 0) .and. (EAmagInterp(1) /= 0d0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NON-MASTER TIME-STEP HAS ZERO AEAmagN= ', AEAmagN(j), &
												', EAmagInterp= ', EAmagInterp(1), &
												' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
												f, ', Q GRID CELL= ', Qind, ', TIME-STEP= ', n, &
												', AND PARTICLE= ', j, &
												' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end if
								end if

								if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. (Qind /= SpecieT(s)%FluxTubeT(f)%NqUBT(1))) then
									if ((SpecieT(s)%FluxTubeT(f)%qGCGT(nn- 1, Qind- 1) > &
	                  qk4(j)) .and. (qk4(j) >= &
	                  SpecieT(s)%FluxTubeT(f)%qGCGT(nn- 1, Qind))) then

										xLinInterp(1)= qk4(j)
										xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%qGCGT(nn- 1, Qind- 1)
										xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%qGCGT(nn- 1, Qind)
										yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn- 1, Qind- 1)
										yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn- 1, Qind)

										call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
											yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
											EAmagInterp(1)= abs(yLinInterp(1))
										end if
										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
											EAmagInterp(1)= (yLinInterp(1))
										end if

										if (EAmagInterp(1) == 0d0) then
											AEAmagN(j)= 0d0
										else
		                  if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
												AEAmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
											end if
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
												AEAmagN(j)= ((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
											end if
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

										if ((AEAmagN(j) == 0) .and. (EAmagInterp(1) /= 0d0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NON-MASTER TIME-STEP HAS ZERO AEAmagN= ', AEAmagN(j), &
												', EAmagInterp= ', EAmagInterp(1), &
												' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
												f, ', Q GRID CELL= ', Qind, ', TIME-STEP= ', n, &
												', AND PARTICLE= ', j, &
												' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end if
								end if

								if (Qind == SpecieT(s)%FluxTubeT(f)%NqUBT(1)) then
									if ((SpecieT(s)%FluxTubeT(f)%qGCGT(nn- 1, Qind- 1) > &
										qk4(j)) .and. (qk4(j) > &
										SpecieT(s)%FluxTubeT(f)%qGCGT(nn- 1, Qind))) then

										xLinInterp(1)= qk4(j)
										xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%qGCGT(nn- 1, Qind- 1)
										xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%qGCGT(nn- 1, Qind)
										yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn- 1, Qind- 1)
										yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EAmagRT(nn- 1, Qind)

										call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
											yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
											EAmagInterp(1)= abs(yLinInterp(1))
										end if
										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
											EAmagInterp(1)= (yLinInterp(1))
										end if

										if (EAmagInterp(1) == 0d0) then
											AEAmagN(j)= 0d0
										else
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
												AEAmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
											end if
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
												AEAmagN(j)= ((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
											end if
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

										if ((AEAmagN(j) == 0) .and. (EAmagInterp(1) /= 0d0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NON-MASTER TIME-STEP HAS ZERO AEAmagN= ', AEAmagN(j), &
												', EAmagInterp= ', EAmagInterp(1), &
												' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
												f, ', Q GRID CELL= ', Qind, ', TIME-STEP= ', n, &
												', AND PARTICLE= ', j, &
												' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end if

									if ((SpecieT(s)%FluxTubeT(f)%qGCGT(nn- 1, Qind) >= &
	                  qk4(j)) .and. (qk4(j) >= &
	                  SpecieT(s)%FluxTubeT(f)%qGHGT(nn- 1, Qind))) then

										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
											EAmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EAmagRT(nn- 1, Qind))
										end if
										if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
											EAmagInterp(1)= (SpecieT(s)%FluxTubeT(f)%EAmagRT(nn- 1, Qind))
										end if

										if (EAmagInterp(1) == 0d0) then
											AEAmagN(j)= 0d0
										else
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
			                  AEAmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
											end if
											if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
												AEAmagN(j)= ((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EAmagInterp(1))
											end if
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

										if ((AEAmagN(j) == 0) .and. (EAmagInterp(1) /= 0d0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NON-MASTER TIME-STEP HAS ZERO AEAmagN= ', AEAmagN(j), &
												', EAmagInterp= ', EAmagInterp(1), &
												' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
												f, ', Q GRID CELL= ', Qind, ', TIME-STEP= ', n, &
												', AND PARTICLE= ', j, &
												' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end if
								end if

								! ----------------------------------------------------

							end do
						end do
					end if

					! ----------------------------------------------------

					! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

					do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
						if ((isnan(real(AEAmagN(j))) .eqv. .true.) .or. &
							(size(AEAmagN(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAmagN= ', &
								AEAmagN(j), ' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
								' IN AMBIPOLAR ELECTRIC FIELD SUBROUTINE' // achar(27) // '[0m.'
						end if
					end do

					! ----------------------------------------------------

				else if (SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 0) then

					AEAmagN(:)= 0d0

				end if

				! ----------------------------------------------------

			end if
		end do

		! ----------------------------------------------------

		!if (SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 1) then
		! do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
		!	 if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
		!		 (n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1))))) then
		!		 do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
		!			 if (rank == 0) then
		!				 if (SpecieT(s)%FluxTubeT(f)%M0phRT(nn- 1, Qind) /= 0) then
		!					 write(*, *) 'EAmagRT s, f, nn, Qind= ', s, f, nn, Qind, &
		!						 SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, Qind)
		!				 end if
		!			 end if
		!		 end do
		!	 end if
		! end do
			!do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			!	if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
			!		 (n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1))))) then
					 !write(*, *) 'AEAmagN s, f= ', s, f, AEAmagN(:)
			!	end if
			!end do
		!end if

		! ----------------------------------------------------

	end subroutine AmbipolarEfieldSub

end module AmbipolarEfield
