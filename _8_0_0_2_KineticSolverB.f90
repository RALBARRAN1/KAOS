module KineticSolverB

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.0.0.2 3D KINETIC SOLVER B:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use RK4DipolePolynomialSolver
use MBVelocityDistribution
use KineticRK4Update
use ParticleCounts
use DistributionFunctions
use ZerothIonMoment
use FirstPerpIonMoment
use FirstParIonMoment
use SecondIonMoment
use SecondPerpIonMoment
use SecondParIonMoment
use AmbipolarEfield
use Gravfield
use PotentialStructure
use ZerothENAMoment
use FirstPENAMoment
use FirstQENAMoment
use FirstPHIENAMoment
use SecondENAMoment
use SecondPENAMoment
use SecondQENAMoment
use SecondPHIENAMoment
use DataExport2

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

	subroutine KineticSolverBSub

		! ----------------------------------------------------

		if (n /= 1) then

			! ----------------------------------------------------

			! UPDATE TIME PARAMETER ON STATISTICAL TIME-STEP:

			nnloopKS0: do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
				if (nn > SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1) then
					exit nnloopKS0
				end if
				if (n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1)) then
					SpecieT(s)%FluxTubeT(f)%TimeT(nn)= Time(1)
				end if
			end do nnloopKS0

			! Set to nnFlag= 1 for injection time-step (0 otherwise)
			do nn= 1, SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1, 1
				if (n /= (nn- 1)*SpecieT(s)%FluxTubeT(f)%Q0ndatfacT(1)) then
					nnFlag(1)= 0d0
				end if
				if (n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%Q0ndatfacT(1)) then
					nnFlag(1)= 1d0
				end if
			end do

			! ----------------------------------------------------

			! UPDATE ALL EXISTING PARTICLES ON INJECTION TIME-STEP:

			! ----------------------------------------------------

			! Update variables on injection time-steps

			nnloopKS1: do nn= 1, SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1, 1
				if (nn > SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1) then
					exit nnloopKS1
				end if
				if (n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%Q0ndatfacT(1)) then

 					nnFlag(1)= 1d0

					! ----------------------------------------------------

					! DIAGNOSTIC FLAGS FOR PROPER INJECTION DENSITY:

					if (SpecieT(s)%FluxTubeT(f)%NsnT(nn- 1) /= (NsTK(1)+ dNsTK1(1)- dNsTK2(1)- dNsTK3(1))) then
					 write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						 ' INCONSISTENT INJECTION DENSITY= ', SpecieT(s)%FluxTubeT(f)%NsnT(nn- 1), &
						 (NsTK(1)+ dNsTK1(1)- dNsTK2(1)- dNsTK3(1)), ' FOR SPECIE= ', s, &
						 ', FLUX TUBE= ', f, ', AND INJECTION TIME-STEP= ', nn, &
						 ' IN KINETIC SOLVER B SUBROUTINE' &
						 // achar(27) // '[0m.'
					end if

					if (nnFlag(1) /= 1) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						 ' INCORRECT nnFlag= ', nnFlag(1), ' FOR SPECIE= ', s, &
						 ', FLUX TUBE= ', f, ', AND INJECTION TIME-STEP= ', nn, &
						 ' IN KINETIC SOLVER B SUBROUTINE' &
						 // achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

					! Update all current particle variables
					do j= 1, (NsTK(1)- dNsTK2(1)- dNsTK3(1)), 1

						! Reset Ambipolar Electric Field Accelerations
						if (ENAflag(j) .eqv. .false.) then
							AEAmag(j)= AEAmagN(j)
							AGmag(j)= AGmagN(j)
							AEPmag(j)= AEPmagN(j)
						end if

						! ----------------------------------------------------

						! Update particle phase-space coordinates
						x(j)= xN(j)
						y(j)= yN(j)
						z(j)= zN(j)
						if (ENAflag(j) .eqv. .false.) then
							Vperp1(j)= Vperp1N(j)
							Vperp2(j)= Vperp2N(j)
							Vperp(j)= VperpN(j)
						else if ((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
							(ENAflag(j) .eqv. .true.)) then
							Vperp1(j)= 0d0
							Vperp2(j)= 0d0
							Vperp(j)= 0d0
						end if
						Vx(j)= VxN(j)
						Vy(j)= VyN(j)
						Vz(j)= VzN(j)

					end do

				end if
			end do nnloopKS1

			if (nnFlag(1) == 1) then

				! ----------------------------------------------------
				! Perform iterative updates of particle positions and velocities
				jloopKS1: do j= 1, (NsTK(1)- dNsTK2(1)- dNsTK3(1)), 1

					! ----------------------------------------------------

					! DIAGNOSTIC FLAG FOR CORRECT ION/ENA PARTICLE TYPE:

					if ((ENAflag(j) .eqv. .true.) .and. (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 0)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ENA PRODUCED WITH', &
							' CHARGE-EXCHANGE TURNED OFF FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
							' IN KINETIC SOLVER B SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

					call KineticRK4UpdateSub

					! ----------------------------------------------------

					if (j < (NsTK(1)- dNsTK2(1)- dNsTK3(1))) then
						cycle jloopKS1
					end if
					if (j >= (NsTK(1)- dNsTK2(1)- dNsTK3(1))) then
						exit jloopKS1
					end if

				end do jloopKS1

				! ----------------------------------------------------

			end if

			! ----------------------------------------------------

			! INITIALIZE LOWER BOUNDARY INJECTED IONS ON INJECTION TIME-STEP:

			! ----------------------------------------------------

			nnloopKS2: do nn= 1, SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1, 1
				if (nn > SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1) then
					exit nnloopKS2
				end if
				if (n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%Q0ndatfacT(1)) then

					! ----------------------------------------------------
					! Initialize all LB injected particle variables
					do j= (NsTK(1)- dNsTK2(1)- dNsTK3(1))+ 1, NsTK(1)- dNsTK3(1), 1

						! ----------------------------------------------------

						! Inject ions at lower altitude boundary
						ENAflag(j)= .false.
						ENAflagN0ind(j)= n
						AEAmag(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*abs(SpecieT(s)%FluxTubeT(f)%EAmagRT(nn- 1, NqLB(1))))
						AGmag(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*abs(SpecieT(s)%FluxTubeT(f)%EGmagRT(NqLB(1))))
						AEPmag(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*abs(SpecieT(s)%FluxTubeT(f)%EPmagRT(NqLB(1))))

						! ----------------------------------------------------

						! Set initial dipole coordinates
						qNp(1)= SpecieT(s)%FluxTubeT(f)%QCellT(NqLB(1))%qGCT(1)
						pNp(1)= SpecieT(s)%FluxTubeT(f)%QCellT(1)%pGCT(1)

						! ----------------------------------------------------

						if (SpecieT(s)%FluxTubeT(f)%QCellT(1)%qGLT(1) <= 0) then ! S. Magnetic Hemisphere
							QloopKSB1: do Qind= NqLB(1), NqUB(1), 1
								if ((Qind == NqLB(1)) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGLT(1) <= qNp(1)) &
									.and. (qNp(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGHT(1))) then
									Qindk1(j)= Qind

									exit QloopKSB1

								else if ((Qind /= NqLB(1)) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGLT(1) < qNp(1)) &
									.and. (qNp(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGHT(1))) then
									Qindk1(j)= Qind

									exit QloopKSB1

								end if
							end do QloopKSB1

							if (SpecieT(s)%FluxTubeT(f)%QCellT(NqLB(1))%qGLT(1) > qNp(1)) then
								! SMH Lower boundary escape
								Qindk1(j)= 0d0
							else if (SpecieT(s)%FluxTubeT(f)%QCellT(NqUB(1))%qGHT(1) < qNp(1)) then
								! SMH Upper boundary escape
								Qindk1(j)= -1d0
							end if
						end if

						! ----------------------------------------------------

						if (SpecieT(s)%FluxTubeT(f)%QCellT(1)%qGLT(1) > 0) then ! N. Magnetic Hemisphere
							QloopKSB2: do Qind= NqLB(1), NqUB(1), 1
								if ((Qind == NqLB(1)) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGLT(1) >= qNp(1)) &
									.and. (qNp(1) >= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGHT(1))) then
									Qindk1(j)= Qind

									exit QloopKSB2

								else if ((Qind /= NqLB(1)) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGLT(1) > qNp(1)) &
									.and. (qNp(1) >= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGHT(1))) then
									Qindk1(j)= Qind

									exit QloopKSB2

								end if
							end do QloopKSB2

							if (SpecieT(s)%FluxTubeT(f)%QCellT(NqLB(1))%qGLT(1) < qNp(1)) then
								! NMH Lower boundary escape
								Qindk1(j)= 0d0
							else if (SpecieT(s)%FluxTubeT(f)%QCellT(NqUB(1))%qGHT(1) > qNp(1)) then
								! NMH Upper boundary escape
								Qindk1(j)= -1d0
							end if
						end if

						! ----------------------------------------------------

						! DIAGNOSTIC FLAG FOR PARTICLE POSITIONS WITHIN BOUNDS:

						if (SpecieT(s)%FluxTubeT(f)%QCellT(1)%qGLT(1) <= 0) then ! S. Magnetic Hemisphere)
							if ((qNp(1) < SpecieT(s)%FluxTubeT(f)%QCellT(NqLB(1))%qGLT(1)) .or. &
								(qNp(1) > SpecieT(s)%FluxTubeT(f)%QCellT(NqLB(1))%qGHT(1))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qniICT IS', &
									' OUTSIDE LOWER BOUNDARY BOUNDS IN SOUTH MAGNETIC HEMISPHERE FOR SPECIE= ', s, ', FLUX TUBE= ', &
									f, ', AND Qind= ', Qind, ' IN KINETIC SOLVER B SUBROUTINE' // achar(27) // '[0m.'
							end if
						end if
						if (SpecieT(s)%FluxTubeT(f)%QCellT(1)%qGLT(1) > 0) then ! N. Magnetic Hemisphere
							if ((qNp(1) > SpecieT(s)%FluxTubeT(f)%QCellT(NqLB(1))%qGLT(1)) .or. &
								(qNp(1) < SpecieT(s)%FluxTubeT(f)%QCellT(NqLB(1))%qGHT(1))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qNp IS', &
									' OUTSIDE LOWER BOUNDARY BOUNDS IN NORTH MAGNETIC HEMISPHERE FOR SPECIE= ', s, ', FLUX TUBE= ', &
									f, ', AND Qind= ', Qind, ' IN KINETIC SOLVER B SUBROUTINE' // achar(27) // '[0m.'
							end if
						end if

						! ----------------------------------------------------

						call RK4DipolePolynomialSolverSub

						! ----------------------------------------------------

						xN(j)= xfinalRK4(1)
						yN(j)= yfinalRK4(1)
						zN(j)= zfinalRK4(1)

						! ----------------------------------------------------

						! Inject ions with MB distributed velocities
						thetaMB(1)= thetafinalRK4(1)
						phiMB(1)= phifinalRK4(1)
						ellMB(1)= 1d0+ 3d0*(cos(thetaMB(1)))**2d0
						TsMB(1)= SpecieT(s)%FluxTubeT(f)%QCellT(NqLB(1))%TsT(1)
						LBflag(1)= 1d0

						! ----------------------------------------------------

						call MBVelocityDistributionSub

						! ----------------------------------------------------

						Vperp1N(j)= Vperp1MB(1)
						Vperp2N(j)= Vperp2MB(1)
						VperpN(j)= VperpMB(1)
						VxN(j)= VxMB(1)
						VyN(j)= VyMB(1)
						VzN(j)= VzMB(1)

						! ----------------------------------------------------

						if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

							VparloopKSB11: do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1

								if ((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, 1, 1)%VparGLT(1) > VparMB(1)) &
									.or. (VparMB(1) > SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1))%VparGHT(1))) then
									Vparindk1(j)= 0d0

									exit VparloopKSB11

								end if
								if ((Vparind == 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, 1, Vparind)%VparGLT(1) <= VparMB(1)) &
									.and. (VparMB(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, 1, Vparind)%VparGHT(1))) then
									Vparindk1(j)= Vparind

									exit VparloopKSB11

								end if
								if ((Vparind /= 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, 1, Vparind)%VparGLT(1) < VparMB(1)) &
									.and. (VparMB(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, 1, Vparind)%VparGHT(1))) then
									Vparindk1(j)= Vparind

									exit VparloopKSB11

								end if
							end do VparloopKSB11

							Vperp1loopKSB1: do Vperp1ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), 1
								if ((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, 1, 1)%Vperp1GLT(1) > Vperp1MB(1)) &
									.or. (Vperp1MB(1) > SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), 1, 1)%Vperp1GHT(1))) then
									Vperp1indk1(j)= 0d0

									exit Vperp1loopKSB1

								end if
								if ((Vperp1ind == 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(Vperp1ind, 1, 1)%Vperp1GLT(1) <= Vperp1MB(1)) &
									.and. (Vperp1MB(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(Vperp1ind, 1, 1)%Vperp1GHT(1))) then
									Vperp1indk1(j)= Vperp1ind

									exit Vperp1loopKSB1

								end if
								if ((Vperp1ind /= 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(Vperp1ind, 1, 1)%Vperp1GLT(1) < Vperp1MB(1)) &
									.and. (Vperp1MB(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(Vperp1ind, 1, 1)%Vperp1GHT(1))) then
									Vperp1indk1(j)= Vperp1ind

									exit Vperp1loopKSB1

								end if
							end do Vperp1loopKSB1

							Vperp2loopKSB1: do Vperp2ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), 1
								if ((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, 1, 1)%Vperp2GLT(1) > Vperp2MB(1)) &
									.or. (Vperp2MB(1) > SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), 1)%Vperp2GHT(1))) then
									Vperp2indk1(j)= 0d0

									exit Vperp2loopKSB1

								end if
								if ((Vperp2ind == 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, Vperp2ind, 1)%Vperp2GLT(1) <= Vperp2MB(1)) &
									.and. (Vperp2MB(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, Vperp2ind, 1)%Vperp2GHT(1))) then
									Vperp2indk1(j)= Vperp2ind

									exit Vperp2loopKSB1

								end if
								if ((Vperp2ind /= 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, Vperp2ind, 1)%Vperp2GLT(1) < Vperp2MB(1)) &
									.and. (Vperp2MB(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, Vperp2ind, 1)%Vperp2GHT(1))) then
									Vperp2indk1(j)= Vperp2ind

									exit Vperp2loopKSB1

								end if
							end do Vperp2loopKSB1

						else

							VparloopKSB12: do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1
								if ((Vparind == 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									VCellT(1, Vparind)%VparGLT(1) <= VparMB(1)) &
									.and. (VparMB(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									VCellT(1, Vparind)%VparGHT(1))) then
									Vparindk1(j)= Vparind

									exit VparloopKSB12

								end if
								if ((Vparind /= 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									VCellT(1, Vparind)%VparGLT(1) < VparMB(1)) &
									.and. (VparMB(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									VCellT(1, Vparind)%VparGHT(1))) then
									Vparindk1(j)= Vparind

									exit VparloopKSB12

								end if
							end do VparloopKSB12

							VperploopKSB1: do Vperpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), 1
								if ((Vperpind == 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									VCellT(Vperpind, 1)%VperpGLT(1) <= VperpMB(1)) &
									.and. (VperpMB(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									VCellT(Vperpind, 1)%VperpGHT(1))) then
									Vperpindk1(j)= Vperpind

									exit VperploopKSB1

								end if
								if ((Vperpind /= 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									VCellT(Vperpind, 1)%VperpGLT(1) < VperpMB(1)) &
									.and. (VperpMB(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									VCellT(Vperpind, 1)%VperpGHT(1))) then
									Vperpindk1(j)= Vperpind

									exit VperploopKSB1

								end if
							end do VperploopKSB1

						end if

						! ----------------------------------------------------

					end do

					! ----------------------------------------------------

					! DIAGNOSTIC FOR PROPER LB INJECTED ION GRID CELL:

					do j= (NsTK(1)- dNsTK2(1)- dNsTK3(1))+ 1, NsTK(1)- dNsTK3(1), 1
						if (Qindk1(j) /= NqLB(1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' LB INJECTED PARTICLE= ', j, &
								', Qindk1= ', Qindk1(j), ' NOT IN LOWER GRID CELL FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', AND STATISTICAL TIME-STEP= ', nn, ' IN KINETIC SOLVER B SUBROUTINE' &
								// achar(27) // '[0m.'
						end if
					end do

					! ----------------------------------------------------

					! INITIALIZE UPPER BOUNDARY INJECTED IONS ON INJECTION TIME-STEP:

					! ----------------------------------------------------

					! Initialize all UB injected particle variables
					do j= (NsTK(1)- dNsTK3(1))+ 1, NsTK(1), 1

						! ----------------------------------------------------

						! Inject ions at upper altitude boundary
						ENAflag(j)= .false.
						ENAflagN0ind(j)= n
						AEAmag(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*abs(SpecieT(s)%FluxTubeT(f)%EAmagRT(nn- 1, NqUB(1))))
						AGmag(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*abs(SpecieT(s)%FluxTubeT(f)%EGmagRT(NqUB(1))))
						AEPmag(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*abs(SpecieT(s)%FluxTubeT(f)%EPmagRT(NqUB(1))))

						! ----------------------------------------------------

						! Set initial dipole coordinates
						qNp(1)= SpecieT(s)%FluxTubeT(f)%QCellT(NqUB(1))%qGCT(1)
						pNp(1)= SpecieT(s)%FluxTubeT(f)%QCellT(1)%pGCT(1)

						! ----------------------------------------------------

						if (SpecieT(s)%FluxTubeT(f)%QCellT(1)%qGLT(1) <= 0) then ! S. Magnetic Hemisphere
							QloopKSB3: do Qind= NqLB(1), NqUB(1), 1
								if ((Qind == NqLB(1)) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGLT(1) <= qNp(1)) &
									.and. (qNp(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGHT(1))) then
									Qindk1(j)= Qind

									exit QloopKSB3

								else if ((Qind /= NqLB(1)) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGLT(1) < qNp(1)) &
									.and. (qNp(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGHT(1))) then
									Qindk1(j)= Qind

									exit QloopKSB3

								end if
							end do QloopKSB3

							if (SpecieT(s)%FluxTubeT(f)%QCellT(NqLB(1))%qGLT(1) > qNp(1)) then
								! SMH Lower boundary escape
								Qindk1(j)= 0d0
							else if (SpecieT(s)%FluxTubeT(f)%QCellT(NqUB(1))%qGHT(1) < qNp(1)) then
								! SMH Upper boundary escape
								Qindk1(j)= -1d0
							end if
						end if

						! ----------------------------------------------------

						if (SpecieT(s)%FluxTubeT(f)%QCellT(1)%qGLT(1) > 0) then ! N. Magnetic Hemisphere
							QloopKSB4: do Qind= NqLB(1), NqUB(1), 1
								if ((Qind == NqLB(1)) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGLT(1) >= qNp(1)) &
									.and. (qNp(1) >= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGHT(1))) then
									Qindk1(j)= Qind

									exit QloopKSB4

								else if ((Qind /= NqLB(1)) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGLT(1) > qNp(1)) &
									.and. (qNp(1) >= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGHT(1))) then
									Qindk1(j)= Qind

									exit QloopKSB4

								end if
							end do QloopKSB4

							if (SpecieT(s)%FluxTubeT(f)%QCellT(NqLB(1))%qGLT(1) < qNp(1)) then
								! NMH Lower boundary escape
								Qindk1(j)= 0d0
							else if (SpecieT(s)%FluxTubeT(f)%QCellT(NqUB(1))%qGHT(1) > qNp(1)) then
								! NMH Upper boundary escape
								Qindk1(j)= -1d0
							end if
						end if

						! ----------------------------------------------------

						! DIAGNOSTIC FLAG FOR PARTICLE POSITIONS WITHIN BOUNDS:

						if (SpecieT(s)%FluxTubeT(f)%QCellT(1)%qGLT(1) <= 0) then ! S. Magnetic Hemisphere)
							if ((qNp(1) < SpecieT(s)%FluxTubeT(f)%QCellT(NqUB(1))%qGLT(1)) .or. &
								(qNp(1) > SpecieT(s)%FluxTubeT(f)%QCellT(NqUB(1))%qGHT(1))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qniICT IS', &
									' OUTSIDE UPPER BOUNDARY BOUNDS IN SOUTH MAGNETIC HEMISPHERE FOR SPECIE= ', s, ', FLUX TUBE= ', &
									f, ', AND Qind= ', Qind, ' IN KINETIC SOLVER B SUBROUTINE' // achar(27) // '[0m.'
							end if
						end if
						if (SpecieT(s)%FluxTubeT(f)%QCellT(1)%qGLT(1) > 0) then ! N. Magnetic Hemisphere
							if ((qNp(1) > SpecieT(s)%FluxTubeT(f)%QCellT(NqUB(1))%qGLT(1)) .or. &
								(qNp(1) < SpecieT(s)%FluxTubeT(f)%QCellT(NqUB(1))%qGHT(1))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qNp IS', &
									' OUTSIDE UPPER BOUNDARY BOUNDS IN NORTH MAGNETIC HEMISPHERE FOR SPECIE= ', s, ', FLUX TUBE= ', &
									f, ', AND Qind= ', Qind, ' IN KINETIC SOLVER B SUBROUTINE' // achar(27) // '[0m.'
							end if
						end if

						! ----------------------------------------------------

						call RK4DipolePolynomialSolverSub

						! ----------------------------------------------------

						xN(j)= xfinalRK4(1)
						yN(j)= yfinalRK4(1)
						zN(j)= zfinalRK4(1)

						! ----------------------------------------------------

						! Inject ions with MB distributed velocities
						thetaMB(1)= thetafinalRK4(1)
						phiMB(1)= phifinalRK4(1)
						ellMB(1)= 1d0+ 3d0*(cos(thetaMB(1)))**2d0
						TsMB(1)= SpecieT(s)%FluxTubeT(f)%QCellT(NqUB(1))%TsT(1)
						LBflag(1)= 0d0

						! ----------------------------------------------------

						call MBVelocityDistributionSub

						! ----------------------------------------------------

						Vperp1N(j)= Vperp1MB(1)
						Vperp2N(j)= Vperp2MB(1)
						VperpN(j)= VperpMB(1)
						VxN(j)= VxMB(1)
						VyN(j)= VyMB(1)
						VzN(j)= VzMB(1)

						! ----------------------------------------------------

						if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

							VparloopKSB21: do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1

								if ((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, 1, 1)%VparGLT(1) > VparMB(1)) &
									.or. (VparMB(1) > SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1))%VparGHT(1))) then
									Vparindk1(j)= 0d0

									exit VparloopKSB21

								end if
								if ((Vparind == 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, 1, Vparind)%VparGLT(1) <= VparMB(1)) &
									.and. (VparMB(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, 1, Vparind)%VparGHT(1))) then
									Vparindk1(j)= Vparind

									exit VparloopKSB21

								end if
								if ((Vparind /= 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, 1, Vparind)%VparGLT(1) < VparMB(1)) &
									.and. (VparMB(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, 1, Vparind)%VparGHT(1))) then
									Vparindk1(j)= Vparind

									exit VparloopKSB21

								end if
							end do VparloopKSB21

							Vperp1loopKSB2: do Vperp1ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), 1
								if ((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, 1, 1)%Vperp1GLT(1) > Vperp1MB(1)) &
									.or. (Vperp1MB(1) > SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), 1, 1)%Vperp1GHT(1))) then
									Vperp1indk1(j)= 0d0

									exit Vperp1loopKSB2

								end if
								if ((Vperp1ind == 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(Vperp1ind, 1, 1)%Vperp1GLT(1) <= Vperp1MB(1)) &
									.and. (Vperp1MB(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(Vperp1ind, 1, 1)%Vperp1GHT(1))) then
									Vperp1indk1(j)= Vperp1ind

									exit Vperp1loopKSB2

								end if
								if ((Vperp1ind /= 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(Vperp1ind, 1, 1)%Vperp1GLT(1) < Vperp1MB(1)) &
									.and. (Vperp1MB(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(Vperp1ind, 1, 1)%Vperp1GHT(1))) then
									Vperp1indk1(j)= Vperp1ind

									exit Vperp1loopKSB2

								end if
							end do Vperp1loopKSB2

							Vperp2loopKSB2: do Vperp2ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), 1
								if ((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, 1, 1)%Vperp2GLT(1) > Vperp2MB(1)) &
									.or. (Vperp2MB(1) > SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), 1)%Vperp2GHT(1))) then
									Vperp2indk1(j)= 0d0

									exit Vperp2loopKSB2

								end if
								if ((Vperp2ind == 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, Vperp2ind, 1)%Vperp2GLT(1) <= Vperp2MB(1)) &
									.and. (Vperp2MB(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, Vperp2ind, 1)%Vperp2GHT(1))) then
									Vperp2indk1(j)= Vperp2ind

									exit Vperp2loopKSB2

								end if
								if ((Vperp2ind /= 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, Vperp2ind, 1)%Vperp2GLT(1) < Vperp2MB(1)) &
									.and. (Vperp2MB(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									V2PerpCellT(1, Vperp2ind, 1)%Vperp2GHT(1))) then
									Vperp2indk1(j)= Vperp2ind

									exit Vperp2loopKSB2

								end if
							end do Vperp2loopKSB2

						else

							VparloopKSB22: do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1
								if ((Vparind == 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									VCellT(1, Vparind)%VparGLT(1) <= VparMB(1)) &
									.and. (VparMB(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									VCellT(1, Vparind)%VparGHT(1))) then
									Vparindk1(j)= Vparind

									exit VparloopKSB22

								end if
								if ((Vparind /= 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									VCellT(1, Vparind)%VparGLT(1) < VparMB(1)) &
									.and. (VparMB(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									VCellT(1, Vparind)%VparGHT(1))) then
									Vparindk1(j)= Vparind

									exit VparloopKSB22

								end if
							end do VparloopKSB22

							VperploopKSB2: do Vperpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), 1
								if ((Vperpind == 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									VCellT(Vperpind, 1)%VperpGLT(1) <= VperpMB(1)) &
									.and. (VperpMB(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									VCellT(Vperpind, 1)%VperpGHT(1))) then
									Vperpindk1(j)= Vperpind

									exit VperploopKSB2

								end if
								if ((Vperpind /= 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									VCellT(Vperpind, 1)%VperpGLT(1) < VperpMB(1)) &
									.and. (VperpMB(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
									VCellT(Vperpind, 1)%VperpGHT(1))) then
									Vperpindk1(j)= Vperpind

									exit VperploopKSB2

								end if
							end do VperploopKSB2

						end if

					end do

					! ----------------------------------------------------

					! DIAGNOSTIC FOR PROPER UB INJECTED ION GRID CELL:

					do j= (NsTK(1)- dNsTK3(1))+ 1, NsTK(1), 1
						if (Qindk1(j) /= NqUB(1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' UB INJECTED PARTICLE= ', j, &
								', Qindk1= ', Qindk1(j), ' NOT IN UPPER GRID CELL FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', AND INJECTION TIME-STEP= ', nn, ' IN KINETIC SOLVER B SUBROUTINE' &
								// achar(27) // '[0m.'
						end if
					end do

					! ----------------------------------------------------

				end if
			end do nnloopKS2

			! ----------------------------------------------------

			! UPDATE ALL PARTICLES ON ALL OTHER TIME-STEPS:

			! ----------------------------------------------------

			! Update variables on all other time-steps (without injection)
			if (nnFlag(1) == 0) then

				do j= 1, NsTK(1), 1

					if (ENAflag(j) .eqv. .false.) then
						AEAmag(j)= AEAmagN(j)
						AGmag(j)= AGmagN(j)
						AEPmag(j)= AEPmagN(j)
					end if

					! ----------------------------------------------------

					! Update particle phase-space coordinates
					x(j)= xN(j)
					y(j)= yN(j)
					z(j)= zN(j)
					if (ENAflag(j) .eqv. .false.) then
						Vperp1(j)= Vperp1N(j)
						Vperp2(j)= Vperp2N(j)
						Vperp(j)= VperpN(j)
					else if ((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
						(ENAflag(j) .eqv. .true.)) then
						Vperp1(j)= 0d0
						Vperp2(j)= 0d0
						Vperp(j)= 0d0
					end if
					Vx(j)= VxN(j)
					Vy(j)= VyN(j)
					Vz(j)= VzN(j)

					! ----------------------------------------------------

					! DIAGNOSTIC FLAGS FOR NAN VALUES:

					if (isnan(x(j)) .eqv. .true.) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' x VALUE IS NAN FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
							' IN KINETIC SOLVER B SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if (isnan(y(j)) .eqv. .true.) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' y VALUE IS NAN FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
							' IN KINETIC SOLVER B SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if (isnan(z(j)) .eqv. .true.) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' z VALUE IS NAN FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
							' IN KINETIC SOLVER B SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if (isnan(Vperp1(j)) .eqv. .true.) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' Vperp1 VALUE IS NAN FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
							' IN KINETIC SOLVER B SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if (isnan(Vperp2(j)) .eqv. .true.) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' Vperp2 VALUE IS NAN FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
							' IN KINETIC SOLVER B SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if (isnan(Vperp(j)) .eqv. .true.) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' Vperp VALUE IS NAN FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
							' IN KINETIC SOLVER B SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if (isnan(Vx(j)) .eqv. .true.) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' Vx VALUE IS NAN FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
							' IN KINETIC SOLVER B SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if (isnan(Vy(j)) .eqv. .true.) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' Vy VALUE IS NAN FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
							' IN KINETIC SOLVER B SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if (isnan(Vz(j)) .eqv. .true.) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' Vz VALUE IS NAN FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
							' IN KINETIC SOLVER B SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

				end do

				! ----------------------------------------------------

				! Perform iterative updates of particle positions and velocities
				jloopKS4: do j= 1, NsTK(1), 1

					! ----------------------------------------------------

					! DIAGNOSTIC FLAG FOR CORRECT ION/ENA PARTICLE TYPE:

					if ((ENAflag(j) .eqv. .true.) .and. (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 0)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ENA PRODUCED WITH', &
							' CHARGE-EXCHANGE TURNED OFF FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
							' IN KINETIC SOLVER B SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

					call KineticRK4UpdateSub

					! ----------------------------------------------------

					! DIAGNOSTIC FLAGS FOR NAN VALUES:

					if (isnan(xN(j)) .eqv. .true.) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' xN VALUE IS NAN FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
							' IN KINETIC SOLVER B SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if (isnan(yN(j)) .eqv. .true.) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' yN VALUE IS NAN FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
							' IN KINETIC SOLVER B SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if (isnan(zN(j)) .eqv. .true.) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' zN VALUE IS NAN FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
							' IN KINETIC SOLVER B SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if (isnan(Vperp1N(j)) .eqv. .true.) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' Vperp1N VALUE IS NAN FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
							' IN KINETIC SOLVER B SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if (isnan(Vperp2N(j)) .eqv. .true.) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' Vperp2N VALUE IS NAN FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
							' IN KINETIC SOLVER B SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if (isnan(VperpN(j)) .eqv. .true.) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' VperpN VALUE IS NAN FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
							' IN KINETIC SOLVER B SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if (isnan(VxN(j)) .eqv. .true.) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' VxN VALUE IS NAN FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
							' IN KINETIC SOLVER B SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if (isnan(VyN(j)) .eqv. .true.) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' VyN VALUE IS NAN FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
							' IN KINETIC SOLVER B SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if (isnan(VzN(j)) .eqv. .true.) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' VzN VALUE IS NAN FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
							' IN KINETIC SOLVER B SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

					if (j < NsTK(1)) then
						cycle jloopKS4
					end if
					if (j >= NsTK(1)) then
						exit jloopKS4
					end if

				end do jloopKS4

				! ----------------------------------------------------

			end if

			! ----------------------------------------------------

			! GATHER STATISTICS AND FORM MACROSCOPIC PARAMETERS:

			! Compute all ion fluid moments
			call ParticleCountsSub
			call DistributionFunctionsSub
			call ZerothIonMomentSub
			call FirstPerpIonMomentSub
			call FirstParIonMomentSub
			call SecondIonMomentSub
			call SecondPerpIonMomentSub
			call SecondParIonMomentSub

			! ----------------------------------------------------

			! DIAGNOSTIC FLAG FOR NON-ZERO GRID CELLS IN SPIN-UP:

			if ((SpecieT(s)%FluxTubeT(f)%DENSITYOUTPUTflagT(1) == 1) .and. (rank == 0)) then
				do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
					if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
						(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1)))) then
						if (nn == SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1) then
							do Qind= NqLB(1), NqUB(1), 1
								if (SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) == 0) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' SPIN-UP SIMULATION HAS ZERO DENSITY AT FINAL TIME FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN KINETIC SOLVER B SUBROUTINE' &
										// achar(27) // '[0m.'
								end if
							end do
						end if
					end if
				end do
			end if

			! ----------------------------------------------------

			! Compute ambipolar and parallel electric fields
			call AmbipolarEfieldSub
			call GravfieldSub
			call PotentialStructureSub

			! ----------------------------------------------------

			! Compute all ENA fluid moments
			!if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
				!call ZerothENAMomentSub
				!call FirstPENAMomentSub
				!call FirstQENAMomentSub
				!call FirstPHIENAMomentSub
				!call SecondENAMomentSub
				!call SecondPENAMomentSub
				!call SecondQENAMomentSub
				!call SecondPHIENAMomentSub
			!end if

			! ----------------------------------------------------

			! EXPORT INTERMEDIARY KINETIC DATA:

			do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
				if (((n /= 1) .and. (nn /= 1)) .and. &
					(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1)) .and. &
					((n /= SpecieT(s)%FluxTubeT(f)%NtT(1)) .and. &
					(nn /= SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))) then

					call DataExport2Sub

				end if
			end do

			! Export injection data on injection time-step
			if ((SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1) == 1)  .and. &
				(rank == 0)) then
				do nn= 1, SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1, 1
					if (((n /= 1) .and. (nn /= 1)) .and. &
						(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%Q0ndatfacT(1)) .and. &
						((n /= SpecieT(s)%FluxTubeT(f)%NtT(1)) .and. &
						(nn /= SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1))) then

						write(nnstring, '(I5)') nn
						write(sstring, '(I5)') s
						write(fstring, '(I5)') f

						expstring= adjustl(adjustr(rankstring) &
							// '_' // adjustl(adjustr(nnstring) &
							// '_' // adjustl(adjustr(sstring) &
							// '_' // adjustl(adjustr(fstring) // '_'))))

						NsnTfile= adjustl(adjustr(expstring) // adjustl(adjustr('NsnTfort.bin')))
						open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
							adjustl(adjustr(NsnTfile))), &
							status= 'replace', form= 'unformatted', access= 'stream')
						write(expint) SpecieT(s)%FluxTubeT(f)%NsnT(nn)
						close(expint)

						NsnRRTfile= adjustl(adjustr(expstring) // adjustl(adjustr('NsnRRTfort.bin')))
						open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
							adjustl(adjustr(NsnRRTfile))), &
							status= 'replace', form= 'unformatted', access= 'stream')
						write(expint) SpecieT(s)%FluxTubeT(f)%NsnRRT(nn)
						close(expint)

						NqLBoutfluxIonRTfile= adjustl(adjustr(expstring) // adjustl(adjustr('NqLBoutfluxIonRTfort.bin')))
						open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
							adjustl(adjustr(NqLBoutfluxIonRTfile))), &
							status= 'replace', form= 'unformatted', access= 'stream')
						write(expint) SpecieT(s)%FluxTubeT(f)%NqLBoutfluxIonRT(nn)
						close(expint)

						LBoutfluxIonRTfile= adjustl(adjustr(expstring) // adjustl(adjustr('LBoutfluxIonRTfort.bin')))
						open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
							adjustl(adjustr(LBoutfluxIonRTfile))), &
							status= 'replace', form= 'unformatted', access= 'stream')
						write(expint) SpecieT(s)%FluxTubeT(f)%LBoutfluxIonRT(nn)
						close(expint)

						NqUBoutfluxIonRTfile= adjustl(adjustr(expstring) // adjustl(adjustr('NqUBoutfluxIonRTfort.bin')))
						open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
							adjustl(adjustr(NqUBoutfluxIonRTfile))), &
							status= 'replace', form= 'unformatted', access= 'stream')
						write(expint) SpecieT(s)%FluxTubeT(f)%NqUBoutfluxIonRT(nn)
						close(expint)

						UBoutfluxIonRTfile= adjustl(adjustr(expstring) // adjustl(adjustr('UBoutfluxIonRTfort.bin')))
						open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
							adjustl(adjustr(UBoutfluxIonRTfile))), &
							status= 'replace', form= 'unformatted', access= 'stream')
						write(expint) SpecieT(s)%FluxTubeT(f)%UBoutfluxIonRT(nn)
						close(expint)

						if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
							NqLBoutfluxENARTfile= adjustl(adjustr(expstring) // adjustl(adjustr('NqLBoutfluxENARTfort.bin')))
							open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
								adjustl(adjustr(NqLBoutfluxENARTfile))), &
								status= 'replace', form= 'unformatted', access= 'stream')
							write(expint) SpecieT(s)%FluxTubeT(f)%NqLBoutfluxENART(nn)
							close(expint)

							LBoutfluxENARTfile= adjustl(adjustr(expstring) // adjustl(adjustr('LBoutfluxENARTfort.bin')))
							open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
								adjustl(adjustr(LBoutfluxENARTfile))), &
								status= 'replace', form= 'unformatted', access= 'stream')
							write(expint) SpecieT(s)%FluxTubeT(f)%LBoutfluxENART(nn)
							close(expint)

							NqUBoutfluxENARTfile= adjustl(adjustr(expstring) // adjustl(adjustr('NqUBoutfluxENARTfort.bin')))
							open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
								adjustl(adjustr(NqUBoutfluxENARTfile))), &
								status= 'replace', form= 'unformatted', access= 'stream')
							write(expint) SpecieT(s)%FluxTubeT(f)%NqUBoutfluxENART(nn)
							close(expint)

							UBoutfluxENARTfile= adjustl(adjustr(expstring) // adjustl(adjustr('UBoutfluxENARTfort.bin')))
							open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
								adjustl(adjustr(UBoutfluxENARTfile))), &
								status= 'replace', form= 'unformatted', access= 'stream')
							write(expint) SpecieT(s)%FluxTubeT(f)%UBoutfluxENART(nn)
							close(expint)
						end if

						if (SpecieT(s)%FluxTubeT(f)%LBCONDITIONflagT(1) == 1) then
							LBNetDensityTfile= adjustl(adjustr(expstring) // adjustl(adjustr('LBNetDensityTfort.bin')))
							open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
								adjustl(adjustr(LBNetDensityTfile))), &
								status= 'replace', form= 'unformatted', access= 'stream')
							write(expint) SpecieT(s)%FluxTubeT(f)%LBNetDensityT(nn)
							close(expint)

							UBNetDensityTfile= adjustl(adjustr(expstring) // adjustl(adjustr('UBNetDensityTfort.bin')))
							open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
								adjustl(adjustr(UBNetDensityTfile))), &
								status= 'replace', form= 'unformatted', access= 'stream')
							write(expint) SpecieT(s)%FluxTubeT(f)%UBNetDensityT(nn)
							close(expint)

						end if

					end if
  			end do
  		end if

			! ----------------------------------------------------

			! EXPORT FINAL KINETIC DATA:

			do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
				if (((n /= 1) .and. (nn /= 1)) .and. &
					(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1)) .and. &
					((n == SpecieT(s)%FluxTubeT(f)%NtT(1)) .and. &
					(nn == SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))) then
					call DataExport2Sub
				end if
			end do

			! Export injection data on injection time-step
			if ((SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1) == 1)  .and. &
				(rank == 0)) then
				do nn= 1, SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1, 1
					if (((n /= 1) .and. (nn /= 1)) .and. &
						(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%Q0ndatfacT(1)) .and. &
						((n == SpecieT(s)%FluxTubeT(f)%NtT(1)) .and. &
						(nn == SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1))) then

						write(nnstring, '(I5)') nn
						write(sstring, '(I5)') s
						write(fstring, '(I5)') f

						expstring= adjustl(adjustr(rankstring) &
							// '_' // adjustl(adjustr(nnstring) &
							// '_' // adjustl(adjustr(sstring) &
							// '_' // adjustl(adjustr(fstring) // '_'))))

						NsnTfile= adjustl(adjustr(expstring) // adjustl(adjustr('NsnTfort.bin')))
						open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
							adjustl(adjustr(NsnTfile))), &
							status= 'replace', form= 'unformatted', access= 'stream')
						write(expint) SpecieT(s)%FluxTubeT(f)%NsnT(nn)
						close(expint)

						NsnRRTfile= adjustl(adjustr(expstring) // adjustl(adjustr('NsnRRTfort.bin')))
						open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
							adjustl(adjustr(NsnRRTfile))), &
							status= 'replace', form= 'unformatted', access= 'stream')
						write(expint) SpecieT(s)%FluxTubeT(f)%NsnRRT(nn)
						close(expint)

						NqLBoutfluxIonRTfile= adjustl(adjustr(expstring) // adjustl(adjustr('NqLBoutfluxIonRTfort.bin')))
						open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
							adjustl(adjustr(NqLBoutfluxIonRTfile))), &
							status= 'replace', form= 'unformatted', access= 'stream')
						write(expint) SpecieT(s)%FluxTubeT(f)%NqLBoutfluxIonRT(nn)
						close(expint)

						LBoutfluxIonRTfile= adjustl(adjustr(expstring) // adjustl(adjustr('LBoutfluxIonRTfort.bin')))
						open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
							adjustl(adjustr(LBoutfluxIonRTfile))), &
							status= 'replace', form= 'unformatted', access= 'stream')
						write(expint) SpecieT(s)%FluxTubeT(f)%LBoutfluxIonRT(nn)
						close(expint)

						NqUBoutfluxIonRTfile= adjustl(adjustr(expstring) // adjustl(adjustr('NqUBoutfluxIonRTfort.bin')))
						open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
							adjustl(adjustr(NqUBoutfluxIonRTfile))), &
							status= 'replace', form= 'unformatted', access= 'stream')
						write(expint) SpecieT(s)%FluxTubeT(f)%NqUBoutfluxIonRT(nn)
						close(expint)

						UBoutfluxIonRTfile= adjustl(adjustr(expstring) // adjustl(adjustr('UBoutfluxIonRTfort.bin')))
						open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
							adjustl(adjustr(UBoutfluxIonRTfile))), &
							status= 'replace', form= 'unformatted', access= 'stream')
						write(expint) SpecieT(s)%FluxTubeT(f)%UBoutfluxIonRT(nn)
						close(expint)

						if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
							NqLBoutfluxENARTfile= adjustl(adjustr(expstring) // adjustl(adjustr('NqLBoutfluxENARTfort.bin')))
							open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
								adjustl(adjustr(NqLBoutfluxENARTfile))), &
								status= 'replace', form= 'unformatted', access= 'stream')
							write(expint) SpecieT(s)%FluxTubeT(f)%NqLBoutfluxENART(nn)
							close(expint)

							LBoutfluxENARTfile= adjustl(adjustr(expstring) // adjustl(adjustr('LBoutfluxENARTfort.bin')))
							open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
								adjustl(adjustr(LBoutfluxENARTfile))), &
								status= 'replace', form= 'unformatted', access= 'stream')
							write(expint) SpecieT(s)%FluxTubeT(f)%LBoutfluxENART(nn)
							close(expint)

							NqUBoutfluxENARTfile= adjustl(adjustr(expstring) // adjustl(adjustr('NqUBoutfluxENARTfort.bin')))
							open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
								adjustl(adjustr(NqUBoutfluxENARTfile))), &
								status= 'replace', form= 'unformatted', access= 'stream')
							write(expint) SpecieT(s)%FluxTubeT(f)%NqUBoutfluxENART(nn)
							close(expint)

							UBoutfluxENARTfile= adjustl(adjustr(expstring) // adjustl(adjustr('UBoutfluxENARTfort.bin')))
							open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
								adjustl(adjustr(UBoutfluxENARTfile))), &
								status= 'replace', form= 'unformatted', access= 'stream')
							write(expint) SpecieT(s)%FluxTubeT(f)%UBoutfluxENART(nn)
							close(expint)
						end if

						if (SpecieT(s)%FluxTubeT(f)%LBCONDITIONflagT(1) == 1) then
							LBNetDensityTfile= adjustl(adjustr(expstring) // adjustl(adjustr('LBNetDensityTfort.bin')))
							open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
								adjustl(adjustr(LBNetDensityTfile))), &
								status= 'replace', form= 'unformatted', access= 'stream')
							write(expint) SpecieT(s)%FluxTubeT(f)%LBNetDensityT(nn)
							close(expint)

							UBNetDensityTfile= adjustl(adjustr(expstring) // adjustl(adjustr('UBNetDensityTfort.bin')))
							open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
								adjustl(adjustr(UBNetDensityTfile))), &
								status= 'replace', form= 'unformatted', access= 'stream')
							write(expint) SpecieT(s)%FluxTubeT(f)%UBNetDensityT(nn)
							close(expint)

						end if

					end if
  			end do
  		end if

			! ----------------------------------------------------

			! PRINT OUT STATISTICAL TIME-STEP:

			do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
				if ((n /= 1) .and. (nn /= 1) .and. &
					(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1))) then
					if (rank == 0) then
						call cpu_time(KS0End)
						write(Nsstring, '(i10)') SpecieT(s)%FluxTubeT(f)%NsnRRT(nn)
						write(nnstring, '(I5)') nn
						write(sstring, '(I5)') s
						write(fstring, '(I5)') f
						write(Timestring, '(I5)') nint(Time(1))
						write(KS0string, '(i10)')  nint(KS0End)
						write(*, *) trim('** COMPLETE: STATISTICAL TIME-STEP= ' &
							// adjustl(nnstring)) // &
							trim(', RANK= ' // adjustl(rankstring)) // &
							trim(', PARTICLE SPECIE= ' // adjustl(sstring)) // &
							trim(', FLUX-TUBE= ' // adjustl(fstring)) // &
							trim(', SIM-TIME= ' // adjustl(Timestring)) // &
							trim(' s., REAL-TIME= ' // adjustl(KS0string)) // &
							trim(' s., TOTAL PARTICLE NUMBER= ' // adjustl(Nsstring))
					end if
				end if
			end do

			! ----------------------------------------------------

			! PRINT OUT INJECTION TIME-STEP:

			!do nn= 1, SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1, 1
			!	if ((n /= 1) .and. (nn /= 1) .and. &
			!		(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%Q0ndatfacT(1))) then
			!		if (rank == 0) then
			!			call cpu_time(KS0End)
			!			write(Nsstring, '(i10)') SpecieT(s)%FluxTubeT(f)%NsnRRT(nn)
			!			write(nnstring, '(I5)') nn
			!			write(sstring, '(I5)') s
			!			write(fstring, '(I5)') f
			!			write(Timestring, '(I5)') nint(Time(1))
			!			write(KS0string, '(i10)')  nint(KS0End)
			!			write(*, *) trim('** COMPLETE: INJECTION TIME-STEP= ' &
			!				// adjustl(nnstring)) // &
			!				trim(', RANK= ' // adjustl(rankstring)) // &
			!				trim(', PARTICLE SPECIE= ' // adjustl(sstring)) // &
			!				trim(', FLUX-TUBE= ' // adjustl(fstring)) // &
			!				trim(', SIM-TIME= ' // adjustl(Timestring)) // &
			!				trim(' s., REAL-TIME= ' // adjustl(KS0string)) // &
			!				trim(' s., TOTAL PARTICLE NUMBER= ' // adjustl(Nsstring))
			!		end if
			!	end if
			!end do

			! ----------------------------------------------------

		end if

		! ----------------------------------------------------

	end subroutine KineticSolverBSub

end module KineticSolverB
