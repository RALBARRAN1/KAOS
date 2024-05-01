module DipolePolynomialSolver

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	5- DIPOLE POLYNOMIAL SOLVER:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! SOLVE QUARTIC DIPOLE POLYNOMIAL BY ANALYTIC ORDER DECOMPOSITION FOR INITIAL PARTICLE
! POSITIONS IN RECTANGULAR CONFIG. SPACE:

	subroutine DipolePolynomialSolverSub

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				do Qind= NqLB(1), NqUB(1), 1

					! ----------------------------------------------------

					allocate(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%rfinalICT &
						(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1)), &
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%thetafinalICT &
						(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1)), &
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%phifinalICT &
						(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1)), &
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%xfinalICT &
						(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1)), &
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%yfinalICT &
						(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1)), &
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%zfinalICT &
						(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1)), &
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qfinalICT &
						(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1)), &
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%pfinalICT &
						(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1)), &
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%ellfinalICT &
						(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1)))

					! ----------------------------------------------------

					! USE q AND p VALUES OF niIC TO GET x, y, z VALUES BY QUARTIC
					! POLYNOMICAL ROOT-FINDER IN FAindIC:

					do FAindIC= 1, SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1), 1

						! ----------------------------------------------------

						ApIC(1)= 0d0  ! From dipole quartic polynomial form
						BpIC(1)= 0d0
						CpIC(1)= 1d0/(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%pniICT(FAindIC)* &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniICT(FAindIC)**2d0)
						DpIC(1)= -1d0/(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniICT(FAindIC)**2d0)

						AbIC(1)= -BpIC(1) ! From resolvent cubic form
						BbIC(1)= ApIC(1)*CpIC(1)- 4d0*DpIC(1)
						CbIC(1)= 4d0*BpIC(1)*DpIC(1)- CpIC(1)**2d0- (ApIC(1)**2d0)*DpIC(1)

						! Resolvent cubic discriminant
						D3IC(1)= (BbIC(1)**2d0)*(AbIC(1)**2d0)- 4d0*CbIC(1)*(AbIC(1)**3d0)- &
							4d0*(BbIC(1)**3d0)+ 18d0*AbIC(1)* BbIC(1)*CbIC(1)- 27d0*(CbIC(1)**2d0)

						D4IC(1)= ((CpIC(1)**2d0)*(BpIC(1)**2d0)*(ApIC(1)**2d0)- 4d0*(CpIC(1)**3d0)* &
							(ApIC(1)**3d0)- 4d0*(CpIC(1)**2d0)*(BpIC(1)**3d0)+ 18d0*(CpIC(1)**3d0)* &
							BpIC(1)*ApIC(1)- 27d0*(CpIC(1)**4d0)+ 256d0*(DpIC(1)**3d0))+ DpIC(1)* &
							(-4d0*(BpIC(1)**3d0)*(ApIC(1)**2d0)+ 18d0*CpIC(1)*BpIC(1)*(ApIC(1)**3d0)+ &
							16d0*(BpIC(1)**4d0)- 80d0*CpIC(1)*(BpIC(1)**2d0)*ApIC(1)- 6d0*(CpIC(1)**2d0)* &
							(ApIC(1)**2d0)+ 144d0*(CpIC(1)**2d0)*BpIC(1))+ (DpIC(1)**2d0)*(-27d0* &
							(ApIC(1)**4d0)+ 144d0*BpIC(1)*(ApIC(1)**2d0)- 128d0*(BpIC(1)**2d0)- &
							192d0*CpIC(1)*ApIC(1)) ! Original quartic disrcriminant

						if (abs(D3IC(1)- D4IC(1)) > 1d-14) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' BAD DISCRIMINANT', &
								' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
								', AND FAindIC= ', FAindIC, ' IN DIPOLE POLYNOMIAL SOLVER', &
								' SUBROUTINE' // achar(27) // '[0m.'
						end if

						AbbIC(1)= BbIC(1)- (AbIC(1)**2d0)/3d0
						BbbIC(1)= 2d0*(AbIC(1)**3d0)/27d0- AbIC(1)*BbIC(1)/3d0+ CbIC(1)
						DeltaIC(1)= -BbbIC(1)/2d0
						EpsilonIC(1)= AbbIC(1)/3d0

						ThetapIC(1)= acos(DeltaIC(1)/(i*sqrt(EpsilonIC(1)**3d0)))

						sigma23IC(1)= 2d0*i*sqrt(EpsilonIC(1))*cos((ThetapIC(1)+ 4d0*pi)/3d0)- &
							AbIC(1)/3d0
						sigmaIC(1)= (abs(real(sigma23IC(1))))

						muIC(1)= sqrt((ApIC(1)**2d0)/4d0- BpIC(1)+ sigmaIC(1))

						if (muIC(1) /= 0) then
							nuIC(1)= sqrt(3d0*(ApIC(1)**2d0)/4d0- (muIC(1)**2d0)- 2d0*BpIC(1)+ &
								(4d0*ApIC(1)*BpIC(1)- 8d0*CpIC(1)- (ApIC(1)**3d0))/(4d0*muIC(1)))
							piiIC(1)= sqrt(3d0*(ApIC(1)**2d0)/4d0- (muIC(1)**2d0)- 2d0*BpIC(1)- &
								(4d0*ApIC(1)*BpIC(1)- 8d0*CpIC(1)- (ApIC(1)**3d0))/(4d0*muIC(1)))
						else if (muIC(1) == 0) then
							nuIC(1)= sqrt(3d0*(ApIC(1)**2d0)/4d0- 2d0*BpIC(1)+ &
								2d0*sqrt((sigmaIC(1)**2d0)- 4d0*DpIC(1)))
							piiIC(1)= sqrt(3d0*(ApIC(1)**2d0)/4d0- 2d0*BpIC(1)- &
								2d0*sqrt((sigmaIC(1)**2d0)- 4d0*DpIC(1)))
						end if

						! Original real quartic root > 0
						gamma3IC(1)= -ApIC(1)/4d0- muIC(1)/2d0+ piiIC(1)/2d0

						! Revert quartic roots with (q, p) into (r, theta)

						r3IC(1)= gamma3IC(1)*RE
						theta3IC(1)= asin(sqrt(r3IC(1) &
							/(RE*SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%pniICT(FAindIC))))

						qtest3IC(1)= (RE**2d0)*cos(theta3IC(1))/(r3IC(1)**2d0)

						ptest3IC(1)= r3IC(1)/(RE*(sin(theta3IC(1))**2d0))

						! Note: gamma3IC is real positive root. For q< 0 (q> 0), phase-shift thetaIC by
						! pi (0) to get correct sign of q (above and below dipole equator).

						rfinalIC(1)= r3IC(1) ! Get final (r, theta) values

						if (SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniICT(FAindIC) &
							<= 0) then ! Phase shift root 3 solution
							thetafinalIC(1)= pi- theta3IC(1)
						else if (SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniICT(FAindIC) &
							> 0) then
							thetafinalIC(1)= theta3IC(1)
						end if

						! Note: Get final (y, z) values with phi= pi/2 s.t. x= 0.

						phifinalIC(1)= pi/2d0 ! Select (arbitrary for dipole) B longitude

						xfinalIC(1)= rfinalIC(1)*sin(thetafinalIC(1))*cos(phifinalIC(1))
						yfinalIC(1)= rfinalIC(1)*sin(thetafinalIC(1))*sin(phifinalIC(1))
						zfinalIC(1)= rfinalIC(1)*cos(thetafinalIC(1))

						if (abs(xfinalIC(1)) < 1d-6) then ! Set Cartesian coords according to B longitude
							xfinalIC(1)= 0d0
						end if

						if (abs(yfinalIC(1)) < 1d-6) then
							yfinalIC(1)= 0d0
						end if

						if (abs(zfinalIC(1)) < 1d-6) then
							zfinalIC(1)= 0d0
						end if

						! Note: Get final (q, p) value and let phid= phi to compare with initial input.

						qfinalIC(1)= (RE**2d0)*cos(thetafinalIC(1))/(rfinalIC(1)**2d0)
						pfinalIC(1)= rfinalIC(1)/(RE*(sin(thetafinalIC(1))**2d0))
						ellfinalIC(1)= 1d0+ 3d0*(cos(thetafinalIC(1)))**2d0

						! ----------------------------------------------------

						! DIAGNOSTIC FLAGS IN ABSOLUTE ERROR FOR QUARTIC DIPOLE POLYNOMIAL ROOTS:

						if (abs(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%pniICT(FAindIC)- &
							pfinalIC(1)) > 1d-12) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' BAD P VALUE', &
								' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
								', AND FAindIC= ', FAindIC, ' IN DIPOLE POLYNOMIAL SOLVER', &
								' SUBROUTINE' // achar(27) // '[0m.'
						end if

						if (abs(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniICT(FAindIC)- &
							qfinalIC(1)) > 1d-12) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' BAD Q VALUE', &
								' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
								', AND FAindIC= ', FAindIC, ' IN DIPOLE POLYNOMIAL SOLVER', &
								' SUBROUTINE' // achar(27) // '[0m.'
						end if

						! ----------------------------------------------------

						! Create nested derived data types
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%rfinalICT(FAindIC)= rfinalIC(1)
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%thetafinalICT(FAindIC)= thetafinalIC(1)
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%phifinalICT(FAindIC)= phifinalIC(1)
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%xfinalICT(FAindIC)= xfinalIC(1)
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%yfinalICT(FAindIC)= yfinalIC(1)
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%zfinalICT(FAindIC)= zfinalIC(1)
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qfinalICT(FAindIC)= qfinalIC(1)
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%pfinalICT(FAindIC)= pfinalIC(1)
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%ellfinalICT(FAindIC)= ellfinalIC(1)

						! ----------------------------------------------------

						! DIAGNOSTIC FLAGS FOR ALL qfinalICT VALUES WITHIN CONFIGURATION-SPACE GRID:

						if (SpecieT(s)%FluxTubeT(f)%qGLT(1, 1) <= 0) then
							if (((Qind == NqLB(1)) .and. &
								(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qfinalICT(FAindIC) < &
								SpecieT(s)%FluxTubeT(f)%qGLT(1, SpecieT(s)%FluxTubeT(f)%NqICAT(1)+ Qind- 1))) &
								.or. ((Qind == NqLB(1)) .and. &
								(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qfinalICT(FAindIC) > &
								SpecieT(s)%FluxTubeT(f)%qGHT(1, SpecieT(s)%FluxTubeT(f)%NqICAT(1)+ Qind- 1)))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qfinalICT= ', &
									SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qfinalICT(FAindIC), &
									' VALUE OUT OF CONFIG-SPACE GRID FOR SPECIE= ', s, ', FLUX TUBE= ', &
									f, ', Qind= ', Qind, ', AND FAindIC= ', FAindIC, ' IN DIPOLE', &
									' POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
							end if
							if (((Qind /= NqLB(1)) .and. &
								(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qfinalICT(FAindIC) <= &
								SpecieT(s)%FluxTubeT(f)%qGLT(1, SpecieT(s)%FluxTubeT(f)%NqICAT(1)+ Qind- 1))) &
								.or. ((Qind /= NqLB(1)) .and. &
								(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qfinalICT(FAindIC) > &
								SpecieT(s)%FluxTubeT(f)%qGHT(1, SpecieT(s)%FluxTubeT(f)%NqICAT(1)+ Qind- 1)))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qfinalICT= ', &
									SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qfinalICT(FAindIC), &
									' VALUE OUT OF CONFIG-SPACE GRID FOR SPECIE= ', s, ', FLUX TUBE= ', &
									f, ', Qind= ', Qind, ', AND FAindIC= ', FAindIC, ' IN DIPOLE', &
									' POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
							end if
						end if
						if (SpecieT(s)%FluxTubeT(f)%qGLT(1, 1) > 0) then
							if (((Qind == NqLB(1)) .and. &
								(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qfinalICT(FAindIC) > &
								SpecieT(s)%FluxTubeT(f)%qGLT(1, SpecieT(s)%FluxTubeT(f)%NqICAT(1)+ Qind- 1))) &
								.or. ((Qind == NqLB(1)) .and. &
								(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qfinalICT(FAindIC) < &
								SpecieT(s)%FluxTubeT(f)%qGHT(1, SpecieT(s)%FluxTubeT(f)%NqICAT(1)+ Qind- 1)))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qfinalICT= ', &
									SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qfinalICT(FAindIC), &
									' VALUE OUT OF CONFIG-SPACE GRID FOR SPECIE= ', s, ', FLUX TUBE= ', &
									f, ', Qind= ', Qind, ', AND FAindIC= ', FAindIC, ' IN DIPOLE', &
									' POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
							end if
							if (((Qind /= NqLB(1)) .and. &
								(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qfinalICT(FAindIC) >= &
								SpecieT(s)%FluxTubeT(f)%qGLT(1, SpecieT(s)%FluxTubeT(f)%NqICAT(1)+ Qind- 1))) &
								.or. ((Qind /= NqLB(1)) .and. &
								(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qfinalICT(FAindIC) < &
								SpecieT(s)%FluxTubeT(f)%qGHT(1, SpecieT(s)%FluxTubeT(f)%NqICAT(1)+ Qind- 1)))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qfinalICT= ', &
									SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qfinalICT(FAindIC), &
									' VALUE OUT OF CONFIG-SPACE GRID FOR SPECIE= ', s, ', FLUX TUBE= ', &
									f, ', Qind= ', Qind, ', AND FAindIC= ', FAindIC, ' IN DIPOLE', &
									' POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
							end if
						end if

						! ----------------------------------------------------

						! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSION, SIZES AND
						! FINITE VALUES:

						if ((isnan(real(rfinalIC(1))) .eqv. .true.) .or. &
							(size(rfinalIC(:)) /= 1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rfinalIC HAS', &
								' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND FAindIC= ', FAindIC, &
								' IN DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((thetafinalIC(1) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%thetafinalICT(FAindIC)) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%thetafinalICT(FAindIC))) &
							.eqv. .true.) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%thetafinalICT(:)) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' thetafinalICT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND FAindIC= ', FAindIC, &
								' IN DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((phifinalIC(1) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%phifinalICT(FAindIC)) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%phifinalICT(FAindIC))) &
							.eqv. .true.) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%phifinalICT(:)) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' phifinalICT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND FAindIC= ', FAindIC, &
								' IN DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((xfinalIC(1) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%xfinalICT(FAindIC)) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%xfinalICT(FAindIC))) &
							.eqv. .true.) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%xfinalICT(:)) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' xfinalICT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND FAindIC= ', FAindIC, &
								' IN DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((yfinalIC(1) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%yfinalICT(FAindIC)) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%yfinalICT(FAindIC))) &
							.eqv. .true.) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%yfinalICT(:)) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' yfinalICT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND FAindIC= ', FAindIC, &
								' IN DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((zfinalIC(1) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%zfinalICT(FAindIC)) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%zfinalICT(FAindIC))) &
							.eqv. .true.) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%zfinalICT(:)) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' zfinalICT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND FAindIC= ', FAindIC, &
								' IN DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((qfinalIC(1) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qfinalICT(FAindIC)) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qfinalICT(FAindIC))) &
							.eqv. .true.) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qfinalICT(:)) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qfinalICT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND FAindIC= ', FAindIC, &
								' IN DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((pfinalIC(1) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%pfinalICT(FAindIC)) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%pfinalICT(FAindIC))) &
							.eqv. .true.) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%pfinalICT(:)) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' pfinalICT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND FAindIC= ', FAindIC, &
								' IN DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((ellfinalIC(1) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%ellfinalICT(FAindIC)) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%ellfinalICT(FAindIC))) &
							.eqv. .true.) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%ellfinalICT(:)) /= &
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ellfinalICT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND FAindIC= ', FAindIC, &
								' IN DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
						end if

						! ----------------------------------------------------

					end do

					! ----------------------------------------------------

					! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES:

					if (abs(size(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%pniICT(:)) &
						- size(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%pfinalICT(:))) /= 0) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' BAD P SIZE FOR ', &
							' SPECIE= ', s, ', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
							' IN DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
					end if

					if (abs(size(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qniICT(:)) &
						- size(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qfinalICT(:))) /= 0) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' BAD Q SIZE FOR ', &
							' SPECIE= ', s, ', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
							' IN DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

				end do
			end do
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

		! PRINT TERMINAL EXPORT PARAMETERS:

		! do s= 1, Stot, 1
! 			do f= 1, SpecieT(s)%NfT(1), 1
! 				do Qind= NqLB(1), NqUB(1), 1
! 					if (rank == 0) then
! 						write(*, *) 'thetafinalICT= ', &
! 						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%thetafinalICT(:)
! 						write(*, *) 'phifinalICT= ', &
! 						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%phifinalICT(:)
! 						write(*, *) 'xfinalICT= ', &
! 						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%xfinalICT(:)
! 						write(*, *) 'yfinalICT= ', &
! 						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%yfinalICT(:)
! 						write(*, *) 'zfinalICT= ', &
! 						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%zfinalICT(:)
! 						write(*, *) 'qfinalICT= ', &
! 						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%qfinalICT(:)
! 						write(*, *) 'pfinalICT= ', &
! 						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%pfinalICT(:)
! 						write(*, *) 'ellfinalICT= ', &
! 						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%ellfinalICT(:)
! 					end if
! 				end do
! 			end do
! 		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				if (rank == 0) then
					call cpu_time(S5End)
					write(S5string, '(i10)')  nint(S5End)
					write(*, *) trim('%% 5- RANK= ' // adjustl(rankstring)) // &
						trim(', REAL-TIME= ' // adjustl(S5string)) // &
						trim(' s. DIPOLE POLYNOMIAL SOLVER COMPLETE %%')
				end if
			end do
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine DipolePolynomialSolverSub

end module DipolePolynomialSolver
