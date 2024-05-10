module Convection

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.11- CONVECTION:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use ConfigGridGenerator
use RK4DipolePolynomialSolver
use VelGridGenerator

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

	subroutine ConvectionSub

  	! ----------------------------------------------------

  	! UPDATE TIME-DEPENDENT PHASE-SPACE GRID WITHOUT CONVECTING FLUX-TUBES:

		if (n == 2) then
	  	if (SpecieT(s)%FluxTubeT(f)%CONVECTIONflagT(1) == 0) then
	  		do nn= 2, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1

	  			! ----------------------------------------------------

	  			SpecieT(s)%FluxTubeT(f)%ndatfacT(nn)= SpecieT(s)%FluxTubeT(f)%ndatfacT(1)

	  			SpecieT(s)%FluxTubeT(f)%d3xCLBT(nn)= SpecieT(s)%FluxTubeT(f)%d3xCLBT(1)
	  			SpecieT(s)%FluxTubeT(f)%sigmaLBT(nn)= SpecieT(s)%FluxTubeT(f)%sigmaLBT(1)
	  			SpecieT(s)%FluxTubeT(f)%d3xCUBT(nn)= SpecieT(s)%FluxTubeT(f)%d3xCUBT(1)
	  			SpecieT(s)%FluxTubeT(f)%sigmaUBT(nn)= SpecieT(s)%FluxTubeT(f)%sigmaUBT(1)
	  			SpecieT(s)%FluxTubeT(f)%LBNominalDensityT(nn)= SpecieT(s)%FluxTubeT(f)%LBNominalDensityT(1)
	  			SpecieT(s)%FluxTubeT(f)%UBNominalDensityT(nn)= SpecieT(s)%FluxTubeT(f)%UBNominalDensityT(1)
	  			SpecieT(s)%FluxTubeT(f)%ns0T(nn)= SpecieT(s)%FluxTubeT(f)%ns0T(1)

	  			! ----------------------------------------------------

	  			SpecieT(s)%FluxTubeT(f)%qGCT(nn, :)= SpecieT(s)%FluxTubeT(f)%qGCT(1, :)
	  			SpecieT(s)%FluxTubeT(f)%hqCT(nn, :)= SpecieT(s)%FluxTubeT(f)%hqCT(1, :)
	  			SpecieT(s)%FluxTubeT(f)%dpCT(nn, :)= SpecieT(s)%FluxTubeT(f)%dpCT(1, :)
	  			SpecieT(s)%FluxTubeT(f)%dqCT(nn, :)= SpecieT(s)%FluxTubeT(f)%dqCT(1, :)
	  			SpecieT(s)%FluxTubeT(f)%dphiCT(nn, :)= SpecieT(s)%FluxTubeT(f)%dphiCT(1, :)
	  			SpecieT(s)%FluxTubeT(f)%rGCT(nn, :)= SpecieT(s)%FluxTubeT(f)%rGCT(1, :)
	  			SpecieT(s)%FluxTubeT(f)%phiGCT(nn, :)= SpecieT(s)%FluxTubeT(f)%phiGCT(1, :)
	  			SpecieT(s)%FluxTubeT(f)%thetaGCT(nn, :)= SpecieT(s)%FluxTubeT(f)%thetaGCT(1, :)
	  			SpecieT(s)%FluxTubeT(f)%ellGCT(nn, :)= SpecieT(s)%FluxTubeT(f)%ellGCT(1, :)
	  			SpecieT(s)%FluxTubeT(f)%qGLT(nn, :)= SpecieT(s)%FluxTubeT(f)%qGLT(1, :)
	  			SpecieT(s)%FluxTubeT(f)%qGHT(nn, :)= SpecieT(s)%FluxTubeT(f)%qGHT(1, :)
	  			SpecieT(s)%FluxTubeT(f)%pGCT(nn, :)= SpecieT(s)%FluxTubeT(f)%pGCT(1, :)
	  			SpecieT(s)%FluxTubeT(f)%d3xCT(nn, :)= SpecieT(s)%FluxTubeT(f)%d3xCT(1, :)
	  			SpecieT(s)%FluxTubeT(f)%TsPerpT(nn, :)= SpecieT(s)%FluxTubeT(f)%TsPerpT(1, :)
	  			SpecieT(s)%FluxTubeT(f)%TsParT(nn, :)= SpecieT(s)%FluxTubeT(f)%TsParT(1, :)
	  			SpecieT(s)%FluxTubeT(f)%TsT(nn, :)= SpecieT(s)%FluxTubeT(f)%TsT(1, :)

	  			if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
	  				SpecieT(s)%FluxTubeT(f)%nsnormCNeutT(nn, :)= SpecieT(s)%FluxTubeT(f)%nsnormCNeutT(1, :)
	  			end if

	  			SpecieT(s)%FluxTubeT(f)%dsICRT(nn, :)= SpecieT(s)%FluxTubeT(f)%dsICRT(1, :)

	  			! ----------------------------------------------------

	  		end do
	  	end if
		end if

  	! ----------------------------------------------------

  	! CONVECT ALL FLUX-TUBES:

		! ----------------------------------------------------

  	if (SpecieT(s)%FluxTubeT(f)%CONVECTIONflagT(1) == 1) then
  		convnloop: do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
  			if ((n /= 1) .and. (nn /= 1) .and. &
  				(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn- 1))) then

  				! ----------------------------------------------------

					! UPDATE PHASE-SPACE GRID:

  				INITIALGRIDflag= 0

  				! Convect L-shells in (p, phi) with constant (q, vperp1, vperp2, vpar) grid:
  				if (SpecieT(s)%FluxTubeT(f)%DYNAMICGRIDflagT(1) == 0) then

  					ns0(1)= ns0IC ! Reference density at lower boundary ghost cell [m^-3]
  					SpecieT(s)%FluxTubeT(f)%nsnormfacT(nn)= nsnormfacIC ! Macroparticle normalization constant
  					zns0Neut(1)= zns0NeutIC ! Reference altitude of neutral density [km]
  					ns0Neut(1)= ns0NeutIC ! Reference neutral density [m^-3]
  					VperpsigmaFac(1)= VperpsigmaFacIC ! Number of MB standard deviations spanned by Vperp grid
  					VparsigmaFac(1)= VparsigmaFacIC ! Number of MB standard deviations spanned by Vpar grid
  					Ti(1)= TiIC ! Thermal ion temperature [K]
  					Te(1)= TeIC ! Electron temperature [K]
  					TNeut(1)= TNeutIC ! Neutral temperature [K]
  					qGA(1)= qGAIC ! Lower boundary q value
  					qGB(1)= qGBIC ! Upper boundary q value

  					Lshell(1)= LshellIC ! L-shell [RE]
  					phiLshell(1)= phiLshellIC ! invariant longitude [rads]

  				end if

  				! Convect L-shells in (p, phi) with dynamic (q, vperp1, vperp2, vpar) grid:
  				! Note: Conserve phase-space density by Louiville mapping
  				if (SpecieT(s)%FluxTubeT(f)%DYNAMICGRIDflagT(1) == 1) then

  					ns0(1)= ns0IC ! Reference density at lower boundary ghost cell [m^-3]
  					SpecieT(s)%FluxTubeT(f)%nsnormfacT(nn)= nsnormfacIC ! Macroparticle normalization constant
  					zns0Neut(1)= zns0NeutIC ! Reference altitude of neutral density [km]
  					ns0Neut(1)= ns0NeutIC ! Reference neutral density [m^-3]
  					VperpsigmaFac(1)= VperpsigmaFacIC ! Number of MB standard deviations spanned by Vperp grid
  					VparsigmaFac(1)= VparsigmaFacIC ! Number of MB standard deviations spanned by Vpar grid
  					Ti(1)= TiIC ! Thermal ion temperature [K]
  					Te(1)= TeIC ! Electron temperature [K]
  					TNeut(1)= TNeutIC ! Neutral temperature [K]
  					qGA(1)= qGAIC ! Lower boundary q value
  					qGB(1)= qGBIC ! Upper boundary q value

  					Lshell(1)= LshellIC ! L-shell [RE]
  					phiLshell(1)= phiLshellIC ! invariant longitude [rads]

  				end if

  				! ----------------------------------------------------

  				call ConfigGridGeneratorSub

  				! ----------------------------------------------------

					! FIXME Keep static velocity grid for now
  				!call VelGridGeneratorSub

  				! ----------------------------------------------------

  				exit convnloop

  			end if
  		end do convnloop
  	end if

		! ----------------------------------------------------

		! UPDATE ION POSITIONS AND VELOCITIES WITH BETRATRON ACCELERATION:

  	if (SpecieT(s)%FluxTubeT(f)%CONVECTIONflagT(1) == 1) then
  		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
  			if ((n /= 1) .and. (nn /= 1) .and. &
  				(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn- 1))) then

write(*, *) rank, 'C1'

					do j= 1, (NsTK(1)- dNsTK2(1)- dNsTK3(1)), 1

						! ----------------------------------------------------

						! UPDATE ION POSITIONS:
						! Note: On statistical time-steps use old values of AEAmagN, AGmagN, AEPmagN
						! (which are overridden on next computational time-step)

						call rsub(rConv(1), xN(j), yN(j), zN(j))
						call thetasub(thetaConv(1), zN(j), rConv(1))
						call qsub(qNp(1), rConv(1), thetaConv(1))
						call ellsub(ellConv(1), thetaConv(1))
						call Bmagsub(BmagConv(1), rConv(1), ellConv(1))

						call musub(muConv(1), SpecieT(s)%msT(1), BmagConv(1), VperpN(j))

						pNp(1)= SpecieT(s)%FluxTubeT(f)%pGCT(nn, 1) ! Update L-shell and longitude
						phiNp(1)= SpecieT(s)%FluxTubeT(f)%phiGCT(nn, 1)

						! ----------------------------------------------------

						call RK4DipolePolynomialSolverSub

						! ----------------------------------------------------

						xN(j)= xfinalRK4(1)
						yN(j)= yfinalRK4(1)
						zN(j)= zfinalRK4(1)

						! ----------------------------------------------------

						! UPDATE ION TRANSVERSE VELOCITIES:
						! Note: first adiabatic invariant is mu= ms*(Vperp**2d0)/(2d0*Bmag) s.t. Vperp= sqrt(2*Bmag*mu/ms)

						if (ENAflag(j) .eqv. .false.) then
							call rsub(rConv(1), xN(j), yN(j), zN(j))
							call thetasub(thetaConv(1), zN(j), rConv(1))
							call ellsub(ellConv(1), thetaConv(1))
							call Bmagsub(BmagConv(1), rConv(1), ellConv(1))

							! Preserve gyro-angle
							if ((Vperp1N(j) > 0d0) .and. (Vperp2N(j) > 0d0)) then
								GAConv(1)= atan(abs(Vperp2N(j))/abs(Vperp1N(j)))
							end if
							if ((Vperp1N(j) < 0d0) .and. (Vperp2N(j) > 0d0)) then
								GAConv(1)= pi- atan(abs(Vperp2N(j))/abs(Vperp1N(j)))
							end if
							if ((Vperp1N(j) < 0d0) .and. (Vperp2N(j) < 0d0)) then
								GAConv(1)= pi+ atan(abs(Vperp2N(j))/abs(Vperp1N(j)))
							end if
							if ((Vperp1N(j) > 0d0) .and. (Vperp2N(j) < 0d0)) then
								GAConv(1)= 2d0*pi- atan(abs(Vperp2N(j))/abs(Vperp1N(j)))
							end if

							! Compute new Vperp value with betatron acceleration
							VperpN(j)= sqrt(2d0*BmagConv(1)*muConv(1)/SpecieT(s)%msT(1))

							! ----------------------------------------------------

							! RECONSTRUCT Vperp1 and Vperp2 COMPONENTS WITH PRESERVED GYRO-ANGLE:

							if ((0d0 .lt. GAConv(1)) .and. (GAConv(1) .le. pi/2d0)) then
								Vperp1N(j)= abs(VperpN(j)*cos(GAConv(1)))
								Vperp2N(j)= abs(VperpN(j)*sin(GAConv(1)))

								! ----------------------------------------------------

								! ENSURE CORRECT SIGNS OF Vperp1 AND Vperp2 VALUES:

								if (Vperp1N(j) < 0d0) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' INCORRECT BETATRON SIGN OF Vperp1N= ', Vperp1N(j), &
										' FOR GYRO-ANGLE [rads]= ', GAConv(1), &
										' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
										', STATISTICAL TIME-STEP= ', nn, &
										', AND PARTICLE= ', j, ' CONVECTION SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								if (Vperp2N(j) < 0d0) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' INCORRECT BETATRON SIGN OF Vperp2N= ', Vperp2N(j), &
										' FOR GYRO-ANGLE [rads]= ', GAConv(1), &
										' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
										', STATISTICAL TIME-STEP= ', nn, &
										', AND PARTICLE= ', j, ' CONVECTION SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								! ----------------------------------------------------

							end if
							if ((pi/2d0 .lt. GAConv(1)) .and. (GAConv(1) .le. pi)) then
								Vperp1N(j)= -abs(VperpN(j)*cos(pi- GAConv(1)))
								Vperp2N(j)= abs(VperpN(j)*sin(pi- GAConv(1)))

								! ----------------------------------------------------

								! ENSURE CORRECT SIGNS OF Vperp1 AND Vperp2 VALUES:

								if (Vperp1N(j) > 0d0) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' INCORRECT BETATRON SIGN OF Vperp1N= ', Vperp1N(j), &
										' FOR GYRO-ANGLE [rads]= ', GAConv(1), &
										' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
										', STATISTICAL TIME-STEP= ', nn, &
										', AND PARTICLE= ', j, ' CONVECTION SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								if (Vperp2N(j) < 0d0) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' INCORRECT BETATRON SIGN OF Vperp2N= ', Vperp2N(j), &
										' FOR GYRO-ANGLE [rads]= ', GAConv(1), &
										' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
										', STATISTICAL TIME-STEP= ', nn, &
										', AND PARTICLE= ', j, ' CONVECTION SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								! ----------------------------------------------------

							end if
							if ((pi .lt. GAConv(1)) .and. (GAConv(1) .le. 3d0*pi/2d0)) then
								Vperp1N(j)= -abs(VperpN(j)*cos(GAConv(1)- pi))
								Vperp2N(j)= -abs(VperpN(j)*sin(GAConv(1)- pi))

								! ----------------------------------------------------

								! ENSURE CORRECT SIGNS OF Vperp1 AND Vperp2 VALUES:

								if (Vperp1N(j) > 0d0) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' INCORRECT BETATRON SIGN OF Vperp1N= ', Vperp1N(j), &
										' FOR GYRO-ANGLE [rads]= ', GAConv(1), &
										' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
										', STATISTICAL TIME-STEP= ', nn, &
										', AND PARTICLE= ', j, ' CONVECTION SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								if (Vperp2N(j) > 0d0) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' INCORRECT BETATRON SIGN OF Vperp2N= ', Vperp2N(j), &
										' FOR GYRO-ANGLE [rads]= ', GAConv(1), &
										' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
										', STATISTICAL TIME-STEP= ', nn, &
										', AND PARTICLE= ', j, ' CONVECTION SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								! ----------------------------------------------------

							end if
							if ((3d0*pi/2d0 .lt. GAConv(1)) .and. (GAConv(1) .le. 2d0*pi)) then
								Vperp1N(j)= abs(VperpN(j)*cos(2d0*pi- GAConv(1)))
								Vperp2N(j)= -abs(VperpN(j)*sin(2d0*pi- GAConv(1)))

								! ----------------------------------------------------

								! ENSURE CORRECT SIGNS OF Vperp1 AND Vperp2 VALUES:

								if (Vperp1N(j) < 0d0) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' INCORRECT BETATRON SIGN OF Vperp1N= ', Vperp1N(j), &
										' FOR GYRO-ANGLE [rads]= ', GAConv(1), &
										' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
										', STATISTICAL TIME-STEP= ', nn, &
										', AND PARTICLE= ', j, ' CONVECTION SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								if (Vperp2N(j) > 0d0) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' INCORRECT BETATRON SIGN OF Vperp2N= ', Vperp2N(j), &
										' FOR GYRO-ANGLE [rads]= ', GAConv(1), &
										' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
										', STATISTICAL TIME-STEP= ', nn, &
										', AND PARTICLE= ', j, ' CONVECTION SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								! ----------------------------------------------------

							end if
						end if

						! ----------------------------------------------------

						! UPDATE ION TRANSLATIONAL VELOCITIES:

						call rsub(rConv(1), xN(j), yN(j), zN(j))
						call thetasub(thetaConv(1), zN(j), rConv(1))

						if (ENAflag(j) .eqv. .false.) then
							phiConv(1)= SpecieT(s)%FluxTubeT(f)%phiGCT(nn, 1)
						else if ((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
							(ENAflag(j) .eqv. .true.)) then
							call phisub(phiConv(1), xN(j), yN(j))
						end if

						! ----------------------------------------------------

						! DIAGNOSTIC FLAGS FOR CONSISTENT PHI VALUE:

						if (ENAflag(j) .eqv. .false.) then
							if (phiConv(1) /= SpecieT(s)%FluxTubeT(f)%phiGCT(nn, 1)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' INCONSISTENT ION phiConv= ', phiConv(1), &
									' AND phiGCT VALUE= ', SpecieT(s)%FluxTubeT(f)%phiGCT(nnind, 1), &
									' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
									', STATISTICAL TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
									' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
							end if
							if (phiNp(1) /= SpecieT(s)%FluxTubeT(f)%phiGCT(nn, 1)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' INCONSISTENT ION phiNp= ', phiNp(1), &
									' AND phiGCT VALUE= ', SpecieT(s)%FluxTubeT(f)%phiGCT(nnind, 1), &
									' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
									', STATISTICAL TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
									' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
							end if
							if (pNp(1) /= SpecieT(s)%FluxTubeT(f)%pGCT(nn, 1)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' INCONSISTENT ION pNp= ', pNp(1), &
									' AND phiGCT VALUE= ', SpecieT(s)%FluxTubeT(f)%phiGCT(nnind, 1), &
									' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
									', STATISTICAL TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
									' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
							end if
						end if

						! ----------------------------------------------------

						call ellsub(ellConv(1), thetaConv(1))

						! ----------------------------------------------------

						! RE-ALIGN VELOCITY VECTOR INTO LOCAL DIPOLE COORDINATES:

						VparConv(1)= 3d0*cos(thetaConv(1))*sin(thetaConv(1))*(VxN(j)* &
							cos(phiConv(1))+ VyN(j)*sin(phiConv(1)))/sqrt(ellConv(1))+ &
							VzN(j)*(3d0*(cos(thetaConv(1))**2d0)- 1d0)/sqrt(ellConv(1))

						! ----------------------------------------------------

						if (ENAflag(j) .eqv. .false.) then
							VpConv(1)= 0d0
							VphiConv(1)= 0d0
						else if ((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
							(ENAflag(j) .eqv. .true.)) then
							VpConv(1)= abs((1d0- 3d0*(cos(thetaConv(1))**2d0))*(VxN(j)* &
								cos(phiConv(1))+ VyN(j)*sin(phiConv(1)))/sqrt(ellConv(1))+ 3d0*VzN(j)* &
								cos(thetaConv(1))*sin(thetaConv(1))/sqrt(ellConv(1)))
							VphiConv(1)= -VxN(j)*sin(phiConv(1))+ VyN(j)*cos(phiConv(1))
						end if

						! ----------------------------------------------------

						! COMPUTE NEW FIELD-ALIGNED CARTESIAN VELOCITY COMPONENTS:

						VxN(j)= cos(phiConv(1))*(3d0*VparConv(1)*cos(thetaConv(1))*sin(thetaConv(1))+ &
							VpConv(1)*(1d0- 3d0*(cos(thetaConv(1))**2d0)))/sqrt(ellConv(1))- &
							VphiConv(1)*sin(phiConv(1))
						VyN(j)= sin(phiConv(1))*(3d0*VparConv(1)*cos(thetaConv(1))*sin(thetaConv(1))+ &
							VpConv(1)*(1d0- 3d0*(cos(thetaConv(1))**2d0)))/sqrt(ellConv(1))+ &
							VphiConv(1)*cos(phiConv(1))
						VzN(j)= VparConv(1)*(3d0*(cos(thetaConv(1))**2d0)- 1d0)/sqrt(ellConv(1))+ &
							3d0*VpConv(1)*cos(thetaConv(1))*sin(thetaConv(1))/sqrt(ellConv(1))

						! ----------------------------------------------------

						! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

						if ((isnan(real(VxN(j))) .eqv. .true.) .or. &
							(size(VxN(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VxN= ', VxN(j), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', STATISTICAL TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(VyN(j))) .eqv. .true.) .or. &
							(size(VyN(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VyN= ', VyN(j), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', STATISTICAL TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(VzN(j))) .eqv. .true.) .or. &
							(size(VzN(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VzN= ', VzN(j), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', STATISTICAL TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						! ----------------------------------------------------

					end do

					write(*, *) rank, 'C2'

				end if
			end do
		end if

		! ----------------------------------------------------

  	!do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
  	!	if ((n /= 1) .and. (nn /= 1) .and. (n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn))) then
  	! update grid with conservation of f and particle number,
  	! update LBNominalDensityT and UBNominalDensityT
  	! update all current particle positions and velocities with convection and with betatron acceleration
  	!do j= 1, (NsTK(1)- dNsTK2(1)- dNsTK3(1)), 1
  	! xN(j), yN(j), zN(j), Vperp1N(j), Vperp2N(j), VperpN(j), VxN(j), VyN(j), VzN(j) for ions and ENAs.
  	! Note: all new LB and UB injected particles are reset on correct Lshell in KineticSolverB given
  	! SpecieT(s)%FluxTubeT(f)%pGCT(nn, Qind) and SpecieT(s)%FluxTubeT(f)%phiGCT(nn, Qind)
  	!1- Conserve f and Update GRID
  	!2- Update LB and UB densities and ensure robust statistics
  	!3- Translate all particle positions in (p, q, phi)
  	!4- Update Vperp values with betatron ACCELERATION
  	!5- On new field line, update AEAmagN(j), AGmagN(j), AEPmagN(j)
  	!6- In KineticUpdateA, KineticUpdateB, KineticUpdateC add parallel centrifugal acceleration

		! ----------------------------------------------------

	end subroutine ConvectionSub

end module Convection
