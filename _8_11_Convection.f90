module Convection

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.11- CONVECTION:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use ConfigGridGenerator
use RK4DipolePolynomialSolver
use ConfigGridGeneratorB
use VelGridGenerator

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

	subroutine ConvectionSub

  	! ----------------------------------------------------

  	! UPDATE GRID PARAMETERS WITHOUT CONVECTION:
		! Note: These are the grid parameters out put from ConfigGridGenerator.f90.
		! Note: Leave all grid parameters in (nn) for BoundaryConditions.f90, and KineticSolverB.f90

  	if ((SpecieT(s)%FluxTubeT(f)%CONVECTIONflagT(1) == 0) .and. (n == 2)) then
  		do nn= 2, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1

  			! ----------------------------------------------------

  			SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)= SpecieT(s)%FluxTubeT(f)%ndatfacGT(1)
				SpecieT(s)%FluxTubeT(f)%nsnormCLBGT(nn)= SpecieT(s)%FluxTubeT(f)%nsnormCLBGT(1)
				SpecieT(s)%FluxTubeT(f)%nsnormCUBGT(nn)= SpecieT(s)%FluxTubeT(f)%nsnormCUBGT(1)

  			! ----------------------------------------------------

  			SpecieT(s)%FluxTubeT(f)%qGCGT(nn, :)= SpecieT(s)%FluxTubeT(f)%qGCGT(1, :)
  			SpecieT(s)%FluxTubeT(f)%hqCGT(nn, :)= SpecieT(s)%FluxTubeT(f)%hqCGT(1, :)
  			SpecieT(s)%FluxTubeT(f)%dpCGT(nn, :)= SpecieT(s)%FluxTubeT(f)%dpCGT(1, :)
  			SpecieT(s)%FluxTubeT(f)%dqCGT(nn, :)= SpecieT(s)%FluxTubeT(f)%dqCGT(1, :)
  			SpecieT(s)%FluxTubeT(f)%dphiCGT(nn, :)= SpecieT(s)%FluxTubeT(f)%dphiCGT(1, :)
  			SpecieT(s)%FluxTubeT(f)%rGCGT(nn, :)= SpecieT(s)%FluxTubeT(f)%rGCGT(1, :)
  			SpecieT(s)%FluxTubeT(f)%phiGCGT(nn, :)= SpecieT(s)%FluxTubeT(f)%phiGCGT(1, :)
  			SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, :)= SpecieT(s)%FluxTubeT(f)%thetaGCGT(1, :)
  			SpecieT(s)%FluxTubeT(f)%ellGCGT(nn, :)= SpecieT(s)%FluxTubeT(f)%ellGCGT(1, :)
  			SpecieT(s)%FluxTubeT(f)%qGLGT(nn, :)= SpecieT(s)%FluxTubeT(f)%qGLGT(1, :)
  			SpecieT(s)%FluxTubeT(f)%qGHGT(nn, :)= SpecieT(s)%FluxTubeT(f)%qGHGT(1, :)
  			SpecieT(s)%FluxTubeT(f)%pGCGT(nn, :)= SpecieT(s)%FluxTubeT(f)%pGCGT(1, :)
  			SpecieT(s)%FluxTubeT(f)%d3xCGT(nn, :)= SpecieT(s)%FluxTubeT(f)%d3xCGT(1, :)
  			SpecieT(s)%FluxTubeT(f)%TsPerpGT(nn, :)= SpecieT(s)%FluxTubeT(f)%TsPerpGT(1, :)
  			SpecieT(s)%FluxTubeT(f)%TsParGT(nn, :)= SpecieT(s)%FluxTubeT(f)%TsParGT(1, :)
  			SpecieT(s)%FluxTubeT(f)%TsGT(nn, :)= SpecieT(s)%FluxTubeT(f)%TsGT(1, :)
				SpecieT(s)%FluxTubeT(f)%TeGT(nn, :)= SpecieT(s)%FluxTubeT(f)%TeGT(1, :)

  			if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
  				SpecieT(s)%FluxTubeT(f)%nsnormCNeutGT(nn, :)= SpecieT(s)%FluxTubeT(f)%nsnormCNeutGT(1, :)
  			end if

  			SpecieT(s)%FluxTubeT(f)%dsICRGT(nn, :)= SpecieT(s)%FluxTubeT(f)%dsICRGT(1, :)

  			! ----------------------------------------------------

  		end do
  	end if

  	! ----------------------------------------------------

  	! UPDATE GRID WITH CONVECTION:

		! ----------------------------------------------------

  	if (SpecieT(s)%FluxTubeT(f)%CONVECTIONflagT(1) == 1) then
  		convnloop1: do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
  			if ((n /= 1) .and. (nn /= 1) .and. &
  				(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1)))) then

  				! ----------------------------------------------------

  				INITIALGRIDflag= 0

  				! Convect L-shells in (p, phi) with constant (q, vperp1, vperp2, vpar) grid:
  				if (SpecieT(s)%FluxTubeT(f)%DYNAMICGRIDflagT(1) == 0) then

  					ns0(1)= ns0IC ! Reference density at lower boundary ghost cell [m^-3]
  					nsnormfac(1)= nsnormfacIC ! Macroparticle normalization constant
  					zns0Neut(1)= zns0NeutIC ! Reference altitude of neutral density [km]
  					ns0Neut(1)= ns0NeutIC ! Reference neutral density [m^-3]
  					VperpsigmaFac(1)= VperpsigmaFacIC ! Number of MB standard deviations spanned by Vperp grid
  					VparsigmaFac(1)= VparsigmaFacIC ! Number of MB standard deviations spanned by Vpar grid
  					Ti(1)= TiIC ! Thermal ion temperature [K]
  					Te(1)= TeIC ! Electron temperature [K]
  					TNeut(1)= TNeutIC ! Neutral temperature [K]
  					qGA(1)= qGAIC ! Lower boundary q value
  					qGB(1)= qGBIC ! Upper boundary q value

						! FIXME insert analytical or empirical or MAGE potential solver here. include a KAOS plasmasphere?

  					Lshell(1)= Lshell(1)+ 0d0 ! L-shell [RE]
  					phiLshell(1)= phiLshell(1)+ 0d0 ! invariant longitude [rads]

  				end if

  				! Convect L-shells in (p, phi) with dynamic (q, vperp1, vperp2, vpar) grid:
  				! Note: Conserve phase-space density by Louiville mapping
  				if (SpecieT(s)%FluxTubeT(f)%DYNAMICGRIDflagT(1) == 1) then

  					ns0(1)= ns0IC ! Reference density at lower boundary ghost cell [m^-3]
  					nsnormfac(1)= nsnormfacIC ! Macroparticle normalization constant
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

					NNtp(1)= SpecieT(s)%FluxTubeT(f)%NtT(1)/SpecieT(s)%FluxTubeT(f)%Q0ndatfacT(1)

					! Time-step interval for moment computation (must be > 2 and excludes initial time-step)
					ndatfacG(1)= SpecieT(s)%FluxTubeT(f)%NtT(1)/NNtp(1)

					! ----------------------------------------------------

					call ConfigGridGeneratorSubB

					! ----------------------------------------------------

					SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)= ndatfacG(1)
					SpecieT(s)%FluxTubeT(f)%nsnormCLBGT(nn)= nsnormCLBG(1)
					SpecieT(s)%FluxTubeT(f)%nsnormCUBGT(nn)= nsnormCUBG(1)

					! ----------------------------------------------------

					! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

					if (isnan(real(SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn))) .eqv. .true.) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ndatfacGT= ', &
							SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn), &
							' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
							' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if (isnan(real(SpecieT(s)%FluxTubeT(f)%nsnormCLBGT(nn))) .eqv. .true.) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' nsnormCLBGT= ', &
							SpecieT(s)%FluxTubeT(f)%nsnormCLBGT(nn), &
							' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
							' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if (isnan(real(SpecieT(s)%FluxTubeT(f)%nsnormCUBGT(nn))) .eqv. .true.) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' nsnormCUBGT= ', &
							SpecieT(s)%FluxTubeT(f)%nsnormCUBGT(nn), &
							' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
							' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

					! RE-INDEX CONFIGURATION SPACE GRID FOR A NON-COMPUTATIONAL LOWER BOUNDARY GHOST CELL:

					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

						SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind)= qGC0(Qind+ 1)
						SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)= hqC0(Qind+ 1)
						SpecieT(s)%FluxTubeT(f)%dpCGT(nn, Qind)= dpC0(Qind+ 1)
						SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind)= dqC0(Qind+ 1)
						SpecieT(s)%FluxTubeT(f)%dphiCGT(nn, Qind)= dphiC0(Qind+ 1)
						SpecieT(s)%FluxTubeT(f)%rGCGT(nn, Qind)= rGC0(Qind+ 1)
						SpecieT(s)%FluxTubeT(f)%phiGCGT(nn, Qind)= phiGC0(Qind+ 1)
						SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind)= thetaGC0(Qind+ 1)
						SpecieT(s)%FluxTubeT(f)%ellGCGT(nn, Qind)= ellGC0(Qind+ 1)
						SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind)= qGL0(Qind+ 1)
						SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)= qGH0(Qind+ 1)
						SpecieT(s)%FluxTubeT(f)%pGCGT(nn, Qind)= pGC0(Qind+ 1)
						SpecieT(s)%FluxTubeT(f)%d3xCGT(nn, Qind)= d3xC0(Qind+ 1)
						SpecieT(s)%FluxTubeT(f)%TsPerpGT(nn, Qind)= TsPerp0(Qind+ 1)
						SpecieT(s)%FluxTubeT(f)%TsParGT(nn, Qind)= TsPar0(Qind+ 1)
						SpecieT(s)%FluxTubeT(f)%TsGT(nn, Qind)= Ts0(Qind+ 1)
						SpecieT(s)%FluxTubeT(f)%TeGT(nn, Qind)= Te0(Qind+ 1)

						if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
							SpecieT(s)%FluxTubeT(f)%nsnormCNeutGT(nn, Qind)= nsnormCNeut0(Qind+ 1)
						end if

						! ----------------------------------------------------

						! Compute field line arc length of BBELF wave-heating region
						SpecieT(s)%FluxTubeT(f)%dsICRGT(nn, Qind)= &
							SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind)

						! ----------------------------------------------------

						! DIAGNOSTIC FLAGS FOR CONSISTENT PHI VALUE:

						if (SpecieT(s)%FluxTubeT(f)%phiGCGT(nn, Qind) /= phiGC0(1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT INITIAL phiGCGT VALUE= ', &
								SpecieT(s)%FluxTubeT(f)%phiGCGT(nn, Qind), ' AND phiGC0T VALUE= ', phiGC0(Qind), &
								' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND Q CELL= ', Qind, &
								' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
						end if

						if (SpecieT(s)%FluxTubeT(f)%pGCGT(nn, Qind) /= pGC0(1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT INITIAL pGCGT VALUE= ', &
								SpecieT(s)%FluxTubeT(f)%pGCGT(nn, Qind), ' AND phiGC0T VALUE= ', pGC0(Qind), &
								' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND Q CELL= ', Qind, &
								' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
						end if

						! ----------------------------------------------------

					end do

					! ----------------------------------------------------

					! FIXME Keep static velocity grid for now but add dynamic grid later
  				!call VelGridGeneratorSub

					! ----------------------------------------------------

					exit convnloop1

				end if
			end do convnloop1

			! ----------------------------------------------------

  	end if

		! ----------------------------------------------------

		! UPDATE ION POSITIONS AND VELOCITIES WITH CONVECTION AND BETATRON ACCELERATION:

		! ----------------------------------------------------

		if (SpecieT(s)%FluxTubeT(f)%CONVECTIONflagT(1) == 1) then
		  convnloopxv: do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
		    if ((n /= 1) .and. (nn /= 1) .and. &
		      (n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1)))) then

					! ----------------------------------------------------

		      do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1

		        ! ----------------------------------------------------

		        ! UPDATE ION POSITIONS:
		        ! Note: AEAmagN, AGmagN, AEPmagN are updated on new grid after Convection.f90 in KineticSolver.f90

						! ----------------------------------------------------

						! Compute magnetic moments before convection
						! with current Vperp value and Bmag
		        call rsub(rConv1(1), xN(j), yN(j), zN(j))
		        call thetasub(thetaConv1(1), zN(j), rConv1(1))
		        call qsub(qNp(1), rConv1(1), thetaConv1(1))
		        call ellsub(ellConv1(1), thetaConv1(1))
		        call Bmagsub(BmagConv1(1), rConv1(1), ellConv1(1))
		        call musub(muConv1(1), SpecieT(s)%msT(1), BmagConv1(1), VperpN(j))

						! Compute new dipole positions after convection
		        pNp(1)= SpecieT(s)%FluxTubeT(f)%pGCGT(nn, 1) ! Update L-shell and longitude
		        phiNp(1)= SpecieT(s)%FluxTubeT(f)%phiGCGT(nn, 1)

		        ! ----------------------------------------------------

		        call RK4DipolePolynomialSolverSub

		        ! ----------------------------------------------------

						! Compute new Cartesian positions after convection
		        xN(j)= xfinalRK4(1)
		        yN(j)= yfinalRK4(1)
		        zN(j)= zfinalRK4(1)

						! ----------------------------------------------------

						! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

						if ((isnan(real(rConv1(1))) .eqv. .true.) .or. &
							(size(rConv1(:)) /= 1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rConv1= ', rConv1(1), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(thetaConv1(1))) .eqv. .true.) .or. &
							(size(thetaConv1(:)) /= 1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' thetaConv1= ', thetaConv1(1), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(qNp(1))) .eqv. .true.) .or. &
							(size(qNp(:)) /= 1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qNp= ', qNp(1), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(ellConv1(1))) .eqv. .true.) .or. &
							(size(ellConv1(:)) /= 1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ellConv1= ', ellConv1(1), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(BmagConv1(1))) .eqv. .true.) .or. &
							(size(BmagConv1(:)) /= 1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' BmagConv1= ', BmagConv1(1), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(muConv1(1))) .eqv. .true.) .or. &
							(size(muConv1(:)) /= 1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' muConv1= ', muConv1(1), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(pNp(1))) .eqv. .true.) .or. &
							(size(pNp(:)) /= 1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' pNp= ', pNp(1), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(phiNp(1))) .eqv. .true.) .or. &
							(size(phiNp(:)) /= 1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' phiNp= ', phiNp(1), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(xN(j))) .eqv. .true.) .or. &
							(size(xN(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' xN= ', xN(j), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(yN(j))) .eqv. .true.) .or. &
							(size(yN(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' yN= ', yN(j), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(zN(j))) .eqv. .true.) .or. &
							(size(zN(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' zN= ', zN(j), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						! ----------------------------------------------------

		        ! DIAGNOSTIC FLAGS FOR CONSISTENT PHI VALUE:

		        if (ENAflag(j) .eqv. .false.) then
		          if (phiNp(1) /= SpecieT(s)%FluxTubeT(f)%phiGCGT(nn, 1)) then
		            write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
		              ' INCONSISTENT ION phiNp= ', phiNp(1), &
		              ' AND phiGCGT VALUE= ', SpecieT(s)%FluxTubeT(f)%phiGCGT(nnind, 1), &
		              ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
		              ', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
		              ' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
		          end if
		          if (pNp(1) /= SpecieT(s)%FluxTubeT(f)%pGCGT(nn, 1)) then
		            write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
		              ' INCONSISTENT ION pNp= ', pNp(1), &
		              ' AND phiGCGT VALUE= ', SpecieT(s)%FluxTubeT(f)%phiGCGT(nnind, 1), &
		              ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
		              ', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
		              ' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
		          end if
		        end if

		        ! ----------------------------------------------------

		        ! UPDATE ION TRANSVERSE VELOCITIES WITH BETATRON ACCELERATION:
		        ! Note: first adiabatic invariant is mu= ms*(Vperp**2d0)/(2d0*Bmag) s.t. Vperp= sqrt(2*Bmag*mu/ms)
						! Note: VxN, VyN, and VzN remain unchanged

						! ----------------------------------------------------

		        if (ENAflag(j) .eqv. .false.) then

							! Compute Vperp values after convection
							! with new Bmag and conserved magnetic moments and gyro-angles.
		          call rsub(rConv2(1), xN(j), yN(j), zN(j))
		          call thetasub(thetaConv2(1), zN(j), rConv2(1))
		          call ellsub(ellConv2(1), thetaConv2(1))
		          call Bmagsub(BmagConv2(1), rConv2(1), ellConv2(1))

		          ! Compute gyro-angle before convection
		          if ((Vperp1N(j) > 0d0) .and. (Vperp2N(j) > 0d0)) then
		            GAConv1(1)= atan(abs(Vperp2N(j))/abs(Vperp1N(j)))
		          end if
		          if ((Vperp1N(j) < 0d0) .and. (Vperp2N(j) > 0d0)) then
		            GAConv1(1)= pi- atan(abs(Vperp2N(j))/abs(Vperp1N(j)))
		          end if
		          if ((Vperp1N(j) < 0d0) .and. (Vperp2N(j) < 0d0)) then
		            GAConv1(1)= pi+ atan(abs(Vperp2N(j))/abs(Vperp1N(j)))
		          end if
		          if ((Vperp1N(j) > 0d0) .and. (Vperp2N(j) < 0d0)) then
		            GAConv1(1)= 2d0*pi- atan(abs(Vperp2N(j))/abs(Vperp1N(j)))
		          end if

							! ----------------------------------------------------

							! COMPUTE NEW Vperp WITH BETATRON ACCELERATION:

		          VperpN(j)= sqrt(2d0*BmagConv2(1)*muConv1(1)/SpecieT(s)%msT(1))

							! ----------------------------------------------------

							! ENSURE NEW MAGNETIC MOMENT IS CONSERVED:

							call musub(muConv2(1), SpecieT(s)%msT(1), BmagConv2(1), VperpN(j))

							if (abs(muConv2(1)- muConv1(1)) .gt. 1e-9) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' ION MAGNETIC MOMENT IS NOT CONSERVED, mu1= ', muConv1(1), &
									', AND mu2= ', muConv2(1), ' FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', MASTER TIME-STEP= ', nn, &
									', AND PARTICLE= ', j, ' CONVECTION SUBROUTINE' &
									// achar(27) // '[0m.'
							end if

		          ! ----------------------------------------------------

		          ! RECONSTRUCT Vperp1 and Vperp2 COMPONENTS WITH PRESERVED GYRO-ANGLE:

		          if ((0d0 .lt. GAConv1(1)) .and. (GAConv1(1) .le. pi/2d0)) then
		            Vperp1N(j)= abs(VperpN(j)*cos(GAConv1(1)))
		            Vperp2N(j)= abs(VperpN(j)*sin(GAConv1(1)))

		            ! ----------------------------------------------------

		            ! ENSURE CORRECT SIGNS OF Vperp1 AND Vperp2 VALUES:

		            if (Vperp1N(j) < 0d0) then
		              write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
		                ' INCORRECT BETATRON SIGN OF Vperp1N= ', Vperp1N(j), &
		                ' FOR GYRO-ANGLE [rads]= ', GAConv1(1), &
		                ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
		                ', MASTER TIME-STEP= ', nn, &
		                ', AND PARTICLE= ', j, ' CONVECTION SUBROUTINE' &
		                // achar(27) // '[0m.'
		            end if

		            if (Vperp2N(j) < 0d0) then
		              write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
		                ' INCORRECT BETATRON SIGN OF Vperp2N= ', Vperp2N(j), &
		                ' FOR GYRO-ANGLE [rads]= ', GAConv1(1), &
		                ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
		                ', MASTER TIME-STEP= ', nn, &
		                ', AND PARTICLE= ', j, ' CONVECTION SUBROUTINE' &
		                // achar(27) // '[0m.'
		            end if

		            ! ----------------------------------------------------

		          end if
		          if ((pi/2d0 .lt. GAConv1(1)) .and. (GAConv1(1) .le. pi)) then
		            Vperp1N(j)= -abs(VperpN(j)*cos(pi- GAConv1(1)))
		            Vperp2N(j)= abs(VperpN(j)*sin(pi- GAConv1(1)))

		            ! ----------------------------------------------------

		            ! ENSURE CORRECT SIGNS OF Vperp1 AND Vperp2 VALUES:

		            if (Vperp1N(j) > 0d0) then
		              write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
		                ' INCORRECT BETATRON SIGN OF Vperp1N= ', Vperp1N(j), &
		                ' FOR GYRO-ANGLE [rads]= ', GAConv1(1), &
		                ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
		                ', MASTER TIME-STEP= ', nn, &
		                ', AND PARTICLE= ', j, ' CONVECTION SUBROUTINE' &
		                // achar(27) // '[0m.'
		            end if

		            if (Vperp2N(j) < 0d0) then
		              write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
		                ' INCORRECT BETATRON SIGN OF Vperp2N= ', Vperp2N(j), &
		                ' FOR GYRO-ANGLE [rads]= ', GAConv1(1), &
		                ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
		                ', MASTER TIME-STEP= ', nn, &
		                ', AND PARTICLE= ', j, ' CONVECTION SUBROUTINE' &
		                // achar(27) // '[0m.'
		            end if

		            ! ----------------------------------------------------

		          end if
		          if ((pi .lt. GAConv1(1)) .and. (GAConv1(1) .le. 3d0*pi/2d0)) then
		            Vperp1N(j)= -abs(VperpN(j)*cos(GAConv1(1)- pi))
		            Vperp2N(j)= -abs(VperpN(j)*sin(GAConv1(1)- pi))

		            ! ----------------------------------------------------

		            ! ENSURE CORRECT SIGNS OF Vperp1 AND Vperp2 VALUES:

		            if (Vperp1N(j) > 0d0) then
		              write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
		                ' INCORRECT BETATRON SIGN OF Vperp1N= ', Vperp1N(j), &
		                ' FOR GYRO-ANGLE [rads]= ', GAConv1(1), &
		                ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
		                ', MASTER TIME-STEP= ', nn, &
		                ', AND PARTICLE= ', j, ' CONVECTION SUBROUTINE' &
		                // achar(27) // '[0m.'
		            end if

		            if (Vperp2N(j) > 0d0) then
		              write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
		                ' INCORRECT BETATRON SIGN OF Vperp2N= ', Vperp2N(j), &
		                ' FOR GYRO-ANGLE [rads]= ', GAConv1(1), &
		                ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
		                ', MASTER TIME-STEP= ', nn, &
		                ', AND PARTICLE= ', j, ' CONVECTION SUBROUTINE' &
		                // achar(27) // '[0m.'
		            end if

		            ! ----------------------------------------------------

		          end if
		          if ((3d0*pi/2d0 .lt. GAConv1(1)) .and. (GAConv1(1) .le. 2d0*pi)) then
		            Vperp1N(j)= abs(VperpN(j)*cos(2d0*pi- GAConv1(1)))
		            Vperp2N(j)= -abs(VperpN(j)*sin(2d0*pi- GAConv1(1)))

		            ! ----------------------------------------------------

		            ! ENSURE CORRECT SIGNS OF Vperp1 AND Vperp2 VALUES:

		            if (Vperp1N(j) < 0d0) then
		              write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
		                ' INCORRECT BETATRON SIGN OF Vperp1N= ', Vperp1N(j), &
		                ' FOR GYRO-ANGLE [rads]= ', GAConv1(1), &
		                ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
		                ', MASTER TIME-STEP= ', nn, &
		                ', AND PARTICLE= ', j, ' CONVECTION SUBROUTINE' &
		                // achar(27) // '[0m.'
		            end if

		            if (Vperp2N(j) > 0d0) then
		              write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
		                ' INCORRECT BETATRON SIGN OF Vperp2N= ', Vperp2N(j), &
		                ' FOR GYRO-ANGLE [rads]= ', GAConv1(1), &
		                ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
		                ', MASTER TIME-STEP= ', nn, &
		                ', AND PARTICLE= ', j, ' CONVECTION SUBROUTINE' &
		                // achar(27) // '[0m.'
		            end if

		            ! ----------------------------------------------------

		          end if
		        end if

						! ----------------------------------------------------

		        ! UPDATE ION TRANSLATIONAL VELOCITIES:
						! Note: Update VxN, VyN, and VzN on new field line

						VparConv2(1)= VparConvSign(j)*abs(sqrt(VxN(j)**2d0+ VyN(j)**2d0+ VzN(j)**2d0))

 						if (ENAflag(j) .eqv. .false.) then
 							VpConv2(1)= 0d0
 							VphiConv2(1)= 0d0
 						else if ((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
 							(ENAflag(j) .eqv. .true.)) then
 							VpConv2(1)= abs((1d0- 3d0*(cos(thetaConv2(1))**2d0))*(VxN(j)* &
 								cos(phiNp(1))+ VyN(j)*sin(phiNp(1)))/sqrt(ellConv2(1))+ 3d0*VzN(j)* &
 								cos(thetaConv2(1))*sin(thetaConv2(1))/sqrt(ellConv2(1)))
 							VphiConv2(1)= -VxN(j)*sin(phiNp(1))+ VyN(j)*cos(phiNp(1))
 						end if

 						! ----------------------------------------------------

 						! COMPUTE NEW FIELD-ALIGNED CARTESIAN VELOCITY COMPONENTS:

 						VxN(j)= cos(phiNp(1))*(3d0*VparConv2(1)*cos(thetaConv2(1))*sin(thetaConv2(1))+ &
 							VpConv2(1)*(1d0- 3d0*(cos(thetaConv2(1))**2d0)))/sqrt(ellConv2(1))- &
 							VphiConv2(1)*sin(phiNp(1))
 						VyN(j)= sin(phiNp(1))*(3d0*VparConv2(1)*cos(thetaConv2(1))*sin(thetaConv2(1))+ &
 							VpConv2(1)*(1d0- 3d0*(cos(thetaConv2(1))**2d0)))/sqrt(ellConv2(1))+ &
 							VphiConv2(1)*cos(phiNp(1))
 						VzN(j)= VparConv2(1)*(3d0*(cos(thetaConv2(1))**2d0)- 1d0)/sqrt(ellConv2(1))+ &
 							3d0*VpConv2(1)*cos(thetaConv2(1))*sin(thetaConv2(1))/sqrt(ellConv2(1))

						! ----------------------------------------------------

		        ! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

						if ((isnan(real(rConv2(1))) .eqv. .true.) .or. &
							(size(rConv2(:)) /= 1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rConv2= ', rConv2(1), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(thetaConv2(1))) .eqv. .true.) .or. &
							(size(thetaConv2(:)) /= 1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' thetaConv2= ', thetaConv2(1), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(ellConv2(1))) .eqv. .true.) .or. &
							(size(ellConv2(:)) /= 1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ellConv2= ', ellConv2(1), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(BmagConv2(1))) .eqv. .true.) .or. &
							(size(BmagConv2(:)) /= 1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' BmagConv2= ', BmagConv2(1), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(GAConv1(1))) .eqv. .true.) .or. &
							(size(GAConv1(:)) /= 1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' GAConv1= ', GAConv1(1), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(VperpN(j))) .eqv. .true.) .or. &
							(size(VperpN(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VperpN= ', VperpN(j), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(Vperp1N(j))) .eqv. .true.) .or. &
							(size(Vperp1N(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp1N= ', Vperp1N(j), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(Vperp2N(j))) .eqv. .true.) .or. &
							(size(Vperp2N(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp2N= ', Vperp2N(j), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(VparConv2(1))) .eqv. .true.) .or. &
							(size(VparConv2(:)) /= 1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VparConv2= ', VparConv2(1), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(VpConv2(1))) .eqv. .true.) .or. &
							(size(VpConv2(:)) /= 1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VpConv2= ', VpConv2(1), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(VphiConv2(1))) .eqv. .true.) .or. &
							(size(VphiConv2(:)) /= 1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VphiConv2= ', VphiConv2(1), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(VxN(j))) .eqv. .true.) .or. &
							(size(VxN(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VxN= ', VxN(j), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(VyN(j))) .eqv. .true.) .or. &
							(size(VyN(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VyN= ', VyN(j), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((isnan(real(VzN(j))) .eqv. .true.) .or. &
							(size(VzN(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VzN= ', VzN(j), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', MASTER TIME-STEP= ', nn, ', AND PARTICLE= ', j, &
								' IN CONVECTION SUBROUTINE' // achar(27) // '[0m.'
						end if

		        ! ----------------------------------------------------

		      end do

					! ----------------------------------------------------

		      exit convnloopxv

		    end if
		  end do convnloopxv
		end if

		! ----------------------------------------------------

  	!do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
  	!	if ((n /= 1) .and. (nn /= 1) .and. (n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1)))) then
  	! update grid with conservation of f and particle number,
  	! update LBNominalDensityGT and UBNominalDensityGT
  	! update all current particle positions and velocities with convection and with betatron acceleration
  	!do j= 1, (NsTK(1)- dNsTK2(1)- dNsTK3(1)), 1
  	! xN(j), yN(j), zN(j), Vperp1N(j), Vperp2N(j), VperpN(j), VxN(j), VyN(j), VzN(j) for ions and ENAs.
  	! Note: all new LB and UB injected particles are reset on correct Lshell in KineticSolverB given
  	! SpecieT(s)%FluxTubeT(f)%pGCGT(nn, Qind) and SpecieT(s)%FluxTubeT(f)%phiGCGT(nn, Qind)
  	!1- Conserve f and Update GRID
  	!2- Update LB and UB densities and ensure robust statistics
  	!3- Translate all particle positions in (p, q, phi)
  	!4- Update Vperp values with betatron ACCELERATION
  	!5- On new field line, update AEAmagN(j), AGmagN(j), AEPmagN(j)
  	!6- In KineticUpdateA, KineticUpdateB, KineticUpdateC add parallel centrifugal acceleration

		! ----------------------------------------------------

	end subroutine ConvectionSub

end module Convection
