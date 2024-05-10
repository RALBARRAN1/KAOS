module Convection

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.11- CONVECTION:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use ConfigGridGenerator
use VelGridGenerator

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

	subroutine ConvectionSub

  	! ----------------------------------------------------

  	! UPDATE TIME-DEPENDENT PHASE-SPACE GRID WITHOUT CONVECTING FLUX-TUBES:

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

  	! ----------------------------------------------------

  	! UPDATE PHASE-SPACE GRID WITH CONVECTING FLUX-TUBES:

  	if (SpecieT(s)%FluxTubeT(f)%CONVECTIONflagT(1) == 1) then
  		convnloop: do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
  			if ((n /= 1) .and. (nn /= 1) .and. &
  				(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn- 1))) then

  				! ----------------------------------------------------

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

  				!call VelGridGeneratorSub

  				! ----------------------------------------------------

  				exit convnloop

  			end if
  		end do convnloop
  	end if

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
