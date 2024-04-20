module ConfigGridGenerator

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	0 - CONFIGURATION-SPACE GRID GENERATOR:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use GridDipolePolynomialSolver

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! CONSTRUCT PHASE-SPACE GRID:

	subroutine ConfigGridGeneratorSub

! What we need for Grid Generator output:
! qGC0T, hqC0T, dpC0T, dqC0T, dphiC0T, rGC0T, phiGC0T, thetaGC0T, ellGC0T, qGL0T, qGH0T, pGC0T, d3xC0T, TsPerp0T, TsPar0T, Ts0T, &
! NqG0T, NqGT, Te0T

	  ! ----------------------------------------------------

    ! SET PRELIMINARY NUMBER OF CONFIGURATION-SPACE GRID CELLS PER PARTICLE SPECIES AND
    ! FLUX TUBE:

    if (qGA <= 0d0) then ! Southern Magnetic Hemisphere
      SMagHemFlag= 1
    end if
    if (qGA > 0d0) then ! Northern Magnetic Hemisphere
      SMagHemFlag= 0
    end if

    ! Set Vperp parameters
    if (ION2VPERPflag == 1) then
      NVperp1GpF= NVperpGpF
      NVperp2GpF= NVperpGpF
      Vperp12sigmaFac= VperpsigmaFac
      Vperp12NlinRange= (NVperp1GpF)/2d0 ! Number of grid cells spanning linear grid
    end if
    if (ION2VPERPflag == 0) then
      VperpNlinRange= (NVperpGpF)/2d0 ! Number of grid cells spanning linear grid
    end if

    ! Set Vpar parameters
    if (SYMVPARflag == 1) then
        VparNlinRange= (NVparGpF)/2d0
    end if
    if (SYMVPARflag == 0) then
        VparNlinRange= NVparGpF
    end if

    if (QEXCHANGEflag == 1) then
      if (ION2VPERPflag == 1) then
        NVpGpF= NVperp1GpF ! Number of Vp Grid Cells (even (div by 2 odd) for +/- log10 Vp grid) (+ 3)
        NVqGpF= NVparGpF ! Number of Vq Grid Cells (even (div by 2 odd) for +/- log10 Vq grid) (+ 3)
        NVphiGpF= NVperp2GpF ! Number of Vphi Grid Cells (even (div by 2 odd) for +/- log10 Vphi grid) (+ 3)
      end if

      if (ION2VPERPflag == 0) then
        write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NO ENA 1D VPERP GRID CONSTRUCTED', &
          ' IN CONFIG. GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
      end if

      if (ICRflag == 0) then
        VpphisigmaFac= Vperp12sigmaFac ! Sigma factor with linear grid to resolve thermal core of MB distrib.
        VqsigmaFac= VparsigmaFac

        VpphiNlinRange= (NVpGpF)/2d0
        if (SYMVPARflag == 1) then
            VqNlinRange= (NVqGpF)/2d0
        end if
        if (SYMVPARflag == 0) then
            VqNlinRange= NVqGpF
        end if
      end if
      if (ICRflag == 1) then
        VpphisigmaFac= Vperp12sigmaFac ! Sigma factor with linear grid to resolve thermal core of MB distrib.
        VqsigmaFac= VparsigmaFac

        VpphiNlinRange= (NVpGpF)/2d0
        if (SYMVPARflag == 1) then
            VqNlinRange= (NVqGpF)/2d0
        end if
        if (SYMVPARflag == 0) then
            VqNlinRange= NVqGpF
        end if
      end if
    end if

    ! Set Velocity space grid spacing by MB sigma value for initial temperature only. FIXME update velocity grid with new base altitude Ti, Te. in time
    TsPerp= Ti ! Initialize as isotropic (in pitch-angle), where Ti= (1/3)*TsPar+ (2/3)*TsPerp
    TsPar= Ti
    Vperp12sigma= sqrt(kB*TsPerp/mO) ! MB sigma for Vperp12
    Vperpsigma= Vperp12sigma
    Vparsigma= sqrt(kB*TsPar/mO) ! MB sigma for Vpar

    if (QEXCHANGEflag == 1) then
      mNeut= mO ! O neutral mass [kg]
      Vpphisigma= sqrt(kB*TNeut/mNeut) ! MB sigma for Vp, Vphi
      Vqsigma= sqrt(kB*TNeut/mNeut) ! MB sigma for Vq
    end if

    ! ----------------------------------------------------

    ! Total number of q grid values (Make even number to avoid equatorial
    ! grid boundary):
    if ((s == 1) .and. (f .lt. SpecieT(s)%NfT(1)/2d0)) then
      SpecieT(s)%FluxTubeT(f)%NqGpT(1)= NqGpF
    end if
    if ((s == 1) .and. (f .ge. SpecieT(s)%NfT(1)/2d0)) then
      SpecieT(s)%FluxTubeT(f)%NqGpT(1)= NqGpF
    end if
    if ((s == 2) .and. (f .lt. SpecieT(s)%NfT(1)/2d0)) then
      SpecieT(s)%FluxTubeT(f)%NqGpT(1)= NqGpF
    end if
    if ((s == 2) .and. (f .ge. SpecieT(s)%NfT(1)/2d0)) then
      SpecieT(s)%FluxTubeT(f)%NqGpT(1)= NqGpF
    end if

		! ---------------------------------------------
		! ENSURE CORRECT CONFIGURATION-SPACE GRID CELL DIMENSION:

		if ((NqUB(1)- NqLB(1)+ 4d0) /= SpecieT(s)%FluxTubeT(f)%NqGpT(1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' (NqUB(1)- NqLB(1)+ 4d0)= ', (NqUB(1)- NqLB(1)+ 3d0), &
				', NqGpT= ', SpecieT(s)%FluxTubeT(f)%NqGpT(1), ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
		end if

		! ---------------------------------------------

    ! SET MAGNETIC L-SHELL VALUES AND PRELIMINARY FIELD-ALIGNED GRID CELLS FOR EACH PARTICLE SPECIES
    ! AND FLUX-TUBE:

    SpecieT(s)%FluxTubeT(f)%NqGT(1)= (SpecieT(s)%FluxTubeT(f)%NqGpT(1)/dq)- 1d0 ! Number of config cells

    allocate(SpecieT(s)%FluxTubeT(f)%qGT(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))
		allocate(SpecieT(s)%FluxTubeT(f)%pGT(1))
    allocate(SpecieT(s)%FluxTubeT(f)%phiGridOutT(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%rGridOutT(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%thetaGridOutT(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%xGridOutT(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%yGridOutT(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%zGridOutT(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))

		SpecieT(s)%FluxTubeT(f)%pGT(1)= Lshell ! Set grid L-shell

    do Qind= 1, SpecieT(s)%FluxTubeT(f)%NqGpT(1), 1
      if (Qind /= SpecieT(s)%FluxTubeT(f)%NqGpT(1)) then
				if (qGA <= 0d0) then ! SMH
	        SpecieT(s)%FluxTubeT(f)%qGT(Qind)= qGA+ (Qind- 1d0)*(abs(qGB- qGA)/SpecieT(s)%FluxTubeT(f)%NqGpT(1))
				end if
				if (qGA > 0d0) then ! NMH
	        SpecieT(s)%FluxTubeT(f)%qGT(Qind)= qGA- (Qind- 1d0)*(abs(qGB- qGA)/SpecieT(s)%FluxTubeT(f)%NqGpT(1))
				end if
      end if
      if (Qind == SpecieT(s)%FluxTubeT(f)%NqGpT(1)) then
        SpecieT(s)%FluxTubeT(f)%qGT(Qind)= qGB
      end if
    end do

    do Qind= 1, SpecieT(s)%FluxTubeT(f)%NqGpT(1), 1

      pGridIn(1)= SpecieT(s)%FluxTubeT(f)%pGT(1)
      qGridIn(1)= SpecieT(s)%FluxTubeT(f)%qGT(Qind)

      call GridDipolePolynomialSolverSub

      SpecieT(s)%FluxTubeT(f)%rGridOutT(Qind)= rGridOut(1)
      SpecieT(s)%FluxTubeT(f)%thetaGridOutT(Qind)= thetaGridOut(1)
			SpecieT(s)%FluxTubeT(f)%phiGridOutT(Qind)= phiGridOut(1)
      SpecieT(s)%FluxTubeT(f)%xGridOutT(Qind)= xGridOut(1)
      SpecieT(s)%FluxTubeT(f)%yGridOutT(Qind)= yGridOut(1)
      SpecieT(s)%FluxTubeT(f)%zGridOutT(Qind)= zGridOut(1)

    end do

    ! ---------------------------------------------

    ! SET CONFIGURATION SPACE GRID LIMITS:

    allocate(SpecieT(s)%FluxTubeT(f)%phiGL0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%phiGH0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%rGL0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%rGH0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%thetaGL0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%thetaGH0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%xGC0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%xGL0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%xGH0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%yGC0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%yGL0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%yGH0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%zGC0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%zGL0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%zGH0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%ellGL0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%ellGH0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%hpC0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%hphiC0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))

    do Qind= NqLB(1), NqUB(1)+ 2, 1

      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGL0T(1)= &
        SpecieT(s)%FluxTubeT(f)%qGT((Qind- 1)*dq+ 1) ! lower q limits of FA cells
      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGH0T(1)= &
        SpecieT(s)%FluxTubeT(f)%qGT((Qind- 1)*dq+ 1+ dq) ! upper q limits of FA cells
      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGC0T(1)= &
        (SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGL0T(1)+ &
	      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGH0T(1))/2d0 ! center q values of FA cells

      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(1)= SpecieT(s)%FluxTubeT(f)%pGT(1) ! center p values of FA cells

      SpecieT(s)%FluxTubeT(f)%phiGL0T(Qind)= &
        SpecieT(s)%FluxTubeT(f)%phiGridOutT((Qind- 1)*dq+ 1) ! lower phid limits of FA cells
      SpecieT(s)%FluxTubeT(f)%phiGH0T(Qind)= &
        SpecieT(s)%FluxTubeT(f)%phiGridOutT((Qind- 1)*dq+ 1+ dq) ! upper phid limits of FA cells

      SpecieT(s)%FluxTubeT(f)%rGL0T(Qind)= &
        SpecieT(s)%FluxTubeT(f)%rGridOutT((Qind- 1)*dq+ 1) ! lower r limits of FA cells
      SpecieT(s)%FluxTubeT(f)%rGH0T(Qind)= &
        SpecieT(s)%FluxTubeT(f)%rGridOutT((Qind- 1)*dq+ 1+ dq) ! upper r limits of FA cells

      SpecieT(s)%FluxTubeT(f)%thetaGL0T(Qind)= &
        SpecieT(s)%FluxTubeT(f)%thetaGridOutT((Qind- 1)*dq+ 1) ! lower theta limits of FA cells
      SpecieT(s)%FluxTubeT(f)%thetaGH0T(Qind)= &
        SpecieT(s)%FluxTubeT(f)%thetaGridOutT((Qind- 1)*dq+ 1+ dq) ! upper theta limits of FA cells

      SpecieT(s)%FluxTubeT(f)%yGL0T(Qind)= &
        SpecieT(s)%FluxTubeT(f)%yGridOutT((Qind- 1)*dq+ 1) ! lower y limits of FA cells
      SpecieT(s)%FluxTubeT(f)%yGH0T(Qind)= &
        SpecieT(s)%FluxTubeT(f)%yGridOutT((Qind- 1)*dq+ 1+ dq) ! upper y limits of FA cells

      SpecieT(s)%FluxTubeT(f)%zGL0T(Qind)= &
        SpecieT(s)%FluxTubeT(f)%zGridOutT((Qind- 1)* dq+ 1) ! lower z limits of FA cells
      SpecieT(s)%FluxTubeT(f)%zGH0T(Qind)= &
        SpecieT(s)%FluxTubeT(f)%zGridOutT((Qind- 1)*dq+ 1+ dq) ! upper z limits of FA cells

      SpecieT(s)%FluxTubeT(f)%xGL0T(Qind)= 0d0 ! Let phi= pi/2 s.t. all FA cells have x= 0
      SpecieT(s)%FluxTubeT(f)%xGH0T(Qind)= 0d0

    end do

    do Qind= NqLB(1), NqUB(1)+ 2, 1

      pGridIn(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(1)
      qGridIn(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGC0T(1)

      call GridDipolePolynomialSolverSub

      ! SET CONFIGURATION SPACE GRID CENTER VALUES, METRIC FACTORS, AND VOLUMES:

      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(1)= rGridOut(1)
      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(1)= thetaGridOut(1)
      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(1)= phiGridOut(1)
      SpecieT(s)%FluxTubeT(f)%xGC0T(Qind)= xGridOut(1)
      SpecieT(s)%FluxTubeT(f)%yGC0T(Qind)= yGridOut(1)
      SpecieT(s)%FluxTubeT(f)%zGC0T(Qind)= zGridOut(1)

      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%ellGC0T(1)= &
        (1d0+ 3d0*cos(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(1))**2d0)
      SpecieT(s)%FluxTubeT(f)%ellGL0T(Qind)= &
        (1d0+ 3d0*cos(SpecieT(s)%FluxTubeT(f)%thetaGL0T(Qind))**2d0)
      SpecieT(s)%FluxTubeT(f)%ellGH0T(Qind)= &
        (1d0+ 3d0*cos(SpecieT(s)%FluxTubeT(f)%thetaGH0T(Qind))**2d0)

      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%hqC0T(1)= &
        abs(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(1)**3d0/ &
        (RE**2d0*sqrt(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%ellGC0T(1)))) ! q metric factor centered in FA cells
      SpecieT(s)%FluxTubeT(f)%hpC0T(Qind)= &
        abs(RE*sin(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(1))**3d0/ &
        (sqrt(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%ellGC0T(1)))) ! p metric factor centered in FA cells
      SpecieT(s)%FluxTubeT(f)%hphiC0T(Qind)= &
        abs(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(1)* &
        sin(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(1))) ! phi metric factor centered in FA cells

      ! Note: Make dpC equal to the Lshell drift limit over entire simulation and dphidC equal to the
      ! phid drift limit in kinetic solver f90 subroutine.

      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dqC0T(1)= &
        abs(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGH0T(1)- SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGL0T(1)) ! dq across FA cells
      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dpC0T(1)= 1d-3 ! dp across FA cells [RE]
      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dphiC0T(1)= pi/180d0 ! dphi across FA cells [1 deg= pi/180 rads]

      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1)= &
        SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%hqC0T(1)* &
        SpecieT(s)%FluxTubeT(f)%hpC0T(Qind)* &
        SpecieT(s)%FluxTubeT(f)%hphiC0T(Qind)* &
        SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dqC0T(1)* &
        SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dpC0T(1)* &
        SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dphiC0T(1)

      if (SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1) == 0d0) then
        SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1)= 1d-9
      end if

      if (SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1) < 0d0) then
        write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
          ' NEGATIVE d3xC0T VALUE IN GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
      end if

    end do

    SpecieT(s)%Qindns0T(1)= 1d0

    ! ---------------------------------------------

    ! Re-index to account for lower and upper ghost cells (non-computational domain)
    SpecieT(s)%FluxTubeT(f)%NqG0T(1)= SpecieT(s)%FluxTubeT(f)%NqGT(1)
    SpecieT(s)%FluxTubeT(f)%NqGT(1)= SpecieT(s)%FluxTubeT(f)%NqGT(1)- 2d0

		allocate(SpecieT(s)%FluxTubeT(f)%Te0T(SpecieT(s)%FluxTubeT(f)%NqG0T(1)))

    do Qind= NqLB(1), NqUB(1)+ 2, 1
      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPerp0T(1)= TsPerp
      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(1)= TsPar
      SpecieT(s)%FluxTubeT(f)%Te0T(Qind)= Te
    end do

    if (SpecieT(s)%FluxTubeT(f)%NqG0T(1) /= (NqUB(1)- NqLB(1))+ 3) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' BAD NqG0T VALUE', &
        ' IN GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
    end if

    ! ---------------------------------------------

  end subroutine ConfigGridGeneratorSub

end module ConfigGridGenerator
