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

		! SET GRID RANGE ALONG FLUX TUBES:

		! ----------------------------------------------------

		if ((INITIALGRIDflag == 1) .and. (nn /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
				' INCONSISTENT INITIAL GRID FLAG FOR SPECIE= ', s, &
				' , FLUX TUBE= ', f, ', AND STATISTICAL TIME-STEP= ', nn, &
				' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((INITIALGRIDflag == 0) .and. (nn == 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
				' INCONSISTENT INITIAL GRID FLAG FOR SPECIE= ', s, &
				' , FLUX TUBE= ', f, ', AND STATISTICAL TIME-STEP= ', nn, &
				' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		NqGpF= (SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 4)
		NqIC(1)= abs(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1d0 ! Q cell range for density initialization
		SpecieT(s)%FluxTubeT(f)%NqICT= NqIC

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

		if ((size(SpecieT(s)%FluxTubeT(f)%NqLBT(:)) /= 1) .or. &
			(isnan(real(SpecieT(s)%FluxTubeT(f)%NqLBT(1))) .eqv. .true.)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqLBT HAS', &
				' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
				' AND FLUX TUBE= ', f, ' IN CONFIGURATION-SPACE GRID GENERATOR', &
				' SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((size(SpecieT(s)%FluxTubeT(f)%NqUBT(:)) /= 1) .or. &
			(isnan(real(SpecieT(s)%FluxTubeT(f)%NqUBT(1))) .eqv. .true.)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqUBT HAS', &
				' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
				' AND FLUX TUBE= ', f, ' IN CONFIGURATION-SPACE GRID GENERATOR', &
				' SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((NqIC(1) /= SpecieT(s)%FluxTubeT(f)%NqICT(1)) .or. &
			(size(SpecieT(s)%FluxTubeT(f)%NqICT(:)) /= 1) .or. &
			(isnan(real(SpecieT(s)%FluxTubeT(f)%NqICT(1))) .eqv. .true.)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqICT HAS', &
				' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
				' AND FLUX TUBE= ', f, ' IN CONFIGURATION-SPACE GRID GENERATOR', &
				' SUBROUTINE' // achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! ALLOCATE INITIALIZATION CONFIG-SPACE DERIVED DATA TYPES NESTED IN FLUX TUBE
		! TYPES NESTED IN PARTICLE SPECIES TYPE:

		! Allocate QCellICT(Qind) derived data type nested in FluxTubeT(f)
		! nested in SpecieT(s)
		if (INITIALGRIDflag == 1) then
			allocate(SpecieT(s)%FluxTubeT(f)%QCellICT(SpecieT(s)%FluxTubeT(f)%NqICT(1)))
		end if

	  ! ----------------------------------------------------

    ! SET PRELIMINARY NUMBER OF CONFIGURATION-SPACE GRID CELLS PER PARTICLE SPECIES AND
    ! FLUX TUBE:

    if (qGA(1) <= 0d0) then ! Southern Magnetic Hemisphere
      SMagHemFlag= 1
    end if
    if (qGA(1) > 0d0) then ! Northern Magnetic Hemisphere
      SMagHemFlag= 0
    end if

    ! Set Vperp parameters
    if (ION2VPERPflag == 1) then
      NVperp1GpF= NVperpGpF
      NVperp2GpF= NVperpGpF
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
        VpphisigmaFac= VperpsigmaFac(1) ! Sigma factor with linear grid to resolve thermal core of MB distrib.
        VqsigmaFac= VparsigmaFac(1)

        VpphiNlinRange= (NVpGpF)/2d0
        if (SYMVPARflag == 1) then
            VqNlinRange= (NVqGpF)/2d0
        end if
        if (SYMVPARflag == 0) then
            VqNlinRange= NVqGpF
        end if
      end if
      if (ICRflag == 1) then
        VpphisigmaFac= VperpsigmaFac(1) ! Sigma factor with linear grid to resolve thermal core of MB distrib.
        VqsigmaFac= VparsigmaFac(1)

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
		TsPerp= Ti(1) ! Initialize as isotropic (in pitch-angle), where Ti= (1/3)*TsPar+ (2/3)*TsPerp
    TsPar= Ti(1)
    Vperp12sigma= sqrt(kB*TsPerp/SpecieT(s)%msT(1)) ! MB sigma for Vperp12
    Vperpsigma= Vperp12sigma
    Vparsigma= sqrt(kB*TsPar/SpecieT(s)%msT(1)) ! MB sigma for Vpar

    if (QEXCHANGEflag == 1) then
      mNeut= SpecieT(s)%msT(1) ! O neutral mass [kg]
      Vpphisigma= sqrt(kB*TNeut(1)/mNeut) ! MB sigma for Vp, Vphi
      Vqsigma= sqrt(kB*TNeut(1)/mNeut) ! MB sigma for Vq
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

		if ((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- &
			SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 4d0) /= SpecieT(s)%FluxTubeT(f)%NqGpT(1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
				' (SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 4d0)= ', &
				(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3d0), &
				', NqGpT= ', SpecieT(s)%FluxTubeT(f)%NqGpT(1), ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
		end if

		! ---------------------------------------------

    ! SET MAGNETIC L-SHELL VALUES AND PRELIMINARY FIELD-ALIGNED GRID CELLS FOR EACH PARTICLE SPECIES
    ! AND FLUX-TUBE:

		SpecieT(s)%FluxTubeT(f)%NqGTp(1)= (SpecieT(s)%FluxTubeT(f)%NqGpT(1)/dq)- 1d0 ! Number of config cells

		if (INITIALGRIDflag == 1) then
	    allocate(qGp(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))
			allocate(pGp(1))
			allocate(phiGp(1))
	    allocate(rGridOutp(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))
	    allocate(thetaGridOutp(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))
	    allocate(xGridOutp(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))
	    allocate(yGridOutp(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))
	    allocate(zGridOutp(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))
		end if

		pGp(1)= Lshell(1) ! Set grid L-shell
		phiGp(1)= phiLshell(1) ! Set grid longitude

    do Qind= 1, SpecieT(s)%FluxTubeT(f)%NqGpT(1), 1
      if (Qind /= SpecieT(s)%FluxTubeT(f)%NqGpT(1)) then
				if (qGA(1) <= 0d0) then ! SMH
	        qGp(Qind)= qGA(1)+ (Qind- 1d0)*(abs(qGB(1)- qGA(1))/SpecieT(s)%FluxTubeT(f)%NqGpT(1))
				end if
				if (qGA(1) > 0d0) then ! NMH
	        qGp(Qind)= qGA(1)- (Qind- 1d0)*(abs(qGB(1)- qGA(1))/SpecieT(s)%FluxTubeT(f)%NqGpT(1))
				end if
      end if
      if (Qind == SpecieT(s)%FluxTubeT(f)%NqGpT(1)) then
        qGp(Qind)= qGB(1)
      end if
    end do

    do Qind= 1, SpecieT(s)%FluxTubeT(f)%NqGpT(1), 1

      pGridIn(1)= pGp(1)
      qGridIn(1)= qGp(Qind)
			phiGridIn(1)= phiGp(1)

      call GridDipolePolynomialSolverSub

      rGridOutp(Qind)= rGridOut(1)
      thetaGridOutp(Qind)= thetaGridOut(1)
      xGridOutp(Qind)= xGridOut(1)
      yGridOutp(Qind)= yGridOut(1)
      zGridOutp(Qind)= zGridOut(1)

			! ---------------------------------------------

			! DIAGNOSTIC FLAGS FOR NAN VALUES:

			if ((isnan(real(rGridOutp(Qind))) .eqv. .true.) .or. &
				(size(rGridOutp(:)) /= SpecieT(s)%FluxTubeT(f)%NqGpT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rGridOutp= ', &
					rGridOutp(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(thetaGridOutp(Qind))) .eqv. .true.) .or. &
				(size(thetaGridOutp(:)) /= SpecieT(s)%FluxTubeT(f)%NqGpT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' thetaGridOutp= ', &
					thetaGridOutp(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(xGridOutp(Qind))) .eqv. .true.) .or. &
				(size(xGridOutp(:)) /= SpecieT(s)%FluxTubeT(f)%NqGpT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' xGridOutp= ', &
					xGridOutp(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(yGridOutp(Qind))) .eqv. .true.) .or. &
				(size(yGridOutp(:)) /= SpecieT(s)%FluxTubeT(f)%NqGpT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' yGridOutp= ', &
					yGridOutp(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(zGridOutp(Qind))) .eqv. .true.) .or. &
				(size(zGridOutp(:)) /= SpecieT(s)%FluxTubeT(f)%NqGpT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' zGridOutp= ', &
					zGridOutp(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			! ---------------------------------------------

    end do

    ! ---------------------------------------------

    ! SET CONFIGURATION SPACE GRID LIMITS:

		if (INITIALGRIDflag == 1) then
	    allocate(SpecieT(s)%FluxTubeT(f)%rGL0T(SpecieT(s)%FluxTubeT(f)%NqGTp(1)))
	    allocate(SpecieT(s)%FluxTubeT(f)%rGH0T(SpecieT(s)%FluxTubeT(f)%NqGTp(1)))
	    allocate(SpecieT(s)%FluxTubeT(f)%thetaGL0T(SpecieT(s)%FluxTubeT(f)%NqGTp(1)))
	    allocate(SpecieT(s)%FluxTubeT(f)%thetaGH0T(SpecieT(s)%FluxTubeT(f)%NqGTp(1)))
	    allocate(SpecieT(s)%FluxTubeT(f)%xGC0T(SpecieT(s)%FluxTubeT(f)%NqGTp(1)))
	    allocate(SpecieT(s)%FluxTubeT(f)%xGL0T(SpecieT(s)%FluxTubeT(f)%NqGTp(1)))
	    allocate(SpecieT(s)%FluxTubeT(f)%xGH0T(SpecieT(s)%FluxTubeT(f)%NqGTp(1)))
	    allocate(SpecieT(s)%FluxTubeT(f)%yGC0T(SpecieT(s)%FluxTubeT(f)%NqGTp(1)))
	    allocate(SpecieT(s)%FluxTubeT(f)%yGL0T(SpecieT(s)%FluxTubeT(f)%NqGTp(1)))
	    allocate(SpecieT(s)%FluxTubeT(f)%yGH0T(SpecieT(s)%FluxTubeT(f)%NqGTp(1)))
	    allocate(SpecieT(s)%FluxTubeT(f)%zGC0T(SpecieT(s)%FluxTubeT(f)%NqGTp(1)))
	    allocate(SpecieT(s)%FluxTubeT(f)%zGL0T(SpecieT(s)%FluxTubeT(f)%NqGTp(1)))
	    allocate(SpecieT(s)%FluxTubeT(f)%zGH0T(SpecieT(s)%FluxTubeT(f)%NqGTp(1)))
	    allocate(SpecieT(s)%FluxTubeT(f)%ellGL0T(SpecieT(s)%FluxTubeT(f)%NqGTp(1)))
	    allocate(SpecieT(s)%FluxTubeT(f)%ellGH0T(SpecieT(s)%FluxTubeT(f)%NqGTp(1)))
	    allocate(SpecieT(s)%FluxTubeT(f)%hpC0T(SpecieT(s)%FluxTubeT(f)%NqGTp(1)))
	    allocate(SpecieT(s)%FluxTubeT(f)%hphiC0T(SpecieT(s)%FluxTubeT(f)%NqGTp(1)))
		end if

    do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 2, 1

      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGL0T(1)= &
        qGp((Qind- 1)*dq+ 1) ! lower q limits of FA cells
      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGH0T(1)= &
        qGp((Qind- 1)*dq+ 1+ dq) ! upper q limits of FA cells
      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGC0T(1)= &
        (SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGL0T(1)+ &
	      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGH0T(1))/2d0 ! center q values of FA cells

      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(1)= pGp(1) ! center p values of FA cells
			SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(1)= phiGp(1) ! center phi values of FA cells

      SpecieT(s)%FluxTubeT(f)%rGL0T(Qind)= &
        rGridOutp((Qind- 1)*dq+ 1) ! lower r limits of FA cells
      SpecieT(s)%FluxTubeT(f)%rGH0T(Qind)= &
        rGridOutp((Qind- 1)*dq+ 1+ dq) ! upper r limits of FA cells

      SpecieT(s)%FluxTubeT(f)%thetaGL0T(Qind)= &
        thetaGridOutp((Qind- 1)*dq+ 1) ! lower theta limits of FA cells
      SpecieT(s)%FluxTubeT(f)%thetaGH0T(Qind)= &
        thetaGridOutp((Qind- 1)*dq+ 1+ dq) ! upper theta limits of FA cells

      SpecieT(s)%FluxTubeT(f)%yGL0T(Qind)= &
        yGridOutp((Qind- 1)*dq+ 1) ! lower y limits of FA cells
      SpecieT(s)%FluxTubeT(f)%yGH0T(Qind)= &
        yGridOutp((Qind- 1)*dq+ 1+ dq) ! upper y limits of FA cells

      SpecieT(s)%FluxTubeT(f)%zGL0T(Qind)= &
        zGridOutp((Qind- 1)* dq+ 1) ! lower z limits of FA cells
      SpecieT(s)%FluxTubeT(f)%zGH0T(Qind)= &
        zGridOutp((Qind- 1)*dq+ 1+ dq) ! upper z limits of FA cells

			SpecieT(s)%FluxTubeT(f)%xGL0T(Qind)= &
        xGridOutp((Qind- 1)* dq+ 1) ! lower x limits of FA cells
      SpecieT(s)%FluxTubeT(f)%xGH0T(Qind)= &
        xGridOutp((Qind- 1)*dq+ 1+ dq) ! upper x limits of FA cells

			! ---------------------------------------------

			! DIAGNOSTIC FLAGS FOR NAN VALUES:

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGL0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGL0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qGL0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGL0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGH0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGH0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qGH0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGH0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGC0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGC0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qGC0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGC0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' pGC0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%rGL0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%rGL0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGTp(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rGL0T= ', &
					SpecieT(s)%FluxTubeT(f)%rGL0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%rGH0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%rGH0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGTp(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rGH0T= ', &
					SpecieT(s)%FluxTubeT(f)%rGH0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%thetaGL0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%thetaGL0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGTp(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' thetaGL0T= ', &
					SpecieT(s)%FluxTubeT(f)%thetaGL0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%thetaGH0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%thetaGH0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGTp(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' thetaGH0T= ', &
					SpecieT(s)%FluxTubeT(f)%thetaGH0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%yGL0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%yGL0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGTp(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' yGL0T= ', &
					SpecieT(s)%FluxTubeT(f)%yGL0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%yGH0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%yGH0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGTp(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' yGH0T= ', &
					SpecieT(s)%FluxTubeT(f)%yGH0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%zGL0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%zGL0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGTp(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' zGL0T= ', &
					SpecieT(s)%FluxTubeT(f)%zGL0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%zGH0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%zGH0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGTp(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' zGH0T= ', &
					SpecieT(s)%FluxTubeT(f)%zGH0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%xGL0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%xGL0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGTp(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' xGL0T= ', &
					SpecieT(s)%FluxTubeT(f)%xGL0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%xGH0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%xGH0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGTp(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' xGH0T= ', &
					SpecieT(s)%FluxTubeT(f)%xGH0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			! ---------------------------------------------

    end do

    do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 2, 1

      pGridIn(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(1)
      qGridIn(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGC0T(1)
			phiGridIn(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(1)

      call GridDipolePolynomialSolverSub

      ! SET CONFIGURATION SPACE GRID CENTER VALUES, METRIC FACTORS, AND VOLUMES:

      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(1)= rGridOut(1)
      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(1)= thetaGridOut(1)
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

			! ---------------------------------------------

			! DIAGNOSTIC FLAGS FOR NAN VALUES:

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rGC0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' thetaGC0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' phiGC0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if (SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(1) /= &
				SpecieT(s)%FluxTubeT(f)%QCell0T(1)%phiGC0T(1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' phiGC0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(1), &
					' DOES NOT EQUAL FLUX-TUBE FOOT PRINT phiCG0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(1)%phiGC0T(1), ' FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if (SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(1) /= &
				SpecieT(s)%FluxTubeT(f)%QCell0T(1)%pGC0T(1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' pGC0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(1), &
					' DOES NOT EQUAL FLUX-TUBE FOOT PRINT phiCG0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(1)%pGC0T(1), ' FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%xGC0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%xGC0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGTp(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' xGC0T= ', &
					SpecieT(s)%FluxTubeT(f)%xGC0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%yGC0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%yGC0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGTp(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' yGC0T= ', &
					SpecieT(s)%FluxTubeT(f)%yGC0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%zGC0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%zGC0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGTp(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' zGC0T= ', &
					SpecieT(s)%FluxTubeT(f)%zGC0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%ellGC0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%ellGC0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ellGC0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%ellGC0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%ellGL0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%ellGL0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGTp(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ellGL0T= ', &
					SpecieT(s)%FluxTubeT(f)%ellGL0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%ellGH0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%ellGH0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGTp(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ellGH0T= ', &
					SpecieT(s)%FluxTubeT(f)%ellGH0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%hqC0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%hqC0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' hqC0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%hqC0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%hpC0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%hpC0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGTp(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' hpC0T= ', &
					SpecieT(s)%FluxTubeT(f)%hpC0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%hphiC0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%hphiC0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGTp(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' hphiC0T= ', &
					SpecieT(s)%FluxTubeT(f)%hphiC0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dqC0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dqC0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dqC0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dqC0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dpC0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dpC0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dpC0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dpC0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dphiC0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dphiC0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dphiC0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dphiC0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' d3xC0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			! ---------------------------------------------

    end do

    SpecieT(s)%Qindns0GT(1)= 1d0

    ! ---------------------------------------------

    ! Re-index to account for lower and upper ghost cells (non-computational domain)
    SpecieT(s)%FluxTubeT(f)%NqG0T(1)= SpecieT(s)%FluxTubeT(f)%NqGTp(1)
    SpecieT(s)%FluxTubeT(f)%NqGT(1)= SpecieT(s)%FluxTubeT(f)%NqGTp(1)- 2d0

		if (INITIALGRIDflag == 1) then
			allocate(SpecieT(s)%FluxTubeT(f)%Te0T(SpecieT(s)%FluxTubeT(f)%NqG0T(1)))
		end if

    do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 2, 1
      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPerp0T(1)= TsPerp
      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(1)= TsPar
      SpecieT(s)%FluxTubeT(f)%Te0T(Qind)= Te(1)
    end do

    if (SpecieT(s)%FluxTubeT(f)%NqG0T(1) /= (SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 3) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' BAD INITIAL NqG0T VALUE', &
        ' IN GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
    end if

		! ---------------------------------------------

		! DIAGNOSTIC FLAGS FOR NAN VALUES:

		if ((isnan(real(SpecieT(s)%Qindns0GT(1))) .eqv. .true.) .or. &
			(size(SpecieT(s)%Qindns0GT(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Qindns0GT= ', &
				SpecieT(s)%Qindns0GT(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
				' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(SpecieT(s)%FluxTubeT(f)%NqG0T(1))) .eqv. .true.) .or. &
			(size(SpecieT(s)%FluxTubeT(f)%NqG0T(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqG0T= ', &
				SpecieT(s)%FluxTubeT(f)%NqG0T(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
				', FLUX TUBE= ', f, ' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(SpecieT(s)%FluxTubeT(f)%NqGT(1))) .eqv. .true.) .or. &
			(size(SpecieT(s)%FluxTubeT(f)%NqGT(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqGT= ', &
				SpecieT(s)%FluxTubeT(f)%NqGT(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
				', FLUX TUBE= ', f, ' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 2, 1

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPerp0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPerp0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' TsPerp0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPerp0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' TsPar0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%Te0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%Te0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Te0T= ', &
					SpecieT(s)%FluxTubeT(f)%Te0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

		end do

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

		if ((size(SpecieT(s)%FluxTubeT(f)%NqG0T(:)) /= 1) .or. &
			(isnan(real(SpecieT(s)%FluxTubeT(f)%NqG0T(1))) .eqv. .true.)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqG0T HAS', &
				' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
				' AND FLUX TUBE= ', f, &
				' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((size(SpecieT(s)%FluxTubeT(f)%NqGT(:)) /= 1) .or. &
			(isnan(real(SpecieT(s)%FluxTubeT(f)%NqGT(1))) .eqv. .true.)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqG0T HAS', &
				' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
				' AND FLUX TUBE= ', f, &
				' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
		end if

		do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 2, 1
			if ((size(SpecieT(s)%FluxTubeT(f)%Te0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1)) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%Te0T(Qind))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Te0T HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if (SpecieT(s)%FluxTubeT(f)%Te0T(Qind) > dNTeEND) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' INITIAL Te0T IS GREATER THAN Te CAP FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if
		end do

		! ----------------------------------------------------

		do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 2, 1

			! ----------------------------------------------------

			! Total initial ion temperature [K]
			! Note: Isotropic temperature from grid data (Ti= Tperp= Tpar= Tperp1= Tperp2)
			SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%Ts0T(1)= (1d0/3d0)*(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(1))+ &
				(2d0/3d0)*(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPerp0T(1))

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGC0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGC0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qGCGT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%hqC0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%hqC0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' hqCGT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dpC0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dpC0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dpCGT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dqC0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dqC0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dqCGT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dphiC0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dphiC0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dphiCGT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPerp0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPerp0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' TsPerpGT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' TsParGT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rGCGT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' phiGCGT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' thetaGCGT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%ellGC0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%ellGC0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ellGCGT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGL0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGL0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qGLGT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGH0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGH0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qGHGT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' pGCGT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' d3xCGT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if (SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1) < 0) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NEGATIVE ', &
					' d3xCGT FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			! ----------------------------------------------------

		end do

		! ----------------------------------------------------

		! SET STATISTICAL TIME-STEP WITH NEW GRID:

		! Injection time-step [s] (must be > 2 and excludes initial time-step)
		! Mean transit time through LB ghost cell

		Q0ndatfac(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(SpecieT(s)%FluxTubeT(f)%NqLBT(1))%hqC0T(1)* &
			SpecieT(s)%FluxTubeT(f)%QCell0T(SpecieT(s)%FluxTubeT(f)%NqLBT(1))%dqC0T(1)* &
			(abs(sqrt(SpecieT(s)%msT(1)/ &
			((kB/1d0)*(SpecieT(s)%FluxTubeT(f)%QCell0T(SpecieT(s)%FluxTubeT(f)%NqLBT(1))%TsPar0T(1))))))/SpecieT(s)%FluxTubeT(f)%hT(1)

		if (INITIALGRIDflag == 1) then
			SpecieT(s)%FluxTubeT(f)%NNtT(1)= SpecieT(s)%FluxTubeT(f)%NtT(1)/Q0ndatfac(1) ! Number of time-steps for injection
			allocate(SpecieT(s)%FluxTubeT(f)%ndatfacGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
		end if

		if (INITIALGRIDflag == 1) then
			! Time-step interval for moment computation (must be > 2 and excludes initial time-step)
			SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)= SpecieT(s)%FluxTubeT(f)%NtT(1)/SpecieT(s)%FluxTubeT(f)%NNtT(1)
		end if

		if (INITIALGRIDflag == 0) then

			NNtp(1)= SpecieT(s)%FluxTubeT(f)%NtT(1)/Q0ndatfac(1)

			! Time-step interval for moment computation (must be > 2 and excludes initial time-step)
			SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)= SpecieT(s)%FluxTubeT(f)%NtT(1)/NNtp(1)

		end if

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR STATISTICAL TIME-STEP SIZE:

		if (SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn) < 3d0) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ndatfacGT= ', &
				SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn), &
				' IS NOT COMPATIBLE WITH KINETIC SOLVER TIME FORMAT FOR SPECIE= ', s, &
				' , FLUX TUBE= ', f, ', AND STATISTICAL TIME-STEP= ', nn, &
				' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! COMPUTE ION DENSITY PROFILE:

		if (rank == 0) then

			if (INITIALGRIDflag == 1) then
				allocate(ICbbp(((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 3)))
			end if

			do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 2, 1

				! ----------------------------------------------------

				! COMPUTE STEADY-STATE ION DENSITY PROFILE WITH NORMALIZATION:
				! NOTE: SET ninormfac SUCH THAT NO FA CELL IS EMPTY.

				ICbbp(Qind)= 2d0*GG*ME*((SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%hqC0T(1)* &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dqC0T(1)* &
					cos(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(1)))/ &
					((SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(1)**2d0)* &
					(sqrt(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%ellGC0T(1)))))

			end do

			IC0bb(1)= sum(ICbbp(1:SpecieT(s)%Qindns0GT(1)))

			do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 2, 1

				!Altitude dependent gravitational acceleration, integration over ds
				ICbb(1)= sum(ICbbp(1:Qind))
        gC(1)= ICbb(1)- IC0bb(1)

				if (qGA(1) <= 0d0) then ! SMH
          argC(1)= ((SpecieT(s)%msT(1)*gC(1))/ &
           	(kB*(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(1)+ SpecieT(s)%FluxTubeT(f)%Te0T(Qind))))
				end if

				if (qGA(1) > 0d0) then ! NMH
          argC(1)= -((SpecieT(s)%msT(1)*gC(1))/ &
           	(kB*(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(1)+ SpecieT(s)%FluxTubeT(f)%Te0T(Qind))))
				end if

				nsC(1)= ns0(1)*exp(argC(1)) ! Number of ions per flux-tube grid cell [m^-3]

				! Normalized number of macroparticles per flux-tube grid cell [unitless]
				nsnormC(1)= (nsC(1)/SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))* &
					(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1))

				if (INITIALGRIDflag == 1) then
					if (SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1) == 1) then
						SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%nsnormCT= nsnormC
					end if
					if (SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1) == 0) then
						if (Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) then
							SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%nsnormCT(1)= SpecieT(s)%FluxTubeT(f)%nsnormCLBInputT(1)
						end if
						if (Qind == SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 2) then
							SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%nsnormCT(1)= SpecieT(s)%FluxTubeT(f)%nsnormCUBInputT(1)
						end if
						if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. (Qind /= SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 2)) then
							SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%nsnormCT(1)= &
								(SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%DensityInputT(1)/SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))* &
								(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1))
						end if
					end if
				end if

				if (INITIALGRIDflag == 0) then
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%nsnormCT= nsnormC
				end if

				! ----------------------------------------------------

				! COMPUTE STEADY-STATE O NEUTRAL DENSITY PROFILE:

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then

					gCNeut(1)= (GG*ME)/(kB*TNeut(1))

					if (qGA(1) <= 0d0) then ! SMH
						argCNeut(1)= mNeut*gCNeut(1)*((1d0/zns0Neut(1))- &
							(1d0/SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(1)))
					end if
					if (qGA(1) > 0d0) then ! NMH
						argCNeut(1)= -mNeut*gCNeut(1)*((1d0/zns0Neut(1))- &
							(1d0/SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(1)))
					end if

					nsCNeut(1)= ns0Neut(1)*exp(argCNeut(1)) ! Number of H neutrals per flux-tube grid cell [m^-3]

					! Normalized number of macroparticles per flux-tube grid cell [unitless]
					nsnormCNeut(1)= (nsCNeut(1)/SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))* &
						(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1))

					! Create nested derived data types
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%nsnormCNeut0T= nsnormCNeut

				end if

				! ----------------------------------------------------

				! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

				if ((size(gC(:)) /= 1) .or. (isnan(real(gC(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' gC HAS', &
						' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
						', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
						' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((size(argC(:)) /= 1) .or. (isnan(real(argC(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' argC HAS', &
						' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
						', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
						' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((size(nsC(:)) /= 1) .or. (isnan(real(nsC(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' nsC HAS', &
						' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
						', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
						' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
				end if

				if (INITIALGRIDflag == 1) then
					if (SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1) == 1) then
						if ((nsnormC(1) /= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%nsnormCT(1)) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%nsnormCT(:)) /= 1) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%nsnormCT(1))) .eqv. &
							.true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' nsnormCT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
								' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
						end if
					end if
				end if

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then

					if ((size(gCNeut(:)) /= 1) .or. (isnan(real(gCNeut(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' gCNeut HAS', &
							' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
							' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((size(argCNeut(:)) /= 1) .or. (isnan(real(argCNeut(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' argCNeut HAS', &
							' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
							' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((size(nsnormCNeut(:)) /= 1) .or. (isnan(real(nsnormCNeut(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' nsnormCNeut HAS', &
							' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
							' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((nsnormCNeut(1) /= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%nsnormCNeut0T(1)) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%nsnormCNeut0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%nsnormCNeut0T(1))) .eqv. &
						.true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' nsnormCNeutGT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
							' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
					end if

				end if

				! ----------------------------------------------------

			end do

			! ----------------------------------------------------

			nsnormCLB(1)= nint(SpecieT(s)%FluxTubeT(f)%QCell0T(SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1)%nsnormCT(1))
			nsnormCUB(1)= nint(SpecieT(s)%FluxTubeT(f)%QCell0T(SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 1)%nsnormCT(1))

			! ----------------------------------------------------

			! DIAGNOSTIC FLAG FOR CONSISTENT LOWER GHOST CELL DENSITIY:

			if (SpecieT(s)%Qindns0GT(1) /= 1d0) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' LOWER GHOST CELL DOES NOT MATCH ALTITUDE OF REFERENCE ', &
					' DENSITY FOR SPECIE= ', s, ', AND FLUX TUBE= ', f, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			!if (nsnormCLB(1) /= (ns0(1)/SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))* &
			!	(SpecieT(s)%FluxTubeT(f)%QCell0T(SpecieT(s)%FluxTubeT(f)%NqLBT(1))%d3xC0T(1))) then
			!	write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
			!		' INCONSISTENT LOWER GHOST CELL DENSITY, nsnormCLB= ', nsnormCLB(1), &
			!		', AND (ns0/nsnormfacT)*d3xCLB= ', &
			!		(ns0(1)/SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))* &
			!		(SpecieT(s)%FluxTubeT(f)%QCell0T(SpecieT(s)%FluxTubeT(f)%NqLBT(1))%d3xC0T(1)), &
			!		' FOR SPECIE= ', s, ', AND FLUX TUBE= ', f, ' IN CONFIGURATION-SPACE GRID', &
			!		' GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			!end if

			! ----------------------------------------------------

		end if

		! ----------------------------------------------------

		! Broadcast boundary ion densities to all ranks
		call mpi_barrier(MPI_COMM_WORLD, ierr)
		call mpi_bcast(nsnormCLB(1), 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
		call mpi_barrier(MPI_COMM_WORLD, ierr)
		call mpi_bcast(nsnormCUB(1), 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

		! ----------------------------------------------------

		! SET BOUNDARY GHOST CELL DENSITIES:

		if (INITIALGRIDflag == 1) then
			allocate(SpecieT(s)%FluxTubeT(f)%d3xCLBGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
			allocate(SpecieT(s)%FluxTubeT(f)%sigmaLBGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
			allocate(SpecieT(s)%FluxTubeT(f)%d3xCUBGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
			allocate(SpecieT(s)%FluxTubeT(f)%sigmaUBGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
			allocate(SpecieT(s)%FluxTubeT(f)%LBNominalDensityGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
			allocate(SpecieT(s)%FluxTubeT(f)%UBNominalDensityGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
			allocate(SpecieT(s)%FluxTubeT(f)%ns0GT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
		end if

		if (SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 0) then
			SpecieT(s)%FluxTubeT(f)%d3xCLBGT(nn)= &
				SpecieT(s)%FluxTubeT(f)%QCell0T(SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1)%d3xC0T(1)
			SpecieT(s)%FluxTubeT(f)%sigmaLBGT(nn)= &
				SpecieT(s)%FluxTubeT(f)%QCell0T(SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1)%d3xC0T(1)/ &
				(SpecieT(s)%FluxTubeT(f)%QCell0T(SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1)%hqC0T(1)* &
				SpecieT(s)%FluxTubeT(f)%QCell0T(SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1)%dqC0T(1))
		end if

		if (SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 1) then
			SpecieT(s)%FluxTubeT(f)%d3xCLBGT(nn)= &
				SpecieT(s)%FluxTubeT(f)%QCell0T(SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1)%d3xC0T(1)
			SpecieT(s)%FluxTubeT(f)%sigmaLBGT(nn)= &
				SpecieT(s)%FluxTubeT(f)%QCell0T(SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1)%d3xC0T(1)/ &
				(SpecieT(s)%FluxTubeT(f)%QCell0T(SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1)%hqC0T(1)* &
				SpecieT(s)%FluxTubeT(f)%QCell0T(SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1)%dqC0T(1))
		end if

		if (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 0) then
			SpecieT(s)%FluxTubeT(f)%d3xCUBGT(nn)= &
				SpecieT(s)%FluxTubeT(f)%QCell0T(SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 1)%d3xC0T(1)
			SpecieT(s)%FluxTubeT(f)%sigmaUBGT(nn)= &
				SpecieT(s)%FluxTubeT(f)%QCell0T(SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 1)%d3xC0T(1)/ &
				(SpecieT(s)%FluxTubeT(f)%QCell0T(SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 1)%hqC0T(1)* &
				SpecieT(s)%FluxTubeT(f)%QCell0T(SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 1)%dqC0T(1))
		end if

		if (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 1) then
			SpecieT(s)%FluxTubeT(f)%d3xCUBGT(nn)= &
				SpecieT(s)%FluxTubeT(f)%QCell0T(SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 1)%d3xC0T(1)
			SpecieT(s)%FluxTubeT(f)%sigmaUBGT(nn)= &
				SpecieT(s)%FluxTubeT(f)%QCell0T(SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 1)%d3xC0T(1)/ &
				(SpecieT(s)%FluxTubeT(f)%QCell0T(SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 1)%hqC0T(1)* &
				SpecieT(s)%FluxTubeT(f)%QCell0T(SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 1)%dqC0T(1))
		end if

		SpecieT(s)%FluxTubeT(f)%LBNominalDensityGT(nn)= nsnormCLB(1)
		SpecieT(s)%FluxTubeT(f)%UBNominalDensityGT(nn)= nsnormCUB(1)
		SpecieT(s)%FluxTubeT(f)%ns0GT(nn)= ns0(1)

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER GHOST CELL DENSITIES:

		if (SpecieT(1)%FluxTubeT(1)%LBCONDITIONflagT(1) == 1) then
			if (SpecieT(s)%FluxTubeT(f)%LBNominalDensityGT(nn) == 0d0) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' ZERO LOWER BOUNDARY DENSITY= ', SpecieT(s)%FluxTubeT(f)%LBNominalDensityGT(nn), &
					' FOR SPECIE= ', s, ', AND STATISTICAL TIME-STEP= ', nn, ', FLUX TUBE= ', f, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if (SpecieT(s)%FluxTubeT(f)%LBNominalDensityGT(nn) == 0d0) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' NOMINAL LOWER BOUNDARY DENSITY= ', SpecieT(s)%FluxTubeT(f)%LBNominalDensityGT(nn), &
					' FOR SPECIE= ', s, ', AND STATISTICAL TIME-STEP= ', nn, ', FLUX TUBE= ', f, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if
		end if

		! ----------------------------------------------------

		! RE-INDEX CONFIGURATION SPACE GRID FOR A NON-COMPUTATIONAL LOWER BOUNDARY GHOST CELL:

		! ----------------------------------------------------
		if (INITIALGRIDflag == 1) then
			allocate(SpecieT(s)%FluxTubeT(f)%qGCGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))
			allocate(SpecieT(s)%FluxTubeT(f)%hqCGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))
			allocate(SpecieT(s)%FluxTubeT(f)%dpCGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))
			allocate(SpecieT(s)%FluxTubeT(f)%dqCGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))
			allocate(SpecieT(s)%FluxTubeT(f)%dphiCGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))
			allocate(SpecieT(s)%FluxTubeT(f)%rGCGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))
			allocate(SpecieT(s)%FluxTubeT(f)%phiGCGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))
			allocate(SpecieT(s)%FluxTubeT(f)%thetaGCGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))
			allocate(SpecieT(s)%FluxTubeT(f)%ellGCGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))
			allocate(SpecieT(s)%FluxTubeT(f)%qGLGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))
			allocate(SpecieT(s)%FluxTubeT(f)%qGHGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))
			allocate(SpecieT(s)%FluxTubeT(f)%pGCGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))
			allocate(SpecieT(s)%FluxTubeT(f)%d3xCGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))
			allocate(SpecieT(s)%FluxTubeT(f)%TsPerpGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))
			allocate(SpecieT(s)%FluxTubeT(f)%TsParGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))
			allocate(SpecieT(s)%FluxTubeT(f)%TsGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))
			allocate(SpecieT(s)%FluxTubeT(f)%nsnormCNeutGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))
			allocate(SpecieT(s)%FluxTubeT(f)%dsICRGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1)))

			do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

				SpecieT(s)%FluxTubeT(f)%qGCGT(1, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%qGC0T(1)
				SpecieT(s)%FluxTubeT(f)%hqCGT(1, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%hqC0T(1)
				SpecieT(s)%FluxTubeT(f)%dpCGT(1, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%dpC0T(1)
				SpecieT(s)%FluxTubeT(f)%dqCGT(1, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%dqC0T(1)
				SpecieT(s)%FluxTubeT(f)%dphiCGT(1, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%dphiC0T(1)
				SpecieT(s)%FluxTubeT(f)%rGCGT(1, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%rGC0T(1)
				SpecieT(s)%FluxTubeT(f)%phiGCGT(1, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%phiGC0T(1)
				SpecieT(s)%FluxTubeT(f)%thetaGCGT(1, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%thetaGC0T(1)
				SpecieT(s)%FluxTubeT(f)%ellGCGT(1, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%ellGC0T(1)
				SpecieT(s)%FluxTubeT(f)%qGLGT(1, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%qGL0T(1)
				SpecieT(s)%FluxTubeT(f)%qGHGT(1, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%qGH0T(1)
				SpecieT(s)%FluxTubeT(f)%pGCGT(1, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%pGC0T(1)
				SpecieT(s)%FluxTubeT(f)%d3xCGT(1, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%d3xC0T(1)
				SpecieT(s)%FluxTubeT(f)%TsPerpGT(1, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%TsPerp0T(1)
				SpecieT(s)%FluxTubeT(f)%TsParGT(1, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%TsPar0T(1)
				SpecieT(s)%FluxTubeT(f)%TsGT(1, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%Ts0T(1)

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
					SpecieT(s)%FluxTubeT(f)%nsnormCNeutGT(1, Qind)= nsnormCNeut0(Qind+ 1)
				end if

				! ----------------------------------------------------

				! Compute field line arc length of BBELF wave-heating region
				SpecieT(s)%FluxTubeT(f)%dsICRGT(1, Qind)= &
					SpecieT(s)%FluxTubeT(f)%hqCGT(1, Qind)*SpecieT(s)%FluxTubeT(f)%dqCGT(1, Qind)

				! ----------------------------------------------------

				! DIAGNOSTIC FLAGS FOR CONSISTENT PHI VALUE:

				if (SpecieT(s)%FluxTubeT(f)%phiGCGT(1, Qind) /= SpecieT(s)%FluxTubeT(f)%QCell0T(1)%phiGC0T(1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT INITIAL phiGCGT VALUE= ', &
						SpecieT(s)%FluxTubeT(f)%phiGCGT(1, Qind), ' AND phiGC0T VALUE= ', &
						SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(1), &
						' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND Q CELL= ', Qind, &
						' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
				end if

				if (SpecieT(s)%FluxTubeT(f)%pGCGT(1, Qind) /= SpecieT(s)%FluxTubeT(f)%QCell0T(1)%pGC0T(1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT INITIAL pGCGT VALUE= ', &
						SpecieT(s)%FluxTubeT(f)%pGCGT(1, Qind), ' AND phiGC0T VALUE= ', &
						SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(1), &
						' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND Q CELL= ', Qind, &
						' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
				end if

				! ----------------------------------------------------

			end do

		end if

		if (INITIALGRIDflag == 0) then
			do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

				SpecieT(s)%FluxTubeT(f)%qGCGT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%qGC0T(1)
				SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%hqC0T(1)
				SpecieT(s)%FluxTubeT(f)%dpCGT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%dpC0T(1)
				SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%dqC0T(1)
				SpecieT(s)%FluxTubeT(f)%dphiCGT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%dphiC0T(1)
				SpecieT(s)%FluxTubeT(f)%rGCGT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%rGC0T(1)
				SpecieT(s)%FluxTubeT(f)%phiGCGT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%phiGC0T(1)
				SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%thetaGC0T(1)
				SpecieT(s)%FluxTubeT(f)%ellGCGT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%ellGC0T(1)
				SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%qGL0T(1)
				SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%qGH0T(1)
				SpecieT(s)%FluxTubeT(f)%pGCGT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%pGC0T(1)
				SpecieT(s)%FluxTubeT(f)%d3xCGT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%d3xC0T(1)
				SpecieT(s)%FluxTubeT(f)%TsPerpGT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%TsPerp0T(1)
				SpecieT(s)%FluxTubeT(f)%TsParGT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%TsPar0T(1)
				SpecieT(s)%FluxTubeT(f)%TsGT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%Ts0T(1)

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
					SpecieT(s)%FluxTubeT(f)%nsnormCNeutGT(nn, Qind)= nsnormCNeut0(Qind+ 1)
				end if

				! ----------------------------------------------------

				! Compute field line arc length of BBELF wave-heating region
				SpecieT(s)%FluxTubeT(f)%dsICRGT(nn, Qind)= &
					SpecieT(s)%FluxTubeT(f)%hqCGT(nn, Qind)*SpecieT(s)%FluxTubeT(f)%dqCGT(nn, Qind)

				! ----------------------------------------------------

				! DIAGNOSTIC FLAGS FOR CONSISTENT PHI VALUE:

				if (SpecieT(s)%FluxTubeT(f)%phiGCGT(nn, Qind) /= SpecieT(s)%FluxTubeT(f)%QCell0T(1)%phiGC0T(1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT INITIAL phiGCGT VALUE= ', &
						SpecieT(s)%FluxTubeT(f)%phiGCGT(nn, Qind), ' AND phiGC0T VALUE= ', &
						SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(1), &
						' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND Q CELL= ', Qind, &
						' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
				end if

				if (SpecieT(s)%FluxTubeT(f)%pGCGT(nn, Qind) /= SpecieT(s)%FluxTubeT(f)%QCell0T(1)%pGC0T(1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT INITIAL pGCGT VALUE= ', &
						SpecieT(s)%FluxTubeT(f)%pGCGT(nn, Qind), ' AND phiGC0T VALUE= ', &
						SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(1), &
						' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND Q CELL= ', Qind, &
						' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
				end if

				! ----------------------------------------------------

			end do

			! ----------------------------------------------------

		end if

		! ----------------------------------------------------

  end subroutine ConfigGridGeneratorSub

end module ConfigGridGenerator
