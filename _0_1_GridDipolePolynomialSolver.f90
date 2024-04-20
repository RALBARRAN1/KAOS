module GridDipolePolynomialSolver

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	0.1 - GRID DIPOLE POLYNOMIAL SOLVER:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! SOLVE QUARTIC DIPOLE POLYNOMIAL BY ANALYTIC ORDER DECOMPOSITION:

  subroutine GridDipolePolynomialSolverSub

  	! ----------------------------------------------------

  	ApGrid(1)= 0d0  ! From dipole quartic polynomial form
  	BpGrid(1)= 0d0
  	CpGrid(1)= 1d0/(pGridIn(1)*qGridIn(1)**2d0)
  	DpGrid(1)= -1d0/(qGridIn(1)**2d0)

  	AbGrid(1)= -BpGrid(1) ! From resolvent cubic form
  	BbGrid(1)= ApGrid(1)*CpGrid(1)- 4d0*DpGrid(1)
  	CbGrid(1)= 4d0*BpGrid(1)*DpGrid(1)- CpGrid(1)**2d0- (ApGrid(1)**2d0)*DpGrid(1)

  	! Resolvent cubic discriminant
  	D3Grid(1)= (BbGrid(1)**2d0)*(AbGrid(1)**2d0)- 4d0*CbGrid(1)*(AbGrid(1)**3d0)- &
  		4d0*(BbGrid(1)**3d0)+ 18d0*AbGrid(1)* BbGrid(1)*CbGrid(1)- 27d0*(CbGrid(1)**2d0)

  	D4Grid(1)= ((CpGrid(1)**2d0)*(BpGrid(1)**2d0)*(ApGrid(1)**2d0)- 4d0*(CpGrid(1)**3d0)* &
  		(ApGrid(1)**3d0)- 4d0*(CpGrid(1)**2d0)*(BpGrid(1)**3d0)+ 18d0*(CpGrid(1)**3d0)* &
  		BpGrid(1)*ApGrid(1)- 27d0*(CpGrid(1)**4d0)+ 256d0*(DpGrid(1)**3d0))+ DpGrid(1)* &
  		(-4d0*(BpGrid(1)**3d0)*(ApGrid(1)**2d0)+ 18d0*CpGrid(1)*BpGrid(1)*(ApGrid(1)**3d0)+ &
  		16d0*(BpGrid(1)**4d0)- 80d0*CpGrid(1)*(BpGrid(1)**2d0)*ApGrid(1)- 6d0*(CpGrid(1)**2d0)* &
  		(ApGrid(1)**2d0)+ 144d0*(CpGrid(1)**2d0)*BpGrid(1))+ (DpGrid(1)**2d0)*(-27d0* &
  		(ApGrid(1)**4d0)+ 144d0*BpGrid(1)*(ApGrid(1)**2d0)- 128d0*(BpGrid(1)**2d0)- &
  		192d0*CpGrid(1)*ApGrid(1)) ! Original quartic disrcriminant

  	if (abs(D3Grid(1)- D4Grid(1)) > 1d-14) then
  		write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' BAD DISCRIMINANT', &
  			' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
        ' IN GRID DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
  	end if

  	AbbGrid(1)= BbGrid(1)- (AbGrid(1)**2d0)/3d0
  	BbbGrid(1)= 2d0*(AbGrid(1)**3d0)/27d0- AbGrid(1)*BbGrid(1)/3d0+ CbGrid(1)
  	DeltaGrid(1)= -BbbGrid(1)/2d0
  	EpsilonGrid(1)= AbbGrid(1)/3d0

  	ThetapGrid(1)= acos(DeltaGrid(1)/(i*sqrt(EpsilonGrid(1)**3d0)))

  	sigma23Grid(1)= 2d0*i*sqrt(EpsilonGrid(1))*cos((ThetapGrid(1)+ 4d0*pi)/3d0)- &
  		AbGrid(1)/3d0
  	sigmaGrid(1)= (abs(real(sigma23Grid(1))))

  	muGrid(1)= sqrt((ApGrid(1)**2d0)/4d0- BpGrid(1)+ sigmaGrid(1))

  	if (muGrid(1) /= 0) then
  		nuGrid(1)= sqrt(3d0*(ApGrid(1)**2d0)/4d0- (muGrid(1)**2d0)- 2d0*BpGrid(1)+ &
  			(4d0*ApGrid(1)*BpGrid(1)- 8d0*CpGrid(1)- (ApGrid(1)**3d0))/(4d0*muGrid(1)))
  		piiGrid(1)= sqrt(3d0*(ApGrid(1)**2d0)/4d0- (muGrid(1)**2d0)- 2d0*BpGrid(1)- &
  			(4d0*ApGrid(1)*BpGrid(1)- 8d0*CpGrid(1)- (ApGrid(1)**3d0))/(4d0*muGrid(1)))
  	else if (muGrid(1) == 0) then
  		nuGrid(1)= sqrt(3d0*(ApGrid(1)**2d0)/4d0- 2d0*BpGrid(1)+ &
  			2d0*sqrt((sigmaGrid(1)**2d0)- 4d0*DpGrid(1)))
  		piiGrid(1)= sqrt(3d0*(ApGrid(1)**2d0)/4d0- 2d0*BpGrid(1)- &
  			2d0*sqrt((sigmaGrid(1)**2d0)- 4d0*DpGrid(1)))
  	end if

  	! Original real quartic root > 0
  	gamma3Grid(1)= -ApGrid(1)/4d0- muGrid(1)/2d0+ piiGrid(1)/2d0

  	! Revert quartic roots with (q, p) into (r, theta)

  	r3Grid(1)= gamma3Grid(1)*RE
  	theta3Grid(1)= asin(sqrt(r3Grid(1) &
  		/(RE*pGridIn(1))))

  	qtest3Grid(1)= (RE**2d0)*cos(theta3Grid(1))/(r3Grid(1)**2d0)

  	ptest3Grid(1)= r3Grid(1)/(RE*(sin(theta3Grid(1))**2d0))

  	! Note: gamma3Grid is real positive root. For q< 0 (q> 0), phase-shift thetaGrid by
  	! pi (0) to get correct sign of q (above and below dipole equator).

  	rGridOut(1)= r3Grid(1) ! Get final (r, theta) values

  	if (qGridIn(1) <= 0) then ! Phase shift root 3 solution
  		thetaGridOut(1)= pi- theta3Grid(1)
  	else if (qGridIn(1) > 0) then
  		thetaGridOut(1)= theta3Grid(1)
  	end if

  	! Note: Get final (y, z) values with phi= pi/2 s.t. x= 0.

    ! FIXME Set phiGridOut according to initial input phi with ExB drift

  	phiGridOut(1)= pi/2d0 ! Select B longitude

  	xGridOut(1)= rGridOut(1)*sin(thetaGridOut(1))*cos(phiGridOut(1))
  	yGridOut(1)= rGridOut(1)*sin(thetaGridOut(1))*sin(phiGridOut(1))
  	zGridOut(1)= rGridOut(1)*cos(thetaGridOut(1))

  	if (abs(xGridOut(1)) < 1d-6) then ! Set Cartesian coords according to B longitude
  		xGridOut(1)= 0d0
  	end if

  	if (abs(yGridOut(1)) < 1d-6) then
  		yGridOut(1)= 0d0
  	end if

  	if (abs(zGridOut(1)) < 1d-6) then
  		zGridOut(1)= 0d0
  	end if

  	! Note: Get final (q, p) value and let phid= phi to compare with initial input.

  	qGridOut(1)= (RE**2d0)*cos(thetaGridOut(1))/(rGridOut(1)**2d0)
  	pGridOut(1)= rGridOut(1)/(RE*(sin(thetaGridOut(1))**2d0))

		if (abs(pGridIn(1)- pGridOut(1)) > 1d-6) then
			pGridOut(1)= pGridIn(1)
		end if

  	! ----------------------------------------------------

  	! DIAGNOSTIC FLAGS IN ABSOLUTE ERROR FOR QUARTIC DIPOLE POLYNOMIAL ROOTS:

  	if (abs(pGridIn(1)- pGridOut(1)) > 1d-6) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCORRECT CONVERSION pGridIn= ', pGridIn(1), &
        ', pGridOut= ', pGridOut(1), ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
        ' IN GRID DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
  	end if

  	if (abs(qGridIn(1)- qGridOut(1)) > 1d0) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCORRECT CONVERSION qGridIn= ', qGridIn(1), &
        ', qGridOut= ', qGridOut(1), ' thetaGridOut= ', thetaGridOut(1), ' rGridOut= ', rGridOut(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
        ' IN GRID DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
  	end if

  	! ----------------------------------------------------

  	! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSION, SIZES AND
  	! FINITE VALUES:

  	if ((isnan(real(xGridOut(1))) .eqv. .true.) .or. &
  		(size(xGridOut(:)) /= 1)) then
  		write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' xGridOut HAS', &
  			' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
  			', FLUX TUBE= ', f, ', AND PARTICLE= ', j, &
        ' IN GRID DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
  	end if

    if ((isnan(real(yGridOut(1))) .eqv. .true.) .or. &
  		(size(yGridOut(:)) /= 1)) then
  		write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' yGridOut HAS', &
  			' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
  			', FLUX TUBE= ', f, ', AND PARTICLE= ', j, &
        ' IN GRID DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
  	end if

    if ((isnan(real(zGridOut(1))) .eqv. .true.) .or. &
  		(size(zGridOut(:)) /= 1)) then
  		write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' zGridOut HAS', &
  			' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
  			', FLUX TUBE= ', f, ', AND PARTICLE= ', j, &
        ' IN GRID DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
  	end if

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine GridDipolePolynomialSolverSub

end module GridDipolePolynomialSolver
