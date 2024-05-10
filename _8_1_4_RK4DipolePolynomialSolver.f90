module RK4DipolePolynomialSolver

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.1.4 - RK4 DIPOLE POLYNOMIAL SOLVER:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! SOLVE QUARTIC DIPOLE POLYNOMIAL BY ANALYTIC ORDER DECOMPOSITION FOR UPDATED PARTICLE
! POSITIONS IN RECTANGULAR CONFIG. SPACE:

	subroutine RK4DipolePolynomialSolverSub

  	! ----------------------------------------------------

  	ApRK4(1)= 0d0  ! From dipole quartic polynomial form
  	BpRK4(1)= 0d0
  	CpRK4(1)= 1d0/(pNp(1)*qNp(1)**2d0)
  	DpRK4(1)= -1d0/(qNp(1)**2d0)

  	AbRK4(1)= -BpRK4(1) ! From resolvent cubic form
  	BbRK4(1)= ApRK4(1)*CpRK4(1)- 4d0*DpRK4(1)
  	CbRK4(1)= 4d0*BpRK4(1)*DpRK4(1)- CpRK4(1)**2d0- (ApRK4(1)**2d0)*DpRK4(1)

  	! Resolvent cubic discriminant
  	D3RK4(1)= (BbRK4(1)**2d0)*(AbRK4(1)**2d0)- 4d0*CbRK4(1)*(AbRK4(1)**3d0)- &
  		4d0*(BbRK4(1)**3d0)+ 18d0*AbRK4(1)* BbRK4(1)*CbRK4(1)- 27d0*(CbRK4(1)**2d0)

  	D4RK4(1)= ((CpRK4(1)**2d0)*(BpRK4(1)**2d0)*(ApRK4(1)**2d0)- 4d0*(CpRK4(1)**3d0)* &
  		(ApRK4(1)**3d0)- 4d0*(CpRK4(1)**2d0)*(BpRK4(1)**3d0)+ 18d0*(CpRK4(1)**3d0)* &
  		BpRK4(1)*ApRK4(1)- 27d0*(CpRK4(1)**4d0)+ 256d0*(DpRK4(1)**3d0))+ DpRK4(1)* &
  		(-4d0*(BpRK4(1)**3d0)*(ApRK4(1)**2d0)+ 18d0*CpRK4(1)*BpRK4(1)*(ApRK4(1)**3d0)+ &
  		16d0*(BpRK4(1)**4d0)- 80d0*CpRK4(1)*(BpRK4(1)**2d0)*ApRK4(1)- 6d0*(CpRK4(1)**2d0)* &
  		(ApRK4(1)**2d0)+ 144d0*(CpRK4(1)**2d0)*BpRK4(1))+ (DpRK4(1)**2d0)*(-27d0* &
  		(ApRK4(1)**4d0)+ 144d0*BpRK4(1)*(ApRK4(1)**2d0)- 128d0*(BpRK4(1)**2d0)- &
  		192d0*CpRK4(1)*ApRK4(1)) ! Original quartic disrcriminant

  	if (abs(D3RK4(1)- D4RK4(1)) > 1d-14) then
  		write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' BAD DISCRIMINANT', &
  			' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND PARTICLE= ', j, &
        ' IN RK4 DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
  	end if

  	AbbRK4(1)= BbRK4(1)- (AbRK4(1)**2d0)/3d0
  	BbbRK4(1)= 2d0*(AbRK4(1)**3d0)/27d0- AbRK4(1)*BbRK4(1)/3d0+ CbRK4(1)
  	DeltaRK4(1)= -BbbRK4(1)/2d0
  	EpsilonRK4(1)= AbbRK4(1)/3d0

  	ThetapRK4(1)= acos(DeltaRK4(1)/(i*sqrt(EpsilonRK4(1)**3d0)))

  	sigma23RK4(1)= 2d0*i*sqrt(EpsilonRK4(1))*cos((ThetapRK4(1)+ 4d0*pi)/3d0)- &
  		AbRK4(1)/3d0
  	sigmaRK4(1)= (abs(real(sigma23RK4(1))))

  	muRK4(1)= sqrt((ApRK4(1)**2d0)/4d0- BpRK4(1)+ sigmaRK4(1))

  	if (muRK4(1) /= 0) then
  		nuRK4(1)= sqrt(3d0*(ApRK4(1)**2d0)/4d0- (muRK4(1)**2d0)- 2d0*BpRK4(1)+ &
  			(4d0*ApRK4(1)*BpRK4(1)- 8d0*CpRK4(1)- (ApRK4(1)**3d0))/(4d0*muRK4(1)))
  		piiRK4(1)= sqrt(3d0*(ApRK4(1)**2d0)/4d0- (muRK4(1)**2d0)- 2d0*BpRK4(1)- &
  			(4d0*ApRK4(1)*BpRK4(1)- 8d0*CpRK4(1)- (ApRK4(1)**3d0))/(4d0*muRK4(1)))
  	else if (muRK4(1) == 0) then
  		nuRK4(1)= sqrt(3d0*(ApRK4(1)**2d0)/4d0- 2d0*BpRK4(1)+ &
  			2d0*sqrt((sigmaRK4(1)**2d0)- 4d0*DpRK4(1)))
  		piiRK4(1)= sqrt(3d0*(ApRK4(1)**2d0)/4d0- 2d0*BpRK4(1)- &
  			2d0*sqrt((sigmaRK4(1)**2d0)- 4d0*DpRK4(1)))
  	end if

  	! Original real quartic root > 0
  	gamma3RK4(1)= -ApRK4(1)/4d0- muRK4(1)/2d0+ piiRK4(1)/2d0

  	! Revert quartic roots with (q, p) into (r, theta)

  	r3RK4(1)= gamma3RK4(1)*RE
  	theta3RK4(1)= asin(sqrt(r3RK4(1) &
  		/(RE*pNp(1))))

  	qtest3RK4(1)= (RE**2d0)*cos(theta3RK4(1))/(r3RK4(1)**2d0)

  	ptest3RK4(1)= r3RK4(1)/(RE*(sin(theta3RK4(1))**2d0))

  	! Note: gamma3RK4 is real positive root. For q< 0 (q> 0), phase-shift thetaRK4 by
  	! pi (0) to get correct sign of q (above and below dipole equator).

  	rfinalRK4(1)= r3RK4(1) ! Get final (r, theta) values

  	if (qNp(1) <= 0) then ! Phase shift root 3 solution
  		thetafinalRK4(1)= pi- theta3RK4(1)
  	else if (qNp(1) > 0) then
  		thetafinalRK4(1)= theta3RK4(1)
  	end if

  	! Note: Get final (y, z) values.

  	phifinalRK4(1)= phiNp(1)

  	xfinalRK4(1)= rfinalRK4(1)*sin(thetafinalRK4(1))*cos(phifinalRK4(1))
  	yfinalRK4(1)= rfinalRK4(1)*sin(thetafinalRK4(1))*sin(phifinalRK4(1))
  	zfinalRK4(1)= rfinalRK4(1)*cos(thetafinalRK4(1))

  	if (abs(xfinalRK4(1)) < 1d-6) then ! Set Cartesian coords according to B longitude
  		xfinalRK4(1)= 0d0
  	end if

  	if (abs(yfinalRK4(1)) < 1d-6) then
  		yfinalRK4(1)= 0d0
  	end if

  	if (abs(zfinalRK4(1)) < 1d-6) then
  		zfinalRK4(1)= 0d0
  	end if

  	! Note: Get final (q, p) values.

  	qfinalRK4(1)= (RE**2d0)*cos(thetafinalRK4(1))/(rfinalRK4(1)**2d0)
  	pfinalRK4(1)= rfinalRK4(1)/(RE*(sin(thetafinalRK4(1))**2d0))

		if (abs(pNp(1)- pfinalRK4(1)) > 1d-6) then
			pfinalRK4(1)= pNp(1)
		end if

  	! ----------------------------------------------------

  	! DIAGNOSTIC FLAGS IN ABSOLUTE ERROR FOR QUARTIC DIPOLE POLYNOMIAL ROOTS:

  	if (abs(pNp(1)- pfinalRK4(1)) > 1d-6) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCORRECT CONVERSION pNp= ', pNp(1), &
				', pfinalRK4= ', pfinalRK4(1), ' FOR SPECIE= ', s, &
				', FLUX TUBE= ', f, ', AND PARTICLE= ', j, &
				' IN RK4 DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
  	end if

  	!if (abs(qNp(1)- qfinalRK4(1)) > 1d0) then
    !  write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCORRECT CONVERSION qNp= ', qNp(1), &
		!		', qfinalRK4= ', qfinalRK4(1), ' FOR SPECIE= ', s, &
		!		', FLUX TUBE= ', f, ', AND PARTICLE= ', j, &
	  !   ' IN RK4 DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
  	!end if

  	! ----------------------------------------------------

  	! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSION, SIZES AND
  	! FINITE VALUES:

		if ((isnan(real(rfinalRK4(1))) .eqv. .true.) .or. &
  		(size(rfinalRK4(:)) /= 1)) then
  		write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rfinalRK4 HAS', &
  			' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
  			', FLUX TUBE= ', f, ', AND PARTICLE= ', j, &
        ' IN RK4 DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
  	end if

		if ((isnan(real(thetafinalRK4(1))) .eqv. .true.) .or. &
  		(size(thetafinalRK4(:)) /= 1)) then
  		write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' thetafinalRK4 HAS', &
  			' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
  			', FLUX TUBE= ', f, ', AND PARTICLE= ', j, &
        ' IN RK4 DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
  	end if

		if ((isnan(real(phifinalRK4(1))) .eqv. .true.) .or. &
  		(size(phifinalRK4(:)) /= 1)) then
  		write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' phifinalRK4 HAS', &
  			' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
  			', FLUX TUBE= ', f, ', AND PARTICLE= ', j, &
        ' IN RK4 DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
  	end if

  	if ((isnan(real(xfinalRK4(1))) .eqv. .true.) .or. &
  		(size(xfinalRK4(:)) /= 1)) then
  		write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' xfinalRK4 HAS', &
  			' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
  			', FLUX TUBE= ', f, ', AND PARTICLE= ', j, &
        ' IN RK4 DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
  	end if

    if ((isnan(real(yfinalRK4(1))) .eqv. .true.) .or. &
  		(size(yfinalRK4(:)) /= 1)) then
  		write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' yfinalRK4 HAS', &
  			' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
  			', FLUX TUBE= ', f, ', AND PARTICLE= ', j, &
        ' IN RK4 DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
  	end if

    if ((isnan(real(zfinalRK4(1))) .eqv. .true.) .or. &
  		(size(zfinalRK4(:)) /= 1)) then
  		write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' zfinalRK4 HAS', &
  			' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
  			', FLUX TUBE= ', f, ', AND PARTICLE= ', j, &
        ' IN RK4 DIPOLE POLYNOMIAL SOLVER SUBROUTINE' // achar(27) // '[0m.'
  	end if

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine RK4DipolePolynomialSolverSub

end module RK4DipolePolynomialSolver
