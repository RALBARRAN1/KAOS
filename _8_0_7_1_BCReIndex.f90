module BCReIndex

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.0.7.1 BOUNDARY CONDITION REINDEX:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

	subroutine BCReIndexSub

    ! ----------------------------------------------------

    ! RE-INDEX TERMS TO ACCOUNT FOR ION ESCAPE:

    ! ----------------------------------------------------

		allocate(BCreal(size(LBoutfluxIonMsk)), BCinteger(size(LBoutfluxIonMsk)), BClogical(size(LBoutfluxIonMsk)))

    qk4(:)= pack(qk4(:), ((LBoutfluxIonMsk(:) .eqv. .false.) .and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    AEAmagN(:)= pack(AEAmagN(:), ((LBoutfluxIonMsk(:) .eqv. .false.) .and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
		AGmagN(:)= pack(AGmagN(:), ((LBoutfluxIonMsk(:) .eqv. .false.) .and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    AEPmagN(:)= pack(AEPmagN(:), ((LBoutfluxIonMsk(:) .eqv. .false.) .and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    Qindk1(:)= pack(Qindk1(:), ((LBoutfluxIonMsk(:) .eqv. .false.) .and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCinteger(:))
    Qindk0(:)= pack(Qindk0(:), ((LBoutfluxIonMsk(:) .eqv. .false.) .and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCinteger(:))
    Vparindk1(:)= pack(Vparindk1(:), ((LBoutfluxIonMsk(:) .eqv. .false.) .and. &
			(UBoutfluxIonMsk(:) .eqv. .false.)), BCinteger(:))

		if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then
      Vperp1indk1(:)= pack(Vperp1indk1(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
				.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCinteger(:))
      Vperp1indk0(:)= pack(Vperp1indk0(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
				.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCinteger(:))
			Vperp2indk1(:)= pack(Vperp2indk1(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
				.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCinteger(:))
			Vperp2indk0(:)= pack(Vperp2indk0(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
				.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCinteger(:))
		else
			Vperpindk1(:)= pack(Vperpindk1(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
				.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCinteger(:))
			Vperpindk0(:)= pack(Vperpindk0(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
				.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCinteger(:))
		end if

		Vparindk0(:)= pack(Vparindk0(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCinteger(:))
		Vpindk1(:)= pack(Vpindk1(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCinteger(:))
    Vqindk1(:)= pack(Vqindk1(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCinteger(:))
    Vphiindk1(:)= pack(Vphiindk1(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCinteger(:))
		VparConvSign(:)= pack(VparConvSign(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
		VxN(:)= pack(VxN(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    VyN(:)= pack(VyN(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    VzN(:)= pack(VzN(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    xN(:)= pack(xN(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    yN(:)= pack(yN(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    zN(:)= pack(zN(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    Vperp1N(:)= pack(Vperp1N(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    Vperp2N(:)= pack(Vperp2N(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    VperpN(:)= pack(VperpN(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    pk4(:)= pack(pk4(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    phik4(:)= pack(phik4(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
		ellk4(:)= pack(ellk4(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    rk4(:)= pack(rk4(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    thetak4(:)= pack(thetak4(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    ENAflag(:)= pack(ENAflag(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BClogical(:))
    ENAflagN0ind(:)= pack(ENAflagN0ind(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCinteger(:))
    ENAflagN1ind(:)= pack(ENAflagN1ind(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCinteger(:))
    AEAmag(:)= pack(AEAmag(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
		AGmag(:)= pack(AGmag(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    AEPmag(:)= pack(AEPmag(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    x(:)= pack(x(:), ((LBoutfluxIonMsk(:) .eqv. .false.) .and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    y(:)= pack(y(:), ((LBoutfluxIonMsk(:) .eqv. .false.) .and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    z(:)= pack(z(:), ((LBoutfluxIonMsk(:) .eqv. .false.) .and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    Vperp1(:)= pack(Vperp1(:), ((LBoutfluxIonMsk(:) .eqv. .false.) .and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    Vperp2(:)= pack(Vperp2(:), ((LBoutfluxIonMsk(:) .eqv. .false.) .and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    Vperp(:)= pack(Vperp(:), ((LBoutfluxIonMsk(:) .eqv. .false.) .and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    Vpar(:)= pack(Vpar(:), ((LBoutfluxIonMsk(:) .eqv. .false.) .and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    Vx(:)= pack(Vx(:), ((LBoutfluxIonMsk(:) .eqv. .false.) .and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    Vy(:)= pack(Vy(:), ((LBoutfluxIonMsk(:) .eqv. .false.) .and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    Vz(:)= pack(Vz(:), ((LBoutfluxIonMsk(:) .eqv. .false.) .and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    Vphi(:)= pack(Vphi(:), ((LBoutfluxIonMsk(:) .eqv. .false.) .and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    Vp(:)= pack(Vp(:), ((LBoutfluxIonMsk(:) .eqv. .false.) .and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    pdriftion(:)= pack(pdriftion(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    qdriftion(:)= pack(qdriftion(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    phidriftion(:)= pack(phidriftion(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    rdriftENA(:)= pack(rdriftENA(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    thetadriftENA(:)= pack(thetadriftENA(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))
    phidriftENA(:)= pack(phidriftENA(:), ((LBoutfluxIonMsk(:) .eqv. .false.) &
			.and. (UBoutfluxIonMsk(:) .eqv. .false.)), BCreal(:))

		deallocate(BCreal, BCinteger, BClogical)

		! ----------------------------------------------------

    ! RE-INDEX TERMS TO ACCOUNT FOR ION AND ENA ESCAPE:

    ! ----------------------------------------------------

		if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then

			allocate(BCreal(size(LBoutfluxENAMsk)), BCinteger(size(LBoutfluxENAMsk)), BClogical(size(LBoutfluxENAMsk)))

      qk4(:)= pack(qk4(:), ((LBoutfluxENAMsk(:) .eqv. .false.) .and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      AEAmagN(:)= pack(AEAmagN(:), ((LBoutfluxENAMsk(:) .eqv. .false.) .and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
			AGmagN(:)= pack(AGmagN(:), ((LBoutfluxENAMsk(:) .eqv. .false.) .and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      AEPmagN(:)= pack(AEPmagN(:), ((LBoutfluxENAMsk(:) .eqv. .false.) .and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      Qindk1(:)= pack(Qindk1(:), ((LBoutfluxENAMsk(:) .eqv. .false.) .and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCinteger(:))
      Qindk0(:)= pack(Qindk0(:), ((LBoutfluxENAMsk(:) .eqv. .false.) .and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCinteger(:))
      Vparindk1(:)= pack(Vparindk1(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCinteger(:))

			if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then
        Vperp1indk1(:)= pack(Vperp1indk1(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
					.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCinteger(:))
        Vperp1indk0(:)= pack(Vperp1indk0(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
					.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCinteger(:))
				Vperp2indk1(:)= pack(Vperp2indk1(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
					.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCinteger(:))
				Vperp2indk0(:)= pack(Vperp2indk0(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
					.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCinteger(:))
			else
				Vperpindk1(:)= pack(Vperpindk1(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
					.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCinteger(:))
				Vperpindk0(:)= pack(Vperpindk0(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
					.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCinteger(:))
			end if

			Vparindk0(:)= pack(Vparindk0(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCinteger(:))
			Vpindk1(:)= pack(Vpindk1(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCinteger(:))
      Vqindk1(:)= pack(Vqindk1(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCinteger(:))
      Vphiindk1(:)= pack(Vphiindk1(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCinteger(:))
			VparConvSign(:)= pack(VparConvSign(:), ((LBoutfluxENAMsk(:) .eqv. .false.) .and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      VxN(:)= pack(VxN(:), ((LBoutfluxENAMsk(:) .eqv. .false.) .and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      VyN(:)= pack(VyN(:), ((LBoutfluxENAMsk(:) .eqv. .false.) .and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      VzN(:)= pack(VzN(:), ((LBoutfluxENAMsk(:) .eqv. .false.) .and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      xN(:)= pack(xN(:), ((LBoutfluxENAMsk(:) .eqv. .false.) .and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      yN(:)= pack(yN(:), ((LBoutfluxENAMsk(:) .eqv. .false.) .and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      zN(:)= pack(zN(:), ((LBoutfluxENAMsk(:) .eqv. .false.) .and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      Vperp1N(:)= pack(Vperp1N(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      Vperp2N(:)= pack(Vperp2N(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      VperpN(:)= pack(VperpN(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      pk4(:)= pack(pk4(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      phik4(:)= pack(phik4(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
			ellk4(:)= pack(ellk4(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      rk4(:)= pack(rk4(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      thetak4(:)= pack(thetak4(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      ENAflag(:)= pack(ENAflag(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BClogical(:))
      ENAflagN0ind(:)= pack(ENAflagN0ind(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCinteger(:))
      ENAflagN1ind(:)= pack(ENAflagN1ind(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCinteger(:))
      AEAmag(:)= pack(AEAmag(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
			AGmag(:)= pack(AGmag(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      AEPmag(:)= pack(AEPmag(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      x(:)= pack(x(:), ((LBoutfluxENAMsk(:) .eqv. .false.) .and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      y(:)= pack(y(:), ((LBoutfluxENAMsk(:) .eqv. .false.) .and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      z(:)= pack(z(:), ((LBoutfluxENAMsk(:) .eqv. .false.) .and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      Vperp1(:)= pack(Vperp1(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      Vperp2(:)= pack(Vperp2(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      Vperp(:)= pack(Vperp(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      Vpar(:)= pack(Vpar(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      Vx(:)= pack(Vx(:), ((LBoutfluxENAMsk(:) .eqv. .false.) .and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      Vy(:)= pack(Vy(:), ((LBoutfluxENAMsk(:) .eqv. .false.) .and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      Vz(:)= pack(Vz(:), ((LBoutfluxENAMsk(:) .eqv. .false.) .and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      Vphi(:)= pack(Vphi(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      Vp(:)= pack(Vp(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      pdriftion(:)= pack(pdriftion(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      qdriftion(:)= pack(qdriftion(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      phidriftion(:)= pack(phidriftion(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      rdriftENA(:)= pack(rdriftENA(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      thetadriftENA(:)= pack(thetadriftENA(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))
      phidriftENA(:)= pack(phidriftENA(:), ((LBoutfluxENAMsk(:) .eqv. .false.) &
				.and. (UBoutfluxENAMsk(:) .eqv. .false.)), BCreal(:))

			deallocate(BCreal, BCinteger, BClogical)

		end if

    ! ----------------------------------------------------

    ! RE-INDEX TERMS TO ACCOUNT FOR LB ION LOSS:

    ! ----------------------------------------------------

    if (SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 1) then

      allocate(BCreal(size(LBreplenishIonMsk)), BCinteger(size(LBreplenishIonMsk)), BClogical(size(LBreplenishIonMsk)))

      qk4(:)= pack(qk4(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      AEAmagN(:)= pack(AEAmagN(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
			AGmagN(:)= pack(AGmagN(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      AEPmagN(:)= pack(AEPmagN(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      Qindk1(:)= pack(Qindk1(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
      Qindk0(:)= pack(Qindk0(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
      Vparindk1(:)= pack(Vparindk1(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))

      if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then
        Vperp1indk1(:)= pack(Vperp1indk1(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
        Vperp1indk0(:)= pack(Vperp1indk0(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
        Vperp2indk1(:)= pack(Vperp2indk1(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
        Vperp2indk0(:)= pack(Vperp2indk0(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
      else
        Vperpindk1(:)= pack(Vperpindk1(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
        Vperpindk0(:)= pack(Vperpindk0(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
      end if

      Vparindk0(:)= pack(Vparindk0(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
      Vpindk1(:)= pack(Vpindk1(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
      Vqindk1(:)= pack(Vqindk1(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
      Vphiindk1(:)= pack(Vphiindk1(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
			VparConvSign(:)= pack(VparConvSign(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      VxN(:)= pack(VxN(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      VyN(:)= pack(VyN(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      VzN(:)= pack(VzN(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      xN(:)= pack(xN(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      yN(:)= pack(yN(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      zN(:)= pack(zN(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      Vperp1N(:)= pack(Vperp1N(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      Vperp2N(:)= pack(Vperp2N(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      VperpN(:)= pack(VperpN(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      pk4(:)= pack(pk4(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      phik4(:)= pack(phik4(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      ellk4(:)= pack(ellk4(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      rk4(:)= pack(rk4(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      thetak4(:)= pack(thetak4(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      ENAflag(:)= pack(ENAflag(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BClogical(:))
      ENAflagN0ind(:)= pack(ENAflagN0ind(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
      ENAflagN1ind(:)= pack(ENAflagN1ind(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
      AEAmag(:)= pack(AEAmag(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
			AGmag(:)= pack(AGmag(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      AEPmag(:)= pack(AEPmag(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      x(:)= pack(x(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      y(:)= pack(y(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      z(:)= pack(z(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      Vperp1(:)= pack(Vperp1(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      Vperp2(:)= pack(Vperp2(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      Vperp(:)= pack(Vperp(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      Vpar(:)= pack(Vpar(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      Vx(:)= pack(Vx(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      Vy(:)= pack(Vy(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      Vz(:)= pack(Vz(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      Vphi(:)= pack(Vphi(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      Vp(:)= pack(Vp(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      pdriftion(:)= pack(pdriftion(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      qdriftion(:)= pack(qdriftion(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      phidriftion(:)= pack(phidriftion(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      rdriftENA(:)= pack(rdriftENA(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      thetadriftENA(:)= pack(thetadriftENA(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      phidriftENA(:)= pack(phidriftENA(:), ((LBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))

      deallocate(BCreal, BCinteger, BClogical)

    end if

    ! ----------------------------------------------------

    ! RE-INDEX TERMS TO ACCOUNT FOR UB ION LOSS:

    ! ----------------------------------------------------

    if (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 1) then

      allocate(BCreal(size(UBreplenishIonMsk)), BCinteger(size(UBreplenishIonMsk)), BClogical(size(UBreplenishIonMsk)))

      qk4(:)= pack(qk4(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      AEAmagN(:)= pack(AEAmagN(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
			AGmagN(:)= pack(AGmagN(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      AEPmagN(:)= pack(AEPmagN(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      Qindk1(:)= pack(Qindk1(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
      Qindk0(:)= pack(Qindk0(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
      Vparindk1(:)= pack(Vparindk1(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))

      if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then
        Vperp1indk1(:)= pack(Vperp1indk1(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
        Vperp1indk0(:)= pack(Vperp1indk0(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
        Vperp2indk1(:)= pack(Vperp2indk1(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
        Vperp2indk0(:)= pack(Vperp2indk0(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
      else
        Vperpindk1(:)= pack(Vperpindk1(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
        Vperpindk0(:)= pack(Vperpindk0(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
      end if

      Vparindk0(:)= pack(Vparindk0(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
      Vpindk1(:)= pack(Vpindk1(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
      Vqindk1(:)= pack(Vqindk1(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
      Vphiindk1(:)= pack(Vphiindk1(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
			VparConvSign(:)= pack(VparConvSign(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      VxN(:)= pack(VxN(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      VyN(:)= pack(VyN(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      VzN(:)= pack(VzN(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      xN(:)= pack(xN(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      yN(:)= pack(yN(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      zN(:)= pack(zN(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      Vperp1N(:)= pack(Vperp1N(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      Vperp2N(:)= pack(Vperp2N(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      VperpN(:)= pack(VperpN(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      pk4(:)= pack(pk4(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      phik4(:)= pack(phik4(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      ellk4(:)= pack(ellk4(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      rk4(:)= pack(rk4(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      thetak4(:)= pack(thetak4(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      ENAflag(:)= pack(ENAflag(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BClogical(:))
      ENAflagN0ind(:)= pack(ENAflagN0ind(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
      ENAflagN1ind(:)= pack(ENAflagN1ind(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCinteger(:))
      AEAmag(:)= pack(AEAmag(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
			AGmag(:)= pack(AGmag(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      AEPmag(:)= pack(AEPmag(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      x(:)= pack(x(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      y(:)= pack(y(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      z(:)= pack(z(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      Vperp1(:)= pack(Vperp1(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      Vperp2(:)= pack(Vperp2(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      Vperp(:)= pack(Vperp(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      Vpar(:)= pack(Vpar(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      Vx(:)= pack(Vx(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      Vy(:)= pack(Vy(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      Vz(:)= pack(Vz(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      Vphi(:)= pack(Vphi(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      Vp(:)= pack(Vp(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      pdriftion(:)= pack(pdriftion(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      qdriftion(:)= pack(qdriftion(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      phidriftion(:)= pack(phidriftion(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      rdriftENA(:)= pack(rdriftENA(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      thetadriftENA(:)= pack(thetadriftENA(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))
      phidriftENA(:)= pack(phidriftENA(:), ((UBreplenishIonMsk(:) .eqv. .false.)), BCreal(:))

      deallocate(BCreal, BCinteger, BClogical)

    end if

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine BCReIndexSub

end module BCReIndex
