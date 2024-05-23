module DataTypeReAllocA

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.0.6- DATA TYPE RE-ALLOCATION A:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

	subroutine DataTypeReAllocASub

    ! ----------------------------------------------------

    ! ALLOCATE TEMPORARY KINETIC ARRAYS:

    ! ----------------------------------------------------

    allocate(LBoutfluxIonMskTMP(NsTK(1)- dNsTK(1)), UBoutfluxIonMskTMP(NsTK(1)- dNsTK(1)), &
			AEAmagNTMP(NsTK(1)- dNsTK(1)), AGmagNTMP(NsTK(1)- dNsTK(1)), AEPmagNTMP(NsTK(1)- dNsTK(1)))

		if (SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 1) then
			allocate(LBreplenishIonMskTMP(NsTK(1)- dNsTK(1)))
		end if
		if (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 1) then
			allocate(UBreplenishIonMskTMP(NsTK(1)- dNsTK(1)))
		end if

		if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
			allocate(LBoutfluxENAMskTMP(NsTK(1)- dNsTK(1)), UBoutfluxENAMskTMP(NsTK(1)- dNsTK(1)))
		end if

		allocate(Qindk1TMP(NsTK(1)- dNsTK(1)), Qindk0TMP(NsTK(1)- dNsTK(1)), &
      Vparindk1TMP(NsTK(1)- dNsTK(1)), Vparindk0TMP(NsTK(1)- dNsTK(1)), &
      Vpindk1TMP(NsTK(1)- dNsTK(1)), Vqindk1TMP(NsTK(1)- dNsTK(1)), &
      Vphiindk1TMP(NsTK(1)- dNsTK(1)))

		if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

			allocate(Vperp1indk1TMP(NsTK(1)- dNsTK(1)), Vperp2indk1TMP(NsTK(1)- dNsTK(1)), &
				Vperp1indk0TMP(NsTK(1)- dNsTK(1)), Vperp2indk0TMP(NsTK(1)- dNsTK(1)))

		else

			allocate(Vperpindk1TMP(NsTK(1)- dNsTK(1)), Vperpindk0TMP(NsTK(1)- dNsTK(1)))

		end if

		allocate(VxNTMP(NsTK(1)- dNsTK(1)), VyNTMP(NsTK(1)- dNsTK(1)), &
      VzNTMP(NsTK(1)- dNsTK(1)), xNTMP(NsTK(1)- dNsTK(1)), &
      yNTMP(NsTK(1)- dNsTK(1)), zNTMP(NsTK(1)- dNsTK(1)), &
			Vperp1NTMP(NsTK(1)- dNsTK(1)), Vperp2NTMP(NsTK(1)- dNsTK(1)), &
      VperpNTMP(NsTK(1)- dNsTK(1)))
		allocate(pk4TMP(NsTK(1)- dNsTK(1)), phik4TMP(NsTK(1)- dNsTK(1)), &
      qk4TMP(NsTK(1)- dNsTK(1)), rk4TMP(NsTK(1)- dNsTK(1)), &
      thetak4TMP(NsTK(1)- dNsTK(1)), ellk4TMP(NsTK(1)- dNsTK(1)))
		allocate(ENAflagTMP(NsTK(1)- dNsTK(1)), ENAflagN0indTMP(NsTK(1)- dNsTK(1)), &
      ENAflagN1indTMP(NsTK(1)- dNsTK(1)))
    allocate(AEAmagTMP(NsTK(1)- dNsTK(1)), AGmagTMP(NsTK(1)- dNsTK(1)), AEPmagTMP(NsTK(1)- dNsTK(1)), &
      xTMP(NsTK(1)- dNsTK(1)), yTMP(NsTK(1)- dNsTK(1)), zTMP(NsTK(1)- dNsTK(1)), &
      Vperp1TMP(NsTK(1)- dNsTK(1)), Vperp2TMP(NsTK(1)- dNsTK(1)), &
      VperpTMP(NsTK(1)- dNsTK(1)), VparTMP(NsTK(1)- dNsTK(1)), &
      VxTMP(NsTK(1)- dNsTK(1)), VyTMP(NsTK(1)- dNsTK(1)), &
      VzTMP(NsTK(1)- dNsTK(1)), VphiTMP(NsTK(1)- dNsTK(1)), VpTMP(NsTK(1)- dNsTK(1)))
    allocate(pdriftionTMP(NsTK(1)- dNsTK(1)), qdriftionTMP(NsTK(1)- dNsTK(1)), &
      phidriftionTMP(NsTK(1)- dNsTK(1)), rdriftENATMP(NsTK(1)- dNsTK(1)), &
      thetadriftENATMP(NsTK(1)- dNsTK(1)), phidriftENATMP(NsTK(1)- dNsTK(1)))

    ! ----------------------------------------------------

    ! COPY ONTO TEMPORARY KINETIC ARRAYS:

    ! ----------------------------------------------------

		LBoutfluxIonMskTMP(:)= LBoutfluxIonMsk(:)
		UBoutfluxIonMskTMP(:)= UBoutfluxIonMsk(:)

		if (SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 1) then
			LBreplenishIonMskTMP(:)= LBreplenishIonMsk(:)
		end if
		if (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 1) then
			UBreplenishIonMskTMP(:)= UBreplenishIonMsk(:)
		end if

		if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
			LBoutfluxENAMskTMP(:)= LBoutfluxENAMsk(:)
			UBoutfluxENAMskTMP(:)= UBoutfluxENAMsk(:)
		end if

    AEAmagNTMP(:)= AEAmagN(:)
		AGmagNTMP(:)= AGmagN(:)
    AEPmagNTMP(:)= AEPmagN(:)
    Qindk1TMP(:)= Qindk1(:)
    Qindk0TMP(:)= Qindk0(:)
    Vparindk1TMP(:)= Vparindk1(:)

		if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

	    Vperp1indk1TMP(:)= Vperp1indk1(:)
			Vperp1indk0TMP(:)= Vperp1indk0(:)
			Vperp2indk1TMP(:)= Vperp2indk1(:)
			Vperp2indk0TMP(:)= Vperp2indk0(:)

		else

			Vperpindk1TMP(:)= Vperpindk1(:)
			Vperpindk0TMP(:)= Vperpindk0(:)

		end if

    Vparindk0TMP(:)= Vparindk0(:)
    Vpindk1TMP(:)= Vpindk1(:)
    Vqindk1TMP(:)= Vqindk1(:)
    Vphiindk1TMP(:)= Vphiindk1(:)
    VxNTMP(:)= VxN(:)
    VyNTMP(:)= VyN(:)
    VzNTMP(:)= VzN(:)
    xNTMP(:)= xN(:)
    yNTMP(:)= yN(:)
    zNTMP(:)= zN(:)
    Vperp1NTMP(:)= Vperp1N(:)
    Vperp2NTMP(:)= Vperp2N(:)
    VperpNTMP(:)= VperpN(:)
    pk4TMP(:)= pk4(:)
    phik4TMP(:)= phik4(:)
    qk4TMP(:)= qk4(:)
    rk4TMP(:)= rk4(:)
    thetak4TMP(:)= thetak4(:)
		ellk4TMP(:)= ellk4(:)
    ENAflagTMP(:)= ENAflag(:)
    ENAflagN0indTMP(:)= ENAflagN0ind(:)
    ENAflagN1indTMP(:)= ENAflagN1ind(:)
    AEAmagTMP(:)= AEAmag(:)
		AGmagTMP(:)= AGmag(:)
    AEPmagTMP(:)= AEPmag(:)
    xTMP(:)= x(:)
    yTMP(:)= y(:)
    zTMP(:)= z(:)
    Vperp1TMP(:)= Vperp1(:)
    Vperp2TMP(:)= Vperp2(:)
    VperpTMP(:)= Vperp(:)
    VparTMP(:)= Vpar(:)
    VxTMP(:)= Vx(:)
    VyTMP(:)= Vy(:)
    VzTMP(:)= Vz(:)
    VphiTMP(:)= Vphi(:)
    VpTMP(:)= Vp(:)
    pdriftionTMP(:)= pdriftion(:)
    qdriftionTMP(:)= qdriftion(:)
    phidriftionTMP(:)= phidriftion(:)
    rdriftENATMP(:)= rdriftENA(:)
    thetadriftENATMP(:)= thetadriftENA(:)
    phidriftENATMP(:)= phidriftENA(:)

    ! ----------------------------------------------------

    ! DE-ALLOCATE KINETIC ARRAYS:

    ! ----------------------------------------------------

    deallocate(LBoutfluxIonMsk, UBoutfluxIonMsk, AEAmagN, AGmagN, AEPmagN)

		if (SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 1) then
			deallocate(LBreplenishIonMsk)
		end if
		if (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 1) then
			deallocate(UBreplenishIonMsk)
		end if

		if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
			deallocate(LBoutfluxENAMsk, UBoutfluxENAMsk)
		end if

		deallocate(Qindk1, Qindk0, Vparindk1, Vparindk0, Vpindk1, Vqindk1, Vphiindk1)

		if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

			deallocate(Vperp1indk1, Vperp1indk0)
			deallocate(Vperp2indk1, Vperp2indk0)

		else

			deallocate(Vperpindk1, Vperpindk0)

		end if

		deallocate(VxN, VyN, VzN, xN, yN, zN, &
			Vperp1N, Vperp2N, VperpN)
		deallocate(pk4, phik4, qk4, rk4, thetak4, ellk4)
		deallocate(ENAflag, ENAflagN0ind, ENAflagN1ind)
    deallocate(AEAmag, AGmag, AEPmag, x, y, z, &
      Vperp1, Vperp2, Vperp, Vpar, &
      Vx, Vy, Vz, Vphi, Vp)
    deallocate(pdriftion, qdriftion, phidriftion, &
      rdriftENA, thetadriftENA, phidriftENA)

    ! ----------------------------------------------------

    ! ALLOCATE KINETIC ARRAYS WITH NEW SIZE:

    ! ----------------------------------------------------

		allocate(LBoutfluxIonMsk(NsTK(1)), UBoutfluxIonMsk(NsTK(1)), AEAmagN(NsTK(1)), AGmagN(NsTK(1)), AEPmagN(NsTK(1)))

		if (SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 1) then
			allocate(LBreplenishIonMsk(NsTK(1)))
		end if
		if (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 1) then
			allocate(UBreplenishIonMsk(NsTK(1)))
		end if

		if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
			allocate(LBoutfluxENAMsk(NsTK(1)), UBoutfluxENAMsk(NsTK(1)))
		end if

		allocate(Qindk1(NsTK(1)), Qindk0(NsTK(1)), Vparindk1(NsTK(1)), &
			Vparindk0(NsTK(1)), Vpindk1(NsTK(1)), Vqindk1(NsTK(1)), Vphiindk1(NsTK(1)))

		if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

			allocate(Vperp1indk1(NsTK(1)), Vperp1indk0(NsTK(1)))
			allocate(Vperp2indk1(NsTK(1)), Vperp2indk0(NsTK(1)))

		else

			allocate(Vperpindk1(NsTK(1)), Vperpindk0(NsTK(1)))

		end if

		allocate(VxN(NsTK(1)), VyN(NsTK(1)), VzN(NsTK(1)), xN(NsTK(1)), yN(NsTK(1)), zN(NsTK(1)), &
			Vperp1N(NsTK(1)), Vperp2N(NsTK(1)), VperpN(NsTK(1)))
		allocate(pk4(NsTK(1)), phik4(NsTK(1)), qk4(NsTK(1)), rk4(NsTK(1)), thetak4(NsTK(1)), ellk4(NsTK(1)))
		allocate(ENAflag(NsTK(1)), ENAflagN0ind(NsTK(1)), ENAflagN1ind(NsTK(1)))
    allocate(AEAmag(NsTK(1)), AGmag(NsTK(1)), AEPmag(NsTK(1)), x(NsTK(1)), y(NsTK(1)), z(NsTK(1)), &
      Vperp1(NsTK(1)), Vperp2(NsTK(1)), Vperp(NsTK(1)), Vpar(NsTK(1)), &
      Vx(NsTK(1)), Vy(NsTK(1)), Vz(NsTK(1)), Vphi(NsTK(1)), Vp(NsTK(1)))
    allocate(pdriftion(NsTK(1)), qdriftion(NsTK(1)), phidriftion(NsTK(1)), &
      rdriftENA(NsTK(1)), thetadriftENA(NsTK(1)), phidriftENA(NsTK(1)))

    ! ----------------------------------------------------

    ! COPY BACK ONTO KINETIC ARRAYS:

    ! ----------------------------------------------------

    do j= 1, NsTK(1), 1
      if ((dNsTK(1) >= 0) .and. (j <= (NsTK(1)- dNsTK(1)))) then
				LBoutfluxIonMsk(j)= LBoutfluxIonMskTMP(j)
				UBoutfluxIonMsk(j)= UBoutfluxIonMskTMP(j)

				if (SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 1) then
					LBreplenishIonMsk(j)= LBreplenishIonMskTMP(j)
				end if
				if (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 1) then
					UBreplenishIonMsk(j)= UBreplenishIonMskTMP(j)
				end if

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
					LBoutfluxENAMsk(j)= LBoutfluxENAMskTMP(j)
					UBoutfluxENAMsk(j)= UBoutfluxENAMskTMP(j)
				end if

        AEAmagN(j)= AEAmagNTMP(j)
				AGmagN(j)= AGmagNTMP(j)
        AEPmagN(j)= AEPmagNTMP(j)
        Qindk1(j)= Qindk1TMP(j)
        Qindk0(j)= Qindk0TMP(j)
        Vparindk1(j)= Vparindk1TMP(j)

				if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

					Vperp1indk1(j)= Vperp1indk1TMP(j)
					Vperp1indk0(j)= Vperp1indk0TMP(j)
					Vperp2indk1(j)= Vperp2indk1TMP(j)
					Vperp2indk0(j)= Vperp2indk0TMP(j)

				else

					Vperpindk1(j)= Vperpindk1TMP(j)
					Vperpindk0(j)= Vperpindk0TMP(j)

				end if

        Vparindk0(j)= Vparindk0TMP(j)
        Vpindk1(j)= Vpindk1TMP(j)
        Vqindk1(j)= Vqindk1TMP(j)
        Vphiindk1(j)= Vphiindk1TMP(j)
        VxN(j)= VxNTMP(j)
        VyN(j)= VyNTMP(j)
        VzN(j)= VzNTMP(j)
        xN(j)= xNTMP(j)
        yN(j)= yNTMP(j)
        zN(j)= zNTMP(j)
        Vperp1N(j)= Vperp1NTMP(j)
        Vperp2N(j)= Vperp2NTMP(j)
        VperpN(j)= VperpNTMP(j)
        pk4(j)= pk4TMP(j)
        phik4(j)= phik4TMP(j)
        qk4(j)= qk4TMP(j)
        rk4(j)= rk4TMP(j)
        thetak4(j)= thetak4TMP(j)
				ellk4(j)= ellk4TMP(j)
        ENAflag(j)= ENAflagTMP(j)
        ENAflagN0ind(j)= ENAflagN0indTMP(j)
        ENAflagN1ind(j)= ENAflagN1indTMP(j)
        AEAmag(j)= AEAmagTMP(j)
				AGmag(j)= AGmagTMP(j)
        AEPmag(j)= AEPmagTMP(j)
        x(j)= xTMP(j)
        y(j)= yTMP(j)
        z(j)= zTMP(j)
        Vperp1(j)= Vperp1TMP(j)
        Vperp2(j)= Vperp2TMP(j)
        Vperp(j)= VperpTMP(j)
        Vpar(j)= VparTMP(j)
        Vx(j)= VxTMP(j)
        Vy(j)= VyTMP(j)
        Vz(j)= VzTMP(j)
        Vphi(j)= VphiTMP(j)
        Vp(j)= VpTMP(j)
        pdriftion(j)= pdriftionTMP(j)
        qdriftion(j)= qdriftionTMP(j)
        phidriftion(j)= phidriftionTMP(j)
        rdriftENA(j)= rdriftENATMP(j)
        thetadriftENA(j)= thetadriftENATMP(j)
        phidriftENA(j)= phidriftENATMP(j)
      end if
			if (dNsTK(1) < 0) then
				LBoutfluxIonMsk(j)= LBoutfluxIonMskTMP(j)
				UBoutfluxIonMsk(j)= UBoutfluxIonMskTMP(j)

				if (SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 1) then
					LBreplenishIonMsk(j)= LBreplenishIonMskTMP(j)
				end if
				if (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 1) then
					UBreplenishIonMsk(j)= UBreplenishIonMskTMP(j)
				end if

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
					LBoutfluxENAMsk(j)= LBoutfluxENAMskTMP(j)
					UBoutfluxENAMsk(j)= UBoutfluxENAMskTMP(j)
				end if

        AEAmagN(j)= AEAmagNTMP(j)
				AGmagN(j)= AGmagNTMP(j)
        AEPmagN(j)= AEPmagNTMP(j)
        Qindk1(j)= Qindk1TMP(j)
        Qindk0(j)= Qindk0TMP(j)
        Vparindk1(j)= Vparindk1TMP(j)

				if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

					Vperp1indk1(j)= Vperp1indk1TMP(j)
					Vperp1indk0(j)= Vperp1indk0TMP(j)
					Vperp2indk1(j)= Vperp2indk1TMP(j)
					Vperp2indk0(j)= Vperp2indk0TMP(j)

				else

					Vperpindk1(j)= Vperpindk1TMP(j)
					Vperpindk0(j)= Vperpindk0TMP(j)

				end if

        Vparindk0(j)= Vparindk0TMP(j)
        Vpindk1(j)= Vpindk1TMP(j)
        Vqindk1(j)= Vqindk1TMP(j)
        Vphiindk1(j)= Vphiindk1TMP(j)
        VxN(j)= VxNTMP(j)
        VyN(j)= VyNTMP(j)
        VzN(j)= VzNTMP(j)
        xN(j)= xNTMP(j)
        yN(j)= yNTMP(j)
        zN(j)= zNTMP(j)
        Vperp1N(j)= Vperp1NTMP(j)
        Vperp2N(j)= Vperp2NTMP(j)
        VperpN(j)= VperpNTMP(j)
        pk4(j)= pk4TMP(j)
        phik4(j)= phik4TMP(j)
        qk4(j)= qk4TMP(j)
        rk4(j)= rk4TMP(j)
        thetak4(j)= thetak4TMP(j)
				ellk4(j)= ellk4TMP(j)
        ENAflag(j)= ENAflagTMP(j)
        ENAflagN0ind(j)= ENAflagN0indTMP(j)
        ENAflagN1ind(j)= ENAflagN1indTMP(j)
        AEAmag(j)= AEAmagTMP(j)
				AGmag(j)= AGmagTMP(j)
        AEPmag(j)= AEPmagTMP(j)
        x(j)= xTMP(j)
        y(j)= yTMP(j)
        z(j)= zTMP(j)
        Vperp1(j)= Vperp1TMP(j)
        Vperp2(j)= Vperp2TMP(j)
        Vperp(j)= VperpTMP(j)
        Vpar(j)= VparTMP(j)
        Vx(j)= VxTMP(j)
        Vy(j)= VyTMP(j)
        Vz(j)= VzTMP(j)
        Vphi(j)= VphiTMP(j)
        Vp(j)= VpTMP(j)
        pdriftion(j)= pdriftionTMP(j)
        qdriftion(j)= qdriftionTMP(j)
        phidriftion(j)= phidriftionTMP(j)
        rdriftENA(j)= rdriftENATMP(j)
        thetadriftENA(j)= thetadriftENATMP(j)
        phidriftENA(j)= phidriftENATMP(j)
      end if
    end do

    ! ----------------------------------------------------

    ! DE-ALLOCATE TEMPORARY KINETIC ARRAYS:

    ! ----------------------------------------------------

    deallocate(LBoutfluxIonMskTMP, UBoutfluxIonMskTMP, AEAmagNTMP, AGmagNTMP, AEPmagNTMP)

		if (SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 1) then
			deallocate(LBreplenishIonMskTMP)
		end if
		if (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 1) then
			deallocate(UBreplenishIonMskTMP)
		end if

		if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
			deallocate(LBoutfluxENAMskTMP, UBoutfluxENAMskTMP)
		end if

		deallocate(Qindk1TMP, Qindk0TMP, Vparindk1TMP, &
			Vparindk0TMP, Vpindk1TMP, Vqindk1TMP, Vphiindk1TMP)

		if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

			deallocate(Vperp1indk1TMP, Vperp1indk0TMP)
			deallocate(Vperp2indk1TMP, Vperp2indk0TMP)

		else

			deallocate(Vperpindk1TMP, Vperpindk0TMP)

		end if

		deallocate(VxNTMP, VyNTMP, VzNTMP, xNTMP, yNTMP, zNTMP, Vperp1NTMP, Vperp2NTMP, VperpNTMP)
		deallocate(pk4TMP, phik4TMP, qk4TMP, rk4TMP, thetak4TMP, ellk4TMP)
		deallocate(ENAflagTMP, ENAflagN0indTMP, ENAflagN1indTMP)
    deallocate(AEAmagTMP, AGmagTMP, AEPmagTMP, xTMP, yTMP, zTMP, Vperp1TMP, Vperp2TMP, VperpTMP, VparTMP, &
      VxTMP, VyTMP, VzTMP, VphiTMP, VpTMP)
    deallocate(pdriftionTMP, qdriftionTMP, phidriftionTMP, rdriftENATMP, &
      thetadriftENATMP, phidriftENATMP)

    ! ----------------------------------------------------

    ! DIAGNOSTICS FOR CONSISTENT ARRAY SIZES:

		if (size(LBoutfluxIonMsk) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' LBoutfluxIonMsk SIZE= ', &
        size(LBoutfluxIonMsk), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

		if (size(UBoutfluxIonMsk) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' UBoutfluxIonMsk SIZE= ', &
        size(UBoutfluxIonMsk), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

		if (SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 1) then
			if (size(LBreplenishIonMsk) /= NsTK(1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' LBreplenishIonMsk SIZE= ', &
					size(LBreplenishIonMsk), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
					' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
					nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
			end if
		end if
		if (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 1) then
			if (size(UBreplenishIonMsk) /= NsTK(1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' UBreplenishIonMsk SIZE= ', &
					size(UBreplenishIonMsk), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
					' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
					nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
			end if
		end if

		if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
			if (size(LBoutfluxENAMsk) /= NsTK(1)) then
	      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' LBoutfluxENAMsk SIZE= ', &
	        size(LBoutfluxENAMsk), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
	        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
	        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
	    end if

			if (size(UBoutfluxENAMsk) /= NsTK(1)) then
	      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' UBoutfluxENAMsk SIZE= ', &
	        size(UBoutfluxENAMsk), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
	        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
	        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
	    end if
		end if

    if (size(AEAmagN) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAmagN SIZE= ', &
        size(AEAmagN), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

		if (size(AGmagN) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGmagN SIZE= ', &
        size(AGmagN), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(AEPmagN) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPmagN SIZE= ', &
        size(AEPmagN), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(Qindk1) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Qindk1 SIZE= ', &
        size(Qindk1), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(Qindk0) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Qindk0 SIZE= ', &
        size(Qindk0), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(Vparindk0) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vparindk0 SIZE= ', &
        size(Vparindk0), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

		if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

			if (size(Vperp1indk0) /= NsTK(1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp1indk0 SIZE= ', &
					size(Vperp1indk0), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
					' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
					nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
			end if

			if (size(Vperp1indk0) /= NsTK(1)) then
	      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp1indk0 SIZE= ', &
	        size(Vperp1indk0), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
	        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
	        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
	    end if

			if (size(Vperp2indk1) /= NsTK(1)) then
	      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp2indk1 SIZE= ', &
	        size(Vperp2indk1), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
	        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
	        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
	    end if

			if (size(Vperp2indk1) /= NsTK(1)) then
	      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp2indk1 SIZE= ', &
	        size(Vperp2indk1), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
	        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
	        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
	    end if

		else

			if (size(Vperpindk0) /= NsTK(1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperpindk0 SIZE= ', &
					size(Vperpindk0), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
					' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
					nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
			end if

			if (size(Vperpindk1) /= NsTK(1)) then
	      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperpindk1 SIZE= ', &
	        size(Vperpindk1), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
	        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
	        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
	    end if

		end if

    if (size(Vparindk1) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vparindk1 SIZE= ', &
        size(Vparindk1), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(Vpindk1) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vpindk1 SIZE= ', &
        size(Vpindk1), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(Vqindk1) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vqindk1 SIZE= ', &
        size(Vqindk1), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(Vphiindk1) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vphiindk1 SIZE= ', &
        size(Vphiindk1), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(VxN) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VxN SIZE= ', &
        size(VxN), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(VyN) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VyN SIZE= ', &
        size(VyN), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(VzN) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VzN SIZE= ', &
        size(VzN), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(xN) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' xN SIZE= ', &
        size(xN), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(yN) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' yN SIZE= ', &
        size(yN), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(zN) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' zN SIZE= ', &
        size(zN), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(Vperp1N) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp1N SIZE= ', &
        size(Vperp1N), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(Vperp2N) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp2N SIZE= ', &
        size(Vperp2N), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(VperpN) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VperpN SIZE= ', &
        size(VperpN), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(pk4) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' pk4 SIZE= ', &
        size(pk4), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(qk4) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qk4 SIZE= ', &
        size(qk4), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(phik4) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' phik4 SIZE= ', &
        size(phik4), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(rk4) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rk4 SIZE= ', &
        size(rk4), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(thetak4) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' thetak4 SIZE= ', &
        size(thetak4), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

		if (size(ellk4) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ellk4 SIZE= ', &
        size(ellk4), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(ENAflag) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ENAflag SIZE= ', &
        size(ENAflag), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(ENAflagN0ind) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ENAflagN0ind SIZE= ', &
        size(ENAflagN0ind), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(ENAflagN1ind) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ENAflagN1ind SIZE= ', &
        size(ENAflagN1ind), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(AEAmag) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAmag SIZE= ', &
        size(AEAmag), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

		if (size(AGmag) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGmag SIZE= ', &
        size(AGmag), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(AEPmag) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPmag SIZE= ', &
        size(AEPmag), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(x) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' x SIZE= ', &
        size(x), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(y) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' y SIZE= ', &
        size(y), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(z) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' z SIZE= ', &
        size(z), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(Vperp1) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp1 SIZE= ', &
        size(Vperp1), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(Vperp2) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp2 SIZE= ', &
        size(Vperp2), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(Vperp) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp SIZE= ', &
        size(Vperp), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(Vx) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vx SIZE= ', &
        size(Vx), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(Vy) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vy SIZE= ', &
        size(Vy), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(Vz) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vz SIZE= ', &
        size(Vz), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(Vphi) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vphi SIZE= ', &
        size(Vphi), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(Vp) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vp SIZE= ', &
        size(Vp), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(pdriftion) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' pdriftion SIZE= ', &
        size(pdriftion), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(qdriftion) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qdriftion SIZE= ', &
        size(qdriftion), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(phidriftion) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' phidriftion SIZE= ', &
        size(phidriftion), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(rdriftENA) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rdriftENA SIZE= ', &
        size(rdriftENA), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(thetadriftENA) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' thetadriftENA SIZE= ', &
        size(thetadriftENA), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

    if (size(phidriftENA) /= NsTK(1)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' phidriftENA SIZE= ', &
        size(phidriftENA), ' DOES NOT EQUAL PRIOR TOTAL PARTICLE NUMBER= ', NsTK(1)- dNsTK(1), &
        ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND MASTER TIME-STEP= ', &
        nn, ' IN DATA TYPE RE-ALLOCATION A SUBROUTINE' // achar(27) // '[0m.'
    end if

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine DataTypeReAllocASub

end module DataTypeReAllocA
