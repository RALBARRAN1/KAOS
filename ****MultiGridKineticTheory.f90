On each time-step set grid space for all time after initial time in kinetic solver for each de-magnetized ion.
Insert this part on ENA production process. 

Do this by using solitary variables with kinetic particle j measurement in n (for computational optimization purposes). Use time-step N_tt for these purposes, where kinetic distribution functions and macroscopic parameters are computed. Use phase-space transform modules in dipole, spherical, and rectangular coordinates, to create
a per-particle computation of metric factors, phase-space limits, and volumes for all Ntt/= 1, (where Ntt= 1 correspondes to the initial conditions). According to these time-dependent phase-space grid transforms, all particles flow. 

do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. & 
				(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn)))) then 
				
				