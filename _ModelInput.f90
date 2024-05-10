module ModelInput

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	MODEL INPUT:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

implicit none

! ----------------------------------------------------

integer, parameter :: dp = kind(1.d0) ! Define maximum compiler specific precision
complex(kind= dp) :: i ! Square-root of -1

! ----------------- I/O PATHS -----------------

! Set (=1) for spin-up simulations (output initial conditions, i.e. spin up simulation (=1))
integer(kind= dp), parameter :: SPINUPflag= 0

! Define I/O paths:
character(*), parameter :: Densitydatadir= '/Users/robertalbarran/Desktop/KAOS_M1/KAOSDensityInput/'

!if (SPINUPflag == 1) then
  !character(*), parameter :: dataexportdir= '/Users/robertalbarran/Desktop/KAOS_M1/KAOSDataSpinUp/'
!end if
!if (SPINUPflag == 0) then
  character(*), parameter :: dataexportdir= '/Users/robertalbarran/Desktop/KAOS_M1/KAOSDataOutput/'
!end if

! ----------------- PHYSICAL CONSTANTS -----------------

real(kind= dp), parameter :: pi= 4.d0*datan(1.d0) ! pi
real(kind= dp), parameter :: kB= 138d-25 ! Boltzmann constant [kg m^2 s^-2 K^-1]
real(kind= dp), parameter :: epsilon0= 8.854d-12 ! Permittivity of Free Space [m^-3 kg^-1 s^4 A^2]
real(kind= dp), parameter :: mu0= 1.257d-6 ! Permeability of Free Space [m kg s^-2 A^-2]
real(kind= dp), parameter :: m= 8.06d15 ! Earth magnetic moment [T m^3]
real(kind= dp), parameter :: GG= 667d-13 ! Universal gravitational constant [N m^2/kg^2]
real(kind= dp), parameter :: ME= 598d22 ! Earth mass [kg]
real(kind= dp), parameter :: RE= 638d4 ! Earth radius [m]
real(kind= dp), parameter :: melec= 9.109d-31 ! Electron mass [kg]

! ----------------- TIME PARAMETERS -----------------

!HERE00
real(kind= dp), parameter :: A= 0d0 ! Start time [s]
real(kind= dp), parameter :: B= 18000d0 ! 25200d0 ! 14400d0 ! (25200 for W0F, 18000d0 for W0 sims, 14400d0 for VV0, VV0F sims) End time [s]
integer(kind= dp), parameter :: Nt= 3d4 ! 2d4 ! 0.2d4 ! Total number of time-steps

! ----------------- ION SIMULATION FLAGS -----------------

! Set lower limit to particle counts per phase space bin (see ion noise limit in Ion limits) (leave =0)
integer(kind= dp), parameter :: IONNOISEflag= 0

! Set ion moment export (leave =1)
integer(kind= dp), parameter :: FLUIDIONEXPORTflag= 1

! Set symmetric +/- parallel velocity space grid for ions and ENAs (leave =1)
integer(kind= dp), parameter :: SYMVPARflag= 1

! Set 2D perpendicular velocity space grid for ions and ENAs (leave =1)
integer(kind= dp), parameter :: ION2VPERPflag= 1

! ----------------- BOUNDARY CONDITIONS FLAGS -----------------

! Set lower boundary ion injection on statistical time-steps (leave =1)
integer(kind= dp), parameter :: LBCONDITIONflag= 1

! Replenish lower boundary ion density on statistical time-steps (leave =0)
integer(kind= dp), parameter :: LBREPLENISHflag= 0

! Set upper boundary ion injection on statistical time-steps (leave =0)
integer(kind= dp), parameter :: UBCONDITIONflag= 0

! Replenish uppwer boundary ion density on statistical time-steps (leave =0)
integer(kind= dp), parameter :: UBREPLENISHflag= 0

! Set intial density to altitude profile (=1) or only on lower boundary (=0) (leave =1)
integer(kind= dp), parameter :: DENSITYPROFILEflag= 1

! Set initial velocities of injected ions equal to zero (leave =0)
integer(kind= dp), parameter :: STATICINJECTIONflag= 0

! ----------------- ION FORCE FLAGS -----------------

! Set wave-particle interactions
integer(kind= dp), parameter :: ICRflag= 0

! Set coherent wave-particle interactions
integer(kind= dp), parameter :: ICRCOHERENCEflag= 0

! Set mirror-force
integer(kind= dp), parameter :: MIRRORflag= 0

! Set gravitational force
integer(kind= dp), parameter :: GRAVflag= 1

! Set ambipolar electric field
integer(kind= dp), parameter :: EAMBflag= 1

! Set self-consistent ambipolar electric field (leave =1 for EAMBflag=1)
integer(kind= dp), parameter :: EAMBSELFCONSISTflag= 1

! Set (=1) for both positive and negative sign of ambipolar electric field (leave =1 for EAMBflag=1)
integer(kind= dp), parameter :: EAMBSIGNflag= 1

! Set inertial component of ambipolar electric field (Note: obsolete, inertial effects are negligible) (leave =0)
integer(kind= dp), parameter :: EAINERTIALflag= 0

! Set electron pressure component of ambipolar electric field (leave =1 for EAMBflag=1)
integer(kind= dp), parameter :: EAPRESSUREflag= 1

! Set parallel electric field
integer(kind= dp), parameter :: EPARflag= 0

! Set flux-tube converction:
integer(kind= dp), parameter :: CONVECTIONflag= 0

! Set dynamic configuration-space and statistical time-step:
integer(kind= dp), parameter :: DYNAMICGRIDflag= 0

! ----------------- ION MOMENT FLAGS -----------------

! Set sliding average filter along flux-tube to moments (for smooth ambipolar electric field, etc) (leave =1)
integer(kind= dp), parameter :: MOMENTFILTERflag= 1

! Set ion distribution function computation (leave =1)
integer(kind= dp), parameter :: PHASEIONDISTRIBflag= 1

! Set ion density computation (leave =1)
integer(kind= dp), parameter :: PHASEDENSITYIONMOMENTflag= 1

! Set ion perp velocity computation (leave =1)
integer(kind= dp), parameter :: PHASEVELPERPIONMOMENTflag= 1

! Set ion par velocity computation (leave =1)
integer(kind= dp), parameter :: PHASEVELPARIONMOMENTflag= 1

! Set ion total energy computation (leave =1)
integer(kind= dp), parameter :: PHASEENERGYIONMOMENTflag= 1

! Set ion perp energy computation (leave =1)
integer(kind= dp), parameter :: PHASEENERGYPERPIONMOMENTflag= 1

! Set ion par energy computation (leave =1)
integer(kind= dp), parameter :: PHASEENERGYPARIONMOMENTflag= 1

! Set ion perp energy centering by perp velocity (leave =0)
integer(kind= dp), parameter :: PHASEENERGYPERPIONMOMENTCENTERflag= 0

! Set ion par energy centering by par velocity (leave =1)
integer(kind= dp), parameter :: PHASEENERGYPARIONMOMENTCENTERflag= 1

! Compute reference ion moments for validation (leave =0)
integer(kind= dp), parameter :: FLUIDIONREFflag= 0

! ----------------- ENA FLAGS -----------------

! Set ENA production
integer(kind= dp), parameter :: QEXCHANGEflag= 0

! Set ENA velocity space count filter
integer(kind= dp), parameter :: ENANOISEflag= 0

! ----------------- ION AND ENA PHASE-SPACE GRID PARAMETERS -----------------

! Define number of particle species and flux tubes
integer(kind= dp), parameter :: Stot= 1d0 ! Number of particle species
integer(kind= dp), parameter :: Nf= 1d0 ! Number of flux tubes per species
real(kind= dp), parameter :: LshellIC= 10d0 ! Initial L-shell
real(kind= dp), parameter :: phiLshellIC= pi/2d0 ! Initial L-shell longitude
real(kind= dp), parameter :: qGAIC= 0.6d0 ! Set lower boundary q value (< for South Magnetic Hemisphere and > for North Magnetic Hemisphere)
real(kind= dp), parameter :: qGBIC= 0.1d0 ! Set upper boundary q value (< for South Magnetic Hemisphere and > for North Magnetic Hemisphere)

integer(kind= dp) :: SMagHemFlag
real(kind= dp), parameter :: mO= (16d0)*(1.66054d-27) ! O+ mass [kg]
real(kind= dp), parameter :: qO= 1.602d-19 ! O+ Charge [C]

integer(kind= dp), parameter :: dq= 1d0 ! dq length of config cells
integer(kind= dp), parameter :: ddVperp1= 1d0 ! dVperp1 length of config cells
integer(kind= dp), parameter :: ddVperp2= 1d0 ! dVperp2 length of config cells
integer(kind= dp), parameter :: ddVperp= 1d0 ! dVperp length of config cells
integer(kind= dp), parameter :: ddVpar= 1d0 ! dVpar length of config cells
integer(kind= dp), parameter :: ddVp= 1d0 ! dVp length of config cells
integer(kind= dp), parameter :: ddVq= 1d0 ! dVq length of config cells
integer(kind= dp), parameter :: ddVphi= 1d0 ! dVphi length of config cells

integer(kind= dp) :: NqGpF ! Preliminary number of Q Grid Cells (including lower and upper boundary ghost cells i.e. +3)

real(kind= dp), parameter :: NVparGpF= 30d0 ! Preliminary number of Vpar Grid Cells (even (div by 2 odd) for +/- log10 Vpar grid) (+ 3)
real(kind= dp), parameter :: NVperpGpF= 30d0 ! Preliminary number of Vperp Grid Cells (even (div by 2 odd) for +/- log10 Vpar grid) (+ 3)
real(kind= dp), parameter :: VperpsigmaFacIC= 8d0 ! (=4 for thermal) Sigma factor with linear grid to resolve thermal core of MB distrib.
real(kind= dp), parameter :: VparsigmaFacIC= 8d0 ! (=4 for thermal) Sigma factor with linear grid to resolve thermal core of MB distrib.

real(kind= dp), parameter :: TiIC= 5d3 ! Ion initialization temperature [K]
real(kind= dp), parameter :: TeIC= 1d4 ! Electron initialization temperature [K]
real(kind= dp), parameter :: TNeutIC= 848d0 ! O neutral temperature [K] from NRLMSISE-00

! ----------------- ION INITIALIZATION PARAMETERS -----------------
! Note: change lower boundary cell (NqLB) and upper boundary cell (NqUB) and macro-particle normalization constant (nsnormfac) such that no cells are empty

! Mind NqLB and NqUB values wrt moment filter sizes below (e.g. M0MAfilterPt, M1Perp1MAfilterPt, ...)
integer(kind= dp), parameter :: NqLB= 1d0 ! Leave =1 (change q ranges to adjust)
integer(kind= dp), parameter :: NqUB= 25d0

real(kind= dp), parameter :: nsnormfacIC= 5d18 !8d14 ! Macro-particle normalization factor (inversely proportional to number of particles)

real(kind= dp), parameter :: dNTe= 0d0 ! Additive increment of Te on statistical time-step [K]
real(kind= dp), parameter :: dNTeEND= 3d10 ! Additive increment of Te Cap on statistical time-step [K]

real(kind= dp), parameter :: ns0IC= 1d8 ! Initial ion reference density at LB ghost cell (Qindns0) [m-3]

! ----------------- NEUTRAL OXYGEN INITIALIZATION PARAMETERS -----------------

! Use from NRLMSISE-00 data at PFISR at VISIONS-1 tof for neutral density profile
real(kind= dp), parameter :: zns0NeutIC= RE+ 719d3 ! Reference O neutral r value [m]
real(kind= dp), parameter :: ns0NeutIC= 1.423d11 ! Reference O neutral number density at zns0Neut [m^-3]

! ----------------- ICR HEATING PARAMETERS -----------------

real(kind= dp), parameter :: EtaLHp= 0.125d0 ! Fraction of LH polarized ELF wave power (Chang, '86)
real(kind= dp), parameter :: XiPerp1p= 0.5d0 ! Fraction of LH Polarized ELF wave power along ePerp1
real(kind= dp), parameter :: XiPerp2p= 0.5d0 ! Fraction of LH Polarized ELF wave power along ePerp2

real(kind= dp), parameter :: lambdaPerpp= 0d0 ! VLF wavelength at local gyrofrequency [m] (set to zero for infinite wavelength limit)

real(kind= dp), parameter :: f0p= 6.5d0 ! 72.5333d0 ! Reference heating gyrofrequency [Hz]
real(kind= dp), parameter :: S0p= 5d-5 ! 1.9052d-10 ! Reference Wave Spectral Energy Density in [(V^2/m^2)/Hz]

real(kind= dp), parameter :: ChiPerp1p= 2.1d0 ! ELF wave spectral index along vperp1
real(kind= dp), parameter :: ChiPerp2p= 2.1d0 ! ELF wave spectral index along vperp2

! ----------------- PARALLEL POTENTIAL PARAMETERS -----------------

integer(kind= dp), parameter :: Eparlim= 1d0 ! Number of computational timesteps before Epar turns on
real(kind= dp), parameter :: PhiPar0p= 5d-6 ! [V] Set to reference Eparallel value (eg. 9.5d-6 V/m for Wu '02)
real(kind= dp), parameter :: dPhiPar0p= 1d0 ! [m] See comment with PhiPar0p above

! ----------------- ION AND ENA LIMITS -----------------

real(kind= dp), parameter :: IonNoiseLimitNph= 2500d0 ! 15d3 ! Minimum number of ion macroparticles per phase space bin
real(kind= dp), parameter :: ENANoiseLimitNph= 5d0 ! Minimum number of ENA macroparticles per phase space bin

integer(kind= dp), parameter :: M0MAfilterPt= 11d0 ! Moving-average filter point for density moment (odd number)
integer(kind= dp), parameter :: M1Perp1MAfilterPt= 11d0 ! Moving-average filter point for perp1 velocity moment (odd number)
integer(kind= dp), parameter :: M1Perp2MAfilterPt= 11d0 ! Moving-average filter point for perp2 velocity moment (odd number)
integer(kind= dp), parameter :: M1ParMAfilterPt= 11d0 ! Moving-average filter point for parallel velocity moment (odd number)
integer(kind= dp), parameter :: M2ParMAfilterPt= 11d0 ! Moving-average filter point for parallel energy moment (odd number)

end module ModelInput
