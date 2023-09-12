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

! Define these I/O paths with local CPU simulations:
! Note: Set (1) Grid directory (created in MATLAB), (2) Export data directory, (3) Spin-up initialization data directory

character(*), parameter :: dataimportdir= '/Users/robertalbarran/Desktop/FULLDESKTOP/KMIOGridW0F2/'
character(*), parameter :: dataexportdir= '/Users/robertalbarran/Desktop/FULLDESKTOP/KMIOData/'
character(*), parameter :: Densitydatadir= '/Users/robertalbarran/Desktop/FULLDESKTOP/KMIODensityInput/'

! Define these I/O paths with M1 Ultra CPU simulations:
!character(*), parameter :: dataimportdir= '/Users/robertalbarran/Desktop/KMIOGridW0F1/'
!character(*), parameter :: dataexportdir= '/Users/robertalbarran/Desktop/KMIODataW0F1ncscOUTPUTa0/'
!character(*), parameter :: Densitydatadir= '/Users/robertalbarran/Desktop/KMIODensityW0F1ncscOUTPUTa0/'

! Define these I/O paths with VEGA CPU simulations:
!character(*), parameter :: dataimportdir= '/home2/albarrar/KMIOGridVV0Ft/'
!character(*), parameter :: dataexportdir= '/home2/albarrar/KMIODataVV0FticrV0epar6scINPUTa4p/'
!character(*), parameter :: Densitydatadir= '/home2/albarrar/KMIODensityVV0FtncscOUTPUTa0/'

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
real(kind= dp) mNeut ! O neutral mass [kg]
real(kind= dp) TNeut ! O neutral temperature [K] from NRLMSISE-00

! ----------------- TIME PARAMETERS -----------------

real(kind= dp), parameter :: A= 0d0 ! Start time [s]
real(kind= dp), parameter :: B= 18000d0 ! 25200d0 ! 14400d0 ! (25200 for W0F, 18000d0 for W0 sims, 14400d0 for VV0, VV0F sims) End time [s]
integer(kind= dp), parameter :: Nt= 3d4 ! 2d4 ! 0.2d4 ! Total number of time-steps

! ----------------- ION SIMULATION FLAGS -----------------

! Set lower limit to particle counts per phase space bin (see ion noise limit in Ion limits) (leave =0)
integer(kind= dp), parameter :: IONNOISEflag= 0

! Set ion moment export (leave =1)
integer(kind= dp), parameter :: FLUIDIONEXPORTflag= 1

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

! Set (=1) for spin-up simulations (output initial conditions)
integer(kind= dp), parameter :: DENSITYOUTPUTflag= 1

! Set (=1) for non spin-up simulations (input initial conditions from spin-up)
integer(kind= dp), parameter :: DENSITYINPUTflag= 0

! ----------------- ION FORCE FLAGS -----------------

! Set wave-particle interactions
integer(kind= dp), parameter :: ICRflag= 0

! Set mirro-force
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

! ----------------- ION MOMENT FLAGS -----------------

! Set sliding average filter along flux-tube to moments (for smooth ambipolar electric field, etc) (leave =1) 
integer(kind= dp), parameter :: MOMENTFILTERflag= 1 

! Set 2D perpendicular velocity space grid (leave =1)
integer(kind= dp), parameter :: ION2VPERPflag= 1

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

! ----------------- ION INITIALIZATION PARAMETERS -----------------
! Note: change lower boundary cell (NqICA) and upper boundary cell (NqICB) and macro-particle normalization constant (nsnormfac) such that no cells are empty

integer(kind= dp), parameter :: NqICA= 1d0
integer(kind= dp), parameter :: NqICB= 25d0
real(kind= dp), parameter :: nsnormfac= 9.5d13 ! (W0F2)

real(kind= dp), parameter :: dNTe= 0d0 ! Additive increment of Te on statistical time-step [K]
real(kind= dp), parameter :: dNTeEND= 3d10 ! Additive increment of Te Cap on statistical time-step [K]
real(kind= dp), parameter :: dns0= 1d0 ! Multiplicative factor of reference ion density

real(kind= dp), parameter :: zns0= RE+ 389d3 ! RE+ 370.78d3 ! Initial ion density profile reference altitude [km] 
real(kind= dp), parameter :: ns0= 1d7*dns0 ! (W0F2) ! Initial ion density profile reference density at zns0 [m-3]

! ----------------- NEUTRAL OXYGEN INITIALIZATION PARAMETERS -----------------

! Use from NRLMSISE-00 data at PFISR at VISIONS-1 tof for neutral density profile
real(kind= dp), parameter :: zns0Neut= RE+ 719d3 ! Reference O neutral r value [m]
real(kind= dp), parameter :: ns0Neut= 1.423d11 ! Reference O neutral number density at zns0Neut [m^-3]

! ----------------- ICR HEATING PARAMETERS -----------------

real(kind= dp), parameter :: EtaLHp= 0.125d0 ! Fraction of LH polarized ELF wave power (Chang, '86)
real(kind= dp), parameter :: XiPerp1p= 0.5d0 ! Fraction of LH Polarized ELF wave power along ePerp1
real(kind= dp), parameter :: XiPerp2p= 0.5d0 ! Fraction of LH Polarized ELF wave power along ePerp2

real(kind= dp), parameter :: lambdaPerpp= 0d0 ! VLF wavelength at local gyrofrequency [m] (set to zero for long wavelength limit)

real(kind= dp), parameter :: f0p= 6.5d0 ! 72.5333d0 ! Reference heating gyrofrequency [Hz]
real(kind= dp), parameter :: S0p= 5d-7 ! 1.9052d-10 ! Reference Wave Spectral Energy Density in [(V^2/m^2)/Hz]

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
