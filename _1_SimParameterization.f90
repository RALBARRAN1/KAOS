module SimParameterization

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	1- SIMULATION PARAMETERIZATION:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use ConfigGridGenerator
use ConfigGridGeneratorB
use VelGridGenerator

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! SET ALL PHASE-SPACE GRID, INITIALIZATION, AND SIMULATION PARAMETERS PER PARTICLE SPECIES,
! FLUX TUBE, AND PHASE-SPACE GRID CELL:

	subroutine SimParameterizationSub

		! ----------------------------------------------------

		! ALLOCATE PARTICLE SPECIES DERIVED DATA TYPE:

		i= (0, 1) ! Define sqrt(-1)
		allocate(SpecieT(Stot))

		! ----------------------------------------------------

		! ALLOCATE FLUX TUBE DERIVED DATA TYPE NESTED IN PARTICLE SPECIES TYPE:

		do s= 1, Stot, 1

			SpecieT(s)%NfT(1)= Nf
			allocate(SpecieT(s)%FluxTubeT(SpecieT(s)%NfT(1)))

		end do

		! ----------------------------------------------------

		! SET SIMULATION FLAGS:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				! ----------------------------------------------------

				! ION POPULATION FLAGS:

				! Create nested derived data types
 				SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT= FLUIDIONEXPORTflag
				SpecieT(s)%FluxTubeT(f)%SYMVPARflagT= SYMVPARflag
				SpecieT(s)%FluxTubeT(f)%LBCONDITIONflagT= LBCONDITIONflag
				SpecieT(s)%FluxTubeT(f)%UBCONDITIONflagT= UBCONDITIONflag
				SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT= LBREPLENISHflag
				SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT= UBREPLENISHflag
				SpecieT(s)%FluxTubeT(f)%STATICINJECTIONflagT= STATICINJECTIONflag
				SpecieT(s)%FluxTubeT(f)%SPINUPflagT= SPINUPflag

				SpecieT(s)%FluxTubeT(f)%IONNOISEflagT= IONNOISEflag
 				SpecieT(s)%FluxTubeT(f)%ICRflagT= ICRflag
				SpecieT(s)%FluxTubeT(f)%ICRCOHERENCEflagT= ICRCOHERENCEflag
				SpecieT(s)%FluxTubeT(f)%MIRRORflagT= MIRRORflag

				SpecieT(s)%FluxTubeT(f)%GRAVflagT= GRAVflag
				SpecieT(s)%FluxTubeT(f)%EAMBflagT= EAMBflag
				SpecieT(s)%FluxTubeT(f)%EAMBSELFCONSISTflagT= EAMBSELFCONSISTflag
				SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT= EAMBSIGNflag
				SpecieT(s)%FluxTubeT(f)%EAINERTIALflagT= EAINERTIALflag
				SpecieT(s)%FluxTubeT(f)%EAPRESSUREflagT= EAPRESSUREflag

				SpecieT(s)%FluxTubeT(f)%EPARflagT= EPARflag
				SpecieT(s)%FluxTubeT(f)%CONVECTIONflagT= CONVECTIONflag
				SpecieT(s)%FluxTubeT(f)%DYNAMICGRIDflagT= DYNAMICGRIDflag

				SpecieT(s)%FluxTubeT(s)%MOMENTFILTERflagT= MOMENTFILTERflag
				SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT= ION2VPERPflag
				SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT= PHASEIONDISTRIBflag
				SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT= PHASEDENSITYIONMOMENTflag
				SpecieT(s)%FluxTubeT(f)%PHASEVELPERPIONMOMENTflagT= PHASEVELPERPIONMOMENTflag
				SpecieT(s)%FluxTubeT(f)%PHASEVELPARIONMOMENTflagT= PHASEVELPARIONMOMENTflag
				SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT= PHASEENERGYIONMOMENTflag
				SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTflagT= &
					PHASEENERGYPERPIONMOMENTflag
				SpecieT(s)%FluxTubeT(f)%PHASEENERGYPARIONMOMENTflagT= &
					PHASEENERGYPARIONMOMENTflag
				SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTCENTERflagT= &
					PHASEENERGYPERPIONMOMENTCENTERflag
				SpecieT(s)%FluxTubeT(f)%PHASEENERGYPARIONMOMENTCENTERflagT= &
					PHASEENERGYPARIONMOMENTCENTERflag
				SpecieT(s)%FluxTubeT(f)%FLUIDIONREFflagT= FLUIDIONREFflag

				SpecieT(s)%FluxTubeT(f)%ENANOISEflagT= ENANOISEflag
				SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT= QEXCHANGEflag

				! ----------------------------------------------------

				! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

				if ((FLUIDIONEXPORTflag /= SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' FLUIDIONEXPORTflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SYMVPARflag /= SpecieT(s)%FluxTubeT(f)%SYMVPARflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%SYMVPARflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%SYMVPARflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' SYMVPARflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((LBCONDITIONflag /= SpecieT(s)%FluxTubeT(f)%LBCONDITIONflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%LBCONDITIONflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%LBCONDITIONflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' LBCONDITIONflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((UBCONDITIONflag /= SpecieT(s)%FluxTubeT(f)%UBCONDITIONflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%UBCONDITIONflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%UBCONDITIONflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' UBCONDITIONflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((LBREPLENISHflag /= SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' LBREPLENISHflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((UBREPLENISHflag /= SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' UBREPLENISHflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((STATICINJECTIONflag /= SpecieT(s)%FluxTubeT(f)%STATICINJECTIONflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%STATICINJECTIONflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%STATICINJECTIONflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' STATICINJECTIONflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SPINUPflag /= SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%SPINUPflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' SPINUPflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((IONNOISEflag /= SpecieT(s)%FluxTubeT(f)%IONNOISEflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%IONNOISEflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%IONNOISEflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' IONNOISEflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((ICRflag /= SpecieT(s)%FluxTubeT(f)%ICRflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%ICRflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%ICRflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ICRflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((ICRCOHERENCEflag /= SpecieT(s)%FluxTubeT(f)%ICRCOHERENCEflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%ICRCOHERENCEflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%ICRCOHERENCEflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ICRCOHERENCEflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((MIRRORflag /= SpecieT(s)%FluxTubeT(f)%MIRRORflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%MIRRORflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%MIRRORflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' MIRRORflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((GRAVflag /= SpecieT(s)%FluxTubeT(f)%GRAVflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%GRAVflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%GRAVflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' GRAVflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((EAMBflag /= SpecieT(s)%FluxTubeT(f)%EAMBflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%EAMBflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%EAMBflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' EAMBflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((EAMBSELFCONSISTflag /= SpecieT(s)%FluxTubeT(f)%EAMBSELFCONSISTflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%EAMBSELFCONSISTflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%EAMBSELFCONSISTflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' EAMBSELFCONSISTflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((EAMBSIGNflag /= SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' EAMBSIGNflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((EAINERTIALflag /= SpecieT(s)%FluxTubeT(f)%EAINERTIALflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%EAINERTIALflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%EAINERTIALflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' EAINERTIALflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((EAPRESSUREflag /= SpecieT(s)%FluxTubeT(f)%EAPRESSUREflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%EAPRESSUREflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%EAPRESSUREflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' EAPRESSUREflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((EPARflag /= SpecieT(s)%FluxTubeT(f)%EPARflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%EPARflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%EPARflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' EPARflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((CONVECTIONflag /= SpecieT(s)%FluxTubeT(f)%CONVECTIONflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%CONVECTIONflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%CONVECTIONflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' CONVECTIONflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((DYNAMICGRIDflag /= SpecieT(s)%FluxTubeT(f)%DYNAMICGRIDflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%DYNAMICGRIDflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%DYNAMICGRIDflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' DYNAMICGRIDflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((MOMENTFILTERflag /= SpecieT(s)%FluxTubeT(f)%MOMENTFILTERflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%MOMENTFILTERflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%MOMENTFILTERflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' MOMENTFILTERflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((ION2VPERPflag /= SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ION2VPERPflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((PHASEIONDISTRIBflag /= SpecieT(s)%FluxTubeT(f)% &
					PHASEIONDISTRIBflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' PHASEIONDISTRIBflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((PHASEDENSITYIONMOMENTflag /= SpecieT(s)%FluxTubeT(f)% &
					PHASEDENSITYIONMOMENTflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT(:)) /= 1) &
					.or. (isnan(real(SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' PHASEDENSITYIONMOMENTflagT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE', &
					' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((PHASEVELPERPIONMOMENTflag /= SpecieT(s)%FluxTubeT(f)% &
					PHASEVELPERPIONMOMENTflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%PHASEVELPERPIONMOMENTflagT(:)) /= 1) &
					.or. (isnan(real(SpecieT(s)%FluxTubeT(f)%PHASEVELPERPIONMOMENTflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' PHASEVELPERPIONMOMENTflagflagT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE', &
					' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((PHASEVELPARIONMOMENTflag /= SpecieT(s)%FluxTubeT(f)% &
					PHASEVELPARIONMOMENTflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%PHASEVELPARIONMOMENTflagT(:)) /= 1) &
					.or. (isnan(real(SpecieT(s)%FluxTubeT(f)%PHASEVELPARIONMOMENTflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' PHASEVELPARIONMOMENTflagT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE', &
					' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((PHASEENERGYIONMOMENTflag /= SpecieT(s)%FluxTubeT(f)% &
					PHASEENERGYIONMOMENTflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(:)) /= 1) &
					.or. (isnan(real(SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' PHASEENERGYIONMOMENTflagT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE', &
					' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((PHASEENERGYPERPIONMOMENTflag /= SpecieT(s)%FluxTubeT(f)% &
					PHASEENERGYPERPIONMOMENTflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTflagT(:)) /= 1) &
					.or. (isnan(real(SpecieT(s)%FluxTubeT(f)% &
					PHASEENERGYPERPIONMOMENTflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' PHASEENERGYPERPIONMOMENTflagT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE', &
					' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((PHASEENERGYPARIONMOMENTflag /= SpecieT(s)%FluxTubeT(f)% &
					PHASEENERGYPARIONMOMENTflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%PHASEENERGYPARIONMOMENTflagT(:)) /= 1) &
					.or. (isnan(real(SpecieT(s)%FluxTubeT(f)% &
					PHASEENERGYPARIONMOMENTflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' PHASEENERGYPARIONMOMENTflagT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE', &
					' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((PHASEENERGYPERPIONMOMENTCENTERflag /= SpecieT(s)%FluxTubeT(f)% &
					PHASEENERGYPERPIONMOMENTCENTERflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTCENTERflagT(:)) /= 1) &
					.or. (isnan(real(SpecieT(s)%FluxTubeT(f)% &
					PHASEENERGYPERPIONMOMENTCENTERflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' PHASEENERGYPERPIONMOMENTCENTERflagT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE', &
					' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((PHASEENERGYPARIONMOMENTCENTERflag /= SpecieT(s)%FluxTubeT(f)% &
					PHASEENERGYPARIONMOMENTCENTERflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%PHASEENERGYPARIONMOMENTCENTERflagT(:)) /= 1) &
					.or. (isnan(real(SpecieT(s)%FluxTubeT(f)% &
					PHASEENERGYPARIONMOMENTCENTERflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' PHASEENERGYPARIONMOMENTCENTERflagT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE', &
					' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((FLUIDIONREFflag /= SpecieT(s)%FluxTubeT(f)% &
					FLUIDIONREFflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%FLUIDIONREFflagT(:)) /= 1) &
					.or. (isnan(real(SpecieT(s)%FluxTubeT(f)%FLUIDIONREFflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' FLUIDIONREFflagT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE', &
					' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((ENANOISEflag /= SpecieT(s)%FluxTubeT(f)%ENANOISEflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%ENANOISEflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%ENANOISEflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ENANOISEflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((QEXCHANGEflag /= SpecieT(s)%FluxTubeT(f)% &
					QEXCHANGEflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(:)) /= 1) &
					.or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' QEXCHANGEflagT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE', &
					' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				! ----------------------------------------------------

				! DIAGNOSTIC FLAGS FOR MOMENT FLAG CONSISTENCY:

				if (((SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT(1) == 1)) .or. &
					((SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEVELPERPIONMOMENTflagT(1) == 1)) .or. &
					((SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEVELPARIONMOMENTflagT(1) == 1)) .or. &
					((SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(1) == 1)) .or. &
					((SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEENERGYPARIONMOMENTflagT(1) == 1)) .or. &
					((SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTflagT(1) == 1))) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT PHASEIONDISTRIBflagT AND ION MOMENT FLAGS ', &
						' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, &
						' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%LBCONDITIONflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT LBREPLENISHflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%UBCONDITIONflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT UBREPLENISHflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEVELPERPIONMOMENTflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT(1) == 0)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT PHASEVELPERPIONMOMENTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEVELPARIONMOMENTflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT(1) == 0)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT PHASEVELPARIONMOMENTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if (((SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT(1) == 0)) .or. &
					((SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEVELPERPIONMOMENTflagT(1) == 0)) .or. &
					((SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEVELPARIONMOMENTflagT(1) == 0))) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT PHASEENERGYIONMOMENTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(1) == 0)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT PHASEENERGYPERPIONMOMENTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEENERGYPARIONMOMENTflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(1) == 0)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT PHASEENERGYPARIONMOMENTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTCENTERflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT PHASEENERGYIONMOMENTflagT AND PHASEENERGYPERPIONMOMENTCENTERflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEENERGYPARIONMOMENTCENTERflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT PHASEENERGYIONMOMENTflagT AND PHASEENERGYPARIONMOMENTCENTERflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%FLUIDIONREFflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT PHASEIONDISTRIBflagT AND FLUIDIONREFflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT PHASEIONDISTRIBflagT AND FLUIDIONEXPORTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if (((SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT(1) == 0)) .or. &
					((SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEVELPARIONMOMENTflagT(1) == 0))) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT EAMBflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%EAINERTIALflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%EAPRESSUREflagT(1) == 0)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT EAINETIALflagT AND EAPRESSUREflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%EAINERTIALflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT EAINETIALflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%EAPRESSUREflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT EAPRESSUREflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%EAMBSELFCONSISTflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT EAMBSELFCONSISTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT EAMBSIGNflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%MOMENTFILTERflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT MOMENTFILTERflagT AND PHASEDENSITYIONMOMENTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEVELPERPIONMOMENTflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%MOMENTFILTERflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT MOMENTFILTERflagT AND PHASEVELPERPIONMOMENTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEVELPARIONMOMENTflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%MOMENTFILTERflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT MOMENTFILTERflagT AND PHASEVELPARIONMOMENTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEENERGYPARIONMOMENTflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%MOMENTFILTERflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT MOMENTFILTERflagT AND PHASEENERGYPARIONMOMENTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT PHASEENERGYIONMOMENTflagT AND QEXCHANGEflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%ENANOISEflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 0)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT ENANOISEflagT AND QEXCHANGEflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%DYNAMICGRIDflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%CONVECTIONflagT(1) == 0)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT DYNAMICGRIDflagT AND CONVECTIONflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1) == 0)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT SPINUPflagT AND FLUIDIONEXPORTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				! ----------------------------------------------------

			end do
		end do

		! ----------------------------------------------------

		! PARTICLE MASS, CHARGE, AND ATOMIC RADIUS PER PARTICLE SPECIES:

		do s= 1, Stot, 1

			SpecieT(s)%msT= mO
			SpecieT(s)%qsT= qO

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

			if ((size(SpecieT(s)%Qindns0GT(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%Qindns0GT(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Qindns0GT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((mO /= SpecieT(s)%msT(1)) .or. (size(SpecieT(s)%msT(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%msT(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' msT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((qO /= SpecieT(s)%qsT(1)) .or. (size(SpecieT(s)%qsT(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%qsT(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qsT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
			end if

			! ----------------------------------------------------

		end do

		! ----------------------------------------------------

		! SET SIMULATION DURATION AND COMPUTATIONAL TIME-STEP:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				h(1)= (B- A)/Nt ! Time-step size [s]

				SpecieT(s)%AT= A ! Create nested derived data types
				SpecieT(s)%BT= B
				SpecieT(s)%FluxTubeT(f)%NtT= Nt
				SpecieT(s)%hT= h(1)

				! ----------------------------------------------------

				! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

				if ((A /= SpecieT(s)%AT) .or. &
					(isnan(real(SpecieT(s)%AT)) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((B /= SpecieT(s)%BT) .or. &
					(isnan(real(SpecieT(s)%BT)) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' BT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((Nt /= SpecieT(s)%FluxTubeT(f)%NtT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%NtT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%NtT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NtT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((h(1) /= SpecieT(s)%hT) .or. &
					(isnan(real(SpecieT(s)%hT)) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' hT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				! ----------------------------------------------------

			end do
		end do

		! ----------------------------------------------------

		! DIAGNOSTIC FLAG FOR CORRECT CONFIGURATION-SPACE GRID WITH MAGNETIC HEMSIPHERE:

		if (abs(qGAIC) <= abs(qGBIC)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
				' INCONSISTENT qGAIC AND qGBIC WITH MAGNETIC HEMISPHERE', &
				' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! DEFINE CONFIGURATION-SPACE BOUNDARY GRID CELLS:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				if (SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1) == 1) then

					SpecieT(s)%FluxTubeT(f)%NqLBT(1)= NqLB
					SpecieT(s)%FluxTubeT(f)%NqUBT(1)= NqUB

				end if
				if (SpecieT(1)%FluxTubeT(1)%SPINUPflagT(1) == 0) then
					write(sstring, '(I5)') s
					write(fstring, '(I5)') f

					open (unit= 0, file= Densitydatadir &
						// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
						// '_' // 'NqLBTfort.bin')), &
						status= 'old', form= 'unformatted', access= 'stream')
					read(0) NqLBp
					close(0)

					open (unit= 0, file= Densitydatadir &
						// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
						// '_' // 'NqUBTfort.bin')), &
						status= 'old', form= 'unformatted', access= 'stream')
					read(0) NqUBp
					close(0)

					SpecieT(s)%FluxTubeT(f)%NqLBT(1)= NqLBp
					SpecieT(s)%FluxTubeT(f)%NqUBT(1)= NqUBp

				end if
			end do
		end do

		! ----------------------------------------------------

		! ALLOCATE CONFIGURATION-SPACE DERIVED DATA TYPES NESTED IN FLUX TUBE TYPES NESTED IN
		! PARTICLE SPECIES TYPE:

		! Allocate QCellT(Qind) derived data type nested in FluxTubeT(f)
		! nested in SpecieT(s)
		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				allocate(SpecieT(s)%FluxTubeT(f)%QCellT(((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1)))
			end do
		end do

		! ----------------------------------------------------

		! IMPORT EQUILIBRIUM DENSITY PARAMETERS FROM PREVIOUS SIMULATION:

		if (SpecieT(1)%FluxTubeT(1)%SPINUPflagT(1) == 0) then

			do s= 1, Stot, 1
				do f= 1, SpecieT(s)%NfT(1), 1

					allocate(DensityInput(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1), &
						TemperatureInput(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1), &
						EAInertialInput(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1), &
						EAPressureInput(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1), &
						EAmagInput(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1))

				end do
			end do

			! ----------------------------------------------------

			do s= 1, Stot, 1
				write(sstring, '(I5)') s
				do f= 1, SpecieT(s)%NfT(1), 1
					write(fstring, '(I5)') f

					open (unit= 0, file= Densitydatadir &
						// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
						// '_' // 'nsnormCLBOutputfort.bin')), &
						status= 'old', form= 'unformatted', access= 'stream')
					read(0) nsnormCLBInput
					close(0)
					open (unit= 0, file= Densitydatadir &
						// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
						// '_' // 'nsnormCUBOutputfort.bin')), &
						status= 'old', form= 'unformatted', access= 'stream')
					read(0) nsnormCUBInput
					close(0)
					open (unit= 0, file= Densitydatadir &
						// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
						// '_' // 'DensityOutputRTfort.bin')), &
						status= 'old', form= 'unformatted', access= 'stream')
					read(0) DensityInput(:)
					close(0)
					open (unit= 0, file= Densitydatadir &
						// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
						// '_' // 'TemperatureOutputRTfort.bin')), &
						status= 'old', form= 'unformatted', access= 'stream')
					read(0) TemperatureInput(:)
					close(0)
					open (unit= 0, file= Densitydatadir &
						// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
						// '_' // 'EAInertialOutputRTfort.bin')), &
						status= 'old', form= 'unformatted', access= 'stream')
					read(0) EAInertialInput(:)
					close(0)
					open (unit= 0, file= Densitydatadir &
						// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
						// '_' // 'EAPressureOutputRTfort.bin')), &
						status= 'old', form= 'unformatted', access= 'stream')
					read(0) EAPressureInput(:)
					close(0)
					open (unit= 0, file= Densitydatadir &
						// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
						// '_' // 'EAmagOutputRTfort.bin')), &
						status= 'old', form= 'unformatted', access= 'stream')
					read(0) EAmagInput(:)
					close(0)
					open (unit= 0, file= Densitydatadir &
						// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
						// '_' // 'LBREPLENISHflagTOutputfort.bin')), &
						status= 'old', form= 'unformatted', access= 'stream')
					read(0) LBREPLENISHflagInput
					close(0)
					open (unit= 0, file= Densitydatadir &
						// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
						// '_' // 'UBREPLENISHflagTOutputfort.bin')), &
						status= 'old', form= 'unformatted', access= 'stream')
					read(0) UBREPLENISHflagInput
					close(0)

					! ----------------------------------------------------

					SpecieT(s)%FluxTubeT(f)%nsnormCLBInputT(1)= nsnormCLBInput
					SpecieT(s)%FluxTubeT(f)%nsnormCUBInputT(1)= nsnormCUBInput

				end do
			end do

			do s= 1, Stot, 1
				do f= 1, SpecieT(s)%NfT(1), 1

					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%DensityInputT(1)= DensityInput(Qind)
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TemperatureInputT(1)= TemperatureInput(Qind)
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%EAInertialInputT(1)= EAInertialInput(Qind)
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%EAPressureInputT(1)= EAPressureInput(Qind)
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%EAmagInputT(1)= EAmagInput(Qind)
					end do

					! ----------------------------------------------------

					if (LBREPLENISHflagInput /= SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' INCONSISTENT LBREPLENISHflagInputT FOR SPECIE= ', s, &
							' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
							' SUBROUTINE' // achar(27) // '[0m.'
					end if
					if (UBREPLENISHflagInput /= SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' INCONSISTENT UBREPLENISHflagInputT FOR SPECIE= ', s, &
							' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
							' SUBROUTINE' // achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

				end do
			end do
		end if

		! ----------------------------------------------------

		! CONSTRUCT INITIAL CONFIGURATION-SPACE GRID:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				! ----------------------------------------------------

				INITIALGRIDflag= 1
				nn= 1

				! ----------------------------------------------------

				if ((SpecieT(1)%FluxTubeT(1)%SPINUPflagT(1) == 1) .or. &
					((SpecieT(1)%FluxTubeT(1)%SPINUPflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%DYNAMICGRIDflagT(1) == 0))) then

					! Set initial grid and density profile:
					ns0(1)= ns0IC ! Reference density at lower boundary ghost cell [m^-3]
					SpecieT(s)%FluxTubeT(f)%nsnormfacT(1)= nsnormfacIC ! Macroparticle normalization constant
					zns0Neut(1)= zns0NeutIC ! Reference altitude of neutral density [km]
					ns0Neut(1)= ns0NeutIC ! Reference neutral density [m^-3]
					VperpsigmaFac(1)= VperpsigmaFacIC ! Number of MB standard deviations spanned by Vperp grid
					VparsigmaFac(1)= VparsigmaFacIC ! Number of MB standard deviations spanned by Vpar grid
					Ti(1)= TiIC ! Thermal ion temperature [K]
					Te(1)= TeIC ! Electron temperature [K]
					TNeut(1)= TNeutIC ! Neutral temperature [K]
					qGA(1)= qGAIC ! Lower boundary q value
					qGB(1)= qGBIC ! Upper boundary q value

				end if
				if ((SpecieT(1)%FluxTubeT(1)%SPINUPflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%DYNAMICGRIDflagT(1) == 1)) then

					! Set initial grid and density profile:
					ns0(1)= ns0IC ! Reference density at lower boundary ghost cell [m^-3]
					SpecieT(s)%FluxTubeT(f)%nsnormfacT(1)= nsnormfacIC ! Macroparticle normalization constant
					zns0Neut(1)= zns0NeutIC ! Reference altitude of neutral density [km]
					ns0Neut(1)= ns0NeutIC ! Reference neutral density [m^-3]
					VperpsigmaFac(1)= VperpsigmaFacIC ! Number of MB standard deviations spanned by Vperp grid
					VparsigmaFac(1)= VparsigmaFacIC ! Number of MB standard deviations spanned by Vpar grid
					Ti(1)= SpecieT(s)%FluxTubeT(f)%QCellT(1)%TemperatureInputT(1) ! Thermal ion temperature [K]
					Te(1)= TeIC ! Electron temperature [K]
					TNeut(1)= TNeutIC ! Neutral temperature [K]
					qGA(1)= qGAIC ! Lower boundary q value
					qGB(1)= qGBIC ! Upper boundary q value

				end if

				Lshell(1)= LshellIC ! L-shell [RE]
				phiLshell(1)= phiLshellIC ! invariant longitude [rads]

				! ----------------------------------------------------

				! ALLOCATE INITIALIZATION CONFIG-SPACE GRID PARAMETERS:

				allocate(qGC0(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3), &
					hqC0(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3), &
					dpC0(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3), &
					dqC0(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3), &
					dphiC0(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3), &
					rGC0(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3), &
					phiGC0(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3), &
					thetaGC0(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3), &
					ellGC0(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3), &
					qGL0(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3), &
					qGH0(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3), &
					pGC0(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3), &
					d3xC0(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3), &
					TsPerp0(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3), &
					TsPar0(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3), &
					Ts0(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3), &
					hpC0(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3), &
					hphiC0(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3), &
					Te0(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3), &
					nsnormC0(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3), &
					nsnormCNeut0(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3))

				! ----------------------------------------------------

				call ConfigGridGeneratorSub

				! ----------------------------------------------------

				allocate(SpecieT(s)%FluxTubeT(f)%QCellICT(SpecieT(s)%FluxTubeT(f)%NqICT(1)))

				! Number of time-steps for injection
				SpecieT(s)%FluxTubeT(f)%NNtT(1)= &
					SpecieT(s)%FluxTubeT(f)%NtT(1)/SpecieT(s)%FluxTubeT(f)%Q0ndatfacT(1)

				! ----------------------------------------------------

				! ALLOCATE ALL GRID PARAMETERS:

	      allocate(SpecieT(s)%FluxTubeT(f)%ndatfacGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
	      allocate(SpecieT(s)%FluxTubeT(f)%nsnormCLBGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
	      allocate(SpecieT(s)%FluxTubeT(f)%nsnormCUBGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
	      allocate(SpecieT(s)%FluxTubeT(f)%nsnormCGT( &
	        (SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))
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
	      allocate(SpecieT(s)%FluxTubeT(f)%TeGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
	        (SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))
	      allocate(SpecieT(s)%FluxTubeT(f)%nsnormCNeutGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
	        (SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))
	      allocate(SpecieT(s)%FluxTubeT(f)%dsICRGT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
	        ((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1)))
				allocate(ICbbp(((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 3)))

				! Time-step interval for moment computation (must be > 2 and excludes initial time-step)
				ndatfacG(1)= SpecieT(s)%FluxTubeT(f)%NtT(1)/SpecieT(s)%FluxTubeT(f)%NNtT(1)

				! ----------------------------------------------------

				call ConfigGridGeneratorSubB

				! ----------------------------------------------------

				SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)= ndatfacG(1)
				SpecieT(s)%FluxTubeT(f)%nsnormCLBGT(nn)= nsnormCLBG(1)
				SpecieT(s)%FluxTubeT(f)%nsnormCUBGT(nn)= nsnormCUBG(1)

				! ----------------------------------------------------

				! RE-INDEX CONFIGURATION SPACE GRID FOR A NON-COMPUTATIONAL LOWER BOUNDARY GHOST CELL:

				do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

					SpecieT(s)%FluxTubeT(f)%nsnormCGT(Qind)= nsnormC0(Qind+ 1)
					SpecieT(s)%FluxTubeT(f)%qGCGT(1, Qind)= qGC0(Qind+ 1)
					SpecieT(s)%FluxTubeT(f)%hqCGT(1, Qind)= hqC0(Qind+ 1)
					SpecieT(s)%FluxTubeT(f)%dpCGT(1, Qind)= dpC0(Qind+ 1)
					SpecieT(s)%FluxTubeT(f)%dqCGT(1, Qind)= dqC0(Qind+ 1)
					SpecieT(s)%FluxTubeT(f)%dphiCGT(1, Qind)= dphiC0(Qind+ 1)
					SpecieT(s)%FluxTubeT(f)%rGCGT(1, Qind)= rGC0(Qind+ 1)
					SpecieT(s)%FluxTubeT(f)%phiGCGT(1, Qind)= phiGC0(Qind+ 1)
					SpecieT(s)%FluxTubeT(f)%thetaGCGT(1, Qind)= thetaGC0(Qind+ 1)
					SpecieT(s)%FluxTubeT(f)%ellGCGT(1, Qind)= ellGC0(Qind+ 1)
					SpecieT(s)%FluxTubeT(f)%qGLGT(1, Qind)= qGL0(Qind+ 1)
					SpecieT(s)%FluxTubeT(f)%qGHGT(1, Qind)= qGH0(Qind+ 1)
					SpecieT(s)%FluxTubeT(f)%pGCGT(1, Qind)= pGC0(Qind+ 1)
					SpecieT(s)%FluxTubeT(f)%d3xCGT(1, Qind)= d3xC0(Qind+ 1)
					SpecieT(s)%FluxTubeT(f)%TsPerpGT(1, Qind)= TsPerp0(Qind+ 1)
					SpecieT(s)%FluxTubeT(f)%TsParGT(1, Qind)= TsPar0(Qind+ 1)
					SpecieT(s)%FluxTubeT(f)%TsGT(1, Qind)= Ts0(Qind+ 1)
					SpecieT(s)%FluxTubeT(f)%TeGT(1, Qind)= Te0(Qind+ 1)

					if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
						SpecieT(s)%FluxTubeT(f)%nsnormCNeutGT(1, Qind)= nsnormCNeut0(Qind+ 1)
					end if

					! ----------------------------------------------------

					! Compute field line arc length of BBELF wave-heating region
					SpecieT(s)%FluxTubeT(f)%dsICRGT(1, Qind)= &
						SpecieT(s)%FluxTubeT(f)%hqCGT(1, Qind)*SpecieT(s)%FluxTubeT(f)%dqCGT(1, Qind)

					! ----------------------------------------------------

					! DIAGNOSTIC FLAGS FOR CONSISTENT PHI VALUE:

					if (SpecieT(s)%FluxTubeT(f)%phiGCGT(1, Qind) /= phiGC0(1)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT INITIAL phiGCGT VALUE= ', &
							SpecieT(s)%FluxTubeT(f)%phiGCGT(1, Qind), ' AND phiGC0T VALUE= ', phiGC0(Qind), &
							' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND Q CELL= ', Qind, &
							' IN SIMULATION PARAMAETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if (SpecieT(s)%FluxTubeT(f)%pGCGT(1, Qind) /= pGC0(1)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT INITIAL pGCGT VALUE= ', &
							SpecieT(s)%FluxTubeT(f)%pGCGT(1, Qind), ' AND phiGC0T VALUE= ', pGC0(Qind), &
							' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND Q CELL= ', Qind, &
							' IN SIMULATION PARAMAETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

				end do

				! ----------------------------------------------------

				call VelGridGeneratorSub

				! ----------------------------------------------------

			end do
		end do

		! ----------------------------------------------------

		! SET NOISE FILTERS FOR ION MOMENTS:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				if (SpecieT(1)%FluxTubeT(1)%MOMENTFILTERflagT(1) == 1) then

					! ----------------------------------------------------

					! DIAGNOSTIC FLAGS CONSISTENT MOMENT FILTER MOVING AVERAGE POINTS:

					if (M0MAfilterPt > ((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ION DENSITY MOMENT FILTER', &
							' MOVING AVERAGE POINT IS GREATER THAN NUMBER OF CONFIGURATION-SPACE GRID CELLS', &
							' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
							' SUBROUTINE' // achar(27) // '[0m.'
					end if

					if (M1Perp1MAfilterPt > ((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ION PERP1 VELOCITY MOMENT FILTER', &
							' MOVING AVERAGE POINT IS GREATER THAN NUMBER OF CONFIGURATION-SPACE GRID CELLS', &
							' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
							' SUBROUTINE' // achar(27) // '[0m.'
					end if

					if (M1Perp2MAfilterPt > ((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ION PERP2 VELOCITY MOMENT FILTER', &
							' MOVING AVERAGE POINT IS GREATER THAN NUMBER OF CONFIGURATION-SPACE GRID CELLS', &
							' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
							' SUBROUTINE' // achar(27) // '[0m.'
					end if

					if (M1ParMAfilterPt > ((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ION PARALLEL VELOCITY MOMENT FILTER', &
							' MOVING AVERAGE POINT IS GREATER THAN NUMBER OF CONFIGURATION-SPACE GRID CELLS', &
							' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
							' SUBROUTINE' // achar(27) // '[0m.'
					end if

					if (M2ParMAfilterPt > ((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ION PARALLEL ENERGY MOMENT FILTER', &
							' MOVING AVERAGE POINT IS GREATER THAN NUMBER OF CONFIGURATION-SPACE GRID CELLS', &
							' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
							' SUBROUTINE' // achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

				end if

			end do
		end do

		! ----------------------------------------------------

		! SET NOISE LIMIT FOR ION POPULATION STATISTICS:

		if (SpecieT(1)%FluxTubeT(1)%IONNOISEflagT(1) == 1) then
			do s= 1, Stot, 1
				do f= 1, SpecieT(s)%NfT(1), 1 ! 0.2d0 limits N= 25 and 0.141 limits N= 50
					IonNoiseLimit(1)= 1d0/(sqrt(IonNoiseLimitNph)) ! Noise limit for ion distribution functions
				end do
			end do
		end if
		if (SpecieT(1)%FluxTubeT(1)%IONNOISEflagT(1) == 0) then
			do s= 1, Stot, 1
				do f= 1, SpecieT(s)%NfT(1), 1
					IonNoiseLimit(1)= 0d0 ! Noise limit for ion distribution functions
				end do
			end do
		end if
		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				SpecieT(s)%FluxTubeT(f)%IonNoiseLimitT= IonNoiseLimit

				! ----------------------------------------------------

				! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

				if ((IonNoiseLimit(1) /= SpecieT(s)%FluxTubeT(f)%IonNoiseLimitT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%IonNoiseLimitT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%IonNoiseLimitT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' IonNoiseLimitT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				! ----------------------------------------------------

			end do
		end do

		! ----------------------------------------------------

		! SET NOISE LIMIT FOR ENA POPULATION STATISTICS:

		if (SpecieT(1)%FluxTubeT(1)%ENANOISEflagT(1) == 1) then
			do s= 1, Stot, 1
				do f= 1, SpecieT(s)%NfT(1), 1 ! 0.2d0 limits N= 25 and 0.141 limits N= 50
						ENANoiseLimit(1)= 1d0/(sqrt(ENANoiseLimitNph)) ! Noise limit for ENA distribution functions
				end do
			end do
		end if
		if (SpecieT(1)%FluxTubeT(1)%ENANOISEflagT(1) == 0) then
			do s= 1, Stot, 1
				do f= 1, SpecieT(s)%NfT(1), 1
					ENANoiseLimit(1)= 0d0 ! Noise limit for ENA distribution functions
				end do
			end do
		end if
		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				SpecieT(s)%FluxTubeT(f)%ENANoiseLimitT= ENANoiseLimit

				! ----------------------------------------------------

				! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

				if ((ENANoiseLimit(1) /= SpecieT(s)%FluxTubeT(f)%ENANoiseLimitT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%ENANoiseLimitT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%ENANoiseLimitT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ENANoiseLimitT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				! ----------------------------------------------------

			end do
		end do

		! ----------------------------------------------------

		! COMPUTE NUMBER OF MPI RANKS:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				call MPI_COMM_SIZE(MPI_COMM_WORLD, ranksize(1), ierr)
			end do
		end do

		allocate(RankT(ranksize(1)))

		! ----------------------------------------------------

		! SET REFERENCE PARALLEL POTENTIAL DROPS:

		! Consider making this time and space dependent
		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				if (SpecieT(s)%FluxTubeT(f)%EPARflagT(1) == 1) then

					! ----------------------------------------------------

					EPar0p(1)= PhiPar0p/dPhiPar0p ! [V/m]

					! Create nested derived data types
					SpecieT(s)%FluxTubeT(f)%PhiPar0BT= EPar0p(1)

					! ----------------------------------------------------

					! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

					if ((EPar0p(1) /= SpecieT(s)%FluxTubeT(f)%PhiPar0BT(1)) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%PhiPar0BT(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%PhiPar0BT(1))) &
						.eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' PhiPar0BT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', AND FLUX TUBE= ', f, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

				end if
			end do
		end do

		! ----------------------------------------------------

		! SET ICR HEATING PARAMETERS:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				if (SpecieT(s)%FluxTubeT(f)%ICRflagT(1) == 1) then

					allocate(SpecieT(s)%FluxTubeT(f)%lambdaPerppT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
	          ((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1)))
					allocate(SpecieT(s)%FluxTubeT(f)%EtaLHpT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
	          ((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1)))
					allocate(SpecieT(s)%FluxTubeT(f)%XiPerp1pT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
	          ((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1)))
					allocate(SpecieT(s)%FluxTubeT(f)%XiPerp2pT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
	          ((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1)))
					allocate(SpecieT(s)%FluxTubeT(f)%S0pT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
	          ((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1)))
					allocate(SpecieT(s)%FluxTubeT(f)%OmegaG0pT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
	          ((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1)))
					allocate(SpecieT(s)%FluxTubeT(f)%ChiPerp1pT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
	          ((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1)))
					allocate(SpecieT(s)%FluxTubeT(f)%ChiPerp2pT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
						((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1)))

					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
						do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1

							! ----------------------------------------------------

							OmegaG0p(1)= 2d0*pi*f0p ! Ref. gyrofreq. corresponding to EPerp0 [rads/s]

							! Create nested derived data types
							SpecieT(s)%FluxTubeT(f)%lambdaPerppT(nn, Qind)= lambdaPerpp
							SpecieT(s)%FluxTubeT(f)%EtaLHpT(nn, Qind)= EtaLHp
							SpecieT(s)%FluxTubeT(f)%XiPerp1pT(nn, Qind)= XiPerp1p
							SpecieT(s)%FluxTubeT(f)%XiPerp2pT(nn, Qind)= XiPerp2p
							SpecieT(s)%FluxTubeT(f)%S0pT(nn, Qind)= S0p
							SpecieT(s)%FluxTubeT(f)%OmegaG0pT(nn, Qind)= OmegaG0p(1)
							SpecieT(s)%FluxTubeT(f)%ChiPerp1pT(nn, Qind)= ChiPerp1p
							SpecieT(s)%FluxTubeT(f)%ChiPerp2pT(nn, Qind)= ChiPerp2p

							! ----------------------------------------------------

							! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

							if ((lambdaPerpp /= SpecieT(s)%FluxTubeT(f)%lambdaPerppT(nn, Qind)) .or. &
								(size(SpecieT(s)%FluxTubeT(f)%lambdaPerppT(:, :)) /= &
									((SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)*((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) .or. &
								(isnan(real(SpecieT(s)%FluxTubeT(f)%lambdaPerppT(nn, Qind))) &
								.eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' lambdaPerppT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE ', &
									' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
									', STATISTICAL TIME-STEP= ', nn, &
									', AND Qind= ', Qind, ' IN SIMULATION', &
									' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((EtaLHp /= SpecieT(s)%FluxTubeT(f)%EtaLHpT(nn, Qind)) .or. &
								(size(SpecieT(s)%FluxTubeT(f)%EtaLHpT(:, :)) /= &
									((SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)*((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) .or. &
								(isnan(real(SpecieT(s)%FluxTubeT(f)%EtaLHpT(nn, Qind))) &
								.eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' EtaLHpT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE ', &
									' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
									', STATISTICAL TIME-STEP= ', nn, &
									', AND Qind= ', Qind, ' IN SIMULATION', &
									' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((XiPerp1p /= SpecieT(s)%FluxTubeT(f)%XiPerp1pT(nn, Qind)) .or. &
								(size(SpecieT(s)%FluxTubeT(f)%XiPerp1pT(:, :)) /= &
									((SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)*((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) .or. &
								(isnan(real(SpecieT(s)%FluxTubeT(f)%XiPerp1pT(nn, Qind))) &
								.eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' XiPerp1pT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE ', &
									' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
									', STATISTICAL TIME-STEP= ', nn, &
									', AND Qind= ', Qind, ' IN SIMULATION', &
									' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((XiPerp2p /= SpecieT(s)%FluxTubeT(f)%XiPerp2pT(nn, Qind)) .or. &
								(size(SpecieT(s)%FluxTubeT(f)%XiPerp2pT(:, :)) /= &
									((SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)*((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) .or. &
								(isnan(real(SpecieT(s)%FluxTubeT(f)%XiPerp2pT(nn, Qind))) &
								.eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' XiPerp2pT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE ', &
									' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
									', STATISTICAL TIME-STEP= ', nn, &
									', AND Qind= ', Qind, ' IN SIMULATION', &
									' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((S0p /= SpecieT(s)%FluxTubeT(f)%S0pT(nn, Qind)) .or. &
								(size(SpecieT(s)%FluxTubeT(f)%S0pT(:, :)) /= &
									((SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)*((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) .or. &
								(isnan(real(SpecieT(s)%FluxTubeT(f)%S0pT(nn, Qind))) &
								.eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' S0pT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE ', &
									' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
									', STATISTICAL TIME-STEP= ', nn, &
									', AND Qind= ', Qind, ' IN SIMULATION', &
									' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((OmegaG0p(1) /= SpecieT(s)%FluxTubeT(f)%OmegaG0pT(nn, Qind)) .or. &
								(size(SpecieT(s)%FluxTubeT(f)%OmegaG0pT(:, :)) /= &
									((SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)*((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) .or. &
								(isnan(real(SpecieT(s)%FluxTubeT(f)%OmegaG0pT(nn, Qind))) &
								.eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' OmegaG0pT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE ', &
									' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
									', STATISTICAL TIME-STEP= ', nn, &
									', AND Qind= ', Qind, ' IN SIMULATION', &
									' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((ChiPerp1p /= SpecieT(s)%FluxTubeT(f)%ChiPerp1pT(nn, Qind)) .or. &
								(size(SpecieT(s)%FluxTubeT(f)%ChiPerp1pT(:, :)) /= &
									((SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)*((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) .or. &
								(isnan(real(SpecieT(s)%FluxTubeT(f)%ChiPerp1pT(nn, Qind))) &
								.eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ChiPerp1pT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE ', &
									' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
									', STATISTICAL TIME-STEP= ', nn, &
									', AND Qind= ', Qind, ' IN SIMULATION', &
									' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((ChiPerp2p /= SpecieT(s)%FluxTubeT(f)%ChiPerp2pT(nn, Qind)) .or. &
								(size(SpecieT(s)%FluxTubeT(f)%ChiPerp2pT(:, :)) /= &
									((SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)*((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) .or. &
								(isnan(real(SpecieT(s)%FluxTubeT(f)%ChiPerp2pT(nn, Qind))) &
								.eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ChiPerp2pT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE ', &
									' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
									', STATISTICAL TIME-STEP= ', nn, &
									', AND Qind= ', Qind, ' IN SIMULATION', &
									' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
							end if

							if (SpecieT(s)%FluxTubeT(f)%ICRCOHERENCEflagT(1) == 1) then
								if (SpecieT(s)%FluxTubeT(f)%XiPerp1pT(nn, Qind) /= &
									SpecieT(s)%FluxTubeT(f)%XiPerp2pT(nn, Qind)) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' INCONSISTENT XiPerp1 and XiPerp2 VALUES FOR COHERENT ICR HEATING', &
										' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
										', STATISTICAL TIME-STEP= ', nn, &
										', AND Qind= ', Qind, ' IN SIMULATION', &
										' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
								end if
								if (SpecieT(s)%FluxTubeT(f)%ChiPerp1pT(nn, Qind) /= &
									SpecieT(s)%FluxTubeT(f)%ChiPerp2pT(nn, Qind)) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' INCONSISTENT ChiPerp1 and ChiPerp2 VALUES FOR COHERENT ICR HEATING', &
										' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
										', STATISTICAL TIME-STEP= ', nn, &
										', AND Qind= ', Qind, ' IN SIMULATION', &
										' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
								end if
							end if

							! ----------------------------------------------------

						end do
					end do
				end if
			end do
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

		! PRINT TERMINAL EXPORT PARAMETERS:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				if (rank == 0) then

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *) '----------------------------------------------------'
					write(*, *) 'Kinetic model of Auroral ion OutflowS (KAOS)'
					write(*, *) 'Robert M. Albarran II, Ph.D. contact: albarran1@atmos.ucla.edu'
					write(*, *) 'Department of Atmospheric and Oceanic Sciences, University of California, Los Angeles'

		      call date_and_time(datechar(1), datechar(2), datechar(3), dateint)
					write(monthstring, '(i10)') dateint(2)
					write(daystring, '(i10)') dateint(3)
					write(yearstring, '(i10)') dateint(1)
					write(hourstring, '(i10)') dateint(5)
					write(minutestring, '(i10)') dateint(6)
					write(secondstring, '(i10)') dateint(7)
	 				write(*, *) trim('Month= ' // adjustl(monthstring))
					write(*, *) trim('Day= ' // adjustl(daystring))
					write(*, *) trim('Year= ' // adjustl(yearstring))
					write(*, *) trim('Hour= ' // adjustl(hourstring))
					write(*, *) trim('Minute= ' // adjustl(minutestring))
					write(*, *) trim('Second= ' // adjustl(secondstring))

					write(*, *)
					write(*, *) '----------------------------------------------------'
					if (SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1) == 1) then
						write(*, *) trim('KAOS SPIN-UP SIMULATION')
					end if
					if (SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1) == 0) then
						write(*, *) trim('KAOS RESTART SIMULATION')
					end if

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *) 'SIMULATION PARAMETERS:'
					write(*, *)

					write(paramstring, '(i10)') ranksize(1)
					write(*, *) trim('Number of MPI ranks= ' // adjustl(paramstring))

					if (qGAIC <= 0d0) then
						write(*, *) 'Southern Magnetic Hemisphere'
					else if (qGAIC > 0d0) then
						write(*, *) 'Northern Magnetic Hemisphere'
					end if

					write(*, *) trim('Data Export Path= ' // adjustl(dataexportdir))

					if (SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1) == 1) then
						write(*, *) trim('Export Spin-up Data Directory Path= ' // adjustl(Densitydatadir))
					end if
					if (SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1) == 0) then
						write(*, *) trim('Import Spin-up Data Directory Path= ' // adjustl(Densitydatadir))
					end if

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *) 'CONFIGURATION-SPACE PARAMETERS:'
					write(*, *)

					write(paramstring, '(D10.2)') LshellIC
					write(*, *) trim('Initial L-shell [RE]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') phiLshellIC
					write(*, *) trim('Initial Invariant Longitude [rads]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') (SpecieT(s)%FluxTubeT(f)%rGCGT(1, SpecieT(s)%FluxTubeT(f)%NqLBT(1))- RE)*1d-3
					write(*, *) trim('Lower Grid Altitude [km]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') (SpecieT(s)%FluxTubeT(f)%rGCGT(1, SpecieT(s)%FluxTubeT(f)%NqUBT(1))- RE)*1d-3
					write(*, *) trim('Upper Grid Altitude [km]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') (SpecieT(s)%FluxTubeT(f)%rGCGT(1, SpecieT(s)%FluxTubeT(f)%NqLBT(1))- RE)*1d-3
					write(*, *) trim('Lower Ion Injection Altitude [km]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') (SpecieT(s)%FluxTubeT(f)%rGCGT(1, SpecieT(s)%FluxTubeT(f)%NqUBT(1))- RE)*1d-3
					write(*, *) trim('Upper Ion Initialization Altitude [km]= ' // adjustl(paramstring))

					write(paramstring, '(i10)') Stot
					write(*, *) trim('Total Number of Particle Species= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(1)%NfT(1)
					write(*, *) trim('Total Number of Flux Tubes= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%NqGT(1)
					write(*, *) trim('Total Number of Ion/ENA Config-Space Grid Cells= ' // adjustl(paramstring))

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *) 'ION VELOCITY-SPACE PARAMETERS:'
					write(*, *)

					if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%V2PerpCellT(1, 1, 1)%Vperp1GLT(1)*1d-3
						write(*, *) trim('Lower Vperp1 Limit [km/s]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%V2PerpCellT(1, 1, 1)%Vperp2GLT(1)*1d-3
						write(*, *) trim('Lower Vperp2 Limit [km/s]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%V2PerpCellT( &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), 1, 1)%Vperp1GHT(1)*1d-3
						write(*, *) trim('Upper Vperp1 Limit [km/s]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%V2PerpCellT( &
							1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), 1)%Vperp2GHT(1)*1d-3
						write(*, *) trim('Upper Vperp2 Limit [km/s]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%V2PerpCellT(1, 1, 1)%VparGLT(1)*1d-3
						write(*, *) trim('Lower Vpar Limit [km/s]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%V2PerpCellT( &
							1, 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1))%VparGHT(1)*1d-3
						write(*, *) trim('Upper Vpar Limit [km/s]= ' // adjustl(paramstring))

						write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1)
						write(*, *) trim('Total Number of Ion Vperp1 Grid Cells= ' // adjustl(paramstring))
						write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1)
						write(*, *) trim('Total Number of Ion Vperp2 Grid Cells= ' // adjustl(paramstring))
					else
						write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1)
						write(*, *) trim('Total Number of Ion Vperp Grid Cells= ' // adjustl(paramstring))
					end if
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)
					write(*, *) trim('Total Number of Ion Vpar Grid Cells= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') Vperp12sigma*1d-3
					write(*, *) trim('Vperp1 and Vperp2 Standard Deviation [km/s]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') VperpsigmaFac
					write(*, *) trim('Number of Linear Vperp1 and Vperp2 Standard Deviations &
						Spanned by Linear Grid= ' // adjustl(paramstring))
					write(paramstring, '(i10)') Vperp12NlinRange
					write(*, *) trim('Number of Vperp1 and Vperp2 Linear Grid Cells Spanning to &
						Last Standard Deviation= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') Vparsigma*1d-3
					write(*, *) trim('Vpar Standard Deviation [km/s]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') VparsigmaFac
					write(*, *) trim('Number of Linear Vpar Standard Deviations &
						Spanned by Linear Grid= ' // adjustl(paramstring))
					write(paramstring, '(i10)') VparNlinRange
					write(*, *) trim('Number of Vpar Linear Grid Cells Spanning to &
						Last Standard Deviation= ' // adjustl(paramstring))

					if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then

						write(*, *)
						write(*, *) '----------------------------------------------------'
						write(*, *) 'ENA VELOCITY-SPACE PARAMETERS:'
						write(*, *)

						write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1)
						write(*, *) trim('Total Number of ENA Vp Grid Cells= ' // adjustl(paramstring))
						write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1)
						write(*, *) trim('Total Number of ENA Vq Grid Cells= ' // adjustl(paramstring))
						write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)
						write(*, *) trim('Total Number of ENA Vphi Grid Cells= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') Vpphisigma*1d-3
						write(*, *) trim('Vp and Vphi Standard Deviation [km/s]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') VpphisigmaFac
						write(*, *) trim('Number of Linear Vp and Vpphi Standard Deviations &
							Spanned by Linear Grid= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') VpphiNlinRange
						write(*, *) trim('Number of Vp and Vphi Linear Grid Cells Spanning to &
							Last Standard Deviation= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') Vqsigma*1d-3
						write(*, *) trim('Vq Standard Deviation [km/s]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') VqsigmaFac
						write(*, *) trim('Number of Linear Vq Standard Deviations &
							Spanned by Linear Grid= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') VqNlinRange
						write(*, *) trim('Number of Vq Linear Grid Cells Spanning to &
							Last Standard Deviation= ' // adjustl(paramstring))
					end if

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *) 'TIME PARAMETERS:'
					write(*, *)

					write(paramstring, '(D10.2)') SpecieT(s)%AT
					write(*, *) trim('Beginning Simulation Time [s]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%BT
					write(*, *) trim('End Simulation Time [s]= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%NtT(1)
					write(*, *) trim('Total Number of Base Time-Steps= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%hT
					write(*, *) trim('Base Time-Step Duration [s]= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%NNtT(1)
					write(*, *) trim('Total Number of Statistical Time-Steps (constant)= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%ndatfacGT(1)
					write(*, *) trim('Initial Number of Base Time-Steps per Statistical Time-Step= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%ndatfacGT(1)*SpecieT(s)%hT
					write(*, *) trim('Initial Statistical Time-Step Duration [s]= ' // adjustl(paramstring))

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *) 'ION SIMULATION FLAGS:'
					write(*, *)

					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%IONNOISEflagT(1)
					write(*, *) trim('Ion Statistical Noise Reduction Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1)
					write(*, *) trim('Ion Fluid Data Export Flag= ' // adjustl(paramstring))

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *) 'BOUNDARY CONDITIONS FLAGS:'
					write(*, *)

					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%LBCONDITIONflagT(1)
					write(*, *) trim('Lower Boundary Condition Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%UBCONDITIONflagT(1)
					write(*, *) trim('Upper Boundary Condition Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1)
					write(*, *) trim('Lower Boundary Density Replenish Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1)
					write(*, *) trim('Upper Boundary Density Replenish Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%STATICINJECTIONflagT(1)
					write(*, *) trim('Static Injection Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1)
					write(*, *) trim('Spin up Flag= ' // adjustl(paramstring))

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *) 'ION FORCE FLAGS:'
					write(*, *)

					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%ICRflagT(1)
					write(*, *) trim('Ion Cyclotron Wave Heating Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%ICRCOHERENCEflagT(1)
					write(*, *) trim('Coherent Ion Cyclotron Wave Heating Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%MIRRORflagT(1)
					write(*, *) trim('Mirror Force Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%GRAVflagT(1)
					write(*, *) trim('Gravitational Force Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%EAMBflagT(1)
					write(*, *) trim('Ambipolar Electric Field Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%EAMBSELFCONSISTflagT(1)
					write(*, *) trim('Self-Consistent Ambipolar Electric Field Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1)
					write(*, *) trim('Self-Consistent Sign of Ambipolar Electric Field Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%EAINERTIALflagT(1)
					write(*, *) trim('Inertial Term of Ambipolar Electric Field Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%EAPRESSUREflagT(1)
					write(*, *) trim('Pressure Term of Ambipolar Electric Field Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%EPARflagT(1)
					write(*, *) trim('Parallel Electric Field Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%CONVECTIONflagT(1)
					write(*, *) trim('Flux-tube Convection Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%DYNAMICGRIDflagT(1)
					write(*, *) trim('Dynamic Resolution Flag= ' // adjustl(paramstring))

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *) 'ION MOMENT FLAGS:'
					write(*, *)

					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%MOMENTFILTERflagT(1)
					write(*, *) trim('Moment Filter Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(1)
					write(*, *) trim('2D Perpendicular Velocity Space Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1)
					write(*, *) trim('Ion Phase-Space Distribution Function Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT(1)
					write(*, *) trim('Ion Phase-Space Density Moment Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%PHASEVELPERPIONMOMENTflagT(1)
					write(*, *) trim('Ion Phase-Space Perpendicular Velocity Moment Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%PHASEVELPARIONMOMENTflagT(1)
					write(*, *) trim('Ion Phase-Space Parallel Velocity Moment Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(1)
					write(*, *) trim('Ion Phase-Space Total Thermal Energy Moment Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTflagT(1)
					write(*, *) trim('Ion Phase-Space Perpendicular Thermal Energy Moment Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%PHASEENERGYPARIONMOMENTflagT(1)
					write(*, *) trim('Ion Phase-Space Parallel Thermal Energy Moment Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTCENTERflagT(1)
					write(*, *) trim('Ion Phase-Space Perpendicular Thermal Energy Moment Centering Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%PHASEENERGYPARIONMOMENTCENTERflagT(1)
					write(*, *) trim('Ion Phase-Space Parallel Thermal Energy Moment Centering Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%FLUIDIONREFflagT(1)
					write(*, *) trim('Ion Reference Moments Flag= ' // adjustl(paramstring))

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *) 'ENA FLAGS:'
					write(*, *)

					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%ENANOISEflagT(1)
					write(*, *) trim('ENA Statistical Noise Reduction Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1)
					write(*, *) trim('Charge-Exchange Flag= ' // adjustl(paramstring))

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *) 'ION INITIALIZATION PARAMETERS:'
					write(*, *)

					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%NqLBT(1)
					write(*, *) trim('Ion Initialization LB Grid Cell= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%NqUBT(1)
					write(*, *) trim('Ion Initialization UB Grid Cell= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') (SpecieT(s)%FluxTubeT(f)%rGCGT(1, SpecieT(s)%FluxTubeT(f)%NqLBT(1))- RE)*1d-3
					write(*, *) trim('Ion Initialization Lower Altitude Boundary [km]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') (SpecieT(s)%FluxTubeT(f)%rGCGT(1, SpecieT(s)%FluxTubeT(f)%NqUBT(1))- RE)*1d-3
					write(*, *) trim('Ion Initialization Upper Altitude Boundary [km]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') (SpecieT(s)%FluxTubeT(f)%rGCGT(1, SpecieT(s)%FluxTubeT(f)%NqLBT(1))- RE)*1d-3
					write(*, *) trim('Ion Initialization Reference Altitude [km]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') ns0IC
					write(*, *) trim('Ion Initialization Reference Density [m^-3]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%TsGT(1, SpecieT(s)%FluxTubeT(f)%NqLBT(1))
					write(*, *) trim('Ion Initialization LB Temperature [K]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%TsPerpGT(1, SpecieT(s)%FluxTubeT(f)%NqLBT(1))
					write(*, *) trim('Ion Initialization LB Perpendicular Temperature [K]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%TsParGT(1, SpecieT(s)%FluxTubeT(f)%NqLBT(1))
					write(*, *) trim('Ion Initialization LB Parallel Temperature [K]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%TeGT(1, SpecieT(s)%FluxTubeT(f)%NqLBT(1))
					write(*, *) trim('Electron Initialization LB Temperature [K]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') dNTe
					write(*, *) trim('Additive Increment of Electron Temperature on Statistical Time-steps [K]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') dNTeEND
					write(*, *) trim('Additive Increment Cap of Electron Temperature on Statistical Time-steps [K]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') ns0(1)
					write(*, *) trim('Reference Density Grid Cell= ' // adjustl(paramstring))

					if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then

						write(*, *)
						write(*, *) '----------------------------------------------------'
						write(*, *) 'NEUTRAL OXYGEN INITIALIZATION PARAMETERS:'
						write(*, *)

						write(paramstring, '(D10.2)') (zns0neutIC- RE)*1d-3
						write(*, *) trim('Neutral Oxygen Initialization Reference Altitude [km]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') ns0neutIC
						write(*, *) trim('Neutral Oxygen Initialization Reference Density [m^-3]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') TNeut
						write(*, *) trim('Neutral Oxygen Initialization Temperature [K]= ' // adjustl(paramstring))

					end if

					if (SpecieT(s)%FluxTubeT(f)%ICRflagT(1) == 1) then

						write(*, *)
						write(*, *) '----------------------------------------------------'
						write(*, *) 'ICR HEATING PARAMETERS:'
						write(*, *)

						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%lambdaPerppT(1, 1)
						write(*, *) trim('BBELF Wavelength [m]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%EtaLHpT(1, 1)
						write(*, *) trim('BBELF Wave Power LHP Fraction= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%XiPerp1pT(1, 1)
						write(*, *) trim('BBELF Wave Power Fraction Along Vperp1= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%XiPerp2pT(1, 1)
						write(*, *) trim('BBELF Wave Power Fraction Along Vperp2= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%S0pT(1, 1)
						write(*, *) trim('Wave Spectral Energy Density [(V^2/m^2)/Hz]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%OmegaG0pT(1, 1)/(2d0*pi)
						write(*, *) trim('Reference Ion Cyclotron Frequency [Hz]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%ChiPerp1pT(1, 1)
						write(*, *) trim('Wave Spectral Index Along Vperp1= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%ChiPerp2pT(1, 1)
						write(*, *) trim('Wave Spectral Index Along Vperp2= ' // adjustl(paramstring))

					end if

					if (SpecieT(s)%FluxTubeT(f)%EPARflagT(1) == 1) then

						write(*, *)
						write(*, *) '----------------------------------------------------'
						write(*, *) 'PARALLEL POTENTIAL PARAMETERS:'
						write(*, *)

						write(paramstring, '(i10)') Eparlim
						write(*, *) trim('Computational Time-step of Activated Potential Drop= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') PhiPar0p
						write(*, *) trim('Reference Parallel Potential Drop [V]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') dPhiPar0p
						write(*, *) trim('Reference Parallel Potential Distance [m]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%PhiPar0BT
						write(*, *) trim('Reference Parallel Electric Field [V/m]= ' // adjustl(paramstring))

					end if

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *) 'ION LIMITS:'
					write(*, *)

					if (SpecieT(1)%FluxTubeT(1)%MOMENTFILTERflagT(1) == 1) then
						write(paramstring, '(i10)') M0MAfilterPt
						write(*, *) trim('Density Moment Moving Average Point= ' // adjustl(paramstring))
						write(paramstring, '(i10)') M1Perp1MAfilterPt
						write(*, *) trim('Perp1 Velocity Moment Moving Average Point= ' // adjustl(paramstring))
						write(paramstring, '(i10)') M1Perp2MAfilterPt
						write(*, *) trim('Perp2 Velocity Moment Moving Average Point= ' // adjustl(paramstring))
						write(paramstring, '(i10)') M1ParMAfilterPt
						write(*, *) trim('Parallel Velocity Moment Moving Average Point= ' // adjustl(paramstring))
						write(paramstring, '(i10)') M2ParMAfilterPt
						write(*, *) trim('Parallel Energy Moment Moving Average Point= ' // adjustl(paramstring))
					end if
					write(paramstring, '(D10.2)') IonNoiseLimitNph
					write(*, *) trim('Minimum Phase-Space Ion Macroparticle Number= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%IonNoiseLimitT(1)
					write(*, *) trim('Ion Noise Limit= ' // adjustl(paramstring))

					if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then

						write(*, *)
						write(*, *) '----------------------------------------------------'
						write(*, *) 'ENA LIMITS:'
						write(*, *)

						write(paramstring, '(D10.2)') ENANoiseLimitNph
						write(*, *) trim('Minimum Phase-Space ENA Macroparticle Number= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%ENANoiseLimitT(1)
						write(*, *) trim('ENA Noise Limit= ' // adjustl(paramstring))

					end if

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *)

				end if
			end do
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				if (rank == 0) then
					call cpu_time(S1End)
					write(S1string, '(i10)')  nint(S1End)
					write(*, *) trim('%% 1- RANK= ' // adjustl(rankstring)) // &
						trim(', REAL-TIME= ' // adjustl(S1string)) // &
						trim(' s. SIMULATION PARAMETERIZATION COMPLETE %%')
				end if
			end do
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine SimParameterizationSub

end module SimParameterization
