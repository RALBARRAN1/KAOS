module DataExport2

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.0.2- DATA EXPORT 2:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! EXPORT DATA AS BINARY FILES:

	subroutine DataExport2Sub

		! ----------------------------------------------------

		if ((SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1) == 1)  .and. &
			(rank == 0)) then

			write(nnstring, '(I5)') nn
			write(sstring, '(I5)') s
			write(fstring, '(I5)') f

			expstring= adjustl(adjustr(rankstring) &
				// '_' // adjustl(adjustr(nnstring) &
				// '_' // adjustl(adjustr(sstring) &
				// '_' // adjustl(adjustr(fstring) // '_'))))

			TimeTfile= adjustl(adjustr(expstring) // adjustl(adjustr('TimeTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
				adjustl(adjustr(TimeTfile))), &
				status= 'replace', form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%TimeT(nn)
			close(expint)

			ndatfacGTfile= adjustl(adjustr(expstring) // adjustl(adjustr('ndatfacGTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
				adjustl(adjustr(ndatfacGTfile))), &
				status= 'replace', form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%ndatfacGT(nn)
			close(expint)

			ns0GTfile= adjustl(adjustr(expstring) // adjustl(adjustr('ns0GTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
				adjustl(adjustr(ns0GTfile))), &
				status= 'replace', form= 'unformatted', access= 'stream')
			write(expint) ns0(1)
			close(expint)

			qGCGTfile= adjustl(adjustr(expstring) // adjustl(adjustr('qGCGTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
				adjustl(adjustr(qGCGTfile))), &
				status= 'replace', form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%qGCGT(nn, :)
			close(expint)

			pGCGTfile= adjustl(adjustr(expstring) // adjustl(adjustr('pGCGTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
				adjustl(adjustr(pGCGTfile))), &
				status= 'replace', form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%pGCGT(nn, :)
			close(expint)

			rGCGTfile= adjustl(adjustr(expstring) // adjustl(adjustr('rGCGTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
				adjustl(adjustr(rGCGTfile))), &
				status= 'replace', form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%rGCGT(nn, :)
			close(expint)

			phiGCGTfile= adjustl(adjustr(expstring) // adjustl(adjustr('phiGCGTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
				adjustl(adjustr(phiGCGTfile))), &
				status= 'replace', form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%phiGCGT(nn, :)
			close(expint)

			thetaGCGTfile= adjustl(adjustr(expstring) // adjustl(adjustr('thetaGCGTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
				adjustl(adjustr(thetaGCGTfile))), &
				status= 'replace', form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, :)
			close(expint)

			TsGTfile= adjustl(adjustr(expstring) // adjustl(adjustr('TsGTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
				adjustl(adjustr(TsGTfile))), &
				status= 'replace', form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%TsGT(nn, :)
			close(expint)

			!do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
			!	TeNTfile= adjustl(adjustr(expstring) // adjustl(adjustr('TeNTfort.bin')))
			!	open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
			!		adjustl(adjustr(TeNTfile))), &
			!		status= 'replace', form= 'unformatted', access= 'stream')
			!	write(expint) SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn)
			!	close(expint)
			!end do

		end if

		! ----------------------------------------------------

		if ((SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1) == 1)  .and. &
			(rank == 0)) then
			if (SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1) == 1) then
				if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then
					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

						write(nnstring, '(I5)') nn
						write(sstring, '(I5)') s
						write(fstring, '(I5)') f
						write(Qindstring, '(I5)') Qind

						expstring= adjustl(adjustr(rankstring) // '_' // &
							adjustl(adjustr(nnstring) // '_' // &
							adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // &
							adjustl(adjustr(Qindstring) // '_')))))

						TeNTfile= adjustl(adjustr(expstring) // adjustl(adjustr('TeNTfort.bin')))
						open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
							adjustl(adjustr(TeNTfile))), &
							status= 'replace', form= 'unformatted', access= 'stream')
						write(expint) SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(nn)
						close(expint)

						N2PerpphRTfile= adjustl(adjustr(expstring) // &
							adjustl(adjustr('N2PerpphRTfort.bin')))
						open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
							adjustl(adjustr(N2PerpphRTfile))), status= 'replace', &
							form= 'unformatted', access= 'stream')
						write(expint) SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%N2PerpphRTp(:, :, :, nn)
						close(expint)

					end do
				end if
			end if
		end if

		if ((SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1) == 1)  .and. &
			(rank == 0)) then
			if (SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1) == 1) then
				if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 0) then
					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

						write(nnstring, '(I5)') nn
						write(sstring, '(I5)') s
						write(fstring, '(I5)') f
						write(Qindstring, '(I5)') Qind

						expstring= adjustl(adjustr(rankstring) // '_' // &
							adjustl(adjustr(nnstring) // '_' // &
							adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // &
							adjustl(adjustr(Qindstring) // '_')))))

						NphRTfile= adjustl(adjustr(expstring) // &
							adjustl(adjustr('NphRTfort.bin')))
						open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
							adjustl(adjustr(NphRTfile))), status= 'replace', &
							form= 'unformatted', access= 'stream')
						write(expint) SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NphRTp(:, :, nn)
						close(expint)

					end do
				end if
			end if
		end if

		if ((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1)  .and. &
			(rank == 0)) then
			do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

				write(nnstring, '(I5)') nn
				write(sstring, '(I5)') s
				write(fstring, '(I5)') f
				write(Qindstring, '(I5)') Qind

				expstring= adjustl(adjustr(rankstring) // '_' // &
					adjustl(adjustr(nnstring) // '_' // &
					adjustl(adjustr(sstring) // '_' // &
					adjustl(adjustr(fstring) // '_' // &
					adjustl(adjustr(Qindstring) // '_')))))

				NphENARTfile= adjustl(adjustr(expstring) // &
					adjustl(adjustr('NphENARTfort.bin')))
				open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
					adjustl(adjustr(NphENARTfile))), &
					status= 'replace', form= 'unformatted', &
					access= 'stream')
				write(expint) SpecieT(s)%FluxTubeT(f)% &
					QCellT(Qind)%NphENARTp(:, :, :, nn)
				close(expint)

			end do
		end if

		! ----------------------------------------------------

		if ((SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1) == 1)  .and. &
			(rank == 0)) then
			if (SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT(1) == 1) then

				write(nnstring, '(I5)') nn
				write(sstring, '(I5)') s
				write(fstring, '(I5)') f

				expstring= adjustl(adjustr(rankstring) // '_' // &
					adjustl(adjustr(nnstring) // '_' // &
					adjustl(adjustr(sstring) // '_' // &
					adjustl(adjustr(fstring) // '_'))))

				M0phRTfile= adjustl(adjustr(expstring) // &
					adjustl(adjustr('M0phRTfort.bin')))
				open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
					adjustl(adjustr(M0phRTfile))), status= 'replace', &
					form= 'unformatted', access= 'stream')
				write(expint) SpecieT(s)%FluxTubeT(f)%M0phRT(nn, :)
				close(expint)

			end if
		end if

		! ----------------------------------------------------

		!if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
		!	do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
		!		if (rank == 0) then

		!			write(nnstring, '(I5)') nn
		!			write(sstring, '(I5)') s
		!			write(fstring, '(I5)') f

		!			expstring= adjustl(adjustr(rankstring) // '_' // &
		!				adjustl(adjustr(nnstring) // '_' // &
		!				adjustl(adjustr(sstring) // '_' // &
		!				adjustl(adjustr(fstring) // '_'))))

		!			M0phENARTfile= adjustl(adjustr(expstring) // &
		!				adjustl(adjustr('M0phENARTfort.bin')))
		!			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
		!				adjustl(adjustr(M0phENARTfile))), status= 'replace', &
		!				form= 'unformatted', access= 'stream')
		!			write(expint) SpecieT(s)%FluxTubeT(f)%M0phENART(nn, :)
		!			close(expint)

		!		end if
		!	end do
		!end if

		! ----------------------------------------------------

		if ((SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1) == 1)  .and. &
			(rank == 0)) then
			if (SpecieT(s)%FluxTubeT(f)%PHASEVELPERPIONMOMENTflagT(1) == 1) then
				if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

					write(nnstring, '(I5)') nn
					write(sstring, '(I5)') s
					write(fstring, '(I5)') f

					expstring= adjustl(adjustr(rankstring) // '_' // &
						adjustl(adjustr(nnstring) // '_' // &
						adjustl(adjustr(sstring) // '_' // &
						adjustl(adjustr(fstring) // '_'))))

					M1Perp1phRTfile= adjustl(adjustr(expstring) // &
						adjustl(adjustr('M1Perp1phRTfort.bin')))
					open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
						adjustl(adjustr(M1Perp1phRTfile))), status= 'replace', &
						form= 'unformatted', access= 'stream')
					write(expint) SpecieT(s)%FluxTubeT(f)%M1Perp1phRT(nn, :)
					close(expint)

					M1Perp2phRTfile= adjustl(adjustr(expstring) // &
						adjustl(adjustr('M1Perp2phRTfort.bin')))
					open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
						adjustl(adjustr(M1Perp2phRTfile))), status= 'replace', &
						form= 'unformatted', access= 'stream')
					write(expint) SpecieT(s)%FluxTubeT(f)%M1Perp2phRT(nn, :)
					close(expint)

				end if
			end if
		end if

		if ((SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1) == 1)  .and. &
			(rank == 0)) then
			if (SpecieT(s)%FluxTubeT(f)%PHASEVELPERPIONMOMENTflagT(1) == 1) then
				if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 0) then

					write(nnstring, '(I5)') nn
					write(sstring, '(I5)') s
					write(fstring, '(I5)') f

					expstring= adjustl(adjustr(rankstring) // '_' // &
						adjustl(adjustr(nnstring) // '_' // &
						adjustl(adjustr(sstring) // '_' // &
						adjustl(adjustr(fstring) // '_'))))

					M1PerpphRTfile= adjustl(adjustr(expstring) // &
						adjustl(adjustr('M1PerpphRTfort.bin')))
					open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
						adjustl(adjustr(M1PerpphRTfile))), status= 'replace', &
						form= 'unformatted', access= 'stream')
					write(expint) SpecieT(s)%FluxTubeT(f)%M1PerpphRT(nn, :)
					close(expint)

				end if
			end if
		end if

		! ----------------------------------------------------

		if ((SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1) == 1)  .and. &
			(rank == 0)) then
			if (SpecieT(s)%FluxTubeT(f)%PHASEVELPARIONMOMENTflagT(1) == 1) then

				write(nnstring, '(I5)') nn
				write(sstring, '(I5)') s
				write(fstring, '(I5)') f

				expstring= adjustl(adjustr(rankstring) // '_' // &
					adjustl(adjustr(nnstring) // '_' // &
					adjustl(adjustr(sstring) // '_' // &
					adjustl(adjustr(fstring) // '_'))))

				M1ParphRTfile= adjustl(adjustr(expstring) // &
					adjustl(adjustr('M1ParphRTfort.bin')))
				open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
					adjustl(adjustr(M1ParphRTfile))), status= 'replace', &
					form= 'unformatted', access= 'stream')
				write(expint) SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn, :)
				close(expint)

			end if
		end if

		if ((SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1) == 1)  .and. &
			(rank == 0)) then
			if (SpecieT(1)%FluxTubeT(1)%MOMENTFILTERflagT(1) == 1) then

				write(nnstring, '(I5)') nn
				write(sstring, '(I5)') s
				write(fstring, '(I5)') f

				expstring= adjustl(adjustr(rankstring) // '_' // &
					adjustl(adjustr(nnstring) // '_' // &
					adjustl(adjustr(sstring) // '_' // &
					adjustl(adjustr(fstring) // '_'))))

				M0FiltAvrgRTfile= adjustl(adjustr(expstring) // &
					adjustl(adjustr('M0FiltAvrgRTfort.bin')))
				open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
					adjustl(adjustr(M0FiltAvrgRTfile))), status= 'replace', &
					form= 'unformatted', access= 'stream')
				write(expint) SpecieT(s)%FluxTubeT(f)%M0FiltAvrgRT(nn, :)
				close(expint)

				M1Perp1FiltAvrgRTfile= adjustl(adjustr(expstring) // &
					adjustl(adjustr('M1Perp1FiltAvrgRTfort.bin')))
				open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
					adjustl(adjustr(M1Perp1FiltAvrgRTfile))), status= 'replace', &
					form= 'unformatted', access= 'stream')
				write(expint) SpecieT(s)%FluxTubeT(f)%M1Perp1FiltAvrgRT(nn, :)
				close(expint)

				M1Perp2FiltAvrgRTfile= adjustl(adjustr(expstring) // &
					adjustl(adjustr('M1Perp2FiltAvrgRTfort.bin')))
				open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
					adjustl(adjustr(M1Perp2FiltAvrgRTfile))), status= 'replace', &
					form= 'unformatted', access= 'stream')
				write(expint) SpecieT(s)%FluxTubeT(f)%M1Perp2FiltAvrgRT(nn, :)
				close(expint)

				M1ParFiltAvrgRTfile= adjustl(adjustr(expstring) // &
					adjustl(adjustr('M1ParFiltAvrgRTfort.bin')))
				open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
					adjustl(adjustr(M1ParFiltAvrgRTfile))), status= 'replace', &
					form= 'unformatted', access= 'stream')
				write(expint) SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(nn, :)
				close(expint)

				M2ParFiltAvrgRTfile= adjustl(adjustr(expstring) // &
					adjustl(adjustr('M2ParFiltAvrgRTfort.bin')))
				open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
					adjustl(adjustr(M2ParFiltAvrgRTfile))), status= 'replace', &
					form= 'unformatted', access= 'stream')
				write(expint) SpecieT(s)%FluxTubeT(f)%M2ParFiltAvrgRT(nn, :)
				close(expint)

			end if
		end if

		! ----------------------------------------------------

		if ((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. (rank == 0)) then

			write(nnstring, '(I5)') nn
			write(sstring, '(I5)') s
			write(fstring, '(I5)') f

			expstring= adjustl(adjustr(rankstring) // '_' // &
				adjustl(adjustr(nnstring) // '_' // &
				adjustl(adjustr(sstring) // '_' // &
				adjustl(adjustr(fstring) // '_'))))

			sigmaIonNeutRTfile= adjustl(adjustr(expstring) // &
				adjustl(adjustr('sigmaIonNeutRTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
				adjustl(adjustr(sigmaIonNeutRTfile))), status= 'replace', &
				form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%sigmaIonNeutRT(nn, :)
			close(expint)

			nuIonNeutRTfile= adjustl(adjustr(expstring) // &
				adjustl(adjustr('nuIonNeutRTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
				adjustl(adjustr(nuIonNeutRTfile))), status= 'replace', &
				form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, :)
			close(expint)

		end if

		! ----------------------------------------------------

		if ((SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1) == 1)  .and. &
			(rank == 0)) then
			if (SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(1) == 1) then

				write(nnstring, '(I5)') nn
				write(sstring, '(I5)') s
				write(fstring, '(I5)') f

				expstring= adjustl(adjustr(rankstring) // '_' // &
					adjustl(adjustr(nnstring) // '_' // &
					adjustl(adjustr(sstring) // '_' // &
					adjustl(adjustr(fstring) // '_'))))

				M2phRTfile= adjustl(adjustr(expstring) // &
					adjustl(adjustr('M2phRTfort.bin')))
				open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
					adjustl(adjustr(M2phRTfile))), status= 'replace', &
					form= 'unformatted', access= 'stream')
				write(expint) SpecieT(s)%FluxTubeT(f)%M2phRT(nn, :)
				close(expint)

				if (SpecieT(s)%FluxTubeT(f)% &
					PHASEENERGYPERPIONMOMENTflagT(1) == 1) then
					if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

						M2Perp1phRTfile= adjustl(adjustr(expstring) // &
							adjustl(adjustr('M2Perp1phRTfort.bin')))
						open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
							adjustl(adjustr(M2Perp1phRTfile))), status= 'replace', &
							form= 'unformatted', access= 'stream')
						write(expint) SpecieT(s)%FluxTubeT(f)%M2Perp1phRT(nn, :)
						close(expint)

						M2Perp2phRTfile= adjustl(adjustr(expstring) // &
							adjustl(adjustr('M2Perp2phRTfort.bin')))
						open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
							adjustl(adjustr(M2Perp2phRTfile))), status= 'replace', &
							form= 'unformatted', access= 'stream')
						write(expint) SpecieT(s)%FluxTubeT(f)%M2Perp2phRT(nn, :)
						close(expint)

					else

						M2PerpphRTfile= adjustl(adjustr(expstring) // &
							adjustl(adjustr('M2PerpphRTfort.bin')))
						open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
							adjustl(adjustr(M2PerpphRTfile))), status= 'replace', &
							form= 'unformatted', access= 'stream')
						write(expint) SpecieT(s)%FluxTubeT(f)%M2PerpphRT(nn, :)
						close(expint)

					end if
				end if

				if (SpecieT(s)%FluxTubeT(f)% &
					PHASEENERGYPARIONMOMENTflagT(1) == 1) then

					M2ParphRTfile= adjustl(adjustr(expstring) // &
						adjustl(adjustr('M2ParphRTfort.bin')))
					open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
						adjustl(adjustr(M2ParphRTfile))), status= 'replace', &
						form= 'unformatted', access= 'stream')
					write(expint) SpecieT(s)%FluxTubeT(f)%M2ParphRT(nn, :)
					close(expint)

				end if

			end if
		end if

		! ----------------------------------------------------

		!if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
		!	do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
		!		if (rank == 0) then

		!			write(nnstring, '(I5)') nn
		!			write(sstring, '(I5)') s
		!			write(fstring, '(I5)') f

		!			expstring= adjustl(adjustr(rankstring) // '_' // &
		!				adjustl(adjustr(nnstring) // '_' // &
		!				adjustl(adjustr(sstring) // '_' // &
		!				adjustl(adjustr(fstring) // '_'))))

		!			M1PphENARTfile= adjustl(adjustr(expstring) // &
		!				adjustl(adjustr('M1PphENARTfort.bin')))
		!			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
		!				adjustl(adjustr(M1PphENARTfile))), status= 'replace', &
		!				form= 'unformatted', access= 'stream')
		!			write(expint) SpecieT(s)%FluxTubeT(f)%M1PphENART(nn, :)
		!			close(expint)

		!		end if
		!	end do
		!end if

		! ----------------------------------------------------

		!if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
		!	do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
		!		if (rank == 0) then

		!			write(nnstring, '(I5)') nn
		!			write(sstring, '(I5)') s
		!			write(fstring, '(I5)') f

		!			expstring= adjustl(adjustr(rankstring) // '_' // &
		!				adjustl(adjustr(nnstring) // '_' // &
		!				adjustl(adjustr(sstring) // '_' // &
		!				adjustl(adjustr(fstring) // '_'))))

		!			M1QphENARTfile= adjustl(adjustr(expstring) // &
		!				adjustl(adjustr('M1QphENARTfort.bin')))
		!			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
		!				adjustl(adjustr(M1QphENARTfile))), status= 'replace', &
		!				form= 'unformatted', access= 'stream')
		!			write(expint) SpecieT(s)%FluxTubeT(f)%M1QphENART(nn, :)
		!			close(expint)

		!		end if
		!	end do
		!end if

		! ----------------------------------------------------

		!if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
		!	do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
		!		if (rank == 0) then

		!			write(nnstring, '(I5)') nn
		!			write(sstring, '(I5)') s
		!			write(fstring, '(I5)') f

		!			expstring= adjustl(adjustr(rankstring) // '_' // &
		!				adjustl(adjustr(nnstring) // '_' // &
		!				adjustl(adjustr(sstring) // '_' // &
		!				adjustl(adjustr(fstring) // '_'))))

		!			M1PHIphENARTfile= adjustl(adjustr(expstring) // &
		!				adjustl(adjustr('M1PHIphENARTfort.bin')))
		!			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
		!				adjustl(adjustr(M1PHIphENARTfile))), status= 'replace', &
		!				form= 'unformatted', access= 'stream')
		!			write(expint) SpecieT(s)%FluxTubeT(f)%M1PHIphENART(nn, :)
		!			close(expint)

		!		end if
		!	end do
		!end if

		! ----------------------------------------------------

		!if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
		!	do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
		!		if (rank == 0) then

		!			write(nnstring, '(I5)') nn
		!			write(sstring, '(I5)') s
		!			write(fstring, '(I5)') f

		!			expstring= adjustl(adjustr(rankstring) // '_' // &
		!				adjustl(adjustr(nnstring) // '_' // &
		!				adjustl(adjustr(sstring) // '_' // &
		!				adjustl(adjustr(fstring) // '_'))))

		!			M2phENARTfile= adjustl(adjustr(expstring) // &
		!				adjustl(adjustr('M2phENARTfort.bin')))
		!			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
		!				adjustl(adjustr(M2phENARTfile))), status= 'replace', &
		!				form= 'unformatted', access= 'stream')
		!			write(expint) SpecieT(s)%FluxTubeT(f)%M2phENART(nn, :)
		!			close(expint)

		!			M2PphENARTfile= adjustl(adjustr(expstring) // &
		!				adjustl(adjustr('M2PphENARTfort.bin')))
		!			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
		!				adjustl(adjustr(M2PphENARTfile))), status= 'replace', &
		!				form= 'unformatted', access= 'stream')
		!			write(expint) SpecieT(s)%FluxTubeT(f)%M2PphENART(nn, :)
		!			close(expint)

		!			M2QphENARTfile= adjustl(adjustr(expstring) // &
		!				adjustl(adjustr('M2QphENARTfort.bin')))
		!			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
		!				adjustl(adjustr(M2QphENARTfile))), status= 'replace', &
		!				form= 'unformatted', access= 'stream')
		!			write(expint) SpecieT(s)%FluxTubeT(f)%M2QphENART(nn, :)
		!			close(expint)

		!			M2PHIphENARTfile= adjustl(adjustr(expstring) // &
		!				adjustl(adjustr('M2PHIphENARTfort.bin')))
		!			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
		!				adjustl(adjustr(M2PHIphENARTfile))), status= 'replace', &
		!				form= 'unformatted', access= 'stream')
		!			write(expint) SpecieT(s)%FluxTubeT(f)%M2PHIphENART(nn, :)
		!			close(expint)

		!		end if
		!	end do
		!end if

		! ----------------------------------------------------

		if ((SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1) == 1) .and. &
			(rank == 0)) then

			write(nnstring, '(I5)') nn
			write(sstring, '(I5)') s
			write(fstring, '(I5)') f

			expstring= adjustl(adjustr(rankstring) // '_' // &
				adjustl(adjustr(nnstring) // '_' // &
				adjustl(adjustr(sstring) // '_' // &
				adjustl(adjustr(fstring) // '_'))))

			EAInertialRTfile= adjustl(adjustr(expstring) // &
				adjustl(adjustr('EAInertialRTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
				adjustl(adjustr(EAInertialRTfile))), status= 'replace', &
				form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%EAInertialRT(nn, :)
			close(expint)

			EAPressureRTfile= adjustl(adjustr(expstring) // &
				adjustl(adjustr('EAPressureRTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
				adjustl(adjustr(EAPressureRTfile))), status= 'replace', &
				form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%EAPressureRT(nn, :)
			close(expint)

			EAmagRTfile= adjustl(adjustr(expstring) // &
				adjustl(adjustr('EAmagRTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
				adjustl(adjustr(EAmagRTfile))), status= 'replace', &
				form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%EAmagRT(nn, :)
			close(expint)

			EGmagRTfile= adjustl(adjustr(expstring) // &
				adjustl(adjustr('EGmagRTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
				adjustl(adjustr(EGmagRTfile))), status= 'replace', &
				form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%EGmagRT(nn, :)
			close(expint)

			EPmagRTfile= adjustl(adjustr(expstring) // &
				adjustl(adjustr('EPmagRTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
				adjustl(adjustr(EPmagRTfile))), status= 'replace', &
				form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, :)
			close(expint)

			PhiParRTfile= adjustl(adjustr(expstring) // &
				adjustl(adjustr('PhiParRTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
				adjustl(adjustr(PhiParRTfile))), status= 'replace', &
				form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%PhiParRT(nn, :)
			close(expint)

		end if

		if ((SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1) == 1) .and. &
			(rank == 0)) then

			write(sstring, '(I5)') s
			write(fstring, '(I5)') f

			expstring= adjustl(adjustr(sstring) // '_' // &
				adjustl(adjustr(fstring) // '_'))

			NqLBTfile= adjustl(adjustr(expstring) // &
				adjustl(adjustr('NqLBTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(Densitydatadir) // &
				adjustl(adjustr(NqLBTfile))), status= 'replace', &
				form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%NqLBT(1)
			close(expint)

			NqUBTfile= adjustl(adjustr(expstring) // &
				adjustl(adjustr('NqUBTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(Densitydatadir) // &
				adjustl(adjustr(NqUBTfile))), status= 'replace', &
				form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%NqUBT(1)
			close(expint)

			DensityOutputRTfile= adjustl(adjustr(expstring) // &
				adjustl(adjustr('DensityOutputRTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(Densitydatadir) // &
				adjustl(adjustr(DensityOutputRTfile))), status= 'replace', &
				form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%DensityOutputRT(:)
			close(expint)

			TemperatureOutputRTfile= adjustl(adjustr(expstring) // &
				adjustl(adjustr('TemperatureOutputRTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(Densitydatadir) // &
				adjustl(adjustr(TemperatureOutputRTfile))), status= 'replace', &
				form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%TemperatureOutputRT(:)
			close(expint)

			EAInertialOutputRTfile= adjustl(adjustr(expstring) // &
				adjustl(adjustr('EAInertialOutputRTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(Densitydatadir) // &
				adjustl(adjustr(EAInertialOutputRTfile))), status= 'replace', &
				form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%EAInertialOutputRT(:)
			close(expint)

			EAPressureOutputRTfile= adjustl(adjustr(expstring) // &
				adjustl(adjustr('EAPressureOutputRTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(Densitydatadir) // &
				adjustl(adjustr(EAPressureOutputRTfile))), status= 'replace', &
				form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%EAPressureOutputRT(:)
			close(expint)

			EAmagOutputRTfile= adjustl(adjustr(expstring) // &
				adjustl(adjustr('EAmagOutputRTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(Densitydatadir) // &
				adjustl(adjustr(EAmagOutputRTfile))), status= 'replace', &
				form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%EAmagOutputRT(:)
			close(expint)

			nsnormCLBGTfile= adjustl(adjustr(expstring) // &
				adjustl(adjustr('nsnormCLBOutputfort.bin')))
			open(unit= expint, file= adjustl(adjustr(Densitydatadir) // &
				adjustl(adjustr(nsnormCLBGTfile))), status= 'replace', &
				form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%nsnormCLBGT(1)
			close(expint)

			nsnormCUBGTfile= adjustl(adjustr(expstring) // &
				adjustl(adjustr('nsnormCUBOutputfort.bin')))
			open(unit= expint, file= adjustl(adjustr(Densitydatadir) // &
				adjustl(adjustr(nsnormCUBGTfile))), status= 'replace', &
				form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%nsnormCUBGT(1)
			close(expint)

			LBREPLENISHflagTfile= adjustl(adjustr(expstring) // &
				adjustl(adjustr('LBREPLENISHflagTOutputfort.bin')))
			open(unit= expint, file= adjustl(adjustr(Densitydatadir) // &
				adjustl(adjustr(LBREPLENISHflagTfile))), status= 'replace', &
				form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1)
			close(expint)

			UBREPLENISHflagTfile= adjustl(adjustr(expstring) // &
				adjustl(adjustr('UBREPLENISHflagTOutputfort.bin')))
			open(unit= expint, file= adjustl(adjustr(Densitydatadir) // &
				adjustl(adjustr(UBREPLENISHflagTfile))), status= 'replace', &
				form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1)
			close(expint)

		end if

		! ----------------------------------------------------
		! ----------------------------------------------------

		if (((n == 1) .and. (nn == 1)) .and. &
			(n /= sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1))) .and. &
			((n /= SpecieT(s)%FluxTubeT(f)%NtT(1)) .and. &
			(nn /= SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))) then
			if (SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1) == 1) then
				if (rank == 0) then
					call cpu_time(SDE2End)
					write(nnstring, '(I5)') nn
					write(sstring, '(I5)') s
					write(fstring, '(I5)') f
					write(Timestring, '(I5)') nint(Time(1))
					write(SDE2string, '(i10)')  nint(SDE2End)
					write(*, *) trim('-- INITIAL DATA EXPORT COMPLETE: STATISTICAL TIME-STEP= ' &
						// adjustl(nnstring)) // &
						trim(', RANK= ' // adjustl(rankstring)) // &
						trim(', PARTICLE SPECIE= ' // adjustl(sstring)) // &
						trim(', FLUX-TUBE= ' // adjustl(fstring)) // &
						trim(', SIM-TIME= ' // adjustl(Timestring)) // &
						trim(', REAL-TIME= ' // adjustl(SDE2string)) // &
						trim(' s.')
				end if
			end if
		end if

		if (((n /= 1) .and. (nn /= 1)) .and. &
			(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1))) .and. &
			((n /= SpecieT(s)%FluxTubeT(f)%NtT(1)) .and. &
			(nn /= SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))) then
			if (SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1) == 1) then
				if (rank == 0) then
					call cpu_time(SDE2End)
					write(nnstring, '(I5)') nn
					write(sstring, '(I5)') s
					write(fstring, '(I5)') f
					write(Timestring, '(I5)') nint(Time(1))
					write(SDE2string, '(i10)')  nint(SDE2End)
					write(*, *) trim('-- DATA EXPORT COMPLETE: STATISTICAL TIME-STEP= ' &
						// adjustl(nnstring)) // &
						trim(', RANK= ' // adjustl(rankstring)) // &
						trim(', PARTICLE SPECIE= ' // adjustl(sstring)) // &
						trim(', FLUX-TUBE= ' // adjustl(fstring)) // &
						trim(', SIM-TIME= ' // adjustl(Timestring)) // &
						trim(' s., REAL-TIME= ' // adjustl(SDE2string)) // &
						trim(' s.')
				end if
			end if
		end if

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine DataExport2Sub

end module DataExport2
