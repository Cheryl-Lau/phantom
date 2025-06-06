!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module readwrite_dumps_common
!
! readwrite_dumps_common
!
! :References: None
!
! :Owner: Daniel Mentiplay
!
! :Runtime parameters: None
!
! :Dependencies: dim, dump_utils, eos, gitinfo, io, options, part,
!   sphNGutils
!
 use dump_utils, only:lenid
 implicit none

contains

!--------------------------------------------------------------------
!+
!  contruct header string based on compile-time options
!  these are for information only (ie. not important for restarting)
!+
!--------------------------------------------------------------------
character(len=lenid) function fileident(firstchar,codestring)
 use part,    only:h2chemistry,mhd,npartoftype,idust,gravity,lightcurve
 use options, only:use_dustfrac
 use dim,     only:use_dustgrowth,phantom_version_string,use_krome,store_dust_temperature
 use gitinfo, only:gitsha
 character(len=2), intent(in) :: firstchar
 character(len=*), intent(in), optional :: codestring
 character(len=10) :: datestring, timestring
 character(len=30) :: string
!
!--print date and time stamp in file header
!
 call date_and_time(datestring,timestring)
 datestring = datestring(7:8)//'/'//datestring(5:6)//'/'//datestring(1:4)
 timestring = timestring(1:2)//':'//timestring(3:4)//':'//timestring(5:)

 string = ' '
 if (gravity) string = trim(string)//'+grav'
 if (npartoftype(idust) > 0) string = trim(string)//'+dust'
 if (use_dustfrac) string = trim(string)//'+1dust'
 if (h2chemistry) string = trim(string)//'+H2chem'
 if (lightcurve) string = trim(string)//'+lightcurve'
 if (use_dustgrowth) string = trim(string)//'+dustgrowth'
 if (use_krome) string = trim(string)//'+krome'
 if (store_dust_temperature) string = trim(string)//'+Tdust'
#ifdef NUCLEATION
 string = trim(string)//'+nucleation'
#endif
 if (present(codestring)) then
    fileident = firstchar//':'//trim(codestring)//':'//trim(phantom_version_string)//':'//gitsha
 else
    fileident = firstchar//':Phantom'//':'//trim(phantom_version_string)//':'//gitsha
 endif

 if (mhd) then
    fileident = trim(fileident)//' (mhd+clean'//trim(string)//')  : '//trim(datestring)//' '//trim(timestring)
 else
    fileident = trim(fileident)//' (hydro'//trim(string)//'): '//trim(datestring)//' '//trim(timestring)
 endif

end function fileident

!--------------------------------------------------------------------
!+
!  extract various options used in Phantom from the fileid string
!+
!--------------------------------------------------------------------
subroutine get_options_from_fileid(fileid,tagged,phantomdump,smalldump,&
                                   use_onefluiddust,ierr)
 character(len=lenid), intent(in)  :: fileid
 logical,              intent(out) :: tagged,phantomdump,smalldump,use_onefluiddust
 integer,              intent(out) :: ierr
!
!--if file is a small dump, return an error code but still read what
!  can be read from a small dump
!
 ierr = 0
 tagged      = .false.
 smalldump   = .false.
 phantomdump = .false.
 if (fileid(2:2)=='T') tagged = .true.
 if (fileid(1:1) /= 'F') then
    !write(*,*) 'ERROR! file header indicates file is not a full dump'
    ierr = 1
    if (fileid(1:1)=='S') smalldump = .true.
 endif
 if (index(fileid,'Phantom') /= 0) then
    phantomdump = .true.
 elseif (index(fileid,'sphNG') /= 0) then
    phantomdump = .false.
    write(*,*) 'reading dump in sphNG format'
 else
    write(*,*) 'WARNING: could not determine Phantom/sphNG from fileident'
    write(*,*) '(assuming sphNG...)'
    phantomdump = .false.
 endif
 if (index(fileid,'+1dust') /= 0) then
    use_onefluiddust = .true.
 else
    use_onefluiddust = .false.
 endif

end subroutine get_options_from_fileid

!---------------------------------------------------------------
!+
!  make sure required arrays have been read from Phantom file
!  and perform basic sanity checks
!+
!---------------------------------------------------------------
subroutine check_arrays(i1,i2,npartoftype,npartread,nptmass,nsinkproperties,massoftype,&
                        alphafile,tfile,phantomdump,got_iphase,got_xyzh,got_vxyzu,got_alpha, &
                        got_krome_mols,got_krome_gamma,got_krome_mu,got_krome_T, &
                        got_abund,got_dustfrac,got_sink_data,got_sink_vels,got_Bxyz,got_psi,got_dustprop,got_pxyzu,got_VrelVf, &
                        got_dustgasprop,got_temp,got_raden,got_kappa,got_Tdust,got_iorig,iphase,&
                        xyzh,vxyzu,pxyzu,alphaind,xyzmh_ptmass,Bevol,iorig,iprint,ierr)
 use dim,  only:maxp,maxvxyzu,maxalpha,maxBevol,mhd,h2chemistry,store_temperature,&
                use_dustgrowth,gr,do_radiation,store_dust_temperature
 use eos,  only:polyk,gamma
 use part, only:maxphase,isetphase,set_particle_type,igas,ihacc,ihsoft,imacc,&
                xyzmh_ptmass_label,vxyz_ptmass_label,get_pmass,rhoh,dustfrac,ndusttypes,norig
 use io,   only:warning,id,master
 use options,    only:alpha,use_dustfrac
 use sphNGutils, only:itype_from_sphNG_iphase,isphNG_accreted
 integer,         intent(in)    :: i1,i2,npartoftype(:),npartread,nptmass,nsinkproperties
 real,            intent(in)    :: massoftype(:),alphafile,tfile
 logical,         intent(in)    :: phantomdump,got_iphase,got_xyzh(:),got_vxyzu(:),got_alpha,got_dustprop(:)
 logical,         intent(in)    :: got_VrelVf,got_dustgasprop(:)
 logical,         intent(in)    :: got_abund(:),got_dustfrac(:),got_sink_data(:),got_sink_vels(:),got_Bxyz(:)
 logical,         intent(in)    :: got_krome_mols(:),got_krome_gamma,got_krome_mu,got_krome_T
 logical,         intent(in)    :: got_psi,got_temp,got_Tdust,got_pxyzu(:),got_raden(:),got_kappa,got_iorig
 integer(kind=1), intent(inout) :: iphase(:)
 integer(kind=8), intent(inout) :: iorig(:)
 real,            intent(inout) :: vxyzu(:,:),Bevol(:,:),pxyzu(:,:)
 real(kind=4),    intent(inout) :: alphaind(:,:)
 real,            intent(inout) :: xyzh(:,:),xyzmh_ptmass(:,:)
 integer,         intent(in)    :: iprint
 integer,         intent(out)   :: ierr
 logical :: use_gas
 integer :: i,itype,nread
 !
 ! particle type information
 !
 if (maxphase==maxp) then
    if (got_iphase) then
       if (phantomdump) then
          do i=i1,i2
             itype = int(iphase(i))
             iphase(i) = isetphase(itype,iactive=.true.)
          enddo
       else
          ! convert from sphNG
          do i=i1,i2
             itype = itype_from_sphNG_iphase(iphase(i))
             iphase(i) = isetphase(itype,iactive=.true.)
             if (itype==isphNG_accreted) then
                !  mark accreted/unknown particle types as dead according
                !  to Phantom (i.e., give negative smoothing length)
                xyzh(4,i) = -abs(xyzh(4,i))
                call set_particle_type(i,igas) ! to give an allowed particle type
             endif
          enddo
       endif
    elseif (any(npartoftype(2:) > 0)) then
       !
       !--iphase is required if there is more than one particle type
       !
       write(*,*) 'error in rdump: need type information but iamtype not present in dump file'
       ierr = 8
       return
    else
       !
       !--iphase does not need to be read if there is only one particle type
       !  but to start the code it should be set such that all particles are active
       !
       do i=i1,i2
          iphase(i) = isetphase(igas,iactive=.true.)
       enddo
    endif
 endif
 if (maxphase==maxp) then
    use_gas = .false.
 else
    use_gas = .true.
 endif

 !
 ! hydrodynamics arrays
 !
 if (any(.not.got_xyzh)) then
    if (id==master .and. i1==1) write(*,*) 'ERROR: x, y, z or h not found in file'
    ierr = 9
    return
 endif
 if (any(.not.got_vxyzu(1:3))) then
    if (id==master .and. i1==1) write(*,*) 'ERROR: missing velocity information from file'
 endif
 if (maxvxyzu==4 .and. .not.got_vxyzu(4)) then
    if (gamma < 1.01) then
       do i=i1,i2
          vxyzu(4,i) = 1.5*polyk
          !print*,'u = ',vxyzu(4,i)
       enddo
       if (id==master .and. i1==1) write(*,*) 'WARNING: u not in file but setting u = 3/2 * cs^2'
    else
       do i=i1,i2
          vxyzu(4,i) = (1.0/(gamma-1.0))*polyk*rhoh(xyzh(4,i),get_pmass(i,use_gas))**(gamma - 1.)
          !print*,'u = ',vxyzu(4,i)
       enddo
       if (id==master .and. i1==1) write(*,*) 'WARNING: u not in file but setting u = (K*rho**(gamma-1))/(gamma-1)'
    endif
 endif
 if (h2chemistry .and. .not.all(got_abund)) then
    if (id==master) write(*,*) 'error in rdump: using H2 chemistry, but abundances not found in dump file'
    ierr = 9
    return
 endif
#ifdef KROME
 if (.not.all(got_krome_mols).and. npartread > 0) then
    if (id==master) write(*,*) 'error in rdump: using KROME chemistry, but abundances not found in dump file'
!     ierr = 9
    return
 endif
 if (.not.got_krome_gamma .and. npartread > 0) then
    if (id==master) write(*,*) 'error in rdump: using KROME chemistry, but gamma not found in dump file'
!     ierr = 9
    return
 endif
 if (.not.got_krome_mu .and. npartread > 0) then
    if (id==master) write(*,*) 'error in rdump: using KROME chemistry, but mu not found in dump file'
!     ierr = 9
    return
 endif
 if (.not.got_krome_T .and. npartread > 0) then
    if (id==master) write(*,*) 'error in rdump: using KROME chemistry, but temperature not found in dump file'
!     ierr = 9
    return
 endif
#endif
 if (store_temperature .and. .not.got_temp) then
    if (id==master .and. i1==1) write(*,*) 'WARNING: missing temperature information from file'
 endif
 if (store_dust_temperature .and. .not.got_Tdust) then
    if (id==master .and. i1==1) write(*,*) 'WARNING: missing dust temperature information from file'
 endif
 if (maxalpha==maxp) then
    if (got_alpha) then
       if (alphafile < 0.99 .and. tfile > 0.) then
          if (any(alphaind(1,i1:i2) > 1.0 .or. alphaind(1,i1:i2) < 0.)) then
             if (id==master) write(iprint,*) 'ERROR! AV alpha < 0 or alpha > 1 in dump file: using alpha'
             alphaind(1,i1:i2) = real(alpha,kind=4)
          endif
       endif
    else
       if (id==master .and. i1==1) write(*,*) 'WARNING: alpha not found in file'
       alphaind(1,i1:i2) = real(alpha,kind=4)
    endif
 endif
 if (npartread > 0) then
    do i = 1, size(massoftype)
       if (npartoftype(i) > 0) then
          if (massoftype(i) <= 0.0) then
             if (id==master .and. i1==1) write(*,*) 'ERROR! mass not set in read_dump (Phantom)'
             ierr = 12
             return
          endif
       endif
    enddo
 endif
 if (use_dustfrac .and. .not. all(got_dustfrac(1:ndusttypes))) then
    if (id==master .and. i1==1) write(*,*) 'ERROR! using one-fluid dust, but no dust fraction found in dump file'
    if (id==master .and. i1==1) write(*,*) ' Setting dustfrac = 0'
    dustfrac = 0.
    !ierr = 13
    return
 endif
 if (use_dustgrowth .and. .not.got_dustprop(1)) then
    write(*,*) 'ERROR! using dustgrowth, but no grain size found in dump file'
    return
 endif
 if (use_dustgrowth .and. .not.got_dustprop(1)) then
    write(*,*) 'ERROR! using dustgrowth, but no grain density found in dump file'
    return
 endif
 if (use_dustgrowth .and. .not.got_VrelVf) then
    write(*,*) 'ERROR! using dustgrowth, but no Vrel/Vfrag found in dump file'
    return
 endif
 if (use_dustgrowth .and. .not.got_dustgasprop(3)) then
    write(*,*) 'ERROR! using dustgrowth, but no St found in dump file'
    return
 endif
 !
 ! sink particle arrays
 !
 nread = 0
 if (nptmass > 0) then
    do i=1,nsinkproperties
       if (.not.got_sink_data(i)) then
          if (i <= 5) then
             if (id==master) write(*,*) 'ERROR! sink particle '//trim(xyzmh_ptmass_label(i))//' not found'
             ierr = 10
             return
          else
             if (id==master) write(*,*) 'WARNING! sink particle '//trim(xyzmh_ptmass_label(i))//' not found'
          endif
       endif
    enddo
    if (.not.all(got_sink_vels(1:3))) then
       if (id==master .and. i1==1) write(*,*) 'WARNING! sink particle velocities not found'
    endif
    if (id==master .and. i1==1) then
       print "(2(a,i4),a)",' got ',nsinkproperties,' sink properties from ',nptmass,' sink particles'
       if (nptmass > 0) print "(1x,47('-'),/,1x,a,'|',4(a9,1x,'|'),/,1x,47('-'))",&
                              'ID',' Mass    ',' Racc    ',' Macc    ',' hsoft   '
       do i=1,min(nptmass,999)
          if (xyzmh_ptmass(4,i) > 0.) print "(i3,'|',4(1pg9.2,1x,'|'))",i,xyzmh_ptmass(4,i), &
                                            xyzmh_ptmass(ihacc,i),xyzmh_ptmass(imacc,i),xyzmh_ptmass(ihsoft,i)
       enddo
       if (nptmass > 0) print "(1x,47('-'))"
    endif
 endif
 !
 ! radiation arrays
 !
 if (do_radiation) then
    if (.not.all(got_raden)) then
       if (id==master .and. i1==1) write(*,*) 'ERROR: RADIATION=yes but radiation arrays not found in Phantom dump file'
       ierr = ierr + 1
    endif
    if (.not.got_kappa) then
       if (id==master .and. i1==1) write(*,*) 'WARNING: RADIATION=yes but opacity not found in Phantom dump file'
    endif
 endif

 !
 ! MHD arrays
 !
 if (mhd) then
    if (.not.all(got_Bxyz(1:3))) then
       if (id==master .and. i1==1) write(*,*) 'WARNING: MHD but magnetic field arrays not found in Phantom dump file'
    endif
    if (.not.got_psi) then
       if (id==master .and. i1==1) write(*,*) 'WARNING! div B cleaning field (Psi) not found in Phantom dump file: assuming psi=0'
       Bevol(maxBevol,i1:i2) = 0.
    endif
 endif

 !
 ! GR arrays
 !
 if (gr) then
    if (.not.all(got_pxyzu(1:3))) then
       write(*,*) 'WARNING: GR but momentum arrays not found in Phantom dump file'
       pxyzu(:,i1:i2) = 0.
    endif
 endif

!
! Particle IDs
!
 if (.not.got_iorig) then
    do i=i1,i2
       iorig(i) = i
    enddo
    norig = i2
    if (id==master .and. i1==1) write(*,*) 'WARNING: Particle IDs not in dump; resetting IDs'
 else
    norig = 0
    do i=i1,i2
       norig = max(norig,iorig(i))
    enddo
 endif

end subroutine check_arrays

end module readwrite_dumps_common
