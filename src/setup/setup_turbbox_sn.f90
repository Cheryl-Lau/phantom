!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setup
!
! Setup routine for a turbulent free-field case resembling a sub-grid model 
!
! :References: None
!
! :Owner: Cheryl Lau
!
! :Runtime parameters:
!   - dist_unit      : *distance unit (e.g. au)*
!   - mass_unit      : *mass unit (e.g. solarm)*
!   - np_req         : *requested number of particles*
!   - pmass          : *particle mass in code units*
!   - rhozero        : *density of box in code units*
!   - cs0            : *initial sound speed in code units*
!   - tmax           : *end time of simulation in code units*
!   - dtmax          : *maximum timestep of simulation in code units*
!
! :Dependencies: boundary, dim, domain, eos, infile_utils, io, options, part,
!   physcon, prompting, setup_params, timestep, unifdis, units
!
 use setup_params,  only:rhozero
 use timestep,      only:dtmax,tmax
 implicit none

 public  :: setpart

 private

 integer :: np_req
 real    :: pmass,cs0_cgs,rms_mach,rhozero_cgs 
 real(kind=8)       :: udist,umass
 character(len=20)  :: dist_unit,mass_unit

contains
!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:maxvxyzu
 use io,           only:master,fatal,warning
 use boundary,     only:xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound,set_boundary
 use part,         only:set_particle_type,igas
 use physcon,      only:pi,mass_proton_cgs,kboltz,years,pc,solarm,micron
 use units,        only:set_units,unit_density,unit_velocity,unit_ergg,utime,select_unit
 use eos,          only:gmw,ieos
 use part,         only:periodic
 use unifdis,      only:set_unifdis
 use options,      only:alphau,nfulldump,nmaxdumps,icooling,ipdv_heating,ishock_heating
 use cooling,      only:ufloor,Tfloor
 use timestep,     only:nout
 use prompting,    only:prompt
 use velfield,     only:set_velfield_from_cubes
 use datafiles,    only:find_phantom_datafile
 use inject,       only:inject_sn,sink_progenitor,frackin,fractherm 
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma
 real,              intent(inout) :: hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 character(len=100) :: filename,temp_ans,dens_ans
 character(len=20), parameter     :: filevx = 'cube_v1.dat'
 character(len=20), parameter     :: filevy = 'cube_v2.dat'
 character(len=20), parameter     :: filevz = 'cube_v3.dat'
 integer, parameter :: maxrow = 1E7
 real    :: xyzh_raw(4,maxrow)
 real    :: vol_box,boxlength,boxsize_fromfile,boxsize_toscale,deltax
 real    :: boxsize_sample,rho_sample,rho_sample_cgs,rho_fracdiff
 real    :: totmass,totmass_req,temp,pmass_cgs,u0,u0_cgs,tmax_cgs,dtmax_cgs
real     :: xmini,xmaxi,ymini,ymaxi,zmini,zmaxi
 real    :: t_ff,rmsmach,v2i,turbfac,turbboxsize,mean_v,sigma_v,cs0
 integer :: i,ierr,io_file,ix,npmax,npart_sample
 logical :: iexist
 character(len=120) :: filex,filey,filez

 !
 !--general parameters
 !
 time = 0.
 gamma = 5./3.  ! adiabatic eos index

 !
 ! Interactive setup
 !
 filename=trim(fileprefix)//'.setup'
 inquire(file=filename,exist=iexist)
 if (iexist) then
    !--read from setup file
    call read_setupfile(filename,ierr)
    if (ierr /= 0) then
       if (id==master) call write_setupfile(filename)
       stop
    endif
 elseif (id==master) then
    mass_unit = 'solarm'
    dist_unit = 'pc'
    ierr = 1
    do while (ierr /= 0)
       call prompt('Enter mass unit (e.g. solarm,jupiterm,earthm)',mass_unit)
       call select_unit(mass_unit,umass,ierr)
       if (ierr /= 0) print "(a)",' ERROR: mass unit not recognised'
    enddo
    ierr = 1
    do while (ierr /= 0)
       call prompt('Enter distance unit (e.g. au,pc,kpc,0.1pc)',dist_unit)
       call select_unit(dist_unit,udist,ierr)
       if (ierr /= 0) print "(a)",' ERROR: length unit not recognised'
    enddo
    call set_units(dist=udist,mass=umass,G=1.d0)

    !
    ! Set number of particles (to be updated after set_unifdis)
    !
    npmax = int(size(xyzh(1,:)))
    np_req = 1E6
    call prompt('Enter total number of particles',np_req,1)
    if (np_req > npmax) call fatal('setup_unifdis_cmi','number of particles exceeded limit')
    !
    ! Particle mass
    !
    pmass_cgs = 1.9891E30 ! 10^-3 solarm
    pmass = pmass_cgs/umass
    call prompt('Enter particle mass in units of '//mass_unit,pmass,0.)
    !
    ! Density
    !
    rhozero_cgs = 6E-21
    call prompt('Enter initial density in g cm^-3',rhozero_cgs,0.)

    !
    ! set initial sound speed
    !
    cs0_cgs = 0.91E5
    call prompt('Enter initial sound speed in cm s^-1',cs0_cgs,0.)
    !
    ! Set timestep and end-time
    !
    dtmax_cgs = 3.15360E9   ! 1E-4 Myr
    tmax_cgs  = 1e3*dtmax_cgs
    dtmax = dtmax_cgs/utime
    tmax  = tmax_cgs/utime
    call prompt('Enter timestep in code units',dtmax,0.)
    call prompt('Enter simulation end-time in code units',tmax,0.)
    !
    ! Set turbulence 
    !
    rms_mach = 10.
    call prompt('Enter the Mach number of the cloud turbulence',rms_mach,0.)

    print*,'*************'
    print*,rhozero_cgs,cs0_cgs 
    if (id==master) call write_setupfile(filename)
    stop 'rerun phantomsetup after editing .setup file'
 else
    stop
 endif
 !
 ! set units
 !
 call set_units(dist=udist,mass=umass,G=1.d0)

 !
 ! Convert units 
 !
 rhozero = rhozero_cgs/unit_density
 cs0 = cs0_cgs/unit_velocity 
 !
 ! Calculate total mass required with npart_requested
 !
 totmass_req = np_req*pmass
 !
 ! Calculate volume needed
 !
 vol_box   = totmass_req/rhozero
 boxlength = vol_box**(1./3.)
 xmini = -0.5*boxlength; xmaxi = 0.5*boxlength
 ymini = -0.5*boxlength; ymaxi = 0.5*boxlength
 zmini = -0.5*boxlength; zmaxi = 0.5*boxlength
!
! Set boundaries
!
xmini = -0.5*boxlength; xmaxi = 0.5*boxlength
ymini = -0.5*boxlength; ymaxi = 0.5*boxlength
zmini = -0.5*boxlength; zmaxi = 0.5*boxlength
call set_boundary(xmini,xmaxi,ymini,ymaxi,zmini,zmaxi)
 print*,'boundaries before: ',xmini,xmaxi,ymini,ymaxi,zmini,zmaxi

 !
 ! Set particle distribution
 !
 npart = 0
 deltax = (xmaxi-xmini)/np_req**(1./3.)
 call set_unifdis('cubic',id,master,xmini,xmaxi,ymini,ymaxi,zmini,zmaxi,deltax,hfact,&
                npart,xyzh,periodic)

 !
 ! Check real boundaries after setup
 !
 xmini = huge(xmini); xmaxi = tiny(xmaxi)
 ymini = huge(ymini); ymaxi = tiny(ymaxi)
 zmini = huge(zmini); zmaxi = tiny(zmaxi)
 do i = 1,npart
    xmini = min(xmini,xyzh(1,i)); xmaxi = max(xmaxi,xyzh(1,i))
    ymini = min(ymini,xyzh(2,i)); ymaxi = max(ymaxi,xyzh(2,i))
    zmini = min(zmini,xyzh(3,i)); zmaxi = max(zmaxi,xyzh(3,i))
 enddo 
 boxlength = max(xmaxi-xmini,ymaxi-ymini,zmaxi-zmini)
 print*,'boundaries after: ',xmini,xmaxi,ymini,ymaxi,zmini,zmaxi

 !
 ! Set particle properties
 !
 npartoftype(:) = 0
 npartoftype(igas) = npart
 massoftype(igas)  = pmass
 do i = 1,npartoftype(igas)
    call set_particle_type(i,igas)
 enddo


 !
 ! Check that the density of xyzh matches the input rhozero_cgs
 !
 boxsize_sample = 1.5*xmaxi
 npart_sample = 0
 do i = 1,npart
    if (xyzh(1,i) >= -0.5*boxsize_sample .and. xyzh(1,i) <= 0.5*boxsize_sample .and. &
        xyzh(2,i) >= -0.5*boxsize_sample .and. xyzh(2,i) <= 0.5*boxsize_sample .and. &
        xyzh(3,i) >= -0.5*boxsize_sample .and. xyzh(3,i) <= 0.5*boxsize_sample) then
       npart_sample = npart_sample + 1
    endif
 enddo
 if (npart_sample == 0) call fatal('setup_unifdis_cmi','no particles within sampling box')

 rho_sample     = npart_sample*massoftype(igas) / boxsize_sample**3.
 rho_sample_cgs = rho_sample*unit_density
 rhozero_cgs    = rhozero*unit_density   !- re-define since it's not written to .setup file
 print*,'required density:  ',rhozero_cgs,'g cm^-3'
 print*,'density from xyzh: ',rho_sample_cgs,'g cm^-3'

 rho_fracdiff = abs(rho_sample_cgs - rhozero_cgs)/rhozero_cgs
 if (rho_fracdiff > 1E-1) then
    call fatal('setup_unifdis_cmi','box density does not match the input value')
 elseif (rho_fracdiff > 1E-3) then
    write(*,'(a15,f3.1,a58)',advance='no') 'Box density is ',rho_fracdiff*100.,'% off from &
             &input value, would you like to proceed ([y]/n)?'
    read(*,'(a)') dens_ans
    if (len(trim(adjustl(dens_ans))) /= 0) then
       if (trim(adjustl(dens_ans)) == 'n') then
          print*,'stopping program - try adjust boxsize_sample'
          stop
       elseif (trim(adjustl(dens_ans)) == 'y') then
          print*,'y - proceeding...'
       else
          print*,'Invalid input'
       endif
    else
       print*,'y - proceeding...'
    endif
 endif
 
 !
 ! Calculate temperature from input sound speed
 !
 temp = (cs0*unit_velocity)**2*gmw*mass_proton_cgs/(gamma*kboltz)
 print*,'temperature: ',temp,'K'
 write(*,'(a)',advance='no') ' Temperature acceptable ([y]/n)?'
 read(*,'(a)') temp_ans
 if (len(trim(adjustl(temp_ans))) /= 0) then
    if (trim(adjustl(temp_ans)) == 'n') then
       print*,'stopping program - adjust gamma, gmw or cs0'
       stop
    elseif (trim(adjustl(temp_ans)) == 'y') then
       print*,'y - proceeding...'
    else
       print*,'Invalid input'
    endif
 else
    print*,'y - proceeding...'
 endif

 !
 ! Polytropic constant
 !
 if (maxvxyzu < 4 .or. gamma <= 1.) then
    polyk = cs0**2
 else
    polyk = 0.
 endif

 !
 ! Turbulent velocity field
 !
 if (rms_mach > 0.) then
    filex  = find_phantom_datafile(filevx,'velfield_sphng_small')
    filey  = find_phantom_datafile(filevy,'velfield_sphng_small')
    filez  = find_phantom_datafile(filevz,'velfield_sphng_small')
    !- To convert endian for different vfield files: setenv GFORTRAN_CONVERT_UNIT big/small_endian

    turbboxsize = boxlength
    call set_velfield_from_cubes(xyzh(:,1:npart),vxyzu(:,1:npart),npart, &
                                 filex,filey,filez,1.,turbboxsize,.false.,ierr)
    if (ierr /= 0) call fatal('setup','error setting up velocity field on clouds')

    rmsmach = 0.0
    print*, 'Turbulence being set by user'
    do i = 1,npart
       v2i     = mag(vxyzu(1:3,i))**2
       rmsmach = rmsmach + v2i/cs0**2
    enddo
    rmsmach = sqrt(rmsmach/npart)
    if (rmsmach > 0.) then
       turbfac = rms_mach/rmsmach ! normalise the energy to the desired mach number
     else
       turbfac = 0.
    endif

    !- Set velocities and get mean 
    open(2040,file='turbbox_velocities.dat',status='replace')
    mean_v = 0.
    do i = 1,npart
       vxyzu(1:3,i) = turbfac*vxyzu(1:3,i)
       write(2040,*) mag(vxyzu(1:3,i))*unit_velocity 
       mean_v = mean_v + mag(vxyzu(1:3,i))
    enddo
    close(2040)
    mean_v = mean_v/real(npart) 
    !- Get initial dispersion 
    sigma_v = 0. 
    do i = 1,npart
       sigma_v = sigma_v + (mag(vxyzu(1:3,i)) - mean_v)**2 
    enddo 
    sigma_v = sqrt(sigma_v/real(npart))
 endif

 !
 ! Set particle energies using temp
 !
 u0_cgs = kboltz * temp / (gmw*mass_proton_cgs*(gamma-1.))
 u0     = u0_cgs/unit_ergg
 do i = 1,npart
   if (size(vxyzu(:,i)) >= 4) vxyzu(4,i) = u0
 enddo


 !
 ! Set runtime parameters
 !
 ieos      = 2     ! adiabatic eos
 nout      = 1
 nmaxdumps = 100
 nfulldump = 1
 Tfloor    = 3.
 ufloor    = kboltz*Tfloor/(gmw*mass_proton_cgs*(gamma-1.))/unit_ergg
 icooling  = 7
 alphau    = 1.0
 ipdv_heating   = 1
 ishock_heating = 1

 !
 ! Set SN injection parameters 
 ! 
 inject_sn = .false.
 sink_progenitor = .false.
 frackin = 1.0
 fractherm = 0.0

 !- Print summary - for checking
 open(unit=2050,file='turbbox_properties.txt',status='replace')
 write(2050,*) ' -SUMMARY- '
 write(2050,*) 'totmass:    ',totmass,totmass*umass,'g'
 write(2050,*) 'rhozero:    ',rhozero,rhozero*unit_density,'g cm^-3'
 write(2050,*) 'box vol:    ',vol_box,vol_box*udist**3,'cm^3'
 write(2050,*) 'box length: ',boxlength,boxlength*udist,'cm'
 write(2050,*) 'temp:       ',temp,'K'
 write(2050,*) 'internal E: ',u0,u0_cgs,'erg g^-1'
 write(2050,*) 'cs0:        ',cs0,cs0*unit_velocity,'cm s^-1'
 if (rms_mach > 0.) then
    write(2050,*) 'rms_mach:   ',rms_mach
    write(2050,*) 'mean_v:     ',mean_v,mean_v*unit_velocity,'cm s^-1'
    write(2050,*) 'sigma_v:    ',sigma_v,sigma_v*unit_velocity,'cm s^-1'
 endif 
 write(2050,*) 'npart:      ',npart
 write(2050,*) 'part mass:  ',massoftype(igas),massoftype(igas)*umass,'g'
 close(2050) 

end subroutine setpart


real function mag(vec)
 real,   intent(in) :: vec(3)

 mag = sqrt(dot_product(vec,vec))

end function mag

!------------------------------------------------------------------------
!
! write setup file
!
!------------------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer :: iunit

 print "(/,a)",' writing setup options file '//trim(filename)
 open(newunit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for uniform setup routine'

 write(iunit,"(/,a)") '# units'
 call write_inopt(dist_unit,'dist_unit','distance unit (e.g. au)',iunit)
 call write_inopt(mass_unit,'mass_unit','mass unit (e.g. solarm)',iunit)
 !
 ! other parameters
 !
 write(iunit,"(/,a)") '# setup'
 call write_inopt(np_req,'np_req','total number of particles requested',iunit)
 call write_inopt(pmass,'pmass','particle mass',iunit)
 call write_inopt(rhozero_cgs,'rhozero_cgs','initial gas density',iunit)
 call write_inopt(cs0_cgs,'cs0_cgs','initial sound speed',iunit)
 call write_inopt(dtmax,'dtmax','timestep in code units',iunit)
 call write_inopt(tmax,'tmax','end-time in code units',iunit)
 call write_inopt(rms_mach,'rms_mach','turbulent rms mach number',iunit)

 close(iunit)

end subroutine write_setupfile

!------------------------------------------------------------------------
!
! read setup file
!
!------------------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use units,        only:select_unit
 use io,           only:error
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",' reading setup options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 !
 ! units
 !
 call read_inopt(mass_unit,'mass_unit',db,errcount=nerr)
 call read_inopt(dist_unit,'dist_unit',db,errcount=nerr)
 !
 ! other parameters
 !
 call read_inopt(np_req,'np_req',db,min=8,errcount=nerr)
 call read_inopt(pmass,'pmass',db,min=0.,errcount=nerr)
 call read_inopt(rhozero_cgs,'rhozero_cgs',db,min=0.,errcount=nerr)
 call read_inopt(cs0_cgs,'cs0_cgs',db,min=0.,errcount=nerr)
 call read_inopt(dtmax,'dtmax',db,min=0.,errcount=nerr)
 call read_inopt(tmax,'tmax',db,min=0.,errcount=nerr)
 call read_inopt(rms_mach,'rms_mach',db,nerr)

 call close_db(db)
 !
 ! parse units
 !
 call select_unit(mass_unit,umass,nerr)
 if (nerr /= 0) then
    call error('setup_unifdis','mass unit not recognised')
    ierr = ierr + 1
 endif
 call select_unit(dist_unit,udist,nerr)
 if (nerr /= 0) then
    call error('setup_unifdis','length unit not recognised')
    ierr = ierr + 1
 endif

 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
