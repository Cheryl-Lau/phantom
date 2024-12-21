!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setup
!
! Setup routine for uniform distribution in StarBench test
!- Particle positions can be done through importing a glass cube (specifying xyzh)
!- or calling setup_unifdis
!
! :References: None
!
! :Owner: Cheryl Lau
!
! :Runtime parameters:
!   - dist_unit      : *distance unit (e.g. au)*
!   - mass_unit      : *mass unit (e.g. solarm)*
!   - np_req         : *requested number of particles (if not using glass)*
!   - pmass          : *particle mass in code units*
!   - rhozero        : *density of box in code units*
!   - cs0            : *initial sound speed in code units*
!   - tmax           : *end time of simulation in code units*
!   - dtmax          : *maximum timestep of simulation in code units*
!   - glass_filename : *filename of glass cube to be imported*
!
! :Dependencies: boundary, dim, domain, eos, infile_utils, io, options, part,
!   physcon, prompting, setup_params, timestep, unifdis, units
!
 use setup_params,  only:rhozero
 use timestep,      only:dtmax,tmax
 use photoionize_cmi, only:tree_accuracy_cmi
 implicit none

 public  :: setpart

 private

 integer :: np_req
 real    :: pmass,cs0,rhozero_cgs 
 real    :: rcut_opennode,rcut_leafpart
 real(kind=8)       :: udist,umass
 character(len=20)  :: dist_unit,mass_unit
 character(len=100) :: glass_filename

 logical :: use_glass = .true.

contains
!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:maxvxyzu
 use io,           only:master,fatal,warning
 use boundary,     only:set_boundary
 use part,         only:set_particle_type,igas
 use physcon,      only:pi,mass_proton_cgs,kboltz,years,pc,solarm,micron
 use units,        only:set_units,unit_density,unit_velocity,unit_ergg,utime,select_unit
 use eos,          only:gmw,ieos
 use part,         only:periodic
 use unifdis,      only:set_unifdis
 use options,      only:alphau,nfulldump,icooling,ipdv_heating,ishock_heating
 use cooling,      only:ufloor,Tfloor
 use timestep,     only:nout
 use prompting,    only:prompt
 use photoionize_cmi, only:monochrom_source,fix_temp_hii,implicit_cmi,treat_Rtype_phase
 use photoionize_cmi, only:photoionize_tree,nHlimit_fac
 use photoionize_cmi, only:limit_voronoi,hlimit_fac,extradist_fac
 use photoionize_cmi, only:rcut_opennode_cgs,rcut_leafpart_cgs,delta_rcut_cgs
 use photoionize_cmi, only:nsetphotosrc,xyztq_setphotosrc_cgs
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
 integer, parameter :: maxrow = 1E7
 real    :: xyzh_raw(4,maxrow)
 real    :: vol_box,boxlength,boxsize_fromfile,boxsize_toscale,deltax
 real    :: boxsize_setunifdis,boxsize_sample,rho_sample,rho_sample_cgs,rho_fracdiff
 real    :: totmass,totmass_req,temp,pmass_cgs,cs0_cgs,u0,u0_cgs,tmax_cgs,dtmax_cgs
 real    :: xmini,xmaxi,ymini,ymaxi,zmini,zmaxi
 integer :: i,ierr,io_file,ix,npmax,npart_sample
 logical :: iexist

 !
 !--general parameters
 !
 time = 0.
 gamma = 1.00011  ! polytropic eos index
 gmw = 1.         ! Mean molecular weight for pure hydrogen

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

    if (.not.use_glass) then
       !
       ! Set number of particles (to be updated after set_unifdis)
       !
       npmax = int(size(xyzh(1,:)))
       np_req = 1E6
       call prompt('Enter total number of particles',np_req,1)
       if (np_req > npmax) call fatal('setup_unifdis_cmi','number of particles exceeded limit')
    else
       !
       ! Get glass file
       !
       glass_filename = 'glassCube_128.dat'
       call prompt('Enter filename of glass cube',glass_filename)
       glass_filename = trim(adjustl(glass_filename))
    endif
    !
    ! Particle mass
    !
    pmass_cgs = 1.9891E30 ! 10^-3 solarm
    pmass = pmass_cgs/umass
    call prompt('Enter particle mass in units of '//mass_unit,pmass,0.)
    !
    ! Density
    !
    rhozero_cgs = 5.21E-21
    call prompt('Enter initial density in g/cm^3',rhozero_cgs,0.)
    print*,'read density:',rhozero_cgs 
    !
    ! set initial sound speed
    !
    cs0_cgs = 0.91E5
    cs0 = cs0_cgs/unit_velocity
    call prompt('Enter initial sound speed in code units',cs0,0.)
    !
    ! Set timestep and end-time
    !
    dtmax_cgs = 3.15360E9   ! 1E-4 Myr
    tmax_cgs  = 4.41504E12  ! 0.14 Myr
    dtmax = dtmax_cgs/utime
    tmax  = tmax_cgs/utime
    call prompt('Enter timestep in code units',dtmax,0.)
    call prompt('Enter simulation end-time in code units',tmax,0.)
    !
    ! Set tree accuracy parameters 
    ! 
    tree_accuracy_cmi = 0.2
    rcut_opennode_cgs = 0.05*pc
    rcut_leafpart_cgs = 0.02*pc
    rcut_opennode = rcut_opennode_cgs/udist 
    rcut_leafpart = rcut_leafpart_cgs/udist 
    call prompt('Enter the tree accuracy for tree-walk',tree_accuracy_cmi,0.)
    call prompt('Enter the radius for leaves in code units',rcut_opennode,0.)
    call prompt('Enter the radius for individual particles in code units',rcut_leafpart,0.)

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
 ! setup particles
 !
 if (.not.use_glass) then
    !
    ! Calculate total mass required with npart_requested
    !
    totmass_req = np_req*pmass
    !
    ! Box density 
    !
    rhozero = rhozero_cgs/unit_density
    !
    ! Calculate volume needed
    !
    vol_box   = totmass_req/rhozero
    boxlength = vol_box**(1./3.)
    xmini = -0.5*boxlength; xmaxi = 0.5*boxlength
    ymini = -0.5*boxlength; ymaxi = 0.5*boxlength
    zmini = -0.5*boxlength; zmaxi = 0.5*boxlength
    !
    ! Set particle distribution
    !
    npart = 0
    deltax = (xmaxi-xmini)/np_req**(1./3.)
    call set_unifdis('cubic',id,master,xmini,xmaxi,ymini,ymaxi,zmini,zmaxi,deltax,hfact,&
                     npart,xyzh_raw,periodic)
    !
    ! Update totmass with the new npart
    !
    totmass = npart*pmass
    !
    ! Update box volume
    !
    vol_box   = totmass/rhozero
    boxlength = vol_box**(1./3.)
    !
    ! Find box size of xyzh_raw from set_unifdis
    !
    boxsize_setunifdis = epsilon(boxsize_setunifdis)
    do i = 1,npart
       do ix = 1,3
          boxsize_setunifdis = max(boxsize_setunifdis,abs(xyzh_raw(ix,i)))
       enddo
    enddo
    boxsize_setunifdis = 2.*boxsize_setunifdis
    boxsize_toscale = boxlength/boxsize_setunifdis
    !
    ! Scale xyzh to the required boxlength
    !
    do i = 1,npart
       do ix = 1,4
          xyzh(ix,i) = xyzh_raw(ix,i) * boxsize_toscale
       enddo
    enddo
    !
    ! Set boundaries
    !
    xmini = -0.5*boxlength; xmaxi = 0.5*boxlength
    ymini = -0.5*boxlength; ymaxi = 0.5*boxlength
    zmini = -0.5*boxlength; zmaxi = 0.5*boxlength
    call set_boundary(xmini,xmaxi,ymini,ymaxi,zmini,zmaxi)
    !
    ! Set particle properties
    !
    npartoftype(:) = 0
    npartoftype(igas) = npart
    massoftype(igas)  = pmass
    do i = 1,npartoftype(igas)
       call set_particle_type(i,igas)
    enddo

 elseif (use_glass) then
    !
    ! Find total number of particles (entries in glass file)
    !
    open(3000,file=glass_filename,status='old',iostat=io_file)
    if (io_file < 0) call fatal('setup_unifdis_cmi','error opening glassfile')
    npart = 0
    do
       read(3000,*,iostat=io_file)
       if (io_file /= 0) then
          exit
       else
          npart = npart + 1
          if (npart > maxrow) call fatal('setup_unifdis_cmi','increase maxrow')
       endif
    enddo
    if (npart == 0) call fatal('setup_unifdis_cmi','glassfile import failed')
    !
    ! Set particle type properties
    !
    npartoftype(:) = 0
    npartoftype(igas) = npart
    massoftype(:) = 0.
    massoftype(igas) = pmass
    do i = 1,npartoftype(igas)
       call set_particle_type(i,igas)
    enddo
    !
    ! Calculate total mass
    !
    totmass = npart * massoftype(igas)
    !
    ! Box density 
    ! 
    rhozero = rhozero_cgs/unit_density
    !
    ! Calculate box size needed
    !
    vol_box   = totmass/rhozero
    boxlength = vol_box**(1./3.)
    !
    ! Read in glass file and find its box size
    !
    rewind(3000)
    boxsize_fromfile = epsilon(boxsize_fromfile)
    do i = 1,npart
       read(3000,*) xyzh_raw(1,i),xyzh_raw(2,i),xyzh_raw(3,i),xyzh_raw(4,i)
       boxsize_fromfile = max(boxsize_fromfile,abs(xyzh_raw(1,i)),abs(xyzh_raw(2,i)),abs(xyzh_raw(3,i)))
    enddo
    close(3000)
    !- catch NaN or 0.
    if (xyzh_raw(1,10) /= xyzh_raw(1,10) .or. xyzh_raw(1,100) == 0.) then
       call fatal('setup_unifdis_cmi','invalid xyzh imported from glass file')
    endif

    boxsize_fromfile = 2.*(boxsize_fromfile + epsilon(boxsize_fromfile))
    boxsize_toscale = boxlength/boxsize_fromfile
    !
    ! Scale positions to box size needed and assign them to particles
    !
    do i = 1,npart
       do ix = 1,4
          xyzh(ix,i) = xyzh_raw(ix,i) * boxsize_toscale
       enddo
    enddo
    !
    ! Set boundaries
    !
    xmini = -0.5*boxlength; xmaxi = 0.5*boxlength
    ymini = -0.5*boxlength; ymaxi = 0.5*boxlength
    zmini = -0.5*boxlength; zmaxi = 0.5*boxlength
    call set_boundary(xmini,xmaxi,ymini,ymaxi,zmini,zmaxi)
 endif

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
 ! Set velocities to 0 and set particle energies using temp
 !
 u0_cgs = kboltz * temp / (gmw*mass_proton_cgs*(gamma-1.))
 u0     = u0_cgs/unit_ergg
 do i = 1,npart
   vxyzu(1:3,i) = 0.
   if (size(vxyzu(:,i)) >= 4) then
     vxyzu(4,i) = u0
   end if
 enddo
 !
 ! Polytropic constant
 !
 if (maxvxyzu < 4 .or. gamma <= 1.) then
    polyk = cs0**2
 else
    polyk = 0.
 endif

 !
 ! Set runtime parameters
 !
 ieos      = 2     ! adiabatic eos
 nout      = 1
 nfulldump = 1
 Tfloor    = 3.
 ufloor    = kboltz*Tfloor/(gmw*mass_proton_cgs*(gamma-1.))/unit_ergg
 icooling  = 0
 alphau    = 1.0
 ipdv_heating   = 1
 ishock_heating = 1

 !
 ! Photoionization settings
 !
 monochrom_source  = .false.
 fix_temp_hii      = .false. 
 implicit_cmi      = .true. 
 treat_Rtype_phase = .false.

 photoionize_tree  = .true.
 nHlimit_fac       = 500
 rcut_opennode_cgs = rcut_opennode*udist 
 rcut_leafpart_cgs = rcut_leafpart*udist 
 delta_rcut_cgs    = 0.05*pc

 limit_voronoi = .true. 
 hlimit_fac = 1E-4
 extradist_fac = 1.0

 !- Print summary - for checking
 print*,' -SUMMARY- '
 print*,'totmass:    ',totmass,totmass*umass,'g'
 print*,'rhozero:    ',rhozero,rhozero*unit_density,'g cm^-3'
 print*,'box vol:    ',vol_box,vol_box*udist**3,'cm^3'
 print*,'box length: ',boxlength,boxlength*udist,'cm'
 print*,'temp:       ',temp,'K'
 print*,'internal E: ',u0,u0_cgs,'erg g^-1'
 print*,'ufloor:     ',ufloor,ufloor*unit_ergg,'erg g^-1'
 print*,'cs0:        ',cs0,cs0*unit_velocity,'cm s^-1'
 print*,'npart:      ',npart
 print*,'part mass:  ',massoftype(igas),massoftype(igas)*umass,'g'
 print*,'dtmax:      ',dtmax,dtmax*utime,'s'
 print*,'tmax:       ',tmax,tmax*utime,'s'

end subroutine setpart

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
 if (.not.use_glass) then
    call write_inopt(np_req,'np_req','total number of particles requested',iunit)
 else
    call write_inopt(glass_filename,'glass_filename','filename of glass cube',iunit)
 endif
 call write_inopt(pmass,'pmass','particle mass',iunit)
 call write_inopt(rhozero_cgs,'rhozero_cgs','initial gas density in g/cm^3',iunit)
 call write_inopt(cs0,'cs0','initial sound speed in code units',iunit)
 call write_inopt(dtmax,'dtmax','timestep in code units',iunit)
 call write_inopt(tmax,'tmax','end-time in code units',iunit)

 call write_inopt(tree_accuracy_cmi,'tree_accuracy_cmi','tree accuracy for cmi nodes',iunit)
 call write_inopt(rcut_opennode,'rcut_opennode','rcut_leaf in code units',iunit)
 call write_inopt(rcut_leafpart,'rcut_leafpart','rcut_part in code units',iunit)

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
 if (.not.use_glass) then
    call read_inopt(np_req,'np_req',db,min=8,errcount=nerr)
 else
    call read_inopt(glass_filename,'glass_filename',db,errcount=nerr)
 endif
 call read_inopt(pmass,'pmass',db,min=0.,errcount=nerr)
 call read_inopt(rhozero_cgs,'rhozero_cgs',db,min=0.,errcount=nerr)
 call read_inopt(cs0,'cs0',db,min=0.,errcount=nerr)
 call read_inopt(dtmax,'dtmax',db,min=0.,errcount=nerr)
 call read_inopt(tmax,'tmax',db,min=0.,errcount=nerr)

 call read_inopt(tree_accuracy_cmi,'tree_accuracy_cmi',db,min=0.,errcount=nerr)
 call read_inopt(rcut_opennode,'rcut_opennode',db,min=0.,errcount=nerr)
 call read_inopt(rcut_leafpart,'rcut_leafpart',db,min=0.,errcount=nerr)

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
