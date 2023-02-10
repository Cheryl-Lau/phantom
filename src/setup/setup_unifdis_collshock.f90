!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setup
!
! Setup routine for uniform distribution injecting two colliding shocks
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - cs0         : *initial sound speed in code units*
!   - dist_unit   : *distance unit (e.g. au)*
!   - ilattice    : *lattice type (1=cubic, 2=closepacked)*
!   - mass_unit   : *mass unit (e.g. solarm)*
!   - npart       : *number of particles*
!   - totmass     : *total mass of particles*
!   - xmax        : *xmax boundary*
!   - xmin        : *xmin boundary*
!   - ymax        : *ymax boundary*
!   - ymin        : *ymin boundary*
!   - zmax        : *zmax boundary*
!   - zmin        : *zmin boundary*
!
! :Dependencies: boundary, cooling, dim, domain, eos, h2cooling,
!   infile_utils, io, mpiutils, options, part, physcon, prompting,
!   set_dust, setup_params, timestep, unifdis, units
!

 implicit none
 public :: setpart

 integer           :: np,npart,ilattice
 real              :: totmass,cs0,xmini,xmaxi,ymini,ymaxi,zmini,zmaxi,rms_mach
 character(len=20) :: dist_unit,mass_unit
 real(kind=8)      :: udist,umass

 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:maxvxyzu
 use setup_params, only:npart_total,rhozero
 use io,           only:master,fatal
 use unifdis,      only:set_unifdis
 use boundary,     only:xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound,set_boundary
 use part,         only:periodic,abundance,iHI
 use part,         only:set_particle_type,igas
 use physcon,      only:pi,mass_proton_cgs,kboltz,years,pc,solarm,micron
 use units,        only:set_units,unit_density,unit_velocity,select_unit
 use domain,       only:i_belong
 use eos,          only:gmw,ieos
 use options,      only:icooling,alpha,alphau,nfulldump
 use timestep,     only:dtmax,tmax,C_cour,C_force,C_cool,tolv,nout
 use cooling,      only:Tfloor
 use prompting,    only:prompt
 use velfield,     only:set_velfield_from_cubes
 use datafiles,    only:find_phantom_datafile
 use inject,       only:init_inject,inject_particles,vel_sk,collaxis,time_sk
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
 character(len=20), parameter     :: filevx = 'cube_v1.dat'
 character(len=20), parameter     :: filevy = 'cube_v2.dat'
 character(len=20), parameter     :: filevz = 'cube_v3.dat'
 character(len=120) :: filex,filey,filez
 character(len=100) :: filename,infilename,cwd
 real    :: vol_box,boxlength,cs0_cgs,deltax,rmsmach,v2i,turbfac,turbboxsize,amplitude,temp
 real    :: vel_sk_cgs,vxyz_avg,vxyz_avg_cgs,vxyz_min,vxyz_min_cgs,vxyz_max,vxyz_max_cgs
 integer :: i,ierr
 logical :: iexist,in_iexist

 !
 !--general parameters
 !
 time = 0.
 if (maxvxyzu < 4) then
    gamma = 1.
 else
    gamma = 5./3.
 endif
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
    ! set boundaries
    !
    boxlength = 0.5 ! pc
    xmini = -0.5*boxlength; xmaxi = 0.5*boxlength
    ymini = -0.5*boxlength; ymaxi = 0.5*boxlength
    zmini = -0.5*boxlength; zmaxi = 0.5*boxlength
    call prompt('enter xmin boundary',xmini)
    call prompt('enter xmax boundary',xmaxi,xmini)
    call prompt('enter ymin boundary',ymini)
    call prompt('enter ymax boundary',ymaxi,ymini)
    call prompt('enter zmin boundary',zmini)
    call prompt('enter zmax boundary',zmaxi,zmini)
    !
    ! number of particles
    !
    np = 1E4
    call prompt(' enter total number of particles',np,1)
    ilattice = 2
    !
    ! total mass
    !
    totmass = 1E0
    call prompt(' enter the total mass of particles',totmass,0.)
    !
    ! sound speed in code units
    !
    cs0_cgs = 7.6E4  ! for around 100K
    cs0 = cs0_cgs/unit_velocity
    call prompt(' enter sound speed in code units',cs0,0.)
    !
    ! set turbulence mach number
    !
    rms_mach = 2.2
    call prompt('Enter the Mach number of the cloud turbulence',rms_mach,0.)

    if (id==master) call write_setupfile(filename)
    stop 'rerun phantomsetup after editing .setup file'
 else
    stop
 endif
 !
 ! set units and boundaries
 !
 call set_units(dist=udist,mass=umass,G=1.d0)
 call set_boundary(xmini,xmaxi,ymini,ymaxi,zmini,zmaxi)
 !
 ! setup particles
 !
 deltax = dxbound/(np**(1./3.))
 npart = 0
 npart_total = 0

 call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,&
                  npart,xyzh,periodic,nptot=npart_total,mask=i_belong)

 npartoftype(:) = 0
 npartoftype(igas) = npart
 print*,' npart = ',npart,npart_total

 if (massoftype(igas) < epsilon(massoftype(igas))) massoftype(igas) = totmass/npart_total

 do i = 1,npartoftype(igas)
    call set_particle_type(i,igas)
 enddo

 print*,' particle mass = ',massoftype(igas)

 temp = (cs0*unit_velocity)**2*gmw*mass_proton_cgs/(gamma*kboltz)
 print*,'gas temp (K) = ',temp

 vol_box = (xmaxi-xmini)*(ymaxi-ymini)*(zmaxi-zmini)
 rhozero = totmass / vol_box

 print*,'density = ',rhozero
 print*,'density_cgs = ',rhozero*unit_density
 print*,'nrho_cgs = ',rhozero*unit_density/mass_proton_cgs

 if (maxvxyzu < 4 .or. gamma <= 1.) then
    polyk = cs0**2
 else
    polyk = 0.
 endif
 do i=1,npart
    vxyzu(1:3,i) = 0.
    if (maxvxyzu >= 4 .and. gamma > 1.) vxyzu(4,i) = cs0**2/(gamma*(gamma-1.))
 enddo
 !
 ! Turbulent velocity field
 !
 vxyzu = 0.
 if (rms_mach > 0.) then
    call getcwd(cwd)

    ! Kennedy or Dial
    filex  = find_phantom_datafile(filevx,'velfield_sphng_small')
    filey  = find_phantom_datafile(filevy,'velfield_sphng_small')
    filez  = find_phantom_datafile(filevz,'velfield_sphng_small')

    turbboxsize = max(abs(xmini),abs(xmaxi),abs(ymini),abs(ymaxi),abs(zmini),abs(zmaxi))
    amplitude   = 1.
    call set_velfield_from_cubes(xyzh(:,1:npart),vxyzu(:,:npart),npart, &
                                 filex,filey,filez,amplitude,turbboxsize,.false.,ierr)
    if (ierr /= 0) call fatal('setup','error setting up velocity field on clouds')

    rmsmach = 0.0
    print*, 'Turbulence being set by user'
    do i = 1,npart
       v2i     = dot_product(vxyzu(1:3,i),vxyzu(1:3,i))
       rmsmach = rmsmach + v2i/cs0**2
    enddo
    rmsmach = sqrt(rmsmach/npart)
    if (rmsmach > 0.) then
       turbfac = rms_mach/rmsmach ! normalise the energy to the desired mach number
     else
       turbfac = 0.
    endif
    do i = 1,npart
       vxyzu(1:3,i) = turbfac*vxyzu(1:3,i)
    enddo

    vxyz_avg = 0.
    vxyz_min = huge(vxyz_min)
    vxyz_max = epsilon(vxyz_max)
    open(2022,file='velfield.txt')
    do i = 1,npart
       vxyz_avg = vxyz_avg + mag(vxyzu(1:3,i))
       if (mag(vxyzu(1:3,i)) < vxyz_min) vxyz_min = mag(vxyzu(1:3,i))
       if (mag(vxyzu(1:3,i)) > vxyz_max) vxyz_max = mag(vxyzu(1:3,i))
       write(2022,*) mag(vxyzu(1:3,i))*unit_velocity*1E-5
    enddo
    close(2022)
    vxyz_avg = vxyz_avg/npart
    vxyz_avg_cgs = vxyz_avg*unit_velocity
    vxyz_min_cgs = vxyz_min*unit_velocity
    vxyz_max_cgs = vxyz_max*unit_velocity
    print*,'vturb (km/s) mean',vxyz_avg_cgs*1E-5,'min',vxyz_min_cgs*1E-5,'max',vxyz_max_cgs*1E-5

 endif
 !
 ! Set default runtime parameters if .in file does not exist
 !
 infilename=trim(fileprefix)//'.in'
 inquire(file=infilename,exist=in_iexist)
 if (.not. in_iexist) then
    tmax      = 1E-2
    dtmax     = 1E-8
    nout      = 1E0
    nfulldump = 1E0
    ieos      = 2
    icooling  = 7   ! JML06 implicit
    Tfloor    = 3.
    ! Set velocity of shock
    vel_sk_cgs = 1E7
    call prompt('enter velocity of shock in cm/s',vel_sk_cgs)
    vel_sk = vel_sk_cgs/unit_velocity
    ! Set axis of collision
    call prompt('enter axis of collision (1-x; 2-y; 3-z)',collaxis)
    ! Set time of shock
    time_sk = 1E-7
    call prompt('enter time of injecting shock',time_sk)
 endif

end subroutine setpart


real function mag(vec)
 real, intent(in) :: vec(3)

 mag = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)

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
 ! boundaries
 !
 write(iunit,"(/,a)") '# boundaries'
 call write_inopt(xmini,'xmin','xmin boundary',iunit)
 call write_inopt(xmaxi,'xmax','xmax boundary',iunit)
 call write_inopt(ymini,'ymin','ymin boundary',iunit)
 call write_inopt(ymaxi,'ymax','ymax boundary',iunit)
 call write_inopt(zmini,'zmin','zmin boundary',iunit)
 call write_inopt(zmaxi,'zmax','zmax boundary',iunit)
 !
 ! other parameters
 !
 write(iunit,"(/,a)") '# setup'
 call write_inopt(np,'np','total number of particles',iunit)
 call write_inopt(totmass,'totmass','total mass of particles',iunit)
 call write_inopt(cs0,'cs0','initial sound speed in code units',iunit)
 call write_inopt(ilattice,'ilattice','lattice type (1=cubic, 2=closepacked)',iunit)
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
 ! boundaries
 !
 call read_inopt(xmini,'xmin',db,errcount=nerr)
 call read_inopt(xmaxi,'xmax',db,min=xmini,errcount=nerr)
 call read_inopt(ymini,'ymin',db,errcount=nerr)
 call read_inopt(ymaxi,'ymax',db,min=ymini,errcount=nerr)
 call read_inopt(zmini,'zmin',db,errcount=nerr)
 call read_inopt(zmaxi,'zmax',db,min=zmini,errcount=nerr)
 !
 ! other parameters
 !
 call read_inopt(np,'np',db,min=8,errcount=nerr)
 call read_inopt(totmass,'totmass',db,min=0.,errcount=nerr)
 call read_inopt(cs0,'cs0',db,min=0.,errcount=nerr)
 call read_inopt(ilattice,'ilattice',db,min=1,max=2,errcount=nerr)
 call read_inopt(rms_mach,'rms_mach',db,errcount=nerr)
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
