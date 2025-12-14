!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setup
!
! This module sets up an ellipsoid with an envalope
!
! :References: None
!
! :Owner: Cheryl Lau
!
! :Runtime parameters:
!   - dist_unit        : *distance unit (e.g. au)*
!   - mass_unit        : *mass unit (e.g. solarm)*
!   - mpart_solarm     : *mass of particle in solarm units*
!   - np_cloud         : *requested number of particles in cloud*
!   - np_envelope      : *requested number of particles within envelope boundaries*
!   - rho_cloud_cgs    : *density of cloud*
!   - rho_envelope_cgs : *density of envelope*
!   - r1_ellipsoid     : *ratio of semi-axis a of both cloud and envelope*
!   - r2_ellipsoid     : *ratio of semi-axis b of both cloud and envelope*
!   - r3_ellipsoid     : *ratio of semi-axis c of both cloud and envelope*
!   - rms_mach         : *Mach number for turbulence*
!   - gmw_in           : *mean molecular weight*
!   - make_sinks       : *flag to dynamically create sinks*
!
! :Dependencies: dim, io, physcon, setup_params, spherical, boundary, prompting, units,
!                eos, part, ptmass, timestep, kernel, options, datafiles
!                photoionize_cmi, inject
!
 implicit none
 public :: setpart

 private

 integer      :: np_cloud,np_envelope
 real(kind=8) :: udist,umass
 real         :: mpart_solarm,rho_cloud_cgs,rho_envelope_cgs,rms_mach,gmw_in
 real         :: r1_ellipsoid,r2_ellipsoid,r3_ellipsoid
 logical      :: make_sinks
 character(len=20) :: dist_unit,mass_unit

contains

!----------------------------------------------------------------
!+
!  setup for an ellipsoid with an outer envelope
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use physcon,      only:pi,solarm,years,au,pc,kboltz,mass_proton_cgs,gg,Rg
 use dim,          only:maxvxyzu
 use setup_params, only:rhozero,npart_total
 use io,           only:master,fatal,iprint,iverbose
 use spherical,    only:set_ellipse
 use boundary,     only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use prompting,    only:prompt
 use units,        only:set_units,select_unit,utime,unit_density,unit_velocity,unit_energ,unit_ergg
 use eos,          only:polyk2,ieos,gmw
 use part,         only:igas,abundance,set_particle_type
 use part,         only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft
 use ptmass,       only:icreate_sinks,rho_crit_cgs,r_crit,h_acc,h_soft_sinksink,h_soft_sinkgas
 use ptmass,       only:r_merge_cond,r_merge_uncond
 use timestep,     only:dtmax,tmax,dtwallmax,C_cour,C_force,C_cool,tolv,nout
 use options,      only:nfulldump,nmaxdumps,icooling,alpha,alphau
 use kernel,       only:hfact_default
 use domain,       only:i_belong
 use cooling,      only:Tfloor,ufloor
 use velfield,     only:set_velfield_from_cubes
 use datafiles,    only:find_phantom_datafile
 use options,      only:ipdv_heating,ishock_heating
 use photoionize_cmi, only:monochrom_source,fix_temp_hii,treat_Rtype_phase
 use photoionize_cmi, only:photoionize_tree,tree_accuracy_cmi,nHlimit_fac
 use photoionize_cmi, only:rcut_opennode_cgs,rcut_leafpart_cgs,delta_rcut_cgs
 use photoionize_cmi, only:sink_ionsrc,inject_rad
 use photoionize_cmi, only:nsetphotosrc,xyztq_setphotosrc_cgs
 use inject,          only:inject_sn,sink_progenitor
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 character(len=20), parameter     :: filevx = 'cube_v1.dat'
 character(len=20), parameter     :: filevy = 'cube_v2.dat'
 character(len=20), parameter     :: filevz = 'cube_v3.dat'
 real,  allocatable :: xyzh_outer(:,:)
 real(kind=8)       :: h_acc_in
 integer            :: ip,npmax,ierr,i
 integer            :: npart_cloud,npart_outer,npmax_env,temp_xyzh_size
 integer            :: npart_temp
 integer(kind=8)    :: npart_total_outer
 real               :: r_cloud(3),r_outer(3)
 real               :: psep_cloud,psep_outer,x,y,z,h,rmax,scale_param
 real               :: cs_cloud,cs_cloud_cgs,rho_cloud,rho_outer,temp_cloud,vol_cloud,vol_outer
 real               :: t_ff,jeans_mass,jeans_mass_cgs,totmass_cloud,totmass_outer
 real               :: temp_envelope,cs_envelope_cgs,cs_envelope
 real               :: rmsmach,v2i,turbfac,turbboxsize,v2_sum,v_rms,v_rms_kms
 real               :: r_sn_cgs,engsn_cgs,pmsncrit_cgs
 real               :: h_acc_cgs,h_soft_sinksink_cgs,h_soft_sinkgas_cgs,rho_crit_cgs_recomm
 logical            :: iexist,in_iexist
 logical            :: remove_envelope     = .false.  ! temporarily remove envelope to see virial
 logical            :: place_sink_in_setup = .false.
 character(len=120) :: filex,filey,filez
 character(len=100) :: filename,infilename
 character(len=40)  :: fmt,lattice
 character(len=9)   :: proceed

 print "(a)", 'Setup for an elongated cloud (ellipsoid) with an envelope'

 filename = trim(fileprefix)//'.setup'
 inquire(file = filename,exist=iexist)
 if (iexist) then
    call read_setupfile(filename,ierr)
    if (ierr /= 0) then
       if (id == master) call write_setupfile(filename)
       stop
    endif
 elseif (id == master) then
    print "(a,/)",trim(filename)//' not found: using interactive setup'
    dist_unit = 'pc'
    mass_unit = 'solarm'
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

    !- Units
    call set_units(dist=udist,mass=umass,G=1.d0)

    !- Limit the total number of particles
    npmax = int(0.9*size(xyzh(1,:))) ! approx max number allowed in sphere given size(xyzh(1,:))

    !- Mass of particles
    mpart_solarm = 1E-2
    call prompt('Enter the mass of particles in solarm units',mpart_solarm,0.)

    !- Settings for the cloud
    np_cloud = 1E6
    call prompt('Enter the approximate number of particles in the cloud',np_cloud,0,npmax)
    rho_cloud_cgs = 1E-21
    call prompt('Enter the density of the cloud in g/cm^3',rho_cloud_cgs,0.)

    !- Settings for the envelope
    np_envelope = 4E3
    npmax_env = npmax - np_cloud
    call prompt('Enter the approximate number of particles within the envelope boundaries',np_envelope,0,npmax_env)
    rho_envelope_cgs = 5E-25
    call prompt('Enter the density of the envelope in g/cm^3',rho_envelope_cgs,0.)

    !- Ratio of semi-axes of ellipsoid
    r1_ellipsoid = 5
    r2_ellipsoid = 2.5
    r3_ellipsoid = 2.5
    call prompt('Enter the ratio of semi-axis a of the ellipsoid',r1_ellipsoid,0.)
    call prompt('Enter the ratio of semi-axis b of the ellipsoid',r2_ellipsoid,0.)
    call prompt('Enter the ratio of semi-axis c of the ellipsoid',r3_ellipsoid,0.)

    !- Turbulence (Virial ratio)
    rms_mach = 7.4
    call prompt('Enter the Mach number of the cloud turbulence',rms_mach,0.)

    !- Mean molecular weight
    gmw_in = 1.29
    call prompt('Enter the mean molecular weight',gmw_in,0.)

    make_sinks = .true.
    call prompt('Do you wish to dynamically create sink particles? ',make_sinks)

    if (id == master) call write_setupfile(filename)
    stop 'please edit .setup file and rerun phantomsetup'
 else
    stop ! MPI, stop on other threads, interactive on master
 endif

 !
 ! units
 !
 call set_units(dist=udist,mass=umass,G=1.d0)

 !
 ! set gamma
 !
 if (maxvxyzu >= 4) then
    gamma = 5./3.
 else
    gamma = 1.
 endif

 !
 ! general parameters
 !
 time        = 0.
 hfact       = hfact_default

 !
 ! Set up the cloud
 !
 massoftype(igas) = mpart_solarm*solarm/umass
 npart_cloud   = np_cloud
 totmass_cloud = npart_cloud*massoftype(igas)
 rho_cloud = rho_cloud_cgs/unit_density
 vol_cloud = totmass_cloud/rho_cloud
 scale_param = (vol_cloud/(4./3.*pi*r1_ellipsoid*r2_ellipsoid*r3_ellipsoid))**(1./3.)

 r_cloud(1) = r1_ellipsoid*scale_param
 r_cloud(2) = r2_ellipsoid*scale_param
 r_cloud(3) = r3_ellipsoid*scale_param
 t_ff       = sqrt(3.*pi/(32.*rho_cloud))
 lattice    = 'closepacked'

 call set_ellipse(trim(lattice),id,master,r_cloud,vol_cloud,psep_cloud,hfact,xyzh,npart,&
                  nptot=npart_total,exactN=.true.,np_requested=npart_cloud,mask=i_belong)

 !- Update parameters with new npart_cloud
 npart_cloud = npart_total
 totmass_cloud = npart_cloud*massoftype(igas)
 rho_cloud = totmass_cloud/vol_cloud

 !
 ! Temperature and Sound speed of cloud
 !
 temp_cloud   = get_eqtemp_from_rho(rho_cloud_cgs)
 cs_cloud_cgs = sqrt(temp_cloud*(gamma*kboltz)/(gmw*mass_proton_cgs))
 cs_cloud     = cs_cloud_cgs/unit_velocity

 !
 ! Apply turbulence to the cloud
 !
 vxyzu = 0.
 if (rms_mach > 0.) then
    filex  = find_phantom_datafile(filevx,'velfield_sphng_small')
    filey  = find_phantom_datafile(filevy,'velfield_sphng_small')
    filez  = find_phantom_datafile(filevz,'velfield_sphng_small')

    rmax = max(r_cloud(1),r_cloud(2),r_cloud(3))
    turbboxsize = 1.1*rmax
    call set_velfield_from_cubes(xyzh(:,1:npart_cloud),vxyzu(:,:npart_cloud),npart_cloud, &
                                 filex,filey,filez,1.,turbboxsize,.false.,ierr)
    if (ierr /= 0) call fatal('setup','error setting up velocity field on clouds')

    rmsmach = 0.0
    do ip = 1,npart_cloud
       v2i     = dot_product(vxyzu(1:3,ip),vxyzu(1:3,ip))
       rmsmach = rmsmach + v2i/cs_cloud**2
    enddo
    rmsmach = sqrt(rmsmach/npart_cloud)
    if (rmsmach > 0.) then
       turbfac = rms_mach/rmsmach ! normalise the energy to the desired mach number
    else
       turbfac = 0.
    endif
    do ip = 1,npart_cloud
       vxyzu(1:3,ip) = turbfac*vxyzu(1:3,ip)
    enddo

    !- Calculate rms-velocity
    v2_sum = 0.
    do ip = 1,npart_cloud
       v2_sum = v2_sum + mag2(vxyzu(1:3,ip))
    enddo
    v_rms = sqrt(1./npart_cloud*v2_sum)
    v_rms_kms = v_rms*unit_velocity*1E-5
 endif

 !
 ! Set initial temperature
 !
 do ip = 1,npart_cloud
    vxyzu(4,ip) = kboltz * temp_cloud / (gmw*mass_proton_cgs*(gamma-1.)) /unit_ergg
 enddo


 !
 ! Set up the outer ellipsoid then remove the cloud to create the envelope
 !
 if (.not.remove_envelope) then
    npart_outer = np_envelope
    totmass_outer = npart_outer*massoftype(igas)
    rho_outer = rho_envelope_cgs/unit_density
    vol_outer = totmass_outer/rho_outer
    scale_param = (vol_outer/(4./3.*pi*r1_ellipsoid*r2_ellipsoid*r3_ellipsoid))**(1./3.)
    r_outer(1) = r1_ellipsoid*scale_param
    r_outer(2) = r2_ellipsoid*scale_param
    r_outer(3) = r3_ellipsoid*scale_param

    temp_xyzh_size = npart_outer*1.2
    npart_temp = 0
    npart_total_outer = 0
    allocate(xyzh_outer(4,temp_xyzh_size))

    call set_ellipse(trim(lattice),id,master,r_outer,vol_outer,psep_outer,hfact,xyzh_outer,npart_temp,&
                     nptot=npart_total_outer,exactN=.true.,np_requested=npart_outer,mask=i_belong)

    !- Update parameters
    totmass_outer = npart_total_outer*massoftype(igas)
    rho_outer = totmass_outer/vol_outer
    rho_envelope_cgs = rho_outer*unit_density

    !
    ! Temperature and Sound speed of envelope
    !
    temp_envelope   = get_eqtemp_from_rho(rho_envelope_cgs)
    cs_envelope_cgs = sqrt(temp_envelope*(gamma*kboltz)/(gmw*mass_proton_cgs))
    cs_envelope     = cs_envelope_cgs/unit_velocity

    !
    ! Add to real particles if particle is outside cloud i.e. is envelope
    !
    npart = npart_cloud
    do ip = 1,npart_total_outer
       x = xyzh_outer(1,ip)
       y = xyzh_outer(2,ip)
       z = xyzh_outer(3,ip)
       if ((x/r_cloud(1))**2 + (y/r_cloud(2))**2 + (z/r_cloud(3))**2 > 0.99) then
          npart = npart + 1
          h = xyzh_outer(4,ip)
          xyzh(1:4,npart) = (/ x,y,z,h /)
          vxyzu(4,npart) = kboltz * temp_envelope / (gmw*mass_proton_cgs*(gamma-1.)) /unit_ergg
       endif
    enddo
    deallocate(xyzh_outer)

    if (npart /= npart_cloud+npart_total_outer) then
       print*,npart_cloud+npart_total_outer-npart,'particles removed'
       print*,'actual number of particles in the envelope ',npart-npart_cloud
    else
       stop 'no particles removed!'
    endif
 endif

 !
 ! set particle properties
 !
 npartoftype(:)    = 0
 npartoftype(igas) = npart
 do i = 1,npartoftype(igas)
    call set_particle_type(i,igas)
 enddo

 !
 ! Manually place a sink as feedback source
 !
 if (place_sink_in_setup) then
    nptmass                      = 1
    xyzmh_ptmass(:,:)            = 0.
    xyzmh_ptmass(1:3,nptmass)    = 0.
    xyzmh_ptmass(4,nptmass)      = 40.*solarm/umass
    xyzmh_ptmass(ihacc,nptmass)  = 0.005*pc/udist
    xyzmh_ptmass(ihsoft,nptmass) = 0.005*pc/udist
    vxyz_ptmass                  = 0.
 endif

 !
 ! Set default runtime parameters if .in file does not exist
 !
 infilename=trim(fileprefix)//'.in'
 inquire(file=infilename,exist=in_iexist)
 if (.not.in_iexist) then
    tmax      = 1.0*t_ff
    dtmax     = 1d-4*t_ff
    nout      = 10
    nfulldump = 1
    nmaxdumps = 1000
    dtwallmax = 0.d0
    iverbose  = 1

    ieos      = 2    ! adiabatic eos with P = (gamma-1)*rho*u
    gmw       = gmw_in
    icooling  = 7
    Tfloor    = 3.
    ufloor    = kboltz*Tfloor/(gmw*mass_proton_cgs*(gamma-1.))/unit_ergg
    ipdv_heating   = 1
    ishock_heating = 1

    !
    ! Sinks settings
    !
    if (make_sinks) then
       icreate_sinks    = 1
       h_acc            = 5.d0*au/udist
       r_crit           = 5.d0*h_acc
       rho_crit_cgs     = 1.d-15 
       rho_crit         = rho_crit_cgs/unit_density
       h_soft_sinkgas   = h_acc
       h_soft_sinksink  = 0.d0
    else
       icreate_sinks = 0
    endif

    !
    ! Photoionization settings
    !
    inject_rad = .false.
    sink_ionsrc = .true.

    monochrom_source  = .false.
    fix_temp_hii      = .false.
    treat_Rtype_phase = .true.

    photoionize_tree  = .true.
    tree_accuracy_cmi = 0.2
    nHlimit_fac       = 100.
    rcut_opennode_cgs = 3.0*pc
    rcut_leafpart_cgs = 1.0*pc
    delta_rcut_cgs    = 0.5*pc

    !
    ! Supernova settings
    !
    inject_sn = .false.
    sink_progenitor = .false.

 endif

 !
 ! Calculate the approx number of stars that will form
 !
 jeans_mass_cgs = (5.*Rg*temp_cloud/(2.*gg*gmw))**(3./2.) * (4./3.*pi*rho_cloud_cgs)**(-1./2.)
 jeans_mass = jeans_mass_cgs/umass


 print*,'-Cloud-'
 print*,'total mass        ',totmass_cloud,mass_unit
 print*,'Jeans mass        ',jeans_mass,mass_unit
 print*,'free-fall time    ',t_ff*utime/(1E6*years),'Myr'
 print*,'temperature       ',temp_cloud,'K'
 print*,'sound speed       ',cs_cloud_cgs*1E-5,'km/s'
 print*,'v_rms             ',v_rms_kms,'km/s'

 if (.not.remove_envelope) then
    print*,'-Envelope-'
    print*,'temperature       ',temp_envelope,'K'
    print*,'sound speed       ',cs_envelope_cgs*1E-5,'km/s'
 endif

 proceed = 'y'
 call prompt('Do you wish to continue?',proceed)
 if (trim(adjustl(proceed)) == 'n') then
    stop
 elseif (trim(adjustl(proceed)) /= 'y') then
    stop 'Invalid input'
 else
    open(2022,file='cloud_envlp_info.txt')
    write(2022,*) '-Simulation-'
    write(2022,*) 'tmax              ',tmax*utime/(1E6*years),'Myr / ',tmax*utime,'s'
    write(2022,*) 'dtmax             ',dtmax*utime/(1E6*years),'Myr / ',dtmax*utime,'s'
    write(2022,*) '-Cloud-'
    write(2022,*) 'density           ',rho_cloud_cgs,'g/cm^3'
    write(2022,*) 'total mass        ',totmass_cloud,mass_unit
    write(2022,*) 'Jeans mass        ',jeans_mass,mass_unit
    write(2022,*) 'semi-axis a       ',r_cloud(1),dist_unit
    write(2022,*) 'semi-axis b       ',r_cloud(2),dist_unit
    write(2022,*) 'semi-axis c       ',r_cloud(3),dist_unit
    write(2022,*) 'volume            ',vol_cloud*udist**3/pc**3,'pc^3'
    write(2022,*) 'free-fall time    ',t_ff,'/ ',t_ff*utime/(1E6*years),'Myr / ',t_ff*utime,'s'
    write(2022,*) 'temperature       ',temp_cloud,'K'
    write(2022,*) 'sound speed       ',cs_cloud_cgs*1E-5,'km/s'
    write(2022,*) 'Mach number       ',rms_mach
    write(2022,*) 'v_rms             ',v_rms_kms,'km/s'
    if (.not.remove_envelope) then
       write(2022,*) '-Envelope-'
       write(2022,*) 'density           ',rho_envelope_cgs,'g/cm^3'
       write(2022,*) 'temperature       ',temp_envelope,'K'
       write(2022,*) 'sound speed       ',cs_envelope_cgs*1E-5,'km/s'
    endif
    close(2022)
 endif

end subroutine setpart



real function mag2(vec)
 real, intent(in) :: vec(3)

 mag2 = vec(1)**2 + vec(2)**2 + vec(3)**2

end function mag2

!
! Estimate temperature with given rho from cooling curve
!
real function get_eqtemp_from_rho(rho_cgs)
 use datafiles, only:find_phantom_datafile
 real,   intent(in) :: rho_cgs
 integer, parameter :: maxrow = 1000
 integer, parameter :: maxcol = 4
 integer :: irow,icol,iclosest
 real    :: rhoteq(maxcol,maxrow)

 open(3000,file=find_phantom_datafile('rhoteq_table.txt','coolcurve_rho_temp'))
 rewind(3000)
 do irow = 1,maxrow
    read(3000,*) (rhoteq(icol,irow), icol=1,maxcol)
 enddo
 close(3000)

 iclosest = minloc(abs(rhoteq(1,:)-rho_cgs),1)
 get_eqtemp_from_rho = rhoteq(3,iclosest)

end function get_eqtemp_from_rho

!----------------------------------------------------------------
!+
!  write parameters to setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only: write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter           :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for sphere setup routines'
 write(iunit,"(/,a)") '# units'
 call write_inopt(dist_unit,'dist_unit','distance unit (e.g. au)',iunit)
 call write_inopt(mass_unit,'mass_unit','mass unit (e.g. solarm)',iunit)

 write(iunit,"(/,a)") '# general parameters'
 call write_inopt(mpart_solarm,'mpart_solarm','mass of particles in solarm units',iunit)
 call write_inopt(rms_mach,'rms_mach','turbulent rms mach number',iunit)
 call write_inopt(gmw_in,'gmw_in','mean molecular weight',iunit)
 call write_inopt(make_sinks,'make_sinks','dynamically create sink particles',iunit)
 call write_inopt(r1_ellipsoid,'r1_ellipsoid','ratio of semi-axis a of ellipsoid',iunit)
 call write_inopt(r2_ellipsoid,'r2_ellipsoid','ratio of semi-axis b of ellipsoid',iunit)
 call write_inopt(r3_ellipsoid,'r3_ellipsoid','ratio of semi-axis c of ellipsoid',iunit)

 write(iunit,"(/,a)") '# cloud settings'
 call write_inopt(np_cloud,'np_cloud','requested number of particles in cloud',iunit)
 call write_inopt(rho_cloud_cgs,'rho_cloud_cgs','density of cloud in g/cm^3',iunit)

 write(iunit,"(/,a)") '# envelope settings'
 call write_inopt(np_envelope,'np_envelope','requested number of particles within envelope boundaries',iunit)
 call write_inopt(rho_envelope_cgs,'rho_envelope_cgs','density of envelope in g/cm^3',iunit)

 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  Read parameters from setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only: open_db_from_file,inopts,read_inopt,close_db
 use io,           only: error
 use units,        only: select_unit
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit = 21
 type(inopts), allocatable     :: db(:)
 integer :: nerr

 !--Read values
 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(mass_unit,'mass_unit',db,ierr)
 call read_inopt(dist_unit,'dist_unit',db,ierr)
 call read_inopt(mpart_solarm,'mpart_solarm',db,ierr)
 call read_inopt(np_cloud,'np_cloud',db,ierr)
 call read_inopt(rho_cloud_cgs,'rho_cloud_cgs',db,ierr)
 call read_inopt(np_envelope,'np_envelope',db,ierr)
 call read_inopt(rho_envelope_cgs,'rho_envelope_cgs',db,ierr)
 call read_inopt(r1_ellipsoid,'r1_ellipsoid',db,ierr)
 call read_inopt(r2_ellipsoid,'r2_ellipsoid',db,ierr)
 call read_inopt(r3_ellipsoid,'r3_ellipsoid',db,ierr)
 call read_inopt(rms_mach,'rms_mach',db,ierr)
 call read_inopt(gmw_in,'gmw_in',db,ierr)
 call read_inopt(make_sinks,'make_sinks',db,ierr)

 call close_db(db)
 !
 ! parse units
 !
 call select_unit(mass_unit,umass,nerr)
 if (nerr /= 0) then
    call error('setup_elongated_envelope','mass unit not recognised')
    ierr = ierr + 1
 endif
 call select_unit(dist_unit,udist,nerr)
 if (nerr /= 0) then
    call error('setup_elongated_envelope','length unit not recognised')
    ierr = ierr + 1
 endif
 if (ierr > 0) then
    print "(1x,a,i2,a)",'setup_elongated_envelope: ',nerr,' error(s) during read of setup file.  Re-writing.'
 endif

end subroutine read_setupfile

end module setup
