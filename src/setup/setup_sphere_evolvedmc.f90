!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setup
!
! This module sets up a sphere with no surrounding medium
! Used for evolving a cloud before injecting SN
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - BEfac            : *over-density factor of the BE sphere [code units]*
!   - BEmass           : *mass radius of the BE sphere [code units]*
!   - BErad_norm       : *normalised radius of the BE sphere*
!   - BErad_phys       : *physical radius of the BE sphere [code units]*
!   - BErho_cen        : *central density of the BE sphere [code units]*
!   - Bzero            : *Magnetic field strength in Gauss*
!   - ang_Bomega       : *Angle (degrees) between B and rotation axis*
!   - angvel           : *angular velocity in rad/s*
!   - cs_sphere_cgs    : *sound speed in sphere/ellipsoid in cm/s*
!   - dist_unit        : *distance unit (e.g. au)*
!   - h_acc            : *accretion radius (code units)*
!   - h_soft_sinksink  : *sink-sink softening radius (code units)*
!   - iBE_options      : *The set of parameters to define the BE sphere*
!   - icreate_sinks    : *1: create sinks.  0: do not create sinks*
!   - is_sphere        : *set up a sphere; otherwise an ellipsoid*
!   - lattice          : *either 'random' or 'closepacked'*
!   - lbox             : *length of a box side in terms of spherical radii*
!   - lbox_x/y/z       : *dimensions of rectangular box in terms of elliptical semiaxes*
!   - mass_unit        : *mass unit (e.g. solarm)*
!   - masstoflux       : *mass-to-magnetic flux ratio in units of critical value*
!   - mc_method        : *applying Monte Carlo method to randomize particles; equivalent to lattice='random'*
!   - np               : *requested number of particles in sphere*
!   - pos_ranh         : *the intent is to slightly randomize particle positions on closepacked lattice*
!   - r_crit           : *critical radius (code units)*
!   - r_sphere         : *radius of sphere in code units*
!   - r_ellipsoid      : *semi-axes of ellipsoid a, b and c in code units*
!   - rho_pert_amp     : *amplitude of density perturbation*
!   - totmass_sphere   : *mass of sphere/ellipsoid in code units*
!   - use_BE_sphere    : *centrally condense as a BE sphere*
!
! :Dependencies: boundary, centreofmass, dim, domain, eos, infile_utils,
!   io, kernel, options, part, physcon, prompting, ptmass, rho_profile,
!   setup_params, spherical, timestep, unifdis, units
!
 use part,    only:periodic
 use dim,     only:use_dust,maxvxyzu,periodic
 implicit none
 public :: setpart

 private
 !--private module variables
 real :: xmini(3),xmaxi(3)
 real :: totmass_sphere,r_sphere,r_ellipsoid(3),cs_sphere,cs_sphere_cgs
 real :: angvel,Bzero_G,masstoflux,ang_Bomega,rms_mach
 real :: rho_pert_amp,lbox,lbox_x,lbox_y,lbox_z
 real :: BErho_cen,BErad_phys,BErad_norm,BEmass,BEfac
 real :: rho_crit_cgs_setup,r_crit_setup,h_acc_setup,h_soft_sinksink_setup,rhofinal_setup
 real(kind=8)     :: udist,umass
 integer          :: np,iBEparam,icreate_sinks_setup
 logical          :: BEsphere,mu_not_B,cs_in_code,pos_ranh,mc_method
 logical          :: is_sphere
 character(len=20)            :: dist_unit,mass_unit,lattice
 character(len= 1), parameter :: labelx(3) = (/'x','y','z'/)

contains

!----------------------------------------------------------------
!+
!  setup for a sphere
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use physcon,      only:pi,solarm,hours,years,au,kboltz,mass_proton_cgs
 use dim,          only:maxvxyzu,h2chemistry,gr
 use setup_params, only:rhozero,npart_total,ihavesetupB  !,rmax
 use io,           only:master,fatal,iprint,iverbose
 use unifdis,      only:set_unifdis
 use spherical,    only:set_sphere,set_ellipse,set_unifdis_sphereN
 use rho_profile,  only:rho_bonnorebert,prompt_BEparameters
 use boundary,     only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use prompting,    only:prompt
 use units,        only:set_units,select_unit,utime,unit_density,unit_Bfield,unit_velocity,unit_pressure,unit_energ,unit_ergg
 use eos,          only:polyk2,ieos,gmw
 use part,         only:Bxyz,Bextx,Bexty,Bextz,igas,idust,abundance,iHI,set_particle_type
 use timestep,     only:dtmax,tmax,dtwallmax,C_cour,C_force,C_cool,tolv,nout
 use centreofmass, only:reset_centreofmass
 use options,      only:nfulldump,nmaxdumps,icooling,alpha,alphau
 use kernel,       only:hfact_default
 use domain,       only:i_belong
 use ptmass,       only:icreate_sinks,rho_crit,rho_crit_cgs,r_crit,h_acc,h_soft_sinksink
 use cooling,      only:Tfloor,ufloor
 use h2cooling,    only:abundc,abundo,abundsi,abunde,dust_to_gas_ratio,iphoto
 use velfield,     only:set_velfield_from_cubes
 use datafiles,    only:find_phantom_datafile
 use random,       only:ran2
 use options,      only:ipdv_heating,ishock_heating
 use stretchmap,   only:rho_func
 use photoionize_cmi, only:monochrom_source,fix_temp_hii,treat_Rtype_phase
 use photoionize_cmi, only:photoionize_tree,tree_accuracy_cmi,nHlimit_fac
 use photoionize_cmi, only:rcut_opennode_cgs,rcut_leafpart_cgs,delta_rcut_cgs
 use photoionize_cmi, only:old_sources_exist
 procedure(rho_func), pointer :: density_func
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
 real(kind=8)       :: h_acc_in
 integer            :: i,nx,np_in,npartsphere,npmax,iBElast,iseed,ierr
 integer            :: iBE
 real               :: rmin,rmax
 real               :: totmass,vol_box,psep,psep_box
 real               :: vol_sphere,dens_sphere,dens_medium,cs_medium,angvel_code,przero
 real               :: totmass_box,t_ff,r2,area,Bzero,rmasstoflux_crit
 real               :: rxy2,rxyz2,phi,dphi,central_density,edge_density,rmsmach,v2i,turbfac,turbboxsize
 real               :: r_sn_cgs,engsn_cgs,pmsncrit_cgs
 real               :: vxyz_avg,vxyz_min,vxyz_max,vxyz_avg_cgs,vxyz_min_cgs,vxyz_max_cgs
 real, allocatable  :: rtab(:),rhotab(:)
 logical            :: iexist,in_iexist
 logical            :: make_sinks = .true.
 character(len=120) :: filex,filey,filez
 character(len=100) :: filename,infilename,cwd
 character(len=40)  :: fmt
 character(len=10)  :: h_acc_char
 character(len=9)   :: c_shape

 print "(a)", 'Sphere or Ellipsoid setup for evolving a cloud'

 filename = trim(fileprefix)//'.setup'
 inquire(file=filename,exist=iexist)
 if (iexist) then
    call read_setupfile(filename,ierr)
    np_in = np
    if (ierr /= 0) then
       if (id==master) call write_setupfile(filename)
       stop
    endif
 elseif (id==master) then
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
    !
    ! units
    !
    call set_units(dist=udist,mass=umass,G=1.d0)
    !
    ! Set number of particles
    !
    npmax = int(2.0/3.0*size(xyzh(1,:))) ! approx max number allowed in sphere given size(xyzh(1,:))
    print*, 'npmax', npmax
    np = 1E6
    !
    ! prompt user for settings
    !
    c_shape = 'sphere' ! 'ellipsoid'
    call prompt('Is the cloud initially a [sphere] or an [ellipsoid]?',c_shape)

    if (c_shape=='sphere') then
       is_sphere = .true.

       BEsphere = .false.
       call prompt('Centrally condense the sphere as a BE sphere?',BEsphere)
       if (BEsphere) then
          print "(a)",' ERROR: Cannot set BEsphere if the box is not a cube'
          stop
       endif

       if (.not. BEsphere) then
          call prompt('Enter the approximate number of particles in the sphere',np,0,npmax)
          np_in    = np
          r_sphere = 5.
          call prompt('Enter radius of sphere in units of '//dist_unit,r_sphere,0.)
          totmass_sphere = 1E4
          call prompt('Enter total mass in sphere in units of '//mass_unit,totmass_sphere,0.)

       else
          call prompt_BEparameters(iBEparam,BErho_cen,BErad_phys,BErad_norm,BEmass,BEfac,umass,udist,au,solarm)
       endif

    elseif (c_shape=='ellipsoid') then
       is_sphere = .false.

       call prompt('Enter the approximate number of particles in the ellipsoid',np,0,npmax)
       np_in          = np
       r_ellipsoid(1) = 5.
       r_ellipsoid(2) = 4.
       r_ellipsoid(3) = 4.
       call prompt('Enter the semi-axis {a} of ellipsoid in units of '//dist_unit,r_ellipsoid(1),0.)
       call prompt('Enter the semi-axis {b} of ellipsoid in units of '//dist_unit,r_ellipsoid(2),0.)
       call prompt('Enter the semi-axis {c} of ellipsoid in units of '//dist_unit,r_ellipsoid(3),0.)
       totmass_sphere = 1E4
       call prompt('Enter total mass in ellipsoid in units of '//mass_unit,totmass_sphere,0.)

    else
       print "(a)",' ERROR: shape not recognised - Cloud has to be either a sphere or an ellipsoid.'
       stop
    endif

    cs_sphere_cgs = 21888.0  ! cm/s ~ 8K assuming mu = 2.31 & gamma = 5/3
    call prompt('Enter sound speed in sphere in units of cm/s',cs_sphere_cgs,0.)

    angvel = 0.
    call prompt('Enter angular rotation speed in rad/s ',angvel,0.)

    mc_method = .false.
    call prompt('Do you intend to use Monte Carlo method for assigning particle positions?',mc_method)
    if (mc_method) then
       lattice = 'random'
    else
       lattice = 'closepacked'
       pos_ranh = .false.
       call prompt('Do to intend to slightly randomize particle positions on the lattice?',pos_ranh)
    endif

    rms_mach = 14.2
    call prompt('Enter the Mach number of the cloud turbulence',rms_mach,0.)

    !
    ! sink particle details to go to .in file
    !
    call prompt('Do you wish to dynamically create sink particles? ',make_sinks)
    if (make_sinks) then
       h_acc_char  = '1.d-3pc'
       call prompt('Enter the accretion radius of the sink (with units; e.g. au,pc,kpc,0.1pc) ',h_acc_char)
       call select_unit(h_acc_char,h_acc_in,ierr)
       h_acc_setup = h_acc_in
       if (ierr==0 ) h_acc_setup = h_acc_setup/udist
       r_crit_setup        = 5.0*h_acc_setup
       rho_crit_cgs_setup  = 1.d-15   ! 10^5 times the initial MC density (Bate etal 1995); default 1.e-10
       icreate_sinks_setup = 1
    else
       icreate_sinks_setup = 0
       rhofinal_setup = 0
       call prompt('Enter final maximum density in g/cm^3 (ignored for <= 0) ',rhofinal_setup)
    endif
    if (id==master) call write_setupfile(filename)
    stop 'please edit .setup file and rerun phantomsetup'
 else
    stop ! MPI, stop on other threads, interactive on master
 endif
 !
 ! units
 !
 call set_units(dist=udist,mass=umass,G=1.d0)
 !
 ! Sphere/Ellipsoid shared parameters
 !
 if (is_sphere) then
    rmax = r_sphere
    vol_sphere  = 4./3.*pi*r_sphere**3
 else
    rmax = max(r_ellipsoid(1),r_ellipsoid(2),r_ellipsoid(3))
    vol_sphere = 4./3.*pi*r_ellipsoid(1)*r_ellipsoid(2)*r_ellipsoid(3)
 endif
 !
 ! boundaries
 !
 do i = 1,3
    xmini(i) = 1.2*rmax
    xmaxi(i) = -xmini(i)
 enddo
 print*,'xmini(1)',xmini(1)
 print*,'xmini(2)',xmini(2)
 print*,'xmini(3)',xmini(3)
 call set_boundary(xmini(1),xmaxi(1),xmini(2),xmaxi(2),xmini(3),xmaxi(3))
 !
 ! convert units of sound speed
 !
 if (cs_in_code) then
    cs_sphere_cgs = cs_sphere*unit_velocity
 else
    cs_sphere     = cs_sphere_cgs/unit_velocity
 endif
 !
 ! Bonnor-Ebert profile (if requested)
 !
 if (BEsphere) then
    iBE = 8192
    allocate(rtab(iBE),rhotab(iBE))
    call rho_bonnorebert(iBEparam,BErho_cen,edge_density,BErad_phys,BErad_norm,BEmass,BEfac,cs_sphere, &
                         iBE,iBElast,rtab,rhotab,ierr)
    central_density = BErho_cen
    r_sphere        = BErad_phys
    totmass_sphere  = BEmass
    if (ierr > 0) call fatal('setup_sphere_sphngsne','Error in calculating Bonnor-Ebert profile')
 endif
 !
 ! general parameters
 !
 time        = 0.
 hfact       = hfact_default
 if (maxvxyzu >=4 ) then
    gamma    = 5./3.
 else
    gamma    = 1.
 endif
 angvel_code = angvel*utime
 rhozero     = totmass_sphere / vol_sphere
 dens_sphere = rhozero
 totmass     = totmass_sphere
 t_ff        = sqrt(3.*pi/(32.*dens_sphere))
 !
 ! setup particles in the sphere; use this routine to get N_sphere as close to np as possible
 !
 if (is_sphere) then
    if (BEsphere) then
       call set_sphere(trim(lattice),id,master,0.,r_sphere,psep,hfact,npart,xyzh, &
                       rhotab=rhotab(1:iBElast),rtab=rtab(1:iBElast),nptot=npart_total,&
                       exactN=.true.,np_requested=np,mask=i_belong)
    else
       density_func => gauss_density_func
!       call set_sphere(trim(lattice),id,master,0.,r_sphere,psep,hfact,npart,xyzh,nptot=npart_total,&
!                       rhofunc=density_func,exactN=.true.,np_requested=np,mask=i_belong)
       call set_sphere(trim(lattice),id,master,0.,r_sphere,psep,hfact,npart,xyzh,nptot=npart_total,&
                       exactN=.true.,np_requested=np,mask=i_belong)
       if (trim(lattice)/='random') print "(a,es10.3)",' Particle separation in sphere = ',psep
    endif
    print "(a)",' Initialised sphere'
 else
    call set_ellipse(trim(lattice),id,master,r_ellipsoid,vol_sphere,psep,hfact,xyzh,npart,&
                     nptot=npart_total,exactN=.true.,np_requested=np,mask=i_belong)
    if (trim(lattice)/='random') print "(a,es10.3)",' Particle separation in ellipsoid = ',psep
    print "(a)",' Initialised ellipsoid'
 endif
 !
 ! Slightly randomize particle positions by some factor of smoothing length h
 !
 if (pos_ranh) then
    iseed = -4385
    do i = 1,size(xyzh(1,:))
       xyzh(1,i) = xyzh(1,i) + (-1)**nint(ran2(iseed))*ran2(iseed)*xyzh(4,i)
       xyzh(2,i) = xyzh(2,i) + (-1)**nint(ran2(iseed))*ran2(iseed)*xyzh(4,i)
       xyzh(3,i) = xyzh(3,i) + (-1)**nint(ran2(iseed))*ran2(iseed)*xyzh(4,i)
    enddo
 endif

 npartsphere = npart
 if (np_in /= npartsphere) np = npartsphere
 massoftype(igas) = totmass_sphere/npartsphere
 !
 ! set particle properties
 !
 npartoftype(:)    = 0
 npartoftype(igas) = npart
 if (massoftype(igas) < epsilon(massoftype(igas))) massoftype(igas) = totmass/npart_total
 do i = 1,npartoftype(igas)
    call set_particle_type(i,igas)
 enddo

 if (BEsphere) deallocate(rtab,rhotab)
 !
 ! reset to centre of mass
 ! (if random, recentering may shift particles outside of the defined range)
 !
 if (trim(lattice)/='random') call reset_centreofmass(npart,xyzh,vxyzu)

 !
 ! temperature set to give a pressure equilibrium
 !
 polyk  = cs_sphere**2
 polyk2 = cs_medium**2

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
    ! Convert endian for different vfield files:
    ! setenv GFORTRAN_CONVERT_UNIT big/small_endian

    turbboxsize = max(abs(xmini(1)),abs(xmaxi(1)),abs(xmini(2)),abs(xmaxi(2)),abs(xmini(3)),abs(xmaxi(3)))
    call set_velfield_from_cubes(xyzh(:,1:npartsphere),vxyzu(:,:npartsphere),npartsphere, &
                                 filex,filey,filez,1.,turbboxsize,.false.,ierr)
    if (ierr /= 0) call fatal('setup','error setting up velocity field on clouds')

    rmsmach = 0.0
    print*, 'Turbulence being set by user'
    do i = 1,npartsphere
       v2i     = dot_product(vxyzu(1:3,i),vxyzu(1:3,i))
       rmsmach = rmsmach + v2i/cs_sphere**2
    enddo
    rmsmach = sqrt(rmsmach/npartsphere)
    if (rmsmach > 0.) then
       turbfac = rms_mach/rmsmach ! normalise the energy to the desired mach number
     else
       turbfac = 0.
    endif
    do i = 1,npartsphere
       vxyzu(1:3,i) = turbfac*vxyzu(1:3,i)
    enddo

    vxyz_avg = 0.
    vxyz_min = huge(vxyz_min)
    vxyz_max = epsilon(vxyz_max)
!    open(2022,file='velfield.txt')
    do i = 1,npart
       vxyz_avg = vxyz_avg + mag(vxyzu(1:3,i))
       if (mag(vxyzu(1:3,i)) < vxyz_min) vxyz_min = mag(vxyzu(1:3,i))
       if (mag(vxyzu(1:3,i)) > vxyz_max) vxyz_max = mag(vxyzu(1:3,i))
!       write(2022,*) mag(vxyzu(1:3,i))*unit_velocity*1E-5
    enddo
!    close(2022)
    vxyz_avg = vxyz_avg/npart
    vxyz_avg_cgs = vxyz_avg*unit_velocity
    vxyz_min_cgs = vxyz_min*unit_velocity
    vxyz_max_cgs = vxyz_max*unit_velocity
    print*,'vturb (km/s) mean',vxyz_avg_cgs*1E-5,'min',vxyz_min_cgs*1E-5,'max',vxyz_max_cgs*1E-5

 endif
 !
 ! velocity field corresponding to uniform rotation
 !
 do i=1,npart
    r2 = dot_product(xyzh(1:3,i),xyzh(1:3,i))
    if (r2 < rmax**2) then
       vxyzu(1,i) = vxyzu(1,i) - angvel_code*xyzh(2,i)
       vxyzu(2,i) = vxyzu(2,i) + angvel_code*xyzh(1,i)
       if (maxvxyzu >= 4) vxyzu(4,i) = 1.5*polyk
    else
       if (maxvxyzu >= 4) vxyzu(4,i) = 1.5*polyk2
    endif
 enddo
 !
 ! Set default runtime parameters if .in file does not exist
 !
 infilename=trim(fileprefix)//'.in'
 inquire(file=infilename,exist=in_iexist)
 if (.not. in_iexist) then
    tmax      = 3.15360E14/utime ! 1E1 Myr
    dtmax     = 1.57680E10/utime ! 5E-4 Myr
    nout      = 1
    nfulldump = 1
    nmaxdumps = -1
    dtwallmax = 1800.  ! s
    iverbose  = 1

    ieos      = 2    ! adiabatic eos with P = (gamma-1)*rho*u
    icooling  = 7
    Tfloor    = 3.
    ufloor    = kboltz*Tfloor/(gmw*mass_proton_cgs*(gamma-1.))/unit_ergg
    ipdv_heating   = 1
    ishock_heating = 1

    icreate_sinks   = icreate_sinks_setup
    rho_crit_cgs    = rho_crit_cgs_setup
    rho_crit        = rho_crit_cgs/unit_density
    r_crit          = r_crit_setup
    h_acc           = h_acc_setup
    h_soft_sinksink = h_soft_sinksink_setup

    !
    ! Photoionization settings
    !
    monochrom_source  = .false.
    fix_temp_hii      = .false.
    treat_Rtype_phase = .true.
    old_sources_exist = .true.

    photoionize_tree  = .true.
    tree_accuracy_cmi = 0.2
    nHlimit_fac       = 50
    rcut_opennode_cgs = 4.6E18   ! 1.5 pc
    rcut_leafpart_cgs = 3.1E18   ! 1.0 pc
    delta_rcut_cgs    = 3.1E17   ! 0.1 pc

 endif
 !
 !--Summarise the sphere
 !
 if (is_sphere) then
    print "(a)", ' Sphere'
    print "(a,i10)",' Input number of particles in sphere = ',np
 else
    print "(a)", ' Ellipsoid'
    print "(a,i10)",' Input number of particles in ellipsoid = ',np
 endif
 print "(1x,50('-'))"
 print "(a)",'  Quantity         (code units)  (physical units)'
 print "(1x,50('-'))"
 fmt = "((a,1pg10.3,3x,1pg10.3),a)"
 print fmt,' Total mass       : ',totmass,totmass*umass,' g'
 if (is_sphere) then
    print fmt,' Mass in sphere   : ',totmass_sphere,totmass_sphere*umass,' g'
    print fmt,' Radius of sphere : ',r_sphere,r_sphere*udist,' cm'
 else
    print fmt,' Mass in ellipsoid   : ',totmass_sphere,totmass_sphere*umass,' g'
    print fmt,' Semi-axis {a} of ellipsoid : ',r_ellipsoid(1),r_ellipsoid(1)*udist,' cm'
    print fmt,' Semi-axis {b} of ellipsoid : ',r_ellipsoid(2),r_ellipsoid(2)*udist,' cm'
    print fmt,' Semi-axis {c} of ellipsoid : ',r_ellipsoid(3),r_ellipsoid(3)*udist,' cm'
 endif
 if (BEsphere) then
    print fmt,' Mean rho sphere  : ',dens_sphere,dens_sphere*unit_density,' g/cm^3'
    print fmt,' central density  : ',central_density,central_density*unit_density,' g/cm^3'
    print fmt,' edge density     : ',edge_density,edge_density*unit_density,' g/cm^3'
    print fmt,' Mean rho medium  : ',dens_medium,dens_medium*unit_density,' g/cm^3'
 else
    if (is_sphere) then
       print fmt,' Density sphere   : ',dens_sphere,dens_sphere*unit_density,' g/cm^3'
    else
       print fmt,' Density ellipsoid   : ',dens_sphere,dens_sphere*unit_density,' g/cm^3'
    endif
    print fmt,' Density medium   : ',dens_medium,dens_medium*unit_density,' g/cm^3'
 endif
 if (is_sphere) then
    print fmt,' cs in sphere     : ',cs_sphere,cs_sphere_cgs,' cm/s'
 else
    print fmt,' cs in ellipsoid     : ',cs_sphere,cs_sphere_cgs,' cm/s'
 endif
 print fmt,' cs in medium     : ',cs_medium,cs_medium*unit_velocity,' cm/s'
 print fmt,' Free fall time   : ',t_ff,t_ff*utime/years,' yrs'
 print fmt,' Angular velocity : ',angvel_code,angvel,' rad/s'
 print fmt,' Turbulent Mach no: ',rms_mach
 print fmt,' Omega*t_ff       : ',angvel_code*t_ff
 print "(1x,50('-'))"

end subroutine setpart

!
! Set Gaussian density profile
!
real function gauss_density_func(r)
 use physcon,  only:pi
 real, intent(in) :: r
 real :: amp,sigma

 sigma = 0.67*r_sphere  !- rho at center is 3x of edge
 gauss_density_func = exp(-r**2/(2*sigma**2))

end function gauss_density_func


real function mag(vec)
 real, intent(in) :: vec(3)

 mag = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)

end function mag

!----------------------------------------------------------------
!+
!  write parameters to setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only: write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter           :: iunit = 20
 integer                      :: i

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for sphere setup routines'
 write(iunit,"(/,a)") '# units'
 call write_inopt(dist_unit,'dist_unit','distance unit (e.g. au)',iunit)
 call write_inopt(mass_unit,'mass_unit','mass unit (e.g. solarm)',iunit)
 write(iunit,"(/,a)") '# resolution'

 write(iunit,"(/,a)") '# options for spherical shape'
 call write_inopt(is_sphere,'is_sphere','the intent is to set up a sphere',iunit)

 if (is_sphere) then
    write(iunit,"(/,a)") '# options for sphere'
    call write_inopt(np,'np','requested number of particles in sphere',iunit)
    call write_inopt(BEsphere,'use_BE_sphere','centrally condense as a BE sphere',iunit)
    if (.not. BEsphere) then
       call write_inopt(r_sphere,'r_sphere','radius of sphere in code units',iunit)
       call write_inopt(totmass_sphere,'totmass_sphere','mass of sphere in code units',iunit)
    else
       call write_inopt(iBEparam,'iBE_options','The set of parameters to define the BE sphere',iunit)
       if (iBEparam==1 .or. iBEparam==2 .or. iBEparam==3) &
          call write_inopt(BErho_cen,'BErho_cen','central density of the BE sphere [code units]',iunit)
       if (iBEparam==1 .or. iBEparam==4 .or. iBEparam==6) &
           call write_inopt(BErad_phys,'BErad_phys','physical radius of the BE sphere [code units]',iunit)
       if (iBEparam==2 .or. iBEparam==4 .or. iBEparam==5) &
           call write_inopt(BErad_norm,'BErad_norm','normalised radius of the BE sphere',iunit)
       if (iBEparam==3 .or. iBEparam==5 .or. iBEparam==6) &
           call write_inopt(BEmass,'BEmass','mass radius of the BE sphere [code units]',iunit)
       if (iBEparam==4 .or. iBEparam==5)                  &
           call write_inopt(BEfac,'BEfac','over-density factor of the BE sphere [code units]',iunit)
    endif
 else
    write(iunit,"(/,a)") '# options for ellipsoid'
    call write_inopt(np,'np','requested number of particles in ellipsoid',iunit)
    call write_inopt(r_ellipsoid(1),'r_ellipsoid(1)','semi-axis {a} of ellipsoid in code units',iunit)
    call write_inopt(r_ellipsoid(2),'r_ellipsoid(2)','semi-axis {b} of ellipsoid in code units',iunit)
    call write_inopt(r_ellipsoid(3),'r_ellipsoid(3)','semi-axis {c} of ellipsoid in code units',iunit)
    call write_inopt(totmass_sphere,'totmass_sphere','mass of ellipsoid in code units',iunit)
 endif

 write(iunit,"(/,a)") '# intended result'
 call write_inopt(lattice,'lattice','particle lattice',iunit)
 call write_inopt(mc_method,'mc_method','the intent is to place particles with Monte Carlo method',iunit)
 if (.not.mc_method .and. lattice/='random') then
    call write_inopt(pos_ranh,'pos_ranh','the intent is to slightly randomize particle positions on lattice',iunit)
 endif
 if (is_sphere) then
    call write_inopt(cs_sphere_cgs,'cs_sphere_cgs','sound speed in sphere in cm/s',iunit)
 else
    call write_inopt(cs_sphere_cgs,'cs_sphere_cgs','sound speed in ellipsoid in cm/s',iunit)
 endif
 call write_inopt(angvel,'angvel','angular velocity in rad/s',iunit)
 call write_inopt(rms_mach,'rms_mach','turbulent rms mach number',iunit)
 write(iunit,"(/,a)") '# Sink properties (values in .in file, if present, will take precedence)'
 call write_inopt(icreate_sinks_setup,'icreate_sinks','1: create sinks.  0: do not create sinks',iunit)
 if (icreate_sinks_setup==1) then
    call write_inopt(rho_crit_cgs_setup,'rho_crit_cgs','critical density in cgs units',iunit)
    call write_inopt(h_acc_setup,'h_acc','accretion radius (code units)',iunit)
    call write_inopt(r_crit_setup,'r_crit','critical radius (code units)',iunit)
 else
    call write_inopt(rhofinal_setup,'rho_final','final maximum density (<=0 to ignore) (cgs units)',iunit)
 endif
 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  Read parameters from setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only: open_db_from_file,inopts,read_inopt,close_db
 use unifdis,      only: is_valid_lattice
 use io,           only: error
 use units,        only: select_unit
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit = 21
 integer                       :: i,nerr,jerr,kerr
 type(inopts), allocatable     :: db(:)

 !--Read values
 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(mass_unit,'mass_unit',db,ierr)
 call read_inopt(dist_unit,'dist_unit',db,ierr)
 call read_inopt(np,'np',db,ierr)
 call read_inopt(lattice,'lattice',db,ierr)
 call read_inopt(mc_method,'mc_method',db,ierr)
 if (.not. mc_method .and. trim(lattice)/='random') then
    call read_inopt(pos_ranh,'randomize_pos',db,ierr)
 endif
 if ((trim(lattice)=='random' .and. .not.mc_method) .or. &
     (trim(lattice)/='random' .and. mc_method)) then
    print*, ' Monte Carlo method will be used if lattice is set to random - Check input consistency'
    stop
 endif
 if (ierr/=0 .or. .not. is_valid_lattice(trim(lattice))) then
    print*, ' invalid lattice.  Setting to closepacked'
    lattice = 'closepacked'
    mc_method = .false.
 endif

 ! Sphere/Ellipsoid params
 call read_inopt(is_sphere,'is_sphere',db,ierr)
 if (is_sphere) then
    call read_inopt(BEsphere,'use_BE_sphere',db,ierr)
 endif
 if (is_sphere) then
    if (.not. BEsphere) then
       call read_inopt(r_sphere,'r_sphere',db,ierr)
       call read_inopt(totmass_sphere,'totmass_sphere',db,ierr)
    else
       call read_inopt(iBEparam,'iBE_options',db,ierr)
       if (iBEparam==1 .or. iBEparam==2 .or. iBEparam==3) call read_inopt(BErho_cen,'BErho_cen',db,ierr)
       if (iBEparam==1 .or. iBEparam==4 .or. iBEparam==6) call read_inopt(BErad_phys,'BErad_phys',db,ierr)
       if (iBEparam==2 .or. iBEparam==4 .or. iBEparam==5) call read_inopt(BErad_norm,'BErad_norm',db,ierr)
       if (iBEparam==3 .or. iBEparam==5 .or. iBEparam==6) call read_inopt(BEmass,'BEmass',db,ierr)
       if (iBEparam==4 .or. iBEparam==5)                  call read_inopt(BEfac,'BEfac',db,ierr)
    endif
 else
    call read_inopt(r_ellipsoid(1),'r_ellipsoid(1)',db,ierr)
    call read_inopt(r_ellipsoid(2),'r_ellipsoid(2)',db,ierr)
    call read_inopt(r_ellipsoid(3),'r_ellipsoid(3)',db,ierr)
    call read_inopt(totmass_sphere,'totmass_sphere',db,ierr)
 endif

 call read_inopt(cs_sphere,'cs_sphere',db,jerr)
 call read_inopt(cs_sphere_cgs,'cs_sphere_cgs',db,kerr)
 cs_in_code = .false.  ! for backwards compatibility
 if (jerr /= 0 .and. kerr == 0) then
    cs_in_code = .false.
 elseif (jerr == 0 .and. kerr /= 0) then
    cs_in_code = .true.
 else
    ierr = ierr + 1
 endif
 call read_inopt(angvel,'angvel',db,ierr)
 call read_inopt(rms_mach,'rms_mach',db,ierr)
 mu_not_B = .true.
 call read_inopt(icreate_sinks_setup,'icreate_sinks',db,ierr)
 if (icreate_sinks_setup==1) then
    call read_inopt(rho_crit_cgs_setup,'rho_crit_cgs',db,ierr)
    call read_inopt(h_acc_setup,'h_acc',db,ierr)
    call read_inopt(r_crit_setup,'r_crit',db,ierr)
 else
    call read_inopt(rhofinal_setup,'rho_final',db,ierr)
 endif
 call close_db(db)
 !
 ! parse units
 !
 call select_unit(mass_unit,umass,nerr)
 if (nerr /= 0) then
    call error('setup_sphere_sphngsne','mass unit not recognised')
    ierr = ierr + 1
 endif
 call select_unit(dist_unit,udist,nerr)
 if (nerr /= 0) then
    call error('setup_sphere_sphngsne','length unit not recognised')
    ierr = ierr + 1
 endif

 if (ierr > 0) then
    print "(1x,a,i2,a)",'setup_sphere_sphngsne: ',nerr,' error(s) during read of setup file.  Re-writing.'
 endif

end subroutine read_setupfile

end module setup
