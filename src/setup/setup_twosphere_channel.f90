!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setup
!
! This module sets up an ellipsoidal envalope that contains two spherical clouds
! with the option to manually create/mimic a radiation-carved channel
!
! :References: None
!
! :Owner: Cheryl Lau
!
! :Runtime parameters:
!   - dist_unit        : *distance unit (e.g. au)*
!   - mass_unit        : *mass unit (e.g. solarm)*
!   - mpart_solarm     : *mass of particle in solarm units*
!   - np_cloud1        : *requested number of particles in cloud1*
!   - np_cloud2        : *requested number of particles in cloud2*
!   - np_envelope      : *requested number of particles within envelope boundaries*
!   - rho_cloud1_cgs   : *density of cloud1*
!   - rho_cloud2_cgs   : *density of cloud2*
!   - rho_envelope_cgs : *density of envelope*
!   - rms_mach_cloud1  : *Mach number for turbulence in cloud1*
!   - rms_mach_cloud2  : *Mach number for turbulence in cloud2*
!   - cloud_sep_pc     : *separation between the two clouds in pc*
!   - r1_envelope      : *ratio of semi-axis a of envelope*
!   - r2_envelope      : *ratio of semi-axis b of envelope*
!   - r3_envelope      : *ratio of semi-axis c of envelope*
!   - omega_channel    : *solid angle of the channel to carve upon request*
!   - pfrac_channel    : *fraction of particles to remove in the channel*
!   - gmw_in           : *mean molecular weight*
!   - make_sinks       : *flag to dynamically create sinks*
!   - create_channel   : *flag to manually open a channel*
!
! :Dependencies: dim, io, physcon, setup_params, spherical, boundary, prompting, units,
!                eos, part, ptmass, timestep, kernel, options, datafiles,
!                photoionize_cmi, inject
!
 implicit none
 public :: setpart

 private

 integer      :: np_cloud1,np_cloud2,np_envelope
 integer      :: iseed = -12345
 real(kind=8) :: udist,umass
 real         :: mpart_solarm,rho_cloud1_cgs,rho_cloud2_cgs,rho_envelope_cgs,gmw_in
 real         :: rms_mach_cloud1,rms_mach_cloud2,r1_envelope,r2_envelope,r3_envelope
 real         :: cloud_sep_pc,omega_channel,pfrac_channel
 logical      :: make_sinks,create_channel
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
 use io,           only:master,fatal,warning,iprint,iverbose
 use spherical,    only:set_sphere,set_ellipse
 use boundary,     only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use prompting,    only:prompt
 use units,        only:set_units,select_unit,utime,unit_density,unit_velocity,unit_energ,unit_ergg,unit_pressure
 use eos,          only:polyk2,ieos,gmw
 use part,         only:igas,abundance,set_particle_type
 use part,         only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft
 use ptmass,       only:icreate_sinks,rho_crit_cgs,r_crit,h_acc,h_soft_sinksink,h_soft_sinkgas
 use timestep,     only:dtmax,tmax,dtwallmax,C_cour,C_force,C_cool,tolv,nout
 use options,      only:nfulldump,nmaxdumps,icooling,alpha,alphau
 use kernel,       only:hfact_default
 use domain,       only:i_belong
 use random,       only:ran2
 use cooling,      only:Tfloor,ufloor
 use velfield,     only:set_velfield_from_cubes
 use datafiles,    only:find_phantom_datafile
 use options,      only:ipdv_heating,ishock_heating
 use photoionize_cmi, only:monochrom_source,fix_temp_hii,treat_Rtype_phase
 use photoionize_cmi, only:photoionize_tree,tree_accuracy_cmi,nHlimit_fac
 use photoionize_cmi, only:rcut_opennode_cgs,rcut_leafpart_cgs,delta_rcut_cgs
 use photoionize_cmi, only:sink_ionsrc,sink_as_cluster,inject_rad,one_sink_ionsrc,isink_ionsrc
 use inject,          only:inject_sn,sink_progenitor,frackin,fractherm 
 use inject,          only:one_sink_progenitor,isink_progenitor 
 use inject,          only:delay_sn_injection,delay_by_mslifetime
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
 integer            :: ip,npmax,ierr,i
 integer            :: npart_cloud1,npart_cloud2,np_outer,npart_outer,npmax_env,xyzh_size_tmp
 integer            :: npart_tmp_cloud1,npart_tmp_cloud2,npart_tmp_outer,npart_diff
 integer(kind=8)    :: npart_total_cloud1,npart_total_cloud2,npart_total_outer
 real,  allocatable :: xyzh_tmp_cloud1(:,:),xyzh_tmp_cloud2(:,:),xyzh_tmp_outer(:,:)
 real,  allocatable :: vxyzu_turb(:,:)
 real(kind=8)       :: h_accs_in
 real               :: r_cloud1,cen_cloud1(3),psep_cloud1,cs_cloud1,cs_cloud1_cgs,rho_cloud1,p_cloud1,p_cloud1_cgs
 real               :: temp_cloud1,vol_cloud1,t_ff_cloud1,totmass_cloud1,u_cloud1
 real               :: rmsmach,v2i,turbfac,turbboxsize,v2_sum,v_rms_cloud1,v_rms_kms_cloud1
 real               :: r_cloud2,cen_cloud2(3),psep_cloud2,cs_cloud2,cs_cloud2_cgs,rho_cloud2
 real               :: temp_cloud2,vol_cloud2,t_ff_cloud2,totmass_cloud2,u_cloud2
 real               :: v_rms_cloud2,v_rms_kms_cloud2
 real               :: mjeans_cloud2,mjeans_cgs_cloud2
 real               :: cloud_sep,minx_clouds,minyz_clouds
 real               :: omega,area,radius,rad_circ,x,y,z,rad_strom,rad_stag,rad_strom_cgs,u_hii,cs_0,cs_i
 real               :: r_outer(3),psep_outer,scale_param,rho_outer,vol_outer,totmass_outer
 real               :: temp_envelope,cs_envelope_cgs,cs_envelope,u_envelope,p_envelope,p_envelope_cgs
 real               :: r_sn_cgs,engsn_cgs,pmsncrit_cgs
 real               :: h_acc_cgs,h_soft_sinksink_cgs,h_soft_sinkgas_cgs,rho_crit_cgs_recomm
 logical            :: iexist,in_iexist,add_particle
 logical            :: remove_cloud1       = .true. ! temporarily removing objects to check virial ratio
 logical            :: remove_cloud2       = .true.
 logical            :: remove_envelope     = .false.
 logical            :: place_sink_in_setup = .true.
 character(len=120) :: filex,filey,filez
 character(len=100) :: filename,infilename
 character(len=40)  :: fmt,lattice_cloud1,lattice_cloud2,lattice_envelope
 character(len=9)   :: proceed

 print "(a)", 'Setup for two spheres with an ellipsoidal envelope'

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
    npmax = int(size(xyzh(1,:))) ! approx max number allowed in sphere given size(xyzh(1,:))

    !- Mass of particles
    mpart_solarm = 1E-2
    call prompt('Enter the mass of particles in solarm units',mpart_solarm,0.)

    !- Settings for the cloud to inject feedback (cloud1)
    np_cloud1 = 1E7
    call prompt('Enter the approximate number of particles in cloud1 (with feedback injected)',np_cloud1,0,npmax)
    rho_cloud1_cgs = 1E-21
    call prompt('Enter the density of cloud1 in g/cm^3',rho_cloud1_cgs,0.)
    rms_mach_cloud1 = 0.
    call prompt('Enter the Mach number of the cloud1 turbulence',rms_mach_cloud1,0.)

    !- Settings for the neighbouring cloud (cloud2)
    np_cloud2 = 1E5
    call prompt('Enter the approximate number of particles in cloud2',np_cloud2,0,npmax)
    rho_cloud2_cgs = 1E-21
    call prompt('Enter the density of cloud2 in g/cm^3',rho_cloud2_cgs,0.)
    rms_mach_cloud2 = 6.
    call prompt('Enter the Mach number of the cloud2 turbulence',rms_mach_cloud2,0.)

    !- Settings for cloud locations
    cloud_sep_pc = 30.
    call prompt('Enter the separation between the centre of clouds in pc',cloud_sep_pc,0.)

    !- Settings for the envelope
    np_envelope = 1E6 !3E5
    call prompt('Enter the approximate number of particles within the envelope boundaries',np_envelope,0,npmax)
    rho_envelope_cgs = 4E-25  ! for 1000K  
    call prompt('Enter the density of the envelope in g/cm^3',rho_envelope_cgs,0.)

    !- Ratio of semi-axes of ellipsoidal envelope
    r1_envelope = 1.
    r2_envelope = 1.
    r3_envelope = 1.
    call prompt('Enter the ratio of semi-axis a of the ellipsoidal envelope ',r1_envelope,0.)
    call prompt('Enter the ratio of semi-axis b of the ellipsoidal envelope ',r2_envelope,0.)
    call prompt('Enter the ratio of semi-axis c of the ellipsoidal envelope ',r3_envelope,0.)

    !- Mean molecular weight
    gmw_in = 1.29
    call prompt('Enter the mean molecular weight',gmw_in,0.)

    make_sinks = .false.
    call prompt('Do you wish to dynamically create sink particles? ',make_sinks)

    create_channel = .false.
    call prompt('Do you wish to manually carve a channel? ',create_channel)

    if (create_channel) then
       omega_channel = 0.1 * 4.*pi
       call prompt('Enter the solid angle of the channel in steradian ',omega_channel,0.,4.*pi)
       pfrac_channel = 0.8
       call prompt('Enter the fraction of particles to remove in the channel ',pfrac_channel,0.,1.)
    endif

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
 massoftype(igas) = mpart_solarm*solarm/umass
 time  = 0.
 hfact = hfact_default

 !
 ! Set up cloud1 - cloud at origin with feedback injected
 !
 totmass_cloud1 = np_cloud1*massoftype(igas)
 rho_cloud1     = rho_cloud1_cgs/unit_density
 vol_cloud1     = totmass_cloud1/rho_cloud1
 r_cloud1       = (3./4.*vol_cloud1/pi)**(1./3.)
 t_ff_cloud1    = sqrt(3.*pi/(32.*rho_cloud1))
 lattice_cloud1 = 'random'
 cen_cloud1 = (/ 0.,0.,0. /)

 if (.not.remove_cloud1) then

    !- Set sphere particles
    npart_cloud1 = 0
    npart_tmp_cloud1 = 0
    npart_total_cloud1 = 0

    xyzh_size_tmp = int(np_cloud1*1.5)
    allocate(xyzh_tmp_cloud1(4,xyzh_size_tmp))

    call set_sphere(trim(lattice_cloud1),id,master,0.,r_cloud1,psep_cloud1,hfact,npart_tmp_cloud1,xyzh_tmp_cloud1,&
                    nptot=npart_total_cloud1,exactN=.true.,np_requested=np_cloud1,mask=i_belong)

    print*,'Number of particles in cloud1: ',npart_tmp_cloud1

    !- Manually clear away some particles in a channel that points towards the neighbouring cloud
    if (create_channel) then
       check_particles: do ip = 1,npart_tmp_cloud1
          x = xyzh_tmp_cloud1(1,ip)
          y = xyzh_tmp_cloud1(2,ip)
          z = xyzh_tmp_cloud1(3,ip)
          !- Remove a portion of particles that fall within the user-defined solid angle omega_channel
          add_particle = .true.
          if (x > 0.) then
             radius   = x
             rad_circ = (abs(y)**2 + abs(z)**2)**(1./2.)
             area     = 4.*pi*rad_circ**2
             omega    = area/radius**2
             if (omega < omega_channel) then
                if (ran2(iseed) < pfrac_channel) add_particle = .false.
             endif
          endif
          if (add_particle) then  !- store particle
             npart = npart + 1
             xyzh(1:4,npart) = xyzh_tmp_cloud1(1:4,ip)
             npart_cloud1 = npart_cloud1 + 1   !- bookkeeping
          endif
       enddo check_particles

       print*,'Number of particles in cloud1 after carving channel: ',npart
       if (npart == npart_tmp_cloud1) call fatal('setup_twosphere_channel','no particles removed from channel')
    else
       npart = npart_tmp_cloud1
       do ip = 1,npart
          xyzh(:,ip) = xyzh_tmp_cloud1(:,ip)
       enddo
       npart_cloud1 = npart  !- bookkeeping
    endif

    deallocate(xyzh_tmp_cloud1)

    !- Temperature and sound speed
    temp_cloud1   = get_eqtemp_from_rho(rho_cloud1_cgs)
    cs_cloud1_cgs = sqrt(temp_cloud1*(gamma*kboltz)/(gmw*mass_proton_cgs))
    cs_cloud1     = cs_cloud1_cgs/unit_velocity

    !- Pressure
    p_cloud1     = cs_cloud1**2*rho_cloud1/gamma
    p_cloud1_cgs = p_cloud1*unit_pressure

    !- Apply turbulence
    vxyzu = 0.
    if (rms_mach_cloud1 > 0.) then

       turbboxsize = 1.1*r_cloud1
       filex  = find_phantom_datafile(filevx,'velfield_sphng_small')
       filey  = find_phantom_datafile(filevy,'velfield_sphng_small')
       filez  = find_phantom_datafile(filevz,'velfield_sphng_small')
       call set_velfield_from_cubes(xyzh(:,1:npart),vxyzu(:,1:npart),npart_cloud1, &
                                    filex,filey,filez,1.,turbboxsize,.false.,ierr)
       if (ierr /= 0) call fatal('setup','error setting up velocity field on clouds')

       rmsmach = 0.0
       do ip = 1,npart
          v2i     = dot_product(vxyzu(1:3,ip),vxyzu(1:3,ip))
          rmsmach = rmsmach + v2i/cs_cloud1**2
       enddo
       rmsmach = sqrt(rmsmach/npart_cloud1)
       if (rmsmach > 0.) then
          turbfac = rms_mach_cloud1/rmsmach ! normalise the energy to the desired mach number
       else
          turbfac = 0.
       endif
       do ip = 1,npart
          vxyzu(1:3,ip) = turbfac*vxyzu(1:3,ip)
       enddo

       !- Calculate rms-velocity
       v2_sum = 0.
       do ip = 1,npart
          v2_sum = v2_sum + mag2(vxyzu(1:3,ip))
       enddo
       v_rms_cloud1 = sqrt(1./npart_cloud1*v2_sum)
       v_rms_kms_cloud1 = v_rms_cloud1*unit_velocity*1E-5
    endif

    !- Set initial temperature
    u_cloud1 = kboltz * temp_cloud1 / (gmw*mass_proton_cgs*(gamma-1.)) /unit_ergg
    do ip = 1,npart
      vxyzu(4,ip) = u_cloud1
    enddo

    !- Estimate Stromgren radius
    rad_strom_cgs = (3.*1E49*mass_proton_cgs**2/(4.*pi*2.7E-13*rho_cloud1_cgs**2))**(1./3.)
    rad_strom = rad_strom_cgs/udist

    !- Stagnation radius of HII region
    cs_i = sqrt(1E4*(gamma*kboltz)/(gmw*mass_proton_cgs)) / unit_velocity  ! cs of ionized medium
    cs_0 = cs_cloud1                                                       ! cs of neutral medium

    rad_stag = (8./3.)**(2./3.) * (cs_i/cs_0)**(4./3.) * rad_strom
    if (rad_stag < 0.2*r_cloud1) call warning('setup_twosphere_channel','HII region could be too small')
    if (rad_stag > r_cloud1) call warning('setup_twosphere_channel','HII region could be beyond the sphere boundaries')
 endif

 !
 ! Set up cloud2 - neighbouring cloud
 !
 totmass_cloud2 = np_cloud2*massoftype(igas)
 rho_cloud2     = rho_cloud2_cgs/unit_density
 vol_cloud2     = totmass_cloud2/rho_cloud2
 r_cloud2       = (3./4.*vol_cloud2/pi)**(1./3.)
 t_ff_cloud2    = sqrt(3.*pi/(32.*rho_cloud2))
 lattice_cloud2 = 'closepacked'
 cloud_sep = cloud_sep_pc*pc/udist
 cen_cloud2     = (/ cloud_sep,0.,0. /)

 if (cloud_sep < r_cloud1+r_cloud2) call warning('setup_twosphere_channel','spheres could overlap')

 if (.not.remove_cloud2) then

    !- Set sphere particles
    npart_cloud2 = 0
    npart_tmp_cloud2 = 0
    npart_total_cloud2 = 0

    xyzh_size_tmp = int(np_cloud2*1.5)
    allocate(xyzh_tmp_cloud2(4,xyzh_size_tmp))

    call set_sphere(trim(lattice_cloud2),id,master,0.,r_cloud2,psep_cloud2,hfact,npart_tmp_cloud2,xyzh_tmp_cloud2,&
                    nptot=npart_total_cloud2,exactN=.true.,np_requested=np_cloud2,mask=i_belong)

    npart_cloud2 = npart_tmp_cloud2

    !- Shift position of sphere down x-axis
    do ip = 1,npart_cloud2
       npart = npart + 1
       xyzh(:,npart) = xyzh_tmp_cloud2(:,ip)
       xyzh(1,npart) = xyzh(1,npart) + cloud_sep
    enddo

    print*,'Number of particles in cloud2: ',npart_cloud2
    print*,'Total number of particles in both clouds: ',npart

    !- Temperature and sound speed
    temp_cloud2   = get_eqtemp_from_rho(rho_cloud2_cgs)
    cs_cloud2_cgs = sqrt(temp_cloud2*(gamma*kboltz)/(gmw*mass_proton_cgs))
    cs_cloud2     = cs_cloud2_cgs/unit_velocity

    !- Apply turbulence
    allocate(vxyzu_turb(4,npart_cloud2))
    do ip = npart_cloud1+1,npart
      vxyzu(1:3,ip) = 0.  ! init
    enddo
    if (rms_mach_cloud2 > 0.) then
       turbboxsize = 1.1*r_cloud2
       filex  = find_phantom_datafile(filevx,'velfield_sphng_small')
       filey  = find_phantom_datafile(filevy,'velfield_sphng_small')
       filez  = find_phantom_datafile(filevz,'velfield_sphng_small')

       call set_velfield_from_cubes(xyzh_tmp_cloud2(:,:),vxyzu_turb(:,:),npart_cloud2, &
            filex,filey,filez,1.,turbboxsize,.false.,ierr)
       if (ierr /= 0) call fatal('setup','error setting up velocity field on clouds')

       rmsmach = 0.
       do ip = 1,npart_cloud2
          v2i     = dot_product(vxyzu_turb(1:3,ip),vxyzu_turb(1:3,ip))
          rmsmach = rmsmach + v2i/cs_cloud2**2
       enddo
       rmsmach = sqrt(rmsmach/npart_cloud2)
       if (rmsmach > 0.) then
          turbfac = rms_mach_cloud2/rmsmach ! normalise the energy to the desired mach number
       else
          turbfac = 0.
       endif
       do ip = 1,npart_cloud2
          vxyzu(1:3,npart_cloud1+ip) = turbfac*vxyzu_turb(1:3,ip)
       enddo

       !- Calculate rms-velocity
       v2_sum = 0.
       do ip = npart_cloud1+1,npart
          v2_sum = v2_sum + mag2(vxyzu(1:3,ip))
       enddo
       v_rms_cloud2 = sqrt(1./npart_cloud2*v2_sum)
       v_rms_kms_cloud2 = v_rms_cloud2*unit_velocity*1E-5
    endif
    deallocate(vxyzu_turb)

    !- Set initial temperature
    u_cloud2 = kboltz * temp_cloud2 / (gmw*mass_proton_cgs*(gamma-1.)) /unit_ergg
    do ip = npart_cloud1+1,npart
      vxyzu(4,ip) = u_cloud2
    enddo

    deallocate(xyzh_tmp_cloud2)
 endif

 !
 ! Set up an outer ellipsoid then remove those within the clouds to create the envelope
 !
 if (.not.remove_envelope) then
    np_outer = np_envelope
    totmass_outer = np_outer*massoftype(igas)
    rho_outer = rho_envelope_cgs/unit_density
    vol_outer = totmass_outer/rho_outer
    scale_param = (vol_outer/(4./3.*pi*r1_envelope*r2_envelope*r3_envelope))**(1./3.)
    r_outer(1) = r1_envelope*scale_param
    r_outer(2) = r2_envelope*scale_param
    r_outer(3) = r3_envelope*scale_param
    lattice_envelope = 'closepacked'

    !- Check that the size of envelope can cover the two clouds
    minx_clouds = r_cloud1
    if (.not.remove_cloud2) minx_clouds = minx_clouds + cloud_sep + r_cloud2
    minyz_clouds = r_cloud1
    if (.not.remove_cloud2) minyz_clouds = 2.*max(r_cloud1,r_cloud2)
    if (2.*r_outer(1) < 1.5*minx_clouds .or. 2.*max(r_outer(2),r_outer(3)) < 1.5*minyz_clouds) then
       call fatal('setup_twosphere_channel','require more particles in envelope')
    endif

    xyzh_size_tmp = int(np_outer*1.5)
    npart_tmp_outer = 0
    npart_total_outer = 0
    allocate(xyzh_tmp_outer(4,xyzh_size_tmp))

    call set_ellipse(trim(lattice_envelope),id,master,r_outer,vol_outer,psep_outer,hfact,xyzh_tmp_outer,npart_tmp_outer,&
                     nptot=npart_total_outer,exactN=.true.,np_requested=np_outer,mask=i_belong)

    !- Update parameters
    totmass_outer = npart_tmp_outer*massoftype(igas)
    rho_outer = totmass_outer/vol_outer
    rho_envelope_cgs = rho_outer*unit_density

    !- Store particle only if it is outside both clouds
    do ip = 1,npart_tmp_outer
       if (.not.remove_cloud1 .and. .not.remove_cloud2) then
          if (mag2(xyzh_tmp_outer(1:3,ip)-cen_cloud1) > (1.01*r_cloud1)**2 .and. &
              mag2(xyzh_tmp_outer(1:3,ip)-cen_cloud2) > (1.01*r_cloud2)**2) then
             npart = npart + 1
             xyzh(:,npart) = xyzh_tmp_outer(:,ip)
          endif
       elseif (.not.remove_cloud1) then
          if (mag2(xyzh_tmp_outer(1:3,ip)-cen_cloud1) > (1.01*r_cloud1)**2) then
             npart = npart + 1
             xyzh(:,npart) = xyzh_tmp_outer(:,ip)
          endif
       elseif (.not.remove_cloud2) then
          if (mag2(xyzh_tmp_outer(1:3,ip)-cen_cloud2) > (1.01*r_cloud2)**2) then
             npart = npart + 1
             xyzh(:,npart) = xyzh_tmp_outer(:,ip)
          endif
       else
          npart = npart + 1
          xyzh(:,npart) = xyzh_tmp_outer(:,ip)
       endif
    enddo
    deallocate(xyzh_tmp_outer)

    npart_outer = npart - npart_cloud1 - npart_cloud2

    !- Temperature and sound speed
    temp_envelope   = get_eqtemp_from_rho(rho_envelope_cgs)
    cs_envelope_cgs = sqrt(temp_envelope*(gamma*kboltz)/(gmw*mass_proton_cgs))
    cs_envelope     = cs_envelope_cgs/unit_velocity

    !- Set initial temp
    u_envelope = kboltz * temp_envelope / (gmw*mass_proton_cgs*(gamma-1.)) /unit_ergg
    do ip = npart_cloud1+npart_cloud2+1,npart
       vxyzu(4,ip) = u_envelope
    enddo

    !- Pressure
    p_envelope     = cs_envelope**2*(rho_envelope_cgs/unit_density)/gamma
    p_envelope_cgs = p_envelope*unit_pressure

    npart_diff = npart_cloud1 + npart_cloud2 + npart_tmp_outer - npart
    if (npart_diff > 0) then
       print*,'Removed ',npart_diff,' particles from envelope that overlap the spheres'
       print*,'Number of particles in the envelope: ',npart_outer
    else
       if (.not.remove_cloud1 .and. .not.remove_cloud2) then
          call fatal('setup_twosphere_channel','no particles removed in evelope')
       endif
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
 npart_total = npart

 !
 ! Manually place a sink as feedback source
 !
 if (place_sink_in_setup) then
    nptmass                      = 1
    xyzmh_ptmass(:,:)            = 0.
    xyzmh_ptmass(1:3,nptmass)    = 0.
    xyzmh_ptmass(4,nptmass)      = 50*solarm/umass
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
    tmax      = 2.0*min(t_ff_cloud1,t_ff_cloud2)
    dtmax     = 0.001*min(t_ff_cloud1,t_ff_cloud2)
    nout      = 1
    nfulldump = 1
    nmaxdumps = 1000
    dtwallmax = 1800.  ! s
    iverbose  = 0

    ieos      = 2    ! adiabatic eos with P = (gamma-1)*rho*u
    gmw       = gmw_in
    icooling  = 7
    Tfloor    = 3.
    ufloor    = kboltz*Tfloor/(gmw*mass_proton_cgs*(gamma-1.))/unit_ergg
    ipdv_heating   = 1
    ishock_heating = 1
    alphau = 1.

    !
    ! Sinks settings
    !
    if (make_sinks) then
       icreate_sinks = 1
       rho_crit_cgs = 1.E-16           ! density above which sink particles are created
       h_acc_cgs    = 0.005*pc         ! accretion radius for new sink particles
       h_soft_sinksink_cgs = 0.005*pc  ! softening length between sink particles
       h_soft_sinkgas_cgs  = 0.        ! softening length for new sink particles

       !- Check if rho_crit is sensible [10^5 times the initial MC density (Bate et al 1995)]
       rho_crit_cgs_recomm = 1.d5*max(rho_cloud1_cgs,rho_cloud2_cgs)
       if (abs(log10(rho_crit_cgs_recomm)-log10(rho_crit_cgs))/log10(rho_crit_cgs_recomm) > 0.1) then
          print*,'Recommend setting rho_crit_cgs to ',rho_crit_cgs_recomm,' instead'
          proceed = 'n'
          call prompt('Do you wish to continue?',proceed)
          if (trim(adjustl(proceed)) == 'n') then
             stop 'Edit rho_crit_cgs and rerun'
          elseif (trim(adjustl(proceed)) /= 'y') then
             stop 'Invalid input'
          endif
       endif

       !- convert to code units
       h_acc           = h_acc_cgs/udist
       r_crit          = 2.*h_acc
       h_soft_sinksink = h_soft_sinksink_cgs/udist
       h_soft_sinkgas  = h_soft_sinkgas_cgs/udist

    else
       icreate_sinks = 0
    endif

    !
    ! Photoionization settings
    !
    inject_rad  = .false.
    sink_ionsrc = .true. 
    one_sink_ionsrc = .true. 
    isink_ionsrc    = 1
    sink_as_cluster = .false. 

    monochrom_source  = .false.
    fix_temp_hii      = .false.
    treat_Rtype_phase = .false.

    photoionize_tree  = .true.
    tree_accuracy_cmi = 0.3
    nHlimit_fac       = 300.
    rcut_opennode_cgs = 3.*pc
    rcut_leafpart_cgs = 2.*pc
    delta_rcut_cgs    = 0.5*pc

    !
    ! Supernova settings
    !
    inject_sn = .false.
    sink_progenitor = .true.
    one_sink_progenitor = .true.
    isink_progenitor = 1
    delay_sn_injection = .false. 
    delay_by_mslifetime = .false. 
    frackin = 1.
    fractherm = 0.

 endif

 !
 ! Calculate the approx number of stars that will spontaneously form in cloud2
 !
 mjeans_cgs_cloud2 = (5.*Rg*temp_cloud2/(2.*gg*gmw))**(3./2.) * (4./3.*pi*rho_cloud2_cgs)**(-1./2.)
 mjeans_cloud2 = mjeans_cgs_cloud2/umass

 if (.not.remove_cloud1) then
    print*,'-Cloud 1-'
    print*,'total mass        ',totmass_cloud1,mass_unit
    print*,'radius            ',r_cloud1,dist_unit
    print*,'free-fall time    ',t_ff_cloud1*utime/(1E6*years),'Myr'
    print*,'temperature       ',temp_cloud1,'K'
    print*,'pressure          ',p_cloud1_cgs*0.1,'Pa'
    print*,'sound speed       ',cs_cloud1_cgs*1E-5,'km/s'
    print*,'v_rms             ',v_rms_kms_cloud1,'km/s'
    print*,'expected R_stag   ',rad_stag,dist_unit
 endif

 if (.not.remove_cloud2) then
    print*,'-Cloud 2-'
    print*,'total mass        ',totmass_cloud2,mass_unit
    print*,'Jeans mass        ',mjeans_cloud2,mass_unit
    print*,'radius            ',r_cloud2,dist_unit
    print*,'free-fall time    ',t_ff_cloud2*utime/(1E6*years),'Myr'
    print*,'temperature       ',temp_cloud2,'K'
    print*,'sound speed       ',cs_cloud2_cgs*1E-5,'km/s'
    print*,'v_rms             ',v_rms_kms_cloud2,'km/s'
 endif

 if (.not.remove_envelope) then
    print*,'-Envelope-'
    print*,'temperature       ',temp_envelope,'K'
    print*,'pressure          ',p_envelope_cgs*0.1,'Pa'  ! g/cm/s^2
    print*,'sound speed       ',cs_envelope_cgs*1E-5,'km/s'
 endif

 proceed = 'y'
 call prompt('Do you wish to continue?',proceed)
 if (trim(adjustl(proceed)) == 'n') then
    stop
 elseif (trim(adjustl(proceed)) /= 'y') then
    stop 'Invalid input'
 else
    open(2022,file='clouds_envlp_info.txt')

    write(2022,*) '-Simulation-'
    write(2022,*) 'tmax              ',tmax*utime/(1E6*years),'Myr / ',tmax*utime,'s'
    write(2022,*) 'dtmax             ',dtmax*utime/(1E6*years),'Myr / ',dtmax*utime,'s'

    if (.not.remove_cloud1) then
       write(2022,*) '-Cloud 1-'
       write(2022,*) 'density           ',rho_cloud1_cgs,'g/cm^3'
       write(2022,*) 'total mass        ',totmass_cloud1,mass_unit
       write(2022,*) 'radius            ',r_cloud1,dist_unit
       write(2022,*) 'volume            ',vol_cloud1*udist**3/pc**3,'pc^3'
       write(2022,*) 'free-fall time    ',t_ff_cloud1,'/ ',t_ff_cloud1*utime/(1E6*years),'Myr / ',t_ff_cloud1*utime,'s'
       write(2022,*) 'temperature       ',temp_cloud1,'K'
       write(2022,*) 'sound speed       ',cs_cloud1_cgs*1E-5,'km/s'
       write(2022,*) 'Mach number       ',rms_mach_cloud1
       write(2022,*) 'v_rms             ',v_rms_kms_cloud1,'km/s'
       write(2022,*) 'Stagnation radius ',rad_stag,dist_unit
    endif

    if (.not.remove_cloud2) then
       write(2022,*) '-Cloud 2-'
       write(2022,*) 'density           ',rho_cloud2_cgs,'g/cm^3'
       write(2022,*) 'total mass        ',totmass_cloud2,mass_unit
       write(2022,*) 'Jeans mass        ',mjeans_cloud2,mass_unit
       write(2022,*) 'radius            ',r_cloud2,dist_unit
       write(2022,*) 'volume            ',vol_cloud2*udist**3/pc**3,'pc^3'
       write(2022,*) 'free-fall time    ',t_ff_cloud2,'/ ',t_ff_cloud2*utime/(1E6*years),'Myr / ',t_ff_cloud2*utime,'s'
       write(2022,*) 'temperature       ',temp_cloud2,'K'
       write(2022,*) 'sound speed       ',cs_cloud2_cgs*1E-5,'km/s'
       write(2022,*) 'Mach number       ',rms_mach_cloud2
       write(2022,*) 'v_rms             ',v_rms_kms_cloud2,'km/s'
    endif

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
 use datafiles,  only:find_phantom_datafile
 real,   intent(in) :: rho_cgs
 integer, parameter :: maxrow = 1000
 integer, parameter :: maxcol = 4
 integer :: irow,icol,iclosest,io
 real    :: rhoteq(maxcol,maxrow)

 open(3000,file=find_phantom_datafile('rhoteq_table.txt','coolcurve_rho_temp'))
 rewind(3000)
 do irow = 1,maxrow
    read(3000,*,iostat=io) (rhoteq(icol,irow), icol=1,maxcol)
    if (io /= 0) stop 'error reading file'
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
 call write_inopt(gmw_in,'gmw_in','mean molecular weight',iunit)
 call write_inopt(make_sinks,'make_sinks','dynamically create sink particles',iunit)
 call write_inopt(create_channel,'create_channel','manually open a channel',iunit)

 write(iunit,"(/,a)") '# cloud1 settings'
 call write_inopt(np_cloud1,'np_cloud1','requested number of particles in cloud1',iunit)
 call write_inopt(rho_cloud1_cgs,'rho_cloud1_cgs','density of cloud1 in g/cm^3',iunit)
 call write_inopt(rms_mach_cloud1,'rms_mach_cloud1','turbulent rms mach number of cloud1',iunit)
 call write_inopt(omega_channel,'omega_channel','solid angle of the channel to carve',iunit)
 call write_inopt(pfrac_channel,'pfrac_channel','fraction of particles to remove in the channel',iunit)

 write(iunit,"(/,a)") '# cloud2 settings'
 call write_inopt(np_cloud2,'np_cloud2','requested number of particles in cloud2',iunit)
 call write_inopt(rho_cloud2_cgs,'rho_cloud2_cgs','density of cloud2 in g/cm^3',iunit)
 call write_inopt(rms_mach_cloud2,'rms_mach_cloud2','turbulent rms mach number of cloud2',iunit)
 call write_inopt(cloud_sep_pc,'cloud_sep_pc','separation from cloud1 in pc',iunit)

 write(iunit,"(/,a)") '# envelope settings'
 call write_inopt(np_envelope,'np_envelope','requested number of particles within envelope boundaries',iunit)
 call write_inopt(rho_envelope_cgs,'rho_envelope_cgs','density of envelope in g/cm^3',iunit)
 call write_inopt(r1_envelope,'r1_envelope','ratio of semi-axis a of ellipsoidal envelope',iunit)
 call write_inopt(r2_envelope,'r2_envelope','ratio of semi-axis b of ellipsoidal envelope',iunit)
 call write_inopt(r3_envelope,'r3_envelope','ratio of semi-axis c of ellipsoidal envelope',iunit)

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
 call read_inopt(np_cloud1,'np_cloud1',db,ierr)
 call read_inopt(np_cloud2,'np_cloud2',db,ierr)
 call read_inopt(np_envelope,'np_envelope',db,ierr)
 call read_inopt(rho_cloud1_cgs,'rho_cloud1_cgs',db,ierr)
 call read_inopt(rho_cloud2_cgs,'rho_cloud2_cgs',db,ierr)
 call read_inopt(rms_mach_cloud1,'rms_mach_cloud1',db,ierr)
 call read_inopt(rms_mach_cloud2,'rms_mach_cloud2',db,ierr)
 call read_inopt(cloud_sep_pc,'cloud_sep_pc',db,ierr)
 call read_inopt(rho_envelope_cgs,'rho_envelope_cgs',db,ierr)
 call read_inopt(r1_envelope,'r1_envelope',db,ierr)
 call read_inopt(r2_envelope,'r2_envelope',db,ierr)
 call read_inopt(r3_envelope,'r3_envelope',db,ierr)
 call read_inopt(omega_channel,'omega_channel',db,ierr)
 call read_inopt(pfrac_channel,'pfrac_channel',db,ierr)
 call read_inopt(gmw_in,'gmw_in',db,ierr)
 call read_inopt(make_sinks,'make_sinks',db,ierr)
 call read_inopt(create_channel,'create_channel',db,ierr)

 call close_db(db)
 !
 ! parse units
 !
 call select_unit(mass_unit,umass,nerr)
 if (nerr /= 0) then
    call error('setup_twosphere_channel','mass unit not recognised')
    ierr = ierr + 1
 endif
 call select_unit(dist_unit,udist,nerr)
 if (nerr /= 0) then
    call error('setup_twosphere_channel','length unit not recognised')
    ierr = ierr + 1
 endif
 if (ierr > 0) then
    print "(1x,a,i2,a)",'setup_twosphere_channel: ',nerr,' error(s) during read of setup file.  Re-writing.'
 endif

end subroutine read_setupfile

end module setup
