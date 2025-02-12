!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module photoionize_cmi
!
! The CMI suite: *photoionize_cmi.F90* kdtree_cmi.f90 hnode_cmi.f90 heating_cooling_cmi.f90 
!                utils_cmi.f90
! This module contains all the subroutines necessary for doing photoionization
! using the Monte Carlo Radiative Transfer code CMacIonize.
!
! Brief instructions on how to compile Phantom with CMacIonize can be found on 
! https://github.com/Cheryl-Lau/phantom/blob/master/phantom_cmi_compilation_instructions.txt
! Note that the compilation method may vary on different machines, especially the ld flags. 
!
!
! Two options are available -
!  1. Pass all particles to CMI for density-mapping and construct grid cell around each
!     individual particle; or
!  2. Extract higher-level tree nodes relative to the given source location(s) and pass
!     them as pseudo-particles to CMI for both density-mapping and grid-construction.
!     Individual particles can be retained in regions which are close to the source.
!     This code is capable of automatically adjusting the tree-walk to ensure that the
!     ionization front is resolved. Activate auto_opennode and tweak parameters including
!     nHlimit_fac, delta_rcut_cgs and tree_accuracy_cmi for best results.
!
! Photoionization heating can be done in the following modes -
!  1. sinks as sources        / user-defined sources
!  2. monochromatic source    / blackbody source         (blackbody only if using sinks)
!  3. instantly heat to 1E4 K / compute heating and cooling
!  4. update u implicitly     / update u explicitly      (implicit only if doing instant-heat)
!
! :References: Petkova,et.al,2021,MNRAS,507,858
!              Bisbas,et.al,2015,MNRAS,453,1324
!              Diaz-Miller,et.al,1998,ApJ,501,192
!
! :Owner: Cheryl Lau (adapted from Maya Petkova's analysis mod)
!
! :Runtime parameters:
!   - sink_ionsrc         : *Using sinks as ionizing sources*
!   - masscrit_ionize_cgs : *Minimum mass of sink that emits ionizing radiation*
!   - niter_mcrt          : *Number of iterations to release photon packets in MCRT simulation*
!   - nphoton             : *Number of photon packets per iteration*
!   - photon_eV           : *Energy of ionizing photons in eV*
!   - monochrom_source    : *Use source with monochromatic frequency which corresponds to energy photon_eV,
!                            else follows a planck distribution at temperature temp_star*
!   - tol_vsite           : *Tolerence in nodes' fractional position change above which the Voronoi
!                            generating-sites will be updated*
!   - lloyd               : *Do Lloyd iterations for Voronoi grid construction*
!   - photoionize_tree    : *Pass tree nodes as pseudo-particles to CMI, else pass all particles*
!   - tree_accuracy_cmi   : *threshold angle to open tree node during tree-walk*
!   - rcut_opennode_cgs   : *Radius within which we must get leaves*
!   - rcut_leafpart_cgs   : *Radius within which we extract individual particles*
!   - auto_opennode       : *Automatically adjust tree-walk according to CMI neutral frac output*
!   - auto_tree_acc       : *Automatically adjust tree_accuracy_cmi to avoid jumps in sparseness of nodes*
!   - delta_rcut_cgs      : *Increase in rcut_opennode_cgs and rcut_leafpart_cgs in an iteration step*
!   - nHlimit_fac         : *Free parameter which controls the resolution of the ionization front*
!   - min_nodesize_toflag : *Minimum node size (relative to root node) to check neutral frac*
!   - crop_domain         : *Automatically crop the simulation domain being exposed to CMI*
!   - crop_fac            : *Cropping out a spherical region of radius crop_fac*rcut_opennode*
!   - temp_hii            : *Presumed temperature of ionized gas*
!   - fix_temp_hii        : *Heats ionized particles to temp_hii, else computes heating and cooling*
!   - implicit_cmi        : *Update internal energy of particles with implicit method, else explicit*
!   - treat_Rtype_phase   : *Check that dt is greater than the duration of R-type phase of HII region*
!                            P.S. if anyone knows how to calculate heating rate under non-ionization
!                            equilibrium without using MC estimator pls let me know....
!
! :Dependencies: infile_utils, physcon, units, io, dim, boundaries, eos, part, kdtree, linklist,
!                kdtree_cmi, hnode_cmi, heatcool_cmi
!
! :Warning: CMacIonize easily runs into seg-fault during grid initialization. Bringing up hlimit_fac 
!           and extradist_fac should resolve the issue, but make sure it is not merging too many cells.
!
 use cmi_fortran_library

 implicit none

 public :: init_ionizing_radiation_cmi,set_ionizing_source_cmi,compute_ionization_cmi
 public :: energy_checks_cmi,energ_implicit_cmi,energ_explicit_cmi
 public :: read_options_photoionize,write_options_photoionize

 logical, public :: inject_rad = .true.  ! switch on/off radiation for testing

 ! Position of sources emitting radiation at current time
 integer, public, parameter :: maxphotosrc = 10
 integer, public :: nphotosrc                      !- Current number of sources
 real   , public :: xyz_photosrc(3,maxphotosrc)    !- Locations of current sources [code units]

 ! Using sinks as source
 logical, public :: sink_ionsrc = .false.
 logical, public :: one_sink_ionsrc = .true.       !- Set a specific sink to be the source
 integer, public :: isink_ionsrc = 5               !- Index of this sink
 real,    public :: masscrit_ionize_cgs = 1.989E34 ! 10 M_sun
 logical, public :: sink_as_cluster = .false.  
 !- or
 ! Manually set location, starting/ending time and ionizing photon flux [cgs units] of sources
 integer, public, parameter :: nsetphotosrc = 1
 real,    public :: xyztq_setphotosrc_cgs(6,nsetphotosrc) = reshape((/ 0.,0.,0. , 0.,1E60,1E49 /),&
                                                                    shape=(/6,nsetphotosrc/))

 ! Monte Carlo simulation settings
 integer, public :: nphoton    = 1E6
 integer, public :: niter_mcrt = 10
 real,    public :: photon_eV  = 13.6   ! used only if sink_ionsrc=F and monochrom_source=T
 real,    public :: tol_vsite  = 1E-4
 logical, public :: lloyd      = .true.
 logical, public :: monochrom_source = .false.   ! else blackbody spec
 logical, public :: force_large_Q = .false. 

 ! Move particles up the tree
 logical, public :: photoionize_tree = .true.

 ! Options for extracting cmi-nodes from kdtree
 real,    public :: tree_accuracy_cmi = 0.3
 real,    public :: rcut_opennode_cgs = 1.2E18   ! 0.4 pc
 real,    public :: rcut_leafpart_cgs = 9.3E17   ! 0.3 pc
 real,    public :: delta_rcut_cgs    = 3.1E16   ! 0.01 pc
 real,    public :: nHlimit_fac       = 100      ! ionization front resolution; recommend 60-100
 real,    public :: min_nodesize_toflag = 1E-5   ! min node size as a fraction of root node
 logical, public :: auto_opennode = .true.
 logical, public :: auto_tree_acc = .false.

 ! Options for cropping simulation domain being passed to CMI (to avoid segfault; use with caution)
 real,    public :: crop_fac     = 1.2          ! bounds = crop_fac*rcut_opennode (>1)
 logical, public :: crop_domain  = .false.

 ! Options for modifying the voronoi grid to merge the smallest cells (to avoid segfault; use with caution)
 logical, public :: limit_voronoi = .true.
 real,    public :: hlimit_fac    = 1E-3    ! threshold in h to merge 
 real,    public :: extradist_fac = 2.0     ! merging distance factor (>= 1.)

 ! Options for heating/cooling
 real,    public :: temp_hii     = 1E4          ! K
 logical, public :: fix_temp_hii = .false.      ! else computes heating and cooling
 logical, public :: implicit_cmi = .true.       ! else updates u explicitly
 logical, public :: treat_Rtype_phase = .false.

 ! Global storages required for updating u
 real,    public,   allocatable :: vxyzu_beforepred(:,:)  ! u before predictor step
 real,    public,   allocatable :: du_cmi(:)              ! to heat implicitly in step_leapfrog
 real,    public,   allocatable :: dudt_cmi(:)            ! to heat explicitly in force

 private

 ! Arrays to store properties of nodes
 integer, parameter   :: maxcminode   = 1E8
 integer, parameter   :: maxleafparts = 1E8
 real,    allocatable :: nxyzm_treetocmi(:,:),ixyzhm_leafparts(:,:)
 integer, allocatable :: nnode_toreplace(:)
 real,    allocatable :: h_solvertocmi(:)

 ! Array to store properties of all sites
 real,    allocatable :: nixyzhmf_cminode(:,:)

 ! Arrays to control subsequent tree-walks
 integer, parameter   :: maxnode_open = 1E6
 integer, allocatable :: nnode_needopen(:)    ! for next iteration
 real,    allocatable :: nxyzrs_nextopen(:,:) ! for next timestep
 real,    allocatable :: nxyzrs_nextopen_updatewalk(:,:)
 integer :: nnextopen,nnextopen_updatewalk
 integer :: ncminode_previter
 real    :: tree_acc_tol = 0.1

 ! Arrays to store CMI outputs for computing du_cmi/dudt_cmi
 real,    allocatable :: nH_allparts(:)

 ! Arrays to store Voronoi grid
 real,    allocatable :: x_old(:),y_old(:),z_old(:)
 integer :: nsite_lastgrid

 integer, parameter :: maxoutfile_ult = 99999
 integer :: maxoutfile = 5000         ! max number of (ni)xyzhmnH output files
 integer :: ncall_writefile = 100     ! interval to write (ni)xyzhmnH output file
 integer :: icall,iunit,ifile,iruncmi
 integer :: ncall_checktreewalk = 50  ! interval to check for unnecessarily-opened nodes
 integer :: nphotosrc_old
 real    :: xyz_photosrc_si(3,maxphotosrc),ionflux_photosrc(maxphotosrc),mass_photosrc(maxphotosrc),masscrit_ionize
 real    :: tree_accuracy_cmi_old
 real    :: rcut_opennode,rcut_leafpart,delta_rcut
 real    :: cen_crop(3),r_crop2,rcrop_min
 real    :: temp_star,totq,freq_photon,u_hii,udist_si,umass_si
 real    :: time0_wall,time0_cpu,time_now_wall,time_now_cpu
 real    :: time_ellapsed_wall,time_ellapsed_cpu
 logical :: is_Rtype_phase,old_sources_exist
 logical :: first_call,first_step,warned

 ! Switches for plotting/debugging
 logical :: write_gamma = .false.          ! write heating rates vs nH (from both phantom and CMI)
 logical :: print_cmi   = .false.          ! show CMI shell outputs
 logical :: write_nH_u_distri  = .false.   ! write u of particles vs nH
 logical :: write_node_prop    = .true.    ! write properties of the current set of cmi-nodes
 logical :: plot_cropped_sites = .false.   ! write properties of the current cropped set of cmi-nodes
 logical :: catch_noroot_parts = .false.   ! write particles with no therm-equil roots

contains

!-----------------------------------------------------------------------------
!+
! Initializations
!+
!-----------------------------------------------------------------------------
subroutine init_ionizing_radiation_cmi(time,npart,xyzh,nptmass,dt)
 use physcon,  only:mass_proton_cgs,kboltz,c,planckh,solarl,eV
 use eos,      only:gmw,gamma
 use io,       only:warning,fatal
 use units,    only:udist,umass,utime,unit_ergg
 use dim,      only:maxvxyzu
 use heatcool_cmi, only:init_ueq_table,precompute_uterms
 use omp_lib
 integer, intent(in) :: npart,nptmass
 real,    intent(in) :: time,dt
 real,    intent(in) :: xyzh(:,:)
 integer :: ip,io_file,isrc
 real    :: hmean,psep,wavelength_cgs,energ_photon_cgs,time_src0,time_src
 real    :: gmw0,csi_cgs,temp_hii_fromcsi,csi_cgs_req
 logical :: compilecond_ok

 print*,'Radiation-hydrodynamics: Phantom is coupled to photoionization code CMacIonize'
 print*,'Injecting ionizing radiation with MCRT method'

 !
 ! User input checking
 !
 if (.not.photoionize_tree .and. maxcminode < npart) call fatal('photoionize_cmi',&
                                                   & 'maxcminode has to be greater than npart')
 if (maxvxyzu < 4) call fatal('photoionize_cmi','Not suitable for isothermal simulations')

 !- Check that the tol_vsite is sensible by estimating mean particle separation
 hmean = 0.
 do ip = 1,npart
    hmean = hmean + xyzh(4,ip)
 enddo
 hmean = hmean/npart
 psep = 2.*hmean / 57.9**(1./3.)
 if (tol_vsite > 0.5*psep) call warning('photoionize_cmi','tol_vsite might be too large')

 !- Check pick-nodes settings
 if (photoionize_tree) then
    if (tree_accuracy_cmi == 0) call warning('photoionize_cmi','extracting leaf nodes only')
    rcut_opennode = rcut_opennode_cgs/udist
    rcut_leafpart = rcut_leafpart_cgs/udist
    delta_rcut    = delta_rcut_cgs/udist
    if (rcut_opennode < rcut_leafpart) call fatal('photoionize_cmi','rcut_leafpart must be &
                                                 & smaller than rcut_opennode')
    !- Check compile conditions
    compilecond_ok = .false.
#ifndef PERIODIC
#ifndef MPI
    compilecond_ok = .true.
#endif
#endif
    if (.not.compilecond_ok) call fatal('photoionize_cmi','current version does not support PERIODIC or MPI')
 endif

 !- Check settings for cropping simulation 
 if (crop_domain) then 
    if (crop_fac < 1.0) call fatal('photoionize_cmi','we need larger nodes on the side to detect the ionization front!')
    if (crop_fac < 2.0) call warning('photoionize_cmi','recommend setting a larger crop_fac')
 endif 

 !- Check options for ionizing source
 if (monochrom_source) then
    if (sink_ionsrc) then
       call warning('photoionize_cmi','Use monochromatic spectrum for emissions from sinks?')
    endif 
    write(*,'(2x,a42,f5.2,a3)') 'Source: monochromatic photons with energy ',photon_eV,' eV'
    if (.not.fix_temp_hii) then ! do heating/cooling
       temp_star = (photon_eV-13.6)*eV /(3./2.*kboltz)
       if (temp_star < tiny(temp_star)) call fatal('photoionize_cmi','photon_eV needs to be greater than 13.6 eV')
       write(*,'(3x,a33,es10.4,a2)') 'setting source temperature to be ',temp_star,' K'
    endif
 else ! planck spectrum
    if (sink_ionsrc) then
       write(*,'(2x,a37)') 'Source: sinks with blackbody spectrum'
       write(*,'(3x,a70)') 'temp_star will be computed using mass of sink particles during runtime'
    else ! user-set
       write(*,'(2x,a45)') 'Source: point sources with blackbody spectrum'
       write(*,'(3x,a68)') 'temp_star will be computed using ionizing photon flux during runtime'
    endif
 endif
 ! Note: temp_star is required if using blackbody source and/or doing heating, however,
 !       both can only take one value regardless of the number of ionizing sources.
 !       Here, we will average it over all present sources.

 !
 ! Determine whether or not ionizing sources exist in starting dump
 !- controls the dt > t_recomb check during source injection
 !
 if (treat_Rtype_phase) then
    old_sources_exist = .false.
    if (sink_ionsrc) then
       if (nptmass > 0) old_sources_exist = .true.
    else
       time_src0 = huge(time_src0)
       do isrc = 1,nsetphotosrc
          time_src = xyztq_setphotosrc_cgs(4,isrc)/utime
          time_src0 = min(time_src0,time_src)
       enddo
       if (time > 0. .and. time_src0 < time-dt) old_sources_exist = .true.
    endif
    if (old_sources_exist) call warning('photoionize_cmi','existing ionizing sources in starting dumpfile')
 endif

 !
 ! Convert cgs units to SI units for CMI
 !
 udist_si = udist/1E2
 umass_si = umass/1E3

 !
 ! Memory allocations and initialization for adaptive tree-walk system
 !
 call allocate_cminode
 call reset_cminode
 call allocate_cmi_inputs_history
 call reset_cmi_inputs_history
 nnextopen = 0
 nnextopen_updatewalk = 0
 nxyzrs_nextopen = 0.
 nxyzrs_nextopen_updatewalk = 0.
 tree_accuracy_cmi_old = tree_accuracy_cmi

 !
 ! Internal energy of ionized particles
 !
 if (fix_temp_hii) then
    print*,'Photoionization heating/cooling disabled - Ionized particles will be heated to temp_hii'
    if (.not.implicit_cmi) then
       print*,' Implicit method is required if using a fixed u_hii. Setting implicit_cmi flag to T.'
       implicit_cmi = .true.
    endif
    gmw0 = gmw
    gmw  = 0.5  !- temporarily change the mean molecular weight of ionized particles
    u_hii = kboltz * temp_hii / (gmw*mass_proton_cgs*(gamma-1.)) /unit_ergg
    csi_cgs = sqrt(temp_hii*(gamma*kboltz)/(gmw*mass_proton_cgs))
    csi_cgs_req = 12.85E5  !- ref: Bisbas15
    temp_hii_fromcsi = (csi_cgs_req)**2*gmw*mass_proton_cgs/(gamma*kboltz)
    print*,' -Ionized gas properties- '
    print*,' internal energy: ',u_hii*unit_ergg,'erg/g'
    print*,' sound speed:     ',csi_cgs,'cm/s'
    print*,' temperature:     ',temp_hii_fromcsi,'K'
    gmw = gmw0
 else
    print*,'Photoionization heating/cooling activated'
 endif

 !- Calculate ionizing photon frequency [Hz]
 energ_photon_cgs = photon_eV*eV
 wavelength_cgs   = planckh*c/energ_photon_cgs
 freq_photon = c/wavelength_cgs

 !- Critical mass of sinks
 if (sink_ionsrc) masscrit_ionize = masscrit_ionize_cgs/umass

 !- Init for photoionization heating-cooling routines
 call reset_energ
 call precompute_uterms
 if (implicit_cmi .and. .not.fix_temp_hii) call init_ueq_table

 !- For indicating R-type phase for heating
 nphotosrc_old = 0
 first_step = .true.

 !- For writing snapshots of nH map
 if (maxoutfile > maxoutfile_ult) call fatal('photoionize_cmi','maxoutfile must be within 5 digits')
 call init_write_snapshot

 !- For writing Voronoi sites files / checking openned nodes
 icall = 0
 first_call = .true.

 !- Time the simulations
 call cpu_time(time0_cpu)
 time0_wall = omp_get_wtime()
 open(2050,file='cpu_wall_time_record.txt',status='replace',iostat=io_file)
 if (io_file /= 0) call fatal('photoionize_cmi','unable to open time-record file')

 open(2060,file='CMI_cpu_wall_time_record.txt',status='replace',iostat=io_file)
 if (io_file /= 0) call fatal('photoionize_cmi','unable to open CMI time-record file')

 open(2070,file='cpu_wall_time_indv_process.txt',status='replace',iostat=io_file)
 if (io_file /= 0) call fatal('photoionize_cmi','unable to open time-record file for ind. processes')

end subroutine init_ionizing_radiation_cmi


!-----------------------------------------------------------------------------
!+
! Routines to prepare for CMI-call and energy-updates at the beginning of each step
!+
!-----------------------------------------------------------------------------
!
! Dynamically set the locations of the sources at current timestep
! Updates nphotosrc and xyz_photosrc
!
subroutine set_ionizing_source_cmi(time,nptmass,xyzmh_ptmass)
 use io,        only:fatal
 use physcon,   only:solarm
 use units,     only:utime,udist
 use utils_cmi, only:mag2
 integer, intent(in) :: nptmass
 real,    intent(in) :: time
 real,    intent(in) :: xyzmh_ptmass(:,:)
 integer :: isrc,ix,isink
 real    :: time_startsrc,time_endsrc,mptmass,fluxq
 real    :: xmin_allsrc,xmax_allsrc,ymin_allsrc,ymax_allsrc,zmin_allsrc,zmax_allsrc
 real    :: disttomin2,disttomax2
 logical :: target_sink_found

 !- Init/Reset
 nphotosrc = 0
 xyz_photosrc   = 0.
 ionflux_photosrc = 0.
 mass_photosrc  = 0.

 !
 ! Extract sources which currently emit radiation
 !
 if (sink_ionsrc .and. nptmass > 0) then !- Check sink
    if (one_sink_ionsrc) then 
       nphotosrc = 1
       target_sink_found = .false. 
       over_sinks: do isink = 1,nptmass 
          if (isink == isink_ionsrc) then 
             mptmass = xyzmh_ptmass(4,isink)
             if (mptmass < 0.) call fatal('photoionize_cmi','user-specified sink has been accreted')
             target_sink_found = .true. 
             xyz_photosrc(1:3,nphotosrc) = xyzmh_ptmass(1:3,isink)
             fluxq = get_ionflux_star(mptmass)
             ionflux_photosrc(nphotosrc) = fluxq    ! [nphoton/s]
             mass_photosrc(nphotosrc)    = mptmass
             exit over_sinks
          endif 
       enddo over_sinks 
       if (.not.target_sink_found) call fatal('photoionize_cmi','target sink not found')
    else
       do isink = 1,nptmass
          mptmass = xyzmh_ptmass(4,isink)
          if (mptmass >= masscrit_ionize) then
             nphotosrc = nphotosrc + 1
             if (nphotosrc > maxphotosrc) call fatal('photoionize_cmi','number of sources &
                                                    &exceeded maxphotosrc')
             xyz_photosrc(1:3,nphotosrc) = xyzmh_ptmass(1:3,isink)
             fluxq = get_ionflux_star(mptmass)
             ionflux_photosrc(nphotosrc) = fluxq    ! [nphoton/s]
             mass_photosrc(nphotosrc)    = mptmass
          endif
       enddo
    endif 
 elseif (.not.sink_ionsrc) then !- Check time
    do isrc = 1,nsetphotosrc
       time_startsrc = xyztq_setphotosrc_cgs(4,isrc)/utime
       time_endsrc   = xyztq_setphotosrc_cgs(5,isrc)/utime
       if (time > time_startsrc .and. time < time_endsrc) then
          nphotosrc = nphotosrc + 1
          if (nphotosrc > maxphotosrc) call fatal('photoionize_cmi','number of sources &
                                                  &exceeded maxphotosrc')
          xyz_photosrc(1:3,nphotosrc) = xyztq_setphotosrc_cgs(1:3,isrc)/udist
          ionflux_photosrc(nphotosrc) = xyztq_setphotosrc_cgs(6,isrc)  ! [nphoton/s]
       endif
    enddo
 endif
 !
 ! Locate the boundary that encapsulates all sources (to be used for cropping)
 ! 
 if (crop_domain) then 
    xmax_allsrc = -huge(xmax_allsrc)
    ymax_allsrc = -huge(ymax_allsrc)
    zmax_allsrc = -huge(zmax_allsrc)
    xmin_allsrc = huge(xmin_allsrc)
    ymin_allsrc = huge(ymin_allsrc)
    zmin_allsrc = huge(zmin_allsrc)
    do isrc = 1,nphotosrc
       xmin_allsrc = min(xmin_allsrc,xyz_photosrc(1,isrc))
       ymin_allsrc = min(ymin_allsrc,xyz_photosrc(2,isrc))
       zmin_allsrc = min(zmin_allsrc,xyz_photosrc(3,isrc))
       xmax_allsrc = max(xmax_allsrc,xyz_photosrc(1,isrc))
       ymax_allsrc = max(ymax_allsrc,xyz_photosrc(2,isrc))
       zmax_allsrc = max(zmax_allsrc,xyz_photosrc(3,isrc))
    enddo 
    !- centre of all sources 
    cen_crop = (/ (xmin_allsrc + xmax_allsrc) /2., &
                  (ymin_allsrc + ymax_allsrc) /2., &
                  (zmin_allsrc + zmax_allsrc) /2.  /)
    !- estimate minimum range to crop out 
    disttomin2 = mag2((/xmin_allsrc,ymin_allsrc,zmin_allsrc/) - cen_crop)
    disttomax2 = mag2((/xmax_allsrc,ymax_allsrc,zmax_allsrc/) - cen_crop)
    rcrop_min  = max(disttomin2,disttomax2)
 endif 

 !- Convert to SI units for CMI param file
 if (nphotosrc >= 1) then
    do isrc = 1,maxphotosrc
       do ix = 1,3
          xyz_photosrc_si(ix,isrc) = xyz_photosrc(ix,isrc)*udist_si
       enddo
    enddo
 endif

 if (.not.inject_rad) nphotosrc = 0

end subroutine set_ionizing_source_cmi

!
! Prepare for energy injection
!
subroutine energy_checks_cmi(xyzh,dt)
 use physcon,      only:solarr,steboltz,pi
 use heatcool_cmi, only:check_to_stop_cooling,compute_Rtype_time
 use io,           only:fatal
 real, intent(in) :: xyzh(:,:)
 real, intent(in) :: dt
 integer :: isrc
 real    :: massmean,lumin_star_cgs,rstar_cgs,t_recomb,qmean
 real    :: rstar_solarr = 10.  !- radius of a typical O-star [R_sun]
 logical :: new_source_injected

 if (.not.fix_temp_hii) call check_to_stop_cooling(nphotosrc)

 if (nphotosrc >= 1) then
    !- Ready to rerun CMI
    call reset_energ

    if (sink_ionsrc .and. .not.monochrom_source) then
       !- Estimate temperature of star with mean mass of current sources
       massmean = 0.
       do isrc = 1,nphotosrc
          massmean = massmean + mass_photosrc(isrc)
       enddo
       massmean = massmean/nphotosrc
       lumin_star_cgs = get_lumin_star(massmean)
       rstar_cgs = rstar_solarr*solarr
       temp_star = (lumin_star_cgs/(4.*pi*rstar_cgs**2*steboltz))**(1./4.)

    elseif (.not.sink_ionsrc .and. .not.monochrom_source) then
       !- Estimate temperature of star with mean ionizing photon flux of current sources
       qmean = 0.
       do isrc = 1,nphotosrc
          qmean = qmean + ionflux_photosrc(isrc)
       enddo
       qmean = qmean/nphotosrc
       temp_star = get_temp_star(qmean)
    endif
 endif

 if (inject_rad .and. treat_Rtype_phase) then
    new_source_injected = nphotosrc > nphotosrc_old
    if (first_step .and. old_sources_exist) new_source_injected = .false.
    if (new_source_injected) then
       is_Rtype_phase = .true.
       !- Check that dt is greater than recombination timescale
       call compute_Rtype_time(nphotosrc,xyz_photosrc,xyzh,t_recomb)
       if (t_recomb > dt) call fatal('photoionize_cmi','dt < t_recomb - require larger timestep')
    else
       if (is_Rtype_phase) is_Rtype_phase = .false.  !- switch back off
    endif
    nphotosrc_old = nphotosrc
 endif
 !- Note: CMI does not model the initial R-type phase during the expansion of HII region,
 !        hence this phase should be an instantaneous process wrt the hydrodynamics -
 !        dt should cover the time duration of R-type phase

end subroutine energy_checks_cmi


!
! Calculate luminosity [cgs units] from mass of sink particle [code units]
!
real function get_lumin_star(mass_star)
 use physcon, only:solarm,solarl
 use units,   only:umass
 real, intent(in)  :: mass_star
 real :: mass_star_cgs

 mass_star_cgs  = mass_star*umass
 get_lumin_star = solarl * (mass_star_cgs/solarm)**(3.5)

end function get_lumin_star

!
! Get ionizing photon flux Q [cgs units] from mass of sink particle [code units]
! obtained from fitting mass-flux in table 1 of Diaz-Miller et al. 1998
!
real function get_ionflux_star(mass_star)
 use physcon,   only:solarm
 use units,     only:umass
 use utils_cmi, only:imf_constant,cluster_totflux
 real, intent(in) :: mass_star
 real :: mass_star_cgs,mass_star_solarm
 real :: mass_sink_solarm,constant,mass_star_lowerlim,totflux

 if (.not.sink_as_cluster) then 
    mass_star_cgs = mass_star*umass
    mass_star_solarm = mass_star_cgs/solarm 

!    print*,'mass of ionizing star',mass_star_solarm
    get_ionflux_star = 10**(48.1 + 0.02*(mass_star_solarm - 20.d0))

!    print*,'forcefully set flux to 1E51'
!   get_ionflux_star = 1E51 ! testing 
   if (force_large_Q) get_ionflux_star = 1E51

   ! if (mass_star_cgs > 3.65E34) then
   !    get_ionflux_star = 10**(2.817*log10(mass_star_cgs) - 49.561)
   ! else
   !    get_ionflux_star = 10**(12.548*log10(mass_star_cgs) - 385.885)
   ! endif
 else 
    mass_sink_solarm = mass_star*umass/solarm
    constant = imf_constant(mass_sink_solarm)
    mass_star_lowerlim = 30 
    get_ionflux_star = cluster_totflux(constant,mass_star_lowerlim)
 endif 

end function get_ionflux_star

!
! Get temperature from ionizing photon flux Q of source
! obtained from fitting flux-temp in table 1 of Diaz-Miller et al. 1998
!
real function get_temp_star(ionflux_src)
 real, intent(in) :: ionflux_src

 get_temp_star =  0.887 * exp(0.209*log10(ionflux_src) + 0.278) + 7613.222

end function get_temp_star


!-----------------------------------------------------------------------------
!+
! Wrapper for computing nH of all particles at current timestep
!+
!-----------------------------------------------------------------------------
subroutine compute_ionization_cmi(time,npart,xyzh,vxyzu)
 use part,     only:massoftype,igas,isdead_or_accreted
 use kdtree,   only:inodeparts,inoderange
 use io,       only:fatal,warning
 use omp_lib
 integer, intent(inout) :: npart
 real,    intent(in)    :: time
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,    allocatable   :: x(:),y(:),z(:),h(:),m(:),nH(:)
 integer :: ip,ip_cmi,npart_cmi,i,n,ipnode,isite,ncminode,isrc
 real    :: nH_part,nH_site

 if (nphotosrc == 0) return
 print*,'Injecting radiation at'
 do isrc = 1,nphotosrc 
    print*,xyz_photosrc(:,isrc), ionflux_photosrc(isrc)
 enddo 

 if (photoionize_tree) then
    !
    ! Pass tree nodes as pseudo-particles to CMI
    !
    call treewalk_run_cmi_iterate(time,xyzh,ncminode)
    call write_nH_snapshot(time,ncminode)

    !- Assign nH to all particles beneath each node using the final nixyzhmf_cminode
    over_entries: do isite = 1,ncminode
       nH_site = nixyzhmf_cminode(8,isite)
       if (nH_site < 0. .or. nH_site > 1.) call fatal('photoionize_cmi','invalid nH')
       n = int(nixyzhmf_cminode(1,isite))
       i = int(nixyzhmf_cminode(2,isite))
       if (n /= 0 .and. i == 0) then !- is node
          over_parts: do ipnode = inoderange(1,n),inoderange(2,n)
             ip = abs(inodeparts(ipnode))
             nH_allparts(ip) = nH_site
          enddo over_parts
       elseif (n == 0 .and. i /= 0) then !- is particle
          nH_allparts(i) = nH_site
       else
          call fatal('photoionize_cmi','unidentified site')
       endif
    enddo over_entries
    call reset_cminode
 else
    !
    ! Pass all particles to grid-construction and density mapping
    !
    call allocate_cmi_inoutputs(x,y,z,h,m,nH)
    npart_cmi = 0
    do ip = 1,npart
       if (.not.isdead_or_accreted(xyzh(4,ip))) then
          npart_cmi = npart_cmi + 1
          x(npart_cmi) = xyzh(1,ip)
          y(npart_cmi) = xyzh(2,ip)
          z(npart_cmi) = xyzh(3,ip)
          h(npart_cmi) = xyzh(4,ip)
          m(npart_cmi) = massoftype(igas)
       endif
    enddo
    call run_cmacionize(npart_cmi,x,y,z,h,m,nH)

    !- Collect nH of all alive particles
    ip_cmi = 0
    do ip = 1,npart
       if (.not.isdead_or_accreted(xyzh(4,ip))) then
          ip_cmi = ip_cmi + 1
          nH_part = nH(ip_cmi)
          if (nH_part < 0.) call fatal('photoionize_cmi','invalid nH')
          nH_allparts(ip) = nH_part
       endif
    enddo
    if (ip_cmi /= npart_cmi) call fatal('photoionize_cmi','number of particles &
                                       & passed to and from CMI do not match')

    call write_nH_snapshot(time,npart,xyzh_parts=xyzh,x_in=x,y_in=y,z_in=z,h_in=h,m_in=m,nH_in=nH)
    call deallocate_cmi_inoutputs(x,y,z,h,m,nH)
 endif

 iruncmi = iruncmi + 1
 first_step = .false.

 !
 ! Time the simulations
 !
 call cpu_time(time_now_cpu)
 time_now_wall = omp_get_wtime()
 time_ellapsed_wall = time_now_wall - time0_wall
 time_ellapsed_cpu  = time_now_cpu  - time0_cpu
 open(2050,file='cpu_wall_time_record.txt',position='append')
 write(2050,*) iruncmi, time, time_ellapsed_cpu, time_ellapsed_wall
 close(2050)

end subroutine compute_ionization_cmi

!-----------------------------------------------------------------------------
!+
! Routines for updating the internal energy of particles using the final nH_parts(:)
! returned from CMI
!- Notes:* u-update needs to be separated from nH computation routines since CMI-call and
!-         implicit update are only done during 1st call of deriv, whereas explicit update
!-         needs to be done in both calls.
!-       * nH_allparts(i) = -1 denotes dead/non-existing particles, if nH >= -epsilon
!          means particle has been dealt by CMI (but not necessarily ionized).
!-       * As the original cooling is being switched off, particles which are not heated
!          by photoionization would still go through these routines to cool.
!        * In situations where no thermal equilibrium is found, it probably means that the
!          heating rate is greater than cooling rate at all temperatures, i.e. particle
!          will heat forever. Here we'll drift it to Tmax, however the sim is likely
!          problematic and we recommend checking the heating rates.
!+
!-----------------------------------------------------------------------------
!
! Computes du_cmi(:); to be heated in step_leapfrog
! also updates the vxyzu (i.e. vpred from derivs)
!
subroutine energ_implicit_cmi(time,npart,xyzh,vxyzu,dt)
 use heatcool_cmi, only:heating_term,cooling_term,get_ueq,compute_du
 use part,         only:rhoh,massoftype,igas
 use physcon,      only:kboltz,mass_proton_cgs,years
 use eos,          only:gamma,gmw
 use units,        only:unit_ergg,unit_energ,utime
 use io,           only:warning
 integer, intent(in)    :: npart
 real,    intent(in)    :: time,dt
 real,    intent(in)    :: xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:)
 integer, parameter :: nbuc = 100
 integer :: ip,npart_heated,ibuc,inoroot,numroots
 integer :: iunit_unH
 real    :: nH,ui,ueq,pmass,gammaheat,lambda,rhoi,du
 real    :: ueq_mean,u_mean,temp_ueq,temp_u,gmw0,tau,tau_mean
 real    :: time_cgs
 real    :: pos_noroot(3,npart),nH_noroot(npart),temp_noroot(npart)
 real    :: nH_buc(nbuc),gamma_buc(nbuc),nentry_buc(nbuc)
 character(len=50) :: ifile_char,filename_unH

 if (nphotosrc == 0) return

 gmw0 = gmw
 gmw = 0.5

 pmass = massoftype(igas)

 !- Check variation of gamma with nH
 if (write_gamma) then
    do ibuc = 1,nbuc
       nH_buc(ibuc) = 1./nbuc * (ibuc-1)
    enddo
    gamma_buc  = 0.
    nentry_buc = 0
 endif

 inoroot  = 0
 u_mean   = 0.
 ueq_mean = 0.
 tau_mean = 0.
 npart_heated = 0
 !$omp parallel do default(none) shared(npart,nH_allparts,xyzh,vxyzu) &
 !$omp shared(vxyzu_beforepred,pmass,temp_star,dt,du_cmi) &
 !$omp shared(fix_temp_hii,u_hii) &
 !$omp shared(write_gamma,nH_buc,unit_energ,utime) &
 !$omp shared(catch_noroot_parts,pos_noroot,nH_noroot,temp_noroot) &
 !$omp shared(gmw,gamma,unit_ergg,inoroot) &
 !$omp private(ibuc) &
 !$omp private(ip,nH,rhoi,ui,ueq,gammaheat,lambda,du,tau,numroots) &
 !$omp reduction(+:ueq_mean,u_mean,tau_mean,npart_heated) &
 !$omp reduction(+:gamma_buc,nentry_buc) &
 !$omp schedule(runtime)
 do ip = 1,npart
    nH = nH_allparts(ip)
    ui = vxyzu_beforepred(4,ip)  !- take vxyzu from before predictor as current u
    if (nH > -epsilon(nH)) then
       du = 0.
       if (fix_temp_hii) then
          if (nH < 0.5 .and. ui < u_hii) then
             npart_heated = npart_heated + 1
             du = u_hii - ui     !- instantly heat particle to 10^4 K
          endif
       else
          if (nH > -epsilon(nH)) then
             rhoi = rhoh(xyzh(4,ip),pmass)
             call heating_term(nH,rhoi,ui,temp_star,gammaheat)
             call cooling_term(ui,lambda)  ! for calculating timescale
             call get_ueq(rhoi,gammaheat,ui,numroots,ueq)
             call compute_du(dt,rhoi,ui,ueq,gammaheat,lambda,tau,du)

             !- check current T and Teq of HII region
             if (nH < 0.8) then
                npart_heated = npart_heated + 1
                u_mean   = u_mean + ui
                ueq_mean = ueq_mean + ueq
                tau_mean = tau_mean + tau
             endif

             !- check gamma distribution
             if (write_gamma) then
                ibuc = minloc(abs(nH_buc(:)-nH),1)
                gamma_buc(ibuc)  = gamma_buc(ibuc) + gammaheat*(unit_energ/utime)
                nentry_buc(ibuc) = nentry_buc(ibuc) + 1
             endif

             !- store properties of particles with no roots for trouble-shooting
             if (catch_noroot_parts) then
                if (numroots == 0) then
                   !$omp critical
                   inoroot = inoroot + 1
                   pos_noroot(1:3,inoroot) = xyzh(1:3,ip)
                   nH_noroot(inoroot)   = nH
                   temp_noroot(inoroot) = ui/kboltz*(gmw*mass_proton_cgs*(gamma-1.))*unit_ergg
                   !$omp end critical
                endif
             endif
          endif
       endif
       !- update vpred
       vxyzu(4,ip) = vxyzu(4,ip) + du
       !- store du into global array
       du_cmi(ip) = du
    endif
 enddo
 !$omp end parallel do

 if (npart_heated == 0) then
    print*,' no ionized particles'
    call warning('photoionize_cmi','no ionized particles')
 else 
    print*,'Number of particles heated:   ',npart_heated
    if (.not.fix_temp_hii) then
       ueq_mean = ueq_mean/npart_heated
       temp_ueq = ueq_mean/kboltz*(gmw*mass_proton_cgs*(gamma-1.))*unit_ergg
       print*,'Drifting HII region to temp [K]:   ',temp_ueq
       u_mean = u_mean/npart_heated
       temp_u = u_mean/kboltz*(gmw*mass_proton_cgs*(gamma-1.))*unit_ergg
       print*,'Current temp of HII region  [K]:   ',temp_u
       tau_mean = tau_mean/npart_heated
       print*,'Time remaining to reach Teq [Myr]: ',tau_mean*utime/(1E6*years)
    endif 
 endif 

 if (write_gamma) then
    do ibuc = 1,nbuc
       if (nentry_buc(ibuc) > 0) then
          gamma_buc(ibuc) = gamma_buc(ibuc)/nentry_buc(ibuc)
       endif
    enddo
    open(3012,file='gamma_nH_distri.txt',status='replace')
    do ibuc = 1,nbuc
       write(3012,*) nH_buc(ibuc), gamma_buc(ibuc)
    enddo
    close(3012)
 endif

 if (catch_noroot_parts .and. inoroot > 1) then
    open(5014, file='noroots.txt',status='replace')
    write(5014,'(5a25)') 'x','y','z','nH','temp'
    do ip = 1,inoroot
       write(5014,*) pos_noroot(1:3,ip), nH_noroot(ip), temp_noroot(ip)
    enddo
    close(5014)
 endif

  if (write_nH_u_distri) then
     if (mod(iruncmi,ncall_writefile) == 0) then
        write(ifile_char,'(i5.5)') ifile-1
        iunit_unH = iunit+5000
        filename_unH = 'nH_u_parts_'//trim(adjustl(ifile_char))//'.txt'
        open(iunit_unH,file=filename_unH,status='replace')
        time_cgs = time*utime
        write(iunit_unH,*) time_cgs
        do ip = 1,npart
           write(iunit_unH,*) nH_allparts(ip), vxyzu(4,ip)
        enddo
        close(iunit_unH)
     endif
  endif

 gmw = gmw0

end subroutine energ_implicit_cmi

!
! Computes dudt_cmi(:); to be heated in force
!
subroutine energ_explicit_cmi(npart,xyzh,vxyzu,dt)
 use heatcool_cmi, only:heating_term,cooling_term,compute_dudt
 use part,         only:rhoh,massoftype,igas
 use physcon,      only:kboltz,mass_proton_cgs
 use eos,          only:gamma,gmw
 use units,        only:unit_ergg
 use io,           only:warning
 integer, intent(in) :: npart
 real,    intent(in) :: dt
 real,    intent(in) :: xyzh(:,:)
 real,    intent(in) :: vxyzu(:,:)
 integer :: ip,npart_heated
 real    :: nH,ui,pmass,gammaheat,lambda,rhoi,dudt
 real    :: u_ionized,uhii_mean,temp_ionized,gmw0

 if (nphotosrc == 0) return

 pmass = massoftype(igas)
 uhii_mean = 0.
 npart_heated = 0

 gmw0 = gmw
 gmw = 0.5

 !$omp parallel do default(none) shared(npart,nH_allparts,xyzh,vxyzu) &
 !$omp shared(temp_star,pmass,dt,dudt_cmi) &
 !$omp private(ip,nH,rhoi,ui,gammaheat,lambda,dudt) &
 !$omp private(u_ionized) &
 !$omp reduction(+:npart_heated,uhii_mean) &
 !$omp schedule(runtime)
 do ip = 1,npart
    nH = nH_allparts(ip)
    if (nH > -epsilon(nH)) then
       rhoi = rhoh(xyzh(4,ip),pmass)
       ui   = vxyzu(4,ip)
       call heating_term(nH,rhoi,ui,temp_star,gammaheat)
       call cooling_term(ui,lambda)
       call compute_dudt(dt,rhoi,ui,gammaheat,lambda,dudt)
       !- store dudt into global array
       dudt_cmi(ip) = dudt
       !- checking
       if (nH < 0.8) then
          u_ionized = vxyzu(4,ip)
          npart_heated = npart_heated + 1
          uhii_mean    = uhii_mean + u_ionized
       endif
    endif
 enddo
 !$omp end parallel do

 if (npart_heated == 0) then
    call warning('photoionize_cmi','No ionized particles')
 else
    print*,'Number of particles heated:   ',npart_heated
    uhii_mean = uhii_mean/npart_heated
    temp_ionized = uhii_mean/kboltz*(gmw*mass_proton_cgs*(gamma-1.))*unit_ergg
    print*,'Current temp of HII region [K]: ',temp_ionized
 endif

 gmw = gmw0

end subroutine energ_explicit_cmi

!-----------------------------------------------------------------------------
!+
! Iterate tree-walk and CMI-call until the ionization front is resolved
! Gives the final set of nixyzhmf_cminode
!+
!-----------------------------------------------------------------------------
subroutine treewalk_run_cmi_iterate(time,xyzh,ncminode)
 use linklist,   only:node,ifirstincell
 use kdtree_cmi, only:extract_cminodes_from_tree
 use hnode_cmi,  only:hnode_iterate
 use utils_cmi,  only:mag2 
 use io,         only:fatal,warning
 real,    intent(in)  :: time
 real,    intent(in)  :: xyzh(:,:)
 integer, intent(out) :: ncminode   !- total number of CMI sites (nodes+particles)
 real,    allocatable :: x(:),y(:),z(:),h(:),m(:),nH(:)
 integer :: maxiter = 100
 integer :: nneedopen    !- number of nodes at ionization front that needs to be opened at current iteration
 integer :: ncloseleaf   !- number of leaves to be replaced
 integer :: nleafparts   !- total number of particles in leaf nodes to be replaced
 integer :: niter,inode,isite,n,inextopen,n_nextopen
 integer :: ncminode_in,ncminode_out  !- local var that controls the actual number of sites passed to CMI
 real    :: size_node,size_root,nH_node,nH_part,nH_limit,mintheta
 logical :: node_checks_passed,must_iter

 !- Init
 node_checks_passed = .false.
 niter = 0
 nneedopen = 0
 do inode = 1,maxnode_open
    nnode_needopen(inode) = 0
 enddo

 resolve_ionfront: do while (.not.node_checks_passed)
    write(*,'(1x,a35,i2)') 'Resolve ionization front iteration ',niter
    ncminode_previter = ncminode
    !
    ! Walk tree to pick cmi-nodes
    !
    call extract_cminodes_from_tree(iruncmi,xyz_photosrc,nphotosrc,&
                                    node,ifirstincell,xyzh,&
                                    nxyzm_treetocmi,ixyzhm_leafparts,nnode_toreplace,&
                                    ncminode,nleafparts,ncloseleaf,&
                                    maxcminode,maxleafparts,&
                                    tree_accuracy_cmi,mintheta,rcut_opennode,rcut_leafpart,&
                                    nnode_needopen,nneedopen,&
                                    nxyzrs_nextopen_updatewalk,nnextopen_updatewalk)
    if (ncminode == 0) call fatal('photoionize_cmi','no nodes found')
    !
    ! Solve for the smoothing lengths of nodes
    !
    call hnode_iterate(iruncmi,node,nxyzm_treetocmi,ncminode,h_solvertocmi)
    !
    ! Combine nodes and individual particles near the source into a single set of sites
    !
    call collect_and_combine_cminodes(ncminode,nleafparts,ncloseleaf)   !- fill nixyzhmf_cminode(1:7,inode)
    if (ncminode > maxcminode) call fatal('photoionize_cmi','number of nodes exceeded maxcminode')
    !
    ! Stop and proceed to next iteration with larger rcut_opennode and rcut_leafpart
    ! if ncminode remained the same after opening the stored nodes
    !
    if (ncminode == ncminode_previter) then
       print*,' ** Likely that the nodes that need to be opened are already leaves **'
       print*,' ** Adjusting rcut_opennode and rcut_leafpart **'
       rcut_opennode = rcut_opennode + delta_rcut
       rcut_leafpart = rcut_leafpart + delta_rcut
       niter = niter + 1
       cycle resolve_ionfront
    endif
    !
    ! Pack final set of CMI input - fill x,y,z,h,m
    !
    call allocate_cmi_inoutputs(x,y,z,h,m,nH)
    if (crop_domain) then
       ncminode_in = 0
       r_crop2 = (rcrop_min + rcut_opennode*crop_fac)**2
       do inode = 1,ncminode
          if (mag2(nixyzhmf_cminode(3:5,inode)-cen_crop) < r_crop2) then
             ncminode_in = ncminode_in + 1 
             x(ncminode_in) = nixyzhmf_cminode(3,inode)
             y(ncminode_in) = nixyzhmf_cminode(4,inode)
             z(ncminode_in) = nixyzhmf_cminode(5,inode)
             h(ncminode_in) = nixyzhmf_cminode(6,inode)
             m(ncminode_in) = nixyzhmf_cminode(7,inode)
          endif 
       enddo 
       write(*,'(2x,a42,i8)') 'Number of pseudo-particles after cropping:',ncminode_in
       if (ncminode_in == 0) then 
          print*,' No particles left!'
          call deallocate_cmi_inoutputs(x,y,z,h,m,nH)
          print*,' ** Adjusting rcut_opennode and rcut_leafpart **'
          rcut_opennode = rcut_opennode + delta_rcut
          rcut_leafpart = rcut_leafpart + delta_rcut
          niter = niter + 1
          cycle resolve_ionfront
       endif 
    else 
       x(1:ncminode) = nixyzhmf_cminode(3,1:ncminode)
       y(1:ncminode) = nixyzhmf_cminode(4,1:ncminode)
       z(1:ncminode) = nixyzhmf_cminode(5,1:ncminode)
       h(1:ncminode) = nixyzhmf_cminode(6,1:ncminode)
       m(1:ncminode) = nixyzhmf_cminode(7,1:ncminode)
       ncminode_in = ncminode
    endif 
    !
    ! Run CMI
    !
    call run_cmacionize(ncminode_in,x,y,z,h,m,nH)
    !
    ! Store CMI output - fill nH
    !
    if (crop_domain) then
       !- Check that the whole ionized region is within the cropped domain 
       call check_cropped_space_ionization(ncminode_in,x,y,z,nH,must_iter)
       if (must_iter) then 
          print*,' Boundary of cropped region is ionized too!'
          call deallocate_cmi_inoutputs(x,y,z,h,m,nH)
          print*,' ** Adjusting rcut_opennode and rcut_leafpart **'
          rcut_opennode = rcut_opennode + delta_rcut
          rcut_leafpart = rcut_leafpart + delta_rcut
          niter = niter + 1
          cycle resolve_ionfront
       endif 
       !- If okay, store nH outputs, assuming that all nodes beyond cropped bounds are neutral
       ncminode_out = 0
       do inode = 1,ncminode
          if (mag2(nixyzhmf_cminode(3:5,inode)-cen_crop) < r_crop2) then
             ncminode_out = ncminode_out + 1 
             nixyzhmf_cminode(8,inode) = nH(ncminode_out) ! assume same order as input
          else 
             nixyzhmf_cminode(8,inode) = 1.  
          endif 
       enddo 
       if (ncminode_out /= ncminode_in) call fatal('photoionize_cmi','number of nodes coming in and out do not match!')
    else
       nixyzhmf_cminode(8,1:ncminode) = nH(1:ncminode)      !- fill nixyzhmf_cminode(8,inode)
    endif 
    call deallocate_cmi_inoutputs(x,y,z,h,m,nH)

    !- For plotting current ionization structure of nodes
    if (write_node_prop) then
       open(2029,file='nixyzhmf_cminode.txt')
       write(2029,'(a25,a25,a25,a25,a25,a25,a25,a25)') 'n','i','x','y','z','h','m','nH'
       do isite = 1,ncminode
          write(2029,*) nixyzhmf_cminode(1:8,isite)
       enddo
       close(2029)
    endif

    !
    ! Check the outputs to adjust tree-walk in the next iteration -
    !- If node is big but partially ionized (nH < 1), open node.
    !
    ! Also record the location & radius of regions (i.e. uppermost nodes) that have been iterated,
    ! and record the minimum node size amongst the descendent nodes within each region -
    !- Effectively 'saves' the final distribution of nodes by the end of the iteration to
    !  give to the next timestep without using node indices (which can change when tree rebuilds)
    !  but only their physical properties.
    !
    node_checks_passed = .true.
    if (auto_opennode) then
       size_root = node(1)%size
       over_cminode: do isite = 1,ncminode
          n = int(nixyzhmf_cminode(1,isite))
          if (n /= 0) then
             size_node = node(n)%size/size_root
             nH_limit = 1/nHlimit_fac * (-1/size_node + nHlimit_fac) !- nH limit increases with cell size
             nH_limit = max(nH_limit,tiny(nH_limit))

             if (size_node > min_nodesize_toflag) then
                nH_node = nixyzhmf_cminode(8,isite)
                if (nH_node < 0. .or. nH_node > 1.) write(*,'(2x,a5,i7,a15)') 'node ',n,' has invalid nH'

                if (nH_node < nH_limit) then
                   !- Store to nnode_needopen for next iteration
                   nneedopen = nneedopen + 1
                   if (nneedopen > maxnode_open) call fatal('photoionize_cmi','too many nodes &
                                                           & need to be resolved')
                   nnode_needopen(nneedopen) = n
                   node_checks_passed = .false.

                   !- Store to nxyzrs_nextopen for next timestep
                   if (niter == 0) then
                      nnextopen = nnextopen + 1
                      nxyzrs_nextopen(1,nnextopen)   = n
                      nxyzrs_nextopen(2:4,nnextopen) = node(n)%xcen(1:3)
                      nxyzrs_nextopen(5,nnextopen)   = node(n)%size
                      nxyzrs_nextopen(6,nnextopen)   = node(n)%size*0.95 !- slightly smaller to allow for minor changes
                   else
                      !- find which node stored in nxyzrs_nextopen is n a descendent of
                      do inextopen = 1,nnextopen
                         n_nextopen = int(nxyzrs_nextopen(1,inextopen))
                         if (int(real(n)/real(2**niter)) == n_nextopen) then
                            !- store minimum descendent-node size
                            nxyzrs_nextopen(6,inextopen) = min(nxyzrs_nextopen(6,inextopen),node(n)%size*0.95)
                         endif
                      enddo
                   endif
                   ! Note: The updated nxyzrs_nextopen is NOT given to the tree until the next timestep /
                   !       when iterations at current timestep have completed.
                endif
             endif
          elseif (n == 0) then !- is particle
             nH_part = nixyzhmf_cminode(8,isite)
             if (nH_part < 0. .or. nH_part > 1.) then
                write(*,'(2x,a9,i9,a15)') 'particle ',int(nixyzhmf_cminode(2,isite)),' has invalid nH'
             endif
          endif
       enddo over_cminode
       niter = niter + 1
       if (niter > maxiter) call fatal('photoionize_cmi','ionization front is far away - increase delta_rcut')
    endif

    ! Note: nnode_needopen accumulates ALL nodes not passing check throughout the iterations;
    !      otherwise during the next tree-walks, upper nodes will not be opened (based on original
    !      tree-opening criteria) and hence not reaching the ones beneath that need to be opened.

 enddo resolve_ionfront

 !
 ! Now sync the completed set of nxyzrs_nextopen with nxyzrs_nextopen_updatewalk
 ! to give to the tree during the next timestep
 !
 if (auto_opennode) then
    nnextopen_updatewalk = nnextopen
    nxyzrs_nextopen_updatewalk(:,:) = nxyzrs_nextopen(:,:)
 endif

 !
 ! Adjusting tree_accuracy_cmi to avoid discontinuities in sparseness of nodes
 ! between the leaves region and the higher-level nodes region
 !
 if (auto_tree_acc) then
    if (abs(mintheta-tree_accuracy_cmi) > tree_acc_tol) then
       tree_accuracy_cmi = mintheta
       write(*,'(2x,a33,f5.2)')'Adjusting tree-opening criteria: ',tree_accuracy_cmi
    endif
 endif

 print*,'treewalk_iterate done'

 !
 ! Shift tree back up whenever possible
 !
 if (auto_opennode) then
    icall = icall + 1
    if (icall == ncall_checktreewalk) then
       if (nnextopen_updatewalk > 0) then
          call remove_unnecessary_opennode(ncminode)
       endif
       icall = 0
    endif
 endif

end subroutine treewalk_run_cmi_iterate

!-----------------------------------------------------------------------------
!+
! Packing the final set of nodes to be then passed to CMI:
!  - Combine xyzm_nodes and their h obtained from hnode_iterate solver.
!  - Remove leaf nodes which were flagged (i.e. close to source) and replace them
!    with individual particles within those leaves.
! Fills the array nixyzhmf_cminode(:,:) -
!   n: node index;      n = 0 if is particle
!   i: particle index;  i = 0 if is node
!+
!-----------------------------------------------------------------------------
subroutine collect_and_combine_cminodes(ncminode,nleafparts,ncloseleaf)
 use utils_cmi, only:record_time 
 use io,        only:fatal
 integer, intent(in)    :: nleafparts,ncloseleaf
 integer, intent(inout) :: ncminode
 integer, parameter :: maxstore = 1E6
 integer :: ip,inode,istore,inode_toreplace
 integer :: n_fromtree,n_toreplace,ncminode_old
 real    :: nixyzhmf_cminode_store(7,maxstore)

 call record_time(iruncmi,'before_collectnodes')
 print*,' Collecting final set of pseduo-particles'

 ncminode_old = ncminode

 !- Combine outputs from kdtree_cmi and hnode_cmi to fill properties of nodes
 do inode = 1,ncminode
    nixyzhmf_cminode(1,inode) = nxyzm_treetocmi(1,inode)
    nixyzhmf_cminode(2,inode) = 0
    nixyzhmf_cminode(3,inode) = nxyzm_treetocmi(2,inode)
    nixyzhmf_cminode(4,inode) = nxyzm_treetocmi(3,inode)
    nixyzhmf_cminode(5,inode) = nxyzm_treetocmi(4,inode)
    nixyzhmf_cminode(6,inode) = h_solvertocmi(inode)
    nixyzhmf_cminode(7,inode) = nxyzm_treetocmi(5,inode)
 enddo

 !- Move entries into new array nixyzhmf_cminode_store if n is in nnode_toreplace
 istore = 0
 inode_toreplace = 1 
 do inode = 1,ncminode 
    n_fromtree  = int(nixyzhmf_cminode(1,inode))
    n_toreplace = nnode_toreplace(inode_toreplace)
    if (n_fromtree /= n_toreplace) then 
       istore = istore + 1 
       nixyzhmf_cminode_store(1:7,istore) = nixyzhmf_cminode(1:7,inode)
    else 
       inode_toreplace = inode_toreplace + 1 
    endif 
 enddo 
 ! check that istore = ncminode - ncloseleaf
 if (istore /= ncminode-ncloseleaf) call fatal('photoionize_cmi','wrong entries')
 
 !- Replace nixyzhmf_cminode with nixyzhmf_cminode_store
 nixyzhmf_cminode = 0. 
 ncminode = istore
 nixyzhmf_cminode(1:7,1:ncminode) = nixyzhmf_cminode_store(1:7,1:ncminode)

 !- Add leaf-particles onto the list and extend ncminode
 over_parts: do ip = 1,nleafparts
    ncminode = ncminode + 1
    if (ncminode > maxcminode) call fatal('photoionize_cmi','number of nodes+particles exceeded maxcminode')
    nixyzhmf_cminode(1,ncminode)   = 0
    nixyzhmf_cminode(2:7,ncminode) = ixyzhm_leafparts(1:6,ip)
 enddo over_parts

 write(*,'(2x,a33,i8)') 'Total number of pseudo-particles:',ncminode

 if (ncminode /= ncminode_old-ncloseleaf+nleafparts) call fatal('photoionize_cmi','final value of &
   & ncminode is inconsistent with inputs from kdtree_cmi')

 call record_time(iruncmi,'after_collectnodes')

end subroutine collect_and_combine_cminodes

!-----------------------------------------------------------------------------
!+
! Check all nodes that were previously opened to resolve into the ionization front;
! If they are no longer heavily ionized [determine this with the current nixyzhmf_cminode],
! move the tree back up by removing entries in nxyzrs_nextopen(_updatewalk).
! Note: rcut_opennode and rcut_leafpart will remain.
!+
!-----------------------------------------------------------------------------
subroutine remove_unnecessary_opennode(ncminode)
 use allocutils, only:allocate_array
 use utils_cmi,  only:mag2 
 integer, intent(in)  :: ncminode
 real,    allocatable :: nHmin(:)
 integer :: nnextopen_old,inextopen,isite,irowafter
 real    :: pos_opennode(3),rad_opennode2,pos_site(3),dist2,nH_site
 real    :: nH,size,nH_limit,nHlimit_fac_rev

 nnextopen_old = nnextopen
 call allocate_array('nHmin',nHmin,nnextopen_old)

 !
 ! Find the lowest nH in each of the stored regions in nxyzrs_nextopen (which can overlap)
 !
 all_opened: do inextopen = 1,nnextopen_old
    nHmin(inextopen) = 1. !- init
    pos_opennode  = nxyzrs_nextopen(2:4,inextopen)
    rad_opennode2 = nxyzrs_nextopen(5,inextopen)**2
    current_nodes: do isite = 1,ncminode
       pos_site = nixyzhmf_cminode(3:5,isite)
       dist2 = mag2(pos_site(1:3)-pos_opennode(1:3))
       if (dist2 < rad_opennode2) then
          nH_site = nixyzhmf_cminode(8,isite)
          nHmin(inextopen) = min(nHmin(inextopen),nH_site)
       endif
    enddo current_nodes
 enddo all_opened

 !
 ! Remove entry in nxyzrs_nextopen(_updatewalk) if nH is now sufficiently high
 !
 nHlimit_fac_rev = nHlimit_fac + 5  !- slightly higher than the opening-criterion
 inextopen = 1
 do while (inextopen <= nnextopen)
    nH = nHmin(inextopen)
    size = nxyzrs_nextopen(5,inextopen)
    nH_limit = 1./nHlimit_fac_rev * (-1./size + nHlimit_fac_rev)
    if (nH > nH_limit .and. size > min_nodesize_toflag) then  !- can undo the resolve-iteration
       do irowafter = inextopen,nnextopen-1
          nxyzrs_nextopen(1:6,irowafter) = nxyzrs_nextopen(1:6,irowafter+1)
          nHmin(irowafter) = nHmin(irowafter+1)
       enddo
       nnextopen = nnextopen - 1
    else
       inextopen = inextopen + 1
    endif
 enddo
 nnextopen_updatewalk = nnextopen
 nxyzrs_nextopen_updatewalk(:,:) = nxyzrs_nextopen(:,:)

 if (nnextopen < nnextopen_old) then
    write(*,'(1x,a8,i5,a37)') 'Removed ', nnextopen_old-nnextopen,' nodes stored in nxyzrs_nextopen'
 endif

 if (allocated(nHmin)) deallocate(nHmin)

end subroutine remove_unnecessary_opennode

!-----------------------------------------------------------------------------
!+
! If only passing a cropped spherical region of the whole simulation domain to CMI,
! check that the boundary is beyond ionization front i.e. not ionized.
! Note: Could happen during the initial iterations of resolve-ionfront iter, so 
!       kill the simulation only if the warning persisted across timesteps
!+
!-----------------------------------------------------------------------------
subroutine check_cropped_space_ionization(nsite,x,y,z,nH,must_iter)
 use utils_cmi,  only:quick_sort,mag2 
 integer, intent(in)  :: nsite 
 real,    intent(in)  :: x(nsite),y(nsite),z(nsite),nH(nsite)
 logical, intent(out) :: must_iter   !- flag to immediately jump to next resolve_ionfront iteration
 integer :: isite,inbound,outbound,ifurthest
 integer :: isite_all(nsite)         !- array of site indices to be sorted 
 real    :: dist2site_all(nsite)     !- array of sites' distances to centre to be sorted 
 real    :: boundary_frac = 0.05
 real    :: disttocen2,nH_site
 logical :: bounds_ok

 if (plot_cropped_sites) then 
    open(2039,file='xyzf_cminode_cropping.txt')
    write(2039,'(a25,a25,a25,a25)') 'x','y','z','nH'
    do isite = 1,nsite
       write(2039,*) x(isite), y(isite), z(isite), nH(isite) 
    enddo
    close(2039)
 endif 

 must_iter = .false. 

 !- Compute distances between each node and centre of sources and store in array 
 !$omp parallel do default(none) shared(nsite,cen_crop,x,y,z,isite_all,dist2site_all) &
 !$omp private(isite,disttocen2) &
 !$omp schedule(runtime)
 do isite = 1,nsite
    disttocen2 = mag2((/ x(isite),y(isite),z(isite) /) - cen_crop)
    isite_all(isite) = isite
    dist2site_all(isite) = disttocen2
 enddo 
 !$omp end parallel do

 !- Sort dist2_site_all array in ascending order along with isite_all 
 call quick_sort(nsite,dist2site_all,isite_all,1,nsite)

 !- Get the furthest sites
 bounds_ok = .false.
 do while (.not.bounds_ok)
    inbound  = int((1.-boundary_frac)*nsite)
    outbound = nsite
    if (outbound-inbound > 5) then 
       bounds_ok = .true.
    else 
       boundary_frac = boundary_frac + 0.05
    endif 
 enddo

 !- Check nH of these sites 
 check_nH: do ifurthest = 1,outbound-inbound+1 
    isite = isite_all(inbound+ifurthest-1) 
    nH_site = nH(isite)
    if (nH_site < 1.) then 
       must_iter = .true. 
       exit check_nH
    endif 
 enddo check_nH

end subroutine check_cropped_space_ionization

!-----------------------------------------------------------------------------
!+
! Pass the positions, masses and smoothing lengths of given sites (nodes/particles)
! to CMI and run MCRT simulation
!+
!-----------------------------------------------------------------------------
subroutine run_cmacionize(nsite,x,y,z,h,m,nH)
 use utils_cmi, only:record_time
 use omp_lib
 integer, intent(in)  :: nsite
 real,    intent(in)  :: x(nsite),y(nsite),z(nsite),h(nsite),m(nsite)
 real,    intent(out) :: nH(nsite)
 integer :: talk,numthreads
 real    :: cputime1,cputime2,walltime1,walltime2

 call cpu_time(cputime1)
 walltime1 = omp_get_wtime()

 call record_time(iruncmi,'beforeCMI')

 !
 ! Initialize CMI
 !
 call write_cmi_infiles(nsite,x,y,z,h,m)
 talk = 0
 if (print_cmi) talk = 1
 !$omp parallel default(none) shared(numthreads)
 numthreads = omp_get_num_threads()
 !$omp end parallel
 call cmi_init("phantom_cmi.param",numthreads,udist_si,umass_si,"Petkova",talk)

 !
 ! Run CMI
 !
 if (.not.print_cmi) write(*,'(2x,a15,i3,a8)') 'Running CMI on ',numthreads,' threads'
 call cmi_compute_neutral_fraction_dp(x(1:nsite),y(1:nsite),z(1:nsite),h(1:nsite),&
                                      m(1:nsite),nH(1:nsite),int8(nsite))
 call cmi_destroy


 call cpu_time(cputime2)
 walltime2 = omp_get_wtime()

 open(2060,file='CMI_cpu_wall_time_record.txt',position='append')
 write(2060,*) iruncmi, cputime2-cputime1, walltime2-walltime1
 close(2060)

 call record_time(iruncmi,'afterCMI')

end subroutine run_cmacionize

!-----------------------------------------------------------------------------
!+
! Writes the Voronoi sites and the .param input files to be read by CMI
!+
!-----------------------------------------------------------------------------
subroutine write_cmi_infiles(nsite_in,x_in,y_in,z_in,h_in,m_in)
 use utils_cmi, only:modify_grid,set_bounds
 use io,        only:fatal
 integer, intent(in) :: nsite_in
 real,    intent(in) :: x_in(:),y_in(:),z_in(:),h_in(:),m_in(:)
 integer :: i,isrc,nsite
 real    :: xmin,xmax,ymin,ymax,zmin,zmax,dx,dy,dz,space_fac
 real    :: xmin_si,ymin_si,zmin_si,dx_si,dy_si,dz_si
 real,   allocatable :: x(:),y(:),z(:)
 logical :: redo_grid

 !- Modify these properties only within this subroutine; avoid intent(in) error 
 nsite = nsite_in   
 allocate(x(nsite))
 allocate(y(nsite))
 allocate(z(nsite))
 x = x_in
 y = y_in
 z = z_in
 if (limit_voronoi) call modify_grid(iruncmi,nsite,x,y,z,h_in,hlimit_fac,extradist_fac)

 !- Set boundaries to only around the given sites
 call set_bounds(nsite,x,y,z,h_in,m_in,xmin,xmax,ymin,ymax,zmin,zmax,dx,dy,dz)

 !- Add some space to the edges and convert to SI units
 space_fac = 1E-2
 xmin_si = (xmin - space_fac*dx) *udist_si
 ymin_si = (ymin - space_fac*dy) *udist_si
 zmin_si = (zmin - space_fac*dz) *udist_si
 dx_si   = dx*(1. + 2*space_fac) *udist_si
 dy_si   = dy*(1. + 2*space_fac) *udist_si
 dz_si   = dz*(1. + 2*space_fac) *udist_si

 !
 ! Rewrite the Voronoi grid generating-sites only if their locations have shifted
 ! beyond the tolerance
 ! (But cmi-node labelling could change across timesteps even if it was the same particle....)
 !
 redo_grid = .false.
 check_poschange: do i = 1,maxcminode
    if (x_old(i) /= 0. .and. x(i) /= 0. .and. y_old(i) /= 0. .and. y(i) /= 0. .and. z_old(i) /= 0. .and. z(i) /= 0.) then
       if (abs(x_old(i)-x(i))/(xmax-xmin) > tol_vsite .or. abs(y_old(i)-y(i))/(ymax-ymin) > tol_vsite .or. &
           abs(z_old(i)-z(i))/(zmax-zmin) > tol_vsite) then
          redo_grid = .true.
          exit check_poschange
       endif
    endif
 enddo check_poschange

 if (redo_grid .or. first_call) then
    write(*,'(2x,a33)') 'Updating voronoi generating sites'
    open(2024,file='phantom_voronoi_sites.txt')
    nsite_lastgrid = nsite  ! for cases where grids from prev steps are retained
    do i = 1,nsite
       write(2024,*) x(i)*udist_si, y(i)*udist_si, z(i)*udist_si
    enddo
    close(2024)
    !- update xyz history
    call reset_cmi_inputs_history
    do i = 1,nsite
       x_old(i) = x(i)
       y_old(i) = y(i)
       z_old(i) = z(i)
    enddo
    first_call = .false.
 endif

 !
 ! Write file containing position and flux of sources
 !
 open(2025,file='source_positions.txt')
 write(2025,*) nphotosrc
 !- total ionizing photon flux Q of all sources
 totq = 0.
 do isrc = 1,nphotosrc
    totq = totq + ionflux_photosrc(isrc)
 enddo
 write(2025,*) totq
 do isrc = 1,nphotosrc
    write(2025,*) xyz_photosrc_si(1:3,isrc), ionflux_photosrc(isrc)/totq
 enddo
 close(2025)

 !
 ! Write .param file - Note: Sub-entries have to be indented by 2 spaces
 !
 open(2026,file='phantom_cmi.param')

 write(2026,'(a15)') "AbundanceModel:"
 write(2026,'(2x,a6)') "C: 0.0"
 write(2026,'(2x,a7)') "He: 0.0"
 write(2026,'(2x,a7)') "Ne: 0.0"
 write(2026,'(2x,a6)') "N: 0.0"
 write(2026,'(2x,a6)') "O: 0.0"
 write(2026,'(2x,a6)') "S: 0.0"

 write(2026,'(a22)') "TemperatureCalculator:"
 write(2026,'(2x,a32)') "do temperature calculation: true"

 write(2026,'(a25)') "DiffuseReemissionHandler:"
 write(2026,'(2x,a10)') "type: None"

 write(2026,'(a14)') "SimulationBox:"
 write(2026,'(2x,a9,e12.3,a4,e12.3,a4,e12.3,a3)') "anchor: [",xmin_si," m, ",ymin_si," m, ",zmin_si," m]"
 write(2026,'(2x,a8,e12.3,a4,e12.3,a4,e12.3,a3)') "sides: [",dx_si," m, ",dy_si," m, ",dz_si," m]"
 write(2026,'(2x,a34)') "periodicity: [false, false, false]"

 write(2026,'(a12)') "DensityGrid:"
 write(2026,'(2x,a13)') "type: Voronoi"
 write(2026,'(2x,a14)') "grid type: Old"

 if (lloyd) write(2026,'(2x,a29)') "number of Lloyd iterations: 5"

 write(2026,'(2x,a29)') "VoronoiGeneratorDistribution:"
 write(2026,'(4x,a21,i10)') "number of positions: ",nsite_lastgrid
 write(2026,'(4x,a35)') "filename: phantom_voronoi_sites.txt"
 write(2026,'(4x,a9)') "type: SPH"

 if (write_gamma) then
    write(2026,'(a24)') "DensityGridWriterFields:"
    write(2026,'(2x,a15)') "HeatingRateH: 1"
 endif

 write(2026,'(a16)') "DensityFunction:"
 write(2026,'(2x,a17)') "type: Homogeneous"

 write(2026,'(a25)') "PhotonSourceDistribution:"
 write(2026,'(2x,a20)') "type: AsciiFileTable"
 write(2026,'(2x,a30)') "filename: source_positions.txt"

 write(2026,'(a21)') "PhotonSourceSpectrum:"
 if (monochrom_source) then
    write(2026,'(2x,a19)') "type: Monochromatic"
    write(2026,'(2x,a11,e12.3,a3)') "frequency: ",freq_photon," Hz"
 else
    write(2026,'(2x,a12)') "type: Planck"
    write(2026,'(2x,a13,e10.3,a2)') "temperature: ",temp_star," K"
    write(2026,'(2x,a28)') "ionizing flux: -1. m^-2 s^-1"
 endif

 write(2026,'(a21)') "IonizationSimulation:"
 write(2026,'(2x,a22,i2)') "number of iterations: ",niter_mcrt
 write(2026,'(2x,a19,i7)') "number of photons: ",nphoton

 write(2026,'(a14)') "CrossSections:"
 write(2026,'(2x,a16)') "type: FixedValue"
 write(2026,'(2x,a24)') "hydrogen_0: 6.3e-18 cm^2"
 write(2026,'(2x,a16)') "helium_0: 0. m^2"
 write(2026,'(2x,a16)') "carbon_1: 0. m^2"
 write(2026,'(2x,a16)') "carbon_2: 0. m^2"
 write(2026,'(2x,a18)') "nitrogen_0: 0. m^2"
 write(2026,'(2x,a18)') "nitrogen_1: 0. m^2"
 write(2026,'(2x,a18)') "nitrogen_2: 0. m^2"
 write(2026,'(2x,a16)') "oxygen_0: 0. m^2"
 write(2026,'(2x,a16)') "oxygen_1: 0. m^2"
 write(2026,'(2x,a14)') "neon_0: 0. m^2"
 write(2026,'(2x,a14)') "neon_1: 0. m^2"
 write(2026,'(2x,a17)') "sulphur_1: 0. m^2"
 write(2026,'(2x,a17)') "sulphur_2: 0. m^2"
 write(2026,'(2x,a17)') "sulphur_3: 0. m^2"

 write(2026,'(a19)') "RecombinationRates:"
 write(2026,'(2x,a16)') "type: FixedValue"
 write(2026,'(2x,a29)') "hydrogen_1: 2.7e-13 cm^3 s^-1"
 write(2026,'(2x,a21)') "helium_1: 0. m^3 s^-1"
 write(2026,'(2x,a21)') "carbon_2: 0. m^3 s^-1"
 write(2026,'(2x,a21)') "carbon_3: 0. m^3 s^-1"
 write(2026,'(2x,a23)') "nitrogen_1: 0. m^3 s^-1"
 write(2026,'(2x,a23)') "nitrogen_2: 0. m^3 s^-1"
 write(2026,'(2x,a23)') "nitrogen_3: 0. m^3 s^-1"
 write(2026,'(2x,a21)') "oxygen_1: 0. m^3 s^-1"
 write(2026,'(2x,a21)') "oxygen_2: 0. m^3 s^-1"
 write(2026,'(2x,a19)') "neon_1: 0. m^3 s^-1"
 write(2026,'(2x,a19)') "neon_2: 0. m^3 s^-1"
 write(2026,'(2x,a22)') "sulphur_2: 0. m^3 s^-1"
 write(2026,'(2x,a22)') "sulphur_3: 0. m^3 s^-1"
 write(2026,'(2x,a22)') "sulphur_4: 0. m^3 s^-1"

 close(2026)
 deallocate(x)
 deallocate(y)
 deallocate(z)

end subroutine write_cmi_infiles

!-----------------------------------------------------------------------------
!+
! Routines for writing node/particle properties to snapshot files
!+
!-----------------------------------------------------------------------------
subroutine init_write_snapshot
 use utils_cmi,  only:gen_filename 
 integer :: ifile_search
 logical :: lastfile_found,iexist
 character(len=50) :: filename_search,filename

 iunit = 4000
 iruncmi = 0
 warned = .false.
 !- Set starting ifile by finding the last snapshot file saved
 ifile_search = maxoutfile
 lastfile_found = .false.
 do while (.not.lastfile_found)
    call gen_filename(photoionize_tree,ifile_search,filename_search)
    inquire(file=trim(filename_search),exist=iexist)
    if (iexist) then
        lastfile_found = .true.
    else
        ifile_search = ifile_search - 1
        if (ifile_search < 0) then
           exit  !- start from _00000
        endif
    endif
 enddo
 ifile = ifile_search + 1
 call gen_filename(photoionize_tree,ifile,filename)
 print*,'Beginning with snapshot file ',trim(filename)

end subroutine init_write_snapshot


subroutine write_nH_snapshot(time,nsite,xyzh_parts,x_in,y_in,z_in,h_in,m_in,nH_in)
 use part,      only:isdead_or_accreted
 use utils_cmi, only:gen_filename 
 use units,     only:utime
 use io,        only:fatal,warning
 integer, intent(in) :: nsite
 real,    intent(in) :: time
 real,    intent(in), optional :: xyzh_parts(:,:)
 real,    intent(in), optional :: x_in(nsite),y_in(nsite),z_in(nsite),h_in(nsite)
 real,    intent(in), optional :: m_in(nsite),nH_in(nsite)
 integer :: isite,ip,ip_cmi
 logical :: got_all
 character(len=50) :: allsites_filename

 !- Checking inputs
 if (.not.photoionize_tree) then
    got_all = .false.
    if (present(xyzh_parts) .and. present(x_in) .and. present(y_in).and. present(z_in) &
        .and. present(h_in) .and. present(m_in) .and. present(nH_in)) then
       got_all = .true.
    endif
    if (.not.got_all) call fatal('photoionize_cmi','Missing inputs for writing snapshots')
 endif

 if (ifile < maxoutfile) then
    if (mod(iruncmi,ncall_writefile) == 0) then

       call gen_filename(photoionize_tree,ifile,allsites_filename)
       print*,'Writing outputs to ',trim(allsites_filename)

       open(unit=iunit,file=trim(allsites_filename),status='replace')

       write(iunit,*) 'time_cgs'
       write(iunit,*) time*utime

       if (photoionize_tree) then
          write(iunit,'(a20,a20,a25,a25,a25,a25,a25,a25)') 'n','i','x','y','z','h','m','nH'
          do isite = 1,nsite
             write(iunit,*) nixyzhmf_cminode(1:8,isite)
          enddo
       else
          write(iunit,'(a20,a25,a25,a25,a25,a25)') 'x','y','z','h','m','nH'
          ip_cmi = 0
          do ip = 1,nsite
             if (.not.isdead_or_accreted(xyzh_parts(4,ip))) then
                ip_cmi = ip_cmi + 1
                write(iunit,*) x_in(ip_cmi),y_in(ip_cmi),z_in(ip_cmi),h_in(ip_cmi),m_in(ip_cmi),nH_in(ip_cmi)
             endif
          enddo
       endif
       close(iunit)

       iunit = iunit + 1
       ifile = ifile + 1
    endif
 else
    if (.not.warned) then
       call warning('photoionize_cmi','maxoutfile reached - no longer writing node info to file')
       warned = .true.
    endif
 endif

end subroutine write_nH_snapshot

!-----------------------------------------------------------------------
!+
!  Routines to init/clear storage arrays during runtime
!+
!-----------------------------------------------------------------------
subroutine allocate_cminode
 use allocutils, only:allocate_array
 use dim,        only:maxp
 call allocate_array('nxyzm_treetocmi',nxyzm_treetocmi,5,maxcminode+1)
 call allocate_array('h_solvertocmi',h_solvertocmi,maxcminode+1)
 call allocate_array('ixyzhm_leafparts',ixyzhm_leafparts,6,maxleafparts+1)
 call allocate_array('nnode_toreplace',nnode_toreplace,maxcminode+1)
 call allocate_array('nixyzhmf_cminode',nixyzhmf_cminode,8,maxcminode+1)
 call allocate_array('nnode_needopen',nnode_needopen,maxnode_open)
 call allocate_array('nxyzrs_nextopen',nxyzrs_nextopen,6,maxnode_open)
 call allocate_array('nxyzrs_nextopen_updatewalk',nxyzrs_nextopen_updatewalk,6,maxnode_open)
 call allocate_array('vxyzu_beforepred',vxyzu_beforepred,4,maxp)
 call allocate_array('du_cmi',du_cmi,maxp)
 call allocate_array('dudt_cmi',dudt_cmi,maxp)
 call allocate_array('nH_allparts',nH_allparts,maxp)
end subroutine allocate_cminode

subroutine reset_cminode
 nxyzm_treetocmi  = 0.
 h_solvertocmi    = 0.
 nnode_toreplace  = 0
 nixyzhmf_cminode = 0.
 ixyzhm_leafparts = 0.
end subroutine reset_cminode

subroutine allocate_cmi_inputs_history
 use allocutils,   only:allocate_array
 call allocate_array('x_old',x_old,maxcminode+1)
 call allocate_array('y_old',y_old,maxcminode+1)
 call allocate_array('z_old',z_old,maxcminode+1)
end subroutine allocate_cmi_inputs_history

subroutine reset_cmi_inputs_history
 x_old = 0.
 y_old = 0.
 z_old = 0.
end subroutine reset_cmi_inputs_history

subroutine reset_energ
 vxyzu_beforepred = 0.
 du_cmi   = 0.
 dudt_cmi = 0.
 nH_allparts = -1.  !- -1 denotes those NOT to be heated
end subroutine reset_energ

subroutine allocate_cmi_inoutputs(x,y,z,h,m,nH)
 use allocutils, only:allocate_array
 real, allocatable, intent(inout) :: x(:),y(:),z(:),h(:),m(:),nH(:)
 call allocate_array('x',x,maxcminode+1)
 call allocate_array('y',y,maxcminode+1)
 call allocate_array('z',z,maxcminode+1)
 call allocate_array('h',h,maxcminode+1)
 call allocate_array('m',m,maxcminode+1)
 call allocate_array('nH',nH,maxcminode+1)
end subroutine allocate_cmi_inoutputs

subroutine deallocate_cmi_inoutputs(x,y,z,h,m,nH)
 real, allocatable, intent(inout) :: x(:),y(:),z(:),h(:),m(:),nH(:)
 if (allocated(x))  deallocate(x)
 if (allocated(y))  deallocate(y)
 if (allocated(z))  deallocate(z)
 if (allocated(h))  deallocate(h)
 if (allocated(m))  deallocate(m)
 if (allocated(nH)) deallocate(nH)
end subroutine deallocate_cmi_inoutputs

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_photoionize(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling photoionization'
 call write_inopt(inject_rad,'inject_rad','Inject radiation',iunit)
 call write_inopt(sink_ionsrc,'sink_ionsrc','Using sinks as ionizing sources',iunit)
 call write_inopt(one_sink_ionsrc,'one_sink_ionsrc','Using one specific sink as ionizing source',iunit)
 call write_inopt(isink_ionsrc,'isink_ionsrc','Sink label for ionizing source',iunit)
 call write_inopt(sink_as_cluster,'sink_as_cluster','Sink represent clusters',iunit)
 call write_inopt(masscrit_ionize_cgs,'masscrit_ionize_cgs','Critical sink mass to begin emitting radiation',iunit)
 call write_inopt(niter_mcrt,'niter_mcrt','Number of photon-release iterations',iunit)
 call write_inopt(nphoton,'nphoton','Number of photons per iteration',iunit)
 call write_inopt(photon_eV,'photon_eV','Energy of ionizing photons in eV',iunit)
 call write_inopt(force_large_Q,'force_large_Q','Forcefully set a bigger Q',iunit)
 call write_inopt(tol_vsite,'tol_vsite','Threshold to update Voronoi gen-sites',iunit)
 call write_inopt(lloyd,'lloyd','Apply Lloyd iteration to construct Voronoi grid',iunit)
 call write_inopt(limit_voronoi,'limit_voronoi','Merge tightly packed cells in Voronoi grid',iunit)
 call write_inopt(hlimit_fac,'hlimit_fac','Threshold fraction in h to merge',iunit)
 call write_inopt(extradist_fac,'extradist_fac','Radius factor to merge',iunit)
 call write_inopt(monochrom_source,'monochrom_source','Use monochromatic source',iunit)
 call write_inopt(photoionize_tree,'photoionize_tree','Pass tree nodes as pseudo-particles to CMI',iunit)
 call write_inopt(tree_accuracy_cmi,'tree_accuracy_cmi','Tree-opening criteria for cmi-nodes',iunit)
 call write_inopt(rcut_opennode_cgs,'rcut_opennode_cgs','Radius within which nodes must be leaves',iunit)
 call write_inopt(rcut_leafpart_cgs,'rcut_leafpart_cgs','Radius within which particles are extracted',iunit)
 call write_inopt(auto_opennode,'auto_opennode','Automatically adjust tree-walk with iterative approach',iunit)
 call write_inopt(auto_tree_acc,'auto_tree_acc','Automatically adjust tree_accuracy_cmi',iunit)
 call write_inopt(delta_rcut_cgs,'delta_rcut_cgs','Increase in rcut_opennode_cgs and rcut_leafpart_cgs per iter step',iunit)
 call write_inopt(nHlimit_fac,'nHlimit_fac','Paramter controlling resolution of ionization front',iunit)
 call write_inopt(min_nodesize_toflag,'min_nodesize_toflag','Minimum node size to check nH',iunit)
 call write_inopt(crop_domain,'crop_domain','Crop simulation domain to pass to CMI',iunit)
 call write_inopt(crop_fac,'crop_fac','Factor of rcut_leaf to crop',iunit)
 call write_inopt(temp_hii,'temp_hii','Temperature of ionized gas',iunit)
 call write_inopt(fix_temp_hii,'fix_temp_hii','Heat ionized particles to specified temp_hii',iunit)
 call write_inopt(implicit_cmi,'implicit_cmi','Use implicit method to update internal energies',iunit)
 call write_inopt(treat_Rtype_phase,'treat_Rtype_phase','Check that dt covers duration of R-type phase',iunit)

end subroutine write_options_photoionize

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_photoionize(name,valstring,imatch,igotall,ierr)
 use io,    only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer,          save        :: ngot = 0
 character(len=30), parameter  :: label = 'read_options_photoionize'

 imatch  = .true.
 igotall = .false.

 select case(trim(name))
 case('inject_rad')
    read(valstring,*,iostat=ierr) inject_rad
    ngot = ngot + 1
 case('sink_ionsrc')
    read(valstring,*,iostat=ierr) sink_ionsrc
    ngot = ngot + 1
 case('one_sink_ionsrc')
    read(valstring,*,iostat=ierr) one_sink_ionsrc
    ngot = ngot + 1
 case('isink_ionsrc')
    read(valstring,*,iostat=ierr) isink_ionsrc
    ngot = ngot + 1
    if (isink_ionsrc <= 0.) call fatal(label,'invalid setting for isink_ionsrc (<=0)')
 case('sink_as_cluster')
    read(valstring,*,iostat=ierr) sink_as_cluster
    ngot = ngot + 1
 case('masscrit_ionize_cgs')
    read(valstring,*,iostat=ierr) masscrit_ionize_cgs
    ngot = ngot + 1
    if (masscrit_ionize_cgs <= 0.) call fatal(label,'invalid setting for masscrit_ionize_cgs (<=0)')
 case('niter_mcrt')
    read(valstring,*,iostat=ierr) niter_mcrt
    ngot = ngot + 1
    if (niter_mcrt <= 0.) call fatal(label,'invalid setting for niter_mcrt (<=0)')
 case('nphoton')
    read(valstring,*,iostat=ierr) nphoton
    ngot = ngot + 1
    if (nphoton <= 0.) call fatal(label,'invalid setting for nphoton (<=0)')
 case('photon_eV')
    read(valstring,*,iostat=ierr) photon_eV
    ngot = ngot + 1
    if (photon_eV <= 0.) call fatal(label,'invalid setting for photon_eV (<=0)')
 case('force_large_Q')
    read(valstring,*,iostat=ierr) force_large_Q
    ngot = ngot + 1
 case('tol_vsite')
    read(valstring,*,iostat=ierr) tol_vsite
    ngot = ngot + 1
    if (tol_vsite < 0.) call fatal(label,'invalid setting for tol_vsite (<0)')
 case('lloyd')
    read(valstring,*,iostat=ierr) lloyd
    ngot = ngot + 1
 case('limit_voronoi')
    read(valstring,*,iostat=ierr) limit_voronoi
    ngot = ngot + 1
 case('hlimit_fac')
    read(valstring,*,iostat=ierr) hlimit_fac
    ngot = ngot + 1
    if (hlimit_fac < 0.) call fatal(label,'invalid setting for hlimit_fac (<0)')
 case('extradist_fac')
    read(valstring,*,iostat=ierr) extradist_fac
    ngot = ngot + 1
    if (extradist_fac < 1.) call fatal(label,'invalid setting for extradist_fac (<1)')
 case('monochrom_source')
    read(valstring,*,iostat=ierr) monochrom_source
    ngot = ngot + 1
 case('photoionize_tree')
    read(valstring,*,iostat=ierr) photoionize_tree
    ngot = ngot + 1
 case('tree_accuracy_cmi')
    read(valstring,*,iostat=ierr) tree_accuracy_cmi
    ngot = ngot + 1
    if (tree_accuracy_cmi < 0. .or. tree_accuracy_cmi > 1.) call fatal(label,'invalid setting for tree_accuracy_cmi')
 case('rcut_opennode_cgs')
    read(valstring,*,iostat=ierr) rcut_opennode_cgs
    ngot = ngot + 1
    if (rcut_opennode_cgs < 0.) call fatal(label,'invalid setting for rcut_opennode_cgs')
 case('rcut_leafpart_cgs')
    read(valstring,*,iostat=ierr) rcut_leafpart_cgs
    if (rcut_leafpart_cgs < 0.) call fatal(label,'invalid setting for rcut_leafpart_cgs')
    ngot = ngot + 1
 case('auto_opennode')
    read(valstring,*,iostat=ierr) auto_opennode
    ngot = ngot + 1
 case('auto_tree_acc')
    read(valstring,*,iostat=ierr) auto_tree_acc
    ngot = ngot + 1
 case('delta_rcut_cgs')
    read(valstring,*,iostat=ierr) delta_rcut_cgs
    ngot = ngot + 1
    if (delta_rcut_cgs < 0.) call fatal(label,'invalid setting for delta_rcut_cgs')
 case('nHlimit_fac')
    read(valstring,*,iostat=ierr) nHlimit_fac
    ngot = ngot + 1
    if (nHlimit_fac < 0.) call fatal(label,'invalid setting for nHlimit_fac')
 case('min_nodesize_toflag')
    read(valstring,*,iostat=ierr) min_nodesize_toflag
    ngot = ngot + 1
    if (min_nodesize_toflag < 0.) call fatal(label,'invalid setting for min_nodesize_toflag')
 case('crop_domain')
    read(valstring,*,iostat=ierr) crop_domain
    ngot = ngot + 1
 case('crop_fac')
    read(valstring,*,iostat=ierr) crop_fac
    ngot = ngot + 1
    if (crop_fac <= 0.) call fatal(label,'invalid setting for crop_fac')
 case('temp_hii')
    read(valstring,*,iostat=ierr) temp_hii
    ngot = ngot + 1
    if (temp_hii <= 0.) call fatal(label,'invalid setting for temp_hii (<=0)')
 case('fix_temp_hii')
    read(valstring,*,iostat=ierr) fix_temp_hii
    ngot = ngot + 1
 case('implicit_cmi')
    read(valstring,*,iostat=ierr) implicit_cmi
    ngot = ngot + 1
 case('treat_Rtype_phase')
    read(valstring,*,iostat=ierr) treat_Rtype_phase
    ngot = ngot + 1
 case default
    imatch = .false.
 end select
 igotall = ( ngot >= 31 )

end subroutine read_options_photoionize

end module photoionize_cmi
