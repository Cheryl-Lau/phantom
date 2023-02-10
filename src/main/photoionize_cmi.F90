!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module photoionize_cmi
!
! The CMI suite: *photoionize_cmi.f90* kdtree_cmi.f90 hnode_cmi.f90
! This module contains all the subroutines necessary for doing photoionization
! using the Monte Carlo Radiative Transfer code CMacIonize
!
! Two options are available -
!  1. Pass all particles to CMI for density-mapping and construct grid cell around each
!     individual particle; or
!  2. Extract higher-level tree nodes relative to the given source location(s) and pass
!     them as pseudo-particles to CMI for both density-mapping and grid-construction.
!     Individual particles can be retained in regions which are close to the source.
!     This code is capable of automatically adjusting the tree-walk to ensure that the
!     ionization front is resolved. Activate auto_opennode and tweak parameters including
!     nHlimit_fac, min_nodesize_toflag and tree_accuracy_cmi for best results.
!
! :References: Petkova,et.al,2021,MNRAS,507,858
!
! :Owner: Cheryl Lau (adapted from Maya Petkova's analysis mod)
!
! :Runtime parameters:
!   - sink_ionsrc         : *Using sinks as ionizing sources*
!   - masscrit_ionize     : *Minimum mass of sink that emit ionizing radiation*
!   - niter_mcrt          : *Number of iterations to release photon packets in MCRT simulation*
!   - nphoton             : *Number of photons per iteration*
!   - wavelength          : *Wavelength of ionizing photons in nm*
!   - temp_hii            : *Temperature of ionized gas*
!   - tol_vsite           : *Tolerence in nodes' position change above which the Voronoi generating-
!                            sites will be updated*
!   - lloyd               : *Do Lloyd iterations for Voronoi grid construction*
!   - photoionize_tree    : *Construct Voronoi grid with tree nodes; else with particles*
!   - tree_accuracy_cmi   : *threshold angle to open tree node during kdtree-walk*
!   - rcut_opennode       : *Radius within which we must get leaves*
!   - rcut_leafpart       : *Radius within which we extract individual particles*
!   - auto_opennode       : *Automatically adjust tree-walk according to CMI neutral frac output*
!   - auto_tree_acc       : *Automatically adjust tree_accuracy_cmi based on runtime outputs*
!   - delta_rcut          : *Increase in rcut_opennode and rcut_leafpart in an iteration step*
!   - nHlimit_fac         : *Parameter which controls the resolution of the ionization front*
!   - min_nodesize_toflag : *Minimum node size (relative to root node) to check neutral frac*
!
! :Dependencies: infile_utils, physcon, units, io, dim, boundaries, eos, part, kdtree
!
!
 use cmi_fortran_library

 implicit none

 public :: init_ionizing_radiation_cmi,set_ionizing_source_cmi,release_ionizing_radiation_cmi
 public :: read_options_photoionize,write_options_photoionize

 ! Position of sources emitting radiation at current time
 integer, public, parameter :: maxphotosrc = 10
 integer, public :: nphotosrc                    !- Current number of sources
 real   , public :: xyz_photosrc(3,maxphotosrc)

 ! Using sinks as source
 logical, public :: sink_ionsrc = .false.
 real,    public :: masscrit_ionize = 7.  ! in code units
 !- or
 ! Manually set location and starting/ending time of ionizing sources
 integer, public, parameter :: nsetphotosrc = 1
 real,    public :: xyzt_setphotosrc(5,nsetphotosrc) = reshape((/0.,0.,0.,1E-50,1E50 /),&
                                                               shape=(/5,nsetphotosrc/))
 ! Monte Carlo simulation settings
 integer, public :: nphoton    = 1E6
 integer, public :: niter_mcrt = 10
 real,    public :: wavelength = 500     ! nm
 real,    public :: temp_hii   = 1E4     ! K
 real,    public :: tol_vsite  = 1E-9    ! in code units
 logical, public :: lloyd      = .true.

 ! Move grid-construction up the tree
 logical, public :: photoionize_tree = .true.

 ! Options for extracting cmi-nodes from kdtree
 real,    public :: tree_accuracy_cmi = 0.1
 real,    public :: rcut_opennode = 0.2 !0.10        ! in code units
 real,    public :: rcut_leafpart = 0.1 !0.05        ! in code units
 real,    public :: delta_rcut    = 0.05 ! 0.01        ! in code units
 real,    public :: nHlimit_fac   = 100 !80          ! ionization front resolution; recommend 40-80
 real,    public :: min_nodesize_toflag = 0.005  ! min node size as a fraction of root node
 logical, public :: auto_opennode = .true.
 logical, public :: auto_tree_acc = .false.

 private

 ! Init arrays to store properties of nodes
 integer, parameter   :: maxcminode   = 1E8
 integer, parameter   :: maxleafparts = 1E8
 real,    allocatable :: nxyzm_treetocmi(:,:),ixyzhm_leafparts(:,:)
 integer, allocatable :: nnode_toreplace(:)
 real,    allocatable :: h_solvertocmi(:)

 ! Init array to store properties of all sites
 real,    allocatable :: nixyzhmf_cminode(:,:)

 ! Init arrays to control subsequent tree-walks
 integer, parameter   :: maxnode_open = 1E6
 integer, allocatable :: nnode_needopen(:)    !- for next iteration
 real,    allocatable :: nxyzrs_nextopen(:,:) !- for next timestep
 real,    allocatable :: nxyzrs_nextopen_updatewalk(:,:)
 integer :: nnextopen,nnextopen_updatewalk
 integer :: ncminode_previter

 real, allocatable :: x_old(:),y_old(:),z_old(:)
 integer :: nsite_lastgrid

 integer, parameter :: maxoutfile_ult = 99999
 integer :: maxoutfile = 5000          !- max number of xyzhmnH output files
 integer :: ncall_writefile  = 10     !- interval to write xyzhmnH output file
 integer :: icall,iunit,ifile,iruncmi
 integer :: ncall_checktreewalk = 100
 real    :: xyz_photosrc_si(3,maxphotosrc),lumin_photosrc(maxphotosrc)
 real    :: tree_accuracy_cmi_old
 real    :: solarl_photonsec,freq_photon,u_hii,udist_si,umass_si
 real    :: time0_wall,time0_cpu,time_now_wall,time_now_cpu
 real    :: time_ellapsed_wall,time_ellapsed_cpu
 logical :: first_call,warned
 logical :: print_cmi = .false.   ! Print CMI outputs to shell

 real :: u_hi ! testing

contains

!----------------------------------------------------------------
!+
! Check params and initialize variables to prepare for CMI call
!+
!----------------------------------------------------------------
subroutine init_ionizing_radiation_cmi(npart,xyzh)
 use physcon,  only:mass_proton_cgs,kboltz,c,planckh,solarl
 use eos,      only:gmw,gamma
 use io,       only:warning,fatal
 use units,    only:udist,umass,unit_ergg
 use dim,      only:maxvxyzu
 use omp_lib
 integer, intent(in) :: npart
 real,    intent(in) :: xyzh(:,:)
 integer :: inode,ip,io_file,ifile_search
 real    :: h_avg,psep,wavelength_cgs,energ_photon
 real    :: gmw0,csi_cgs,temp_hii_fromcsi
 logical :: compilecond_ok,lastfile_found,iexist
 character(len=20) :: ifile_search_char,nixyzhmf_search_filename
 character(len=20) :: xyzhmf_search_filename

 print*,'Radiation-hydrodynamics: Phantom is coupled to photoionization code CMacIonize'
 print*,'Injecting ionizing radiation with MCRT method'

 if (.not.photoionize_tree .and. maxcminode < npart) call fatal('photoionize_cmi',&
                                                   & 'maxcminode has to be greater than npart')
 if (maxvxyzu < 4) call fatal('photoionize_cmi','Not suitable for isothermal simulations')

 !- Check that the tol_vsite is sensible by estimating mean particle separation
 h_avg = 0.
 do ip = 1,npart
    h_avg = h_avg + xyzh(4,ip)
 enddo
 h_avg = h_avg/npart
 psep = 2*h_avg / (57.9)**(1./3.)
 if (tol_vsite > 0.1*psep) call warning('photoionize_cmi','tol_vsite might be too large')

 !- Check pick-nodes settings
 if (photoionize_tree) then
    if (tree_accuracy_cmi == 0) call warning('photoionize_cmi','extracting only leaf nodes')
    if (rcut_opennode < rcut_leafpart) call fatal('photoionize_cmi','rcut_leafpart must be &
                                                 & smaller than rcut_opennode')
    compilecond_ok = .false.
#ifndef PERIODIC
#ifndef MPI
    compilecond_ok = .true.
#endif
#endif
    if (.not.compilecond_ok) call fatal('photoionize_cmi','current version does not support &
                                       & MPI or PERIODIC')
 endif

 !- Convert cgs units to SI units
 udist_si = udist/1E2
 umass_si = umass/1E3

 !- Init memory
 call allocate_cminode
 call reset_cminode
 call allocate_cmi_inputs_history
 call reset_cmi_inputs_history

 ! Init array to store memory for tree-walk
 nnextopen = 0
 nnextopen_updatewalk = 0
 do inode = 1,maxnode_open
    nxyzrs_nextopen(1:6,inode) = (/ 0.,0.,0.,0.,0.,0. /)
    nxyzrs_nextopen_updatewalk(1:6,inode) = (/ 0.,0.,0.,0.,0.,0. /)
 enddo
 tree_accuracy_cmi_old = tree_accuracy_cmi

 !- Internal energy u of ionized particles
 gmw0 = gmw
 gmw = 0.5  !- temporarily change the mean molecular weight of ionized particles
 u_hii = kboltz * temp_hii / (gmw*mass_proton_cgs*(gamma-1.)) /unit_ergg
 csi_cgs = sqrt(temp_hii*(gamma*kboltz)/(gmw*mass_proton_cgs))
 temp_hii_fromcsi = (12.85E5)**2*gmw*mass_proton_cgs/(gamma*kboltz)
 print*,' -Ionized gas properties- '
 print*,' internal energy: ',u_hii*unit_ergg,'erg/g'
 print*,' sound speed:     ',csi_cgs,'cm/s'
 print*,' temperature:     ',temp_hii_fromcsi,'K'
 gmw = gmw0  !- reset

 u_hi = kboltz * 1E2 / (gmw*mass_proton_cgs*(gamma-1.))/unit_ergg  ! testing

 !- Calculate solar luminosity [photon/sec]
 wavelength_cgs = wavelength*1E-7
 energ_photon   = planckh*c/wavelength_cgs
 solarl_photonsec = solarl/energ_photon

 !- Calculate photon frequency [Hz]
 freq_photon = c/wavelength_cgs

 !- For writing xyzhmf output files
 if (maxoutfile > maxoutfile_ult) then
    call fatal('photoionize_cmi','maxoutfile must be less than 5 digits, or modify line ')
 endif
 iunit = 4000
 iruncmi = 0
 warned = .false.

 if (photoionize_tree) then
    !- Set starting val of ifile by finding the last nixyzhmf_* saved
    ifile_search = 1000
    lastfile_found = .false.
    do while (.not.lastfile_found)
       write(ifile_search_char,'(i5.5)') ifile_search
       nixyzhmf_search_filename = trim('nixyzhmf_')//trim(ifile_search_char)//trim('.txt')
       inquire(file=nixyzhmf_search_filename,exist=iexist)
       if (iexist) then
          lastfile_found = .true.
       else
          ifile_search = ifile_search - 1
          if (ifile_search < 0) then
             exit
          endif
       endif
    enddo
    ifile = ifile_search + 1
    print*,'Start with file nixyzhmf_',ifile
 else
    !- Set starting val of ifile by finding the last xyzhmf_* saved
    ifile_search = 1000
    lastfile_found = .false.
    do while (.not.lastfile_found)
       write(ifile_search_char,'(i5.5)') ifile_search
       xyzhmf_search_filename = trim('xyzhmf_')//trim(ifile_search_char)//trim('.txt')
       inquire(file=xyzhmf_search_filename,exist=iexist)
       if (iexist) then
          lastfile_found = .true.
       else
          ifile_search = ifile_search - 1
          if (ifile_search < 0) then
             exit
          endif
       endif
    enddo
    ifile = ifile_search + 1
    print*,'Start with file xyzhmf_',ifile
 endif

 !- For writing Voronoi sites files / checking openned nodes
 icall = 0
 first_call = .true.

 !- Time the simulations
 call cpu_time(time0_cpu)
 time0_wall = omp_get_wtime()
 open(2050,file='cpu_wall_time_record.txt',status='replace',iostat=io_file)
 if (io_file /= 0) call fatal('photoionize_cmi','unable to open time-record file')

end subroutine init_ionizing_radiation_cmi

!----------------------------------------------------------------
!+
! Dynamically set the locations of the sources at current timestep
! Updates nphotosrc and xyz_photosrc
!+
!----------------------------------------------------------------
subroutine set_ionizing_source_cmi(time,nptmass,xyzmh_ptmass)
 use io,       only:fatal
 use physcon,  only:solarm
 use units,    only:umass
 integer, intent(in) :: nptmass
 real,    intent(in) :: time
 real,    intent(in) :: xyzmh_ptmass(:,:)
 integer :: isrc,ix,isink
 real    :: time_startsrc,time_endsrc,mptmass,lumin,mass_8solarm

 !- Init/Reset
 nphotosrc = 0
 do isrc = 1,maxphotosrc
    xyz_photosrc(1:3,isrc) = (/ 0.,0.,0. /)
    lumin_photosrc(isrc)   = 0.
 enddo
 !
 ! Extract sources which currently emit radiation
 !
 if (sink_ionsrc .and. nptmass > 0.) then !- Check sink
    do isink = 1,nptmass
       mptmass = xyzmh_ptmass(4,isink)
       if (mptmass >= masscrit_ionize) then
          nphotosrc = nphotosrc + 1
          if (nphotosrc > maxphotosrc) call fatal('photoionize_cmi','number of sources &
                                                 & exceeded maxphotosrc')
          xyz_photosrc(1:3,nphotosrc) = xyzmh_ptmass(1:3,isink)
          !- Compute luminosity with sink mass
          lumin_photosrc(nphotosrc) = get_lumin(mptmass)
       endif
    enddo
 elseif (.not.sink_ionsrc) then !- Check time
    do isrc = 1,nsetphotosrc
       time_startsrc = xyzt_setphotosrc(4,isrc)
       time_endsrc   = xyzt_setphotosrc(5,isrc)
       if (time > time_startsrc .and. time < time_endsrc) then
          nphotosrc = nphotosrc + 1
          if (nphotosrc > maxphotosrc) call fatal('photoionize_cmi','number of sources &
                                                 & exceeded maxphotosrc')
          xyz_photosrc(1:3,nphotosrc) = xyzt_setphotosrc(1:3,isrc)
          !- Assume fixed luminosity
          mass_8solarm = (8.*solarm)/umass
          lumin_photosrc(nphotosrc) = get_lumin(mass_8solarm)
       endif
    enddo
 endif

 !- Convert to SI units for CMI param file
 if (nphotosrc >= 1) then
    do isrc = 1,maxphotosrc
       do ix = 1,3
          xyz_photosrc_si(ix,isrc) = xyz_photosrc(ix,isrc)*udist_si
       enddo
    enddo
 endif

end subroutine set_ionizing_source_cmi

!
! Calculate luminosity [s^-1] from mass of sink particle
!
real function get_lumin(mptmass)
 use physcon,  only:solarm
 use units,    only:umass
 real, intent(in)  :: mptmass
 real    :: mptmass_cgs

 mptmass_cgs = mptmass*umass
 get_lumin   = solarl_photonsec * (mptmass_cgs/solarm)**(3.5)

end function get_lumin

!----------------------------------------------------------------
!+
! Wrapper for updating particle energies
! (assumes particle is ionized if neutral fraction nH falls below 0.5)
!+
!----------------------------------------------------------------
subroutine release_ionizing_radiation_cmi(time,npart,xyzh,vxyzu)
 use part,     only:massoftype,igas,isdead_or_accreted
 use kdtree,   only:inodeparts,inoderange
 use units,    only:utime
 use io,       only:fatal,warning
 use omp_lib
 integer, intent(inout) :: npart
 real,    intent(in)    :: time
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,    allocatable   :: x(:),y(:),z(:),h(:),m(:),nH(:)
 integer :: ip,ip_cmi,npart_cmi,i,n,ipnode,isite,ncminode
 real    :: nH_part,nH_site
 real    :: u_ionized  ! testing
 character(len=50) :: xyzhmf_parts_filename,ifile_char

 if (nphotosrc >= 1) then
    if (photoionize_tree) then
       !- Walk tree and run CMI
       call treewalk_run_cmi_iterate(time,xyzh,ncminode)

       !- Update u of all particles beneath each node using the final nixyzhmf_cminode
       over_entries: do isite = 1,ncminode
          nH_site = nixyzhmf_cminode(8,isite)
          if (nH_site < 0. .or. nH_site > 1.) call fatal('photoionize_cmi','invalid nH')
          if (nH_site < 0.5) then
             n = int(nixyzhmf_cminode(1,isite))
             i = int(nixyzhmf_cminode(2,isite))
             if (n /= 0 .and. i == 0) then !- is node
                over_parts: do ipnode = inoderange(1,n),inoderange(2,n)
                   ip = abs(inodeparts(ipnode))
                   if (vxyzu(4,ip) < u_hii) then
                      vxyzu(4,ip) = u_hii
                   endif
!                   u_ionized = u_hii*(1.0-nH_site) + u_hi*nH_site  ! testing
!                   if (vxyzu(4,ip) < u_ionized) vxyzu(4,ip) = u_ionized
                enddo over_parts
             elseif (n == 0 .and. i /= 0) then !- is particle
                if (vxyzu(4,i) < u_hii) then
                   vxyzu(4,i) = u_hii
                endif
!                u_ionized = u_hii*(1.0-nH_site) + u_hi*nH_site  ! testing
!                if (vxyzu(4,i) < u_ionized) vxyzu(4,i) = u_ionized
             else
                call fatal('photoionize_cmi','unidentified site')
             endif
          endif
       enddo over_entries
       call reset_cminode

    else ! Pass all particles to grid-construction and density mapping
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

       !- Update u of particles using the computed neutral frac
       ip_cmi = 0
       do ip = 1,npart
          if (.not.isdead_or_accreted(xyzh(4,ip))) then
             ip_cmi = ip_cmi + 1
             nH_part = nH(ip_cmi)
             if (nH_part < 0.) call fatal('photoionize_cmi','invalid nH')
             if (nH_part < 0.5) then
                if (vxyzu(4,ip) < u_hii) then
                   vxyzu(4,ip) = u_hii
                endif
             endif
          endif
       enddo 
       if (ip_cmi /= npart_cmi) call fatal('photoionize_cmi','number of particles &
         & passed to and from CMI do not match')

       !- Write xyzhmf of particle to a snapshot file
       if (iunit < maxoutfile) then
          if (mod(iruncmi,ncall_writefile) == 0) then
             write(ifile_char,'(i5.5)') ifile
             xyzhmf_parts_filename = trim('xyzhmf_')//trim(ifile_char)//trim('.txt')
             print*,'Writing outputs to ',xyzhmf_parts_filename
             open(unit=iunit,file=xyzhmf_parts_filename,status='replace')
             write(iunit,*) 'time_cgs'
             write(iunit,*) time*utime
             write(iunit,'(a20,a25,a25,a25,a25,a25)') 'x','y','z','h','m','nH'
             ip_cmi = 0
             do ip = 1,npart
                if (.not.isdead_or_accreted(xyzh(4,ip))) then
                   ip_cmi = ip_cmi + 1
                   write(iunit,*) x(ip_cmi),y(ip_cmi),z(ip_cmi),h(ip_cmi),m(ip_cmi),nH(ip_cmi)
                endif
             enddo
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

       call deallocate_cmi_inoutputs(x,y,z,h,m,nH)
    endif

    iruncmi = iruncmi + 1

    !- Time the simulations
    call cpu_time(time_now_cpu)
    time_now_wall = omp_get_wtime()
    time_ellapsed_wall = time_now_wall - time0_wall
    time_ellapsed_cpu  = time_now_cpu  - time0_cpu
    open(2050,file='cpu_wall_time_record.txt',position='append')
    write(2050,*) iruncmi, time_ellapsed_cpu, time_ellapsed_wall
    close(2050)
 endif

end subroutine release_ionizing_radiation_cmi

!----------------------------------------------------------------
!+
! Iterate tree-walk and CMI-call until the ionization front is resolved
! Gives the final set of nixyzhmf_cminode
!+
!----------------------------------------------------------------
subroutine treewalk_run_cmi_iterate(time,xyzh,ncminode)
 use linklist,   only:node,ifirstincell
 use kdtree_cmi, only:extract_cminodes_from_tree
 use hnode_cmi,  only:hnode_iterate
 use units,      only:utime
 use io,         only:fatal,warning
 real,    intent(in)  :: time
 real,    intent(in)  :: xyzh(:,:)
 integer, intent(out) :: ncminode   !- total number of CMI sites (nodes+particles)
 real,    allocatable :: x(:),y(:),z(:),h(:),m(:),nH(:)
 integer :: maxiter = 50
 integer :: nneedopen    !- number of nodes at ionization front that needs to be opened at current iteration
 integer :: ncloseleaf   !- number of leaves to be replaced
 integer :: nleafparts   !- total number of particles in leaf nodes to be replaced
 integer :: niter,inode,isite,n,i,inextopen,n_nextopen,maxnextopen
 real    :: size_node,size_root,nH_node,nH_part,nH_limit,tree_accuracy_cmi_new
 logical :: node_checks_passed
 character(len=50) :: nixyzhmf_cminode_filename,ifile_char

 !- Init
 node_checks_passed = .false.
 niter = 0
 nneedopen = 0
 do inode = 1,maxnode_open
    nnode_needopen(inode) = 0
 enddo
 !
 ! Shift tree back up whenever possible
 !
 if (auto_opennode) then
    icall = icall + 1
    if (icall == ncall_checktreewalk) then
       if (nnextopen_updatewalk > 0) then
          call remove_unnecessary_opennode(nHlimit_fac,ncminode,min_nodesize_toflag)
       endif
       icall = 0
    endif
 endif

 resolve_ionfront: do while (.not.node_checks_passed)
    write(*,'(1x,a35,i2)') 'Resolve ionization front iteration ',niter
    ncminode_previter = ncminode
    !
    ! Walk tree to pick cmi-nodes
    !
    call extract_cminodes_from_tree(xyz_photosrc,nphotosrc,&
                                    node,ifirstincell,xyzh,&
                                    nxyzm_treetocmi,ixyzhm_leafparts,nnode_toreplace,&
                                    ncminode,nleafparts,ncloseleaf,&
                                    maxcminode,maxleafparts,&
                                    tree_accuracy_cmi,rcut_opennode,rcut_leafpart,&
                                    nnode_needopen,nneedopen,&
                                    nxyzrs_nextopen_updatewalk,nnextopen_updatewalk)
    if (ncminode == 0) call fatal('photoionize_cmi','no nodes found')

    ! testing
!    open(2033,file='nxyzm_treetocmi.txt')
!    write(2033,*) 'ncminode',ncminode
!    do i = 1,ncminode
!       write(2033,*) nxyzm_treetocmi(1:5,i)
!    enddo
!    close(2033)

    !
    ! Solve for the smoothing lengths of nodes
    !
    call hnode_iterate(node,nxyzm_treetocmi,ncminode,h_solvertocmi)
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
    ! Run CMI
    !
    call allocate_cmi_inoutputs(x,y,z,h,m,nH)
    x(1:ncminode) = nixyzhmf_cminode(3,1:ncminode)
    y(1:ncminode) = nixyzhmf_cminode(4,1:ncminode)
    z(1:ncminode) = nixyzhmf_cminode(5,1:ncminode)
    h(1:ncminode) = nixyzhmf_cminode(6,1:ncminode)
    m(1:ncminode) = nixyzhmf_cminode(7,1:ncminode)
    call run_cmacionize(ncminode,x,y,z,h,m,nH)
    !
    ! Store CMI output
    !
    nixyzhmf_cminode(8,1:ncminode) = nH(1:ncminode) !- fill nixyzhmf_cminode(8,inode)
    call deallocate_cmi_inoutputs(x,y,z,h,m,nH)

    ! testing
    open(2029,file='nixyzhmf_cminode.txt')
    write(2029,'(a25,a25,a25,a25,a25,a25,a25,a25)') 'n','i','x','y','z','h','m','nH'
    do isite = 1,ncminode
       write(2029,*) nixyzhmf_cminode(1:8,isite)
    enddo
    close(2029)

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
                if (nH_node < 0. .or. nH_node > 1.) then
                   write(*,'(2x,a5,i7,a15)') 'node ',n,' has invalid nH'
                endif

                if (nH_node < nH_limit) then
                   !- Store to nnode_needopen for next iteration
                   nneedopen = nneedopen + 1
                   if (nneedopen > maxnode_open) call fatal('photoionize_cmi','too many nodes &
                                                 need to be resolved')
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

    !- testing
!    open(2031,file='nnode_needopen.txt')
!    write(2031,*) nneedopen
!    do i = 1,nneedopen
!       write(2031,*) nnode_needopen(i)
!    enddo
!    close(2031)

    ! Note: nnode_needopen accumulates ALL nodes not passing check throughout the iterations;
    !      otherwise during the next tree-walks, upper nodes will not be opened (based on original
    !      tree-opening criteria) and hence not reaching the ones beneath that need to be opened.

 enddo resolve_ionfront

 !- testing
! open(2032,file='nxyzrs_nextopen.txt')
! write(2032,*) nnextopen
! do i = 1,nnextopen
!    write(2032,*) nxyzrs_nextopen(1:6,i)
! enddo
! close(2032)

 !
 ! Now sync the completed set of nxyzrs_nextopen with nxyzrs_nextopen_updatewalk
 ! to give to the tree during the next timestep
 !
 if (auto_opennode) then
    nnextopen_updatewalk = nnextopen
    nxyzrs_nextopen_updatewalk(:,:) = nxyzrs_nextopen(:,:)
 endif
 !
 ! Lower tree_accuracy_cmi if too many nodes need to be opened
 ! likely that the radiation is escaping in all directions and reaching large distances
 !
 if (auto_tree_acc) then
    maxnextopen = 1E3   !- define a max value for nnextopen to parametrize tree_accuracy
    if (nnextopen >= maxnextopen) call fatal('photoionize_cmi','require a larger maxnextopen')
    if (nnextopen > 300) then
       tree_accuracy_cmi_new = tree_accuracy_cmi_old*(-0.4/maxnextopen*nnextopen+0.9)
       if (tree_accuracy_cmi_new < tree_accuracy_cmi) then
          tree_accuracy_cmi = tree_accuracy_cmi_new
          write(*,'(2x,a33,f5.2)')'Adjusting tree-opening criteria: ',tree_accuracy_cmi
       endif
    endif
 endif

 print*,'treewalk_iterate done'

 !
 ! Write time and nixyzhmf_cminode(1:8,1:ncminode) to a snapshot file
 !
 if (iunit < maxoutfile) then
    if (mod(iruncmi,ncall_writefile) == 0) then
       write(ifile_char,'(i5.5)') ifile
       nixyzhmf_cminode_filename = trim('nixyzhmf_')//trim(ifile_char)//trim('.txt')
       print*,'Writing outputs to ',nixyzhmf_cminode_filename
       open(unit=iunit,file=nixyzhmf_cminode_filename,status='replace')
       write(iunit,*) 'time_cgs'
       write(iunit,*) time*utime
       write(iunit,'(a20,a20,a25,a25,a25,a25,a25,a25)') 'n','i','x','y','z','h','m','nH'
       do isite = 1,ncminode
          write(iunit,*) nixyzhmf_cminode(1:8,isite)
       enddo
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

end subroutine treewalk_run_cmi_iterate

!----------------------------------------------------------------
!+
! Packing the final set of nodes to be then passed to CMI:
!  - Combine xyzm_nodes and their h obtained from hnode_iterate solver.
!  - Remove leaf nodes which were flagged (i.e. close to source) and replace them
!    with individual particles within those leaves.
! Fills the array nixyzhmf_cminode(:,:) -
!   n: node index;      n = 0 if is particle
!   i: particle index;  i = 0 if is node
!+
!----------------------------------------------------------------
subroutine collect_and_combine_cminodes(ncminode,nleafparts,ncloseleaf)
 use io,  only:fatal
 integer, intent(in)    :: nleafparts,ncloseleaf
 integer, intent(inout) :: ncminode
 integer :: ip,inode,ientry,irepentry,irowafter,ncminode_old,ncminode_curr
 integer :: n_fromtree,n_toreplace,i

 ! testing
! open(2028,file='nodespicked_cmi.txt')
! write(2028,*) 'ncminode',ncminode
! write(2028,*) 'ncloseleaf',ncloseleaf
! write(2028,*) 'nleafparts',nleafparts
! write(2028,'(a20,a20)') 'n','nleaf'
! do i = 1,ncminode
!    write(2028,*) nxyzm_treetocmi(1,i), nnode_toreplace(i)
! enddo
! close(2028)

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

 !- Match n against those in nnode_toreplace and remove entries by shifting up the table
 !  (assumes both arrays are sorted in the same order, which should be.)
 ncminode_old = ncminode
 ientry = 1
 irepentry = 1
 over_entries: do while (ncminode > ncminode_old-ncloseleaf)
    n_fromtree  = int(nixyzhmf_cminode(1,ientry))
    n_toreplace = nnode_toreplace(irepentry)
    if (n_fromtree == n_toreplace) then
       do irowafter = ientry,ncminode-1
          nixyzhmf_cminode(1:7,irowafter) = nixyzhmf_cminode(1:7,irowafter+1)
       enddo
       ncminode = ncminode - 1
       irepentry = irepentry + 1
       if (irepentry > ncloseleaf+1) call fatal('photoionize_cmi','erroneous irepentry')
    else
       ientry = ientry + 1
    endif
    if (ientry > ncminode_old) call fatal('photoionize_cmi','erroneous ientry')
 enddo over_entries

 ! testing
! open(2030,file='nodesremoved_cmi.txt')
! write(2030,*) 'ncminode',ncminode
! write(2030,'(a20,a20)') 'n','nleaf'
! do i = 1,ncminode
!    write(2030,*) nixyzhmf_cminode(1,i), nnode_toreplace(i)
! enddo
! close(2030)

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

 ! testing
! open(2031,file='partsadded_cmi.txt')
! write(2031,*) 'ncminode',ncminode
! write(2031,'(a20,a20)') 'inodes','nleaf'
! do i = 1,ncminode
!    write(2031,*) nixyzhmf_cminode(2,i), nnode_toreplace(i)
! enddo
! write(2031,*) 'particle indices'
! do i = 1,nleafparts
!    write(2031,*) ixyzhm_leafparts(1,i)
! enddo
! close(2031)

end subroutine collect_and_combine_cminodes

!----------------------------------------------------------------
!+
! Check all nodes which were previously opened to resolve into the ionization front;
! If they are no longer heavily ionized [determine this with the currrent nixyzhmf_cminode],
! move the tree back up by removing entries in nxyzrs_nextopen(_updatewalk).
!+
!----------------------------------------------------------------
subroutine remove_unnecessary_opennode(nHlimit_fac,ncminode,min_nodesize_toflag)
 use allocutils,  only:allocate_array
 use linklist,    only:node
 integer, intent(in)    :: ncminode
 real,    intent(in)    :: nHlimit_fac,min_nodesize_toflag
 real,    allocatable   :: nHmin(:)
 integer :: nnextopen_old,inextopen,icminode,n_node,irowafter
 real    :: pos_opennode(3),rad_opennode2,pos_node(3),dist2,nH_node
 real    :: nH,size,nH_limit,nHlimit_fac_rev

 nnextopen_old = nnextopen
 call allocate_array('nHmin',nHmin,nnextopen_old )
 !
 ! Find the lowest nH in each of the stored regions in nxyzrs_nextopen (which can overlap)
 !
 all_opened: do inextopen = 1,nnextopen_old
    nHmin(inextopen) = 1. !- init
    pos_opennode = nxyzrs_nextopen(2:4,inextopen)
    rad_opennode2 = nxyzrs_nextopen(5,inextopen)**2.
    current_nodes: do icminode = 1,ncminode
       n_node = int(nixyzhmf_cminode(1,icminode))
       if (n_node /= 0) then
          pos_node = nixyzhmf_cminode(3:5,icminode)
          dist2 = mag2(pos_node(1:3)-pos_opennode(1:3))
          if (dist2 < rad_opennode2) then
             nH_node = nixyzhmf_cminode(8,icminode)
             nHmin(inextopen) = min(nHmin(inextopen),nH_node)
          endif
       endif
    enddo current_nodes
 enddo all_opened
 !
 ! Remove entry in nxyzrs_nextopen(_updatewalk) if nH is sufficiently high
 !
 nHlimit_fac_rev = nHlimit_fac + 5.  !- slightly higher than the opening-criterion
 inextopen = 1
 do while (inextopen <= nnextopen)
    nH = nHmin(inextopen)
    size = nxyzrs_nextopen(5,inextopen)
    nH_limit = 1/nHlimit_fac_rev * (-1/size + nHlimit_fac)
    if (nH > nH_limit .and. size > min_nodesize_toflag) then !- can undo the resolve-iteration
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

 print*,'nnextopen_old/nnextopen: ',nnextopen_old,nnextopen

 if (allocated(nHmin)) deallocate(nHmin)

end subroutine remove_unnecessary_opennode


real function mag2(vec)
 real, intent(in) :: vec(3)

 mag2 = vec(1)**2 + vec(2)**2 + vec(3)**2

end function mag2

!----------------------------------------------------------------
!+
! Pass the positions, masses and smoothing lengths of given sites
! to CMI and run MCRT simulation;
! This applies to both tree nodes and particles.
!+
!----------------------------------------------------------------
subroutine run_cmacionize(nsite,x,y,z,h,m,nH)
 use omp_lib
 integer, intent(in)  :: nsite
 real,    intent(in)  :: x(nsite),y(nsite),z(nsite),h(nsite),m(nsite)
 real,    intent(out) :: nH(nsite)
 integer :: talk,numthreads,isite

 !- Testing
! open(2027,file='xyzhm_cmi.txt')
! write(2027,*) nsite
! do isite = 1,nsite
!    write(2027,*) x(isite),y(isite),z(isite),h(isite),m(isite)
! enddo
! close(2027)
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
 if (.not.print_cmi) write(*,'(2x,a11)') 'Running CMI'
 call cmi_compute_neutral_fraction_dp(x(1:nsite),y(1:nsite),z(1:nsite),h(1:nsite),&
                                      m(1:nsite),nH(1:nsite),int8(nsite))
 call cmi_destroy

 !- Testing
! open(2077,file='xyzhmf_cmi.txt')
! write(2077,*) nsite
! do isite = 1,nsite
!    write(2077,*) x(isite),y(isite),z(isite),h(isite),m(isite),nH(isite)
! enddo
! close(2077)

end subroutine run_cmacionize

!-----------------------------------------------------------------------
!+
! Writes the Voronoi sites and the .param input files to be read by CMI
!+
!-----------------------------------------------------------------------
subroutine write_cmi_infiles(nsite,x,y,z,h,m)
 use io,  only:fatal
 integer, intent(in) :: nsite
 real,    intent(in) :: x(:),y(:),z(:),h(:),m(:)
 integer :: i,isrc
 real    :: xmin,xmax,ymin,ymax,zmin,zmax,dx,dy,dz,space_fac
 real    :: xmin_si,ymin_si,zmin_si,dx_si,dy_si,dz_si
 real    :: totlumin
 logical :: redo_grid
 !
 ! Set boundaries only around the given sites
 ! also check that the input xyzhm are sensible
 !
 xmax = -huge(xmax)
 ymax = -huge(ymax)
 zmax = -huge(zmax)
 xmin =  huge(xmin)
 ymin =  huge(ymin)
 zmin =  huge(zmin)
 !$omp parallel do default(none) shared(nsite,x,y,z,h,m) private(i) &
 !$omp reduction(min:xmin,ymin,zmin) &
 !$omp reduction(max:xmax,ymax,zmax)
 do i = 1,nsite
    if (x(i) /= x(i) .or. y(i) /= y(i) .or. z(i) /= z(i) .or. &
        h(i) < tiny(h) .or. m(i) < tiny(m)) then
       print*,x(i),y(i),z(i),h(i),m(i)
       call fatal('photoionize_cmi','invalid xyzhm input')
    endif
    xmin = min(xmin,x(i))
    ymin = min(ymin,y(i))
    zmin = min(zmin,z(i))
    xmax = max(xmax,x(i))
    ymax = max(ymax,y(i))
    zmax = max(zmax,z(i))
 enddo
 !$omp end parallel do
 dx = abs(xmax - xmin)
 dy = abs(ymax - ymin)
 dz = abs(zmax - zmin)

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
 !
 redo_grid = .false.
 check_poschange: do i = 1,maxcminode
    if (abs(x_old(i)-x(i)) > tol_vsite .or. abs(y_old(i)-y(i)) > tol_vsite .or. &
        abs(z_old(i)-z(i)) > tol_vsite) then
       redo_grid = .true.
       exit check_poschange
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
 ! Write file containing position and luminosity of sources
 !
 open(2025,file='source_positions.txt')
 write(2025,*) nphotosrc
 totlumin = 0.
 do isrc = 1,nphotosrc
    totlumin = totlumin + lumin_photosrc(isrc)
 enddo
 write(2025,*) totlumin
 do isrc = 1,nphotosrc
    write(2025,*) xyz_photosrc_si(1:3,isrc), lumin_photosrc(isrc)/totlumin
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
 write(2026,'(2x,a33)') "do temperature calculation: false"

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

 write(2026,'(a16)') "DensityFunction:"
 write(2026,'(2x,a17)') "type: Homogeneous"

 write(2026,'(a25)') "PhotonSourceDistribution:"
 write(2026,'(2x,a20)') "type: AsciiFileTable"
 write(2026,'(2x,a30)') "filename: source_positions.txt"

 write(2026,'(a21)') "PhotonSourceSpectrum:"
 write(2026,'(2x,a19)') "type: Monochromatic"
 write(2026,'(2x,a11,e12.3,a3)') "frequency: ",freq_photon," Hz"

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

end subroutine write_cmi_infiles

!-----------------------------------------------------------------------
!+
!  Subroutines to init/clear storage arrays during runtime
!+
!-----------------------------------------------------------------------
subroutine allocate_cminode
 use allocutils, only:allocate_array

 call allocate_array('nxyzm_treetocmi',nxyzm_treetocmi,5,maxcminode+1)
 call allocate_array('h_solvertocmi',h_solvertocmi,maxcminode+1)
 call allocate_array('ixyzhm_leafparts',ixyzhm_leafparts,6,maxleafparts+1)
 call allocate_array('nnode_toreplace',nnode_toreplace,maxcminode+1)
 call allocate_array('nixyzhmf_cminode',nixyzhmf_cminode,8,maxcminode+1)
 call allocate_array('nnode_needopen',nnode_needopen,maxnode_open)
 call allocate_array('nxyzrs_nextopen',nxyzrs_nextopen,6,maxnode_open)
 call allocate_array('nxyzrs_nextopen_updatewalk',nxyzrs_nextopen_updatewalk,6,maxnode_open)

end subroutine allocate_cminode


subroutine reset_cminode
 integer :: inode,ip

 !$omp parallel do default(none) shared(nxyzm_treetocmi,ixyzhm_leafparts) &
 !$omp shared(nnode_toreplace,h_solvertocmi,nixyzhmf_cminode) &
 !$omp private(inode)
 do inode = 1,maxcminode
    nxyzm_treetocmi(1:5,inode)  = (/ 0.,0.,0.,0.,0./)
    h_solvertocmi(inode)        = 0.
    nnode_toreplace(inode)      = 0
    nixyzhmf_cminode(1:8,inode) = (/ 0.,0.,0.,0.,0.,0.,0.,0. /)
 enddo
 !$omp end parallel do

 !$omp parallel do default(none) shared(ixyzhm_leafparts) &
 !$omp private(ip)
 do ip = 1,maxleafparts
    ixyzhm_leafparts(1:6,ip) = (/ 0.,0.,0.,0.,0.,0. /)
 enddo
 !$omp end parallel do

end subroutine reset_cminode


subroutine allocate_cmi_inputs_history
 use allocutils,   only:allocate_array

 call allocate_array('x_old',x_old,maxcminode+1)
 call allocate_array('y_old',y_old,maxcminode+1)
 call allocate_array('z_old',z_old,maxcminode+1)

end subroutine allocate_cmi_inputs_history


subroutine reset_cmi_inputs_history
 integer :: inode

 !$omp parallel do default(none) shared(x_old,y_old,z_old) &
 !$omp private(inode)
 do inode = 1,maxcminode
   x_old(inode) = 0.
   y_old(inode) = 0.
   z_old(inode) = 0.
 enddo
 !$omp end parallel do

end subroutine reset_cmi_inputs_history


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

 if (allocated(x)) deallocate(x)
 if (allocated(y)) deallocate(y)
 if (allocated(z)) deallocate(z)
 if (allocated(h)) deallocate(h)
 if (allocated(m)) deallocate(m)
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
 call write_inopt(sink_ionsrc,'sink_ionsrc','Using sinks as ionizing sources',iunit)
 call write_inopt(masscrit_ionize,'masscrit_ionize','Critical sink mass to begin emitting radiation',iunit)
 call write_inopt(niter_mcrt,'niter_mcrt','number of photon-release iterations',iunit)
 call write_inopt(nphoton,'nphoton','number of photons per iteration',iunit)
 call write_inopt(wavelength,'wavelength','wavelength of photons in nm',iunit)
 call write_inopt(temp_hii,'temp_hii','Temperature of ionized gas',iunit)
 call write_inopt(tol_vsite,'tol_vsite','Threshold to update Voronoi gen-sites',iunit)
 call write_inopt(lloyd,'lloyd','Apply Lloyd iteration to construct Voronoi grid',iunit)
 call write_inopt(photoionize_tree,'photoionize_tree','Construct Voronoi grid around tree nodes',iunit)
 call write_inopt(tree_accuracy_cmi,'tree_accuracy_cmi','Tree-opening criteria for cmi-nodes',iunit)
 call write_inopt(rcut_opennode,'rcut_opennode','Radius within which nodes must be leaves',iunit)
 call write_inopt(rcut_leafpart,'rcut_leafpart','Radius within which particles are extracted',iunit)
 call write_inopt(auto_opennode,'auto_opennode','Automatically adjust tree-walk with iterative approach',iunit)
 call write_inopt(auto_tree_acc,'auto_tree_acc','Automatically adjust tree_accuracy_cmi',iunit)
 call write_inopt(delta_rcut,'delta_rcut','Increase in rcut_opennode and rcut_leafpart per iter step',iunit)
 call write_inopt(nHlimit_fac,'nHlimit_fac','Paramter controlling resolution of ionization front',iunit)
 call write_inopt(min_nodesize_toflag,'min_nodesize_toflag','Minimum node size to check nH',iunit)

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
 case('sink_ionsrc')
    read(valstring,*,iostat=ierr) sink_ionsrc
    ngot = ngot + 1
 case('masscrit_ionize')
    read(valstring,*,iostat=ierr) masscrit_ionize
    ngot = ngot + 1
    if (masscrit_ionize <= 0.) call fatal(label,'invalid setting for masscrit_ionize (<=0)')
 case('niter_mcrt')
    read(valstring,*,iostat=ierr) niter_mcrt
    ngot = ngot + 1
    if (niter_mcrt <= 0.) call fatal(label,'invalid setting for niter_mcrt (<=0)')
 case('nphoton')
    read(valstring,*,iostat=ierr) nphoton
    ngot = ngot + 1
    if (nphoton <= 0.) call fatal(label,'invalid setting for nphoton (<=0)')
case('wavelength')
    read(valstring,*,iostat=ierr) wavelength
    ngot = ngot + 1
    if (wavelength <= 0.) call fatal(label,'invalid setting for wavelength (<=0)')
 case('temp_hii')
    read(valstring,*,iostat=ierr) temp_hii
    ngot = ngot + 1
    if (temp_hii <= 0.) call fatal(label,'invalid setting for temp_hii (<=0)')
 case('tol_vsite')
    read(valstring,*,iostat=ierr) tol_vsite
    ngot = ngot + 1
    if (tol_vsite < 0.) call fatal(label,'invalid setting for tol_vsite (<0)')
 case('lloyd')
    read(valstring,*,iostat=ierr) lloyd
    ngot = ngot + 1
 case('photoionize_tree')
    read(valstring,*,iostat=ierr) photoionize_tree
    ngot = ngot + 1
 case('tree_accuracy_cmi')
    read(valstring,*,iostat=ierr) tree_accuracy_cmi
    ngot = ngot + 1
    if (tree_accuracy_cmi < 0. .or. tree_accuracy_cmi > 1.) call fatal(label,'invalid setting for tree_accuracy_cmi')
 case('rcut_opennode')
    read(valstring,*,iostat=ierr) rcut_opennode
    ngot = ngot + 1
    if (rcut_opennode < 0.) call fatal(label,'invalid setting for rcut_opennode')
 case('rcut_leafpart')
    read(valstring,*,iostat=ierr) rcut_leafpart
    if (rcut_leafpart < 0.) call fatal(label,'invalid setting for rcut_leafpart')
    ngot = ngot + 1
case('auto_opennode')
    read(valstring,*,iostat=ierr) auto_opennode
    ngot = ngot + 1
case('auto_tree_acc')
    read(valstring,*,iostat=ierr) auto_tree_acc
    ngot = ngot + 1
case('delta_rcut')
    read(valstring,*,iostat=ierr) delta_rcut
    ngot = ngot + 1
    if (delta_rcut < 0.) call fatal(label,'invalid setting for delta_rcut')
case('nHlimit_fac')
    read(valstring,*,iostat=ierr) nHlimit_fac
    ngot = ngot + 1
    if (nHlimit_fac < 0.) call fatal(label,'invalid setting for nHlimit_fac')
case('min_nodesize_toflag')
    read(valstring,*,iostat=ierr) min_nodesize_toflag
    ngot = ngot + 1
    if (min_nodesize_toflag < 0.) call fatal(label,'invalid setting for min_nodesize_toflag')
 case default
    imatch = .false.
 end select
 igotall = ( ngot >= 17 )

end subroutine read_options_photoionize

end module photoionize_cmi
