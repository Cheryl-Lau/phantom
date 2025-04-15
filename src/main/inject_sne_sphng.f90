!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module inject
!
! Routine for injecting Type II supernovae
!
! :References: Lucas,et.al,2020,MNRAS,493,4700
!              Maeder,2009,Springer
!
! :Owner: Cheryl Lau (based on origional SPHNG addsupernova.f by Ian Bonnell)
!
! :Runtime parameters:
!   - maxsn           : *maximum number of supernovae at a time*
!   - r_sn_cgs        : *supernova blast radius*
!   - engsn_cgs       : *total energy [erg] released from suepernova*
!   - frackin         : *fraction of kinetic energy in eng_sn*
!   - fractherm       : *fraction of thermal energy in eng_sn*
!   - pmsncrit_cgs    : *critical mass of sink particles*
!   - sink_progenitor : *if true, inits sne when sink exceeds critical mass;
!                        otherwise inits sn(e) with user-given time and location*
!
! :Dependencies: eos, infile_utils, io, part, partinject, physcon, dim, datafiles,
!                     timestep, cooling, units, random, ptmass
!
 use dim,  only:maxptmass
 implicit none
 character(len=*), parameter, public :: inject_type = 'supernovae'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject

 logical, public :: inject_sn = .true.  ! switch on/off sn for testing

 ! Runtime parameters for supernovae injection to read from input file
 integer, public :: maxsn     = 100
 real,    public :: r_sn_cgs  = 3.086E17    ! 0.1 pc
 real,    public :: engsn_cgs = 1E51        ! 1E51 erg
 real,    public :: frackin   = 1.0
 real,    public :: fractherm = 0.0

 real,    public :: pmsncrit_cgs = 1.59E34  ! 8 solarm
 logical, public :: sink_progenitor = .false.
 logical, public :: one_sink_progenitor = .false.
 integer, public :: isink_progenitor = 5

 logical, public :: gaussian_vprofile = .true.     ! Set up Gasussian SN particle radial velocity profile 
 logical, public :: uniform_vprofile = .false.     ! Set up uniform velocity profile 

 logical, public :: delay_sn_injection  = .true.
 logical, public :: delay_by_mslifetime = .false.   ! delay SN inject by MS lifetime, else 10 dtmax

 private

 ! Set sne properties if not using sinks as progenitors
 integer, parameter :: maxsn_insert = 1
 real    :: mstar_cgs = 1.59E34  ! 8 solarm
 real    :: xyzt_sn_insert_cgs(4,maxsn_insert) = reshape((/ 0., 0., 0., 0. /), &
                                                           shape=(/4,maxsn_insert/))

 ! Global storage for all sne (also used for switching-off cooling)
 integer, parameter :: maxallsn = 500
 integer :: nallsn,iprog_allsn(maxallsn),iinsert_allsn(maxallsn)
 real    :: xyzht_allsn(5,maxallsn),vxyz_allsn(3,maxallsn),m_allsn(maxallsn)

 logical :: snflag_time(maxsn_insert) ! Flag when sn is injected
 logical :: snflag_sink(maxptmass)    ! Flag when sink particle has boomed
 logical :: queueflag(maxallsn)

 integer :: maxvprofile = 1E6
 integer :: maxnpartsn  = 3000
 real,   allocatable :: vprofile(:,:) ! Tabulated velocity profiles of sn particles
 logical :: first_sn

 integer :: iseed = -432
 real    :: r_sn,engsn,pmsncrit,mstar
 real    :: xyzt_sn_insert(4,maxsn_insert)

contains
!-----------------------------------------------------------------------
!+
!  Initialize global variables or arrays needed for injection routine
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 use datafiles,  only:find_phantom_datafile
 use units,      only:umass,udist,utime,unit_energ
 use physcon,    only:solarm
 use io,         only:warning
 integer, intent(out) :: ierr
 integer  :: i,ip,isn,ientry
 !
 ! return without error
 !
 ierr = 0
 !
 ! Init flags
 !
 if (sink_progenitor) then
    do ip = 1,maxptmass
       snflag_sink(ip) = .false.
    enddo
 else
    do isn = 1,maxsn_insert
       snflag_time(isn) = .false.
    enddo
 endif
 !
 ! Calculate mass of progenitor
 !
 if (.not.sink_progenitor) then
    if (mstar_cgs < pmsncrit_cgs) call warning('inject_sne_sphng','recommend setting a larger mstar_cgs')
    mstar = mstar_cgs/umass
    print*,'mstar_cgs mstar',mstar_cgs,mstar
 endif
 !
 ! Import tabulated velocity profile (Precomputed solutions to
 ! Ek = sum_{i}{Nsn} 1/2 mi vi(a)^2 with vi(a) defined by vrfunc)
 !
 open(1723,file=find_phantom_datafile('vprofile_sn.txt','sne_sphng'))
 allocate(vprofile(4,maxvprofile))
 do ientry = 1,maxvprofile
   read(1723,*) (vprofile(i,ientry), i=1,4)  ! a,Nsn,m,Ek
 enddo
 print*,'Done loading vprofile'
 close(1723)
 !
 ! Extract a_vrad from velocity profile table only for the first sn
 !
 first_sn = .true.
 !
 ! Init counters/table for perm-storage of all sne taken place in simulation
 !
 nallsn = 0
 xyzht_allsn = 0.
 vxyz_allsn  = 0.
 if (sink_progenitor) then
    m_allsn = 0.
    iprog_allsn = 0.
 else
    iinsert_allsn = 0.
 endif
 queueflag = .false.
 !
 ! Convert input params to code units
 !
 r_sn  = r_sn_cgs/udist
 engsn = engsn_cgs/unit_energ
 pmsncrit = pmsncrit_cgs/umass
 if (.not.sink_progenitor) then
    xyzt_sn_insert(1:3,:) = xyzt_sn_insert_cgs(1:3,:)/udist
    xyzt_sn_insert(4,:)   = xyzt_sn_insert_cgs(4,:)  /utime
 endif

 !
 ! Remove threshold requirement if specifying sink 
 !
 if (sink_progenitor .and. one_sink_progenitor) then 
     print*,' Setting sink',isink_progenitor,'to detonate; ignore mass threshold'
     pmsncrit = tiny(pmsncrit)  
 endif 

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine handling supernovae injection
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npartoftype,dtinject)
 use io,         only:master,fatal,warning
 use part,       only:nptmass,massoftype,igas,kill_particle,hfact
 use partinject, only:add_or_update_particle
 use physcon,    only:pi,pc 
 use units,      only:unit_energ,udist,unit_velocity 
 use timestep,   only:dtmax
 use cooling,    only:snecoolingoff,maxcoolingoff,xyzh_coolingoff,dt_cooloff,ncoolingoff
 real,    intent(in)    :: time,dtlast
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 integer  :: ipart,isn,iadd,ia,ip,ii,icf,iallsn,ncand
 integer  :: npartold,nsn,nkilled,npartsn,ipartsn,progindex,iprog_sn(maxsn),iinsert_sn(maxsn)
 integer  :: i_rightcand,iEk
 integer  :: indexallsn(maxsn)
 integer, allocatable :: i_cand(:)
 real,    allocatable :: Ek_cand(:) 
 real     :: xyz_sn(3,maxsn),vxyz_sn(3,maxsn),m_sn(maxsn),h_sn(maxsn),mradsn,engkin,engtherm,numden,t_sn
 real     :: xyz_partsn(3),vxyz_partsn(3),xyz_ref(3,6),r_ref,hnew,unew,dist,dist_nearby
 real     :: vrad,a_vrad,r_vrad,aN_scalefactor
 real     :: ekintot,ethermtot,etot_allparts_old,etot_allparts,endiff_cgs,ekin_allparts_old,etherm_allparts_old
 real     :: vrad_min,vrad_max,r_ref_max,r_vrad_max
 logical  :: timematch_inject,queueflag_iallsn

 !
 ! Init/Reset variables of sne properties injected at current timestep (temp-storage)
 !
 nsn = 0              ! number of sne going off now
 if (sink_progenitor) then
    iprog_sn = 0      ! sink index of progenitors
    m_sn  = 0.        ! mass of progenitors
 else
    iinsert_sn = 0.   ! index of user-defined sne
 endif
 h_sn  = 0.           ! smoothing length of progenitors
 xyz_sn  = 0.
 vxyz_sn = 0.
 !
 ! Detect and put all sne into a dynamically-created queue (perm-storage)
 !
 if (sink_progenitor) then
    if (nptmass > 0) call check_sink(xyzmh_ptmass,vxyz_ptmass,nptmass,pmsncrit,time)
 else
    call check_time(time,npart,xyzh)
 endif
 !
 ! Loop through perm-storage and move those to inject now into temp-storage
 !
 if (nallsn >= 1) then
    do iallsn = 1,nallsn
       timematch_inject = abs(xyzht_allsn(5,iallsn) - time) < dtmax+tiny(dtmax)
       queueflag_iallsn = queueflag(iallsn)
       if (timematch_inject .and. .not.queueflag_iallsn) then
          nsn = nsn + 1
          if (nsn > maxsn) call fatal('inject_sne_sphng','number of sn exceeded maxsn')
          xyz_sn(1:3,nsn)  = xyzht_allsn(1:3,iallsn)
          vxyz_sn(1:3,nsn) = vxyz_allsn(1:3,iallsn)
          h_sn(nsn) = xyzht_allsn(4,iallsn)
          if (sink_progenitor) then
             m_sn(nsn) = m_allsn(iallsn)
             iprog_sn(nsn) = iprog_allsn(iallsn)
          else
             iinsert_sn(nsn) = iinsert_allsn(iallsn)
          endif
          indexallsn(nsn) = iallsn  !- store the iallsn index associated with isn
       endif
    enddo
 endif

 !
 ! Inject supernovae from temp-storage
 !
 if (inject_sn .and. nsn >= 1) then
    print*,' *****************'
    write(*,'(1x,a10,i3,a17)') 'Injecting ',nsn,' supernova(e) at:'
    print*,' *****************'
    do isn = 1,nsn
       print*, xyz_sn(1:3,isn)
    enddo

    !- Checking
    ekintot = 0.
    ethermtot = 0.
    do ipart = 1,npart
       ekintot   = ekintot + 1./2.*massoftype(igas)*dot_product(vxyzu(1:3,ipart),vxyzu(1:3,ipart))
       ethermtot = ethermtot + vxyzu(4,ipart)*massoftype(igas)
    enddo
    ekin_allparts_old = ekintot
    etherm_allparts_old = ethermtot
    etot_allparts_old =  ekintot + ethermtot

    over_sne: do isn = 1,nsn
       !
       ! Kill all particles within blast radius
       !
       mradsn  = 0.
       nkilled = 0
       do ipart = 1,npart
          dist = mag(xyzh(1:3,ipart) - xyz_sn(1:3,isn))
          if (dist < r_sn) then
             call kill_particle(ipart,npartoftype)
             nkilled = nkilled + 1
             mradsn  = mradsn + massoftype(igas)
          endif
       enddo
       if (nkilled == 0) then
          if (sink_progenitor) then
             call warning('inject_sne_sphng','no particles found within blast radius')
          else
             call warning('inject_sne_sphng','no particles found within blast radius')
          endif
       endif

       !
       ! Adding in 1/4 mass of progenitor into supernova
       !
       if (sink_progenitor) then
          mradsn = mradsn + 0.25*m_sn(isn)
          progindex = iprog_sn(isn)
!          xyzmh_ptmass(4,progindex) = 0.75*xyzmh_ptmass(4,progindex)
       else
          mradsn = mradsn + 0.25*mstar
       endif

       !
       ! Energies of supernova
       !
       engkin   = engsn * frackin
       engtherm = engsn * fractherm

       !
       ! Set number of new particles assuming mass is conserved
       !
       npartsn = nint(mradsn/massoftype(igas))
       print*,'Number of SN particles: ',npartsn
       if (npartsn > nkilled) print*,'adding ',npartsn-nkilled,' extra particles'
       if (npartsn < 6) call fatal('inject_sne_sphng','too few sn particles')
       if (npartsn > maxnpartsn .and. engkin > 0.d0) then 
          call warning('inject_sne_sphng','number of sn particles is beyond pre-computed table range')
          npartsn = maxnpartsn
          print*,'Setting npartsn to ',npartsn 
       endif 

       !
       ! Internal energy of sn particles
       !
       unew = (engtherm/npartsn)/massoftype(igas)

       !
       ! Set up velocity profile of sn particles
       !
       if (engkin > 0.d0) then 
          if (gaussian_vprofile) then 
             aN_scalefactor = 1. ! to avoid compiler warning
             if (first_sn) then
                call interp_from_vrprofile(npartsn,massoftype(igas),engkin,a_vrad) 
                aN_scalefactor = a_vrad*npartsn**(1./2.)  ! Store scale factor for forthcoming sne
                first_sn = .false.
             else
                a_vrad = aN_scalefactor/npartsn**(1./2.)
             endif
          elseif (uniform_vprofile) then 
             vrad = sqrt((2*engkin)/(npartsn*massoftype(igas)))
          endif 
       else 
          a_vrad = 0.d0
          vrad = 0.d0 
       endif 

       !
       ! Estimate new smoothing length
       !
       numden = npartsn / (4./3.*pi*r_sn**3)
       hnew = hfact * (numden)**(-1./3.)

       !
       ! Add npartsn particles back-in in a symmetrical manner -
       !    Generates a random set of 6 refpoints at same distance and right-angles to
       !    each other around the origin;
       !    Uses the 6 refpoints to place 6 particles as a group at a time.
       !
       vrad_max = tiny(vrad_max)
       r_ref_max = tiny(r_ref_max)
       open(unit=2040,file='sn_velocity_profile.dat',status='replace')
       write(2040,'(2a20)') 'r [cm]','vel [cm/s]'

       npartold = npart
       over_partsngroup: do iadd = npartold+1,npartold+npartsn,6
          ipartsn = iadd
          call gen_refpoints(xyz_ref,r_ref)
          if (gaussian_vprofile) then 
             r_vrad = 2.d0*r_ref   ! 2*r_dist/r_sn; scaling max=r_sn to max=2 of vrad func
             vrad = vrfunc(a_vrad,r_vrad) 
          endif 
          over_partsn: do ia = 1,6
             xyz_partsn  = xyz_sn(1:3,isn) + xyz_ref(1:3,ia)*r_sn
             vxyz_partsn = vxyz_sn(1:3,isn) + vrad*xyz_ref(1:3,ia)/r_ref
             call add_or_update_particle(igas,xyz_partsn,vxyz_partsn,hnew,unew,ipartsn,npart,&
                                         npartoftype,xyzh,vxyzu)
             ipartsn = ipartsn + 1
          enddo over_partsn
          ! Record max vel in mid-regions and vel on the outermost 
          if (vrad > vrad_max) then 
             vrad_max = vrad 
             r_vrad_max = r_ref 
          endif 
          if (r_ref > r_ref_max) then 
             r_ref_max = r_ref
             vrad_min  = vrad 
          endif 
          write(2040,'(2e20.10)') r_ref*r_sn*udist, vrad*unit_velocity 
       enddo over_partsngroup

       print*,'vrad_min [cm/s] ',vrad_min*unit_velocity,'at radius [pc]',r_ref_max*r_sn*udist/pc 
       print*,'vrad_max [cm/s] ',vrad_max*unit_velocity,'at radius [pc]',r_vrad_max*r_sn*udist/pc 
       print*,'approx distance to form shell [pc] ',(r_ref_max*r_sn-r_vrad_max*r_sn)*vrad_max/(vrad_max-vrad_min) *udist/pc 
       close(2040)

 !      !
 !      ! Put nearby particles around r_sn to smallest timestep to avoid being 'shocked'
 !      !
 !      over_nearbyparts: do ip = 1,npart
 !         dist_nearby = mag(xyzh(1:3,in) - xyz_sn(1:3,isn))
 !         if (dist_nearby > r_sn .and. dist_nearby <= r_sn+2.0) then
 !            call add_or_update_particle(igas,xyzh(1:3,in),vxyzu(1:3,in),xyzh(4,in),vxyzu(4,in),in,&
 !                                        npart,npartoftype,xyzh,vxyzu)
 !         endif
 !      enddo over_nearbyparts

       !
       ! Flag to avoid further detonation
       !
       if (sink_progenitor) then
          ip = iprog_sn(isn)
          snflag_sink(ip) = .true.
       else
          ii = iinsert_sn(isn)
          snflag_time(ii) = .true.
       endif
       ! to avoid picking it up again from perm storage queue
       iallsn = indexallsn(isn)
       queueflag(iallsn) = .true.

    enddo over_sne

    !- Checking
    ekintot = 0.
    ethermtot = 0.
    do ipart = 1,npart
       ekintot  = ekintot + 1./2.*massoftype(igas)*dot_product(vxyzu(1:3,ipart),vxyzu(1:3,ipart))
       ethermtot = ethermtot + vxyzu(4,ipart)*massoftype(igas)
    enddo
    etot_allparts = ekintot + ethermtot

    print*,'Added KE [erg]: ',(ekintot - ekin_allparts_old) *unit_energ 
    print*,'Added TE [erg]: ',(ethermtot - etherm_allparts_old) *unit_energ

    endiff_cgs = (etot_allparts - etot_allparts_old) *unit_energ
    print*,'Energy added [erg]: ',endiff_cgs

    !- check if they are in the same order of magnitude
    if (abs(log10(engsn_cgs*nsn) - log10(endiff_cgs)) > 0.5) then
       call fatal('inject_sne_sphng','Amount of energy injected does not match the required value')
    endif

 endif
 !
 ! Store coordinates of all sne which, at the current time, needs to have
 ! its nearby particles' cooling switched off
 !
 if (snecoolingoff) then
    ! Refresh table
    ncoolingoff = 0
    xyzh_coolingoff = 0.
    ! Store sn properties into xyzh_coolingoff table (to be then passed to cooling)
    over_allsne: do icf = 1,maxallsn
       t_sn = xyzht_allsn(5,icf)
       if (t_sn > 0. .and. time >= t_sn .and. time < (t_sn+dt_cooloff)) then
          ncoolingoff = ncoolingoff + 1
          xyzh_coolingoff(1:4,ncoolingoff) = xyzht_allsn(1:4,icf)
       endif
    enddo over_allsne
 endif
 !
 ! Timestep constraint
 !
 dtinject = huge(dtinject)

end subroutine inject_particles

! ----------------------------------------------------------------------------------
!
! Subroutines for dynamically initiating supernovae to put into the queue (perm-storage)
! during runtime
!
! ----------------------------------------------------------------------------------
subroutine check_sink(xyzmh_ptmass,vxyz_ptmass,nptmass,pmsncrit,time)
 use units,    only:utime,umass
 use physcon,  only:solarm
 use io,       only:fatal
 use timestep, only:dtmax
 integer, intent(in) :: nptmass
 real,    intent(in) :: pmsncrit,time
 real,    intent(in) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer  :: ip
 real     :: mptmass,mptmass_solarm,t_ms_cgs,t_sn
 logical  :: snflag_sinkip

 over_sinks: do ip = 1,nptmass

    if (one_sink_progenitor .and. ip /= isink_progenitor) cycle over_sinks 

    mptmass = xyzmh_ptmass(4,ip)
    snflag_sinkip = snflag_sink(ip)
    if (mptmass >= pmsncrit .and. .not.snflag_sinkip) then
       !
       ! Time to inject sn - estimated using star MS lifetime (obtained by fitting Table 25.6, P.628, Maeder 2009)
       !
       if (delay_sn_injection) then 
          if (delay_by_mslifetime) then 
             mptmass_solarm = mptmass*umass/solarm
             t_ms_cgs = (1.165E8/(mptmass_solarm-4.242) + 1.143E6) *365*24*3600
          else 
             t_ms_cgs = (10.*dtmax)*utime   ! any number of dtmax
          endif 
          write(*,'(1x,a5,i4,a38,f5.2,a4)') 'Sink ',ip,' set as progenitor - detonating after ',t_ms_cgs/(1E6*365*24*3600),' Myr'
          t_sn = time + t_ms_cgs/utime
       else 
          write(*,'(1x,a5,i4,a35)') 'Sink ',ip,' set as progenitor - detonating now'
          t_sn = time
       endif 
       !
       ! Store sn into queue
       !
       nallsn = nallsn + 1
       if (nallsn > maxallsn) call fatal('inject_sne_sphng','number of sn exceeded maxallsn')
       xyzht_allsn(1:3,nallsn) = xyzmh_ptmass(1:3,ip)
       xyzht_allsn(4,nallsn)   = xyzmh_ptmass(5,ip)
       xyzht_allsn(5,nallsn)   = t_sn
       vxyz_allsn(1:3,nallsn)  = vxyz_ptmass(1:3,ip)
       m_allsn(nallsn)     = mptmass
       iprog_allsn(nallsn) = ip
       !
       ! Flag sink to avoid being placed into queue again
       !
       snflag_sink(ip) = .true.
    endif
 enddo over_sinks

end subroutine check_sink


subroutine check_time(time,npart,xyzh)
 use timestep, only:dtmax
 use io,       only:fatal
 integer, intent(in)  :: npart
 real,    intent(in)  :: time,xyzh(:,:)
 real,    allocatable :: sep(:)
 integer  :: isn,i,iclosest
 logical  :: timematch,snflag_time_isn

 do isn = 1,maxsn_insert
    timematch = abs(xyzt_sn_insert(4,isn) - time) < dtmax
    snflag_time_isn = snflag_time(isn)

    if (timematch .and. .not.snflag_time_isn) then
       print*, 'timematch ok'
       !
       ! Estimate h of progenitor using h of the current closest particle (??)
       !
       allocate (sep(npart))
       do i = 1,npart
          sep(i) = mag(xyzh(1:3,i) - xyzt_sn_insert(1:3,isn))
       enddo
       iclosest = minloc(sep(:),1)
       deallocate (sep)
       !
       ! Store sn into queue
       !
       nallsn = nallsn + 1
       if (nallsn > maxallsn) call fatal('inject_sne_sphng','number of sn exceeded maxallsn')
       xyzht_allsn(1:3,nallsn) = xyzt_sn_insert(1:3,isn)
       xyzht_allsn(4,nallsn)   = xyzh(4,iclosest)
       xyzht_allsn(5,nallsn)   = time
       vxyz_allsn(1:3,nallsn)  = (/ 0.,0.,0. /)  ! assuming source is stationary
       iinsert_allsn(nallsn)   = isn             ! to flag the sn
    endif
 enddo

 ! Note: Do not directly transfer the user-given sn table into the perm-storage
 !       since we need to extract the closest h at the point of injection.

end subroutine check_time

! ----------------------------------------------------------------------------------
!
! Subroutine for generating the 6 reference points
!
! ----------------------------------------------------------------------------------
subroutine gen_refpoints(xyz_ref,r_ref)
 use random, only:ran2
 use ptmass, only:h_acc
 real, intent(out) :: xyz_ref(3,6),r_ref
 real  :: x1,y1,z1,x2,y2,z2,dot_ab
 real  :: a_vec(3),b_vec(3),c_vec(3)

 ! Random vector point a within sphere of unit radius
 ! and beyond a small accretion radius
 a_vec = (/ 1.,1.,1. /)
 do while (mag(a_vec) > 1.d0 .or. mag(a_vec) < h_acc/r_sn)   
    x1 = 2.*(ran2(iseed)-0.5)
    y1 = 2.*(ran2(iseed)-0.5)
    z1 = 2.*(ran2(iseed)-0.5)
    a_vec = (/ x1,y1,z1 /)
 enddo
 r_ref = mag(a_vec)
 ! Random vector point b at right-angles to a with mag of r_ref
 x2 = 2.*(ran2(iseed)-0.5)
 y2 = 2.*(ran2(iseed)-0.5)
 z2 = 2.*(ran2(iseed)-0.5)
 b_vec = (/ x2,y2,z2 /)
 dot_ab = dot_product(a_vec,b_vec)
 b_vec = b_vec - dot_ab*a_vec(1:3)/(mag(a_vec))**2
 b_vec = r_ref * b_vec(1:3)/mag(b_vec)
 ! Vector point c at right-angles to a & b with mag of r_ref
 c_vec = cross_product(a_vec,b_vec)
 c_vec = r_ref * c_vec(1:3)/mag(c_vec)
 ! Generate the 6 points
 xyz_ref(:,1) = a_vec(:)
 xyz_ref(:,2) = b_vec(:)
 xyz_ref(:,3) = c_vec(:)
 xyz_ref(:,4) = -1*a_vec(:)
 xyz_ref(:,5) = -1*b_vec(:)
 xyz_ref(:,6) = -1*c_vec(:)

end subroutine gen_refpoints

! ----------------------------------------------------------------------------------
!
! Subroutine for generating velocity profile of sn particles
!
! ----------------------------------------------------------------------------------
real function vrfunc(a,r)
 real, intent(in) :: a,r

 if (r >= 0 .and. r < 1) then
    vrfunc = a * (3.*r-9./4.*r**2)
 elseif (r >= 1 .and. r < 2) then
    vrfunc = a * 3./4.*(2-r)**2
 else
    vrfunc = 0.
 endif

end function vrfunc

!
! Extract indicies of all entries in vprofile table (and its corresponding Ek)
! where both npartsn and m are closest to the required value 
!  vprofile(4,maxvprofile)  a,Nsn,m,Ek
!
subroutine interp_from_vrprofile(npartsn,m,engkin,a_vrad) 
 use io,    only:fatal
 use units, only:umass,udist,unit_energ,unit_velocity 
 integer, intent(in)  :: npartsn
 real,    intent(in)  :: m,engkin
 real,    intent(out) :: a_vrad
 integer :: i,ifill,npartsn_in,npartsn_closest,mindiff_npartsn,ncand
 real    :: m_in,mindiff_m,m_closest,engkin_extracted,engkin_in,mindiff_engkin,a_vrad_cgs
 real    :: tol_m = 1E+29

 ! Use the max/min value if parameters are beyond vprofile range
 m_in = m
 npartsn_in = npartsn
 if (m > 1E-1) then
    print*,'Warning - particle mass is beyond vprofile range. Setting mi to 1E-1.'
    m_in = 1E-1
 elseif (m < 1E-3) then
    print*,'Warning - particle mass is smaller than vprofile range. Setting mi to 1E-3.'
    m_in = 1E-3
 endif
 if (npartsn > maxnpartsn) then
    print*,'Warning - number of supernova particles is beyond vprofile range. Setting npartsn to ',maxnpartsn
    npartsn_in = maxnpartsn
 endif

 ! Convert units 
 m_in = m_in*umass 
 engkin_in = engkin*unit_energ 


 ! Find value of closest npart_in, closest m_in, and closest engkin_in
 mindiff_m = huge(mindiff_m)
 mindiff_npartsn = huge(mindiff_npartsn)
 mindiff_engkin = huge(mindiff_engkin)
 do i = 1,maxvprofile
    if (abs(m_in - vprofile(3,i)) < mindiff_m) then 
       mindiff_m = abs(m_in - vprofile(3,i)) 
       m_closest = vprofile(3,i)
    endif 
    if (abs(npartsn_in - vprofile(2,i)) < mindiff_npartsn) then 
       mindiff_npartsn = abs(npartsn_in - vprofile(2,i))
       npartsn_closest = vprofile(2,i) 
    endif 
 enddo 

 print*,'m_in,m_closest,mindiff_m',m_in,m_closest,mindiff_m
 print*,'npartsn_in,npartsn_closest,mindiff_npartsn',npartsn_in,npartsn_closest,mindiff_npartsn

 ! Extract rows whoes m and npartsn match npartsn_closest and m_closest
 mindiff_engkin = huge(mindiff_engkin)
 do i = 1,maxvprofile
    if (vprofile(2,i) == npartsn_closest .and.  &
        vprofile(3,i) > m_closest-tol_m  .and.  &
        vprofile(3,i) < m_closest+tol_m) then 

       if (abs(engkin_in - vprofile(4,i)) < mindiff_engkin) then 
          mindiff_engkin = abs(engkin_in - vprofile(4,i)) 
          engkin_extracted = vprofile(4,i)
          a_vrad_cgs = vprofile(1,i)
       endif 
    endif 
 enddo 

 a_vrad = a_vrad_cgs/unit_velocity
 print*,'engkin_extracted,a_vrad',engkin_extracted,a_vrad

end subroutine interp_from_vrprofile


! ----------------------------------------------------------------------------------
!
! Vector math functions
!
! ----------------------------------------------------------------------------------
real function mag(vec)
 real, intent(in) :: vec(3)

 mag = sqrt(dot_product(vec,vec))

end function mag

function cross_product(vec1,vec2)
 real, intent(in)   :: vec1(3),vec2(3)
 real, dimension(3) :: cross_product

 cross_product(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
 cross_product(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
 cross_product(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)

end function cross_product

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use physcon,      only: au, solarm, years
 use infile_utils, only: write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options for injecting supernova'
 call write_inopt(inject_sn,'inject_sn','inject SNe',iunit)
 call write_inopt(sink_progenitor,'sink_progenitor','init sne with sinks',iunit)
 call write_inopt(one_sink_progenitor,'one_sink_progenitor','only set one sink to be progenitor',iunit)
 call write_inopt(isink_progenitor,'isink_progenitor','sink label to detonate',iunit)
 call write_inopt(engsn_cgs,'engsn_cgs','total energy released from sn',iunit)
 call write_inopt(frackin,'frackin','fraction of kinetic energy',iunit)
 call write_inopt(fractherm,'fractherm','fraction of thermal energy',iunit)
 call write_inopt(r_sn_cgs,'r_sn_cgs','blast radius of supernova',iunit)
 call write_inopt(maxsn,'maxsn','maximum number of supernovae at a time',iunit)
 call write_inopt(pmsncrit_cgs,'pmsncrit_cgs','critical mass of sinks',iunit)
 call write_inopt(delay_sn_injection,'delay_sn_injection','delay injection',iunit)
 call write_inopt(delay_by_mslifetime,'delay_by_mslifetime','delay injection by MS lifetime',iunit)

end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io,  only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical, intent(out) :: imatch,igotall
 integer,intent(out) :: ierr
 integer, save :: ngot = 0
 integer       :: noptions
 character(len=30), parameter :: label = 'read_options_inject'

 imatch  = .true.
 igotall = .false.

 select case(trim(name))
 case('inject_sn')
    read(valstring,*,iostat=ierr) inject_sn
    ngot = ngot + 1
 case('sink_progenitor')
    read(valstring,*,iostat=ierr) sink_progenitor
    ngot = ngot + 1
 case('one_sink_progenitor')
    read(valstring,*,iostat=ierr) one_sink_progenitor
    ngot = ngot + 1
 case('isink_progenitor')
    read(valstring,*,iostat=ierr) isink_progenitor
    ngot = ngot + 1
 case('engsn_cgs')
    read(valstring,*,iostat=ierr) engsn_cgs
    ngot = ngot + 1
    if (engsn_cgs < 0.) call fatal(label,'invalid setting for engsn_cgs (<0)')
 case('frackin')
    read(valstring,*,iostat=ierr) frackin
    ngot = ngot + 1
    if (frackin < 0.) call fatal(label,'invalid setting for frackin (<0)')
 case('fractherm')
    read(valstring,*,iostat=ierr) fractherm
    ngot = ngot + 1
    if (fractherm < 0.) call fatal(label,'invalid setting for fractherm (<0)')
 case('r_sn_cgs')
    read(valstring,*,iostat=ierr) r_sn_cgs
    ngot = ngot + 1
    if (r_sn_cgs < 0.) call fatal(label,'invalid setting for r_sn_cgs (<0)')
 case('maxsn')
    read(valstring,*,iostat=ierr) maxsn
    ngot = ngot + 1
    if (maxsn < 1.) call fatal(label,'invalid setting for maxsn (<1)')
 case('pmsncrit_cgs')
    read(valstring,*,iostat=ierr) pmsncrit_cgs
    ngot = ngot + 1
    if (pmsncrit_cgs < 0.) call fatal(label,'invalid setting for pmsncrit_cgs (<0)')
 case('delay_sn_injection')
    read(valstring,*,iostat=ierr) delay_sn_injection
    ngot = ngot + 1
 case('delay_by_mslifetime')
    read(valstring,*,iostat=ierr) delay_by_mslifetime
    ngot = ngot + 1

 case default
    imatch = .false.
 end select

 noptions = 12
 igotall  = (ngot >= noptions)

end subroutine read_options_inject

end module inject
