!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module testptmass
!
! Unit tests of the ptmass/sink particles module
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, checksetup, deriv, dim, energies, eos,
!   fileutils, io, kdtree, kernel, mpiutils, options, part, physcon,
!   ptmass, random, setbinary, setdisc, spherical, step_lf_global,
!   stretchmap, testutils, timestep, units
!
 implicit none
 public :: test_ptmass

 private

contains

subroutine test_ptmass(ntests,npass)
 use dim,      only:maxp,mhd,periodic,gravity,maxptmass
 use io,       only:id,master,iverbose,iskfile
 use part,     only:init_part,npart,npartoftype,massoftype,xyzh,hfact,vxyzu,fxyzu,fext,&
                    xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,nptmass,epot_sinksink,&
                    ihacc,isdead_or_accreted,igas,divcurlv,iphase,isetphase,maxphase,&
                    Bevol,fxyzu,set_particle_type,ispinx,ispiny,ispinz,poten,rhoh
 use eos,             only:gamma,polyk
 use timestep,        only:dtmax,C_force,tolv
 use testutils,       only:checkval,checkvalf,update_test_scores
 use setbinary,       only:set_binary
 use step_lf_global,  only:step,init_step
 use energies,        only:compute_energies,etot,totmom,epot,angtot !,accretedmass
 use ptmass,          only:get_accel_sink_sink,ptmass_accrete,h_soft_sinksink, &
                           ptmass_create,h_acc,get_accel_sink_gas,f_acc,finish_ptmass, &
                           ipart_rhomax,icreate_sinks, &
                           idxmsi,idymsi,idzmsi,idmsi,idspinxsi,idspinysi,idspinzsi, &
                           idvxmsi,idvymsi,idvzmsi,idfxmsi,idfymsi,idfzmsi, &
                           ndptmass,update_ptmass, &
                           r_merge_uncond,r_merge_cond,r_merge_uncond2,r_merge_cond2,r_merge2
 use physcon,         only:pi
 use setdisc,         only:set_disc
 use spherical,       only:set_sphere
 use boundary,        only:set_boundary
 use deriv,           only:get_derivs_global
 use kdtree,          only:tree_accuracy
 !use readwrite_dumps, only:write_fulldump,write_smalldump
 use fileutils,       only:getnextfilename
 use checksetup,      only:check_setup
 use options,         only:ieos,iexternalforce
 use units,           only:set_units
 use kernel,          only:kernel_softening
#ifdef IND_TIMESTEPS
 use part,            only:ibin,ibin_old
#endif
 use mpiutils,        only:bcast_mpi,reduce_in_place_mpi,reduceloc_mpi
 use stretchmap,      only:rho_func
 use random,          only:ran2
 integer, intent(inout) :: ntests,npass
 integer                :: i,j,nsteps,nbinary_tests,itest,nerr,nwarn,itestp
 integer                :: nparttot,merge_ij(maxptmass),merge_n,nsink0,nsinkF
 logical                :: test_binary,test_accretion,test_createsink,test_softening,test_merger
 logical                :: accreted,merged,merged_expected
 real                   :: m1,m2,a,ecc,hacc1,hacc2,dt,dtext,t,dtnew,dr,v2
 real                   :: etotin,totmomin,dtsinksink,omega,mred,errmax,angmomin
 real                   :: r2,r2min,xcofm(3),totmass,dum,dum2,psep,tolen,dtmp(3)
 real                   :: xyzm_ptmass_old(4,1), vxyz_ptmass_old(3,1)
 real                   :: q,phisoft,fsoft,mu,v_c1,v_c2,r1,omega1,omega2,mv0,angmom0,mtot0,mvF,angmomF,mtotF
 real                   :: dptmass(ndptmass,maxptmass)
 real                   :: dptmass_thread(ndptmass,maxptmass)
 real                   :: fxyz_sinksink(4,maxptmass),rhomax_test,rhomax
 integer                :: norbits,itmp,ierr,iseed
 integer                :: nfailed(56),imin(1)
 integer                :: id_rhomax,ipart_rhomax_global
 integer(kind=1)        :: ibin_wakei
 character(len=20)      :: dumpfile,filename
 procedure(rho_func), pointer :: density_func
 logical                :: print_sink_paths = .false. ! print sink paths in the merger test

 if (id==master) write(*,"(/,a,/)") '--> TESTING PTMASS MODULE'

 density_func => gaussianr

 test_binary = .true.
 test_accretion = .true.
 test_createsink = .true.
 test_softening = .true.
 test_merger = .true.
 nbinary_tests = 3
 !
 !--general settings
 !
 polyk = 0.
 gamma = 1.
 iexternalforce = 0
 call init_part()
!
!  Test 1: orbit of a single binary
!  Test 2: with gas disc around it
!
 testbinary: if (test_binary) then
    !
    !--no gas particles
    !
    npart = 0
    npartoftype(:) = 0
    massoftype = 0.

    xyzh(:,:)  = 0.
    vxyzu(:,:) = 0.
    fxyzu(:,:) = 0.
    fext(:,:)  = 0.
    Bevol(:,:) = 0.

#ifdef IND_TIMESTEPS
    ibin(:) = 0_1
    ibin_old(:) = 0_1
#endif
    iverbose = 0
    tree_accuracy = 0.
    h_soft_sinksink = 0.

    binary_tests: do itest = 1,nbinary_tests
       select case(itest)
       case(2,3)
          if (periodic) then
             if (id==master) write(*,"(/,a)") '--> skipping circumbinary disc test (-DPERIODIC is set)'
             cycle binary_tests
          else
             if (itest==3) then
                if (id==master) write(*,"(/,a)") '--> testing integration of disc around eccentric binary'
             else
                if (id==master) write(*,"(/,a)") '--> testing integration of circumbinary disc'
             endif
          endif
       case default
          if (id==master) write(*,"(/,a)") '--> testing integration of binary orbit'
       end select
       !
       !--setup sink-sink binary (no gas particles)
       !
!       time = 0.
       nptmass = 0
       m1    = 1.
       m2    = 1.
       a     = 1.
       if (itest==3) then
          ecc = 0.5
       else
          ecc   = 0.
       endif
       hacc1  = 0.35
       hacc2  = 0.35
       C_force = 0.25
       t = 0.
       call set_units(mass=1.d0,dist=1.d0,G=1.d0)
       call set_binary(m1,m2,a,ecc,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,verbose=.false.)
       if (ierr /= 0) nerr = nerr + 1
       if (itest==2 .or. itest==3) then
          !  add a circumbinary gas disc around it
          nparttot = 1000
          call set_disc(id,master,nparttot=nparttot,npart=npart,rmin=1.5*a,rmax=15.*a,p_index=1.5,q_index=0.75,&
                        HoverR=0.1,disc_mass=0.01*m1,star_mass=m1+m2,gamma=gamma,&
                        particle_mass=massoftype(igas),hfact=hfact,xyzh=xyzh,vxyzu=vxyzu,&
                        polyk=polyk,verbose=.false.)
          npartoftype(1) = npart

          !
          ! check that no errors occurred when setting up disc
          !
          nfailed = 0
          call check_setup(nerr,nwarn)
          call checkval(nerr,0,0,nfailed(1),'no errors during disc setup')
          !call checkval(nwarn,0,0,nfailed(2),'no warnings during disc setup')
          call update_test_scores(ntests,nfailed,npass)
       endif

       tolv = 1.e3
       iverbose = 0
       ieos = 3
       !print*,'initial v = ',vxyz_ptmass(:,1)
       !
       ! initialise forces
       !
       if (id==master) then
          call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_sinksink,epot_sinksink,dtsinksink,0,0.,merge_ij,merge_n)
       endif
       fxyz_ptmass(:,:) = 0.
       call bcast_mpi(epot_sinksink)
       call bcast_mpi(dtsinksink)

       fext(:,:) = 0.
       do i=1,npart
          call get_accel_sink_gas(nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),xyzmh_ptmass,&
                   fext(1,i),fext(2,i),fext(3,i),dum,massoftype(igas),fxyz_ptmass,dum,dum2)
       enddo
       if (id==master) fxyz_ptmass(:,:) = fxyz_ptmass(:,:) + fxyz_sinksink(:,:)
       call reduce_in_place_mpi('+',fxyz_ptmass)
       !
       !--take the sink-sink timestep specified by the get_forces routine
       !
       if (id==master) print*,' dt for sinks = ',C_force*dtsinksink
       dt      = C_force*dtsinksink !2.0/(nsteps)
       dtmax   = dt  ! required prior to derivs call, as used to set ibin
       !
       !--compute SPH forces
       !
       if (itest==2 .or. itest==3) then
          fxyzu(:,:) = 0.
          call get_derivs_global()
       endif
       !
       !--evolve this for a number of orbits
       !
       call compute_energies(t)
       etotin   = etot
       totmomin = totmom
       angmomin = angtot
       !
       !--check that initial potential on the two sinks is correct
       !
       nfailed(:) = 0
       if (itest==1) then
          call checkval(epot_sinksink,-m1*m2/a,epsilon(0.),nfailed(1),'potential energy')
          call update_test_scores(ntests,nfailed,npass)
          !
          !--check initial angular momentum on the two sinks is correct
          !
          call checkval(angtot,m1*m2*sqrt(a/(m1 + m2)),1.e6*epsilon(0.),nfailed(1),'angular momentum')
          call update_test_scores(ntests,nfailed,npass)
       endif
       !
       !--determine number of steps per orbit for information
       !
       mred    = m1*m2/(m1 + m2)
       omega   = sqrt(mred/a**3)
       nsteps  = int(2.*pi/omega/dt) + 1
       if (itest==2 .or. itest==3) then
          norbits = 10
       else
          norbits = 100
       endif
       if (id==master) print*,' nsteps per orbit = ',nsteps,' norbits = ',norbits
       nsteps = nsteps*norbits
       errmax = 0.
       dumpfile='test_00000'
       f_acc = 1.
       call init_step(npart,t,dtmax)
       do i=1,nsteps
          t = t + dt
          dtext = dt
          if (id==master .and. iverbose > 2) write(*,*) ' t = ',t,' dt = ',dt
          call step(npart,npart,t,dt,dtext,dtnew)
          call compute_energies(t)
          errmax = max(errmax,abs(etot - etotin))
          if (itest==2) then
             !   write(1,*) i,t,angtot,totmom,etot,accretedmass,epot
             !   print*,i,t,angtot,totmom,etot,accretedmass,epot
          endif
       enddo
       call compute_energies(t)
       nfailed(:) = 0
       select case(itest)
       case(3)
#ifdef IND_TIMESTEPS
          call checkval(angtot,angmomin,2.1e-6,nfailed(3),'angular momentum')
          call checkval(totmom,totmomin,5.e-6,nfailed(2),'linear momentum')
#else
          call checkval(angtot,angmomin,1.2e-6,nfailed(3),'angular momentum')
          call checkval(totmom,totmomin,4.e-14,nfailed(2),'linear momentum')
#endif
          call checkval(etotin+errmax,etotin,1.2e-2,nfailed(1),'total energy')
       case(2)
          call checkval(angtot,angmomin,3.e-7,nfailed(3),'angular momentum')
          call checkval(totmom,totmomin,6.e-14,nfailed(2),'linear momentum')
          tolen = 2.e-3
          if (gravity) tolen = 3.1e-3
          call checkval(etotin+errmax,etotin,tolen,nfailed(1),'total energy')
       case default
          call checkval(angtot,angmomin,3.e-14,nfailed(3),'angular momentum')
          call checkval(totmom,totmomin,epsilon(0.),nfailed(2),'linear momentum')
          call checkval(etotin+errmax,etotin,3.e-8,nfailed(1),'total energy')
       end select
       !
       !--check energy conservation
       !
       do i=1,3
          call update_test_scores(ntests,nfailed(i:i),npass)
       enddo
    enddo binary_tests

 endif testbinary

!
!  Test of softening between sinks
!
 testsoftening: if (test_softening) then
    if (id==master) write(*,"(/,a)") '--> testing softening in sink particle binary'
    nptmass = 0
    npart = 0
    npartoftype = 0
    m1    = 1.
    m2    = 1.
    a     = 1.
    ecc   = 0.
    hacc1 = 0.
    hacc2 = 0.
    t     = 0.

    h_soft_sinksink = 0.8*a

    call set_units(mass=1.d0,dist=1.d0,G=1.d0)
    call set_binary(m1,m2,a,ecc,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,verbose=.false.)

    q   = a/h_soft_sinksink
    call kernel_softening(q*q,q,phisoft,fsoft)

    ! Test energy and momentum conservation
    mu = (m1*m2)/(m1+m2)
    r1 = sqrt(xyzmh_ptmass(1,1)**2+xyzmh_ptmass(2,1)**2+xyzmh_ptmass(3,1)**2)
    r2 = sqrt(xyzmh_ptmass(1,2)**2+xyzmh_ptmass(2,2)**2+xyzmh_ptmass(3,2)**2)
    omega1 = sqrt(m2*fsoft/(r1*h_soft_sinksink**2))
    omega2 = sqrt(m1*fsoft/(r2*h_soft_sinksink**2))
    v_c1 = omega1*r1
    v_c2 = omega2*r2
    vxyz_ptmass(1,1) = 0.
    vxyz_ptmass(2,1) = v_c1
    vxyz_ptmass(3,1) = 0.
    vxyz_ptmass(1,2) = 0.
    vxyz_ptmass(2,2) = -v_c2
    vxyz_ptmass(3,2) = 0.
    call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,epot_sinksink,dtsinksink,0,0.,merge_ij,merge_n)
    call compute_energies(t)
    etotin   = etot
    totmomin = totmom
    angmomin = angtot

    call checkval(epot,m1*m2*(phisoft)/h_soft_sinksink,2.*epsilon(0.),nfailed(1),'potential energy')
    call update_test_scores(ntests,nfailed(1:1),npass)

    C_force = 0.25
    dt      = 0.3*C_force*dtsinksink
    !if (id==master) print*,' dt for sinks = ',dt
    dtmax   = dt
    omega   = omega1
    nsteps  = int(2.*pi/omega/dt) + 1
    norbits = 10
    if (id==master) print*,' nsteps per orbit = ',nsteps,' norbits = ',norbits
    nsteps = nsteps*norbits
    errmax = 0.
    iverbose = 0
    do i=1,nsteps
       t = t + dt
       dtext = dt
       if (id==master .and. iverbose > 2) write(*,*) ' t = ',t,' dt = ',dt
       call step(npart,npart,t,dt,dtext,dtnew)
!      write(1,'(10(es10.3,1x))') xyzmh_ptmass(:,1)
!      write(1,'(10(es10.3,1x))') xyzmh_ptmass(:,2)
       call compute_energies(t)
       errmax = max(errmax,abs(etot - etotin))
    enddo
    call compute_energies(t)
    nfailed(:) = 0
    call checkval(angtot,angmomin,2.e-14,nfailed(1),'angular momentum')
    call checkval(totmom,totmomin,tiny(0.),nfailed(2),'linear momentum')
    call checkval(etotin+errmax,etotin,2.e-9,nfailed(3),'total energy')
!    call checkval(      ,r_max,1.e-10,nfailed(4),'radius')

    call update_test_scores(ntests,nfailed(1:4),npass)

    ! Reset sink softening
    h_soft_sinksink = 0.0

 endif testsoftening

!
!  Tests of accrete_particle routine
!
 xyzmh_ptmass(:,:) = 0.
 vxyz_ptmass(:,:)  = 0.

 testaccretion: if (test_accretion) then
    if (id==master) write(*,"(/,a)") '--> testing accretion onto sink particles'
    nptmass = 1
    !--setup 1 point mass at (-5,-5,-5)
    xyzmh_ptmass(1:3,1)   = 1.
    xyzmh_ptmass(4,1)     = 10. ! mass of sink
    xyzmh_ptmass(ihacc,1) = 20. ! accretion radius
    vxyz_ptmass(1:3,1)    = -40.
    fxyz_ptmass(1:3,1)    = 40.
    massoftype(1)   = 10.
    !--setup 1 SPH particle at (5,5,5)
    if (id==master) then
       call set_particle_type(1,igas)
       npartoftype(igas) = 1
       npart        = 1
       xyzh(1:3,1)  = 5.
       xyzh(4,1)    = 0.01
       vxyzu(1:3,1) = 80.
       fxyzu(1:3,1) = 20.
    else
       npartoftype(igas) = 0
       npart        = 0
    endif
    xyzm_ptmass_old = xyzmh_ptmass(1:4,1:nptmass)
    vxyz_ptmass_old = vxyz_ptmass (1:3,1:nptmass)
    dr = sqrt(dot_product(xyzh(1:3,1) - xyzmh_ptmass(1:3,1),xyzh(1:3,1) - xyzmh_ptmass(1:3,1)))
    !--perform a test of the accretion of the SPH particle by the point mass
    nfailed(:)  = 0
    !--check energies before accretion event
    t=0.
    call compute_energies(t)
    etotin   = etot
    totmomin = totmom
    angmomin = angtot
    ibin_wakei = 0

    dptmass(:,1:nptmass) = 0.
    !$omp parallel default(shared) private(i) firstprivate(dptmass_thread)
    dptmass_thread(:,1:nptmass) = 0.
    !$omp do
    do i=1,npart
       call ptmass_accrete(1,nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),&
                           vxyzu(1,i),vxyzu(2,i),vxyzu(3,i),fxyzu(1,i),fxyzu(2,i),fxyzu(3,i), &
                           igas,massoftype(igas),xyzmh_ptmass,vxyz_ptmass, &
                           accreted,dptmass_thread,t,1.0,ibin_wakei,ibin_wakei)
    enddo
    !$omp enddo
    !$omp critical(dptmassadd)
    dptmass(:,1:nptmass) = dptmass(:,1:nptmass) + dptmass_thread(:,1:nptmass)
    !$omp end critical(dptmassadd)
    !$omp end parallel

    call reduce_in_place_mpi('+',dptmass(:,1:nptmass))

    if (id==master) call update_ptmass(dptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,nptmass)

    call bcast_mpi(xyzmh_ptmass(:,1:nptmass))
    call bcast_mpi(vxyz_ptmass(:,1:nptmass))
    call bcast_mpi(fxyz_ptmass(:,1:nptmass))

    if (id==master) then
       call checkval(accreted,.true.,nfailed(1),'accretion flag')
       !--check that h has been changed to indicate particle has been accreted
       call checkval(isdead_or_accreted(xyzh(4,1)),.true.,nfailed(2),'isdead_or_accreted flag')
    endif
    call checkval(xyzmh_ptmass(1,1),3.,tiny(0.),nfailed(3),'x(ptmass) after accretion')
    call checkval(xyzmh_ptmass(2,1),3.,tiny(0.),nfailed(4),'y(ptmass) after accretion')
    call checkval(xyzmh_ptmass(3,1),3.,tiny(0.),nfailed(5),'z(ptmass) after accretion')
    call checkval(vxyz_ptmass(1,1),20.,tiny(0.),nfailed(6),'vx(ptmass) after accretion')
    call checkval(vxyz_ptmass(2,1),20.,tiny(0.),nfailed(7),'vy(ptmass) after accretion')
    call checkval(vxyz_ptmass(3,1),20.,tiny(0.),nfailed(8),'vz(ptmass) after accretion')
    call checkval(fxyz_ptmass(1,1),30.,tiny(0.),nfailed(9), 'fx(ptmass) after accretion')
    call checkval(fxyz_ptmass(2,1),30.,tiny(0.),nfailed(10),'fy(ptmass) after accretion')
    call checkval(fxyz_ptmass(3,1),30.,tiny(0.),nfailed(11),'fz(ptmass) after accretion')

    call update_test_scores(ntests,nfailed(1:2),npass)
    call update_test_scores(ntests,nfailed(3:5),npass)
    call update_test_scores(ntests,nfailed(6:8),npass)
    call update_test_scores(ntests,nfailed(9:11),npass)

    !--compute energies after accretion event
    nfailed(:) = 0
    call compute_energies(t)
    call checkval(angtot,angmomin,1.e-10,nfailed(3),'angular momentum')
    call checkval(totmom,totmomin,epsilon(0.),nfailed(2),'linear momentum')
    !call checkval(etot,etotin,1.e-6,'total energy',nfailed(1))
    call update_test_scores(ntests,nfailed(3:3),npass)
    call update_test_scores(ntests,nfailed(2:2),npass)

 endif testaccretion
!
!  Test sink particle creation
!
 testcreatesink: if (test_createsink) then
    do itest=1,2
       select case(itest)
       case(2)
          if (id==master) write(*,"(/,a)") '--> testing sink particle creation (sin)'
       case default
          if (id==master) write(*,"(/,a)") '--> testing sink particle creation (uniform density)'
       end select
       xyzmh_ptmass(:,:) = 0.
       vxyz_ptmass(:,:) = 0.
       nptmass = 0

       npart = 0
       npartoftype(:) = 0
       massoftype = 0.

       xyzh(:,:)  = 0.
       vxyzu(:,:) = 0.
       fxyzu(:,:) = 0.
       fext(:,:)  = 0.
       Bevol(:,:) = 0.
       iverbose   = 1
       call set_boundary(-1.,1.,-1.,1.,-1.,1.)
       !
       ! set up gas particles in a uniform sphere with radius R=0.2
       !
       psep = 0.05  ! required as a variable since this may change under conditions not requested here
       if (itest==2) then
          ! use random so particle with maximum density is unique
          call set_sphere('cubic',id,master,0.,0.2,psep,hfact,npartoftype(igas),xyzh,rhofunc=density_func)
       else
          call set_sphere('cubic',id,master,0.,0.2,psep,hfact,npartoftype(igas),xyzh)
       endif
       totmass = 1.0
       massoftype(igas) = totmass/real(npartoftype(igas))
       npart = npartoftype(igas)
       !
       ! give inward radial velocities
       !
       itestp = npart
       r2min = huge(r2min)
       xcofm(:) = 0.
       do i=1,npart
          r2 = dot_product(xyzh(1:3,i),xyzh(1:3,i))
          if (r2 < r2min) then
             itestp = i
             r2min = r2
          endif
          xcofm = xcofm + xyzh(1:3,i)
       enddo

       xcofm = massoftype(igas)*xcofm/totmass
       if (maxphase==maxp) iphase(1:npart) = isetphase(igas,iactive=.true.)
       !
       ! set up tree for neighbour finding
       ! and make sure that gravitational potential energy has been computed
       !
       tree_accuracy = 0.
       icreate_sinks = 1
       iverbose = 1
       call get_derivs_global()
       !
       ! check that particle being tested is at the maximum density
       !
       if (itest==2 .and. gravity) then
          imin = minloc(xyzh(4,1:npart))
          itestp = imin(1)
          rhomax_test = rhoh(xyzh(4,itestp),massoftype(igas))
          !
          ! only check on the thread that has rhomax
          !
          ipart_rhomax_global = ipart_rhomax
          call reduceloc_mpi('max',ipart_rhomax_global,id_rhomax)
          if (id == id_rhomax) then
             rhomax = rhoh(xyzh(4,ipart_rhomax),massoftype(igas))
             call checkval(rhomax,rhomax_test,epsilon(0.),nfailed(1),'rhomax')
          else
             itestp = -1 ! set itest = -1 on other threads
             call checkval(ipart_rhomax,-1,0,nfailed(1),'ipart_rhomax')
          endif
          call update_test_scores(ntests,nfailed(1:1),npass)
       endif
       !
       ! check energies before insertion of sink
       !
       call compute_energies(t)
       !print*,' got ETOT = ',etot,totmom,epot
       etotin   = etot
       totmomin = totmom
       angmomin = angtot
       !
       ! now create point mass by accreting these particles
       !
       h_acc = 0.15

       !
       ! if gravity is not enabled, then need to choose a particle to create ptmass from
       !
       if (.not. gravity) then
          ipart_rhomax_global = itestp
          call reduceloc_mpi('max',ipart_rhomax_global,id_rhomax)
       endif
       call ptmass_create(nptmass,npart,itestp,xyzh,vxyzu,fxyzu,fext,divcurlv,poten,&
                          massoftype,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,0.)
       !
       ! check that creation succeeded
       !
       nfailed(:) = 0
       call checkval(nptmass,1,0,nfailed(1),'nptmass=1')
       call update_test_scores(ntests,nfailed,npass)
       !
       ! check that linear and angular momentum and energy is conserved
       !
       nfailed(:) = 0
       call compute_energies(t)
       !print*,' got ETOT = ',etot,totmom,epot
       call checkval(angtot,angmomin,1.e-10,nfailed(3),'angular momentum')
       call checkval(totmom,totmomin,epsilon(0.),nfailed(2),'linear momentum')
       !call checkval(etot,etotin,1.e-6,nfailed(1),'total energy')
       call update_test_scores(ntests,nfailed,npass)

       iverbose = 0
       call finish_ptmass(nptmass)
    enddo
 endif testcreatesink
!
!  Test sink particle creation
!
 testsinkmerger: if (test_merger) then
    iseed           = -74205
    nfailed(:)      = 0
    iverbose        = 0
    nptmass         = 2
    npart           = 0
    h_acc           = 0.1
    h_soft_sinksink =    h_acc
    r_merge_uncond  = 2.*h_acc    ! sinks will unconditionally merge if they touch
    r_merge_cond    = 4.*h_acc    ! sinks will merge if bound within this radius
    r_merge_uncond2 = r_merge_uncond**2
    r_merge_cond2   = r_merge_cond**2
    r_merge2        = max(r_merge_uncond2,r_merge_cond2)
    do itest=1,8
       t                 = 0.
       xyzmh_ptmass(:,:) = 0.
       xyzmh_ptmass(4,:) = 1.
       xyzmh_ptmass(ihacc,:) = h_acc
       vxyz_ptmass(:,:)  = 0.
       select case(itest)
       case(1)
          if (id==master) write(*,"(/,a)") '--> testing fast flyby: no merger'
          ! fast flyby within r_merge_uncond < r < r_merge_cond
          xyzmh_ptmass(1,1) =  1.
          xyzmh_ptmass(2,1) =  1.5*h_acc
          vxyz_ptmass(1,1)  = -10.
          merged_expected   = .false.
       case(2)
          if (id==master) write(*,"(/,a)") '--> testing fast flyby: impact so merger'
          ! fast flyby within r < r_merge_uncond
          xyzmh_ptmass(1,1) =  1.
          xyzmh_ptmass(2,1) =  0.5*h_acc
          vxyz_ptmass(1,1)  = -10.
          merged_expected   = .true.
       case(3)
          if (id==master) write(*,"(/,a)") '--> testing slow flyby: capture and merger'
          ! slow flyby within r_merge_uncond < r < r_merge_cond
          xyzmh_ptmass(1,1) =  1.
          xyzmh_ptmass(2,1) =  1.5*h_acc
          vxyz_ptmass(1,1)  = -1.
          merged_expected   = .true.
       case(4)
          if (id==master) write(*,"(/,a)") '--> testing slow flyby: impact and merger'
          ! slow flyby within r < r_merge_cond
          xyzmh_ptmass(1,1) =  1.
          xyzmh_ptmass(2,1) =  0.5*h_acc
          vxyz_ptmass(1,1)  = -1.
          merged_expected   = .true.
       case(5)
          if (id==master) write(*,"(/,a)") '--> testing flyby: slingshot & no merger'
          ! flyby within r_merge_uncond < r < r_merge_cond
          xyzmh_ptmass(1,1) =  1.
          xyzmh_ptmass(2,1) =  1.5*h_acc
          vxyz_ptmass(1,1)  = -5.
          merged_expected   = .false.
       case(6)
          if (id==master) write(*,"(/,a)") '--> testing orbit: stable & no merger'
          ! stable orbit within r >  r_merge_cond
          xyzmh_ptmass(1,1) = 2.5*h_acc
          vxyz_ptmass(2,1)  = sqrt(0.25*xyzmh_ptmass(4,1)/xyzmh_ptmass(1,1))
          merged_expected   = .false.
       case(7)
          if (id==master) write(*,"(/,a)") '--> testing orbit: decaying & merger'
          ! decaying orbit within r >  r_merge_cond
          xyzmh_ptmass(1,1) = 2.5*h_acc
          vxyz_ptmass(2,1)  = 0.9*sqrt(0.25*xyzmh_ptmass(4,1)/xyzmh_ptmass(1,1))
          merged_expected   = .true.
       case(8)
          if (id==master) write(*,"(/,a)") '--> testing multiple sink interations'
          nptmass = 100
          do i = 1,nptmass
             xyzmh_ptmass(1:3,i) = ( (/ran2(iseed),ran2(iseed),ran2(iseed)/) - 0.5) * 2.  ! in range (-1,1)
             vxyz_ptmass(1:3,i)  = ( (/ran2(iseed),ran2(iseed),ran2(iseed)/) - 0.5) * 6.  ! in range (-3,3)
          enddo
          merged_expected   = .true. ! this logical does not have meaning here
       end select
       if (itest /= 8) then
          xyzmh_ptmass(1:3,2) = -xyzmh_ptmass(1:3,1)
          vxyz_ptmass(1:3,2)  = -vxyz_ptmass(1:3,1)
       endif
       !
       ! get initial values
       call get_momenta(nptmass,nsink0,xyzmh_ptmass,vxyz_ptmass,angmom0,mv0,mtot0,periodic)
       !
       ! initialise forces
       !
       if (id==master) then
          call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_sinksink,epot_sinksink,dtsinksink,0,0.,merge_ij,merge_n)
       endif
       fxyz_ptmass(:,:) = 0.
       call bcast_mpi(epot_sinksink)
       call bcast_mpi(dtsinksink)

       if (id==master) fxyz_ptmass(:,:) = fxyz_ptmass(:,:) + fxyz_sinksink(:,:)
       call reduce_in_place_mpi('+',fxyz_ptmass)
       !
       ! integrate
       !
       nsteps = 1000
       dt     = huge(dt)
       do i = 1,nptmass-1
          do j = i+1,nptmass
             dtmp = xyzmh_ptmass(1:3,i)-xyzmh_ptmass(1:3,j)
             r2   = dot_product(dtmp,dtmp)
             dtmp = vxyz_ptmass(1:3,i)-vxyz_ptmass(1:3,j)
             v2   = dot_product(dtmp,dtmp)+epsilon(v2)
             dt   = min(dt,sqrt(r2/v2))
          enddo
       enddo
       dt     = 2.*dt/nsteps
       dtmax  = dt*nsteps
       t      = 0.
       print*, itest,dt
       call init_step(npart,t,dtmax)
       if (print_sink_paths) then
          write(333,*) itest,0,xyzmh_ptmass(1:4,1),xyzmh_ptmass(1:4,2)
          if (itest==8) then
             do j = 1,nptmass
                if (xyzmh_ptmass(4,j) > 0.) write(334,*) t,j,xyzmh_ptmass(1:4,j)
             enddo
          endif
       endif
       do i=1,nsteps
          t = t + dt
          dtext = dt
          if (id==master .and. iverbose > 2) write(*,*) ' t = ',t,' dt = ',dt
          call step(npart,npart,t,dt,dtext,dtnew)
          if (print_sink_paths) then
             write(333,*) itest,i,xyzmh_ptmass(1:4,1),xyzmh_ptmass(1:4,2)
             if (itest==8) then
                do j = 1,nptmass
                   if (xyzmh_ptmass(4,j) > 0.) write(334,*) t,j,xyzmh_ptmass(1:4,j)
                enddo
             endif
          endif
       enddo
       !
       ! check results
       call get_momenta(nptmass,nsinkF,xyzmh_ptmass,vxyz_ptmass,angmomF,mvF,mtotF,periodic)
       !
       if (xyzmh_ptmass(4,2) < 0.) then
          merged = .true.
       else
          merged = .false.
       endif
       if (itest==8) then
          call checkval(nsinkF,41,0,nfailed(itest),'final number of sinks')
       else
          call checkval(merged,merged_expected,nfailed(itest),'merger')
          if (merged_expected) then
             call checkval(xyzmh_ptmass(1,1),0.,epsilon(0.),nfailed(2*itest),'final x-position')
             call checkval(xyzmh_ptmass(2,1),0.,epsilon(0.),nfailed(3*itest),'final y-position')
             v2 = dot_product(vxyz_ptmass(1:2,1),vxyz_ptmass(1:2,1))
             call checkval(sqrt(v2),0.,epsilon(0.),nfailed(4*itest),'final velocity')
          endif
       endif
       call checkval(mvF,    mv0,    1.e-13,nfailed(5*itest),'conservation of linear momentum')
       call checkval(angmomF,angmom0,1.e-13,nfailed(6*itest),'conservation of angular momentum')
       call checkval(mtotF,  mtot0,  1.e-13,nfailed(7*itest),'conservation of mass')
    enddo
    call update_test_scores(ntests,nfailed(1:56),npass)

 endif testsinkmerger

 !--reset stuff, turn off sink creation & clean up temporary files
 itmp    = 201
 nptmass = 0
 icreate_sinks  = 0
 r_merge_uncond = 0.
 r_merge_cond   = 0.
 close(iskfile,iostat=ierr)
 write(filename,"(i3)") iskfile
 filename = 'fort.'//trim(adjustl(filename))
 open(unit=itmp,file=filename,status='old',iostat=ierr)
 close(itmp,status='delete',iostat=ierr)
 open(unit=itmp,file='Sink00.sink',status='old',iostat=ierr)
 close(itmp,status='delete',iostat=ierr)

 if (id==master) write(*,"(/,a)") '<-- PTMASS TEST COMPLETE'

end subroutine test_ptmass

subroutine get_momenta(nptmass,nsnk,xyzmh_ptmass,vxyz_ptmass,ang,mv,mtot,periodic)
 use part, only:ispinx,ispiny,ispinz
 integer, intent(in)  :: nptmass
 integer, intent(out) :: nsnk
 real,    intent(in)  :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real,    intent(out) :: ang,mv,mtot
 logical, intent(in)  :: periodic
 integer              :: i
 real                 :: angx,angy,angz,xmom,ymom,zmom

 angx = 0.
 angy = 0.
 angz = 0.
 xmom = 0.
 ymom = 0.
 zmom = 0.
 mv   = 0.
 mtot = 0.
 nsnk = 0
 do i = 1,nptmass
    if (xyzmh_ptmass(4,i) > 0.) then
       nsnk = nsnk + 1
       angx = angx + xyzmh_ptmass(4,i)*(xyzmh_ptmass(2,i)*vxyz_ptmass(3,i) - xyzmh_ptmass(3,i)*vxyz_ptmass(2,i))
       angy = angy + xyzmh_ptmass(4,i)*(xyzmh_ptmass(3,i)*vxyz_ptmass(1,i) - xyzmh_ptmass(1,i)*vxyz_ptmass(3,i))
       angz = angz + xyzmh_ptmass(4,i)*(xyzmh_ptmass(1,i)*vxyz_ptmass(2,i) - xyzmh_ptmass(2,i)*vxyz_ptmass(1,i))
       angx = angx + xyzmh_ptmass(ispinx,i)
       angy = angy + xyzmh_ptmass(ispiny,i)
       angz = angz + xyzmh_ptmass(ispinz,i)
       xmom = xmom + xyzmh_ptmass(4,i)*vxyz_ptmass(1,i)
       ymom = ymom + xyzmh_ptmass(4,i)*vxyz_ptmass(2,i)
       zmom = zmom + xyzmh_ptmass(4,i)*vxyz_ptmass(3,i)
       mtot = mtot + xyzmh_ptmass(4,i)
    endif
 enddo
 mv  = sqrt(xmom*xmom + ymom*ymom + zmom*zmom)
 ang = sqrt(angx*angx + angy*angy + angz*angz)
 if (periodic) ang = 0. ! since this will not be conserved with periodic boundaries

end subroutine get_momenta

real function gaussianr(r)
 real, intent(in) :: r

 gaussianr = exp(-(r/0.05)**2) !1./(r**2 + 0.0001**2)

end function gaussianr

end module testptmass
