!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module inject
!
! Routine for faking an HII region
!
! :References: Stromgren B. 1939, ApJ, 89, 526
!
! :Owner: Cheryl Lau
!
! :Runtime parameters:
!   - lumin_LyC : *Rate of Lyman continuum photon emission from star*
!
! :Dependencies: eos, infile_utils, io, part, partinject, physcon, dim, datafiles,
!                     timestep, cooling, units
!
! Note: Only works with code units of pc;M_sun
!
 implicit none
 character(len=*), parameter, public :: inject_type = 'HII_region'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject

 ! Runtime parameters
 real, public :: lumin_LyC = 1E39  ! s^-1
 real, public :: temp_hii  = 1E4   ! K

 private

 ! Setting location and time
 real    :: xyzt_hii(4) = (/ -3.,0.,0.,1E-4 /)
 logical :: hii_flag

contains
!-----------------------------------------------------------------------
!+
!  Initialize global variables or arrays needed for injection routine
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 integer, intent(out) :: ierr
 !
 ! Flag when HII region has been injected
 !
 hii_flag = .false.
 !
 ! return without error
 !
 ierr = 0

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine handling HII region injection
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npartoftype,dtinject)
 use io,         only:id,master,fatal
 use part,       only:nptmass,massoftype,iphase,igas,kill_particle,hfact,ibin
 use partinject, only:add_or_update_particle
 use physcon,    only:mass_proton_cgs,kboltz
 use eos,        only:gamma,gmw
 use timestep,   only:dtmax
 use units,      only:unit_density,unit_ergg
 use cooling,    only:ipart_nocooling,maxpart_nocooling,npartnocool
 real,    intent(in)    :: time,dtlast
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 integer  :: npart_nearby,ip,ientry
 real     :: totdens,totu,dist,rhoi,ui,dens_bg,u_bg,rho_cgs,R_st,dR_st,u_hii,unew
 logical  :: timematch
 logical  :: inject_hii = .true.  ! switch on/off HII region for testing
 !
 ! Check time
 !
 timematch = abs(xyzt_hii(4) - time) < dtmax
 !
 ! Inject HII region by heating particles within pre-defined Stromgren sphere
 !
 if (timematch .and. .not.hii_flag) then
    !
    ! Get current average density and temperature around this region
    !
    totdens = 0.
    totu    = 0.
    npart_nearby = 0
    do ip = 1,npart
       dist = mag(xyzh(1:3,ip) - xyzt_hii(1:3))
       if (dist < 10.) then
          rhoi = massoftype(igas) * (xyzh(4,ip)/hfact)**3
          ui   = vxyzu(4,ip)
          totdens = totdens + rhoi
          totu    = totu + ui
          npart_nearby = npart_nearby + 1
       endif
    enddo
    if (npart_nearby < 10) call fatal('inject_test_hiiregion','not enough particles to estimate &
                                      Stromgren sphere radius')
    dens_bg = totdens/npart_nearby
    u_bg    = totu/npart_nearby
    !
    ! Compute Stromgren sphere and transition radius
    !
    rho_cgs = dens_bg*unit_density ! g cm^-3
    R_st  = 0.3 * (lumin_LyC/1E49)**(1./3.) * (rho_cgs/1E-20)**(-2./3.)  ! pc
    dR_st = 0.0002 * (rho_cgs/1E-20)**(-1)  ! pc
    print*,'Stromgren sphere',R_st
    open(2023,file='Stromgren sphere.txt')
    write(2023,*) R_st
    write(2023,*) dR_st
    !
    ! Set up storage of particle indicies to stop cooling
    !
    npartnocool = 0
    do ientry = 1,maxpart_nocooling
       ipart_nocooling(ientry) = 0
    enddo
    !
    ! Heat particles within Stromgren sphere, and smooth temp gradient across ionization front
    !
    u_hii = kboltz*temp_hii/(gmw*mass_proton_cgs*(gamma-1.0))/unit_ergg  ! T -> u (code units)
    over_parts: do ip = 1,npart
       dist = mag(xyzh(1:3,ip) - xyzt_hii(1:3))
       if (dist < (R_st+dR_st)) then
          if (dist < R_st) then
             unew = u_hii
          else
             unew = u_hii + (dist-R_st)*(u_bg-u_hii)/dR_st  ! Interpolate to smooth transition
          endif
          vxyzu(4,ip) = unew
          ! Store particle indicies to stop cooling
          npartnocool = npartnocool + 1
          if (npartnocool > maxpart_nocooling) call fatal('inject_test_hiiregion','Number of particles to &
                                                           stop cooling exceeded limit')
          ipart_nocooling(npartnocool) = ip
       endif
    enddo over_parts
    print*,'number of particles heated',npartnocool
    write(2023,*) npartnocool
    close(2023)
    !
    ! Flag to avoid further injections
    !
    hii_flag = .true.
 endif
 !
 ! Timestep constraint
 !
 dtinject = huge(dtinject)

end subroutine inject_particles

real function mag(vec)
 real, intent(in) :: vec(3)

 mag = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)

end function mag

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
 call write_inopt(lumin_LyC,'lumin_LyC','Rate of LyC photon emission',iunit)
 call write_inopt(temp_hii,'temp_hii','Temperature of particles in HII region',iunit)


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
 case('lumin_LyC')
    read(valstring,*,iostat=ierr) lumin_LyC
    ngot = ngot + 1
    if (lumin_LyC <= 0.) call fatal(label,'invalid setting for lumin_LyC (<=0)')
case('temp_hii')
    read(valstring,*,iostat=ierr) temp_hii
    ngot = ngot + 1
    if (temp_hii <= 0.) call fatal(label,'invalid setting for temp_hii (<=0)')
 case default
    imatch = .false.
 end select

 noptions = 2
 igotall  = (ngot >= noptions)

end subroutine read_options_inject

end module inject
