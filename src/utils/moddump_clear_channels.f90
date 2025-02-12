!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module moddump
!
! Routine to manually carve out conic channels in a spherical HII region
! assumming the HII region is centered at the origin
!
! :References: None
!
! :Owner: Cheryl Lau
!
! :Runtime parameters: 
!   - rad_chnl          : *maximum radius of the channell to clear*
!   - omegafrac_chnl    : *size of channell as a fraction of sphere's total solid angle*
!   - pfrac_chnl        : *fraction of particles to remove within channel*
!
! :Dependencies: None
!
 implicit none

 real   :: rad_chnl = 12.0
 real   :: omegafrac_chnl = 0.05
 real   :: pfrac_chnl = 0.9

 integer, parameter :: nchnl = 4
 real   :: rchnl_vec(3,nchnl) = reshape((/ 1,0,0, &
                                          -1, 0,0, &
                                          0, 1,0, &
                                          0, -1,0 /), &
                                         shape=(/3,nchnl/))

!- Combos 
!  1  0  0
! -1  0  0 
!  0  1  0
!  0 -1  0
!  0  0  1
!  0  0 -1

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use random,  only:ran2
 use part,    only:igas
 use physcon, only:pi
 use io,      only:fatal
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,    allocatable   :: xyzh_store(:,:),vxyzu_store(:,:)
 integer :: maxpart,npart_store,ip,ichnl 
 integer :: iseed = -12345
 real    :: omega,omega_chnl,xp,yp,zp,xc,yc,zc
 real    :: rp(3),rc(3),absrp,absrc,theta_pc 
 logical :: add_particle

 omega_chnl = omegafrac_chnl * 4.*pi

 maxpart = npart
 allocate(xyzh_store(4,maxpart))
 allocate(vxyzu_store(4,maxpart))
 xyzh_store(:,:)  = 0. 
 vxyzu_store(:,:) = 0.
 npart_store = 0

 !- Go through each particle and only store those which are outside the channels
 over_parts: do ip = 1,npart 
    xp = xyzh(1,ip)
    yp = xyzh(2,ip)
    zp = xyzh(3,ip)
    rp = (/ xp,yp,zp /)
    absrp = sqrt(mag2(rp))

    add_particle = .true. 
    if (absrp < rad_chnl) then 

       each_chnl: do ichnl = 1,nchnl 
          xc = rchnl_vec(1,ichnl)
          yc = rchnl_vec(2,ichnl)
          zc = rchnl_vec(3,ichnl)
          rc = (/ xc,yc,zc /) *rad_chnl
          absrc = sqrt(mag2(rc))

          theta_pc = acos(dot_product(rp,rc)/(absrp*absrc))  ! angle between rp and rc 
          if (theta_pc < pi/2.d0) then                       ! within 90 deg 
             omega = 2.d0*pi*(1.d0-cos(theta_pc))            ! solid angle subtended by cone of 2.d0*theta_pc 
             if (omega < omega_chnl) then
                if (ran2(iseed) < pfrac_chnl) add_particle = .false.
             endif
          endif 
       enddo each_chnl 
    endif

    if (add_particle) then 
       npart_store = npart_store + 1
       xyzh_store(:,npart_store)  = xyzh(:,ip)
       vxyzu_store(:,npart_store) = vxyzu(:,ip)
    endif 
 enddo over_parts 
 
 print*,'Number of particles before and after carving: ',npart,npart_store
 if (npart == npart_store) call fatal('moddump_clear_channel','no particles removed')

 npart = npart_store
 xyzh(:,:)  = xyzh_store(:,:)
 vxyzu(:,:) = vxyzu_store(:,:)

 npartoftype(:) = 0
 npartoftype(igas) = npart

 return
end subroutine modify_dump


real function mag2(vec)
 real,   intent(in) :: vec(3)

 mag2 = dot_product(vec,vec)

end function mag2


end module moddump
