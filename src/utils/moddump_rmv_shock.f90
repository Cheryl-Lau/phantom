!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module moddump
!
! Routine to suppress a shock by redistributing particle energies
!
! :References: None
!
! :Owner: Cheryl Lau
!
! :Runtime parameters: 
!
! :Dependencies: None
!
 implicit none

 real :: r_rmshock = 10.0

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,       only:igas,kill_particle
 use partinject, only:add_or_update_particle
 use physcon,    only:pi
 use io,         only:fatal
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: ip,np_rmshock,np_placed,ipartadd
 real    :: ran_x,ran_y,ran_z
 real    :: utot,umean,hmean,r,xyz(3),vxyz(3),x,y,z,maxu

 maxu = tiny(maxu)
 utot = 0
 hmean = 0
 np_rmshock = 0
 do ip = 1,npart 
    r = sqrt(dot_product(xyzh(1:3,ip),xyzh(1:3,ip)))
    if (r < r_rmshock) then 
       utot = utot + vxyzu(4,ip) 
       hmean = hmean + xyzh(4,ip)
       np_rmshock = np_rmshock + 1
       maxu = max(maxu,vxyzu(4,ip))
       call kill_particle(ip) 
    endif 
 enddo 
 print*,'maxu before',maxu

 umean = utot / np_rmshock
 hmean = hmean / np_rmshock

 maxu = tiny(maxu)
 np_placed = 0
 ipartadd = npart
 do while (np_placed <= np_rmshock)
    call random_number(ran_x)
    call random_number(ran_y)
    call random_number(ran_z)
    x = (2.*ran_x-1)*r_rmshock
    y = (2.*ran_y-1)*r_rmshock
    z = (2.*ran_z-1)*r_rmshock
    if (sqrt(x**2 + y**2 + z**2) < r_rmshock) then 
        np_placed = np_placed + 1 
        ipartadd = ipartadd + 1 
        xyz = (/ x, y, z /)
        vxyz = (/ 0., 0., 0. /)
        call add_or_update_particle(igas,xyz,vxyz,hmean,umean,ipartadd,npart,&
                                    npartoftype,xyzh,vxyzu)
        maxu = max(maxu,umean)
    endif 
 enddo
 print*,'maxu after',maxu
 
 npartoftype(:) = 0
 npartoftype(igas) = npart

 return
end subroutine modify_dump

end module moddump
