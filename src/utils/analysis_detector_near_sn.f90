!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine which measures the interpolated gas flow properties (v,rho,P)
! at a certain target point away from the supernova
!
! :References: None
!
! :Owner: Cheryl Lau
!
! :Runtime parameters: None
!
! :Dependencies:
!
 implicit none
 character(len=20), parameter, public :: analysistype = 'sn_output'

 public :: do_analysis

 private

 real    :: xyz_target(3) = (/ 10., 0., 0. /)  ! Position of detector

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use linklist, only:node,ifirstincell,listneigh
 use kdtree,   only:getneigh
 use dim,      only:maxneigh
 use kernel,   only:get_kernel,cnormk,radkern2
 use units,    only:udist,utime,unit_velocity,unit_density,unit_pressure
 use io,       only:fatal
 use part,     only:hfact,rhoh,massoftype,igas
 use eos,      only:gamma
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer, parameter :: neighcachesize = 1E5
 integer :: ip,nneigh,ixyzcachesize,ineigh,ipart
 real    :: xyzcache(3,neighcachesize)
 real    :: hmean,rad_neigh
 real    :: pmass,x,y,z,h,rho,vx,vy,vz,v,u
 real    :: rho_target,v_target,thermpr_target,rampr_target
 real    :: time_si
 real    :: xyz_target_si(3) 
 character(len=70) :: filename

 !- Particle mass
 pmass = massoftype(igas)

 !- Get list of neighbours around detector point 
 hmean = 0.
 do ip = 1,npart 
    hmean = hmean + xyzh(4,ip)
 enddo 
 hmean = hmean/npart 
 rad_neigh = hmean*2.
 call getneigh(node,xyz_target,0.,rad_neigh,3,listneigh,nneigh,xyzh,xyzcache,neighcachesize,ifirstincell,.false.)

 !- Compute properties by interpolating from true neighbours 
 vx_sum  = 0.
 rho_sum = 0.
 thermpr_sum = 0.
 rampr_sum   = 0.
 over_neigh: do ineigh = 1,nneigh
    ip = listneigh(ineigh)
    xyz_b = xyzh(1:3,ip)
    h_b = xyzh(4,ip)
    dr2 = mag2(xyz_b - xyz_target)
    q2  = dr2 / (h_b**2) 
    if (q2 < radkern2) then !- within compact support radius
       q = sqrt(q2) 
       call get_kernel(q2,q,wkern,grkern)
       ! get neigh particle properties 
       vx_b  = vxyzu(1,ip)
       rho_b = rhoh(h_b,pmass)
       rampr_b   = rho_b*mag2(vxyzu(1:3,ip))
       thermpr_b = rho_b*(gamma-1.)*vxyzu(4,ip)
       ! Compute SPH sum
       wkern_norm = cnormk/(h_b**3)*wkern 
       vx_sum  = vx_sum + vx_b*pmass/rho_b*wkern_norm 
       rho_sum = rho_sum + pmass*wkern_norm
       rampr_sum   = rampr_sum + rampr_b*pmass/rho_b*wkern_norm 
       thermpr_sum = thermpr_sum + rampr_b*pmass/rho_b*wkern_norm 
    endif 
 enddo over_neigh

 ! Convert to SI units 
 rho_target = rho_sum*unit_density*1000.
 vx_target  = vx_sum*unit_velocity*0.01
 rampr_target = rampr_sum*unit_pressure*0.1
 thermpr_target = thermpr_target*unit_pressure*0.1
 do i = 1,3
    xyz_target(i) = xyz_target(i)*udist*0.01
 enddo 
 time_si  = time*utime

 filename = 'gasflow_'//TRIM(dumpfile)//'.dat'
 open(unit=2206,file=filename)
 write(2206,'(4a20)') 'time [s]','x [m]','y [m]','z [m]'
 write(2206,*) time_si, xyz_target(1:3)
 write(2206,'(4a25)') 'rho [kg m^-3]','v_x [m s^-1]','therm pr [kg m^-1 s^-2]','ram pr [kg m^-1 s^-2]'
 write(2206,*) rho_target, vx_target, thermpr_target, rampr_target 

 close(2206)

end subroutine do_analysis


real function mag2(vec)
 real,   intent(in) :: vec(3)

 mag2 = dot_product(vec,vec)

end function mag2

!--------------------------------------------------------------------------
end module analysis
