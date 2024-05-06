!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine which measures the interpolated gas flow properties (v,rho,P)
! along a line from the supernova 
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
 character(len=20), parameter, public :: analysistype = 'sn_shock'

 public :: do_analysis

 private

 integer :: npoint = 200
 real    :: maxr = 70. 

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use linklist, only:node,ifirstincell,listneigh,set_linklist
 use kdtree,   only:getneigh
 use dim,      only:maxneigh
 use kernel,   only:get_kernel,cnormk,radkern2
 use units,    only:udist,utime,unit_velocity,unit_density,unit_pressure,unit_ergg
 use io,       only:fatal,warning 
 use part,     only:hfact,rhoh,massoftype,igas
 use eos,      only:gamma
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer, parameter :: neighcachesize = 1E5
 integer :: i,ip,iclosest,ineigh,nneigh,ixyzcachesize,ipoint 
 real    :: xyzcache(3,neighcachesize)
 real    :: pmass,hmean,dist2,dist2_min,rad_neigh
 real    :: vx_sum,u_sum,rho_sum,thermpr_sum,rampr_sum
 real    :: dr2,q2,q,wkern,wkern_norm,grkern
 real    :: xyz_b(3),h_b,vx_b,u_b,rho_b,rampr_b,thermpr_b
 real    :: rho_target,vx_target,u_target,thermpr_target,rampr_target
 real    :: time_cgs,xyz_target_cgs(3),xyz_target(3)
 real    :: xloc
 real,   allocatable :: dumxyzh(:,:)
 character(len=70) :: filename

 filename = 'shock1d_'//TRIM(dumpfile)//'.dat'
 open(unit=2206,file=filename,status='replace')

 time_cgs  = time*utime

 write(2206,'(1a20)') 'time [s]'
 write(2206,'(1e20.10)') time_cgs
 write(2206,'(6a25)') 'x [cm]','rho [g cm^-3]','v_x [cm s^-1]','u [erg g^-1]','therm pr [g cm^-1 s^-2]','ram pr [g cm^-1 s^-2]'
 close(2206)

 !- Particle mass
 pmass = massoftype(igas)

 !- Build tree 
 allocate(dumxyzh(4,npart))
 dumxyzh = xyzh
 call set_linklist(npart,npart,dumxyzh,vxyzu)

 over_points: do ipoint = 1,npoint
   xloc = maxr/real(npoint)*real(ipoint) 
   xyz_target(1:3) = (/ xloc,0.,0./)

   !- Estimate the compact support radius at target point 
   dist2_min = huge(dist2_min)
   do ip = 1,npart 
      if (xyzh(4,ip) > tiny(dist2)) then  ! nearest alive particle 
         dist2 = mag2(xyzh(1:3,ip)-xyz_target(1:3))
         if (dist2 < dist2_min) then 
            dist2_min = dist2 
            iclosest  = ip
         endif 
      endif 
   enddo 
   rad_neigh = xyzh(4,iclosest) * 2. 
   if (rad_neigh < tiny(dist2)) call fatal('analysis_1d_shock','rad_neigh = 0')

   !- Get list of neighbours around detector point 
   call getneigh(node,xyz_target,0.,rad_neigh,3,listneigh,nneigh,xyzh,xyzcache,neighcachesize,ifirstincell,.false.)
   if (nneigh < 50) then 
      print*,'rad_neigh,nneigh',rad_neigh,nneigh
      call warning('analysis_1d_shock','not enough trial neighbours')
   endif 

   !- Compute properties by interpolating from true neighbours 
   vx_sum  = 0.
   u_sum = 0.
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
         vx_b  = sqrt(mag2(vxyzu(1:3,ip)))  ! vxyzu(1,ip)
         u_b   = vxyzu(4,ip)
         rho_b = rhoh(h_b,pmass)
         rampr_b   = rho_b*mag2(vxyzu(1:3,ip))    ! rho*v2 
         thermpr_b = rho_b*(gamma-1.)*vxyzu(4,ip) ! rho*(gamma-1)*u
         ! Compute SPH sum
         wkern_norm = cnormk/(h_b**3)*wkern 
         vx_sum  = vx_sum + vx_b*pmass/rho_b*wkern_norm 
         u_sum   = u_sum + u_b*pmass/rho_b*wkern_norm 
         rho_sum = rho_sum + pmass*wkern_norm
         rampr_sum   = rampr_sum + rampr_b*pmass/rho_b*wkern_norm 
         thermpr_sum = thermpr_sum + thermpr_b*pmass/rho_b*wkern_norm 
      endif 
   enddo over_neigh

   ! Convert to cgs units 
   rho_target = rho_sum*unit_density
   vx_target  = vx_sum*unit_velocity
   u_target   = u_sum*unit_ergg
   rampr_target = rampr_sum*unit_pressure
   thermpr_target = thermpr_sum*unit_pressure
   do i = 1,3
      xyz_target_cgs(i) = xyz_target(i)*udist
   enddo 

   open(unit=2206,file=filename,position='append')
   write(2206,'(6e25.10)') xyz_target_cgs(1), rho_target, vx_target, u_target, thermpr_target, rampr_target 
   close(2206)

 enddo over_points 



 deallocate(dumxyzh)

end subroutine do_analysis


real function mag2(vec)
 real,   intent(in) :: vec(3)

 mag2 = dot_product(vec,vec)

end function mag2

!--------------------------------------------------------------------------
end module analysis
