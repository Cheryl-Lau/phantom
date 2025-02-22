!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Semi-confined SN
! Analysis routine which determines the location of the SN shell by measuring 
! the interpolated u along the opposite direction of the outflow 
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
 character(len=20), parameter, public :: analysistype = 'sn_shell'

 public :: do_analysis

 private

 integer :: npoint = 300
 real    :: maxr = 30. 

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
 integer :: i,ip,iclosest,ineigh,nneigh,ixyzcachesize,ipoint,n
 real    :: xyzcache(3,neighcachesize)
 real    :: pmass,hmean,dist2,dist2_min,rad_neigh
 real    :: u_sum,rho_sum
 real    :: dr2,q2,q,wkern,wkern_norm,grkern
 real    :: xyz_b(3),h_b,u_b,rho_b
 real    :: rho_target,u_target
 real    :: time_cgs,xyz_target_cgs(3),xyz_target(3)
 real    :: xloc
 real,   allocatable :: dumxyzh(:,:)
 character(len=70) :: filename

 filename = 'shell_'//TRIM(dumpfile)//'.dat'
 open(unit=2206,file=filename,status='replace')

 time_cgs  = time*utime

 write(2206,'(1a20)') 'time [s]'
 write(2206,'(1e20.10)') time_cgs
 write(2206,'(3a25)') 'x [cm]','rho [g cm^-3]','u [erg g^-1]'
 close(2206)

 !- Particle mass
 pmass = massoftype(igas)

 !- Build tree 
 allocate(dumxyzh(4,npart))
 dumxyzh = xyzh
 call set_linklist(npart,npart,dumxyzh,vxyzu)

 over_points: do ipoint = 1,npoint
   xloc = -1 * maxr/real(npoint)*real(ipoint) 
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
   if (rad_neigh < tiny(dist2)) call fatal('analysis_shell','rad_neigh = 0')

   nneigh = 0 
   n = 0
   do while (nneigh < 50)
      rad_neigh = rad_neigh * 1.1  ! try increase 
      !- Get list of neighbours around detector point 
      call getneigh(node,xyz_target,0.,rad_neigh,3,listneigh,nneigh,xyzh,xyzcache,neighcachesize,ifirstincell,.false.)
      if (nneigh < 50) then 
         print*,'rad_neigh,nneigh',rad_neigh,nneigh
         call warning('analysis_shell','not enough trial neighbours')
      endif 
      n = n + 1 
      if (n > 1000) call fatal('analysis_shell','cannot find neighbours')
   enddo 

   !- Compute properties by interpolating from true neighbours 
   u_sum = 0.
   rho_sum = 0.
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
         u_b   = vxyzu(4,ip)
         rho_b = rhoh(h_b,pmass)
         ! Compute SPH sum
         wkern_norm = cnormk/(h_b**3)*wkern 
         u_sum   = u_sum + u_b*pmass/rho_b*wkern_norm 
         rho_sum = rho_sum + pmass*wkern_norm
      endif 
   enddo over_neigh

   ! Convert to cgs units 
   rho_target = rho_sum*unit_density
   u_target   = u_sum*unit_ergg
   do i = 1,3
      xyz_target_cgs(i) = xyz_target(i)*udist
   enddo 

   open(unit=2206,file=filename,position='append')
   write(2206,'(3e25.10)') xyz_target_cgs(1), rho_target, u_target
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
