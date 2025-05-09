!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine which generates the velocity field from the particles 
! Make three 3-D cubes - vx, vy, vz 
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
 character(len=20), parameter, public :: analysistype = '3d_velfield'

 public :: do_analysis

 private

 integer, parameter :: nrefpoint = 8E6 
 integer :: npointx = 200   ! number of points on x-axis 
 integer :: npointy = 200   
 integer :: npointz = 200  
 real    :: centre_in(3) = (/ 0., 0., 0. /)

 real    :: maxr = 10.00
 logical :: use_whole_box = .false. 

 integer :: isink_centre = 14 
 logical :: use_sink = .false. 

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use linklist, only:node,ifirstincell,listneigh,set_linklist
 use kdtree,   only:getneigh
 use dim,      only:maxneigh
 use kernel,   only:get_kernel,cnormk,radkern2
 use units,    only:udist,utime,unit_velocity,unit_density,unit_pressure,unit_ergg
 use io,       only:fatal,warning 
  use part,    only:xyzmh_ptmass,vxyz_ptmass,nptmass
 use part,     only:hfact,rhoh,massoftype,igas
 use eos,      only:gamma
 use omp_lib 
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer, parameter :: neighcachesize = 1E5
 integer :: ineigh,nneigh,ixyzcachesize,n,ip,idone
 real    :: xyzcache(3,neighcachesize)
 real    :: vx_sum,vy_sum,vz_sum,vx_target,vy_target,vz_target
 real    :: dr2,q2,q,wkern,wkern_norm,grkern
 real    :: xyz_b(3),h_b,vx_b,vy_b,vz_b,rho_b
 integer :: i,iref,nref,ipointx,ipointy,ipointz
 real    :: time_cgs,pmass,hmean,radneigh 
 real    :: xmin,xmax,ymin,ymax,zmin,zmax,sizex,sizey,sizez,xloc,yloc,zloc,dx,dy,dz
 integer :: ixyz_ref(3,nrefpoint),ixyz_target(3),centre(3) 
 real    :: xyz_ref(3,nrefpoint),vxyz_ref(3,nrefpoint),xyz_target(3),xyz_target_cgs(3) 
 real    :: percentcount,percent 
 real    :: radneighfac,sep
 real,   allocatable :: dumxyzh(:,:)
 integer :: numthreads,ithread
 character(len=70) :: filenamex,filenamey,filenamez

 if (npointx*npointy*npointz > nrefpoint) then 
    print*,'required nrefpoint:',npointx*npointy*npointz
    call fatal('analysis_velfield','set a bigger nrefpoint')
 endif 
 
 !$omp parallel default(none) shared(numthreads)
 numthreads = omp_get_num_threads()
 !$omp end parallel
 print*,'Running on',numthreads,'threads'

 time_cgs  = time*utime

 !- Particle mass
 pmass = massoftype(igas)

 !- Set centre 
 if (use_sink) then 
    if (isink_centre > nptmass) call fatal('analysis_velfield','sink not found.') 
    centre = xyzmh_ptmass(1:3,isink_centre)
 else 
    centre = centre_in
 endif 

 !- Build tree 
 allocate(dumxyzh(4,npart))
 dumxyzh = xyzh
 call set_linklist(npart,npart,dumxyzh,vxyzu)

 if (use_whole_box) then 
   !- Determine box boundaries 
   xmax = -huge(xmax)
   ymax = -huge(ymax)
   zmax = -huge(zmax)
   xmin =  huge(xmin)
   ymin =  huge(ymin)
   zmin =  huge(zmin)
   !$omp parallel do default(none) shared(npart,xyzh) private(i) &
   !$omp reduction(min:xmin,ymin,zmin) &
   !$omp reduction(max:xmax,ymax,zmax) &
   !$omp schedule(runtime)
   do i = 1,npart
      xmin = min(xmin,xyzh(1,i))
      ymin = min(ymin,xyzh(2,i))
      zmin = min(zmin,xyzh(3,i))
      xmax = max(xmax,xyzh(1,i))
      ymax = max(ymax,xyzh(2,i))
      zmax = max(zmax,xyzh(3,i))
   enddo
   !$omp end parallel do
   sizex = abs(xmax - xmin)
   sizey = abs(ymax - ymin)
   sizez = abs(zmax - zmin)
 else 
    sizex = 2.*maxr 
    sizey = 2.*maxr 
    sizez = 2.*maxr 
 endif 

 !- Set up all measuring points in cube 
 nref = 0
 do ipointx = 1,npointx
   dx = sizex/real(npointx)
   xloc = -sizex/2. + dx*real(ipointx) + centre(1)
   do ipointy = 1,npointy
      dy = sizey/real(npointy)
      yloc = -sizey/2. + dy*real(ipointy)  + centre(2)
      do ipointz = 1,npointz
         dz = sizez/real(npointz)
         zloc = -sizez/2. + dz*real(ipointz) + centre(3)
         nref = nref + 1 
         xyz_ref(1:3,nref) = (/xloc,yloc,zloc/)
         ixyz_ref(1:3,nref) = (/ipointx,ipointy,ipointz/)
      enddo 
   enddo 
 enddo 

 !- Estimate compact support radius 
 sep = max(dx,dy,dz)
 hmean = 0.
 do i = 1,npart 
    if (xyzh(4,i) > tiny(hmean)) hmean = hmean + xyzh(4,i)
 enddo 
 hmean = hmean/npart 
 radneighfac = 5.*hmean/sep 


 print*,'Begin interpolation'
 idone = 0
 percentcount = 0. 

 ! Interpolate fluid properties at each ref point 
 !$omp parallel do default(none) shared(nref,ixyz_ref,xyz_ref,vxyz_ref,npart,xyzh,vxyzu,pmass) &
 !$omp shared(radneigh,ifirstincell,node,unit_velocity,radneighfac,sep,numthreads) &
 !$omp private(iref,ixyz_target,xyz_target,xyz_target_cgs,nneigh,n,ineigh,ip,xyzcache) &
 !$omp private(vx_sum,vy_sum,vz_sum,vx_b,vy_b,vz_b,vx_target,vy_target,vz_target) &
 !$omp private(xyz_b,h_b,rho_b,dr2,q2,q,wkern,grkern,wkern_norm,percent,idone,percentcount,ithread) &
 !$omp schedule(runtime)
 do iref = 1,nref 

!    !$omp atomic update 
!    idone = idone + 1 
!    percent = 100.*real(idone)/real(nref/numthreads)
!    if (percent > percentcount) then 
!       print*,'loading ',nint(percent),'%'
!       !$omp atomic update 
!       percentcount = percentcount + 2. 
!    endif 

    ixyz_target(1:3) = ixyz_ref(1:3,iref)
    xyz_target(1:3) = xyz_ref(1:3,iref)

    nneigh = 0 
    n = 0
    radneigh = radneighfac * sep 
    do while (nneigh < 50)
       !- Get list of neighbours around detector point 
       call getneigh(node,xyz_target,0.,radneigh,3,listneigh,nneigh,xyzh,xyzcache,neighcachesize,ifirstincell,.false.)
       if (nneigh < 10) then 
          call warning('analysis_velfield','not enough trial neighbours')
       endif 
       n = n + 1 
       radneigh = radneigh * 1.5  ! try increase 
       if (n > 1e4) call fatal('analysis_velfield','cannot find neighbours')
    enddo 

    !- Compute a component of velocity by interpolating from true neighbours 
    vx_sum = 0.
    vy_sum = 0.
    vz_sum = 0. 

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
         rho_b = rhoh(h_b,pmass)
         vx_b = vxyzu(1,ip)
         vy_b = vxyzu(2,ip)
         vz_b = vxyzu(3,ip)
         ! Compute SPH sum
         wkern_norm = cnormk/(h_b**3)*wkern 
         vx_sum = vx_sum + vx_b*pmass/rho_b*wkern_norm 
         vy_sum = vy_sum + vy_b*pmass/rho_b*wkern_norm 
         vz_sum = vz_sum + vz_b*pmass/rho_b*wkern_norm 
       endif 
    enddo over_neigh

    ! Convert to cgs units 
    vx_target  = vx_sum*unit_velocity
    vy_target  = vy_sum*unit_velocity
    vz_target  = vz_sum*unit_velocity

    vxyz_ref(1:3,iref) = (/vx_target,vy_target,vz_target/)
 enddo
 !$omp end parallel do


 filenamex = 'velfield_x_'//TRIM(dumpfile)//'.dat'
 open(unit=2206,file=filenamex,status='replace')
 write(2206,'(1a20)') 'time [s]'
 write(2206,'(1e20.10)') time_cgs
 do iref = 1,nref
    write(2206,'(3i10,4e25.10)') ixyz_ref(1:3,iref), xyz_ref(1:3,iref), vxyz_ref(1,iref)
 enddo 
 close(2206)

 filenamey = 'velfield_y_'//TRIM(dumpfile)//'.dat'
 open(unit=2207,file=filenamey,status='replace')
 write(2207,'(1a20)') 'time [s]'
 write(2207,'(1e20.10)') time_cgs
 do iref = 1,nref
    write(2207,'(3i10,4e25.10)') ixyz_ref(1:3,iref), xyz_ref(1:3,iref), vxyz_ref(2,iref)
 enddo 
 close(2207)

 filenamez = 'velfield_z_'//TRIM(dumpfile)//'.dat'
 open(unit=2208,file=filenamez,status='replace')
 write(2208,'(1a20)') 'time [s]'
 write(2208,'(1e20.10)') time_cgs
 do iref = 1,nref
    write(2208,'(3i10,4e25.10)') ixyz_ref(1:3,iref), xyz_ref(1:3,iref), vxyz_ref(3,iref)
 enddo 
 close(2208)


 deallocate(dumxyzh)

end subroutine do_analysis


real function mag2(vec)
 real,   intent(in) :: vec(3)

 mag2 = dot_product(vec,vec)

end function mag2

!--------------------------------------------------------------------------
end module analysis
