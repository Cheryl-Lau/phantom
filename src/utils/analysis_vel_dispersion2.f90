!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine which measures the velocity dispersion of all/high-density particles
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
 character(len=20), parameter, public :: analysistype = 'vel_dispersion'

 public :: do_analysis

 private

 real    :: maxrad = 4e1
 real    :: rholimit_cgs  = 1e-21
 real    :: templimit     = 8000 
 logical :: only_highrho  = .true. 
 logical :: only_hightemp = .false. 
 
contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use linklist, only:node,ifirstincell,listneigh,set_linklist
 use kdtree,   only:getneigh
 use dim,      only:maxneigh
 use units,    only:udist,utime,umass,unit_velocity,unit_density,unit_pressure,unit_ergg
 use io,       only:fatal,warning 
 use part,     only:hfact,rhoh,massoftype,igas
 use eos,      only:gamma,gmw 
 use physcon,  only:gg,pi,kboltz,mass_proton_cgs
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer, parameter :: neighcachesize = 1E6
 integer :: ip,ineigh,nneigh,ixyzcachesize,n,ip_neigh,consneigh
 real    :: xyzcache(3,neighcachesize)
 real    :: dist2_allneigh(neighcachesize)
 real    :: pmass,dist2,rho,rad2_limit,mean_v,sigma_v,mean_rho,rad,temp,totmass,virial 
 real    :: rad_cgs,totmass_cgs,mean_rho_cgs,sigma_v_cgs,virial_cgs
 real,   allocatable :: dumxyzh(:,:)
 character(len=70) :: filename

 !- Particle mass
 pmass = massoftype(igas)

 !- Build tree 
 allocate(dumxyzh(4,npart))
 dumxyzh = xyzh
 call set_linklist(npart,npart,dumxyzh,vxyzu)

 open(unit=2026,file='velocity_dispersion2_'//TRIM(dumpfile)//'.dat',status='replace')
 open(unit=2027,file='virial_term2_'//TRIM(dumpfile)//'.dat',status='replace')
 open(unit=2028,file='density_sizescale2_'//TRIM(dumpfile)//'.dat',status='replace')

 !$omp parallel do default(none) shared(npart,pmass,xyzh,vxyzu,rholimit_cgs,templimit) &
 !$omp shared(node,hfact,maxrad,xyzcache,ifirstincell,only_highrho,only_hightemp) &
 !$omp shared(umass,udist,unit_velocity,unit_density,unit_ergg,gamma,gmw) &
 !$omp private(ip,rho,temp,nneigh,rad,rad2_limit,totmass,mean_v,sigma_v,mean_rho,ineigh,ip_neigh,dist2,virial) &
 !$omp private(dist2_allneigh,consneigh) &
 !$omp private (rad_cgs,totmass_cgs,mean_rho_cgs,sigma_v_cgs,virial_cgs) &
 !$omp schedule(runtime)
 over_part: do ip = 1,npart
    print*,'ip/npart',ip,'/',npart

    rho  = pmass*(hfact/abs(xyzh(4,ip)))**3
    temp = vxyzu(4,ip)/kboltz*(gmw*mass_proton_cgs*(gamma-1.))*unit_ergg

    if (only_highrho  .and. rho  < rholimit_cgs/unit_density) cycle over_part
    if (only_hightemp .and. temp < templimit)                 cycle over_part

    ! Get all neigh 
    call getneigh(node,xyzh(1:3,ip),0.,maxrad,3,listneigh,nneigh,xyzh,xyzcache,neighcachesize,ifirstincell,.false.)
    if (nneigh == 0) cycle over_part

    !- Calculate distances to target particle 
    do ineigh = 1,nneigh
       ip_neigh = listneigh(ineigh)
       dist2 = mag2(xyzh(1:3,ip) - xyzh(1:3,ip_neigh))
       dist2_allneigh(ineigh) = dist2  
    enddo 

    !- Sort all neighbours by distance
    call quick_sort(nneigh,dist2_allneigh,listneigh,1,nneigh)

    irad = 0
    consneigh = 1
    more_neigh: do while (consneigh <= nneigh)
       irad = irad + 1 
       if (irad > maxrad) exit more_neigh

       ! consider the first consneigh neighbours
       rad = sqrt(dist2_allneigh(consneigh))
       totmass = pmass*consneigh
       !- Calculte mean 
       mean_v   = 0. 
       mean_rho = 0.
       do ineigh = 1,consneigh 
          ip_neigh = listneigh(ineigh)
          mean_v   = mean_v + sqrt(mag2(vxyzu(1:3,ip_neigh)))
          mean_rho = mean_rho + pmass*(hfact/abs(xyzh(4,ip_neigh)))**3
       enddo 
       mean_v   = mean_v/real(consneigh) 
       mean_rho = mean_rho/real(consneigh) 
       !- Calculate std 
       sigma_v = 0.
       do ineigh = 1,consneigh 
          ip_neigh = listneigh(ineigh)
          sigma_v  = sigma_v + (sqrt(mag2(vxyzu(1:3,ip_neigh))) - mean_v)**2 
       enddo 
       sigma_v = sqrt(sigma_v/real(consneigh)) 
       if (sigma_v > 0) virial = 2.d0*(consneigh*pmass)/(sigma_v**2*rad)

       !- Convert to physical units 
       rad_cgs = rad * udist 
       totmass_cgs  = totmass * umass 
       mean_rho_cgs = mean_rho * unit_density 
       sigma_v_cgs  = sigma_v * unit_velocity 
       if (sigma_v > 0) then 
          virial_cgs = virial * (gg*umass) / (unit_velocity**2 * udist)
       else 
          virial_cgs = -1 
       endif 

       rad_all(irad) = rad_cgs
       mean_rho_all(irad) = mean_rho_cgs
       sigma_v_all(irad) = sigma_v_cgs
       virial_all(irad) = virial_cgs
       
       ! Measure next sigma at radius that corresponds to this many extra neighbours 
       consneigh = consneigh + 50*int(irad**3) 

    enddo more_neigh

    !$omp critical 
    write(2025,*) rad_all(irad)
    write(2026,*) rad_cgs, sigma_v_cgs
    write(2027,*) rad_cgs, virial_cgs
    write(2028,*) rad_cgs, mean_rho_cgs
    !$omp end critical


 enddo over_part
 !$omp end parallel do

 close(2026)
 close(2027)
 close(2028) 
 deallocate(dumxyzh)

end subroutine do_analysis

!-----------------------------------------------------------------------
!+
! Tools 
!+
!-----------------------------------------------------------------------
real function mag2(vec)
 real,   intent(in) :: vec(3)

 mag2 = dot_product(vec,vec)

end function mag2

!
! Routines to sort input array along with iarray (carrying indices)
!
recursive subroutine quick_sort(n,array,iarray,first,last)
 integer, intent(in)    :: n
 integer, intent(in)    :: first,last
 integer, intent(inout) :: iarray(n)
 real,    intent(inout) :: array(n)
 integer :: partition,nleft,nright

 if (first < last .and. n > 0) then
     !- Set pivot
    call partition_pos(n,array,iarray,first,last,partition)

    nleft = partition - first
    call quick_sort(nleft,array,iarray,first,partition-1)

    nright = last - partition
    call quick_sort(nright,array,iarray,partition+1,last)
 endif

end subroutine quick_sort

!
! Function for finding a pivot such that those smaller than pivot
! would be on the left and vice versa
!
subroutine partition_pos(n,array,iarray,first,last,partition)
 integer, intent(in)    :: n
 integer, intent(in)    :: first,last
 integer, intent(inout) :: iarray(n)
 real,    intent(inout) :: array(n)
 integer, intent(out)   :: partition
 integer :: i,j,itemp
 real    :: temp,pivot

 pivot = array(last)
 i = first - 1        !- pointer for greater element

 do j = first,last
    if (array(j) < pivot) then
       i = i + 1
       temp = array(i)
       array(i) = array(j)
       array(j) = temp

       itemp = iarray(i)
       iarray(i) = iarray(j)
       iarray(j) = itemp
    endif
 enddo
 temp = array(i+1)
 array(i+1) = array(last)
 array(last) = temp

 itemp = iarray(i+1)
 iarray(i+1) = iarray(last)
 iarray(last) = itemp

 partition = i + 1

end subroutine partition_pos

!--------------------------------------------------------------------------
end module analysis
