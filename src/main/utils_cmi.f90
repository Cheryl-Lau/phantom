!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module utils_cmi
!
! The CMI suite: photoionize_cmi.F90 kdtree_cmi.f90 hnode_cmi.f90 heating_cooling_cmi.f90 
!                *utils_cmi.f90*
! This module contains all the miscellaneous subroutines required by the above modules 
!
! :References: none
!
! :Owner: Cheryl Lau
!
! :Dependencies: io
!
 implicit none

 public :: modify_grid,set_bounds
 public :: mag2,quick_sort,gen_filename

 private

contains

!-----------------------------------------------------------------------
!+
! Routine to merge the smallest Voronoi grid cells 
! done to avoid seg-fault in CMI grid initialization stage 
!+
!-----------------------------------------------------------------------
subroutine modify_grid(npart,x,y,z,h,hlimit_fac,extradist_fac)
 use io,   only:warning,fatal
 integer, intent(inout) :: npart
 real,    intent(inout) :: x(npart),y(npart),z(npart)
 real,    intent(in)    :: h(npart)
 real,    intent(in)    :: hlimit_fac
 real,    intent(in)    :: extradist_fac
 integer, parameter :: maxgrp = 10000
 integer, parameter :: maxp_per_grp = 5000
 integer :: ip,ip_neigh,npart_smallh,npart_mod,npart_merged
 integer :: k,ngroup,np_in_k,i_in_k
 integer :: parts_grp(maxgrp,maxp_per_grp)  ! stores index of particles in each group 
 integer :: npart_grp(maxgrp)               ! stores number of particles in each group 
 real    :: hmin,hmax,hlimit 
 real    :: xyz_ip(3),xyz_neigh(3),dist2,xmean,ymean,zmean
 real,    allocatable :: x_mod(:),y_mod(:),z_mod(:)
 logical :: flag_particle(npart)
 logical :: check_modgrid = .true.
 
 ! testing 
 if (check_modgrid) then 
    open(2050,file='before_mod_grid.txt')
    do ip = 1,npart
       write(2050,*) x(ip), y(ip), z(ip)
    enddo 
    close(2050)
 endif 

 !- find the max and min of h
 hmin = huge(hmin)
 hmax = tiny(hmax)
 do ip = 1,npart
    hmin = min(hmin,h(ip))
    hmax = max(hmax,h(ip))
 enddo 

 !- set filter - h smaller than this limit will need to be dealt with 
 hlimit = hmin + (hmax-hmin)*hlimit_fac

 !- Init arrays for storing groups
 ngroup = 0
 parts_grp(:,:) = 0. 
 npart_grp(:)   = 0 
 flag_particle(:) = .false. 

 !- Init for the final set of particles
 npart_mod = 0
 allocate(x_mod(npart))
 allocate(y_mod(npart))
 allocate(z_mod(npart))


 k = 0  ! group index
 do ip = 1,npart
    if (h(ip) < hlimit) then 
       if (.not.flag_particle(ip)) then 
          !- start new group
          k = k + 1 
          if (k > maxgrp) call fatal('utils_cmi','number of groups exceeded limit')
          !- put first particle in
          np_in_k = 1 
          parts_grp(k,np_in_k) = ip 
          npart_grp(k)         = np_in_k
          !- loop over all other particles 
          xyz_ip = (/ x(ip), y(ip), z(ip) /)
          over_neigh: do ip_neigh = 1,npart 
             if (ip_neigh == ip) cycle over_neigh 
             if (h(ip_neigh) < hlimit .and. .not.flag_particle(ip_neigh)) then 
                xyz_neigh = (/ x(ip_neigh), y(ip_neigh), z(ip_neigh) /)
                dist2 = mag2(xyz_ip - xyz_neigh)
                !- extract those with overlapping h
                if (dist2 < ((h(ip) + h(ip_neigh))*extradist_fac)**2) then 
                   np_in_k = np_in_k + 1
                   if (np_in_k > maxp_per_grp) call fatal('utils_cmi','number of particles in group exceeded limit')
                   parts_grp(k,np_in_k) = ip_neigh 
                   npart_grp(k)         = np_in_k
                   flag_particle(ip_neigh) = .true.  ! to avoid being placed in group again
                endif 
             endif 
          enddo over_neigh 
       endif 
    elseif (h(ip) >= hlimit) then  
       !- store directly 
       npart_mod = npart_mod + 1 
       x_mod(npart_mod) = x(ip)
       y_mod(npart_mod) = y(ip)
       z_mod(npart_mod) = z(ip)
    endif 
 enddo 

 ngroup = k

 !- Loop through each group to locate its centre 
 npart_merged = 0
 do k = 1,ngroup 
    np_in_k = npart_grp(k)
    npart_merged = npart_merged + np_in_k
    xmean = 0.
    ymean = 0.
    zmean = 0. 
    over_parts: do i_in_k = 1,np_in_k
       ip = parts_grp(k,i_in_k)  ! extract particle index 
       xmean = xmean + x(ip)
       ymean = ymean + y(ip)
       zmean = zmean + z(ip) 
    enddo over_parts
    xmean = xmean/real(np_in_k)
    ymean = ymean/real(np_in_k)
    zmean = zmean/real(np_in_k)
    !- store it as a new particle 
    npart_mod = npart_mod + 1 
    x_mod(npart_mod) = xmean
    y_mod(npart_mod) = ymean
    z_mod(npart_mod) = zmean
 enddo 

 if (npart_mod == npart) then
    call warning('utils_cmi','no cells merged')
 else 
    write(*,'(2x,i6,a38,i5)') npart_merged,' Voronoi generation sites merged into ',ngroup
    write(*,'(2x,a20,i7,a4,i7)') 'Changing nsite from ',npart,' to ',npart_mod
 endif 
 if ((npart-npart_mod)/npart > 0.05) call fatal('utils_cmi','merged too many small cells!') 

 !- Output the modified sites 
 npart = npart_mod
 x(1:npart) = x_mod(1:npart)
 y(1:npart) = y_mod(1:npart)
 z(1:npart) = z_mod(1:npart)

 deallocate(x_mod)
 deallocate(y_mod)
 deallocate(z_mod)

 ! testing 
 if (check_modgrid) then 
    open(2070,file='after_mod_grid.txt')
    do ip = 1,npart
       write(2070,*) x(ip), y(ip), z(ip)
    enddo 
    close(2070)
 endif 

end subroutine modify_grid 

!-----------------------------------------------------------------------
!+
! Find boundaries with given xyz
!+
!-----------------------------------------------------------------------
subroutine set_bounds(nsite,x,y,z,h,m,xmin,xmax,ymin,ymax,zmin,zmax,dx,dy,dz)
 use io,   only:fatal
 integer, intent(in)  :: nsite 
 real,    intent(in)  :: x(nsite),y(nsite),z(nsite),h(nsite),m(nsite) 
 real,    intent(out) :: xmin,xmax,ymin,ymax,zmin,zmax,dx,dy,dz 
 integer :: i 
 
 xmax = -huge(xmax)
 ymax = -huge(ymax)
 zmax = -huge(zmax)
 xmin =  huge(xmin)
 ymin =  huge(ymin)
 zmin =  huge(zmin)
 !$omp parallel do default(none) shared(nsite,x,y,z,h,m) private(i) &
 !$omp reduction(min:xmin,ymin,zmin) &
 !$omp reduction(max:xmax,ymax,zmax) &
 !$omp schedule(runtime)
 do i = 1,nsite
    if (x(i) /= x(i) .or. y(i) /= y(i) .or. z(i) /= z(i) .or. &
        h(i) < tiny(h) .or. m(i) < tiny(m)) then
       print*,x(i),y(i),z(i),h(i),m(i)
       call fatal('utils_cmi','invalid xyzhm input')
    endif
    xmin = min(xmin,x(i))
    ymin = min(ymin,y(i))
    zmin = min(zmin,z(i))
    xmax = max(xmax,x(i))
    ymax = max(ymax,y(i))
    zmax = max(zmax,z(i))
 enddo
 !$omp end parallel do
 dx = abs(xmax - xmin)
 dy = abs(ymax - ymin)
 dz = abs(zmax - zmin)

end subroutine set_bounds 

!-----------------------------------------------------------------------
!+
! Routines to sort input array along with iarray (carrying indices)
!+
!-----------------------------------------------------------------------
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

!-----------------------------------------------------------------------
!+
! Math tools 
!+
!-----------------------------------------------------------------------
real function mag2(vec)
 real,   intent(in) :: vec(3)

 mag2 = dot_product(vec,vec)

end function mag2

!-----------------------------------------------------------------------
!+
! Routine to generate output file name 
!+
!-----------------------------------------------------------------------
subroutine gen_filename(photoionize_tree,ifile,filename)
 integer,           intent(in)  :: ifile
 logical,           intent(in)  :: photoionize_tree
 character(len=50), intent(out) :: filename
 character(len=5)   :: ifile_char

 write(ifile_char,'(i5.5)') ifile  !- convert to str
 if (photoionize_tree) then
    filename = 'nixyzhmf_'//trim(ifile_char)//'.txt'
 else
    filename = 'xyzhmf_'//trim(ifile_char)//'.txt'
 endif

end subroutine gen_filename

end module utils_cmi
