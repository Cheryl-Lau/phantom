!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module kdtree_cmi
!
! The CMI suite: photoionize_cmi.f90 *kdtree_cmi.f90* hnode_cmi.f90 heating_cooling_cmi.f90
! This module contains all the subroutines necessary for extracting nodes from
! the built kdtree to be then passed to the photoionization code CMacIonize.
!
! :References: none
!
! :Owner: Cheryl Lau
!
! :Dependencies: io, part, dtypekdtree, kdtree
!
 implicit none

 public :: extract_cminodes_from_tree

 private

contains
!-----------------------------------------------------------------------
!+
! Walks tree to extract higher-level nodes relative to the ionizing sources;
! Uses nodes instead of individual particles to contruct grid and
! perform density-mapping.
! Includes additional opening criteria that come from CMI outputs.
!+
!-----------------------------------------------------------------------
subroutine extract_cminodes_from_tree(xyz_photosrc,nphotosrc,&
                                      node,ifirstincell,xyzh,&
                                      nxyzm_treetocmi,ixyzhm_leafparts,nnode_toreplace,&
                                      ncminode,nleafparts,ncloseleaf,&
                                      maxcminode,maxleafparts,&
                                      tree_accuracy_cmi,mintheta,rcut_opennode,rcut_leafpart,&
                                      nnode_needopen,nneedopen,&
                                      nxyzrs_nextopen,nnextopen)
 use part,        only:massoftype,igas
 use dtypekdtree, only:kdnode
 use kdtree,      only:irootnode,inodeparts,inoderange
 use kdtree,      only:maxlevel,maxlevel_indexed
 use io,          only:fatal
 type(kdnode), intent(in)    :: node(:)
 integer,      intent(in)    :: ifirstincell(:)
 real,         intent(in)    :: xyzh(:,:)
 integer,      intent(in)    :: nphotosrc
 real,         intent(in)    :: xyz_photosrc(:,:)
 real,         intent(in)    :: tree_accuracy_cmi,rcut_opennode,rcut_leafpart
 integer,      intent(in)    :: nnode_needopen(:)
 integer,      intent(in)    :: nneedopen
 real,         intent(in)    :: nxyzrs_nextopen(:,:)
 integer,      intent(in)    :: nnextopen
 integer,      intent(in)    :: maxcminode,maxleafparts
 integer,      intent(inout) :: ncminode,ncloseleaf,nleafparts
 integer,      intent(inout) :: nnode_toreplace(:)
 real,         intent(inout) :: nxyzm_treetocmi(:,:),ixyzhm_leafparts(:,:)
 real,         intent(out)   :: mintheta
 integer, parameter :: istacksize = 300
 integer :: nstack(istacksize)
 integer :: n,istack,il,ir,isrc,ip,ipnode,inode,n_needopen,iregion
 real    :: dist2,pos_node(3),size_node,mass_node,pos_region(3),rad_region2,minsize_region
 real    :: rcut_opennode2,rcut_leafpart2,rcut_opennode2_inbound
 real    :: size_node2,tree_accuracy_cmi2,theta2,mintheta2
 logical :: open_node,open_node_onesrc,close_to_src

 write(*,'(2x,a60,f5.2)') 'Walking kdtree to extract cmi-nodes with opening criterion: ',tree_accuracy_cmi
 write(*,'(3x,a22,f10.2)') 'radius of leaf nodes: ',rcut_opennode
 write(*,'(3x,a26,f6.2)')  'radius of ind. particles: ',rcut_leafpart

 if (maxlevel > maxlevel_indexed) then
    call fatal('kdtree_cmi','This module only works if maxlevel is within maxlevel_indexed. &
              &Recompile with larger MAXP/NCELLSMAX or use individual particles instead by &
              &setting the variable photoionize_tree (in module photoionize_cmi) to false.')
 endif

 rcut_opennode2 = rcut_opennode**2.
 rcut_leafpart2 = rcut_leafpart**2.
 tree_accuracy_cmi2 = tree_accuracy_cmi**2.

 !- Init counters
 ncminode   = 0  ! number of nodes picked
 ncloseleaf = 0  ! number of leaves to be replaced
 nleafparts = 0  ! number of particles within those leaves

 !- Init params for extracting nodes at boundary of leaf region
 mintheta2 = huge(mintheta2)
 rcut_opennode2_inbound = rcut_opennode2 - 0.01*rcut_opennode2

 istack = 1
 nstack(istack) = irootnode
 over_stack: do while (istack /= 0)
    n = nstack(istack)
    istack = istack - 1
    pos_node   = node(n)%xcen(1:3)
    size_node  = node(n)%size
    size_node2 = size_node**2.
    mass_node  = node(n)%mass
    il = node(n)%leftchild
    ir = node(n)%rightchild
    !
    ! Check if node should be opened based on source locations
    !
    open_node = .false.
    close_to_src = .false.
    over_sources: do isrc = 1,nphotosrc
       dist2 = mag2(pos_node(1:3)-xyz_photosrc(1:3,isrc))
       theta2 = size_node2/dist2
       open_node_onesrc = ( theta2 > tree_accuracy_cmi2 )
       if (open_node_onesrc .or. dist2 < rcut_opennode2) then
          open_node = .true.
          !- extract theta of those at boundary of leaf region
          if (dist2 < rcut_opennode2 .and. dist2 > rcut_opennode2_inbound) then
             mintheta2 = min(mintheta2,theta2)
          endif
          if (dist2 < rcut_leafpart2) then
             close_to_src = .true.
             exit over_sources
          endif
       endif
    enddo over_sources
    !
    ! Check if node has to be opened to resolve into the ionization front
    ! [used within tree-walk iteration]
    !
    if (.not.open_node .and. nneedopen > 0) then
       over_bigcells: do inode = 1,nneedopen
          n_needopen = nnode_needopen(inode)
          if (n == n_needopen) then
             open_node = .true.
             exit over_bigcells
          endif
       enddo over_bigcells
    endif
    !
    ! Check if node should be opened based on iteration outcome from previous step
    ! (i.e. the 'saved' final grid)
    !
    if (.not.open_node .and. nnextopen > 0) then
       over_prevbigcells: do iregion = 1,nnextopen
          pos_region  = nxyzrs_nextopen(2:4,iregion)
          rad_region2 = nxyzrs_nextopen(5,iregion)**2.
          minsize_region = nxyzrs_nextopen(6,iregion)
          dist2 = mag2(pos_node(1:3)-pos_region(1:3))
          if (dist2 < rad_region2) then
             if (size_node > minsize_region) then
                open_node = .true.
                exit over_prevbigcells
             endif
          endif
       enddo over_prevbigcells
    endif
    !
    ! Walk tree to extract nodes and particles within leaves
    !
    if (open_node) then
       is_leaf: if (ifirstincell(n) /= 0) then
          !- Store node
          ncminode = ncminode + 1
          if (ncminode > maxcminode) call fatal('kdtree','ncminode exceeded maxcminode')
          nxyzm_treetocmi(1,ncminode)   = n
          nxyzm_treetocmi(2:4,ncminode) = pos_node(1:3)
          nxyzm_treetocmi(5,ncminode)   = mass_node
          if (close_to_src) then
             !- Store n of leaf [to be removed]
             ncloseleaf = ncloseleaf + 1
             nnode_toreplace(ncloseleaf) = n
             !- Store particles within leaf
             over_parts: do ipnode = inoderange(1,n),inoderange(2,n)
                nleafparts = nleafparts + 1
                if (nleafparts > maxleafparts) call fatal('kdtree','nleafparts exceeded maxleafparts')
                ip = abs(inodeparts(ipnode))  !- include both active and inactive particles
                ixyzhm_leafparts(1,nleafparts)   = ip
                ixyzhm_leafparts(2:5,nleafparts) = xyzh(1:4,ip)
                ixyzhm_leafparts(6,nleafparts)   = massoftype(igas)
             enddo over_parts
          endif
       else ! not leaf
          if (istack+2 > istacksize) call fatal('getneigh','stack overflow in getneigh')
          !- Push children to stack
          if (il /= 0) then
             istack = istack + 1
             nstack(istack) = il
          endif
          if (ir /= 0) then
             istack = istack + 1
             nstack(istack) = ir
          endif
       endif is_leaf

    else ! not open; accept as node
       !- Store node
       ncminode = ncminode + 1
       if (ncminode > maxcminode) call fatal('kdtree','ncminode exceeded maxcminode')
       nxyzm_treetocmi(1,ncminode)   = n
       nxyzm_treetocmi(2:4,ncminode) = pos_node(1:3)
       nxyzm_treetocmi(5,ncminode)   = mass_node
    endif
 enddo over_stack

 mintheta = sqrt(mintheta2)

end subroutine extract_cminodes_from_tree


real function mag2(vec)
 real,   intent(in) :: vec(3)

 mag2 = dot_product(vec,vec)

end function mag2

end module kdtree_cmi
