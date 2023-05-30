!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module hnode_cmi
!
! The CMI suite: photoionize_cmi.f90 kdtree_cmi.f90 *hnode_cmi.f90* heating_cooling_cmi.f90
! This module contains all the subroutines necessary for solving the smoothing lengths
! of cmi-nodes picked from walking the kd-tree.
!
! We pass higher-level nodes (instead of individual particles) to the photoionization
! code CMacIonize to perform density-mapping. For this, we require the xyzhm
! of all nodes. xyzm are collected during tree walk; here, we compute the remaining h.
!
! :References: none
!
! :Owner: Cheryl Lau
!
! :Runtime parameters:
!   - hfact_node : *proportionality fac of hnode to mean local node spacing*
!   - tolh_node  : *tolerance on hnode-iteration*
!   - kmid       : *an arbitrary level midway on the tree to be used in the neighbour-find*
!
! :Dependencies: infile_utils, io, dtypekdtree, kdtree, kernel
!
!
 implicit none

 public :: hnode_iterate

 real,    public :: hfact_node = 1.1
 real,    public :: tolh_node  = 1E-2
 integer, public :: kmid = 10

 private
 integer, parameter :: maxnodeneigh   = 1E7
 integer, parameter :: neighcachesize = 5000
 integer, parameter :: maxkmid = 1E2

contains
!-----------------------------------------------------------------------
!+
!  Solve for h of nodes using their number densities;
!  Since node masses differ, density is no longer proportional to h
!+
!-----------------------------------------------------------------------
subroutine hnode_iterate(node,nxyzm_treetocmi,ncminode,h_solvertocmi)
 use dtypekdtree, only:kdnode
 use io,          only:fatal,warning
 use omp_lib
 type(kdnode), intent(in)    :: node(:)
 integer,      intent(in)    :: ncminode
 real,         intent(in)    :: nxyzm_treetocmi(:,:)
 real,         intent(inout) :: h_solvertocmi(:)
 integer :: listnodeneigh(maxnodeneigh)
 integer :: icminode,na,nnodeneigh,avgneigh
 integer :: maxnode_bruteforce = 1E6
 integer :: n_kmid(maxkmid)
 integer :: nkmid,ierr,i,inosol
 real    :: xyzcache_nodeneigh(3,neighcachesize)
 real    :: avgneigh_est,size_kmid
 real    :: pos_node(3),size_node
 real    :: time_neigh,time1_neigh,time2_neigh,time_iterate,time1_iterate,time2_iterate
 real    :: tottime_neigh,tottime_iterate
 real    :: pos_node_nosol(3,10000)
 logical :: brute_force,node_failed

 write(*,'(2x,a30,i6,a6)') 'Solving smoothing lengths for ',ncminode,' nodes'

 if (ncminode < maxnode_bruteforce) then
    print*,'  brute-force activated'
    brute_force = .true.  !- to avoid overhead in neigh-find
    if (ncminode > maxnodeneigh) call fatal('hnode_cmi','maxnodeneigh has to be greater than ncminode')
 else
    print*,'  brute-force deactivated - switching to neighbour-find algorithm'
    brute_force = .false.
    nkmid = 0
    call init_get_nodeneigh(node,ncminode,nxyzm_treetocmi,n_kmid,nkmid,avgneigh_est,size_kmid)
    if (nkmid == 0) call fatal('hnode_cmi','init_get_nodeneigh failed')
 endif

 node_failed = .false.
 inosol = 0

 avgneigh = 0
 tottime_neigh   = 0.
 tottime_iterate = 0.
 !$omp parallel do default(none) shared(ncminode,nxyzm_treetocmi,node) &
 !$omp shared(hfact_node,h_solvertocmi,tolh_node,brute_force) &
 !$omp shared(n_kmid,nkmid,avgneigh_est,size_kmid) &
 !$omp shared(node_failed,inosol,pos_node_nosol) &
 !$omp private(icminode,na,pos_node,size_node,nnodeneigh,listnodeneigh,xyzcache_nodeneigh) &
 !$omp private(ierr) &
 !$omp private(time_neigh,time1_neigh,time2_neigh,time_iterate,time1_iterate,time2_iterate) &
 !$omp reduction(+:tottime_neigh,tottime_iterate) &
 !$omp reduction(+:avgneigh) &
 !$omp schedule(runtime)
 over_na: do icminode = 1,ncminode
    na = int(nxyzm_treetocmi(1,icminode))
    pos_node  = node(na)%xcen(1:3)
    size_node = node(na)%size

    !- Get trial neighbours nb
    time1_neigh = omp_get_wtime()
    nnodeneigh = 0
    if (brute_force) then
       call get_nodeneigh_bruteforce(node,na,nxyzm_treetocmi,ncminode,nnodeneigh,listnodeneigh,&
                                     xyzcache_nodeneigh)
    else
       call get_nodeneigh(node,na,pos_node,nxyzm_treetocmi,ncminode,nnodeneigh,listnodeneigh,&
                          xyzcache_nodeneigh,n_kmid,nkmid,avgneigh_est,size_kmid)
    endif
    if (nnodeneigh == 0) call fatal('hnode_cmi','no neighbours found')
    time2_neigh = omp_get_wtime()

    !- Solve for its smoothing length
    time1_iterate = omp_get_wtime()
    call hnode_newton_raphson(node,icminode,pos_node,size_node,ncminode,nnodeneigh,listnodeneigh,&
                              xyzcache_nodeneigh,avgneigh,h_solvertocmi,ierr)
    if (ierr /= 0) then
       write(*,'(3x,a30,i7,a25)') 'Newton-Raphson failed for node',na,', trying bisection method'
       call hnode_bisection(node,icminode,pos_node,size_node,ncminode,nnodeneigh,listnodeneigh,&
                            xyzcache_nodeneigh,avgneigh,h_solvertocmi,ierr)
    endif

    !- If h completely failed to converge, simply do h = hfact*node-separation
    !  (acceptable since the h of nodes only affect the resolution in mapping but not the physics)
    if (ierr /= 0) then
       !$omp critical
       inosol = inosol + 1
       pos_node_nosol(1:3,inosol) = pos_node(1:3)
       node_failed = .true.
       !$omp end critical
       h_solvertocmi(icminode) = hfact_node*2.*size_node
       write(*,'(3x,a30,i7,a31)') 'Bisection also failed for node',na,', using node size to estimate h'
    endif

    time2_iterate = omp_get_wtime()

    time_neigh = time2_neigh - time1_neigh
    time_iterate = time2_iterate - time1_iterate
    tottime_neigh = tottime_neigh + time_neigh
    tottime_iterate = tottime_iterate + time_iterate

 enddo over_na
 !$omp end parallel do

 if (node_failed) then
    open(2052,file='nodes_no_sol.txt',status='replace')
    do i = 1,inosol
       write(2052,*) pos_node_nosol(1:3,i)
    enddo
    close(2052)
    if (inosol >= 10) call warning('hnode_cmi','Many nodes have smoothing lengths not converging')
 endif

 avgneigh = nint(real(avgneigh)/real(ncminode))
 write(*,'(3x,a41,i2)') '- Average number of neighbours per node: ',avgneigh

 write(*,'(3x,a18,es18.6)') '- neigh-find time ',tottime_neigh
 write(*,'(3x,a15,es21.6)') '- iterate time ',tottime_iterate

end subroutine hnode_iterate

!
! Solve for h with Newton-Raphson method
!
subroutine hnode_newton_raphson(node,icminode,pos_na,size_na,ncminode,nnodeneigh,listnodeneigh,&
                                xyzcache_nodeneigh,avgneigh,h_solvertocmi,ierr)
 use dtypekdtree, only:kdnode
 type(kdnode), intent(in)    :: node(:)
 integer,      intent(in)    :: icminode,ncminode,nnodeneigh
 integer,      intent(in)    :: listnodeneigh(:)
 real,         intent(in)    :: xyzcache_nodeneigh(:,:)
 real,         intent(in)    :: pos_na(3),size_na
 integer,      intent(inout) :: avgneigh
 real,         intent(out)   :: h_solvertocmi(:)
 integer,      intent(out)   :: ierr
 integer :: nrealneigh,niter
 real    :: h,h0,h_old,func,dfunc
 logical :: converged

 converged = .false.
 ierr  = 0
 niter = 0
 h0 = hfact_node*2.*size_na   !- Init with estimated value of h
 h  = h0

 do while (.not.converged)
    h_old = h
    call getfunc_from_neigh_sum(h,pos_na,node,nnodeneigh,listnodeneigh,xyzcache_nodeneigh,&
                                nrealneigh,func,dfunc)
    !- Update h
    h = h - func/dfunc

    if (abs(h-h_old)/h0 < tolh_node .and. h > 0.) then
       converged = .true.
       h_solvertocmi(icminode) = h
       avgneigh = avgneigh + nrealneigh
    endif
    niter = niter + 1
    if (niter > 100) exit
 enddo
 if (.not.converged) ierr = 1

end subroutine hnode_newton_raphson

!
! Solve for h with bisection method
!
subroutine hnode_bisection(node,icminode,pos_na,size_na,ncminode,nnodeneigh,listnodeneigh,&
                           xyzcache_nodeneigh,avgneigh,h_solvertocmi,ierr)
 use dtypekdtree, only:kdnode
 use io,          only:fatal
 type(kdnode), intent(in)    :: node(:)
 integer,      intent(in)    :: icminode,ncminode,nnodeneigh
 integer,      intent(in)    :: listnodeneigh(:)
 real,         intent(in)    :: xyzcache_nodeneigh(:,:)
 real,         intent(in)    :: pos_na(3),size_na
 integer,      intent(inout) :: avgneigh
 real,         intent(out)   :: h_solvertocmi(:)
 integer,      intent(out)   :: ierr
 integer :: nrealneigh,niter
 real    :: hmin,hmax,hmid,func_mid,func_min,func_max,dfunc
 logical :: converged

 converged = .false.
 ierr  = 0
 niter = 0
 hmin  = hfact_node*1E-2*size_na   !- Init search range
 hmax  = hfact_node*1E+2*size_na

 do while (.not.converged)
    !- Eval func at hmin
    call getfunc_from_neigh_sum(hmin,pos_na,node,nnodeneigh,listnodeneigh,xyzcache_nodeneigh,&
                                nrealneigh,func_min,dfunc)  !- dfunc not used here
    !- Eval func at hmax
    call getfunc_from_neigh_sum(hmax,pos_na,node,nnodeneigh,listnodeneigh,xyzcache_nodeneigh,&
                                nrealneigh,func_max,dfunc)
    if (func_max*func_min > 0) exit

    !- Update h and eval func at h
    hmid = (hmax+hmin)/2.
    call getfunc_from_neigh_sum(hmid,pos_na,node,nnodeneigh,listnodeneigh,xyzcache_nodeneigh,&
                                nrealneigh,func_mid,dfunc)

    if (func_mid*func_min < 0) then
       hmax = hmid
    elseif (func_mid*func_max < 0) then
       hmin = hmid
    else
       call fatal('hnode_cmi','an error occurred in bisection root-finding')
    endif

    if (abs(hmax-hmin)/hmid < tolh_node) then
       hmid = (hmax+hmin)/2.
       if (hmid > 0.) then
          converged = .true.
          h_solvertocmi(icminode) = hmid
          avgneigh = avgneigh + nrealneigh
       endif
    endif
    niter = niter + 1
    if (niter > 100) exit
 enddo
 if (.not.converged) ierr = 1

end subroutine hnode_bisection

!
! Sum over neighbours and compute func for a given h
!
subroutine getfunc_from_neigh_sum(h,pos_node,node,nnodeneigh,listnodeneigh,xyzcache_nodeneigh,&
                                  nrealneigh,func,dfunc)
 use dtypekdtree, only:kdnode
 use kernel,      only:get_kernel,cnormk,radkern2
 type(kdnode), intent(in)  :: node(:)
 integer,      intent(in)  :: nnodeneigh
 integer,      intent(in)  :: listnodeneigh(:)
 real,         intent(in)  :: xyzcache_nodeneigh(:,:)
 real,         intent(in)  :: h
 real,         intent(in)  :: pos_node(3)
 integer,      intent(out) :: nrealneigh
 real,         intent(out) :: func,dfunc
 integer :: ineigh,nb
 real    :: wab,dwdh,wab_sum,dwdh_sum
 real    :: h1,dr2,q,q2,wkern,grkern,pos_neigh(3)

 wab_sum  = 0
 dwdh_sum = 0
 nrealneigh = 0
 h1 = 1./h
 over_neigh: do ineigh = 1,nnodeneigh
    if (ineigh <= neighcachesize) then
       pos_neigh = xyzcache_nodeneigh(1:3,ineigh)
    else
       nb = listnodeneigh(ineigh)
       pos_neigh = node(nb)%xcen(1:3)
    endif
    dr2 = mag2(pos_node(1:3)-pos_neigh(1:3))
    q2 = dr2 * h1**2
    if (q2 < radkern2) then   !- r/h < compact support size
       q = sqrt(q2)
       call get_kernel(q2,q,wkern,grkern)
       wab  = wkern
       dwdh = q*grkern + 3.*wkern
       wab_sum  = wab_sum + wab
       dwdh_sum = dwdh_sum + dwdh
       nrealneigh = nrealneigh + 1
    endif
 enddo over_neigh
 wab_sum  = cnormk*h1**3 * wab_sum
 dwdh_sum = -cnormk*h1**4 * dwdh_sum

 !- Solves the eqn h = hfact*n^(-1/3) where n = sum W_ab(h)
 func  = hfact_node*(wab_sum)**(-1./3.) - h
 dfunc = (-1./3.*hfact_node*(wab_sum)**(-4./3.) * dwdh_sum) -1

end subroutine getfunc_from_neigh_sum

!-----------------------------------------------------------------------
!+
! Get neighbours with brute-force
!+
!-----------------------------------------------------------------------
subroutine get_nodeneigh_bruteforce(node,na,nxyzm_treetocmi,ncminode,nnodeneigh,listnodeneigh,&
                                    xyzcache_nodeneigh)
 use dtypekdtree, only:kdnode
 use io,          only:fatal
 type(kdnode), intent(in)    :: node(:)
 integer,      intent(in)    :: na  !- target node
 real,         intent(in)    :: nxyzm_treetocmi(:,:)
 integer,      intent(in)    :: ncminode
 integer,      intent(inout) :: nnodeneigh
 integer,      intent(out)   :: listnodeneigh(:)
 real,         intent(out)   :: xyzcache_nodeneigh(:,:)
 integer :: inode,nb

 do inode = 1,ncminode
    nb = int(nxyzm_treetocmi(1,inode))
    if (nb /= na) then
       nnodeneigh = nnodeneigh + 1
       if (nnodeneigh > maxnodeneigh) call fatal('hnode_cmi','nnodeneigh exceeded maxnodeneigh')
       if (nnodeneigh <= neighcachesize) then
          listnodeneigh(nnodeneigh) = nb
          xyzcache_nodeneigh(1:3,nnodeneigh) = node(nb)%xcen(1:3)
       else
          listnodeneigh(nnodeneigh) = nb
       endif
    endif
 enddo

end subroutine get_nodeneigh_bruteforce

!-----------------------------------------------------------------------
!+
! Get trial neighbours for a given node na using their indicies
!+
!-----------------------------------------------------------------------
subroutine init_get_nodeneigh(node,ncminode,nxyzm_treetocmi,n_kmid,nkmid,avgneigh_est,size_kmid)
 use dtypekdtree, only:kdnode
 use physcon,     only:pi
 use io,          only:fatal
 type(kdnode), intent(in)    :: node(:)
 integer,      intent(in)    :: ncminode
 real,         intent(in)    :: nxyzm_treetocmi(:,:)
 integer,      intent(inout) :: nkmid
 real,         intent(inout) :: avgneigh_est,size_kmid
 integer,      intent(out)   :: n_kmid(:)
 integer :: max_n,inode,kmid_old,kmid_start,kmid_end,n
 logical :: kmid_ok

 !- Make sure kmid is above all nodes
 max_n = 1
 do inode = 1,ncminode
    max_n = max(max_n,int(nxyzm_treetocmi(1,inode)))
 enddo
 kmid_old = kmid
 kmid_ok = .false.
 do while (.not.kmid_ok)
    kmid_ok = ( 2**kmid+(2**(kmid+1)-1) < max_n/2. )
    if (.not.kmid_ok) kmid = kmid - 1
 enddo
 if (kmid <= 0) call fatal('hnode_cmi','invalid kmid')
 if (kmid /= kmid_old) print*,'changing kmid to ',kmid

 !- Get all nodes on level kmid and its avg node size
 kmid_start = int(2**kmid)
 kmid_end   = int(2**(kmid+1)-1)
 nkmid = 0
 size_kmid = 0.
 do n = kmid_start,kmid_end
    nkmid = nkmid + 1
    if (nkmid > maxkmid) call fatal('hnode_cmi','require a larger maxkmid')
    n_kmid(nkmid) = n
    size_kmid = size_kmid + node(n)%size
 enddo
 size_kmid = size_kmid/nkmid

 !- Estimate number of neighbours
 avgneigh_est = 10.67*pi*hfact_node**3.

end subroutine init_get_nodeneigh


subroutine get_nodeneigh(node,na,pos_na,nxyzm_treetocmi,ncminode,nnodeneigh,listnodeneigh,&
                         xyzcache_nodeneigh,n_kmid,nkmid,avgneigh_est,size_kmid)
 use dtypekdtree, only:kdnode
 use io,          only:fatal
 type(kdnode), intent(in)    :: node(:)
 integer,      intent(in)    :: na  !- target node
 real,         intent(in)    :: pos_na(3)
 real,         intent(in)    :: nxyzm_treetocmi(:,:)
 integer,      intent(in)    :: ncminode,nkmid
 integer,      intent(in)    :: n_kmid(maxkmid)
 real,         intent(in)    :: avgneigh_est,size_kmid
 integer,      intent(inout) :: nnodeneigh
 integer,      intent(out)   :: listnodeneigh(:)
 real,         intent(out)   :: xyzcache_nodeneigh(:,:)
 integer :: ka,icminode,ikmid,n,nb,kb,deltak
 integer :: na_neigh_above,kabove,kabove_start,kabove_end
 real    :: rcut,rcut2,size_neigh,pos_n(3),dist2
 logical :: close_to_na(maxkmid)
 !
 ! Define rcut of na
 !
 ka = int(log(real(na))/log(2.))
 kabove = ka - 1  ! <- can use even higher levels in case there're 'jumps' in tree-walk
 kabove_start = int(2**kabove)
 kabove_end   = int(2**(kabove+1)-1)
 rcut = epsilon(rcut)
 do na_neigh_above = kabove_start,kabove_end
    size_neigh = node(na_neigh_above)%size
    rcut = max(rcut,size_neigh)
 enddo
 rcut  = rcut*avgneigh_est**(1./3.) + size_kmid
 rcut2 = rcut*rcut
 !
 ! Compute distance of na to nodes on kmid to determine in or out
 !
 do ikmid = 1,maxkmid
    if (ikmid <= nkmid) then
       n = n_kmid(ikmid)
       pos_n = node(n)%xcen(1:3)
       dist2 = mag2(pos_n(1:3)-pos_na(1:3))
       if (dist2 < rcut2) then
          close_to_na(ikmid) = .true.
       else
          close_to_na(ikmid) = .false.
       endif
    else
       close_to_na(ikmid) = .false.
    endif
 enddo
 !
 ! Loop over all other cmi-nodes to pick the right nb
 !
 over_nb: do icminode = 1,ncminode
    nb = int(nxyzm_treetocmi(1,icminode))
    kb = int(log(real(nb))/log(2.))
    deltak = kb - kmid
    if (deltak < 0) call fatal('hnode_cmi','kmid should be above kb')
    do ikmid = 1,nkmid
       if (int(real(nb)/real(2**deltak)) == n_kmid(ikmid) .and. nb /= na) then
          if (close_to_na(ikmid)) then
             nnodeneigh = nnodeneigh + 1
             if (nnodeneigh <= neighcachesize) then
                listnodeneigh(nnodeneigh) = nb
                xyzcache_nodeneigh(1:3,nnodeneigh) = node(nb)%xcen(1:3)
             else
                listnodeneigh(nnodeneigh) = nb
             endif
          endif
       endif
    enddo
 enddo over_nb

end subroutine get_nodeneigh


real function mag2(vec)
 real,   intent(in) :: vec(3)

 mag2 = dot_product(vec,vec)

end function mag2

end module hnode_cmi
