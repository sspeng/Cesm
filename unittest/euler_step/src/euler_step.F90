program main
use kinds , only: real_kind
use element_mod, only: element_t
use dimensions_mod, only: np, npdg, nlev, qsize
use hybrid_mod, only: hybrid_t
use derivative_mod, only: derivative_t
use hybvcoord_mod, only:hvcoord_t
implicit none
  integer :: np1_qdp = 1
  integer :: n0_qdp = 2
  real(kind=real_kind) :: dt = 150.0
  type(element_t) :: elem(43)
  type(element_t) :: elem_test(43)
  type(hybrid_t) :: hybrid
  type(derivative_t) :: deriv
  type(hvcoord_t) :: hvcoord
  integer, parameter :: nets = 1
  integer, parameter :: nete = 43
  integer :: DSSopt = 3
  integer :: rhs_multiplier = 0

  real(kind=real_kind), dimension(np, np                        ) :: divdp, dpdiss
  real(kind=real_kind), dimension(np, np                        ) :: div
  real(kind=real_kind), dimension(np, np, nlev                  ) :: dpdissk
  real(kind=real_kind), dimension(np, np, 2                     ) :: gradQ
  real(kind=real_kind), dimension(np, np, 2, nlev               ) :: Vstar
  real(kind=real_kind), dimension(np, np, nlev                  ) :: Qtens
  real(kind=real_kind), dimension(np, np, nlev                  ) :: dp, dp_star
  real(kind=real_kind), dimension(np, np, nlev, qsize, nets:nete) :: Qtens_biharmonic
  real(kind=real_kind), pointer, dimension(:,:,:)                 :: DSSvar
  real(kind=real_kind) :: dp0(nlev), qim_val(nlev), qmax_val(nlev)
  integer :: ie, q, i, j, k, d, l
  integer :: rhs_viss = 0
  integer :: qbeg, qend, kbeg, kend
  integer :: kptr

  interface
    subroutine euler_step( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
      use kinds , only: real_kind
      use element_mod, only: element_t
      use hybvcoord_mod, only:hvcoord_t
      use hybrid_mod, only: hybrid_t
      use derivative_mod, only: derivative_t
      integer :: np1_qpd, n0_qdp, nets, nete, DSSopt, rhs_multiplier
      real(kind=real_kind) :: dt
      type(element_t) :: elem(:)
      type(hybrid_t) :: hybrid
      type(derivative_t) :: deriv
      type(hvcoord_t) :: hvcoord
    end subroutine euler_step
  end interface

  call athread_init()

  do ie=nets,nete
    do q=1,qsize
      do k=1,nlev
        do j=1,np
          do i=1,np
            elem(ie)%state%Qdp(i,j,k,q,n0_qdp) = 1!ie * 1000000 + k * 100 + j * 10 + i + 0.2
            elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = 2!ie * 1000000 + k * 100 + j * 10 + i + 0.1
          enddo
        enddo
      enddo
    enddo
  enddo

  do ie = nets, nete
    do k = 1, nlev
      do j = 1, np
        do i = 1, np
          elem(ie)%derived%dp(i,j,k) = 1 !(ie * 1000000 + k * 100 + j * 10 + i)*2
          elem(ie)%derived%divdp_proj(i,j,k) = 1 ! ie * 1000000 + k * 100 + j * 10 + i
          elem(ie)%derived%divdp(i,j,k) = 1 !ie * 1000000 + k * 100 + j * 10 + i
          elem(ie)%derived%dpdiss_biharmonic(i,j,k) = 1! ie * 1000000 + k * 100 + j * 10 + i
        enddo
      enddo
    enddo
  enddo

  do ie = nets, nete
    do k = 1, nlev
      do j = 1, np
        do i = 1, np
          elem(ie)%derived%vn0(i,j,1,k) = 1!(ie * 1000000 + k * 100 + j * 10 + i)
          elem(ie)%derived%vn0(i,j,2,k) = 1!(ie * 1000000 + k * 100 + j * 10 + i)*2
        enddo
      enddo
    enddo
  enddo

  do l = 1, np
    do i = 1, np
      deriv%Dvv(i, l) = 10*l + i
    enddo
  enddo

  do ie = nets, nete
    do j = 1, np
      do i = 1, np
        elem(ie)%Dinv(i, j, 1, 1) = 1
        elem(ie)%Dinv(i, j, 2, 1) = 2
        elem(ie)%Dinv(i, j, 1, 2) = 3
        elem(ie)%Dinv(i, j, 2, 2) = 4
      enddo
    enddo
  enddo

  do ie = nets, nete
    do j = 1, np
      do i = 1, np
        elem(ie)%metdet(i, j) = 1
        elem(ie)%rmetdet(i, j) = 1
      enddo
    enddo
  enddo

  do ie = nets, nete
    do j = 1, np
      do i = 1, np
        elem(ie)%variable_hyperviscosity(i,j) = 1
        elem(ie)%tensorVisc(i,j,1,1) = 1
        elem(ie)%tensorVisc(i,j,1,2) = 1
        elem(ie)%tensorVisc(i,j,1,2) = 1
        elem(ie)%tensorVisc(i,j,2,2) = 1
      enddo
    enddo
  enddo

  call euler_step( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
end program main

subroutine euler_step( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
use kinds , only: real_kind
use physical_constants, only: rrearth
use element_mod, only: element_t
use dimensions_mod, only: np, npdg, nlev, qsize, max_corner_elem, max_neigh_edges, nelemd
use hybvcoord_mod, only:hvcoord_t
use hybrid_mod, only: hybrid_t
use derivative_mod, only: derivative_t
use edgetype_mod, only       : EdgeBuffer_t
implicit none
  integer, intent(in) :: np1_qdp, n0_qdp
  real(kind=real_kind) :: dt
  type(element_t), intent(inout) :: elem(:)
  type(hybrid_t), intent(in) :: hybrid
  type(derivative_t), intent(in) :: deriv
  type(hvcoord_t), intent(in) :: hvcoord
  integer, intent(in):: nets
  integer, intent(in) :: nete
  integer :: DSSopt
  integer :: rhs_multiplier
  real(kind=real_kind), dimension(np, np                           ) :: divdp, dpdiss
  real(kind=real_kind), dimension(np, np                           ) :: div
  real(kind=real_kind), dimension(np, np,    nlev                  ) :: dpdissk
  real(kind=real_kind), dimension(np, np, 2                        ) :: gradQ
  real(kind=real_kind), dimension(np, np, 2, nlev                  ) :: Vstar
  real(kind=real_kind), dimension(np, np,    nlev                  ) :: Qtens
  real(kind=real_kind), dimension(np, np,    nlev                  ) :: dp, dp_star
  real(kind=real_kind), dimension(np, np,    nlev, qsize, nets:nete) :: Qtens_biharmonic
  real(kind=real_kind), dimension(np, np,    nlev,        nets:nete) :: dp_temp
  real(kind=real_kind), dimension(np, np,    nlev, qsize, nets:nete) :: dp_star_temp
  real(kind=real_kind), dimension(np, np, 2, nlev, qsize, nets:nete) :: gradQ_temp
  real(kind=real_kind), dimension(           nlev, qsize, nets:nete) :: qmax
  real(kind=real_kind), dimension(           nlev, qsize, nets:nete) :: qmin
  real(kind=real_kind), dimension(np, np,    nlev, qsize, nets:nete) :: Qtens_temp
  real(kind=real_kind), pointer, dimension(:,:,:) :: DSSvar
  real(kind=real_kind) :: dp0(nlev), qim_val(nlev), qmax_val(nlev)
  type (EdgeBuffer_t) :: edgeAdv, edgeAdvp1, edgeAdvQminmax, edgeAdv1,  edgeveloc
  integer :: ie, q, i, j, k, d
  integer :: limiter_option = 8
  integer :: rhs_viss = 0
  integer :: qbeg, qend, kbeg, kend
  integer :: kptr
  real(kind=real_kind) :: nu_p = 0.0D5
  real(kind=real_kind) :: nu_q = -1
  type(element_t) :: elem_test(43)

  interface
    subroutine biharmonic_wk_scalar(elem, Qtens_biharmonic, deriv, edgeAdv, hybrid, nets, nete)
      use kinds , only: real_kind
      use element_mod, only: element_t
      use hybrid_mod, only: hybrid_t
      use derivative_mod, only: derivative_t
      use edgetype_mod, only       : EdgeBuffer_t
      type(element_t) :: elem(:)
      real(kind=real_kind) :: Qtens_biharmonic(:,:,:,:,:)
      type(EdgeBuffer_t) :: edgeAdv
      type(derivative_t) :: deriv
      type(hybrid_t) :: hybrid
      integer :: nets, nete
    end subroutine biharmonic_wk_scalar
  end interface
!#define QMAX
!#ifdef QMAX
!  external :: slave_euler_step
!  type param_t
!    integer*8 :: qdp_s_ptr, qdp_leap_ptr, dp_s_ptr, divdp_proj_s_ptr   &
!        , Qtens_biharmonic, qmax, qmin
!    real(kind=real_kind) :: dt
!    integer :: nets, nete, np1_qdp, n0_qdp, DSSopt, rhs_multiplier, qsize
!  end type param_t
!  type(param_t) :: param_s
!  param_s%qdp_s_ptr = loc(elem(nets)%state%Qdp(:,:,:,:,:))
!  param_s%qdp_leap_ptr = loc(elem((nets+1))%state%Qdp(:,:,:,:,:))
!  param_s%dp_s_ptr = loc(elem(nets)%derived%dp(:,:,:))
!  param_s%divdp_proj_s_ptr = loc(elem(nets)%derived%divdp_proj(:,:,:))
!  param_s%Qtens_biharmonic = loc(Qtens_biharmonic)
!  param_s%qmax = loc(qmax)
!  param_s%qmin = loc(qmin)
!  param_s%dt = dt
!  param_s%nets = nets
!  param_s%nete = nete
!  param_s%np1_qdp = np1_qdp
!  param_s%n0_qdp = n0_qdp
!  param_s%DSSopt = DSSopt
!  param_s%rhs_multiplier = rhs_multiplier
!  param_s%qsize = qsize
!  call athread_init()
!  !call athread_spawn(slave_euler_step, param_s)
!  !call athread_join()
!#else
  external :: slave_euler_v
  type param_2d_t
    integer*8 :: qdp_s_ptr, qdp_leap_ptr, divdp_proj, dp, vn0, Dvv, Dinv       \
    , metdet, rmetdet, Qtens_biharmonic, divdp, dpdiss_biharmonic, spheremp    \
    , qmax, qmin, Qtens_temp, dp_star_temp, dp_temp
    real(kind=real_kind) :: dt, rrearth, nu_p, nu_q
    integer :: nets, nete, rhs_multiplier, qsize, n0_qdp, np1_qdp, limiter_option \
        , rhs_viss
  end type param_2d_t
  type(param_2d_t) :: param_2d_s
  param_2d_s%qdp_s_ptr = loc(elem(nets)%state%Qdp(:,:,:,:,:))
  param_2d_s%qdp_leap_ptr = loc(elem((nets+1))%state%Qdp(:,:,:,:,:))
  param_2d_s%divdp_proj = loc(elem(nets)%derived%divdp_proj(:,:,:))
  param_2d_s%dp = loc(elem(nets)%derived%dp(:,:,:))
  param_2d_s%vn0 = loc(elem(nets)%derived%vn0(:,:,:,:))
  param_2d_s%Dvv = loc(deriv%Dvv)
  param_2d_s%Dinv = loc(elem(nets)%Dinv(:,:,:,:))
  param_2d_s%metdet = loc(elem(nets)%metdet(:,:))
  param_2d_s%rmetdet = loc(elem(nets)%rmetdet(:,:))
  param_2d_s%Qtens_biharmonic = loc(Qtens_biharmonic)
  param_2d_s%divdp = loc(elem(nets)%derived%divdp)
  param_2d_s%dpdiss_biharmonic = loc(elem(nets)%derived%dpdiss_biharmonic)
  param_2d_s%spheremp = loc(elem(nets)%spheremp)
  param_2d_s%qmax = loc(qmax)
  param_2d_s%qmin = loc(qmin)
  param_2d_s%Qtens_temp = loc(Qtens_temp(1,1,1,1,nets))
  param_2d_s%dp_star_temp = loc(dp_star_temp(1,1,1,1,nets))
  param_2d_s%dp_temp = loc(dp_temp(1,1,1,nets))
  param_2d_s%dt = dt
  param_2d_s%rrearth = rrearth
  param_2d_s%nu_p = nu_p
  param_2d_s%nu_q = nu_q
  param_2d_s%nets = nets
  param_2d_s%nete = nete
  param_2d_s%rhs_multiplier = rhs_multiplier
  param_2d_s%qsize = qsize
  param_2d_s%n0_qdp = n0_qdp
  param_2d_s%np1_qdp = np1_qdp
  param_2d_s%limiter_option = limiter_option
  param_2d_s%rhs_viss = rhs_viss
  call athread_spawn(slave_euler_v, param_2d_s)
  call athread_join()
!#endif

  allocate(edgeAdv%putmap(max_neigh_edges,nelemd))
  allocate(edgeAdv%getmap(max_neigh_edges,nelemd))
  allocate(edgeAdv%reverse(max_neigh_edges,nelemd))
  allocate(edgeAdv%buf(2000))
  allocate(edgeAdv%receive(2000))
  do ie = nets, nete
    do i = 1, max_neigh_edges
      edgeAdv%putmap(i, ie) = i;
      edgeAdv%getmap(i, ie) = i;
      edgeAdv%reverse(i, ie) = .true.
    enddo
  enddo
  !call biharmonic_wk_scalar(elem, Qtens_biharmonic, deriv, edgeAdv, hybrid, nets, nete)
#if 1
#define PRINT_QTEN
#ifdef PRINT_QTEN
  do ie=nets,nete
    do q=1,qsize
      do k=1,nlev
        do j=1,np
          do i=1,np
            print *, dp_star_temp(i,j,k,q,ie)
          enddo
        enddo
      enddo
    enddo
  enddo
#else
  do ie=nets,nete
    do q=1,qsize
      do k=1, nlev
        print *, qmax(k,q,ie), qmin(k,q,ie)
      enddo
    enddo
  enddo

  qbeg = 1
  qend = qsize
  kbeg = 1
  kend = nlev
  rhs_viss = 0
#endif
#endif

  deallocate(edgeAdv%putmap)
  deallocate(edgeAdv%getmap)
  deallocate(edgeAdv%reverse)
  deallocate(edgeAdv%buf)
  deallocate(edgeAdv%receive)
end subroutine euler_step

subroutine biharmonic_wk_scalar(elem, qtens, deriv, edgeq, hybrid, nets, nete)
use kinds , only: real_kind
use physical_constants, only: rrearth
use element_mod, only: element_t
use control_mod, only : north, south, east, west, neast, nwest, seast, swest
use dimensions_mod, only: np, npdg, nlev, qsize, max_corner_elem, max_neigh_edges, nelemd
use hybvcoord_mod, only:hvcoord_t
use hybrid_mod, only: hybrid_t
use derivative_mod, only: derivative_t
use edgetype_mod, only       : EdgeBuffer_t
use control_mod, only : hypervis_scaling, hypervis_power
implicit none

type(hybrid_t),  intent(in) :: hybrid
type(element_t), intent(inout), target :: elem(:)
integer :: nets, nete
real(kind=real_kind), dimension(np, np, nlev, qsize, nets:nete) :: Qtens
type(EdgeBuffer_t),  intent(inout) :: edgeq
type(derivative_t), intent(in) :: deriv


! local
!integer :: k, kptr, i, j, ie, ic, q
!real(kind=real_kind), dimension(np,np) :: lap_p
!integer :: reverse(4, nets:nete)
!integer :: step_elem
!external :: slave_biharmonic_wk_scalar
!type param_t
!  integer*8 :: Dinv, rspheremp, spheremp, variable_hyperviscosity, tensorVisc, Qtens     \
!      , Dvv, buf
!  integer*8 :: putmap, reverse
!  real(kind=real_kind) :: rrearth, hypervis_scaling, hypervis_power
!  integer :: nets, nete, qsize, step_elem
!end type param_t
!type(param_t) :: param_s
!do ie = nets, nete
!  if (edgeq%reverse(west, ie)) then
!    reverse(1, ie) = 1
!  else
!    reverse(1, ie) = 0
!  endif
!  if (edgeq%reverse(east, ie)) then
!    reverse(2, ie) = 1
!  else
!    reverse(2, ie) = 0
!  endif
!  if (edgeq%reverse(south, ie)) then
!    reverse(3, ie) = 1
!  else
!    reverse(3, ie) = 0
!  endif
!  if (edgeq%reverse(north, ie)) then
!    reverse(4, ie) = 1
!  else
!    reverse(4, ie) = 0
!  endif
!enddo
!
!param_s%Dinv = loc(elem(nets)%Dinv)
!param_s%rspheremp = loc(elem(nets)%rspheremp)
!param_s%spheremp = loc(elem(nets)%spheremp)
!param_s%variable_hyperviscosity = loc(elem(nets)%variable_hyperviscosity)
!param_s%tensorVisc = loc(elem(nets)%tensorVisc)
!param_s%Qtens = loc(Qtens)
!param_s%Dvv = loc(deriv%Dvv)
!param_s%buf = loc(edgeq%buf)
!param_s%putmap = loc(edgeq%putmap(1,nets))
!param_s%reverse = loc(reverse)
!param_s%rrearth = rrearth
!param_s%hypervis_scaling = hypervis_scaling
!param_s%hypervis_power = hypervis_power
!param_s%nets = nets
!param_s%nete = nete
!param_s%qsize = qsize
!param_s%step_elem = (loc(elem(nets+1)%Dinv) - loc(elem(nets)%Dinv))/8
!call athread_spawn(slave_biharmonic_wk_scalar, param_s)
!call athread_join()

integer :: k, kptr, i, j, ie, ic, q

external :: slave_biharmonic_wk_scalar_2
type param_2_t
  integer*8 :: Dinv, rspheremp, spheremp, variable_hyperviscosity, tensorVisc, Qtens     \
      , Dvv, receive
  integer*8 :: getmap
  real(kind=real_kind) :: rrearth, hypervis_scaling, hypervis_power
  integer :: nets, nete, qsize, step_elem
end type param_2_t
type(param_2_t) :: param_s2

param_s2%Dinv = loc(elem(nets)%Dinv)
param_s2%rspheremp = loc(elem(nets)%rspheremp)
param_s2%spheremp = loc(elem(nets)%spheremp)
param_s2%variable_hyperviscosity = loc(elem(nets)%variable_hyperviscosity)
param_s2%tensorVisc = loc(elem(nets)%tensorVisc)
param_s2%Qtens = loc(Qtens)
param_s2%Dvv = loc(deriv%Dvv)
param_s2%receive = loc(edgeq%receive)
param_s2%getmap = loc(edgeq%getmap(1,nets))
param_s2%rrearth = rrearth
param_s2%hypervis_scaling = hypervis_scaling
param_s2%hypervis_power = hypervis_power
param_s2%nets = nets
param_s2%nete = nete
param_s2%qsize = qsize
param_s2%step_elem = (loc(elem(nets+1)%Dinv) - loc(elem(nets)%Dinv))/8
call athread_spawn(slave_biharmonic_wk_scalar_2, param_s2)
call athread_join()


#if 0
  do ie = nets, 1
    do q = 1, qsize
      do k = 1, nlev
        do j = 1, np
          do i = 1, np
            !print *, elem(ie)%variable_hyperviscosity(i,j), elem(ie)%tensorVisc(i,j,1,1) \
            !    , elem(ie)%tensorVisc(i,j,1,2), elem(ie)%tensorVisc(i,j,1,2)       \
            !    , elem(ie)%tensorVisc(i,j,2,2)
            print *, Qtens(i,j,k,q,ie)
          enddo
        enddo
      enddo
    enddo
  enddo

#endif


end subroutine biharmonic_wk_scalar
