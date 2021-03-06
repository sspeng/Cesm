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
  type(element_t) :: elem(42)
  type(element_t) :: elem_test(42)
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
  integer :: ie, q, i, j, k
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

  do ie=nets,nete
    do q=1,qsize
      do k=1,nlev
        do j=1,np
          do i=1,np
            elem(ie)%state%Qdp(i,j,k,q,n0_qdp) = ie * 1000000 + q * 10000 + k * 100 + j * 10 + i + 0.2
            elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = ie * 1000000 + q * 10000 + k * 100 + j * 10 + i + 0.1
          enddo
        enddo
      enddo
    enddo
  enddo

  do ie = nets, nete
    do k = 1, nlev
      do j = 1, np
        do i = 1, np
          elem(ie)%derived%dp(i,j,k) = 8.0
          elem(ie)%derived%divdp_proj(i,j,k) = 2.0
        enddo
      enddo
    enddo
  enddo


  call euler_step( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
end program main

subroutine euler_step( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
use kinds , only: real_kind
use element_mod, only: element_t
use dimensions_mod, only: np, npdg, nlev, qsize
use hybvcoord_mod, only:hvcoord_t
use hybrid_mod, only: hybrid_t
use derivative_mod, only: derivative_t
implicit none
  integer, intent(in) :: np1_qdp, n0_qdp
  real(kind=real_kind) :: dt
  type(element_t), intent(inout), target :: elem(:)
  type(hybrid_t), intent(in) :: hybrid
  type(derivative_t), intent(in) :: deriv
  type(hvcoord_t), intent(in) :: hvcoord
  integer, intent(in):: nets
  integer, intent(in) :: nete
  integer :: DSSopt
  integer :: rhs_multiplier
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
  integer :: ie, q, i, j, k
  integer :: rhs_viss = 0
  integer :: qbeg, qend, kbeg, kend
  integer :: kptr
  type(element_t) :: elem_test(42)

  external :: slave_euler_step
  type param_t
    integer*8 :: qdp_s_ptr, qdp_leap_ptr,dp_s_ptr, dp_leap_ptr, divdp_proj_s_ptr, divdp_proj_leap_ptr, qdp_test_ptr
    real(kind=real_kind) :: dt
    integer :: nets, nete, np1_qdp, n0_qdp, DSSopt, rhs_multiplier, qsize
  end type param_t
  type(param_t) :: param_s
  param_s%qdp_s_ptr = loc(elem(1)%state%Qdp(:,:,:,:,:))
  param_s%qdp_leap_ptr = loc(elem(2)%state%Qdp(:,:,:,:,:))
  param_s%dp_s_ptr = loc(elem(1)%derived%dp(:,:,:))
  param_s%dp_leap_ptr = loc(elem(2)%derived%dp(:,:,:))
  param_s%divdp_proj_s_ptr = loc(elem(1)%derived%divdp_proj(:,:,:))
  param_s%divdp_proj_leap_ptr = loc(elem(2)%derived%divdp_proj(:,:,:))
  param_s%qdp_test_ptr = loc(elem_test(1)%state%Qdp(:,:,:,:,:))
  param_s%dt = dt
  param_s%nets = nets
  param_s%nete = nete
  param_s%np1_qdp = np1_qdp
  param_s%n0_qdp = n0_qdp
  param_s%DSSopt = DSSopt
  param_s%rhs_multiplier = rhs_multiplier
  param_s%qsize = qsize
  call athread_init()
  call athread_spawn(slave_euler_step, param_s)
  call athread_join()

  qbeg = 1
  qend = qsize
  kbeg = 1
  kend = nlev
  rhs_viss = 0


end subroutine euler_step
