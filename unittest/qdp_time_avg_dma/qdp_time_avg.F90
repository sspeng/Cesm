program main
use element_mod, only: element_t
use element_mod, only: real_kind
use element_mod, only: qsize, np, nlev
implicit none
  integer :: nets = 1
  integer :: nete = 42
  integer :: limiter_option = 8
  integer :: rkstage = 3
  integer :: n0_qdp = 1
  integer :: np1_qdp = 2
  real(kind=real_kind) :: nu_p = 0
  type(element_t) :: elem(42)
  integer :: ie, i, j, k, q
  interface
    subroutine qdp_time_avg(elem , rkstage , n0_qdp , np1_qdp , limiter_option , nu_p , nets , nete )
      use element_mod, only: element_t
      use element_mod, only: real_kind
      integer :: rkstage, n0_qdp, np1_qdp, limiter_option, nets, nete
      real(kind=real_kind) :: nu_p
      type(element_t) :: elem(:)
    end subroutine qdp_time_avg
  end interface

  do ie=nets,nete
    do q=1,qsize
      do k=1,nlev
        do j=1,np
          do i=1,np
            elem(ie)%state%Qdp(i,j,k,q,n0_qdp) = ie * 1000000 + q * 10000 + k * 100 + j * 10 + i + 0.1
            elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = ie * 1000000 + q * 10000 + k * 100 + j * 10 + i + 0.2
          enddo
        enddo
      enddo
    enddo
  enddo

  call qdp_time_avg(elem, rkstage, n0_qdp, np1_qdp, limiter_option, nu_p, nets, nete)

end program main

subroutine qdp_time_avg(elem , rkstage , n0_qdp , np1_qdp , limiter_option , nu_p , nets , nete )
use element_mod, only: element_t
use element_mod, only: real_kind
use element_mod, only: qsize, np, nlev
implicit none
  integer, intent(in) :: nets, nete, n0_qdp, np1_qdp, rkstage, limiter_option
  real(kind=real_kind) :: nu_p
  integer :: ie, i, j, k, q
  type(element_t) , intent(inout), target :: elem(:)
  external :: slave_qdp_time_avg

  real(kind=real_kind) :: a


  type param
    !type(element_t), pointer :: elem(:)
    integer*4 :: nets, nete, np, nlev, n0_qdp, np1_qdp, rkstage, limiter_option, qsize
    integer*8 :: qdp_ptr1, qdp_ptr2
  end type param

  type(param) :: param_s
  param_s%qsize = qsize
  param_s%np = np
  param_s%nlev = nlev
  param_s%nets = nets
  param_s%nete = nete
  param_s%n0_qdp = n0_qdp
  param_s%np1_qdp = np1_qdp
  param_s%rkstage = rkstage
  param_s%limiter_option = limiter_option
  param_s%qdp_ptr1 = loc(elem(1)%state%Qdp(:,:,:,:,:))
  param_s%qdp_ptr2 = loc(elem(2)%state%Qdp(:,:,:,:,:))

  print *, elem(1)%state%Qdp(1,1,1,1,1), param_s%np, param_s%nlev

  call athread_init()
  call athread_spawn(slave_qdp_time_avg, param_s)
  call athread_join()
  print *, elem(1)%state%Qdp(1,1,1,1,1), param_s%np, param_s%nlev
  do ie=nets,nete
    do q=1,qsize
      do k=1,nlev
        do j=1,np
          do i=1,np
            elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = ( elem(ie)%state%Qdp(i,j,k,q,n0_qdp) + &
              (rkstage-1)*elem(ie)%state%Qdp(i,j,k,q,np1_qdp) ) / rkstage
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine qdp_time_avg
