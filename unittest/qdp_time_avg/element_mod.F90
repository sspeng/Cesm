module element_mod
  implicit none
  integer, public, parameter :: real_kind = 8
  integer, public, parameter :: timelevels = 3
  integer, public, parameter :: np = 4
  integer, public, parameter :: nlev = 30
  integer, public, parameter :: qsize_d = 25
  integer, public, parameter :: qsize = 25

  type, public :: elem_state_t
    real (kind=real_kind) :: v   (np,np,2,nlev,timelevels)            ! velocity                           1
    real (kind=real_kind) :: T   (np,np,nlev,timelevels)              ! temperature                        2
    real (kind=real_kind) :: dp3d(np,np,nlev,timelevels)              ! delta p on levels                  8
    real (kind=real_kind) :: lnps(np,np,timelevels)                   ! log surface pressure               3
    real (kind=real_kind) :: ps_v(np,np,timelevels)                   ! surface pressure                   4
    real (kind=real_kind) :: phis(np,np)                              ! surface geopotential (prescribed)  5
    real (kind=real_kind) :: Q   (np,np,nlev,qsize_d)                 ! Tracer concentration               6
    real (kind=real_kind) :: Qdp (np,np,nlev,qsize_d,2)               ! Tracer mass                        7
  end type elem_state_t

  type, public :: element_t
     type (elem_state_t) :: state
  end type element_t

end module element_mod
