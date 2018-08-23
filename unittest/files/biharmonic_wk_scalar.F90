subroutine biharmonic_wk_scalar(elem,qtens,deriv,edgeq,hybrid,nets,nete)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute weak biharmonic operator
!    input:  qtens = Q
!    output: qtens = weak biharmonic of Q
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(inout), target :: elem(:)
integer :: nets,nete
real (kind=real_kind), dimension(np,np,nlev,qsize,nets:nete) :: qtens
type (EdgeBuffer_t)  , intent(inout) :: edgeq
type (derivative_t)  , intent(in) :: deriv

! local
integer :: k,kptr,i,j,ie,ic,q
real (kind=real_kind), dimension(np,np) :: lap_p
logical var_coef1

   !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad)
   !so tensor is only used on second call to laplace_sphere_wk
   var_coef1 = .true.
   if(hypervis_scaling > 0)    var_coef1 = .false.



   do ie=nets,nete
      do q=1,qsize
         do k=1,nlev    !  Potential loop inversion (AAM)
            lap_p(:,:)=qtens(:,:,k,q,ie)
! Original use of qtens on left and right hand sides caused OpenMP errors (AAM)
           !qtens(:,:,k,q,ie)=laplace_sphere_wk(lap_p,deriv,elem(ie),var_coef=var_coef1)
           call laplace_sphere_wk(lap_p,deriv,elem(ie),qtens(:,:,k,q,ie),var_coef=var_coef1)
         enddo
         kptr = nlev*(q-1)
         call edgeVpack(edgeq, qtens(:,:,1:nlev,q,ie),nlev,kptr,ie)
      enddo
   enddo

   call bndry_exchangeV(hybrid,edgeq)

   do ie=nets,nete
      ! apply inverse mass matrix, then apply laplace again
      do q=1,qsize
        kptr = nlev*(q-1)
        call edgeVunpack(edgeq, qtens(:,:,1:nlev,q,ie),nlev,kptr,ie)
        do k=1,nlev    !  Potential loop inversion (AAM)
           lap_p(:,:)=elem(ie)%rspheremp(:,:)*qtens(:,:,k,q,ie)
           ! qtens(:,:,k,q,ie)=laplace_sphere_wk(lap_p,deriv,elem(ie),var_coef=.true.)
           call laplace_sphere_wk(lap_p,deriv,elem(ie),qtens(:,:,k,q,ie),var_coef=.true.)
        enddo
      enddo
   enddo
end subroutine biharmonic_wk_scalar


subroutine laplace_sphere_wk(s,deriv,elem,laplace,var_coef)
!
!   input:  s = scalar
!   ouput:  -< grad(PHI), grad(s) >   = weak divergence of grad(s)
!     note: for this form of the operator, grad(s) does not need to be made C0
!
  real(kind=real_kind), intent(in) :: s(np,np)
  logical, intent(in) :: var_coef
  type (derivative_t), intent(in) :: deriv
  type (element_t), intent(in) :: elem
  real(kind=real_kind)             :: laplace(np,np)
  real(kind=real_kind)             :: laplace2(np,np)
  integer i,j,l

  ! Local
  real(kind=real_kind) :: grads(np,np,2), oldgrads(np,np,2)

  !call gradient_sphere(s,deriv,elem%Dinv, grads)

  real(kind=real_kind) :: ds(np,np,2)

  real(kind=real_kind) ::  dsdx00
  real(kind=real_kind) ::  dsdy00
  real(kind=real_kind) ::  v1(np,np),v2(np,np)

  do j=1,np
     do l=1,np
        dsdx00=0.0d0
        dsdy00=0.0d0
!DIR$ UNROLL(NP)
        do i=1,np
           dsdx00 = dsdx00 + deriv%Dvv(i,l  )*s(i,j  )
           dsdy00 = dsdy00 + deriv%Dvv(i,l  )*s(j  ,i)
        end do
        v1(l  ,j  ) = dsdx00*rrearth
        v2(j  ,l  ) = dsdy00*rrearth
     end do
  end do
  ! convert covarient to latlon
  !OMP_COLLAPSE_SIMD
  !DIR_VECTOR_ALIGNED
  do j=1,np
     do i=1,np
        grads(i,j,1)=elem%Dinv(i,j,1,1)*v1(i,j) + elem%Dinv(i,j,2,1)*v2(i,j)
        grads(i,j,2)=elem%Dinv(i,j,1,2)*v1(i,j) + elem%Dinv(i,j,2,2)*v2(i,j)
     enddo
  enddo
!!!!!!!!!!!!!!!! gradient_sphere !!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if (var_coef) then
     if (hypervis_power/=0 ) then
        ! scalar viscosity with variable coefficient
        grads(:,:,1) = grads(:,:,1)*elem%variable_hyperviscosity(:,:)
        grads(:,:,2) = grads(:,:,2)*elem%variable_hyperviscosity(:,:)
     else if (hypervis_scaling /=0 ) then
        ! tensor hv, (3)
        oldgrads=grads
        do j=1,np
           do i=1,np
              grads(i,j,1) = sum(oldgrads(i,j,:)*elem%tensorVisc(i,j,1,:))
              grads(i,j,2) = sum(oldgrads(i,j,:)*elem%tensorVisc(i,j,2,:))
           end do
        end do
     else
        ! do nothing: constant coefficient viscsoity
     endif
  endif

  ! note: divergnece_sphere and divergence_sphere_wk are identical *after* bndry_exchange
  ! if input is C_0.  Here input is not C_0, so we should use divergence_sphere_wk().
  call divergence_sphere_wk(grads,deriv,elem,laplace)

end subroutine laplace_sphere_wk

subroutine gradient_sphere(s,deriv,Dinv,ds)
!
!   input s:  scalar
!   output  ds: spherical gradient of s, lat-lon coordinates
!

  type (derivative_t), intent(in) :: deriv
  real(kind=real_kind), intent(in), dimension(np,np,2,2) :: Dinv
  real(kind=real_kind), intent(in) :: s(np,np)

  real(kind=real_kind) :: ds(np,np,2)

  integer i
  integer j
  integer l

  real(kind=real_kind) ::  dsdx00
  real(kind=real_kind) ::  dsdy00
  real(kind=real_kind) ::  v1(np,np),v2(np,np)

  do j=1,np
     do l=1,np
        dsdx00=0.0d0
        dsdy00=0.0d0
!DIR$ UNROLL(NP)
        do i=1,np
           dsdx00 = dsdx00 + deriv%Dvv(i,l  )*s(i,j  )
           dsdy00 = dsdy00 + deriv%Dvv(i,l  )*s(j  ,i)
        end do
        v1(l  ,j  ) = dsdx00*rrearth
        v2(j  ,l  ) = dsdy00*rrearth
     end do
  end do
  ! convert covarient to latlon
  !OMP_COLLAPSE_SIMD
  !DIR_VECTOR_ALIGNED
  do j=1,np
     do i=1,np
        ds(i,j,1)=Dinv(i,j,1,1)*v1(i,j) + Dinv(i,j,2,1)*v2(i,j)
        ds(i,j,2)=Dinv(i,j,1,2)*v1(i,j) + Dinv(i,j,2,2)*v2(i,j)
     enddo
  enddo

  end subroutine gradient_sphere

  subroutine divergence_sphere_wk(v,deriv,elem,div)
!
!   input:  v = velocity in lat-lon coordinates
!   ouput:  div(v)  spherical divergence of v, integrated by parts
!
!   Computes  -< grad(psi) dot v >
!   (the integrated by parts version of < psi div(v) > )
!
!   note: after DSS, divergence_sphere() and divergence_sphere_wk()
!   are identical to roundoff, as theory predicts.
!
    real(kind=real_kind), intent(in) :: v(np,np,2)  ! in lat-lon coordinates
    type (derivative_t), intent(in) :: deriv
    type (element_t), intent(in) :: elem
    real(kind=real_kind) :: div(np,np)

    ! Local

    integer i,j,m,n

    real(kind=real_kind) ::  vtemp(np,np,2)
    real(kind=real_kind) ::  ggtemp(np,np,2)
    real(kind=real_kind) ::  gtemp(np,np,2)
    real(kind=real_kind) ::  psi(np,np)
    real(kind=real_kind) :: xtmp

    ! latlon- > contra
    !OMP_COLLAPSE_SIMD
    !DIR_VECTOR_ALIGNED
    do j=1,np
       do i=1,np
          vtemp(i,j,1)=(elem%Dinv(i,j,1,1)*v(i,j,1) + elem%Dinv(i,j,1,2)*v(i,j,2))
          vtemp(i,j,2)=(elem%Dinv(i,j,2,1)*v(i,j,1) + elem%Dinv(i,j,2,2)*v(i,j,2))
       enddo
    enddo

    do n=1,np
       do m=1,np

          div(m,n)=0
!DIR$ UNROLL(NP)
          do j=1,np
             div(m,n)=div(m,n)-(elem%spheremp(j,n)*vtemp(j,n,1)*deriv%Dvv(m,j) &
                              +elem%spheremp(m,j)*vtemp(m,j,2)*deriv%Dvv(n,j)) &
                              * rrearth
          enddo

#if 0
! debug the above formula using the N^4 slow formulation:
          psi=0
          psi(m,n)=1
          ggtemp=gradient_sphere(psi,deriv,elem%Dinv)
          ! latlon -> covarient
          do j=1,np
             do i=1,np
                gtemp(i,j,1)=(elem%D(1,1,i,j)*ggtemp(i,j,1) + elem%D(2,1,i,j)*ggtemp(i,j,2))
                gtemp(i,j,2)=(elem%D(1,2,i,j)*ggtemp(i,j,1) + elem%D(2,2,i,j)*ggtemp(i,j,2))
             enddo
          enddo
! grad(psi) dot v:
          xtmp=0
          do j=1,np
          do i=1,np
             xtmp=xtmp-elem%spheremv(i,j)*(vtemp(i,j,1)*gtemp(i,j,1)+vtemp(i,j,2)*gtemp(i,j,2))
          enddo
          enddo
          if (abs(xtmp-div(m,n)) > 3e-17) then
             print *,m,n,xtmp,div(m,n),xtmp-div(m,n)
          endif
#endif
       end do
    end do

  end subroutine divergence_sphere_wk

  subroutine edgeVpack(edge,v,vlyr,kptr,ielem)
    use dimensions_mod, only : np, max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest

    type (EdgeBuffer_t)                :: edge
    integer,              intent(in)   :: vlyr
    real (kind=real_kind),intent(in)   :: v(np,np,vlyr)
    integer,              intent(in)   :: kptr
    integer,              intent(in)   :: ielem
!    type (EdgeDescriptor_t),intent(in) :: desc

    ! Local variables
    integer :: i,k,ir,ll,iptr

    integer :: is,ie,in,iw,edgeptr


    is = edge%putmap(south,ielem)
    ie = edge%putmap(east,ielem)
    in = edge%putmap(north,ielem)
    iw = edge%putmap(west,ielem)
!JMD    call t_adj_detailf(+2)
!dir$ ivdep
    do k=1,vlyr
       iptr = np*(kptr+k-1)
!DIR$ UNROLL(NP)
       do i=1,np
          edge%buf(iptr+ie+i)   = v(np ,i ,k) ! East
          edge%buf(iptr+is+i)   = v(i  ,1 ,k) ! South
          edge%buf(iptr+in+i)   = v(i  ,np,k) ! North
          edge%buf(iptr+iw+i)   = v(1  ,i ,k) ! West
       enddo
    enddo

    !  This is really kludgy way to setup the index reversals
    !  But since it is so a rare event not real need to spend time optimizing

    if(edge%reverse(south,ielem)) then
       do k=1,vlyr
          iptr = np*(kptr+k-1)+is
          do i=1,np
             ir = np-i+1
             edge%buf(iptr+ir)=v(i,1,k)
          enddo
       enddo
    endif

    if(edge%reverse(east,ielem)) then
       do k=1,vlyr
          iptr=np*(kptr+k-1)+ie
          do i=1,np
             ir = np-i+1
             edge%buf(iptr+ir)=v(np,i,k)
          enddo
       enddo
    endif

    if(edge%reverse(north,ielem)) then
       do k=1,vlyr
          iptr=np*(kptr+k-1)+in
          do i=1,np
             ir = np-i+1
             edge%buf(iptr+ir)=v(i,np,k)
          enddo
       enddo
    endif

    if(edge%reverse(west,ielem)) then
       do k=1,vlyr
          iptr=np*(kptr+k-1)+iw
          do i=1,np
             ir = np-i+1
             edge%buf(iptr+ir)=v(1,i,k)
          enddo
       enddo
    endif

! SWEST
    do ll=swest,swest+max_corner_elem-1
        if (edge%putmap(ll,ielem) /= -1) then
            edgeptr = edge%putmap(ll,ielem)+1
            do k=1,vlyr
                iptr = (kptr+k-1)+edgeptr
                edge%buf(iptr) = v(1, 1, k)
            end do
        end if
    end do

! SEAST
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
        if (edge%putmap(ll,ielem) /= -1) then
            edgeptr = edge%putmap(ll,ielem)+1
            do k=1,vlyr
               iptr = (kptr+k-1)+edgeptr
               edge%buf(iptr)=v(np, 1, k)
            end do
        end if
    end do

! NEAST
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
        if (edge%putmap(ll,ielem) /= -1) then
            edgeptr = edge%putmap(ll,ielem)+1
            do k=1,vlyr
               iptr = (kptr+k-1)+edgeptr
               edge%buf(iptr) = v(np, np, k)
            end do
        end if
    end do

! NWEST
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
        if (edge%putmap(ll,ielem) /= -1) then
            edgeptr = edge%putmap(ll,ielem)+1
            do k=1,vlyr
               iptr = (kptr+k-1)+edgeptr
               edge%buf(iptr) = v(1, np, k)
            end do
        end if
    end do

!JMD    call t_adj_detailf(-2)

  end subroutine edgeVpack
