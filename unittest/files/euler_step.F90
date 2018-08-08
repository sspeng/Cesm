subroutine euler_step( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
! ===================================
! This routine is the basic foward
! euler component used to construct RK SSP methods
!
!           u(np1) = u(n0) + dt2*DSS[ RHS(u(n0)) ]
!
! n0 can be the same as np1.
!
! DSSopt = DSSeta or DSSomega:   also DSS eta_dot_dpdn or omega
!
! ===================================
use perf_mod      , only : t_startf, t_stopf            ! _EXTERNAL
use kinds          , only : real_kind
use dimensions_mod , only : np, npdg, nlev
use hybrid_mod     , only : hybrid_t
use element_mod    , only : element_t
use derivative_mod , only : derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
use edge_mod       , only : edgevpack, edgevunpack
use bndry_mod      , only : bndry_exchangev
use hybvcoord_mod  , only : hvcoord_t
#if USE_CUDA_FORTRAN
use cuda_mod, only: euler_step_cuda
#endif
implicit none
integer              , intent(in   )         :: np1_qdp, n0_qdp
real (kind=real_kind), intent(in   )         :: dt
type (element_t)     , intent(inout), target :: elem(:)
type (hvcoord_t)     , intent(in   )         :: hvcoord
type (hybrid_t)      , intent(in   )         :: hybrid
type (derivative_t)  , intent(in   )         :: deriv
integer              , intent(in   )         :: nets
integer              , intent(in   )         :: nete
integer              , intent(in   )         :: DSSopt
integer              , intent(in   )         :: rhs_multiplier

! local
real(kind=real_kind), dimension(np,np                       ) :: divdp, dpdiss
real(kind=real_kind), dimension(np,np                       ) :: div
real(kind=real_kind), dimension(np,np,nlev                  ) :: dpdissk
real(kind=real_kind), dimension(np,np,2                     ) :: gradQ
real(kind=real_kind), dimension(np,np,2,nlev                ) :: Vstar
real(kind=real_kind), dimension(np,np  ,nlev                ) :: Qtens
real(kind=real_kind), dimension(np,np  ,nlev                ) :: dp,dp_star
real(kind=real_kind), dimension(np,np  ,nlev,qsize,nets:nete) :: Qtens_biharmonic
real(kind=real_kind), pointer, dimension(:,:,:)               :: DSSvar
real(kind=real_kind) :: dp0(nlev),qmin_val(nlev),qmax_val(nlev)
integer :: ie,q,i,j,k
integer :: rhs_viss = 0
integer :: qbeg, qend, kbeg, kend
integer :: kptr

!call t_startf('euler_step')
do k = 1 , nlev
  dp0(k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
           ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
enddo
if ( npdg > 0 ) then
  call euler_step_dg( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
  return
endif
#if USE_CUDA_FORTRAN
call euler_step_cuda( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
return
#endif
call t_barrierf('sync_euler_step', hybrid%par%comm)
call t_startf('euler_step')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   compute Q min/max values for lim8
!   compute biharmonic mixing term f
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
qbeg = 1
qend = qsize
kbeg = 1
kend = nlev

rhs_viss = 0
if ( limiter_option == 8  ) then
  ! when running lim8, we also need to limit the biharmonic, so that term needs
  ! to be included in each euler step.  three possible algorithms here:
  ! 1) most expensive:
  !     compute biharmonic (which also computes qmin/qmax) during all 3 stages
  !     be sure to set rhs_viss=1
  !     cost:  3 biharmonic steps with 3 DSS
  !
  ! 2) cheapest:
  !     compute biharmonic (which also computes qmin/qmax) only on first stage
  !     be sure to set rhs_viss=3
  !     reuse qmin/qmax for all following stages (but update based on local qmin/qmax)
  !     cost:  1 biharmonic steps with 1 DSS
  !     main concern:  viscosity
  !
  ! 3)  compromise:
  !     compute biharmonic (which also computes qmin/qmax) only on last stage
  !     be sure to set rhs_viss=3
  !     compute qmin/qmax directly on first stage
  !     reuse qmin/qmax for 2nd stage stage (but update based on local qmin/qmax)
  !     cost:  1 biharmonic steps, 2 DSS
  !
  !  NOTE  when nu_p=0 (no dissipation applied in dynamics to dp equation), we should
  !        apply dissipation to Q (not Qdp) to preserve Q=1
  !        i.e.  laplace(Qdp) ~  dp0 laplace(Q)
  !        for nu_p=nu_q>0, we need to apply dissipation to Q * diffusion_dp
  !
  ! initialize dp, and compute Q from Qdp (and store Q in Qtens_biharmonic)
  call t_startf('euler_step_limiter_option_8')
  call t_startf('qmin_qmax')
  do ie = nets , nete
    ! add hyperviscosity to RHS.  apply to Q at timelevel n0, Qdp(n0)/dp
    do k = 1 , nlev    !  Loop index added with implicit inversion (AAM)
      do j=1,np
        do i=1,np
          dp(i,j,k) = elem(ie)%derived%dp(i,j,k) - rhs_multiplier*dt*elem(ie)%derived%divdp_proj(i,j,k)
          do q = 1, qsize
            Qtens_biharmonic(i,j,k,q,ie) = elem(ie)%state%Qdp(i,j,k,q,n0_qdp)/dp(i,j,k)
          enddo
        enddo
      enddo
    enddo

    do q = 1, qsize
      do k= 1, nlev
        qmin_val(k) = +1.0e+24
        qmax_val(k) = -1.0e+24
        do j=1,np
          do i=1,np
!             Qtens_biharmonic(i,j,k,q,ie) = elem(ie)%state%Qdp(i,j,k,q,n0_qdp)/dp(i,j,k)
            qmin_val(k) = min(qmin_val(k),Qtens_biharmonic(i,j,k,q,ie))
            qmax_val(k) = max(qmax_val(k),Qtens_biharmonic(i,j,k,q,ie))
          enddo
        enddo

        if ( rhs_multiplier == 1 ) then
            qmin(k,q,ie)=min(qmin(k,q,ie),qmin_val(k))
            qmin(k,q,ie)=max(qmin(k,q,ie),0d0)
            qmax(k,q,ie)=max(qmax(k,q,ie),qmax_val(k))
        else
            qmin(k,q,ie)=max(qmin_val(k),0d0)
            qmax(k,q,ie)=qmax_val(k)
        endif
      enddo
    enddo
  enddo
  call t_stopf('qmin_qmax')
  if ( rhs_multiplier == 0 ) then
    ! update qmin/qmax based on neighbor data for lim8
    call t_startf('rhs_multiplier_0')
    call neighbor_minmax(hybrid,edgeAdvQminmax,nets,nete,qmin(:,:,nets:nete),qmax(:,:,nets:nete))
    call t_stopf('rhs_multiplier_0')
  endif

  if ( rhs_multiplier == 2 ) then
    call t_startf('rhs_multiplier_2')
     rhs_viss = 3
    ! two scalings depending on nu_p:
    ! nu_p=0:    qtens_biharmonic *= dp0                   (apply viscsoity only to q)
    ! nu_p>0):   qtens_biharmonc *= elem()%psdiss_ave      (for consistency, if nu_p=nu_q)
    if ( nu_p > 0 ) then
      do ie = nets , nete
        do k = 1 , nlev
          do j=1,np
            do i=1,np
              dpdiss(i,j) = elem(ie)%derived%dpdiss_ave(i,j,k)
            enddo
          enddo
          do q = 1 , qsize
            ! NOTE: divide by dp0 since we multiply by dp0 below
            do j=1,np
              do i=1,np
                Qtens_biharmonic(i,j,k,q,ie)=Qtens_biharmonic(i,j,k,q,ie)*dpdiss(i,j)/dp0(k)
              enddo
            enddo
          enddo
        enddo
      enddo
    endif
    call t_stopf('rhs_multiplier_2')
#ifdef OVERLAP
    call t_startf('euler_step_overlap_1')
    call t_startf('neighbor_minmax_start_and_scalar_1')
    call neighbor_minmax_start(hybrid,edgeAdvQminmax,nets,nete,qmin(:,:,nets:nete),qmax(:,:,nets:nete))
    call biharmonic_wk_scalar(elem,qtens_biharmonic,deriv,edgeAdv,hybrid,nets,nete)
    call t_stopf('neighbor_minmax_start_and_scalar_1')
    do ie = nets, nete
      do q = 1, qsize
        do k = 1, nlev
          !OMP_COLLAPSE_SIMD
          !DIR_VECTOR_ALIGNED
          do j=1,np
            do i=1,np
              ! note: biharmonic_wk() output has mass matrix already applied. Un-apply since we apply again below:
              qtens_biharmonic(i,j,k,q,ie) = &
                   -rhs_viss*dt*nu_q*dp0(k)*Qtens_biharmonic(i,j,k,q,ie) / elem(ie)%spheremp(i,j)
            enddo
          enddo
        enddo
      enddo
    enddo

    call t_startf('neighbor_minmax_finish')
    call neighbor_minmax_finish(hybrid,edgeAdvQminmax,nets,nete,qmin(:,:,nets:nete),qmax(:,:,nets:nete))
    call t_stopf('neighbor_minmax_finish')
    call t_stopf('euler_step_overlap_1')
#else
    call t_startf('neighbor_minmax_start_and_scalar_2')
    call neighbor_minmax(hybrid,edgeAdvQminmax,nets,nete,qmin(:,:,nets:nete),qmax(:,:,nets:nete))
    call biharmonic_wk_scalar(elem,qtens_biharmonic,deriv,edgeAdv,hybrid,nets,nete)
    call t_stopf('neighbor_minmax_start_and_scalar_2')
    call t_startf('euler_step_overlap_2')
    do ie = nets, nete
      do q = 1, qsize
        do k = 1, nlev
          !OMP_COLLAPSE_SIMD
          !DIR_VECTOR_ALIGNED
          do j=1,np
            do i=1,np
          ! note: biharmonic_wk() output has mass matrix already applied. Un-apply since we apply again below:
              qtens_biharmonic(i,j,k,q,ie) = &
                   -rhs_viss*dt*nu_q*dp0(k)*Qtens_biharmonic(i,j,k,q,ie) / elem(ie)%spheremp(i,j)
            enddo
          enddo
        enddo
      enddo
    enddo
    call t_stopf('euler_step_overlap_2')
#endif

!      call biharmonic_wk_scalar_minmax( elem , qtens_biharmonic , deriv , edgeAdvQ3 , hybrid , &
!           nets , nete , qmin(:,:,nets:nete) , qmax(:,:,nets:nete) )
!      do ie = nets , nete
!        do k = 1 , nlev
!          do q = 1 , qsize
!            ! note: biharmonic_wk() output has mass matrix already applied. Un-apply since we apply again below:
!            do j=1,np
!              do i=1,np
!                qtens_biharmonic(i,j,k,q,ie) = &
!                     -rhs_viss*dt*nu_q*dp0(k)*Qtens_biharmonic(i,j,k,q,ie) / elem(ie)%spheremp(i,j)
!              enddo
!            enddo
!          enddo
!        enddo
!      enddo
!
  endif
endif  ! compute biharmonic mixing term and qmin/qmax
call t_stopf('euler_step_limiter_option_8')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   2D Advection step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call t_startf('euler_step_2d_advec')
do ie = nets , nete
  ! note: eta_dot_dpdn is actually dimension nlev+1, but nlev+1 data is
  ! all zero so we only have to DSS 1:nlev
  if ( DSSopt == DSSeta         ) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
  if ( DSSopt == DSSomega       ) DSSvar => elem(ie)%derived%omega_p(:,:,:)
  if ( DSSopt == DSSdiv_vdp_ave ) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)

  ! Compute velocity used to advance Qdp
  do k = 1 , nlev    !  Loop index added (AAM)
    ! derived variable divdp_proj() (DSS'd version of divdp) will only be correct on 2nd and 3rd stage
    ! but that's ok because rhs_multiplier=0 on the first stage:
    do j=1,np
      do i=1,np
        dp(i,j,k) = elem(ie)%derived%dp(i,j,k) - rhs_multiplier * dt * elem(ie)%derived%divdp_proj(i,j,k)
        Vstar(i,j,1,k) = elem(ie)%derived%vn0(i,j,1,k) / dp(i,j,k)
        Vstar(i,j,2,k) = elem(ie)%derived%vn0(i,j,2,k) / dp(i,j,k)
      enddo
    enddo
  enddo

  ! advance Qdp
  do q = 1 , qsize
    do k = 1 , nlev  !  dp_star used as temporary instead of divdp (AAM)
      ! div( U dp Q),

      do j=1,np
        do i=1,np
          gradQ(i,j,1) = Vstar(i,j,1,k) * elem(ie)%state%Qdp(i,j,k,q,n0_qdp)
          gradQ(i,j,2) = Vstar(i,j,2,k) * elem(ie)%state%Qdp(i,j,k,q,n0_qdp)
        enddo
      enddo

      ! dp_star(:,:,k) = divergence_sphere( gradQ , deriv , elem(ie) )
      call t_startf('divergence_sphere')
      call divergence_sphere( gradQ , deriv , elem(ie), dp_star(:,:,k) )
      call t_stopf('divergence_sphere')
      do j=1,np
        do i=1,np
          Qtens(i,j,k) = elem(ie)%state%Qdp(i,j,k,q,n0_qdp) - dt * dp_star(i,j,k)
        enddo
      enddo

      ! optionally add in hyperviscosity computed above:
      if ( rhs_viss /= 0 ) then
        do j=1,np
          do i=1,np
            Qtens(i,j,k) = Qtens(i,j,k) + Qtens_biharmonic(i,j,k,q,ie)
          enddo
        enddo
      endif
    enddo

    if ( limiter_option == 8) then
      do k = 1 , nlev  ! Loop index added (AAM)
        ! UN-DSS'ed dp at timelevel n0+1:
        do j=1,np
          do i=1,np
            dp_star(i,j,k) = dp(i,j,k) - dt * elem(ie)%derived%divdp(i,j,k)
          enddo
        enddo

        if ( nu_p > 0 .and. rhs_viss /= 0 ) then
          ! add contribution from UN-DSS'ed PS dissipation
!            dpdiss(:,:) = ( hvcoord%hybi(k+1) - hvcoord%hybi(k) ) * elem(ie)%derived%psdiss_biharmonic(:,:)
         do j=1,np
            do i=1,np
                dpdiss(i,j) = elem(ie)%derived%dpdiss_biharmonic(i,j,k)
             dp_star(i,j,k) = dp_star(i,j,k) - rhs_viss * dt * nu_q * dpdiss(i,j) / elem(ie)%spheremp(i,j)
            enddo
          enddo
        endif
      enddo
      ! apply limiter to Q = Qtens / dp_star
      call t_startf('limiter_optim_iter_full')
      call limiter_optim_iter_full( Qtens(:,:,:) , elem(ie)%spheremp(:,:) , qmin(:,q,ie) , &
                                    qmax(:,q,ie) , dp_star(:,:,:))
      call t_stopf('limiter_optim_iter_full')
    endif


    ! apply mass matrix, overwrite np1 with solution:
    ! dont do this earlier, since we allow np1_qdp == n0_qdp
    ! and we dont want to overwrite n0_qdp until we are done using it
    do k = 1 , nlev
      do j=1,np
        do i=1,np
          elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = elem(ie)%spheremp(i,j) * Qtens(i,j,k)
        enddo
      enddo
    enddo

    if ( limiter_option == 4 ) then
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      ! sign-preserving limiter, applied after mass matrix
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      call t_startf('limiter2d_zero')
      call limiter2d_zero( elem(ie)%state%Qdp(:,:,:,q,np1_qdp))
      call t_stopf('limiter2d_zero')
    endif

    kptr = nlev*(q-1)
    call t_startf('edgeVpack_1')
    call edgeVpack(edgeAdvp1  , elem(ie)%state%Qdp(:,:,:,q,np1_qdp) , nlev , kptr , ie )
    call t_stopf('edgeVpack_1')
   enddo


   if ( DSSopt == DSSeta         ) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
   if ( DSSopt == DSSomega       ) DSSvar => elem(ie)%derived%omega_p(:,:,:)
   if ( DSSopt == DSSdiv_vdp_ave ) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)
   ! also DSS extra field
   do k = 1 , nlev
     do j=1,np
        do i=1,np
          DSSvar(i,j,k) = elem(ie)%spheremp(i,j) * DSSvar(i,j,k)
        enddo
     enddo
   enddo

   kptr = nlev*qsize
   call t_startf('edgeVpack_2')
   call edgeVpack( edgeAdvp1 , DSSvar(:,:,1:nlev) , nlev , kptr , ie )
   call t_stopf('edgeVpack_2')
enddo

call t_startf('bndry_exchangeV')
call bndry_exchangeV( hybrid , edgeAdvp1    )
call t_stopf('bndry_exchangeV')
call t_stopf('euler_step_2d_advec')

call t_startf('euler_step_2d_advec_2')
do ie = nets , nete
  ! only perform this operation on thread which owns the first tracer
     if ( DSSopt == DSSeta         ) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
     if ( DSSopt == DSSomega       ) DSSvar => elem(ie)%derived%omega_p(:,:,:)
     if ( DSSopt == DSSdiv_vdp_ave ) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)
     kptr = qsize*nlev + kbeg -1
     call t_startf('edgeVunpack_1')
     call edgeVunpack( edgeAdvp1 , DSSvar(:,:,kbeg:kend) , nlev , kptr , ie )
     call t_stopf('edgeVunpack_1')
     do k = 1, nlev
       !OMP_COLLAPSE_SIMD
       !DIR_VECTOR_ALIGNED
       do j=1,np
         do i=1,np
           DSSvar(i,j,k) = DSSvar(i,j,k) * elem(ie)%rspheremp(i,j)
         enddo
       enddo
     enddo
  do q = 1, qsize
    kptr = nlev*(q-1) + kbeg - 1
    call t_startf('edgeVunpack_2')
    call edgeVunpack( edgeAdvp1 , elem(ie)%state%Qdp(:,:,kbeg:kend,q,np1_qdp) , nlev , kptr , ie )
    call t_stopf('edgeVunpack_2')
      do k = 1, nlev
        !OMP_COLLAPSE_SIMD
        !DIR_VECTOR_ALIGNED
        do j=1,np
          do i=1,np
            elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = elem(ie)%rspheremp(i,j) * elem(ie)%state%Qdp(i,j,k,q,np1_qdp)
          enddo
        enddo
      enddo
  enddo

enddo
call t_stopf('euler_step_2d_advec_2')
call t_stopf('euler_step')
end subroutine euler_step







#ifdef LIMITER_REWRITE_OPT
  subroutine limiter_optim_iter_full(ptens,sphweights,minp,maxp,dpmass)
    !
    !The idea here is the following: We need to find a grid field which is closest
    !to the initial field (in terms of weighted sum), but satisfies the min/max constraints.
    !So, first we find values which do not satisfy constraints and bring these values
    !to a closest constraint. This way we introduce some mass change (addmass),
    !so, we redistribute addmass in the way that l2 error is smallest.
    !This redistribution might violate constraints thus, we do a few iterations.
    !
    ! O. Guba ~2012                    Documented in Guba, Taylor & St-Cyr, JCP 2014
    ! I. Demeshko & M. Taylor 7/2015:  Removed indirect addressing.
    ! N. Lopez & M. Taylor 8/2015:     Mass redistributon tweak which is better at
    !                                  linear coorelation preservation
    !
    use kinds         , only : real_kind
    use dimensions_mod, only : np, nlev

    real (kind=real_kind), dimension(nlev), intent(inout)   :: minp, maxp
    real (kind=real_kind), dimension(np*np,nlev), intent(inout)   :: ptens
    real (kind=real_kind), dimension(np*np,nlev), intent(in), optional  :: dpmass
    real (kind=real_kind), dimension(np*np), intent(in)   :: sphweights

    real (kind=real_kind), dimension(np,np) :: ptens_mass
    integer  k1, k, i, j, iter, weightsnum
    real (kind=real_kind) :: addmass, weightssum, mass, sumc
    real (kind=real_kind) :: x(np*np),c(np*np)
    integer :: maxiter = np*np-1
    real (kind=real_kind) :: tol_limiter = 5e-14

    do k = 1, nlev

     do k1=1,np*np
       c(k1)=sphweights(k1)*dpmass(k1,k)
       x(k1)=ptens(k1,k)/dpmass(k1,k)
     enddo

     sumc=sum(c)
     if (sumc <= 0 ) CYCLE   ! this should never happen, but if it does, dont limit
     mass=sum(c*x)

      ! relax constraints to ensure limiter has a solution:
      ! This is only needed if runnign with the SSP CFL>1 or
      ! due to roundoff errors
      if( mass < minp(k)*sumc ) then
        minp(k) = mass / sumc
      endif
      if( mass > maxp(k)*sumc ) then
        maxp(k) = mass / sumc
      endif

      do iter=1,maxiter

      addmass=0.0d0

       do k1=1,np*np
         if((x(k1)>maxp(k))) then
           addmass=addmass+(x(k1)-maxp(k))*c(k1)
           x(k1)=maxp(k)
         endif
         if((x(k1)<minp(k))) then
           addmass=addmass-(minp(k)-x(k1))*c(k1)
           x(k1)=minp(k)
         endif
       enddo !k1

       if(abs(addmass)<=tol_limiter*abs(mass)) exit

       weightssum=0.0d0
!       weightsnum=0
       if(addmass>0)then
        do k1=1,np*np
          if(x(k1)<maxp(k))then
            weightssum=weightssum+c(k1)
!            weightsnum=weightsnum+1
          endif
        enddo !k1
        do k1=1,np*np
          if(x(k1)<maxp(k))then
              x(k1)=x(k1)+addmass/weightssum
!              x(k1)=x(k1)+addmass/(c(k1)*weightsnum)
          endif
        enddo
      else
        do k1=1,np*np
          if(x(k1)>minp(k))then
            weightssum=weightssum+c(k1)
!            weightsnum=weightsnum+1
          endif
        enddo
        do k1=1,np*np
          if(x(k1)>minp(k))then
            x(k1)=x(k1)+addmass/weightssum
!           x(k1)=x(k1)+addmass/(c(k1)*weightsnum)
          endif
        enddo
      endif


   enddo!end of iteration

   do k1=1,np*np
      ptens(k1,k)=x(k1)
   enddo

  enddo

  do k = 1, nlev
    do k1=1,np*np
      ptens(k1,k)=ptens(k1,k)*dpmass(k1,k)
    enddo
  enddo

  end subroutine limiter_optim_iter_full
#endif
