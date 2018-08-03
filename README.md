---
  title: CesmEulerStep
  tags:
    - 2018 08-02
---
#### EulerStep
- ***数据***：`qmin(:,:,nets:nete),qmax(:,:,nets:nete)),elem(ie)%spheremp(i,j),elem(ie)%derived%dpdiss_biharmonic(i,j,k),elem(ie)%derived%divdp(i,j,k),Qtens_biharmonic(i,j,k,q,ie),dp_star(i,j,k),elem(ie)%state%Qdp(i,j,k,q,n0_qdp),elem(ie)%derived%vn0(i,j,2,k),elem(ie)%derived%eta_dot_dpdn(:,:,:),elem(ie)%derived%omega_p(:,:,:),elem(ie)%derived%divdp_proj(:,:,:),`, 需要计算更新的是`qmin(:,:,nets:nete),qmax(:,:,nets:nete)),elem(ie)%state%Qdp(i,j,k,q,n0_qdp),Qtens_biharmonic(i,j,k,q,ie)`.

- ***timing***

  从下面的表格的数据发现运行时间主要花在计算`qmin, qmax, qens_Biharmonic` `divergence_sphere, limiter_optim_iter_full,`
  ***当 ( limiter_option == 4 ) 有`limiter2d_zero`, 但是B算例里没有调用该函数***

|name|ncall|time|
|:---|
|euler_step                          | 555264|960.371|
|euler_step_limiter_option_8         | 555264|329.691|
|qmin_qmax                           | 555264| 68.740|
|rhs_multiplier_0                    | 185088| 22.219|
|euler_step_2d_advec                 | 555264|587.563|
|divergence_sphere                   |1.8e+10|141.416|
|limiter_optim_iter_full             |5.9e+08|160.334|
|edgeVpack_1                         |5.9e+08| 53.602|
|edgeVpack_2                         |2.3e+07|  2.056|
|bndry_exchangeV                     | 555264| 33.990|
|euler_step_2d_advec_2               | 555264| 43.068|
|edgeVunpack_1                       |2.3e+07|  1.643|
|edgeVunpack_2                       |5.9e+08| 33.061|
|rhs_multiplier_2                    | 185088| 16.730|
|euler_step_overlap_1                | 185088|221.948|
|biharmonic_wk_scalar                | 185088|191.713|
|neighbor_minmax_finish              | 185088| 13.390|

- 函数调用关系

<img src="http://yuml.me/diagram/scruffy/class/[note:函数调用关系 {bg:cornsilk}], [euler_step]*-*>[divergence_sphere], [euler_step]-1>[limiter_optim_iter_full], [euler_step]*-*>[biharmonic_wk_scalar],  [biharmonic_wk_scalar]^[National], [biharmonic_wk_scalar]^[laplace_sphere_wk]" >


- 从核数据划分

  完成计算***Qdps_biharmonic***的从核函数的代码，将`element%state%Qdp(i,j,k,q,2)`的按ie和q轴的二维划分数据。

|ie \ q |CID1 |CID2 |CID3|CID4 |CID4 |CID6 |CID7|
|:--- |
|elem(1:6)q(1:3)|SPE1 |SPE2 |SPE3|SPE4 |SPE4 |SPE6 |elem(1:6)q(22:24)|
|SPE8 |SPE9 |SPE10 |SPE11|SPE12 |SPE13 |SPE14 |SPE15|
|SPE16 |SPE17 |SPE18 |SPE19|SPE20 |SPE21 |SPE22 |SPE23|
|SPE24 |SPE25 |SPE26 |SPE27|SPE28 |SPE29 |SPE30 |SPE31|
|SPE32 |SPE33 |SPE34 |SPE35|SPE36 |SPE37 |SPE38 |SPE39|
|SPE40 |SPE41 |SPE42 |SPE43|SPE44 |SPE45 |SPE46 |SPE47|
|SPE48 |SPE49 |SPE50 |SPE51|SPE52 |SPE53 |SPE54 |SPE55|
|elem(43:48)q(1:3) |SPE57 |SPE58 |SPE59|SPE60 |SPE61 |SPE62 |elem(42:48)q(22:24)|

- 循环结构调整

  源码里循环结构中，原本内存连续的循环层次依次是ie, q, k, j, i。
  但源码中把q层放到最下层的循环体中，会造成Qdp数组访存不连续.
``` c
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
```

  现将循环体拆成两块。

 ``` c
  do ie = nets , nete
    ! add hyperviscosity to RHS.  apply to Q at timelevel n0, Qdp(n0)/dp
    do k = 1 , nlev    !  Loop index added with implicit inversion (AAM)
      do j=1,np
        do i=1,np
          dp(i,j,k) = elem(ie)%derived%dp(i,j,k) - rhs_multiplier*dt*elem(ie)%derived%divdp_proj(i,j,k)
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
            Qtens_biharmonic(i,j,k,q,ie) = elem(ie)%state%Qdp(i,j,k,q,n0_qdp)/dp(i,j,k)
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
  ```



- 寄存器通讯

  考虑到`elem%derived%dpdiss_ave(i,j,k)`, `elem(ie)%derived%eta_dot_dpdn(:,:,:)` ,`elem(ie)%derived%omega_p(:,:,:)`,`elem(ie)%derived%divdp_proj(:,:,:)`数组与q轴无关，但可以考虑q轴第一个核取得数据，然后利用寄存器通信将获取的数组传给其他核。

  <img src="http://yuml.me/diagram/scruffy/class/[note:寄存器通信{bg:cornsilk}], [SPE0: athread_get data]*-regcomm*>[SPE7], [SPE0: athread_get data]*--regcomm*>[SPE6], [SPE0: athread_get data]*--regcomm*>[SPE5], [SPE0: athread_get data]*--regcomm*>[SPE4], [SPE0: athread_get data]*--regcomm*>[SPE2], [SPE0: athread_get data]*--regcomm*>[SPE1]" >
