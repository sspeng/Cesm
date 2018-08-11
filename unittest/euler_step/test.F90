program comp
implicit none
  !integer :: a = 5
  !integer :: b = 6
  !integer :: m = 7
  !integer :: c, d
  !c = min(a, b, m)
  !d = max(a, b, m)
  !print *, c, d
  !integer :: i, j, k
  !integer :: a(10, 10)
  !do j = 1, 10
  !  do i = 1, 10
  !    a(i, j) = j * 10 + i
  !  enddo
  !enddo
!
  !do j = 1, 10
  !  do i = 1, 10
  !    if ((i>3) .and. (i<5)) CYCLE
  !    !print *, a(i, j)
  !    if (i>6) exit
  !  enddo
  !enddo
  !print *, 0D0
  integer :: a(16), b(16), c(16)
  integer :: i, j, k, sum_a = 0
  do i = 1, 16
    a(i) = 1
    b(i) = 2
  enddo
  c = a*b
  do i = 1, 16
    print *, c(i)
  enddo
  sum_a = sum(a*b)
  print *, sum_a
end program comp
