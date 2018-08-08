program comp
implicit none
  !integer :: a = 5
  !integer :: b = 6
  !integer :: m = 7
  !integer :: c, d
  !c = min(a, b, m)
  !d = max(a, b, m)
  !print *, c, d
  integer :: i, j, k
  integer :: a(10, 10)
  do j = 1, 10
    do i = 1, 10
      a(i, j) = j * 10 + i
    enddo
  enddo

  do j = 1, 10
    do i = 1, 10
      if ((i>3) .and. (i<5)) CYCLE
      print *, a(i, j)
      if (i>6) exit
    enddo
  enddo
end program comp
