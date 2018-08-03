program comp
implicit none
  integer :: a = 5
  integer :: b = 6
  integer :: m = 7
  integer :: c, d
  c = min(a, b, m)
  d = max(a, b, m)
  print *, c, d
end program comp
