! I am writing this subroutine in an effort to speed up the Python code that I have written for parsing large number of PDB files, in order to calculate the Weighted Contact Number for sets of coordinates of different atoms.
! This is done by calling this Fortran function from my Python wrapper get_res_prop_wcn_bf_f90.py

! Amir Shahmoradi, Thursday 1:56 PM, Aug 14 2014, iCMB, UT Austin

subroutine get_wcn_invsq(natoms,crd,wcn)
implicit none
integer, intent(in) :: natoms
real*8, intent(in)  :: crd(natoms,3)
real*8, intent(out) :: wcn(natoms)
integer :: i,j
real*8 :: distance_sq

wcn = 0.d0

do i=1,natoms
  do j=1,natoms
    if (j/=i) then
      distance_sq = (crd(i,1)-crd(j,1))**2 + (crd(i,2)-crd(j,2))**2 + (crd(i,3)-crd(j,3))**2
      wcn(i) = wcn(i) + 1.d0/distance_sq
    end if
  end do
end do

end subroutine get_wcn_invsq