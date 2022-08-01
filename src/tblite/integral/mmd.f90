! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

!> @file tblite/integral/dipole.f90
!> Provides McMurchie-Davidson recursion relations

!> Implementation of McMurchie-Davidson recursion relations
module tblite_integral_mmd
   use mctc_env, only : wp
   implicit none
   private

   public :: e_function, e_gradient, e_derivative


contains

pure subroutine e_function(lj, li, xij, rpj, rpi, e)
   integer, intent(in) :: lj, li
   real(wp), intent(in) :: xij, rpj(3), rpi(3)
   real(wp), intent(out) :: e(0:, 0:, 0:, :)

   integer :: i, j, m, n

   e(:, :, :, :) = 0.0_wp
   e(0, 0, 0, :) = 1.0_wp
   do m = 1, 3
      do i = 1, li
         e(0, 0, i, m) = &
            &   rpi(m) * e(0, 0, i-1, m) &
            & + e(1, 0, i-1, m)
         do n = 1, i-1
            e(n, 0, i, m) = &
               &   xij * e(n-1, 0, i-1, m) &
               & + rpi(m) * e(n, 0, i-1, m) &
               & + (1 + n) * e(n+1, 0, i-1, m)
         end do
         e(i, 0, i, m) = &
            &   xij * e(i-1, 0, i-1, m) &
            & + rpi(m) * e(i, 0, i-1, m)
      end do
      do j = 1, lj
         do i = 0, li
            e(0, j, i, m) = &
               &   rpj(m) * e(0, j-1, i, m) &
               & + e(1, j-1, i, m)
            do n = 1, i+j-1
               e(n, j, i, m) = &
                  &   xij * e(n-1, j-1, i, m) &
                  & + rpj(m) * e(n, j-1, i, m) &
                  & + (1 + n) * e(n+1, j-1, i, m)
            end do
            e(i+j, j, i, m) = &
               &   xij * e(i+j-1, j-1, i, m) &
               & + rpj(m) * e(i+j, j-1, i, m)
         end do
      end do
   end do
end subroutine e_function

pure subroutine e_gradient(lj, li, xij, rpj, rpi, e, drpj, drpi, de)
   integer, intent(in) :: lj, li
   real(wp), intent(in) :: xij, rpj(3), rpi(3)
   real(wp), intent(in) :: drpj(3), drpi(3)
   real(wp), intent(in) :: e(0:, 0:, 0:, :)
   real(wp), intent(out) :: de(0:, 0:, 0:, :)

   integer :: i, j, m, n

   de(:, :, :, :) = 0.0_wp
   do m = 1, 3
      do i = 1, li
         de(0, 0, i, m) = &
            &   rpi(m) * de(0, 0, i-1, m) &
            & + drpi(m) * e(0, 0, i-1, m) &
            & + de(1, 0, i-1, m)
         do n = 1, i-1
            de(n, 0, i, m) = &
               &   xij * de(n-1, 0, i-1, m) &
               & + rpi(m) * de(n, 0, i-1, m) &
               & + drpi(m) * e(n, 0, i-1, m) &
               & + (1 + n) * de(n+1, 0, i-1, m)
         end do
         de(i, 0, i, m) = &
            &   xij * de(i-1, 0, i-1, m) &
            & + rpi(m) * de(i, 0, i-1, m) &
            & + drpi(m) * e(i, 0, i-1, m)
      end do
      do j = 1, lj
         do i = 0, li
            de(0, j, i, m) = &
               &   rpj(m) * de(0, j-1, i, m) &
               & + drpj(m) * e(0, j-1, i, m) &
               & + de(1, j-1, i, m)
            do n = 1, i+j-1
               de(n, j, i, m) = &
                  &   xij * de(n-1, j-1, i, m) &
                  & + rpj(m) * de(n, j-1, i, m) &
                  & + drpj(m) * e(n, j-1, i, m) &
                  & + (1 + n) * de(n+1, j-1, i, m)
            end do
            de(i+j, j, i, m) = &
               &   xij * de(i+j-1, j-1, i, m) &
               & + rpj(m) * de(i+j, j-1, i, m) &
               & + drpj(m) * e(i+j, j-1, i, m)
         end do
      end do
   end do
end subroutine e_gradient

pure subroutine e_derivative(e, ai, lj, li, de)
   real(wp), intent(in) :: e(0:, 0:, 0:, :)
   real(wp), intent(in) :: ai
   integer, intent(in) :: lj, li
   real(wp), intent(out) :: de(0:, 0:, 0:, :)

   integer :: m, i, j

   de(:, :, :, :) = 0.0_wp
   do m = 1, 3
      do i = 0, li
         do j = 0, lj
            de(0, j, i, m) = 2 * ai * e(0, j, i+1, m)
            if (i > 0) de(0, j, i, m) = de(0, j, i, m) - i * e(0, j, i-1, m)
         end do
      end do
   end do
end subroutine e_derivative

end module tblite_integral_mmd
