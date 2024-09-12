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

!> @file tblite/integral/overlap.f90
!> Provides evaluation of overlap integrals

!> Implementation of overlap integrals
module tblite_integral_overlap
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use tblite_basis_type, only : basis_type, cgto_type
   use tblite_integral_trafo, only : transform0, transform1, transform2
   use tblite_integral_mmd, only : e_function, e_derivative, e_gradient
   implicit none
   private

   public :: overlap_cgto, overlap_grad_cgto
   public :: get_overlap
   public :: maxl, msao

   interface get_overlap
      module procedure :: get_overlap_lat
   end interface get_overlap

   integer, parameter :: maxl = 6
   integer, parameter :: maxl2 = maxl*2
   integer, parameter :: msao(0:maxl) = [1, 3, 5, 7, 9, 11, 13]
   integer, parameter :: mlao(0:maxl) = [1, 3, 6, 10, 15, 21, 28]
   integer, parameter :: lmap(0:maxl) = [0, 1, 4, 10, 20, 35, 56]
   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: sqrtpi3 = sqrtpi**3

   ! x (+1), y (-1), z (0) in [-1, 0, 1] sorting
   integer, parameter :: lx(3, 84) = reshape([&
      & 0, &
      & 0,0,1, &
      & 2,0,0,1,1,0, &
      & 3,0,0,2,2,1,0,1,0,1, &
      & 4,0,0,3,3,1,0,1,0,2,2,0,2,1,1, &
      & 5,0,0,3,3,2,2,0,0,4,4,1,0,0,1,1,3,1,2,2,1, &
      & 6,0,0,3,3,0,5,5,1,0,0,1,4,4,2,0,2,0,3,3,1,2,2,1,4,1,1,2, &
      & 0, &
      & 1,0,0, &
      & 0,2,0,1,0,1, &
      & 0,3,0,1,0,2,2,0,1,1, &
      & 0,4,0,1,0,3,3,0,1,2,0,2,1,2,1, &
      & 0,5,0,2,0,3,0,3,2,1,0,4,4,1,0,1,1,3,2,1,2, &
      & 0,6,0,3,0,3,1,0,0,1,5,5,2,0,0,2,4,4,2,1,3,1,3,2,1,4,1,2, &
      & 0, &
      & 0,1,0, &
      & 0,0,2,0,1,1, &
      & 0,0,3,0,1,0,1,2,2,1, &
      & 0,0,4,0,1,0,1,3,3,0,2,2,1,1,2, &
      & 0,0,5,0,2,0,3,2,3,0,1,0,1,4,4,3,1,1,1,2,2, &
      & 0,0,6,0,3,3,0,1,5,5,1,0,0,2,4,4,0,2,1,2,2,3,1,3,1,1,4,2], &
      & shape(lx), order=[2, 1])


contains


pure subroutine overlap_cgto(cgtoj, cgtoi, r2, vec, intcut, overlap)
   !> Description of contracted Gaussian function on center i
   type(cgto_type), intent(in) :: cgtoi
   !> Description of contracted Gaussian function on center j
   type(cgto_type), intent(in) :: cgtoj
   !> Square distance between center i and j
   real(wp), intent(in) :: r2
   !> Distance vector between center i and j, ri - rj
   real(wp), intent(in) :: vec(3)
   !> Maximum value of integral prefactor to consider
   real(wp), intent(in) :: intcut
   !> Overlap integrals for the given pair i  and j
   real(wp), intent(out) :: overlap(msao(cgtoj%ang), msao(cgtoi%ang))

   integer :: ip, jp, mli, mlj, li(3), lj(3)
   real(wp) :: eab, oab, est, rpi(3), rpj(3), cc, val, pre, e0(3)
   real(wp) :: s3d(mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: et(0:maxl2, 0:maxl, 0:maxl, 3)

   s3d(:, :) = 0.0_wp

   do ip = 1, cgtoi%nprim
      do jp = 1, cgtoj%nprim
         eab = cgtoi%alpha(ip) + cgtoj%alpha(jp)
         oab = 1.0_wp/eab
         est = cgtoi%alpha(ip) * cgtoj%alpha(jp) * r2 * oab
         if (est > intcut) cycle
         pre = exp(-est) * sqrtpi3*sqrt(oab)**3
         rpi = -vec * cgtoj%alpha(jp) * oab
         rpj = +vec * cgtoi%alpha(ip) * oab
         cc = cgtoi%coeff(ip) * cgtoj%coeff(jp) * pre
         call e_function(cgtoj%ang, cgtoi%ang, 0.5_wp*oab, rpj, rpi, et)
         do mli = 1, mlao(cgtoi%ang)
            do mlj = 1, mlao(cgtoj%ang)
               li = lx(:, mli+lmap(cgtoi%ang))
               lj = lx(:, mlj+lmap(cgtoj%ang))
               e0 = [et(0, lj(1), li(1), 1), et(0, lj(2), li(2), 2), et(0, lj(3), li(3), 3)]
               s3d(mlj, mli) = s3d(mlj, mli) + cc * product(e0)
            end do
         end do
      end do
   end do

   call transform0(cgtoj%ang, cgtoi%ang, s3d, overlap)

end subroutine overlap_cgto


subroutine overlap_grad_cgto(cgtoj, cgtoi, r2, vec, intcut, overlap, doverlap)
   !> Description of contracted Gaussian function on center i
   type(cgto_type), intent(in) :: cgtoi
   !> Description of contracted Gaussian function on center j
   type(cgto_type), intent(in) :: cgtoj
   !> Square distance between center i and j
   real(wp), intent(in) :: r2
   !> Distance vector between center i and j, ri - rj
   real(wp), intent(in) :: vec(3)
   !> Maximum value of integral prefactor to consider
   real(wp), intent(in) :: intcut
   !> Overlap integrals for the given pair i  and j
   real(wp), intent(out) :: overlap(msao(cgtoj%ang), msao(cgtoi%ang))
   !> Overlap integral gradient for the given pair i  and j
   real(wp), intent(out) :: doverlap(3, msao(cgtoj%ang), msao(cgtoi%ang))

   integer :: ip, jp, mli, mlj, li(3), lj(3)
   real(wp) :: eab, oab, est, rpi(3), rpj(3), cc, val, grad(3), pre, e0(3), d0(3)
   real(wp) :: s3d(mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: ds3d(3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: et(0:maxl2, 0:maxl, 0:maxl, 3), ed(0:maxl2, 0:maxl, 0:maxl, 3)
   real(wp) :: dd(0:maxl2, 0:maxl, 0:maxl, 3)

   s3d(:, :) = 0.0_wp
   ds3d(:, :, :) = 0.0_wp

   do ip = 1, cgtoi%nprim
      do jp = 1, cgtoj%nprim
         eab = cgtoi%alpha(ip) + cgtoj%alpha(jp)
         oab = 1.0_wp/eab
         est = cgtoi%alpha(ip) * cgtoj%alpha(jp) * r2 * oab
         if (est > intcut) cycle
         pre = exp(-est) * sqrtpi3*sqrt(oab)**3
         rpi = -vec * cgtoj%alpha(jp) * oab
         rpj = +vec * cgtoi%alpha(ip) * oab
         cc = cgtoi%coeff(ip) * cgtoj%coeff(jp) * pre
         call e_function(cgtoj%ang+1, cgtoi%ang+1, 0.5_wp*oab, rpj, rpi, et)
         call e_derivative(et, cgtoi%alpha(ip), cgtoj%ang, cgtoi%ang, ed)

         do mli = 1, mlao(cgtoi%ang)
            do mlj = 1, mlao(cgtoj%ang)
               li = lx(:, mli+lmap(cgtoi%ang))
               lj = lx(:, mlj+lmap(cgtoj%ang))
               e0 = [et(0, lj(1), li(1), 1), et(0, lj(2), li(2), 2), et(0, lj(3), li(3), 3)]
               d0 = [ed(0, lj(1), li(1), 1), ed(0, lj(2), li(2), 2), ed(0, lj(3), li(3), 3)]
               s3d(mlj, mli) = s3d(mlj, mli) + cc * product(e0)

               grad = [d0(1)*e0(2)*e0(3), e0(1)*d0(2)*e0(3), e0(1)*e0(2)*d0(3)]
               ds3d(:, mlj, mli) = ds3d(:, mlj, mli) + cc*grad
            end do
         end do
      end do
   end do

   call transform0(cgtoj%ang, cgtoi%ang, s3d, overlap)
   call transform1(cgtoj%ang, cgtoi%ang, ds3d, doverlap)

end subroutine overlap_grad_cgto


!> Evaluate overlap for a molecular structure
subroutine get_overlap_lat(mol, trans, cutoff, bas, overlap)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Realspace cutoff
   real(wp), intent(in) :: cutoff
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Overlap matrix
   real(wp), intent(out) :: overlap(:, :)

   integer :: iat, jat, izp, jzp, itr, is, js
   integer :: ish, jsh, ii, jj, iao, jao, nao
   real(wp) :: r2, vec(3), cutoff2
   real(wp), allocatable :: stmp(:)

   overlap(:, :) = 0.0_wp

   allocate(stmp(msao(bas%maxl)**2))
   cutoff2 = cutoff**2

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(mol, bas, trans, cutoff2, overlap) private(r2, vec, stmp) &
   !$omp private(iat, jat, izp, jzp, itr, is, js, ish, jsh, ii, jj, iao, jao, nao)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      do jat = 1, mol%nat
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         do itr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            if (r2 > cutoff2) cycle
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is+ish)
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js+jsh)
                  call overlap_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                     & r2, vec, bas%intcut, stmp)

                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  !$omp simd collapse(2)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                           & + stmp(jao + nao*(iao-1))
                     end do
                  end do

               end do
            end do

         end do
      end do
   end do

end subroutine get_overlap_lat


end module tblite_integral_overlap
