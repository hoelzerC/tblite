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

module test_integral_multipole
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use tblite_basis_type
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_cutoff, only : get_lattice_points
   use tblite_integral_dipole
   use tblite_integral_multipole
   implicit none
   private

   public :: collect_integral_multipole

   real(wp), parameter :: thr = 5e+6_wp*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains


!> Collect all exported unit tests
subroutine collect_integral_multipole(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("dipole-trans-ss", test_dipole_ss), &
      new_unittest("dipole-trans-pp", test_dipole_pp), &
      new_unittest("dipole-trans-dd", test_dipole_dd), &
      ! new_unittest("dipole-grad-ss", test_dipole_grad_ss), &
      new_unittest("multipole-grad-pp", test_multipole_grad_pp2), &
      new_unittest("multipole-grad-pp", test_multipole_grad_pp), &
      new_unittest("multipole-h2o-ss", test_multipole_ss_h2), &
      new_unittest("multipole-h2o-pp", test_multipole_pp_h2), &
      new_unittest("multipole-h2o-dd", test_multipole_dd_h2), &
      new_unittest("multipole-01-ss", test_multipole_ss_MB1643_01), &
      new_unittest("multipole-01-pp", test_multipole_pp_MB1643_01), &
      new_unittest("multipole-01-dd", test_multipole_dd_MB1643_01) &
      ]

end subroutine collect_integral_multipole


subroutine make_basis(bas, mol, ng)
   type(basis_type), intent(out) :: bas
   type(structure_type), intent(in) :: mol
   integer, intent(in) :: ng

   integer, parameter :: nsh(20) = [&
      & 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3]
   integer, parameter :: lsh(3, 20) = reshape([&
      & 0, 0, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0, &
      & 0, 1, 0,  0, 1, 0,  0, 1, 2,  0, 1, 0,  0, 1, 2,  0, 1, 2,  0, 1, 2, &
      & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 0,  0, 1, 2], &
      & shape(lsh))
   integer, parameter :: pqn(3, 20) = reshape([&
      & 1, 0, 0,  1, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0, &
      & 2, 2, 0,  2, 2, 0,  2, 2, 3,  3, 3, 0,  3, 3, 3,  3, 3, 3,  3, 3, 3, &
      & 3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  4, 4, 0,  4, 4, 3], &
      & shape(pqn))
   real(wp), parameter :: zeta(3, 20) = reshape([&
      & 1.230000_wp, 0.000000_wp, 0.000000_wp, 1.669667_wp, 1.500000_wp, 0.000000_wp, &
      & 0.750060_wp, 0.557848_wp, 0.000000_wp, 1.034720_wp, 0.949332_wp, 0.000000_wp, &
      & 1.479444_wp, 1.479805_wp, 0.000000_wp, 2.096432_wp, 1.800000_wp, 0.000000_wp, &
      & 2.339881_wp, 2.014332_wp, 0.000000_wp, 2.439742_wp, 2.137023_wp, 0.000000_wp, &
      & 2.416361_wp, 2.308399_wp, 0.000000_wp, 3.084104_wp, 2.312051_wp, 2.815609_wp, &
      & 0.763787_wp, 0.573553_wp, 0.000000_wp, 1.184203_wp, 0.717769_wp, 1.300000_wp, &
      & 1.352531_wp, 1.391201_wp, 1.000000_wp, 1.773917_wp, 1.718996_wp, 1.250000_wp, &
      & 1.816945_wp, 1.903247_wp, 1.167533_wp, 1.981333_wp, 2.025643_wp, 1.702555_wp, &
      & 2.485265_wp, 2.199650_wp, 2.476089_wp, 2.329679_wp, 2.149419_wp, 1.950531_wp, &
      & 0.875961_wp, 0.631694_wp, 0.000000_wp, 1.267130_wp, 0.786247_wp, 1.380000_wp],&
      & shape(zeta))

   integer :: isp, izp, ish, stat
   integer, allocatable :: nshell(:)
   type(cgto_type), allocatable :: cgto(:, :)

   nshell = nsh(mol%num)
   allocate(cgto(maxval(nshell), mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, nshell(isp)
         call slater_to_gauss(ng, pqn(ish, izp), lsh(ish, izp), zeta(ish, izp), &
            & cgto(ish, isp), .true., stat)
      end do
   end do

   call new_basis(bas, mol, nshell, cgto, 1.0_wp)

end subroutine make_basis


subroutine test_dipole_ss(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6
   integer :: stat, i
   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3), r2
   real(wp) :: overlap(1, 1), dipolei(3, 1, 1), dipolej(3, 1, 1)

   call slater_to_gauss(ng, 2, 0, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 1, 0, 1.0_wp, cgtoj, .true., stat)

   call random_number(vec)
   vec = vec - 0.5_wp
   r2 = sum(vec**2)

   call dipole_cgto(cgtoi, cgtoj, r2, vec, 100.0_wp, overlap, dipolej)

   vec(:) = -vec

   call dipole_cgto(cgtoj, cgtoi, r2, vec, 100.0_wp, overlap, dipolei)

   do i = 1, 3
      call check(error, dipolei(i, 1, 1) + vec(i) * overlap(1, 1), dipolej(i, 1, 1), thr=thr)
      if (allocated(error)) return
   end do

end subroutine test_dipole_ss


subroutine test_dipole_pp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6, n = 2

   integer :: stat, i, j, k
   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3), r2
   real(wp) :: overlap(3, 3), dipolei(3, 3, 3), dipolej(3, 3, 3)

   call slater_to_gauss(ng, 2, 1, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 3, 1, 1.0_wp, cgtoj, .true., stat)

   call random_number(vec)
   vec = vec - 0.5_wp
   r2 = sum(vec**2)

   call dipole_cgto(cgtoi, cgtoj, r2, vec, 100.0_wp, overlap, dipolej)

   vec(:) = -vec

   call dipole_cgto(cgtoj, cgtoi, r2, vec, 100.0_wp, overlap, dipolei)

   do i = 1, 3
      do j = 1, 3
         do k = 1, 3
            call check(error, dipolei(k, j, i) + vec(k) * overlap(j, i), dipolej(k, i, j), thr=thr)
            if (allocated(error)) return
         end do
      end do
   end do

end subroutine test_dipole_pp


subroutine test_dipole_dd(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6, n = 2

   integer :: stat, i, j, k
   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3), r2
   real(wp) :: overlap(5, 5), dipolei(3, 5, 5), dipolej(3, 5, 5)

   call slater_to_gauss(ng, 3, 2, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 4, 2, 1.0_wp, cgtoj, .true., stat)

   call random_number(vec)
   vec = vec - 0.5_wp
   r2 = sum(vec**2)

   call dipole_cgto(cgtoi, cgtoj, r2, vec, 100.0_wp, overlap, dipolej)

   vec(:) = -vec

   call dipole_cgto(cgtoj, cgtoi, r2, vec, 100.0_wp, overlap, dipolei)

   do i = 1, 5
      do j = 1, 5
         do k = 1, 3
            call check(error, dipolei(k, j, i) + vec(k) * overlap(j, i), dipolej(k, i, j), thr=thr)
            if (allocated(error)) return
         end do
      end do
   end do

end subroutine test_dipole_dd




subroutine test_multipole_mol(error, cgtoi, cgtoj, vec, refovlp, refdp, refqp)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Input structure
   type(cgto_type), intent(in) :: cgtoi, cgtoj
   real(wp), intent(in) :: vec(:)
   !> Reference values for overlap, dipole and quadrupole integrals
   real(wp), intent(in) :: refovlp(:, :), refdp(:, :, :), refqp(:, :, :)

   integer, parameter :: ao(0:6) = [1, 3, 5, 7, 9, 11, 13]

   integer :: stat, i,j,k
   real(wp) :: r2
   real(wp), allocatable :: overlap(:, :), dipole(:, :, :), quadrupole(:, :, :)

   allocate(overlap(ao(cgtoj%ang), ao(cgtoi%ang)))
   allocate(dipole(3, ao(cgtoj%ang), ao(cgtoi%ang)))
   allocate(quadrupole(6, ao(cgtoj%ang), ao(cgtoi%ang)))

   r2 = sum(vec**2)

   call multipole_cgto(cgtoj, cgtoi, r2, vec, 100.0_wp, overlap, dipole, quadrupole)

   write(*,*) "test_multipole_mol"
   write(*,*) "vec"
   write(*,*) vec
   write(*,*) "overlap"
   write(*,*) overlap
   write(*,*) "dipole"
   write(*,*) dipole
   write(*,*) "quadrupole"
   write(*,*) quadrupole

   ! overlap
   do i = 1, size(overlap, 2)
      do j = 1, size(overlap, 1)
         call check(error, overlap(j, i), refovlp(j, i), thr=thr)
         if (allocated(error)) then
            print '(2es20.13)', overlap(j, i), refovlp(j, i), &
               & overlap(j, i) - refovlp(j, i)
            return
         end if
      end do
   end do
   ! dipole
   do i = 1, size(dipole, 3)
      do j = 1, size(dipole, 2)
         do k = 1, size(dipole, 1)
            call check(error, dipole(k, j, i), refdp(k, j, i), thr=thr)
            if (allocated(error)) then
               print '(2es20.13)', dipole(k, j, i), refdp(k, j, i), &
                  & dipole(k, j, i) - refdp(k, j, i)
               return
            end if
         end do
      end do
   end do
   ! quadrupole
   do i = 1, size(quadrupole, 3)
      do j = 1, size(quadrupole, 2)
         do k = 1, size(quadrupole, 1)
            call check(error, quadrupole(k, j, i), refqp(k, j, i), thr=thr)
            if (allocated(error)) then
               print '(2es20.13)', quadrupole(k, j, i), refqp(k, j, i), &
                  & quadrupole(k, j, i) - refqp(k, j, i)
               return
            end if
         end do
      end do
   end do

end subroutine test_multipole_mol


subroutine test_multipole_ss_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6
   integer :: stat
   real(wp) :: vec(3)
   type(structure_type) :: mol
   type(cgto_type) :: cgtoi, cgtoj

   ! references
   real(wp), parameter :: refovlp(1,1) = reshape([&
      &0.76295311418571510_wp],shape(refovlp))
   real(wp), parameter :: refdp(3, 1, 1) = reshape([&
      &0.0_wp, 0.0_wp, -0.70416920345374734_wp],shape(refdp))
   real(wp), parameter :: refqp(6, 1, 1) = reshape([&
      &-0.34407168144461275_wp, 0.0_wp, -0.34407168144461275_wp,& 
      &0.0_wp, 0.0_wp, 0.68814336288922595_wp],shape(refqp))
   
   call get_structure(mol, "MB16-43", "H2")
   vec(:) = mol%xyz(:, 2) - mol%xyz(:, 1)

   call slater_to_gauss(ng, 2, 0, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 1, 0, 1.0_wp, cgtoj, .true., stat)

   call test_multipole_mol(error, cgtoi, cgtoj, vec, refovlp, refdp, refqp)

end subroutine test_multipole_ss_h2


subroutine test_multipole_ss_MB1643_01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6
   integer :: stat
   real(wp) :: vec(3)
   type(structure_type) :: mol
   type(cgto_type) :: cgtoi, cgtoj

   ! references
   real(wp), parameter :: refovlp(1,1) = reshape([&
      &0.0351230934573249_wp],shape(refovlp))
   real(wp), parameter :: refdp(3, 1, 1) = reshape([&
      &-0.13703277985614506_wp, 0.0780384030559952_wp, 0.0555602485121848_wp],shape(refdp))
   real(wp), parameter :: refqp(6, 1, 1) = reshape([&
      &0.43947308532585705_wp, -0.49680825340407353_wp,-0.14997939849296427_wp,&
      &-0.35370777644218176_wp, 0.20143202269567789_wp, -0.28949368683289278_wp],shape(refqp))
   
   call get_structure(mol, "MB16-43", "01")
   vec(:) = mol%xyz(:, 2) - mol%xyz(:, 1)

   call slater_to_gauss(ng, 2, 0, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 1, 0, 1.0_wp, cgtoj, .true., stat)

   call test_multipole_mol(error, cgtoi, cgtoj, vec, refovlp, refdp, refqp)

end subroutine test_multipole_ss_MB1643_01


subroutine test_multipole_pp_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6
   integer :: stat
   real(wp) :: vec(3)
   type(structure_type) :: mol
   type(cgto_type) :: cgtoi, cgtoj

   ! references
   real(wp), parameter :: refovlp(3, 3) = reshape([&
      &0.81113650577472274_wp,        0.00000000000000000_wp,      0.00000000000000000_wp,&
      &0.00000000000000000_wp,        0.62061321111714240_wp,      0.00000000000000000_wp,& 
      &0.00000000000000000_wp,        0.00000000000000000_wp,      0.81113650577472274_wp],shape(refovlp))
   real(wp), parameter :: refdp(3, 3, 3) = reshape([&
      &0.00000000000000000_wp,       0.00000000000000000_wp,      -0.45805108559310459_wp,&
      &0.00000000000000000_wp,       0.68164325623828437_wp,       0.00000000000000000_wp,&
      &0.00000000000000000_wp,       0.00000000000000000_wp,       0.00000000000000000_wp,&
      &0.00000000000000000_wp,      -0.45805108559310459_wp,       0.00000000000000000_wp,&
      &0.00000000000000000_wp,       0.00000000000000000_wp,      -0.13060587954988973_wp,&
      &-0.45805108559310459_wp,      0.00000000000000000_wp,       0.00000000000000000_wp,& 
      &0.00000000000000000_wp,       0.00000000000000000_wp,       0.00000000000000000_wp,&
      &0.68164325623828437_wp,       0.00000000000000000_wp,       0.00000000000000000_wp,&
      &0.00000000000000000_wp,       0.00000000000000000_wp,      -0.45805108559310459_wp],shape(refdp))
   real(wp), parameter :: refqp(6, 3, 3) = reshape([&
      &-1.8851243601670538_wp,       0.0000000000000000_wp,        3.3195861328492393_wp,&
      &0.0000000000000000_wp,        0.0000000000000000_wp,       -1.4344617726821856_wp,&
      &0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,&
      &0.0000000000000000_wp,        2.0876348926441031_wp,        0.0000000000000000_wp,&
      &0.0000000000000000_wp,        2.6023552465081465_wp,        0.0000000000000000_wp,&
      &0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,&
      &0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,&
      &0.0000000000000000_wp,        3.0530178339930147_wp,        0.0000000000000000_wp,&
      &-1.8097587378293412_wp,       0.0000000000000000_wp,       -1.8097587378293412_wp,&
      &0.0000000000000000_wp,        0.0000000000000000_wp,        3.6195174756586814_wp,&
      &0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,&
      &3.0530178339930147_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,&
      &0.0000000000000000_wp,        2.6023552465081465_wp,        0.0000000000000000_wp,&
      &0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,&
      &0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,&
      &2.0876348926441031_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,&
      &3.3195861328492393_wp,        0.0000000000000000_wp,       -1.8851243601670538_wp,&
      &0.0000000000000000_wp,        0.0000000000000000_wp,       -1.4344617726821856_wp],shape(refqp))

   call get_structure(mol, "MB16-43", "H2")
   vec(:) = mol%xyz(:, 2) - mol%xyz(:, 1)

   call slater_to_gauss(ng, 2, 1, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 3, 1, 1.0_wp, cgtoj, .true., stat)

   call test_multipole_mol(error, cgtoi, cgtoj, vec, refovlp, refdp, refqp)

end subroutine test_multipole_pp_h2

subroutine test_multipole_pp_MB1643_01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6
   integer :: stat
   real(wp) :: vec(3)
   type(structure_type) :: mol
   type(cgto_type) :: cgtoi, cgtoj

   ! references   
   real(wp), parameter :: refovlp(3, 3) = reshape([&
      &-4.1455855514000860E-004_wp, -4.021387048375E-002_wp,      9.9182754014421420E-002_wp,&
      &-4.021387048375E-002_wp,      2.7438061136871601E-002_wp,  7.0614187955769625E-002_wp,&
      & 9.9182754014421420E-002_wp,  7.0614187955769625E-002_wp, -0.11809279820837945_wp],shape(refovlp))
   real(wp), parameter :: refdp(3, 3, 3) = reshape([&
      &0.00041034403144520_wp,      -0.03435490909856620_wp,      -0.00016637490961302_wp,&
      &0.10385031908209739_wp,      -0.14240942129477671_wp,       0.04072838623244988_wp,&
      &-0.17329979446068089_wp,      0.35123598976430753_wp,       0.10385031908209737_wp,&       
      &0.10385031908209738_wp,      -0.00016637490961302_wp,      -0.15906230014800446_wp,&       
      &-0.07151766009976784_wp,      0.04072838623244988_wp,       0.00470403790088692_wp,&
      &-0.12338258178396026_wp,      0.10385031908209737_wp,       0.27930798563309445_wp,&      
      &-0.37309048084113516_wp,      0.00041034403144520_wp,       0.10385031908209738_wp,&      
      &-0.26562562816912394_wp,      0.10385031908209739_wp,      -0.07151766009976784_wp,&  
      &0.36422424822833688_wp,      -0.17329979446068089_wp,      -0.12338258178396026_wp],shape(refdp))
   real(wp), parameter :: refqp(6, 3, 3) = reshape([&
      &-0.15770909098483668_wp,      0.025928944772843666_wp,      0.31624461599443976_wp,&
      &-0.000401000638767978_wp,    -0.010512948921785338_wp,     -0.15853552500960308_wp,&
      &-0.22896018743248364_wp,      0.55311483786829152_wp,      -0.23821931808435831_wp,&
      &-0.18543961238688714_wp,      0.24196759700786297_wp,       0.46717950551684195_wp,& 
      &-0.18625778961468331_wp,     -0.89796306001535808_wp,       0.58753976533016727_wp,& 
      &-0.31647395865265016_wp,      0.55311483786829152_wp,      -0.40128197571548374_wp,&  
      &-0.22896018743248364_wp,     -0.00040100063876797806_wp,    0.39222096268112078_wp,&
      &0.59201357149903289_wp,       0.023642263392206951_wp,     -0.16326077524863719_wp,&
      &0.00087154172855858_wp,      -0.18543961238688714_wp,      -0.21914863581657129_wp,&   
      &-0.11316426961185033_wp,      0.064445593913935839_wp,      0.21827709408801266_wp,&   
      &-0.13260816050908453_wp,     -0.31647395865265016_wp,      -0.15407220451244608_wp,&   
      &-0.67876931649499428_wp,      0.59201357149903289_wp,       0.28668036502153060_wp,& 
      &1.3686485781571567_wp,        0.024468697416973349_wp,     -0.96736660244167294_wp,&   
      &-0.63169409903538976_wp,     -0.00040100063876797806_wp,   -0.40128197571548374_wp,&     
      &0.97442351650503423_wp,      -0.63169409903538976_wp,      -0.15407220451244619_wp,&  
      &0.46198777455299284_wp,      -0.18543961238688714_wp,      -0.82035131199258804_wp,&    
      &-0.41596035842396395_wp,      0.80703775759983309_wp,       0.09838330745180013_wp,&
      &0.57457888187165307_wp,      -0.31647395865265016_wp,       0.31757705097216393_wp,&     
      &-0.040213870483757_wp,        0.040213870483757_wp,        -0.080427740967513_wp],shape(refqp))

   call get_structure(mol, "MB16-43", "01")
   vec(:) = mol%xyz(:, 2) - mol%xyz(:, 1)

   call slater_to_gauss(ng, 2, 1, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 3, 1, 1.0_wp, cgtoj, .true., stat)

   call test_multipole_mol(error, cgtoi, cgtoj, vec, refovlp, refdp, refqp)

end subroutine test_multipole_pp_MB1643_01 

subroutine test_multipole_dd_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6
   integer :: stat
   real(wp) :: vec(3)
   type(structure_type) :: mol
   type(cgto_type) :: cgtoi, cgtoj

   ! references
   real(wp), parameter :: refovlp(5, 5) = reshape([&
      &0.63641740867659413_wp,        0.0000000000000000_wp,       0.0000000000000000_wp,       -4.8074067159589095E-017_wp,   0.0000000000000000_wp,&
      &0.0000000000000000_wp,         0.67912031280053986_wp,      0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,&
      &0.0000000000000000_wp,         0.0000000000000000_wp,       0.67912031280053986_wp,       0.0000000000000000_wp,        0.0000000000000000_wp,&
      &-4.8074067159589095E-017_wp,   0.0000000000000000_wp,       0.0000000000000000_wp,        0.84655682077983663_wp,       0.0000000000000000_wp,&
      &0.0000000000000000_wp,         0.0000000000000000_wp,       0.0000000000000000_wp,        0.0000000000000000_wp,        0.84655682077983707_wp&
      &],shape(refovlp))
   real(wp), parameter :: refdp(3, 5, 5) = reshape([&
      &0.0000000000000000_wp,        0.0000000000000000_wp,      -0.18438637930134960_wp,     -0.92390316157108776_wp,       0.0000000000000000_wp,        0.0000000000000000_wp,       0.0000000000000000_wp,      -0.92390316157108776_wp,       0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        1.0004790555510639_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,      -0.24225382619666419_wp,       0.0000000000000000_wp,        0.0000000000000000_wp,       0.0000000000000000_wp,      -0.51046526539807591_wp,&
      &0.0000000000000000_wp,        0.0000000000000000_wp,       0.0000000000000000_wp,      -0.51046526539807568_wp,       0.0000000000000000_wp,        0.0000000000000000_wp,       1.0004790555510639_wp,       0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,       -0.24225382619666419_wp,       0.0000000000000000_wp,        0.51046526539807591_wp,       0.0000000000000000_wp,       -0.51046526539807568_wp,       0.0000000000000000_wp,        0.0000000000000000_wp,       0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,       0.67899669545898877_wp,      0.0000000000000000_wp,&
      &0.0000000000000000_wp,        0.0000000000000000_wp,      -0.67899669545898877_wp,      0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,      -0.51046526539807591_wp,      0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.67899669545898877_wp,       0.0000000000000000_wp,        0.67899669545898877_wp,       0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,       0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,       0.0000000000000000_wp,      -0.51046526539807568_wp& 
      &],shape(refdp))
   real(wp), parameter :: refqp(6, 5, 5) = reshape([&
      &-2.2806251375435522_wp,       0.0000000000000000_wp,       -2.2806251375435540_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        4.5612502750871062_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,       0.0000000000000000_wp,        2.4332011245606582_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        2.4332011245606582_wp,        0.0000000000000000_wp,       -3.5284839312466305_wp,        0.0000000000000000_wp,        3.5284839312466305_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,&
      &-3.5284839312466332_wp,       0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,       0.89613760533488729_wp,       0.0000000000000000_wp,        0.0000000000000000_wp,        1.3532377313228956_wp,        0.0000000000000000_wp,       -4.1628563798380895_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        2.8096186485151922_wp,        0.0000000000000000_wp,        2.7580470555804926_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,&
      &0.0000000000000000_wp,        3.8338974898794427_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,       3.8338974898794427_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.89613760533488729_wp,       0.0000000000000000_wp,        0.0000000000000000_wp,        2.7580470555804926_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,       -4.1628563798380895_wp,        0.0000000000000000_wp,        1.3532377313228956_wp,&
      &0.0000000000000000_wp,        0.0000000000000000_wp,        2.8096186485151922_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,       -3.8338974898794427_wp,       0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        3.8338974898794427_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,       -3.3234121129047169_wp,        0.0000000000000000_wp,        3.3234121129047169_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        2.7580470555804926_wp,&
      &0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,       -2.7580470555804926_wp,        0.0000000000000000_wp,       2.0371711475365100_wp,        0.0000000000000000_wp,        2.0371711475365100_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,       -4.0743422950730235_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,       -3.3234121129047183_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,&
      &0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        2.7580470555804926_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,       0.0000000000000000_wp,        0.0000000000000000_wp,        2.7580470555804926_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,        2.0371711475365135_wp,        0.0000000000000000_wp,        2.0371711475365135_wp,        0.0000000000000000_wp,        0.0000000000000000_wp,       -4.0743422950730253_wp&
      &],shape(refqp))

   call get_structure(mol, "MB16-43", "H2")
   vec(:) = mol%xyz(:, 2) - mol%xyz(:, 1)

   call slater_to_gauss(ng, 3, 2, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 4, 2, 1.0_wp, cgtoj, .true., stat)

   call test_multipole_mol(error, cgtoi, cgtoj, vec, refovlp, refdp, refqp)

end subroutine test_multipole_dd_h2

subroutine test_multipole_dd_MB1643_01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6
   integer :: stat
   real(wp) :: vec(3)
   type(structure_type) :: mol
   type(cgto_type) :: cgtoi, cgtoj

   ! references   
   real(wp), parameter :: refovlp(5, 5) = reshape([&
      &-1.0617158832898710E-002_wp,   0.13680142879883883_wp,      -7.7906651608812730E-002_wp,   2.0332091769322991E-002_wp,  -3.4272994947879387E-002_wp,&
      &0.13680142879883883_wp,       -6.7939066627836553E-002_wp,   6.3983559997925391E-002_wp,  -4.5812616909048434E-003_wp,   0.11160497018192761_wp,&
      &-7.7906651608812744E-002_wp,   6.3983559997925404E-002_wp,   7.9760738834713674E-003_wp,   0.10861361329288158_wp,      -6.7145920029129549E-004_wp,&
      &2.0332091769322984E-002_wp,   -4.5812616909047810E-003_wp,   0.10861361329288155_wp,      -0.14047505829930806_wp,      -0.13715867112024530_wp,&
      &-3.4272994947879359E-002_wp,   0.11160497018192761_wp,      -6.7145920029129517E-004_wp,  -0.13715867112024527_wp,       9.3599039219173855E-003_wp&
      &],shape(refovlp))
   real(wp), parameter :: refdp(3, 5, 5) = reshape([&
      &1.9417196142526572E-002_wp,  -1.1057843096947462E-002_wp,  -7.3436172251956033E-002_wp, -0.31555071871180973_wp,       0.31936266375870759_wp,       0.28121337436826610_wp,       0.31936266375870764_wp,        6.3366326449384433E-002_wp, -0.16014739449002807_wp,      -0.11482016532534997_wp,      -0.18246619699583916_wp,       0.15130731202790215_wp,      -4.9345192116617559E-002_wp,  -0.11893578868557436_wp,      -0.25505269204685294_wp,      -0.50766434183678966_wp,       0.12789087849018740_wp,        7.6323875540349695E-002_wp,   0.22041568018126589_wp,      -0.10528122904628488_wp,       -2.4817574975612947E-002_wp,  -0.26354197526998441_wp,       -2.1921465053334664E-002_wp,   2.6673704500472774E-002_wp,  -3.9028165930426101E-002_wp,&
      &-0.15034488181912573_wp,      1.0993524284060807E-002_wp,  -0.34518467668827524_wp,      9.6250373087492635E-003_wp,   0.41994231774651591_wp,       0.12789087849018740_wp,      -0.35592451974792833_wp,       -4.3465464018654199E-002_wp, -0.10528122904628488_wp,       0.25597901916826526_wp,       2.6673704500472770E-002_wp,  -2.1921465053334664E-002_wp,  -7.7586697333960376E-003_wp,   6.8302111928364760E-003_wp,  -0.15034488181912570_wp,       0.27020146115134669_wp,       0.44117017970244587_wp,        9.6250373087492635E-003_wp,  -4.4494913050023795E-002_wp,   2.6281878854779873E-002_wp,    1.7795540466728044E-002_wp,   0.21214364356411824_wp,       -7.5869911044430666E-002_wp,   0.13816892630854855_wp,       0.18129397848186929_wp,&
      &-1.7781669493499644E-002_wp, -0.48198374212012074_wp,       9.3004368912372209E-002_wp, -0.12789165276829417_wp,       0.38383751299645374_wp,      -0.29190412710359492_wp,      -0.15461481354088641_wp,        0.42146849078157678_wp,     -0.25094816632110367_wp,      -0.15034488181912575_wp,       0.18783092561150200_wp,       2.4896764814235306E-002_wp,   0.12789087849018740_wp,      -0.34518467668827524_wp,       0.38356337588956824_wp,      -0.10528122904628488_wp,      -9.9191167941207279E-003_wp,   -4.4494913050023774E-002_wp,  -2.1921465053334664E-002_wp,   0.46530677645951601_wp,       -0.17396962863542861_wp,      -0.15034488181912573_wp,       -1.8663458624633353E-002_wp,   8.3942206820167842E-002_wp,   9.6250373087492635E-003_wp&
      &],shape(refdp))
   real(wp), parameter :: refqp(6, 5, 5) = reshape([&
      &-0.12524301014530043_wp,       0.47188281014078509_wp,       0.43463596248251307_wp,      -0.48556048725992618_wp,       0.27652044315699054_wp,      -0.30939295233721265_wp,      -0.18469908533624224_wp,       -1.0130370053802371_wp,       0.35424156618832114_wp,      -0.61575957160181394_wp,       0.45193559725219212_wp,      -0.16954248085207912_wp,       -1.3920622192802321_wp,       -0.31131380991308832_wp,        1.2955099602533231_wp,        0.45193559725219212_wp,      -7.9547545128537678E-002_wp,   9.6552259026908782E-002_wp,   0.51179830588892750_wp,       0.51798131348199239_wp,      -1.1831999587942268_wp,      -0.18287288935824733_wp,       0.28530391879075001_wp,        0.67140165290529952_wp,       1.0838587873565708_wp,&
      &0.28164565838342037_wp,        4.7896160392585374E-002_wp,   0.48579501526449165_wp,      -0.16918265266110888_wp,      -1.1317549477491560_wp,        1.2409985842459266_wp,       -0.68148167800986625_wp,       -0.68032162819119768_wp,     -0.27024808576081571_wp,      -0.34508473543005969_wp,      -0.56067695605472867_wp,       0.29484742383724960_wp,        0.49589892799551827_wp,       0.41617048500354248_wp,        9.5684792980506861E-002_wp,  -0.17308294434260629_wp,      -0.71101790884079197_wp,       0.65733252658569330_wp,       0.16167739182947408_wp,      -1.1791828040455485_wp,       -0.29078219400396482_wp,     -4.1608591827600142E-002_wp,   0.52185027745985479_wp,        0.81228583820728639_wp,       0.74533128511862112_wp,&
      &-0.60112241441034886_wp,      -0.26336204276935515_wp,      -0.85420546458367608_wp,      -0.21116342379693750_wp,       0.38965291789913503_wp,      -4.1169363101337828E-002_wp,  -0.95060068356700045_wp,       -1.9957867107661180_wp,      -7.4413696175467262E-003_wp,   0.56094776566786542_wp,       0.25029078337123267_wp,       0.95583910263305172_wp,       -0.56958908795747898_wp,      -0.34508473543005969_wp,       -0.67968382486832868_wp,       0.31929830458624608_wp,      -0.63804922138063713_wp,      -0.60679412972884650_wp,       0.11619894392078223_wp,      -0.17308294434260629_wp,       0.16506728487512451_wp,      0.52185027745985502_wp,      -0.16076348658699605_wp,       -4.5234654669392874E-002_wp,   0.25261707977366626_wp,&
      &-4.1608591827600142E-002_wp,   0.14228724432296341_wp,      -9.1853593186670246E-002_wp,  -0.58980451703816228_wp,      -0.32465613012271877_wp,       0.26036392840877332_wp,      -0.85420546458367608_wp,        1.6178240438181257_wp,       0.32944058862938963_wp,      -0.35903499518527682_wp,       4.6070747601058004E-002_wp,   0.55439029529748929_wp,       -7.4413696175467262E-003_wp,  -0.28737578159876598_wp,       -0.19535530011221244_wp,      -0.62571551788890600_wp,      -0.78079381597841013_wp,       0.94896179344918496_wp,       0.96066010621856179_wp,       0.40993993566435583_wp,      -0.32324627556027885_wp,     -0.87963607166079505_wp,      -0.19522426985619717_wp,        1.1907936018957344_wp,       -6.0661530583146377E-002_wp,&
      &-0.37712408265070962_wp,      -0.31115753023493914_wp,       2.2548343275174783_wp,       -0.37465318334171960_wp,      -0.70316380619723051_wp,       0.91825766531562103_wp,      -0.20669028318033117_wp,       -1.5516705213202475_wp,      -0.86436861574572021_wp,      -0.17447814628470540_wp,      -0.11423908864504706_wp,       0.70670412630881774_wp,       -0.42508422272346769_wp,       0.97860770439076772_wp,       -0.45996846903312449_wp,       0.72733400834022843_wp,      -0.13814483251492415_wp,       0.67012680073306763_wp,      -0.45671264524138300_wp,       0.59811330154804865_wp,      -1.0532354834391493_wp,       6.5616618424618522E-002_wp,   0.50835214851767063_wp,       -0.68148167800986614_wp,       0.95583910263305150_wp,&
      &0.54488333492147856_wp,        0.38965291789913514_wp,      -1.8151507348690832_wp,        0.78010145761152960_wp,       0.49589892799551827_wp,      -0.60679412972884650_wp,      -1.1697543755106645_wp,        -0.19491607138595007_wp,      0.14424741161779658_wp,       0.55439029529748929_wp,       0.16167739182947408_wp,      -4.5234654669392874E-002_wp,   -0.35947422391153910_wp,      -0.47353479705872958_wp,        0.71928592772139643_wp,      -0.12457850448931884_wp,       0.74533128511862112_wp,      -0.32465613012271877_wp,       0.59811330154804843_wp,       0.56472581165647129_wp,      -0.41932804856994432_wp,     -0.88994220502398558_wp,      -4.1169363101337828E-002_wp,    4.6070747601058004E-002_wp,   0.32521639336751434_wp&
      &],shape(refqp))

   call get_structure(mol, "MB16-43", "01")
   vec(:) = mol%xyz(:, 2) - mol%xyz(:, 1)

   call slater_to_gauss(ng, 3, 2, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 4, 2, 1.0_wp, cgtoj, .true., stat)

   call test_multipole_mol(error, cgtoi, cgtoj, vec, refovlp, refdp, refqp)

end subroutine test_multipole_dd_MB1643_01


end module test_integral_multipole
