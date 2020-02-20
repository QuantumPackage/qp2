BEGIN_PROVIDER [ double precision, slater_bragg_radii, (0:100)]
 implicit none
 BEGIN_DOC
 ! atomic radii in Angstrom defined in table I of JCP 41, 3199 (1964) Slater
 ! execpt for the Hydrogen atom where we took the value of Becke (1988, JCP)
 END_DOC

 slater_bragg_radii = 0.d0

 slater_bragg_radii(1)  = 0.35d0
 slater_bragg_radii(2)  = 0.35d0

 slater_bragg_radii(3)  = 1.45d0
 slater_bragg_radii(4)  = 1.05d0

 slater_bragg_radii(5)  = 0.85d0
 slater_bragg_radii(6)  = 0.70d0
 slater_bragg_radii(7)  = 0.65d0
 slater_bragg_radii(8)  = 0.60d0
 slater_bragg_radii(9)  = 0.50d0
 slater_bragg_radii(10) = 0.45d0

 slater_bragg_radii(11)  = 1.80d0
 slater_bragg_radii(12)  = 1.70d0

 slater_bragg_radii(13)  = 1.50d0
 slater_bragg_radii(14)  = 1.25d0
 slater_bragg_radii(15)  = 1.10d0
 slater_bragg_radii(16)  = 1.00d0
 slater_bragg_radii(17)  = 1.00d0
 slater_bragg_radii(18)  = 1.00d0

 slater_bragg_radii(19)  = 2.20d0
 slater_bragg_radii(20)  = 1.80d0


 slater_bragg_radii(21)  = 1.60d0
 slater_bragg_radii(22)  = 1.40d0
 slater_bragg_radii(23)  = 1.34d0
 slater_bragg_radii(24)  = 1.40d0
 slater_bragg_radii(25)  = 1.40d0
 slater_bragg_radii(26)  = 1.40d0
 slater_bragg_radii(27)  = 1.35d0
 slater_bragg_radii(28)  = 1.35d0
 slater_bragg_radii(29)  = 1.35d0
 slater_bragg_radii(30)  = 1.35d0

 slater_bragg_radii(31)  = 1.30d0
 slater_bragg_radii(32)  = 1.25d0
 slater_bragg_radii(33)  = 1.15d0
 slater_bragg_radii(34)  = 1.15d0
 slater_bragg_radii(35)  = 1.15d0
 slater_bragg_radii(36)  = 1.10d0

 slater_bragg_radii(37)  = 2.35d0
 slater_bragg_radii(38)  = 2.00d0
 slater_bragg_radii(39)  = 1.80d0
 slater_bragg_radii(40)  = 1.55d0
 slater_bragg_radii(41)  = 1.45d0
 slater_bragg_radii(42)  = 1.45d0
 slater_bragg_radii(43)  = 1.35d0
 slater_bragg_radii(44)  = 1.30d0
 slater_bragg_radii(45)  = 1.35d0
 slater_bragg_radii(46)  = 1.40d0
 slater_bragg_radii(47)  = 1.60d0
 slater_bragg_radii(48)  = 1.55d0
 slater_bragg_radii(49)  = 1.55d0
 slater_bragg_radii(50)  = 1.45d0
 slater_bragg_radii(51)  = 1.45d0
 slater_bragg_radii(52)  = 1.40d0
 slater_bragg_radii(53)  = 1.40d0
 slater_bragg_radii(54)  = 1.40d0
 slater_bragg_radii(55)  = 2.60d0
 slater_bragg_radii(56)  = 2.15d0
 slater_bragg_radii(57)  = 1.95d0
 slater_bragg_radii(58)  = 1.85d0
 slater_bragg_radii(59)  = 1.85d0
 slater_bragg_radii(60)  = 1.85d0
 slater_bragg_radii(61)  = 1.85d0
 slater_bragg_radii(62)  = 1.85d0
 slater_bragg_radii(63)  = 1.85d0
 slater_bragg_radii(64)  = 1.80d0
 slater_bragg_radii(65)  = 1.75d0
 slater_bragg_radii(66)  = 1.75d0
 slater_bragg_radii(67)  = 1.75d0
 slater_bragg_radii(68)  = 1.75d0
 slater_bragg_radii(69)  = 1.75d0
 slater_bragg_radii(70)  = 1.75d0
 slater_bragg_radii(71)  = 1.75d0
 slater_bragg_radii(72)  = 1.55d0
 slater_bragg_radii(73)  = 1.45d0
 slater_bragg_radii(74)  = 1.35d0
 slater_bragg_radii(75)  = 1.30d0
 slater_bragg_radii(76)  = 1.30d0
 slater_bragg_radii(77)  = 1.35d0
 slater_bragg_radii(78)  = 1.35d0
 slater_bragg_radii(79)  = 1.35d0
 slater_bragg_radii(80)  = 1.50d0
 slater_bragg_radii(81)  = 1.90d0
 slater_bragg_radii(82)  = 1.75d0
 slater_bragg_radii(83)  = 1.60d0
 slater_bragg_radii(84)  = 1.90d0
 slater_bragg_radii(85)  = 1.50d0
 slater_bragg_radii(86)  = 1.50d0

END_PROVIDER

BEGIN_PROVIDER [double precision, slater_bragg_radii_ua, (0:100)]
 implicit none
 integer :: i
 do i = 0, 100
  slater_bragg_radii_ua(i) = slater_bragg_radii(i) * 1.889725989d0
 enddo
END_PROVIDER

BEGIN_PROVIDER [double precision, slater_bragg_radii_per_atom, (nucl_num)]
 implicit none
 integer :: i
 do i = 1, nucl_num
  slater_bragg_radii_per_atom(i) = slater_bragg_radii(max(1,int(nucl_charge(i))))
 enddo
END_PROVIDER

BEGIN_PROVIDER [double precision, slater_bragg_radii_per_atom_ua, (nucl_num)]
 implicit none
 integer :: i
 do i = 1, nucl_num
  slater_bragg_radii_per_atom_ua(i) = slater_bragg_radii_ua(max(1,int(nucl_charge(i))))
 enddo
END_PROVIDER

BEGIN_PROVIDER [double precision, slater_bragg_type_inter_distance, (nucl_num, nucl_num)]
 implicit none
 integer :: i,j
 double precision :: xhi_tmp,u_ij
 slater_bragg_type_inter_distance = 0.d0
 do i = 1, nucl_num
  do j = i+1, nucl_num
   xhi_tmp = slater_bragg_radii_per_atom(i) / slater_bragg_radii_per_atom(j)
   u_ij = (xhi_tmp - 1.d0 ) / (xhi_tmp +1.d0)
   slater_bragg_type_inter_distance(i,j) = u_ij  / (u_ij * u_ij - 1.d0)
  enddo
 enddo
END_PROVIDER

BEGIN_PROVIDER [double precision, slater_bragg_type_inter_distance_ua, (nucl_num, nucl_num)]
 implicit none
 integer :: i,j
 double precision :: xhi_tmp,u_ij
 slater_bragg_type_inter_distance_ua = 0.d0
 do i = 1, nucl_num
  do j = i+1, nucl_num
   xhi_tmp = slater_bragg_radii_per_atom_ua(i) / slater_bragg_radii_per_atom_ua(j)
   u_ij = (xhi_tmp - 1.d0 ) / (xhi_tmp +1.d0)
   slater_bragg_type_inter_distance_ua(i,j) = u_ij  / (u_ij * u_ij - 1.d0)
   if(slater_bragg_type_inter_distance_ua(i,j).gt.0.5d0)then
    slater_bragg_type_inter_distance_ua(i,j) = 0.5d0
   else if( slater_bragg_type_inter_distance_ua(i,j) .le.-0.5d0)then
    slater_bragg_type_inter_distance_ua(i,j) = -0.5d0
   endif
  enddo
 enddo
END_PROVIDER
