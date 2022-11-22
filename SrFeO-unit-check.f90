program SrFeO_unit_check
implicit none

integer, parameter :: dp=selected_real_kind(15, 300)
real(kind=dp), parameter :: pi=4.0_dp*atan(1.0_dp)

! prim_cell is the primitive unit cell. cell and cell_ortho are the regular unit cell we are using in
! either non-orthogonal or orthogonal bases. n_cell and n_cell_ortho are the same but contain 
! 27 copies of the basis cell to account for neighbour interactions. distances holds the distances

real(kind=dp), dimension(:, :), allocatable :: prim_cell, cell, cell_ortho
real(kind=dp) :: a, b, c, t, id
integer :: N_prim, N, istat, i, j, l, counter, x_prim_mul, y_prim_mul, x_unit, y_unit, prim_in_unit

! This is a cool hcp sim http://lampx.tugraz.at/~hadley/ss1/crystalstructure/structures/hcp/hcp.php

! Define lattice constants in orthogonal and non-orthogonal bases
! a is the characteristic length (nearest neighbour distance) and b and c are the multiples of a
! in their respective directions. Here a and b are for room temp SrFe12O19 from literature
a=5.88121_dp

b=1.0_dp
c=23.023_dp/a

! Define number of atoms in primitive unit cell n
N_prim=11

! Define number of primitive unit cells in the unit cell we will use
! Will probably be something like 4, 2 in x and y, to preserve translational symmetry when converting to an orthogonal basis
! the number of primitive cells used in the unit cell in either the x or the y directions are held in x_prim_mul and y_prim_mul
! This program is for materials/structures using the P63/mmc space group, where the z/c basis vector is the same in the 
! non-orthogonal and orthogonal bases, so there doesn't need to be multiple primitive cells in the z direction

x_prim_mul=2
y_prim_mul=3

prim_in_unit=x_prim_mul*y_prim_mul

N=N_prim*prim_in_unit

! Create primitive cell
allocate (prim_cell(6, N_prim), stat=istat)
if(istat/=0) stop 'Error allocating prim_cell'

! Primitive unit cell is based of data from table 1: https://www.sciencedirect.com/science/article/pii/030488539290253K
! Fe2 is 2b not 4e. The functions at the bottom of the program can be used as a key for the element, material and wyckoff ids

! a, b, c, element_id, material_id, wyckoff_id

prim_cell(:, 1) = [1.0_dp/3.0_dp, 2.0_dp/3.0_dp, 3.0_dp/4.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
prim_cell(:, 2) = [0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 2.0_dp, 2.0_dp]
prim_cell(:, 3) = [0.0_dp, 0.0_dp, 1.0_dp/4.0_dp, 3.0_dp, 2.0_dp, 3.0_dp]
prim_cell(:, 4) = [1.0_dp/3.0_dp, 2.0_dp/3.0_dp, 0.02728_dp, 4.0_dp, 2.0_dp, 4.0_dp]
prim_cell(:, 5) = [1.0_dp/3.0_dp, 2.0_dp/3.0_dp, 0.19080_dp, 5.0_dp, 2.0_dp, 5.0_dp]
prim_cell(:, 6) = [0.16878_dp, 2*(0.16878_dp), -0.10921_dp, 6.0_dp, 2.0_dp, 6.0_dp]
prim_cell(:, 7) = [0.0_dp, 0.0_dp, 0.1516_dp, 7.0_dp, 3.0_dp, 7.0_dp]
prim_cell(:, 8) = [1.0_dp/3.0_dp, 2.0_dp/3.0_dp, -0.0552_dp, 8.0_dp, 3.0_dp, 8.0_dp]
prim_cell(:, 9) = [0.1828_dp, 2*(0.1828_dp), 1.0_dp/4.0_dp, 9.0_dp, 3.0_dp, 9.0_dp]
prim_cell(:, 10) = [0.1563_dp, 2*(0.1563_dp), 0.0524_dp, 10.0_dp, 3.0_dp, 6.0_dp]
prim_cell(:, 11) = [0.5051_dp, 2*(0.5051_dp), 0.1510_dp, 11.0_dp, 3.0_dp, 6.0_dp]

! Implement periodic boundary conditions, if an atom has fraction coodinates that are not 0 =< a, b, c =< 1, add or subtract 
! multiples of 1 until all of the atoms satisfy this condition

! Do while there are still atoms outside of pbc
counter=1
do while (counter > 0)

    ! Add or subtract 1.0_dp from the atoms out of pbc
    do i=1, 3, 1
        do j=1, N_prim, 1
            if (prim_cell(i, j) < 0.0_dp) then
                prim_cell(i, j)=prim_cell(i, j)+1.0_dp
            elseif (prim_cell(i, j) > 1.0_dp) then
                prim_cell(i, j)=prim_cell(i, j)-1.0_dp
            end if
        end do
    end do

    ! Assume there are no more atoms out of pbc
    counter=0

    ! Check if there are any atoms now out of pbc
    do i=1, 3, 1
        do j=1, N_prim, 1
            if (prim_cell(i, j) < 0.0_dp) then
                counter=counter+1
            elseif (prim_cell(i, j) > 1.0_dp) then
                counter=counter+1
            end if
        end do
    end do
end do

! Create unit cell using multiples of the primitive cell at the origin
allocate (cell(6, N), stat=istat)
if(istat/=0) stop 'Error allocating cell'

counter=1
do l=1, N_prim, 1
    do x_unit=0, x_prim_mul-1, 1
        do y_unit=0, y_prim_mul-1, 1
            cell(:, counter)=[prim_cell(1, l)+real(x_unit, kind=dp), prim_cell(2, l)+real(y_unit, kind=dp), prim_cell(3, l)&
            &, prim_cell(4, l), prim_cell(5, l), prim_cell(6, l)]
            counter=counter+1
        end do
    end do
end do

! Put this unit cell into an orthogonal basis
allocate (cell_ortho(6, N), stat=istat)
if(istat/=0) stop 'Error allocating cell_ortho'

cell_ortho=cell

! This also loops the values that are are further than the bounds of the unit cell under the new basis to the opposite side
! The values at which they loop around are found from how many primitive cell are in the x direction of the unit cell
! As the y fractions will stay the same
do i=1, N, 1
    t=cell(1, i) + cell(2, i)*sin(pi*7.0_dp/6.0_dp)
    if (t > real(x_prim_mul, kind=dp)) then
        cell_ortho(1, i) = t-real(x_prim_mul, kind=dp)
    elseif (t < 0.0_dp) then
        cell_ortho(1, i) = t+real(x_prim_mul, kind=dp)
    else
        cell_ortho(1, i) = t
    end if
end do

! Printing Fe1 positions using the element_id function
do i=1, N, 1
    id=cell_ortho(4, i)
    if (element_id(id) == 'Fe1') print *, cell_ortho(1:3,i)
end do

! Deallocate arrays
deallocate (prim_cell, stat=istat)
if(istat/=0) stop 'Error deallocating prim_cell'
deallocate (cell, stat=istat)
if(istat/=0) stop 'Error deallocating cell'
deallocate (cell_ortho, stat=istat)
if(istat/=0) stop 'Error deallocating cell_ortho'

! Subroutines and functions
contains

! Function to print a string for an atom's element corresponding to an index
character(len=3) function element_id(num)
real(kind=dp) :: num

if (num == 1.0_dp) then
    element_id='Sr1' 
elseif (num == 2.0_dp) then
    element_id='Fe1'
elseif (num == 3.0_dp) then
    element_id='Fe2'
elseif (num == 4.0_dp) then
    element_id='Fe3'
elseif (num == 5.0_dp) then 
    element_id='Fe4'
elseif (num == 6.0_dp) then
    element_id='Fe5'
elseif (num == 7.0_dp) then
    element_id='O1 '
elseif (num == 8.0_dp) then
    element_id='O2 '
elseif (num == 9.0_dp) then
    element_id='O3 '
elseif (num == 10.0_dp) then
    element_id='O4 '
elseif (num == 11.0_dp) then
    element_id='O5 '
else
    element_id='N/A'
end if

return
end function element_id

! Function to print a string for an atom's material corresponding to an index
character(len=2) function material_id(num)
real(kind=dp) :: num

if (num == 1.0_dp) then
    material_id='Sr'
elseif (num == 2.0_dp) then
    material_id='Fe'
elseif (num == 3.0_dp) then
    material_id='O '
else
    material_id='NA'
end if

return
end function material_id

! Function to print a string for an atom's wyckoff position corresponding to an index
character(len=3) function wyckoff_id(num)
real(kind=dp) :: num

if (num == 1.0_dp) then
    wyckoff_id='2d '
elseif (num == 2.0_dp) then
    wyckoff_id='2a '
elseif (num == 3.0_dp) then
    wyckoff_id='2b '
elseif (num == 4.0_dp) then
    wyckoff_id='4f1'
elseif (num == 5.0_dp) then
    wyckoff_id='4f2'
elseif (num == 6.0_dp) then
    wyckoff_id='12k'
elseif (num == 7.0_dp) then
    wyckoff_id='4e '
elseif (num == 8.0_dp) then
    wyckoff_id='4f '
elseif (num == 9.0_dp) then
    wyckoff_id='6h '
else
    wyckoff_id='N/A'
end if

return
end function wyckoff_id

end program SrFeO_unit_check
