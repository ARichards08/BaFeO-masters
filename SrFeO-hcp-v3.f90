program BaFeO_hcp
implicit none

integer, parameter :: dp=selected_real_kind(15, 300)
real(kind=dp), parameter :: pi=4.0_dp*atan(1.0_dp)

! prim_cell is the primitive unit cell. cell and cell_ortho are the regular unit cell we are using in
! either non-orthogonal or orthogonal bases. n_cell and n_cell_ortho are the same but contain 
! 27 copies of the basis cell to account for neighbour interactions. distances holds the distances

real(kind=dp), dimension(:, :), allocatable :: atom_data, prim_cell, n_cell, n_cell_ortho, cell, cell_ortho, distances
real(kind=dp) :: a, b, c, x, y, z, t, r, r2, min_r, min_r2, min_r3
integer :: N_prim, N, istat, i, j, k, l, counter, nn, nnn, nnnn, x_prim_mul, y_prim_mul, x_unit, y_unit, prim_in_unit, atom_types
integer :: failcount, x_int, y_int, z_int

! This is a cool hcp sim http://lampx.tugraz.at/~hadley/ss1/crystalstructure/structures/hcp/hcp.php

! Define lattice constants in orthogonal and non-orthogonal bases
! a is the characteristic length (nearest neighbour distance) and b and c are the multiples of a
! in their respective directions. Here a and c are for room temp SrFe12O19 from literature
a=5.88121_dp

b=1.0_dp
c=23.023_dp/a

! Define number of atoms in primitive unit cell n
! Found by summing the multiplicities of the atom's wyckoff position
N_prim=64

! Define number of primitive unit cells in the unit cell we will use
! Will probably be something like 4, 2 in x and y, to preserve translational symmetry when converting to an orthogonal basis
! the number of primitive cells used in the unit cell in either the x or the y directions are held in x_prim_mul and y_prim_mul
! This program is for materials/structures using the P63/mmc space group, where the z/c basis vector is the same in the 
! non-orthogonal and orthogonal bases, so there doesn't need to be multiple primitive cells in the z direction

x_prim_mul=1
y_prim_mul=2

prim_in_unit=x_prim_mul*y_prim_mul

! Atoms in the total unit cell
N=N_prim*prim_in_unit

! Atom data is based of data from table 1: https://www.sciencedirect.com/science/article/pii/030488539290253K
! Fe2 is 2b not 4e. The functions at the bottom of the program can be used as a key for the element, material and wyckoff ids

atom_types=11

allocate (atom_data(6, atom_types), stat=istat)
if(istat/=0) stop 'Error allocating atom_cell'

! a, b, c, element_id, material_id, wyckoff_id

atom_data(:, 1) = [1.0_dp/3.0_dp, 2.0_dp/3.0_dp, 3.0_dp/4.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
atom_data(:, 2) = [0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 2.0_dp, 2.0_dp]
atom_data(:, 3) = [0.0_dp, 0.0_dp, 1.0_dp/4.0_dp, 3.0_dp, 2.0_dp, 3.0_dp]
atom_data(:, 4) = [1.0_dp/3.0_dp, 2.0_dp/3.0_dp, 0.02728_dp, 4.0_dp, 2.0_dp, 4.0_dp]
atom_data(:, 5) = [1.0_dp/3.0_dp, 2.0_dp/3.0_dp, 0.19080_dp, 5.0_dp, 2.0_dp, 4.0_dp]
atom_data(:, 6) = [0.16878_dp, 2*(0.16878_dp), -0.10921_dp, 6.0_dp, 2.0_dp, 5.0_dp]
atom_data(:, 7) = [0.0_dp, 0.0_dp, 0.1516_dp, 7.0_dp, 3.0_dp, 6.0_dp]
atom_data(:, 8) = [1.0_dp/3.0_dp, 2.0_dp/3.0_dp, -0.0552_dp, 8.0_dp, 3.0_dp, 4.0_dp]
atom_data(:, 9) = [0.1828_dp, 2*(0.1828_dp), 1.0_dp/4.0_dp, 9.0_dp, 3.0_dp, 7.0_dp]
atom_data(:, 10) = [0.1563_dp, 2*(0.1563_dp), 0.0524_dp, 10.0_dp, 3.0_dp, 5.0_dp]
atom_data(:, 11) = [0.5051_dp, 2*(0.5051_dp), 0.1510_dp, 11.0_dp, 3.0_dp, 5.0_dp]


! Create primitive cell
allocate (prim_cell(6, N_prim), stat=istat)
if(istat/=0) stop 'Error allocating prim_cell'


call data_to_prim(atom_data, prim_cell)

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

! Create the unit with the 26 surrouding identical neighbours, held in n_cell.
! The surrouing cells are created by copying the initial cell with the corner on the origin,
! and adding or subtracting x_prim_mul or y_prim_mul or -1 or 1 in the z direction
! for at least 1 of each of the axes
allocate (n_cell(6, N*27), stat=istat)
if(istat/=0) stop 'Error allocating n_cell'

counter=1
do l=1, N, 1
    do i=-x_prim_mul, x_prim_mul, x_prim_mul
        do j=-y_prim_mul, y_prim_mul, y_prim_mul
            do k=-1, 1, 1 
                n_cell(:, counter) =[cell(1, l)+real(i, kind=dp), cell(2, l)+real(j, kind=dp), cell(3, l)+real(k, kind=dp)&
                &, cell(4, l), cell(5, l), cell(6, l)]
                counter=counter+1
            end do
        end do
    end do
end do

! Need to now change to an orthogonal basis
! Originally set the positions to the same fractions

! Only the x coordinates change, as the change in length in the y dorection doesn't affect the fractional constants.
! If in non ortho basis an atoms x+ (sin 30)*y > a, need to loop it in ortho basis

! Now do it for the neighbouring cells (27)
allocate (n_cell_ortho(6, N*27), stat=istat)
if(istat/=0) stop 'Error allocating n_cell_ortho'

n_cell_ortho=n_cell

do i=1, N*27, 1
    t=n_cell(1, i) + n_cell(2, i)*sin(pi*7.0_dp/6.0_dp)
    if (t > 2.0_dp*real(x_prim_mul, kind=dp)) then
        n_cell_ortho(1, i) = t-3.0_dp*real(x_prim_mul, kind=dp)
    elseif (t < -real(x_prim_mul, kind=dp)) then
        n_cell_ortho(1, i) = t+3.0_dp*real(x_prim_mul, kind=dp)
    else
        n_cell_ortho(1, i) = t
    end if
end do

! Now need to calculate distances between the atoms.

! Need to only compare the centre unit cell to itself and the surrouding cells
! This is done by rounding the atoms position's to 0.0001 to check that the two
! atoms aren't in the same place

! Distances array holds the number of interactions N*((N*27)-1), and the element, material and wyckoff ids
! of the atom in the center cell, as well as the atom it is interacting with
! for the second index, 1 is the distance, 2-4 is the center atom information, 5-7 is the interacting atom information
allocate (distances((N)*((N*27)-1), 7), stat=istat)
if(istat/=0) stop 'Error allocating distances'


! Loop over every atom twice to compare distances, skipping comparing the same atom to itself
! Also give each distance an atom_id, saying which position in the central unit cell it relates to
counter=1
failcount=0
do i=1, N, 1
    do j=1, N*27, 1
        x=cell_ortho(1, i)-n_cell_ortho(1, j)
        y=cell_ortho(2, i)-n_cell_ortho(2, j)
        z=cell_ortho(3, i)-n_cell_ortho(3, j)
        x_int=nint(x*10000)
        y_int=nint(y*10000)
        z_int=nint(z*10000)



        if (x_int/=0 .or. y_int/=0 .or. z_int/=0) then
            r2=(cell_ortho(1, i)-n_cell_ortho(1, j))**2 + ((cell_ortho(2, i)-n_cell_ortho(2, j))*cos(pi/6.0_dp)*b)**2&
            & + ((cell_ortho(3, i) - n_cell_ortho(3, j))*c)**2
            r=sqrt(r2)

            distances(counter, :)=[r, cell_ortho(4, i), cell_ortho(5, i), cell_ortho(6, i)&
            &, n_cell_ortho(4, i), n_cell_ortho(5, i), n_cell_ortho(6, i)]
            counter=counter+1
        else
            ! Used for bug testing
            failcount=failcount+1
        end if
    end do
end do

!print*, failcount

! Round distances to 4 dp
do i=1, size(distances, 1), 1
    distances(i, 1)=real(nint(distances(i, 1)*10000), kind=dp)*0.0001_dp
end do

! Find the smallest distance
min_r=distances(1, 1)
do i=1, size(distances, 1), 1
    if (distances(i, 1) < min_r) min_r=distances(i, 1)
end do

! Find number of nn 
nn=0
do i=1, size(distances, 1), 1
    if (distances(i, 1) == min_r .and. distances(i, 2)==2) nn=nn+1
end do

print *, 'Nearest:'
print *, min_r, min_r*a, nn

! Find nnn distance and number
min_r2=0.0_dp
do i=1, size(distances, 1), 1
    if (distances(i, 1) > min_r2) min_r2=distances(i, 1)
end do
do i=1, size(distances, 1), 1
    if (distances(i, 1) > min_r .and. distances(i, 1) < min_r2) min_r2=distances(i, 1)
end do

nnn=0
do i=1, size(distances, 1), 1
    if (distances(i, 1) == min_r2 .and. distances(i, 2)==2) nnn=nnn+1
end do

print *, 'Second nearest:'
print *, min_r2, min_r2*a, nnn

! Find nnnn distance and number
min_r3=0.0_dp
do i=1, size(distances, 1), 1
    if (distances(i, 1) > min_r3) min_r3=distances(i, 1)
end do
do i=1, size(distances, 1), 1
    if (distances(i, 1) > min_r2 .and. distances(i, 1) < min_r3) min_r3=distances(i, 1)
end do

nnnn=0
do i=1, size(distances, 1), 1
    if (distances(i, 1) == min_r3 .and. distances(i, 2)==2) nnnn=nnnn+1
end do

print *, 'Third nearest:'
print *, min_r3, min_r3*a, nnnn



! Deallocate arrays
deallocate (atom_data, stat=istat)
if(istat/=0) stop 'Error deallocating atom_data'
deallocate (prim_cell, stat=istat)
if(istat/=0) stop 'Error deallocating prim_cell'
deallocate (cell, stat=istat)
if(istat/=0) stop 'Error deallocating cell'
deallocate (cell_ortho, stat=istat)
if(istat/=0) stop 'Error deallocating cell_ortho'
deallocate (n_cell, stat=istat)
if(istat/=0) stop 'Error deallocating n_cell'
deallocate (n_cell_ortho, stat=istat)
if(istat/=0) stop 'Error deallocating n_cell_ortho'
deallocate (distances, stat=istat)
if(istat/=0) stop 'Error deallocating distances'

! Subroutines and functions
contains


! Subroutine to turn the atom data into a primitive unit cell using the wyckoff positions
subroutine data_to_prim(positions, cell)
real(kind=dp), dimension(:, :), intent(in) :: positions
real(kind=dp), dimension(:, :), intent(inout) :: cell
real(kind=dp) :: x, y, z, atom_element, atom_material, atom_wyckoff
integer :: i, counter

counter=1
do i=1, size(positions, 2)

    x=positions(1, i)
    y=positions(2, i)
    z=positions(3, i)
    atom_element=positions(4, i)
    atom_material=positions(5, i)
    atom_wyckoff=positions(6, i)

    if (trim(wyckoff_id(atom_wyckoff)) == '2d') then
        ! (1/3,2/3,3/4), (2/3,1/3,1/4)   
        cell(:, counter) = [1.0_dp/3.0_dp, 2.0_dp/3.0_dp, 3.0_dp/4.0_dp, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [2.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/4.0_dp, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
    elseif (trim(wyckoff_id(atom_wyckoff)) == '2a') then
        ! (0,0,0), (0,0,1/2)
        cell(:, counter) = [0.0_dp, 0.0_dp, 0.0_dp, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [0.0_dp, 0.0_dp, 1.0_dp/2.0_dp, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
    elseif (trim(wyckoff_id(atom_wyckoff)) == '2b') then
        ! (0,0,1/4), (0,0,3/4)
        cell(:, counter) = [0.0_dp, 0.0_dp, 1.0_dp/4.0_dp, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [0.0_dp, 0.0_dp, 3.0_dp/4.0_dp, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
    elseif (trim(wyckoff_id(atom_wyckoff)) == '4f') then
        ! (1/3,2/3,z), (2/3,1/3,z+1/2), (2/3,1/3,-z), (1/3,2/3,-z+1/2)
        cell(:, counter) = [1.0_dp/3.0_dp, 2.0_dp/3.0_dp, z, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [2.0_dp/3.0_dp, 1.0_dp/3.0_dp, z+1.0_dp/2.0_dp, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [2.0_dp/3.0_dp, 1.0_dp/3.0_dp, -z, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [1.0_dp/3.0_dp, 2.0_dp/3.0_dp, -z+1.0_dp/2.0_dp, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
    elseif (trim(wyckoff_id(atom_wyckoff)) == '12k') then
        ! (x,2x,z), (-2x,-x,z), (x,-x,z), (-x,-2x,z+1/2)
        ! (2x,x,z+1/2), (-x,x,z+1/2), (2x,x,-z), (-x,-2x,-z)
        ! (-x,x,-z), (-2x,-x,-z+1/2), (x,2x,-z+1/2), (x,-x,-z+1/2)
        cell(:, counter) = [x, 2.0_dp*x, z, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [-2.0_dp*x, -x, z, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [x, -x, z, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [-x, -2.0_dp*x, z+1.0_dp/2.0_dp, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [2.0_dp*x, x, z+1.0_dp/2.0_dp, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [-x, x, z+1.0_dp/2.0_dp, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [2.0_dp*x, x, -z, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [-x, -2.0_dp*x, -z, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [-x, x, -z, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [-2.0_dp*x, -x, -z+1.0_dp/2.0_dp, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [x, 2.0_dp*x, -z+1.0_dp/2.0_dp, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [x, -x, -z+1.0_dp/2.0_dp, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
    elseif (trim(wyckoff_id(atom_wyckoff)) == '4e') then
        ! (0,0,z), (0,0,z+1/2), (0,0,-z), (0,0,-z+1/2)
        cell(:, counter) = [0.0_dp, 0.0_dp, z, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [0.0_dp, 0.0_dp, z+1.0_dp/2.0_dp, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [0.0_dp, 0.0_dp, -z, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [0.0_dp, 0.0_dp, -z+1.0_dp/2.0_dp, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
    elseif (trim(wyckoff_id(atom_wyckoff)) == '6h') then
        ! (x,2x,1/4), (-2x,-x,1/4), (x,-x,1/4)
        ! (-x,-2x,3/4), (2x,x,3/4), (-x,x,3/4)
        cell(:, counter) = [x, 2.0_dp*x, 1.0_dp/4.0_dp, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [-2.0_dp*x, -x, 1.0_dp/4.0_dp, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [x, -x, 1.0_dp/4.0_dp, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [-x, -2.0_dp*x, 3.0_dp/4.0_dp, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [2.0_dp*x, x, 3.0_dp/4.0_dp, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
        cell(:, counter) = [-x, x, 3.0_dp/4.0_dp, atom_element, atom_material, atom_wyckoff]
        counter=counter+1
    end if
end do

end subroutine data_to_prim





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
    wyckoff_id='4f '
elseif (num == 5.0_dp) then
    wyckoff_id='12k'
elseif (num == 6.0_dp) then
    wyckoff_id='4e '
elseif (num == 7.0_dp) then
    wyckoff_id='6h '
else
    wyckoff_id='N/A'
end if

return
end function wyckoff_id

end program BaFeO_hcp
