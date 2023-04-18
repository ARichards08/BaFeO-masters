program BaFeO_hcp
implicit none

integer, parameter :: dp=selected_real_kind(15, 300)
real(kind=dp), parameter :: pi=4.0_dp*atan(1.0_dp)
character(len=*), parameter :: interaction_fmt="(6I4.1, E15.6E2)" , test_fmt="(BN, G15.6E2, BN, G15.6E2)"

! prim_cell is the primitive unit cell. cell and cell_ortho are the regular unit cell we are using in
! either non-orthogonal or orthogonal bases. n_cell and n_cell_ortho are the same but contain 
! 27 copies of the basis cell to account for neighbour interactions. distances holds the distances

real(kind=dp), dimension(:, :), allocatable :: atom_data, prim_cell, n_cell_ortho, cell, cell_ortho, distances
real(kind=dp), dimension(:, :), allocatable :: Fe_distances, mag_wyckoff
real(kind=dp), dimension(25) :: correct_Fe_d
real(kind=dp) :: a, b, c, x, y, z, t, r, r2, atom_i, atom_j, min_r, min_2r, U
logical :: anisotropy_calculation


integer, dimension(25) :: correct_nn
integer :: N_prim, N, istat, i, j, k, l, counter, nn, x_prim_mul, y_prim_mul, x_unit, y_unit, prim_in_unit, atom_types
integer :: failcount, x_int, y_int, z_int, mag_atoms, cell_mag_atoms, mag_atom_types, nnn

! File output
character(:), allocatable :: fname, ucf_fname, mat_fname, interaction_type, num_char, mag_mo_char
character(len=132), dimension(:, :), allocatable :: mat_array
character(len=132) :: tnum, tnum2
integer :: num_materials, interactions, O_mat
integer, dimension(:, :), allocatable :: int_out
real(kind=dp), dimension(:, :), allocatable :: exchange_out

! Testing
real(kind=dp), dimension(4) :: angle_out

! Visualisation of the hcp unit cell http://lampx.tugraz.at/~hadley/ss1/crystalstructure/structures/hcp/hcp.php

! Define potential U in eV
U=9.25_dp

! Will this be an anisotropy calculation
anisotropy_calculation=.TRUE.

! ucf and mat file name
fname="SrFeO-test"

! Define lattice constants in orthogonal and non-orthogonal bases
! a is the characteristic length (nearest neighbour distance) and b and c are the multiples of a
! in their respective directions. Here a and c are for room temp SrFe12O19 from literature
!a=5.88121_dp
!a=0.58876_dp
a=5.895_dp

b=a ! non-ortho so same as a
!c=23.023_dp/a
!c=2.31885_dp/a
c=23.199_dp

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

atom_data(:, 1) = [1.0_dp/3.0_dp, 2.0_dp/3.0_dp, 3.0_dp/4.0_dp, 1.0_dp, 5.0_dp, 1.0_dp]
atom_data(:, 2) = [0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 2.0_dp]
atom_data(:, 3) = [0.0_dp, 0.0_dp, 1.0_dp/4.0_dp, 2.0_dp, 1.0_dp, 3.0_dp]
atom_data(:, 4) = [1.0_dp/3.0_dp, 2.0_dp/3.0_dp, 0.02728_dp, 2.0_dp, 2.0_dp, 4.0_dp]
atom_data(:, 5) = [1.0_dp/3.0_dp, 2.0_dp/3.0_dp, 0.19080_dp, 2.0_dp, 3.0_dp, 5.0_dp]
atom_data(:, 6) = [0.16878_dp, 2*(0.16878_dp), -0.10921_dp, 2.0_dp, 4.0_dp, 6.0_dp]
atom_data(:, 7) = [0.0_dp, 0.0_dp, 0.1516_dp, 3.0_dp, 6.0_dp, 7.0_dp]
atom_data(:, 8) = [1.0_dp/3.0_dp, 2.0_dp/3.0_dp, -0.0552_dp, 3.0_dp, 7.0_dp, 8.0_dp]
atom_data(:, 9) = [0.1828_dp, 2*(0.1828_dp), 1.0_dp/4.0_dp, 3.0_dp, 8.0_dp, 9.0_dp]
atom_data(:, 10) = [0.1563_dp, 2*(0.1563_dp), 0.0524_dp, 3.0_dp, 9.0_dp, 6.0_dp]
atom_data(:, 11) = [0.5051_dp, 2*(0.5051_dp), 0.1510_dp, 3.0_dp, 10.0_dp, 6.0_dp]

! From the data, find the number of different element_ids
! that have Fe atoms assigned to them, and
! save the corresponding different wyckoff positions to an array

! mag_wyckoff
! 1st column: central atom wyckoff, 2nd column: interacting atom wyckoff
! 3rd column: nn interaction distance, 4th column: nn number of interactions
! 5th column: nnn interaction distance, 6th column: nnn number of interactions
! because every magnetic wyckoff position can interact with all the others including itself in neighbouring cells,
! mag_wyckoff is 6 x (number of magnetic wyckoff positions)**2 size
mag_atom_types=0
do i=1, size(atom_data, 2), 1
    if (trim(element_id(atom_data(4, i))) == 'Fe') mag_atom_types=mag_atom_types+1
end do

allocate (mag_wyckoff(6, mag_atom_types**2), stat=istat)
if(istat/=0) stop 'Error allocating mag_wyckoff'

mag_wyckoff=0.0_dp

counter=1
do i=1, size(atom_data, 2), 1
    if (trim(element_id(atom_data(4, i))) == 'Fe') then
        do j=1, size(atom_data, 2), 1
            if (trim(element_id(atom_data(4, j))) == 'Fe') then
                mag_wyckoff(1:2, counter)=[atom_data(6, i), atom_data(6, j)]
                counter=counter+1 
            end if
        end do
   end if
end do

! Calculate the number of magnetic atoms in a primitive unit cell using the multiplicity of the magnetic wyckoff positions
mag_atoms=0

do i=1, size(atom_data, 2), 1
    if (trim(element_id(atom_data(4, i))) == 'Fe') then
        mag_atoms = mag_atoms + multiplicity(atom_data(6, i))
    end if
end do

cell_mag_atoms=mag_atoms*prim_in_unit

! Create primitive cell
! Contains all the non-orthogonal position data, element_id, material_id, wyckoff_id
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
! Same data as atomic 
allocate (cell(7, N), stat=istat)
if(istat/=0) stop 'Error allocating cell'

counter=1
do l=1, N_prim, 1
    do x_unit=0, x_prim_mul-1, 1
        do y_unit=0, y_prim_mul-1, 1
            cell(:, counter)=[prim_cell(1, l)+real(x_unit, kind=dp), prim_cell(2, l)+real(y_unit, kind=dp), prim_cell(3, l)&
            &, prim_cell(4, l), prim_cell(5, l), prim_cell(6, l), real(counter-1, kind=dp)]
            counter=counter+1
        end do
    end do
end do

! Create an array to put the unit cell into orthogonal basis
allocate (cell_ortho(7, N), stat=istat)
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

! Create the unit with the 26 surrouding identical neighbours, held in n_cell_ortho.
! The surrounding cells are created by copying the initial cell with the corner on the origin,
! and adding or subtracting x_prim_mul or y_prim_mul or -1 or 1 in the z direction
! for at least 1 of each of the axes
! n_cell ortho will hold this information in indices 8-10

allocate (n_cell_ortho(10, N*27), stat=istat)
if(istat/=0) stop 'Error allocating n_cell_ortho'

counter=1
do l=1, N, 1
    do i=-x_prim_mul, x_prim_mul, x_prim_mul
        do j=-y_prim_mul, y_prim_mul, y_prim_mul
            do k=-1, 1, 1
                n_cell_ortho(:, counter) =[cell_ortho(1, l)+real(i, kind=dp), cell_ortho(2, l)+real(j, kind=dp), cell_ortho(3, l)+&
                &real(k, kind=dp), cell_ortho(4, l), cell_ortho(5, l), cell_ortho(6, l), cell_ortho(7, l),&
                & real(i, kind=dp)/real(x_prim_mul, kind=dp), real(j, kind=dp)/real(y_prim_mul, kind=dp), real(k, kind=dp)]
                counter=counter+1
            end do
        end do
    end do
end do


! Now need to calculate distances between the atoms.

! Need to only compare the centre unit cell to itself and the surrouding cells
! This is done by rounding the atoms position's to 0.0001 to check that the two
! atoms aren't in the same place

! Distances array holds the number of interactions N*((N*27)-1), and the atom_ids of the central atom, and then
! the interacting atom
! Atom_i is the central unit cell atom and atom_j is the interacting atom
! Indices 4-5 hold the wyckoff positions for atom i and atom j
! Indices 6-8 hold the information of which surrounding cell the interacting atom is in
allocate (distances(8, (N)*((N*27)-1)), stat=istat)
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
        x_int=nint(x*100000)
        y_int=nint(y*100000)
        z_int=nint(z*100000)

        ! To prevent comparing an atom with itself
        ! Distances are calculated in fractional 
        if (x_int/=0 .or. y_int/=0 .or. z_int/=0) then
            r2=(x*a)**2 + (y*cos(pi/6.0_dp)*b)**2 + (z*c)**2
            r=sqrt(r2)

            distances(:, counter)=[r, cell_ortho(7, i), n_cell_ortho(7, j), cell_ortho(6, i), n_cell_ortho(6, j)&
            &, n_cell_ortho(8, j), n_cell_ortho(9, j), n_cell_ortho(10, j)]
            counter=counter+1
        else
            ! Used for bug testing
            failcount=failcount+1
        end if
    end do
end do

!print*, failcount



! Calculating the distances only between magnetic (Fe in this case) atoms

! We only care about the distances between the iron atoms for the exchange interaction constants

allocate (Fe_distances(8, (cell_mag_atoms)*((cell_mag_atoms*27)-1)), stat=istat)
if(istat/=0) stop 'Error allocating Fe_distances'

Fe_distances=0.0_dp
k=1

do i=1, size(distances, 2), 1
    atom_i=distances(2, i)
    atom_j=distances(3, i)

    counter=0
    ! Checks to see if the central atom and the interacting atom are both magnetic (Fe here)
    ! If they both are, the distances can be added to Fe_distances
    ! Note that the atom_i+1 is because of the atom_id being -1 to the array index
    if (trim(element_id(cell_ortho(4, nint(atom_i+1)))) == "Fe" .and. trim(element_id(cell_ortho(4, nint(atom_j+1)))) == "Fe") then
        counter=1
    end if

    ! If both are magnetic
    if (counter == 1) then
        Fe_distances(:, k) = distances(:, i)
        k=k+1
    end if
end do

! For each magnetic wyckoff position, calculate the two shortest distances and save 
! the distances and the number of times they occur in mag_wyckoff
do i=1, size(mag_wyckoff, 2), 1

    ! Finds the minimum interaction distance between each combination of two magnetic wyckoff positions, min_r
    counter=0
    do j=1, size(Fe_distances, 2), 1
        ! If the distance in Fe_distances is between the sites we are looking at in this loop 
        if ((Fe_distances(4, j) == mag_wyckoff(1, i)) .and. (Fe_distances(5, j) == mag_wyckoff(2, i))) then
            ! This if statement checks to see if this is the first time the loop has found a distance in the list belonging to i
            ! If so needs to set this distance as the temporary min_r value   
            if (counter == 0) then
                min_r=round6(Fe_distances(1, j))
                counter=counter+1
            else
                ! Finds min_r for that wyckoff position
                if (round6(Fe_distances(1, j)) < min_r) min_r=round6(Fe_distances(1, j))
            end if
        end if
    end do

    ! Saves the min_r value to the mag_wyckoff array
    mag_wyckoff(3, i)=min_r

    ! Finds the number of nn
    nn=0
    do j=1, size(Fe_distances, 2), 1
        if ((Fe_distances(4, j) == mag_wyckoff(1, i)) .and. (Fe_distances(5, j) == mag_wyckoff(2, i))) then
            if (round6(Fe_distances(1, j)) == min_r) nn=nn+1
        end if
    end do
    
    ! nn has a lot of double ups, everything is x2 because of the number of primitive cells in a unit cell
    ! and also the multiplicity of the interacting wyckoff positions

    nn = nn / (prim_in_unit * multiplicity(mag_wyckoff(2, i)))

    ! Saves nn to the mag_wyckoff array
    mag_wyckoff(4, i)=real(nn, kind=dp)


    ! Repeat for nnn
    counter=0
    do j=1, size(Fe_distances, 2), 1
        ! If the distance in Fe_distances is between the sites we are looking at in this loop
        if ((Fe_distances(4, j) == mag_wyckoff(1, i)) .and. (Fe_distances(5, j) == mag_wyckoff(2, i))) then
            ! This if statement checks to see if this is the first time the loop has found a distance in the list belonging to i
            ! If so needs to set this distance as the temporary min_2r value
            if (counter == 0 .and. round6(Fe_distances(1, j)) > min_r) then
                min_2r=round6(Fe_distances(1, j))
                counter=counter+1
            else
                ! Finds min_2r for that wyckoff position
                if (round6(Fe_distances(1, j)) > min_r .and. round6(Fe_distances(1, j)) < min_2r) min_2r=round6(Fe_distances(1, j))
            end if
        end if
    end do

    ! Saves the min_2r value to the mag_wyckoff array
    mag_wyckoff(5, i)=min_2r

    ! Finds the number of nnn
    nnn=0
    do j=1, size(Fe_distances, 2), 1
        if ((Fe_distances(4, j) == mag_wyckoff(1, i)) .and. (Fe_distances(5, j) == mag_wyckoff(2, i))) then
            if (round6(Fe_distances(1, j)) == min_2r) nnn=nnn+1
        end if
    end do

    ! nn has a lot of double ups, everything is x2 because of the number of primitive cells in a unit cell
    ! and also the multiplicity of the interacting wyckoff positions

    nnn = nnn / (prim_in_unit * multiplicity(mag_wyckoff(2, i)))

    ! Saves nn to the mag_wyckoff array
    mag_wyckoff(6, i)=real(nnn, kind=dp)

    
end do


! Correct values for the nearest neighbour number and the nearest neighbour distances for each type of interaction.
! Keepin mind that the distances are in nm here
correct_nn=[6, 2, 3, 3, 1, 2, 6, 3, 3, 1, 6, 6, 3, 1, 2, 6, 6, 1, 1, 2, 6, 6, 6, 6, 2]
correct_Fe_d=[0.589_dp, 0.580_dp, 0.346_dp, 0.557_dp, 0.305_dp, 0.580_dp, 0.589_dp, 0.619_dp, 0.367_dp, 0.371_dp, 0.346_dp,&
& 0.619_dp, 0.363_dp, 0.379_dp, 0.350_dp, 0.557_dp, 0.367_dp, 0.379_dp, 0.277_dp, 0.351_dp, 0.305_dp, 0.371_dp, 0.350_dp,&
& 0.351_dp, 0.291_dp]

! Printing results
interactions=0

!print *, 'Central atom wyckoff, Interacting atom wyckoff, Interaction distance, Number of interactions, Number from Lit,&
!& Relative distance error from Lit, Potential U (eV), Exchange Constant J, NNN distance, NNN number'

do i=1, size(mag_wyckoff, 2), 1
!    print *, wyckoff_id(mag_wyckoff(1, i)), ',', wyckoff_id(mag_wyckoff(2, i)), mag_wyckoff(3, i), nint(mag_wyckoff(4, i)),&
!    & correct_nn(i), (mag_wyckoff(3, i)-correct_Fe_d(i)*10)/(correct_Fe_d(i)*10), U,&
!    & exchange_J_Nov(mag_wyckoff(1, i), mag_wyckoff(2, i), U, 1), mag_wyckoff(5, i), nint(mag_wyckoff(6, i))

    ! Counting interactions of nn and some nnn
    ! Discounting some combinations, 2a+2b, 2a+2a, 2b+2b, 2a+4f2, 2b+4f1
    if ((trim(wyckoff_id(mag_wyckoff(1, i))) == '2a' .and. trim(wyckoff_id(mag_wyckoff(2, i))) == '2b') .or. &
    &(trim(wyckoff_id(mag_wyckoff(1, i))) == '2b' .and. trim(wyckoff_id(mag_wyckoff(2, i))) == '2a') .or. &
    &(trim(wyckoff_id(mag_wyckoff(1, i))) == '2a' .and. trim(wyckoff_id(mag_wyckoff(2, i))) == '4f2') .or. &
    &(trim(wyckoff_id(mag_wyckoff(1, i))) == '4f2' .and. trim(wyckoff_id(mag_wyckoff(2, i))) == '2a') .or. &
    &(trim(wyckoff_id(mag_wyckoff(1, i))) == '2b' .and. trim(wyckoff_id(mag_wyckoff(2, i))) == '4f1') .or. &
    &(trim(wyckoff_id(mag_wyckoff(1, i))) == '4f1' .and. trim(wyckoff_id(mag_wyckoff(2, i))) == '2b') .or. &
    &(trim(wyckoff_id(mag_wyckoff(1, i))) == '2a' .and. trim(wyckoff_id(mag_wyckoff(2, i))) == '2a') .or. &
    &(trim(wyckoff_id(mag_wyckoff(1, i))) == '2b' .and. trim(wyckoff_id(mag_wyckoff(2, i))) == '2b')) then

        continue
    else
        interactions=interactions + multiplicity(mag_wyckoff(2, i))*nint(mag_wyckoff(4, i))*prim_in_unit
    end if

    ! Only counting 4f1-12k and 12k-12k nnn
    if ((trim(wyckoff_id(mag_wyckoff(1, i))) == '4f1' .and. trim(wyckoff_id(mag_wyckoff(2, i))) == '12k') .or. &
    &(trim(wyckoff_id(mag_wyckoff(1, i))) == '12k' .and. trim(wyckoff_id(mag_wyckoff(2, i))) == '4f1') .or. &
    &(trim(wyckoff_id(mag_wyckoff(1, i))) == '12k' .and. trim(wyckoff_id(mag_wyckoff(2, i))) == '12k')) then
        interactions=interactions + multiplicity(mag_wyckoff(2, i))*nint(mag_wyckoff(6, i))*prim_in_unit
    end if 
    
end do
!print *, interactions

! Printing the wyckoff positions of each atom number
!do i=1, N, 1
!    if (trim(element_id(cell_ortho(4, i))) == "Fe") print *, nint(cell_ortho(7, i)), wyckoff_id(cell_ortho(6, i))
!end do


! Allocate arrays for output

! int_out holds atom sites i and j, as well as the x y and z values of which neighbour cell the interactions is happening from
! In addition to an integer expressing if the interaction is first (1) or second nearest neighbour (2)
allocate (int_out(6, interactions), stat=istat)
if(istat/=0) stop 'Error allocating int_out'

! Holds the value of the interaction constant, as well as the bond angle
allocate (exchange_out(2, interactions), stat=istat)
if(istat/=0) stop 'Error allocating exchange_out'

! Saving the output data to arrays to be sorted before being written to a file
counter=0
do i=1, size(mag_wyckoff, 2), 1

    do j=1, size(Fe_distances, 2), 1

        ! If the distance being checked is from the same wyckoff positions being checked
        if ((Fe_distances(4, j) == mag_wyckoff(1, i)) .and. (Fe_distances(5, j) == mag_wyckoff(2, i))) then
            ! Writing out nearest neighbour interactions
            if (round6(Fe_distances(1, j)) == mag_wyckoff(3, i)) then
                ! Check if interaction is 0 before printing
                if (exchange_J_Nov(mag_wyckoff(1, i), mag_wyckoff(2, i), U, 1) /= 0) then

                    int_out(:, counter+1)=[nint(Fe_distances(2:3, j)), nint(Fe_distances(6:8, j)), 1]
                    exchange_out(1, counter+1)=exchange_J_Nov(mag_wyckoff(1, i), mag_wyckoff(2, i), U, 1)
                    exchange_out(2, counter+1)=bond_angle(int_out(1, counter+1), int_out(2, counter+1), int_out(3, counter+1),&
                    & int_out(4, counter+1),int_out(5, counter+1))
                    counter=counter+1
                end if
            end if

            ! Writing out next nearest neighbour interactions for 4f1-12k and 12k-12k
            if ((round6(Fe_distances(1, j)) == mag_wyckoff(5, i)) .and. &
                &((trim(wyckoff_id(mag_wyckoff(1, i))) == '4f1' .and. trim(wyckoff_id(mag_wyckoff(2, i))) == '12k') .or. &
                &(trim(wyckoff_id(mag_wyckoff(1, i))) == '12k' .and. trim(wyckoff_id(mag_wyckoff(2, i))) == '4f1') .or. &
                &(trim(wyckoff_id(mag_wyckoff(1, i))) == '12k' .and. trim(wyckoff_id(mag_wyckoff(2, i))) == '12k'))) then

                ! Check if interaction is 0 before printing
                if (exchange_J_Nov(mag_wyckoff(1, i), mag_wyckoff(2, i), U, 1) /= 0) then

                    int_out(:, counter+1)=[nint(Fe_distances(2:3, j)), nint(Fe_distances(6:8, j)), 2]
                    exchange_out(1, counter+1)=exchange_J_Nov(mag_wyckoff(1, i), mag_wyckoff(2, i), U, 1)
                    exchange_out(2, counter+1)=bond_angle(int_out(1, counter+1), int_out(2, counter+1), int_out(3, counter+1),&
                    & int_out(4, counter+1),int_out(5, counter+1))
                    counter=counter+1
                end if
            end if
        end if

    end do
end do


print *, maxval(exchange_out(2, :))
print *, minval(exchange_out(2, :))
print *, sum(exchange_out, 2)/size(exchange_out, 2)
print *, size(exchange_out, 2)

!!!!!
! Testing
!!!!!
open (unit=30, file="angle-strength.dat", iostat=istat, status='replace')
if (istat/=0) stop "Error opening .dat file"

counter=0
do i=1, size(exchange_out, 2), 1

    angle_out=small_dist_test(int_out(1, i), int_out(2, i), int_out(3, i), int_out(4, i), int_out(5, i))

    if (angle_out(2) < 3.2 .and. angle_out(3) < 3.2) then
        write (unit=30, fmt=*, iostat=istat) trim(material_id(cell(5, (int_out(1, i)+1)))),&
    & " ", trim(material_id(cell(5, (int_out(2, i)+1)))), " ", trim(material_id(angle_out(1))), " ",&
    & exchange_out(2, i), angle_out(2:)
        if (istat/=0) stop "Error writing to .dat file 1"
        counter=counter+1
    end if
end do

close (unit=30, iostat=istat)
if (istat/=0) stop "Error closing .dat file"

print *, counter

! Creating the .ucf file to output

!Specifying constants for the model
num_materials=7
interaction_type=trim("isotropic")

! Writing the .ucf file name to a variable
ucf_fname=trim(fname//".ucf")

! Opening .ucf file
open (unit=10, file=ucf_fname, iostat=istat, status='replace')
if (istat/=0) stop "Error opening .ucf file" 

! Writing to .ucf file

! Header
write (unit=10, fmt=*, iostat=istat) "#unit cell size"
if (istat/=0) stop "Error writing to .ucf file header 1"
write (unit=10, fmt=*, iostat=istat) a*x_prim_mul, b*cos(pi/6.0_dp)*y_prim_mul, c
if (istat/=0) stop "Error writing to .ucf file header 2"
write (unit=10, fmt=*, iostat=istat) "#unit cell vectors"
if (istat/=0) stop "Error writing to .ucf file header 3"
write (unit=10, fmt=*, iostat=istat) "1,0,0"
if (istat/=0) stop "Error writing to .ucf file header 4"
write (unit=10, fmt=*, iostat=istat) "0,1,0"
if (istat/=0) stop "Error writing to .ucf file header 5"
write (unit=10, fmt=*, iostat=istat) "0,0,1"
if (istat/=0) stop "Error writing to .ucf file header 6"

! Atoms positions
write (unit=10, fmt=*, iostat=istat) "#Atoms"
if (istat/=0) stop "Error writing to .ucf file atoms 1"
write (unit=10, fmt=*, iostat=istat) N, num_materials
if (istat/=0) stop "Error writing to .ucf file atoms 2"

do i=1, N, 1
    ! Accounts for the issue of having 5 oxygen materials but only specifying one in the material file    
    if (nint(cell_ortho(5, i)) > 6) then
        O_mat=6
    else
        O_mat=nint(cell_ortho(5, i))
    end if

    write (unit=10, fmt=*, iostat=istat) nint(cell_ortho(7, i)), cell_ortho(1, i)/x_prim_mul, cell_ortho(2, i)/y_prim_mul,&
    & cell_ortho(3, i), O_mat, 0, 0 
    if (istat/=0) stop "Error writing to .ucf file atoms 3"
end do

! Interactions
write (unit=10, fmt=*, iostat=istat) "#Interactions"
if (istat/=0) stop "Error writing to .ucf file interactions 1"
write (unit=10, fmt=*, iostat=istat) interactions, interaction_type
if (istat/=0) stop "Error writing to .ucf file interactions 2"

! Writing interactions out, sorted by atom site i and then by atom site j
counter=0
do i=1, N, 1
    do j=1, N, 1
        do k=1, size(exchange_out, 2), 1
            if ((int_out(1, k) == i-1) .and. (int_out(2, k) == j-1)) then
                write (unit=10, fmt=interaction_fmt, iostat=istat) counter, int_out(:5, k), exchange_out(1, k)
                if (istat/=0) stop "Error writing to .ucf file interactions 3"

                counter=counter+1
            end if
        end do
    end do
end do

! Closing .ucf file
close (unit=10, iostat=istat)
if (istat/=0) stop "Error closing .ucf file"



! Creating the material file

! Writing the .mat file name to a variable
mat_fname=trim(fname//".mat")

! Opening .ucf file
open (unit=20, file=mat_fname, iostat=istat, status='replace')
if (istat/=0) stop "Error opening .mat file"

! Writing to .ucf file

! Header
write(tnum, *) num_materials
num_char=trim(adjustl(tnum))

write (unit=20, fmt=*, iostat=istat) "material:num-materials=", num_char
if (istat/=0) stop "Error writing to .mat file header"

! Material array
! Start with magnetic atoms
! First is names, second is wyckoff positions for the magnetic atoms

allocate (mat_array(num_materials, 2), stat=istat)
if(istat/=0) stop 'Error allocating mat_names_array'

mat_array(:5, 1)=["2a ", "2b ", "4f1", "4f2", "12k"]
mat_array(:, 2)=["a ", "b ", "f1", "g2", "k ", "Ba", "O "]

! Write magnetic atom information
do i=1, mag_atom_types, 1
    write(tnum, *) i
    num_char=trim(adjustl(tnum))

    write(tnum2, *) abs(mag_moment_Nov(trim(mat_array(i, 1)), U))
    mag_mo_char=trim(adjustl(tnum2))

    write (unit=20, fmt=*, iostat=istat) "material["//num_char//"]:material-name="//trim(mat_array(i, 2)) 
    if (istat/=0) stop "Error writing to .mat file 1"
    write (unit=20, fmt=*, iostat=istat) "material["//num_char//"]:material-element="//trim(mat_array(i, 2))
    if (istat/=0) stop "Error writing to .mat file 2"
    write (unit=20, fmt=*, iostat=istat) "material["//num_char//"]:damping-constant=1.0"
    if (istat/=0) stop "Error writing to .mat file 3"
    write (unit=20, fmt=*, iostat=istat) "material["//num_char//"]:atomic-spin-moment="//mag_mo_char//"!muB" 
    if (istat/=0) stop "Error writing to .mat file 4"
    write (unit=20, fmt=*, iostat=istat) "material["//num_char//"]:unit-cell-category="//num_char
    if (istat/=0) stop "Error writing to .mat file 5"

    if ((mat_array(i, 1) == "4f1") .or. (mat_array(i, 1) == "4f2")) then
        write (unit=20, fmt=*, iostat=istat) "material["//num_char//"]:initial-spin-direction=0, 0, -1"
        if (istat/=0) stop "Error writing to .mat file 6"
    else
        write (unit=20, fmt=*, iostat=istat) "material["//num_char//"]:initial-spin-direction=0, 0, +1"
        if (istat/=0) stop "Error writing to .mat file 7"

        ! Anisotropy constraints
        if (anisotropy_calculation .eqv. .TRUE.) then
            write (unit=20, fmt=*, iostat=istat) "material["//num_char//"]:constraint-angle-phi=45.0"
            if (istat/=0) stop "Error writing to .mat file 8"
        end if
    end if

!    write (unit=20, fmt=*, iostat=istat) "material["//num_char//"]:initial-spin-direction=random"
!    if (istat/=0) stop "Error writing to .mat file 8"
    write (unit=20, fmt=*, iostat=istat) "material["//num_char//"]:uniaxial-anisotropy-constant=1e-23"
    if (istat/=0) stop "Error writing to .mat file 9"
    write (unit=20, fmt=*, iostat=istat) ""
    if (istat/=0) stop "Error writing to .mat file 10"
end do

! Write non-magnetic atom information
do i=mag_atom_types+1, num_materials, 1
    write(tnum, *) i
    num_char=trim(adjustl(tnum))

    write (unit=20, fmt=*, iostat=istat) "material["//num_char//"]:material-name="//mat_array(i, 2)
    if (istat/=0) stop "Error writing to .mat file non mag 1"
    write (unit=20, fmt=*, iostat=istat) "material["//num_char//"]:material-element="//mat_array(i, 2)
    if (istat/=0) stop "Error writing to .mat file non mag 2"
    write (unit=20, fmt=*, iostat=istat) "material["//num_char//"]:damping-constant=1.0"
    if (istat/=0) stop "Error writing to .mat file non mag 3"
    write (unit=20, fmt=*, iostat=istat) "material["//num_char//"]:unit-cell-category="//num_char
    if (istat/=0) stop "Error writing to .mat file non mag 4"
    write (unit=20, fmt=*, iostat=istat) "material["//num_char//"]:non-magnetic=keep"
    if (istat/=0) stop "Error writing to .mat file non mag 5"
    write (unit=20, fmt=*, iostat=istat) ""
    if (istat/=0) stop "Error writing to .mat file non mag 6"
end do

! Closing .mat file
close (unit=20, iostat=istat)
if (istat/=0) stop "Error closing .mat file"


! Deallocate arrays
deallocate (atom_data, stat=istat)
if(istat/=0) stop 'Error deallocating atom_data'
deallocate (mag_wyckoff, stat=istat)
if(istat/=0) stop 'Error deallocating mag_wyckoff'
deallocate (prim_cell, stat=istat)
if(istat/=0) stop 'Error deallocating prim_cell'
deallocate (cell, stat=istat)
if(istat/=0) stop 'Error deallocating cell'
deallocate (cell_ortho, stat=istat)
if(istat/=0) stop 'Error deallocating cell_ortho'
deallocate (n_cell_ortho, stat=istat)
if(istat/=0) stop 'Error deallocating n_cell_ortho'
deallocate (distances, stat=istat)
if(istat/=0) stop 'Error deallocating distances'
deallocate (Fe_distances, stat=istat)
if(istat/=0) stop 'Error deallocating Fe_distances'
deallocate (mat_array, stat=istat)
if(istat/=0) stop 'Error deallocating mat_array'
deallocate (int_out, stat=istat)
if(istat/=0) stop 'Error deallocating int_out'
deallocate (exchange_out, stat=istat)
if(istat/=0) stop 'Error deallocating exchange_out'
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
    elseif ((trim(wyckoff_id(atom_wyckoff)) == '4f') .or. (trim(wyckoff_id(atom_wyckoff)) == '4f1') .or.&
        & (trim(wyckoff_id(atom_wyckoff)) == '4f2')) then
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
character(len=2) function element_id(num)
real(kind=dp) :: num

if (num == 1.0_dp) then
    element_id='Sr' 
elseif (num == 2.0_dp) then
    element_id='Fe'
elseif (num == 3.0_dp) then
    element_id='O '
else
    element_id='NA'
end if

return
end function element_id

! Function to print a string for an atom's material corresponding to an index
character(len=3) function material_id(num)
real(kind=dp) :: num

if (num == 0.0_dp) then
    material_id='Fe1'
elseif (num == 1.0_dp) then
    material_id='Fe2'
elseif (num == 2.0_dp) then
    material_id='Fe3'
elseif (num == 3.0_dp) then
    material_id='Fe4'
elseif (num == 4.0_dp) then
    material_id='Fe5'
elseif (num == 5.0_dp) then
    material_id='Sr '
elseif (num == 6.0_dp) then
    material_id='O1'
elseif (num == 7.0_dp) then
    material_id='O2'
elseif (num == 8.0_dp) then
    material_id='O3'
elseif (num == 9.0_dp) then
    material_id='O4'
elseif (num == 10.0_dp) then
    material_id='O5'
else
    material_id='N/A'
end if

return
end function material_id

! Function to return a string for an atom's wyckoff position corresponding to an index
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

! Function to return the multiplicity of an atom type when given it's wyckoff position id
integer function multiplicity(wyckoff)
real(kind=dp) :: wyckoff

if (wyckoff == 1.0_dp) then
    multiplicity=2
elseif (wyckoff == 2.0_dp) then
    multiplicity=2
elseif (wyckoff == 3.0_dp) then
    multiplicity=2
elseif (wyckoff == 4.0_dp) then
    multiplicity=4
elseif (wyckoff == 5.0_dp) then
    multiplicity=4
elseif (wyckoff == 6.0_dp) then
    multiplicity=12
elseif (wyckoff == 7.0_dp) then
    multiplicity=4
elseif (wyckoff == 8.0_dp) then
    multiplicity=4
elseif (wyckoff == 9.0_dp) then
    multiplicity=6
else
    multiplicity=0
end if

return
end function multiplicity

! Function to round a position or distance to 6dp
real(kind=dp) function round6(long_num)
real(kind=dp) :: long_num

round6=real((nint(long_num * 1000000)) * 0.000001, kind=dp)

return
end function round6


! Function to return the exchange integral value J(meV) for an iron atom pair and repulsion potential U(eV)
! Values come from Fig 2. and Fig 3. in https://journals.aps.org/prb/pdf/10.1103/PhysRevB.71.184433
! Then uses a quadratic fit to the data. Only recommended to interpolate between 3.47 and 10.41 eV
! Points are only actually at U= 3.47eV, 6.94eV and 10.41eV
real(kind=dp) function exchange_J_Nov(wyckoff_i, wyckoff_j, U, unit)
real(kind=dp), intent(in) :: wyckoff_i, wyckoff_j, U
integer, intent(in) :: unit

character(:), allocatable :: wyck_i, wyck_j
real(kind=dp) :: a, b, c

wyck_i=trim(wyckoff_id(wyckoff_i))
wyck_j=trim(wyckoff_id(wyckoff_j))

! f1=-, f2=-, k=+, a=+, b=+

if ((wyck_i=="2a" .and. wyck_j=="4f1") .or. (wyck_i=="4f1" .and. wyck_j=="2a")) then
    a=12.0847_dp ; b=-1.40185_dp ; c=0.0506747_dp
!    a=0.056709 ; b=-1.52217 ; c=12.831
!    exchange_J_Nov=-5.01_dp*10.0_dp**(-21.0_dp)
elseif ((wyck_i=="2a" .and. wyck_j=="12k") .or. (wyck_i=="12k" .and. wyck_j=="2a")) then
    a=2.35593_dp ; b=-0.507986_dp ; c=0.0253374_dp
!    a=0.0524558 ; b=-1.42884 ; c=12.3204
!    exchange_J_Nov=0.0_dp
!elseif ((wyck_i=="2b" .and. wyck_j=="4f1") .or. (wyck_i=="4f1" .and. wyck_j=="2b")) then
!    a=0.644068_dp ; b=-0.0903629_dp ; c=0.00351908_dp
!    a=0.0496203 ; b=-1.3895 ; c=12.0483
!    exchange_J_Nov=-1.00_dp*10.0_dp**(-22.0_dp)
elseif ((wyck_i=="2b" .and. wyck_j=="4f2") .or. (wyck_i=="4f2" .and. wyck_j=="2b")) then
    a=12.7797_dp ; b=-1.49954_dp ; c=0.0548976_dp
!    a=0.0311899 ; b=-0.913174 ; c=8.30876
!    exchange_J_Nov=-5.01_dp*10.0_dp**(-21.0_dp)
elseif ((wyck_i=="2b" .and. wyck_j=="12k") .or. (wyck_i=="12k" .and. wyck_j=="2b")) then
    a=4.69492_dp ; b=-0.627656_dp ; c=0.0260412_dp
!    a=0.0241971 ; b=-0.605029 ; c=4.64307
!    exchange_J_Nov=2.00_dp*10.0_dp**(-21.0_dp)
elseif ((wyck_i=="4f1" .and. wyck_j=="4f2") .or. (wyck_i=="4f2" .and. wyck_j=="4f1")) then
    a=3.83051_dp ; b=-0.608118_dp ; c=0.0274488_dp
!    a=0.025519 ; b=-0.579585 ; c=3.72495
!    exchange_J_Nov=1.00_dp*10.0_dp**(-21.0_dp)
elseif ((wyck_i=="4f1" .and. wyck_j=="12k") .or. (wyck_i=="12k" .and. wyck_j=="4f1")) then
    a=8.38983_dp ; b=-0.940263_dp ; c=0.0330793_dp
!    a=0.025519 ; b=-0.510957 ; c=2.36747
!    exchange_J_Nov=-4.01_dp*10.0_dp**(-21.0_dp)
elseif ((wyck_i=="4f2" .and. wyck_j=="12k") .or. (wyck_i=="12k" .and. wyck_j=="4f2")) then
    a=12.3559_dp ; b=-1.44092_dp ; c=0.05349_dp
!    a=0.00567078 ; b=-0.117902 ; c=0.731472
!    exchange_J_Nov=-5.01_dp*10.0_dp**(-21.0_dp)
elseif (wyck_i=="4f1" .and. wyck_j=="4f1") then
    a=5.8549_dp ; b=-0.927841_dp ; c=0.0413623_dp
!    a=0.0415904 ; b=-0.931174 ; c=5.86063
!    exchange_J_Nov=1.00_dp*10.0_dp**(-21.0_dp)
elseif (wyck_i=="4f2" .and. wyck_j=="4f2") then
    a=2.45098_dp ; b=-0.387636_dp ; c=0.0153073_dp
!    a=0.0148096 ; b=-0.380439 ; c=2.42676
!    exchange_J_Nov=6.01_dp*10.0_dp**(-22.0_dp)
elseif (wyck_i=="12k" .and. wyck_j=="12k") then
    a=1.81961_dp ; b=-0.38029_dp ; c=0.0210068_dp
!    a=0.0206911 ; b=-0.375048 ; c=1.79875
!    exchange_J_Nov=1.50_dp*10.0_dp**(-21.0_dp)
else
    a=0.0_dp ; b=0.0_dp ; c=0.0_dp
!    exchange_J_Nov=0.0_dp
end if

!exchange_J_Nov = a*U**2.0_dp + b*U + c
exchange_J_Nov = a + b*U + C*U**2.0_dp

! Incorporate the magnetic moments and a factor of 0.5

exchange_J_Nov = 0.5_dp*exchange_J_Nov*mag_moment_Nov(wyck_i, U)*mag_moment_Nov(wyck_j, U)

if (unit == 1) exchange_J_Nov = exchange_J_Nov * 1.6_dp*10.0_dp**(-22.0_dp)

return

end function exchange_J_Nov



! Function to return the exchange integral value J(meV) for an iron atom pair and repulsion potential U(eV)
! Values come from Fig 4. in https://link.springer.com/content/pdf/10.1038/s41598-021-81028-7.pdf
! Then uses a quadratic fit to the data. Only recommended to interpolate between 2 and 5 eV
! Points are only actually at U= 2eV, 3eV, 4eV and 5eV
real(kind=dp) function exchange_J_Tej(wyckoff_i, wyckoff_j, U, unit)
real(kind=dp), intent(in) :: wyckoff_i, wyckoff_j, U
integer, intent(in) :: unit

character(:), allocatable :: wyck_i, wyck_j
real(kind=dp) :: a, b, c

wyck_i=trim(wyckoff_id(wyckoff_i))
wyck_j=trim(wyckoff_id(wyckoff_j))

! f1=-, f2=-, k=+, a=+, b=+

if ((wyck_i=="2a" .and. wyck_j=="4f1") .or. (wyck_i=="4f1" .and. wyck_j=="2a")) then
    a=14.0361_dp ; b=-2.20987_dp ; c=0.125561_dp
elseif ((wyck_i=="2a" .and. wyck_j=="12k") .or. (wyck_i=="12k" .and. wyck_j=="2a")) then
    a=3.62791_dp ; b=-1.17164_dp ; c=0.113789_dp
elseif ((wyck_i=="2b" .and. wyck_j=="4f2") .or. (wyck_i=="4f2" .and. wyck_j=="2b")) then
    a=17.3831_dp ; b=-2.67679_dp ; c=0.137332_dp
elseif ((wyck_i=="2b" .and. wyck_j=="12k") .or. (wyck_i=="12k" .and. wyck_j=="2b")) then
    a=8.80572_dp ; b=-1.89204_dp ; c=0.137332_dp
elseif ((wyck_i=="4f1" .and. wyck_j=="4f2") .or. (wyck_i=="4f2" .and. wyck_j=="4f1")) then
    a=8.87164_dp ; b=-2.88397_dp ; c=0.357063_dp
elseif ((wyck_i=="4f1" .and. wyck_j=="12k") .or. (wyck_i=="12k" .and. wyck_j=="4f1")) then
    a=10.7511_dp ; b=-1.7343_dp ; c=0.102018_dp
elseif ((wyck_i=="4f2" .and. wyck_j=="12k") .or. (wyck_i=="12k" .and. wyck_j=="4f2")) then
    a=20.136_dp ; b=-4.07209_dp ; c=0.286435_dp
elseif (wyck_i=="4f1" .and. wyck_j=="4f1") then
    a=2.19339_dp ; b=-0.639574_dp ; c=0.0588565_dp
elseif (wyck_i=="12k" .and. wyck_j=="12k") then
    a=17.1202_dp ; b=-4.03049_dp ; c=0.282511_dp
else
    a=0.0_dp ; b=0.0_dp ; c=0.0_dp
end if

exchange_J_Tej = a + b*U + c*U**2.0_dp

! Incorporate the magnetic moments and a factor of 0.5

exchange_J_Tej = 0.5_dp*exchange_J_Tej*mag_moment_Tej(wyck_i, U)*mag_moment_Tej(wyck_j, U)

if (unit == 1) exchange_J_Tej = exchange_J_Tej * 1.6_dp*10.0_dp**(-22.0_dp)

return

end function exchange_J_Tej

! Function to return an iron atom's magnetic moment based on their wyckoff position and the Hubbard parameter U
real(kind=dp) function mag_moment_Nov(wyckoff, U)
real(kind=dp), intent(in) :: U
character(len=*), intent(in):: wyckoff
real(kind=dp) :: a, b, c, s

if (trim(wyckoff)=='2a') then
    a=3.6765_dp; b=0.108501_dp; c=-0.00436014_dp
    s=1.0_dp
elseif (trim(wyckoff)=='2b') then
    a=3.486_dp; b=0.12853_dp; c=-0.00498302_dp
    s=1.0_dp
elseif (trim(wyckoff)=='4f1') then
    a=3.388_dp; b=0.155043_dp; c=-0.00664402_dp
    s=-1.0_dp
elseif (trim(wyckoff)=='4f2') then
    a=3.33_dp; b=0.181556_dp; c=-0.00830503_dp
    s=-1.0_dp
elseif (trim(wyckoff)=='12k') then
    a=3.6855_dp; b=0.105331_dp; c=-0.00394489_dp
    s=1.0_dp
else
    a=0.0_dp; b=0.0_dp; c=0.0_dp
    s=1.0_dp
end if

mag_moment_Nov= s*(a + b*U + c*U**2.0_dp)

return

end function mag_moment_Nov


! Function to return an iron atom's magnetic moment based on their wyckoff position and the Hubbard parameter U
! Based off of Tej data, Only defined between 0eV and 5eV
real(kind=dp) function mag_moment_Tej(wyckoff, U)
real(kind=dp), intent(in) :: U
character(:), allocatable, intent(in):: wyckoff
real(kind=dp) :: a, b, c, s

if (trim(wyckoff)=='2a') then
    a=3.94552_dp; b=0.117854_dp; c=-0.00746542_dp
    s=1.0_dp
elseif (trim(wyckoff)=='2b') then
    a=3.68561_dp; b=0.150722_dp; c=-0.0105912_dp
    s=1.0_dp
elseif (trim(wyckoff)=='4f1') then
    a=3.73244_dp; b=0.126792_dp; c=-0.00666733_dp
    s=-1.0_dp
elseif (trim(wyckoff)=='4f2') then
    a=3.49567_dp; b=0.25712_dp; c=-0.0216979_dp
    s=-1.0_dp
elseif (trim(wyckoff)=='12k') then
    a=3.88846_dp; b=0.133988_dp; c=-0.0100592_dp
    s=1.0_dp
else
    a=0.0_dp; b=0.0_dp; c=0.0_dp
    s=1.0_dp
end if

mag_moment_Tej= s*(a + b*U + c*U**2.0_dp)

return

end function mag_moment_Tej

! Bond Angle Calculation function
! Atom_i and Atom_j are atom id numbers, and x y and z are the relative cell coordinates of j to i
real(kind=dp) function bond_angle(atom_i, atom_j, x, y, z)
integer, intent(in) :: atom_i, atom_j, x, y, z
integer :: O_atom_index, i
real(kind=dp) :: sum_dist, sum_dist_min, dist_x, dist_y, dist_z, O_x, O_y, O_z, i_x, i_y, i_z, j_x, j_y, j_z, i_j, i_O, j_O
character(:), allocatable :: mat_i, mat_j

! Array that will contain the acceptable oxygen wyckoff positions for a pair of iron wyckoff positions
character(len=3), dimension(2) :: O_accepted_wyck

! Set sum_dst_min to a high initial value
sum_dist_min=100000.0_dp

! Set iron atom material
mat_i=trim(material_id(cell_ortho(5, atom_i+1)))
mat_j=trim(material_id(cell_ortho(5, atom_j+1)))

! Check accepted oxygen wyckoff positions for the given iron materials
if ((mat_i=="Fe1" .and. mat_j=="Fe3") .or. (mat_i=="Fe3" .and. mat_j=="Fe1")) then
    do i=1, size(atom_data, 2), 1
        if (trim(material_id(atom_data(5, i))) == "O4") O_accepted_wyck=[wyckoff_id(atom_data(6, i)), "N/A"]
    end do
elseif ((mat_i=="Fe1" .and. mat_j=="Fe5") .or. (mat_i=="Fe5" .and. mat_j=="Fe1")) then
    do i=1, size(atom_data, 2), 1
        if (trim(material_id(atom_data(5, i))) == "O4") O_accepted_wyck=[wyckoff_id(atom_data(6, i)), "N/A"]
    end do
elseif ((mat_i=="Fe2" .and. mat_j=="Fe5") .or. (mat_i=="Fe5" .and. mat_j=="Fe2")) then
    do i=1, size(atom_data, 2), 1
        if (trim(material_id(atom_data(5, i))) == "O1") O_accepted_wyck=[wyckoff_id(atom_data(6, i)), "N/A"]
    end do
elseif ((mat_i=="Fe2" .and. mat_j=="Fe4") .or. (mat_i=="Fe4" .and. mat_j=="Fe2")) then
    do i=1, size(atom_data, 2), 1
        if (trim(material_id(atom_data(5, i))) == "O3") O_accepted_wyck=[wyckoff_id(atom_data(6, i)), "N/A"]
    end do
elseif ((mat_i=="Fe3" .and. mat_j=="Fe5") .or. (mat_i=="Fe5" .and. mat_j=="Fe3")) then
    do i=1, size(atom_data, 2), 1
        if (trim(material_id(atom_data(5, i))) == "O2") O_accepted_wyck=[wyckoff_id(atom_data(6, i)), "N/A"]
    end do
elseif (mat_i=="Fe4" .and. mat_j=="Fe4") then
    do i=1, size(atom_data, 2), 1
        if (trim(material_id(atom_data(5, i))) == "O3") O_accepted_wyck=[wyckoff_id(atom_data(6, i)), "N/A"]
    end do
elseif ((mat_i=="Fe4" .and. mat_j=="Fe5") .or. (mat_i=="Fe5" .and. mat_j=="Fe4")) then
    do i=1, size(atom_data, 2), 1
        if (trim(material_id(atom_data(5, i))) == "O5") O_accepted_wyck=[wyckoff_id(atom_data(6, i)), "N/A"]
    end do
elseif (mat_i=="Fe5" .and. mat_j=="Fe5") then
    do i=1, size(atom_data, 2), 1
        if (trim(material_id(atom_data(5, i))) == "O2") O_accepted_wyck(1)=wyckoff_id(atom_data(6, i))
        if (trim(material_id(atom_data(5, i))) == "O5") O_accepted_wyck(2)=wyckoff_id(atom_data(6, i))
    end do
else
    O_accepted_wyck=["N/A", "N/A"]
end if

!  Iterate through all atoms in n_cell_ortho, checking each oxygen atom
do i=1, size(n_cell_ortho, 2), 1

    ! Check if atom element is oxygen
    if (trim(element_id(n_cell_ortho(4, i))) == "O") then

        ! Check if the oxygen's wyckoff position should be accepted. Also if there are no accepted positions due to
        ! lack of data
        if ( any(O_accepted_wyck==wyckoff_id(n_cell_ortho(6, i))) .or. all(O_accepted_wyck=="N/A") ) then

            ! Atom_i will always be in the central unit cell. Atom_j is dictated by x, y and z
            ! Adding atom_i - O distance to sum_dist
            dist_x=(cell_ortho(1, atom_i+1) - (n_cell_ortho(1, i)+n_cell_ortho(8, i)*x_prim_mul))*a
            dist_y=(cell_ortho(2, atom_i+1) - (n_cell_ortho(2, i)+n_cell_ortho(9, i)*y_prim_mul))*b*cos(pi/6.0_dp)
            dist_z=(cell_ortho(3, atom_i+1) - (n_cell_ortho(3, i)+n_cell_ortho(10, i)))*c

            sum_dist=sqrt(dist_x**2.0_dp + dist_y**2.0_dp + dist_z**2.0_dp)

            ! Adding atom_j - O distance to sum_dist
            dist_x=((cell_ortho(1, atom_j+1)+x*x_prim_mul) - (n_cell_ortho(1, i)+n_cell_ortho(8, i)*x_prim_mul))*a
            dist_y=((cell_ortho(2, atom_j+1)+y*y_prim_mul) - (n_cell_ortho(2, i)+n_cell_ortho(9, i)*y_prim_mul))*b*cos(pi/6.0_dp)
            dist_z=((cell_ortho(3, atom_j+1)+z) - (n_cell_ortho(3, i)+n_cell_ortho(10, i)))*c

            sum_dist=sum_dist + sqrt(dist_x**2.0_dp + dist_y**2.0_dp + dist_z**2.0_dp)
    
            ! If the newly calculated sum_dist is greater than the sum_dist_min value, replace the min value and save O atom info
            if (sum_dist < sum_dist_min) then
                sum_dist_min=sum_dist
                O_atom_index=i
                O_x=n_cell_ortho(1, i)+n_cell_ortho(8, i)*x_prim_mul
                O_y=n_cell_ortho(2, i)+n_cell_ortho(9, i)*y_prim_mul
                O_z=n_cell_ortho(3, i)+n_cell_ortho(10, i)
            end if
        end if
    end if
end do

! Simplifying the equations, calculating the fractional cell coordinates including the relative neighbouring cell
i_x=cell_ortho(1, atom_i+1)
i_y=cell_ortho(2, atom_i+1)
i_z=cell_ortho(3, atom_i+1)

j_x=cell_ortho(1, atom_j+1)+x*x_prim_mul
j_y=cell_ortho(2, atom_j+1)+y*y_prim_mul
j_z=cell_ortho(3, atom_j+1)+z

! Calculate the side lengths, converting into angstrom
dist_x=(i_x - O_x)*a
dist_y=(i_y - O_y)*b*cos(pi/6.0_dp)
dist_z=(i_z - O_z)*c
i_O=sqrt(dist_x**2.0_dp + dist_y**2.0_dp + dist_z**2.0_dp)

dist_x=(j_x - O_x)*a
dist_y=(j_y - O_y)*b*cos(pi/6.0_dp)
dist_z=(j_z - O_z)*c
j_O=sqrt(dist_x**2.0_dp + dist_y**2.0_dp + dist_z**2.0_dp)

dist_x=(i_x - j_x)*a
dist_y=(i_y - j_y)*b*cos(pi/6.0_dp)
dist_z=(i_z - j_z)*c
i_j=sqrt(dist_x**2.0_dp + dist_y**2.0_dp + dist_z**2.0_dp)

! Cosine rule calculation for the bond angle in degrees
bond_angle=acos((i_O**2.0_dp+j_O**2.0_dp-i_j**2.0_dp)/(2.0_dp*i_O*j_O))*180.0_dp/pi

return

end function bond_angle


function small_dist_test(atom_i, atom_j, x, y, z)
real(kind=dp), dimension(4) :: small_dist_test
integer, intent(in) :: atom_i, atom_j, x, y, z
integer :: O_atom_index, i
real(kind=dp) :: sum_dist, sum_dist_min, dist_x, dist_y, dist_z, O_x, O_y, O_z, i_x, i_y, i_z, j_x, j_y, j_z, i_j, i_O, j_O
real(kind=dp) :: bond_angle, O_mat
character(:), allocatable :: mat_i, mat_j

! Array that will contain the acceptable oxygen wyckoff positions for a pair of iron wyckoff positions
character(len=3), dimension(2) :: O_accepted_wyck

! Set sum_dst_min to a high initial value
sum_dist_min=100000.0_dp

! Set iron atom material
mat_i=trim(material_id(cell_ortho(5, atom_i+1)))
mat_j=trim(material_id(cell_ortho(5, atom_j+1)))

! Check accepted oxygen wyckoff positions for the given iron materials
if ((mat_i=="Fe1" .and. mat_j=="Fe3") .or. (mat_i=="Fe3" .and. mat_j=="Fe1")) then
    do i=1, size(atom_data, 2), 1
        if (trim(material_id(atom_data(5, i))) == "O4") O_accepted_wyck=[wyckoff_id(atom_data(6, i)), "N/A"]
    end do
elseif ((mat_i=="Fe1" .and. mat_j=="Fe5") .or. (mat_i=="Fe5" .and. mat_j=="Fe1")) then
    do i=1, size(atom_data, 2), 1
        if (trim(material_id(atom_data(5, i))) == "O4") O_accepted_wyck=[wyckoff_id(atom_data(6, i)), "N/A"]
    end do
elseif ((mat_i=="Fe2" .and. mat_j=="Fe5") .or. (mat_i=="Fe5" .and. mat_j=="Fe2")) then
    do i=1, size(atom_data, 2), 1
        if (trim(material_id(atom_data(5, i))) == "O1") O_accepted_wyck=[wyckoff_id(atom_data(6, i)), "N/A"]
    end do
elseif ((mat_i=="Fe2" .and. mat_j=="Fe4") .or. (mat_i=="Fe4" .and. mat_j=="Fe2")) then
    do i=1, size(atom_data, 2), 1
        if (trim(material_id(atom_data(5, i))) == "O3") O_accepted_wyck=[wyckoff_id(atom_data(6, i)), "N/A"]
    end do
elseif ((mat_i=="Fe3" .and. mat_j=="Fe5") .or. (mat_i=="Fe5" .and. mat_j=="Fe3")) then
    do i=1, size(atom_data, 2), 1
        if (trim(material_id(atom_data(5, i))) == "O2") O_accepted_wyck=[wyckoff_id(atom_data(6, i)), "N/A"]
    end do
elseif (mat_i=="Fe4" .and. mat_j=="Fe4") then
    do i=1, size(atom_data, 2), 1
        if (trim(material_id(atom_data(5, i))) == "O3") O_accepted_wyck=[wyckoff_id(atom_data(6, i)), "N/A"]
    end do
elseif ((mat_i=="Fe4" .and. mat_j=="Fe5") .or. (mat_i=="Fe5" .and. mat_j=="Fe4")) then
    do i=1, size(atom_data, 2), 1
        if (trim(material_id(atom_data(5, i))) == "O5") O_accepted_wyck=[wyckoff_id(atom_data(6, i)), "N/A"]
    end do
elseif (mat_i=="Fe5" .and. mat_j=="Fe5") then
    do i=1, size(atom_data, 2), 1
        if (trim(material_id(atom_data(5, i))) == "O2") O_accepted_wyck(1)=wyckoff_id(atom_data(6, i))
        if (trim(material_id(atom_data(5, i))) == "O5") O_accepted_wyck(2)=wyckoff_id(atom_data(6, i))
    end do
else
    O_accepted_wyck=["N/A", "N/A"]
end if

!  Iterate through all atoms in n_cell_ortho, checking each oxygen atom
do i=1, size(n_cell_ortho, 2), 1

    ! Check if atom element is oxygen
    if (trim(element_id(n_cell_ortho(4, i))) == "O") then

        ! Check if the oxygen's wyckoff position should be accepted. Also if there are no accepted positions due to
        ! lack of data
        if ( any(O_accepted_wyck==wyckoff_id(n_cell_ortho(6, i))) .or. all(O_accepted_wyck=="N/A") ) then

            ! Atom_i will always be in the central unit cell. Atom_j is dictated by x, y and z
            ! Adding atom_i - O distance to sum_dist
            dist_x=(cell_ortho(1, atom_i+1) - (n_cell_ortho(1, i)+n_cell_ortho(8, i)*x_prim_mul))*a
            dist_y=(cell_ortho(2, atom_i+1) - (n_cell_ortho(2, i)+n_cell_ortho(9, i)*y_prim_mul))*b*cos(pi/6.0_dp)
            dist_z=(cell_ortho(3, atom_i+1) - (n_cell_ortho(3, i)+n_cell_ortho(10, i)))*c

            sum_dist=sqrt(dist_x**2.0_dp + dist_y**2.0_dp + dist_z**2.0_dp)

            ! Adding atom_j - O distance to sum_dist
            dist_x=((cell_ortho(1, atom_j+1)+x*x_prim_mul) - (n_cell_ortho(1, i)+n_cell_ortho(8, i)*x_prim_mul))*a
            dist_y=((cell_ortho(2, atom_j+1)+y*y_prim_mul) - (n_cell_ortho(2, i)+n_cell_ortho(9, i)*y_prim_mul))*b*cos(pi/6.0_dp)
            dist_z=((cell_ortho(3, atom_j+1)+z) - (n_cell_ortho(3, i)+n_cell_ortho(10, i)))*c

            sum_dist=sum_dist + sqrt(dist_x**2.0_dp + dist_y**2.0_dp + dist_z**2.0_dp)

            ! If the newly calculated sum_dist is greater than the sum_dist_min value, replace the min value and save O atom info
            if (sum_dist < sum_dist_min) then
                sum_dist_min=sum_dist
                O_atom_index=i
                O_x=n_cell_ortho(1, i)+n_cell_ortho(8, i)*x_prim_mul
                O_y=n_cell_ortho(2, i)+n_cell_ortho(9, i)*y_prim_mul
                O_z=n_cell_ortho(3, i)+n_cell_ortho(10, i)
            end if
        end if
    end if
end do

O_mat=n_cell_ortho(5, O_atom_index)

i_x=cell_ortho(1, atom_i+1)
i_y=cell_ortho(2, atom_i+1)
i_z=cell_ortho(3, atom_i+1)

j_x=cell_ortho(1, atom_j+1)+x*x_prim_mul
j_y=cell_ortho(2, atom_j+1)+y*y_prim_mul
j_z=cell_ortho(3, atom_j+1)+z

! Calculate the side lengths
! Need to take into account not just the fractional cell coordinates but also which relative
! neighbouring cell the atoms might be in
dist_x=(i_x - O_x)*a
dist_y=(i_y - O_y)*b*cos(pi/6.0_dp)
dist_z=(i_z - O_z)*c
i_O=sqrt(dist_x**2.0_dp + dist_y**2.0_dp + dist_z**2.0_dp)

dist_x=(j_x - O_x)*a
dist_y=(j_y - O_y)*b*cos(pi/6.0_dp)
dist_z=(j_z - O_z)*c
j_O=sqrt(dist_x**2.0_dp + dist_y**2.0_dp + dist_z**2.0_dp)

dist_x=(i_x - j_x)*a
dist_y=(i_y - j_y)*b*cos(pi/6.0_dp)
dist_z=(i_z - j_z)*c
i_j=sqrt(dist_x**2.0_dp + dist_y**2.0_dp + dist_z**2.0_dp)

bond_angle=acos((i_O**2.0_dp+j_O**2.0_dp-i_j**2.0_dp)/(2.0_dp*i_O*j_O))*180.0_dp/pi

small_dist_test=(/O_mat, i_O, j_O, i_j/)

return

end function small_dist_test

end program BaFeO_hcp

