program hcp
implicit none

integer, parameter :: dp=selected_real_kind(15, 300)
real(kind=dp), parameter :: pi=4.0_dp*atan(1.0_dp)

real(kind=dp), dimension(:, :), allocatable :: original, ortho, original_single, ortho_single, distances
!real(kind=dp), dimension(:), allocatable ::
real(kind=dp) :: a, b, c, x, y, z, t, r, r2, min_r, min_r2, min_r3
integer :: N, istat, i, j, k, counter, nn, nnn, nnnn, failcount

! This is a cool hcp sim http://lampx.tugraz.at/~hadley/ss1/crystalstructure/structures/hcp/hcp.php

! Define lattice constants in orthogonal and non-orthogonal bases
a=0.229_dp
b=a
c=sqrt(8.0_dp/3.0_dp)*a

!x=a
!y=a*cos(pi/6.0_dp)
!z=c

! First need to set up the rhombus unit cell, and then convert to the orthogonal unti cell
! a=x and b=y coordinate kind of. c=z.  a and b are equal in length, with a 60 deg angle between

! Define number of atoms in unit cell n
N=2*4

! Allocate arrays
! Assign atom positions in terms of lattice constant fractions
! This is in the non-orthogonal form

! Setup the unit cell, plus the 26 neighbouring cells. Do this by adding either -1, 0 or +1 to each coordinate
! In the positions of the unit cell, and looping over -1, 0, +1.
allocate (original(3, N*27), stat=istat)
if(istat/=0) stop 'Error allocating periodic_orig'
allocate (ortho(3, N*27), stat=istat)
if(istat/=0) stop 'Error allocating periodic_ortho'

counter=0
do i=-2, 2, 2
    do j=-2, 2, 2
        do k=-1, 1, 1
            original(:, 1+counter*N) = [0.0_dp+i, 0.0_dp+j, 0.0_dp+k]
            original(:, 2+counter*N) = [(1.0_dp/3.0_dp)+i, (1.0_dp/3.0_dp)+j, 0.5_dp+k]
            original(:, 3+counter*N) = [1.0_dp+i, 0.0_dp+j, 0.0_dp+k]
            original(:, 4+counter*N) = [(4.0_dp/3.0_dp)+i, (1.0_dp/3.0_dp)+j, 0.5_dp+k]
            original(:, 5+counter*N) = [0.0_dp+i, 1.0_dp+j, 0.0_dp+k]
            original(:, 6+counter*N) = [(1.0_dp/3.0_dp)+i, (4.0_dp/3.0_dp)+j, 0.5_dp+k]
            original(:, 7+counter*N) = [1.0_dp+i, 1.0_dp+j, 0.0_dp+k]
            original(:, 8+counter*N) = [(4.0_dp/3.0_dp)+i, (4.0_dp/3.0_dp)+j, 0.5_dp+k]
            counter=counter+1
        end do
    end do
end do

! Need to now change to an orthogonal basis
! Originally set the positions to the same fractions

! Only the x coordinates change, as the change in length in the y dorection doesn't affect the fractional constants.
! If in non ortho basis an atoms x+ (sin 30)*y > a, need to loop it in ortho basis

! Now do it for the neighbouring cells (27)

ortho=original

do i=1, N*27, 1
    t=original(1, i) + original(2, i)*sin(pi/6.0_dp)
    if (t > 4.0_dp) then
        ortho(1, i) = t-6.0_dp
    elseif (t < -2.0_dp) then
        ortho(1, i) = t+6.0_dp
    else
        ortho(1, i) = t
    end if
end do

! Now need to calculate distances between the atoms.

! FIrst create a sinle unit cell at the origin

allocate (original_single(3, N), stat=istat)
if(istat/=0) stop 'Error allocating periodic_orig'
allocate (ortho_single(3, N), stat=istat)
if(istat/=0) stop 'Error allocating periodic_ortho'

original_single(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
original_single(:, 2) = [(1.0_dp/3.0_dp), (1.0_dp/3.0_dp), 0.5_dp]
original_single(:, 3) = [1.0_dp, 0.0_dp, 0.0_dp]
original_single(:, 4) = [(4.0_dp/3.0_dp), (1.0_dp/3.0_dp), 0.5_dp]
original_single(:, 5) = [0.0_dp, 1.0_dp, 0.0_dp]
original_single(:, 6) = [(1.0_dp/3.0_dp), (4.0_dp/3.0_dp), 0.5_dp]
original_single(:, 7) = [1.0_dp, 1.0_dp, 0.0_dp]
original_single(:, 8) = [(4.0_dp/3.0_dp), (4.0_dp/3.0_dp), 0.5_dp]
ortho_single=original_single

do i=1, N, 1
    t=original_single(1, i) + original_single(2, i)*sin(pi/6.0_dp)
    if (t > 1.0_dp) then
        ortho_single(1, i) = t-1.0_dp
    else
        ortho_single(1, i) = t
    end if
end do

! Need to only compare the centre unit cell to itself and the surrouding cells
! Allocate distances array

allocate (distances((N)*((N*27)-1), 2), stat=istat)
if(istat/=0) stop 'Error allocating distances'


! Loop over every atom twice to compare distances, skipping comparing the same atom to itself
! Also give each distance an atom_id, saying which position in the central unit cell it relates to
counter=1
failcount=0
do i=1, N, 1
    do j=1, N*27, 1
        x=ortho_single(1, i)-ortho(1, j)
        y=ortho_single(2, i)-ortho(2, j)
        z=ortho_single(3, i)-ortho(3, j)
        x=nint(x*10000)
        y=nint(y*10000)
        z=nint(z*10000)

        if (x/=0 .or. y/=0 .or. z/=0) then
            r2=(ortho_single(1, i)-ortho(1, j))**2 + ((ortho_single(2, i)-ortho(2, j))*cos(pi/6.0_dp))**2 + ((ortho_single(3, i)&
            & - ortho(3, j))*sqrt(8.0_dp/3.0_dp))**2
            r=sqrt(r2)

            distances(counter, :)=[r, real(i, kind=dp)]
            counter=counter+1
        else
            failcount=failcount+1
        end if
    end do
end do
print *, 'fail=', failcount, 'pass=', counter

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

print *, 'Nearest'
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
deallocate (original, stat=istat)
if(istat/=0) stop 'Error deallocating original'
deallocate (ortho, stat=istat)
if(istat/=0) stop 'Error deallocating ortho'
deallocate (original_single, stat=istat)
if(istat/=0) stop 'Error deallocating original_single'
deallocate (ortho_single, stat=istat)
if(istat/=0) stop 'Error deallocating ortho_single'
deallocate (distances, stat=istat)
if(istat/=0) stop 'Error deallocating distances'

end program hcp
