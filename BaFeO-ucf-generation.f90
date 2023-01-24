program BaFeO_ucf
implicit none

integer :: stat, num_atoms, num_materials, num_interactions
character(len=40) :: filename, interaction_type

!Specifying constants for the model
num_atoms=128
num_materials=7
num_interactions=736
interaction_type=trim("isotropic")

! Writing the .ucf file name to a variable
filename="BaFeO-test.ucf"
filename=trim(filename)

! Opening .ucf file
open (unit=10, file=filename, iostat=stat, status='replace')
if (stat/=0) stop "Error opening .ucf file" 

!!!
! Writing to .ucf file
!!!

! Header
write (unit=10, fmt=*, iostat=stat) "#unit cell size"
if (stat/=0) stop "Error writing to .ucf file header 1"
write (unit=10, fmt=*, iostat=stat) 6.01681046,  10.421420066049103,  23.49536771
if (stat/=0) stop "Error writing to .ucf file header 2"
write (unit=10, fmt=*, iostat=stat) "#unit cell vectors"
if (stat/=0) stop "Error writing to .ucf file header 3"
write (unit=10, fmt=*, iostat=stat) "1,0,0"
if (stat/=0) stop "Error writing to .ucf file header 4"
write (unit=10, fmt=*, iostat=stat) "0,1,0"
if (stat/=0) stop "Error writing to .ucf file header 5"
write (unit=10, fmt=*, iostat=stat) "0,0,1"
if (stat/=0) stop "Error writing to .ucf file header 6"

! Atoms positions
write (unit=10, fmt=*, iostat=stat) "#Atoms"
if (stat/=0) stop "Error writing to .ucf file atoms 1"
write (unit=10, fmt=*, iostat=stat) num_atoms, num_materials
if (stat/=0) stop "Error writing to .ucf file atoms 2"

! Interactions
write (unit=10, fmt=*, iostat=stat) "#Interactions"
if (stat/=0) stop "Error writing to .ucf file interactions 1"
write (unit=10, fmt=*, iostat=stat) num_interactions, interaction_type
if (stat/=0) stop "Error writing to .ucf file interactions 2"

! Closing .ucf file
close (unit=10, iostat=stat)
if (stat/=0) stop "Error closing .ucf file"

end program BaFeO_ucf
