! This Fortran program attempts to find the distance of individual residues from the geometrical center of the protein.
! This code now also outputs the radius of gyration of proteins ina separater file.
! This code uses the output data from get_pdb_prop_CO.f90 and get_res_crd_bf.py

! INPUT:
!       -- Path to the input file pdb_prop_CO.out containing pdb names and number of residues
!       -- Path to the input representative atomic coordinates and Bfactors (../properties/res_crd_bf)
!       -- Path and name of the output file containing distance from COM of each residue for all pdb structures.
!       -- Path and name of the output file containing the radii of gyrations of PDBs.

! Amir Shahmoradi, Saturday 1:27 PM, April 18 2015, iCMB, UT Austin

program get_res_prop_dist_from_COM
implicit none

! pdb file variables:
  integer, parameter                            :: npdb=213                ! number of pdbs in pdb_prop_CO.out file
  integer                                       :: natoms, nres  
  character(len=3), dimension(:)  , allocatable :: res_name
  integer         , dimension(:)  , allocatable :: res_num
  real*8          , dimension(:,:), allocatable :: crd
  real*8          , dimension(:)  , allocatable :: bfactor,bfactor_min,bfactor_max,distance
  real*8          , dimension(3)                :: COM                     ! geometrical center of protein
  real*8                                        :: radius_of_gyration      ! geometrical center of protein
  character(len=6)                              :: pdb,pdb_name            ! The name of the pdb structure and the SINGLE chain contained in it.
  
  
! input files and variables:
  character(len=300) :: co_in                                              ! input pdb file. ATTN: Any non ATOM records will be ignored!
  character(len=300) :: crd_in                                             ! input pdb file. ATTN: Any non ATOM records will be ignored!
  character(len=300) :: output                                             ! output file containing distances of residues from COM of proteins
  character(len=300) :: rogout                                             ! output file containing radii of gyrations
  
! input files variables:
  integer            :: iostat, counter
  integer            :: co_in_unit=10,crd_in_unit=11                       ! input file units
  integer            :: output_unit=21,rogout_unit=22                      ! output file units
  integer            :: ios

! other variables:
  integer                                 :: i,ii,j
  real*8 :: tstart,tend

! First get the command line arguments
  if (command_argument_count()/=4) then
    write(*,*)
    write(*,*) "Incorrect number of input arguments on the command line."
    write(*,*) "Correct use:"
    write(*,'(1A180)') "./a.out <input: pdb_prop_CO.out> <input: crd_bf file> <output file: distances from COM>"
    write(*,*)
    stop
  end if
  call get_command_argument(1,co_in)
  call get_command_argument(2,crd_in)
  call get_command_argument(3,output)
  call get_command_argument(4,rogout)
  
! OPEN INPUT FILES:
  open(unit=co_in_unit,file=trim(adjustl(co_in)),status='old')
  open(unit=crd_in_unit,file=trim(adjustl(crd_in)),status='old')
  read(co_in_unit,*)        ! read file header
  read(crd_in_unit,*)       ! read file header

  open(unit=output_unit,file=trim(adjustl(output)),status='replace')
  write(output_unit,'(8A25)') 'pdb','resnam','resnum','distance_from_COM','distance_normalized','bf','bfmin','bfmax'
  
  open(unit=rogout_unit,file=trim(adjustl(rogout)),status='replace')
  write(rogout_unit,'(3A25)') 'pdb','nres','radius_gyration'
  
call cpu_time(tstart)
do ii = 1,npdb

  read(co_in_unit,*) pdb,natoms,nres
  write(*,*) "processing PDB number ", ii, " : ", pdb, " out of 213 pdbs"
  allocate(crd(3,nres),bfactor(nres),bfactor_min(nres),bfactor_max(nres),res_name(nres),res_num(nres),distance(nres))
  ! Now read the crd file:
  do j = 1,nres
    read(crd_in_unit,*) pdb_name, res_name(j), res_num(j), crd(1:3,j), bfactor(j), bfactor_min(j), bfactor_max(j)
    if (pdb /= pdb_name) then
      write(*,*) 'Something is heavily fishy here!'
      write(*,*) 'pdb names do not match!', pdb, pdb_name; write(*,*) 'Program Aborted'
      STOP
    end if
  end do
  
  ! now calculate COM and distances for the input PDB:
  do i = 1,3
    com(i) = sum(crd(i,1:nres)) / dble(nres)
  end do
  radius_of_gyration = 0.d0
  do i = 1,nres
    distance(i) = sqrt( (crd(1,i)-com(1))**2 + (crd(2,i)-com(2))**2 + (crd(3,i)-com(3))**2 )
    radius_of_gyration = radius_of_gyration + distance(i)*distance(i)
  end do
   radius_of_gyration = sqrt( radius_of_gyration / dble(nres) )

  do i = 1,nres
    ! NOW WRITE OUT THE FIRST FOUR LINES OF THE OUTPUT WCN FILE
    write(output_unit,'(2A25,1I25,5E25.5)') pdb, res_name(i), res_num(i), distance(i), distance(i)/radius_of_gyration, bfactor(i), bfactor_min(i), bfactor_max(i)
  end do
  
  write(rogout_unit,'(1A25,1I25,1E25.5)') pdb, nres, radius_of_gyration
  
  deallocate(crd,bfactor,bfactor_min,bfactor_max,distance,res_name,res_num)
    
end do
call cpu_time(tend)

write(*,*) new_line('a'), 'Overall it took ', tend-tstart, ' seconds for all 213 proteins.' 

close(co_in_unit)
close(crd_in_unit)
close(output_unit)
close(rogout_unit)

end program get_res_prop_dist_from_COM
