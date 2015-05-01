! This Fortran program attempts to find the nearest neighbor distances (nnd) of sites in proteins.

! GOAL:
!       -- To investigate the distribution of nearest neighbor distances of sites in proteins, using different atomic representations of sites in proteins.
! INPUT:
!       -- Path to the input file pdb_prop_CO.out containing pdb names and number of residues
!       -- Path to the input representative atomic coordinates and Bfactors (../../../properties/res_crd_bf)
!       -- Path and name of the output file containing the residue name and number and distance of the nearest neighbor of sites.

! Amir Shahmoradi, Thursday 6:00 PM, April 31 2015, iCMB, UT Austin

program get_res_prop_nnd
implicit none
! pdb file variables:
  integer, parameter                            :: npdb=213                ! number of pdbs in pdb_prop_CO.out file
  integer                                       :: natoms, nres  
  character(len=3), dimension(:)  , allocatable :: res_name
  integer         , dimension(:)  , allocatable :: res_num
  real*8          , dimension(:,:), allocatable :: crd
  real*8          , dimension(:)  , allocatable :: bfactor
  character(len=6)                              :: pdb,pdb_name            ! The name of the pdb structure and the SINGLE chain contained in it.
! input files and variables:
  character(len=300) :: co_in                                              ! input pdb_prop_CO.out file, containing the number of amino acids in each pdb file.
  character(len=300) :: crd_in                                             ! input pdb file. ATTN: Any non ATOM records will be ignored!
  character(len=300) :: nnd_out	                                           ! output file containing on each line, the nnd of each site.
  character(len=300) :: sum_out	                                           ! output summary file containing the pdb names and the best performing cutoff_distances and the corresponding Spearman correlations and the length of the protein on the last column.
! input files variables:
  integer            :: iostat, counter
  integer            :: co_in_unit=10,crd_in_unit=11                        ! input file units
  integer            :: nnd_out_unit=21
  integer            :: ios
! variables and parameters for nnd calculations:
  real*8 , dimension(:)     , allocatable :: nnd                    ! The array of nnd for all sites in a given protein.
  integer, dimension(:)     , allocatable :: nnid                   ! The array of residue IDs of the nearest neighbours of each individual site in a given protein.
! other variables:
  integer                                 :: i,ii,j
  real*8 :: tstart,tend

! First get the command line arguments
  if (command_argument_count()/=3) then
    write(*,*)
    write(*,*) "Incorrect number of input arguments on the command line."
    write(*,*) "Correct use:"
    write(*,'(1A81)') "./a.out <input: pdb_prop_CO.out> <input: crd_bf file> <output: res_prop_nnd file>"
    write(*,*)
    stop
  end if

  call get_command_argument(1,co_in)
  call get_command_argument(2,crd_in)
  call get_command_argument(3,nnd_out)
  
! OPEN INPUT FILES:
  open(unit=co_in_unit,file=trim(adjustl(co_in)),status='old')
  open(unit=crd_in_unit,file=trim(adjustl(crd_in)),status='old')
  read(co_in_unit,*)        ! read file header
  read(crd_in_unit,*)       ! read file header
! OPEN OUTPUT FILE:
  open(unit=nnd_out_unit,file=trim(adjustl(nnd_out)),status='replace')
  write(nnd_out_unit,'(7A20)') 'pdb','aa1_resname','aa2_resname','aa1_resnum','aa2_resnum','seq_distance','nnd'
  
call cpu_time(tstart)

do ii = 1,npdb

  read(co_in_unit,*) pdb,natoms,nres
  write(*,*) "processing PDB number ", ii, " : ", pdb, " out of 213 pdbs"
  allocate(crd(3,nres),bfactor(nres),nnd(nres),nnid(nres),res_name(nres),res_num(nres))
  
  ! NOW READ THE CRD FILE:
  do j = 1,nres
    read(crd_in_unit,*) pdb_name, res_name(j), res_num(j), crd(1:3,j), bfactor(j)
    if (pdb /= pdb_name) then
      write(*,*) 'Something is heavily fishy here!'
      write(*,*) 'pdb names do not match!', pdb, pdb_name; write(*,*) 'Program Aborted'
      STOP
    end if
  end do

  
  ! NOW CALCULATE NND FOR ALL VALUES OF PARAMETERS:
  call nnd_finder(nres,crd,nnid,nnd)
  
  ! NOW WRITE OUT THE FIRST FOUR LINES OF THE OUTPUT nnd FILE
  do j = 1,nres
    write(nnd_out_unit,'(3A20,3I20,1F20.3)') pdb,res_name(j),res_name(nnid(j)),res_num(j),res_num(nnid(j)),abs(res_num(j)-res_num(nnid(j))),nnd(j)
  end do

  deallocate(crd,bfactor,nnd,nnid,res_name,res_num)
  
end do
call cpu_time(tend)

write(*,*) new_line('a'), 'Overall it took ', tend-tstart, ' seconds for all 213 proteins.' 

close(co_in_unit)
close(crd_in_unit)
close(nnd_out_unit)

end program get_res_prop_nnd


! This Fortran subroutine takes in the 3D coordinates of a set of atoms (or sites). On the output it gives the nearest neighbor distance (nnd) and the corresponding ID of the neighboring site for all the atoms that were given in the input.

subroutine nnd_finder(nres,crd,nnid,nnd)
  implicit none
  integer, intent(in)          :: nres
  real*8 , intent(in)          :: crd(3,nres)
  real*8 , intent(out)         :: nnd(nres)
  integer, intent(out)         :: nnid(nres)
  integer                      :: i,j
  real*8                       :: distance_sq(nres)
  nnd = 0.d0
  nnid = 0
  do i=1,nres
    distance_sq(i) = huge(distance_sq)
    do j=1,nres
      if (i /= j) then
        distance_sq(j) = (crd(1,i)-crd(1,j))**2 + (crd(2,i)-crd(2,j))**2 + (crd(3,i)-crd(3,j))**2
      end if
    end do
    nnd(i)    = sqrt( minval(distance_sq) )
    nnid(i:i) = minloc(distance_sq)
  end do
end subroutine nnd_finder