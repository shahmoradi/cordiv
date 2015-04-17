! This Fortran program attempts to calculate the Contact Order of the 209 enzyme proteins in cordiv dataset format according to the original definition CO given by Plaxco et a. 1998, for a given range of contact cutoff distances.

! INPUT:
!       -- CN definition to be used, with weighting options: h (Heaviside step function) as the only option at the moment.
!       -- Path to the input file pdb_prop_CO.out containing pdb names and number of residues
!       -- Path to the input representative atomic coordinates and Bfactors (../../../properties/res_crd_bf)
!       -- Path and name of the output file containing the CO values for each protein in a row, for a range of cutoff values.

! Amir Shahmoradi, Friday 11:31 AM, April 17 2015, iCMB, UT Austin

program CO_cutoff_varied
  implicit none
! pdb file variables:
  integer, parameter                            :: npdb = 213                ! number of pdbs in pdb_prop_CO.out file
  integer                                       :: natoms, nres  
  character(len=3)                              :: res_name
  integer                                       :: res_num
  real*8          , dimension(:,:), allocatable :: crd
  real*8                                        :: bfactor
  character(len=6)                              :: pdb,pdb_name            ! The name of the pdb structure and the SINGLE chain contained in it.

! input files and variables:
  character(len=1)   :: model                                              ! the CN model and definition: h (Heaviside: cutoff)
  character(len=300) :: co_in                                              
  character(len=300) :: crd_in                                             
  character(len=300) :: exp_out	                                           ! output file containing on each line, the Spearman correlations of Seq.Ent of a single pdb file with different exponential definitions of WCN (different cutoff_distances).

! input files variables:
  integer            :: iostat, counter
  integer            :: co_in_unit=10,crd_in_unit=11                        ! input file units
  integer            :: exp_out_unit=21
  integer            :: ios
  character(len=100) :: exp_output_format

! variables and parameters for CO calculations:
  real*8                                  :: contact_order      ! function calculating CO
  real*8                                  :: free_param_min, free_param_max,stride
  integer                                 :: nstride     ! number of strides
  real*8 , dimension(:)     , allocatable :: free_param, CO   ! The array of the values of the free parameter and the corresponding Contact Order.
  character(len=6)                        :: dummy_char

! other variables:
  integer                                 :: i,ii,j
  real*8 :: tstart,tend

! First get the command line arguments
  if (command_argument_count()/=4) then
    write(*,*)
    write(*,*) "Incorrect number of input arguments on the command line."
    write(*,*) "Correct use:"
    write(*,'(1A105)') "./a.out <input: CN model: h> <input: pdb_prop_CO.out> <input: crd_bf file> <output: CO. vs. param. file>"
    write(*,*)
    stop
  end if
  call get_command_argument(1,model)
  call get_command_argument(2,co_in)
  call get_command_argument(3,crd_in)
  call get_command_argument(4,exp_out)
  if ( model == 'h' ) then
    free_param_min = 0.d0
    free_param_max = 50.d0
    stride = 0.1d0
  else
    write (*,*) 'invalid input model! : ', model
    write (*,*) 'Only simple cutoff model is supported at the moment! : ', model
    write (*,*) 'program aborted'
    STOP
  end if

! allocate free_param values:
  nstride = nint ( (free_param_max-free_param_min) / stride )
  allocate(free_param(nstride),CO(nstride))
  do i = 1,nstride
    free_param(i) = free_param_min + dble(i)*stride
  end do
  
! OPEN INPUT FILES:
  open(unit=co_in_unit,file=trim(adjustl(co_in)),status='old')
  open(unit=crd_in_unit,file=trim(adjustl(crd_in)),status='old')
  read(co_in_unit,*)        ! read file header
  read(crd_in_unit,*)       ! read file header
! OPEN OUTPUT FILES:
  write(dummy_char,'(1I6)') nstride
  exp_output_format = trim(adjustl('(1A20,' // trim(adjustl(dummy_char)) // 'F20.6)'))
  open(unit=exp_out_unit,file=trim(adjustl(exp_out)),status='replace')
  write(exp_out_unit,exp_output_format) 'free_parameter',(free_param(i),i=1,nstride)
  
call cpu_time(tstart)
do ii = 1,npdb

  read(co_in_unit,*) pdb,natoms,nres
  write(*,*) "processing PDB number ", ii, " : ", pdb, " out of 213 pdbs"
  allocate(crd(3,nres))

  ! Now read the crd file:
  do j = 1,nres
    read(crd_in_unit,*) pdb_name, res_name, res_num, crd(1:3,j), bfactor
  end do
  
  ! Now calculate CN for all values of parameters:
  do i = 1,nstride
    CO(i) = contact_order(model,free_param(i),nres,crd)
  end do

  deallocate(crd)
  
  ! NOW WRITE OUT THE FIRST FOUR LINES OF THE OUTPUT CO FILE
  write(exp_out_unit,exp_output_format) pdb,(CO(i),i=1,nstride)

end do
call cpu_time(tend)

write(*,*) new_line('a'), 'Overall it took ', tend-tstart, ' seconds for all 213 proteins.' 

close(co_in_unit)
close(crd_in_unit)
close(exp_out_unit)

end program CO_cutoff_varied


! This Fortran subroutine takes in the 3D coordinates of a set of atoms or amino acid representative coordinates. On the output it gives the contact_order of the protein given the free parameter of the model.

function contact_order(model,free_param,nres,crd)
  
  implicit none
  character(len=1), intent(in) :: model
  integer, intent(in)          :: nres
  real*8, intent(in)           :: free_param,crd(3,nres)
  real*8                       :: contact_order     ! function
  integer                      :: i,j
  real*8                       :: distance_sq,free_param_sq
  integer                      :: ncontacts   ! number of contacts detected in the file

  ncontacts = 0.d0
  contact_order = 0.d0
  free_param_sq = free_param*free_param
  do i=1,nres-1
    do j=i+1,nres
      distance_sq = (crd(1,i)-crd(1,j))**2 + (crd(2,i)-crd(2,j))**2 + (crd(3,i)-crd(3,j))**2
      !if (model == 'h') then
        if (distance_sq <= free_param_sq) then
          contact_order = contact_order + dble(j-i)
          ncontacts = ncontacts + 1
        end if
      !end if
    end do
  end do
  
  contact_order = contact_order/dble(nres*ncontacts)
  
end function contact_order
