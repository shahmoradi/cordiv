! This Fortran program attempts to calculate the Contact Order of a given protein in PDB format according to the original definition CO given by Plaxco et a. 1998.

! GOALS:
!       -- To see if the CO of proteins has any predictive power for the observed diversity of the correlations between sequence entropy and structural variables such as WCN and RSA.
!       -- To find which definition of the Contact Order performs the best. The original definition, defines a contact between two residues if any atoms of the two amino acids are closer than 6 Angstroms to each other. An alternative definition could be just taking the center of mass of each of the side chains or the complete residues (with backbone atoms) and define contact between two residues if their centers of mass is closer than 6 Angstroms to each other.

! INPUT:
!       -- Path to the input pdb file
!       -- Path and name of the output file containing different definitions of Contact Order for each pdb structure on each line.

! METHOD:
!       The code calculates the distance matrix of all atoms in the PDB file and based on it, finds the Contact Order of the given protein.

! Amir Shahmoradi, Saturday 5:22 PM, August 23 2014, iCMB, UT Austin

program get_pdb_prob_CO
implicit none
! pdb file variables:
  character(len=6), dimension(:)  , allocatable :: record
  integer         , dimension(:)  , allocatable :: atom_num
  character(len=4), dimension(:)  , allocatable :: atom_name
  character(len=1), dimension(:)  , allocatable :: alt_loc_ind
  character(len=3), dimension(:)  , allocatable :: res_name
  character(len=1), dimension(:)  , allocatable :: chain_id
  integer         , dimension(:)  , allocatable :: res_num
  character(len=1), dimension(:)  , allocatable :: res_code
  real*8          , dimension(:,:), allocatable :: crd
  real*8          , dimension(:)  , allocatable :: occupancy
  real*8          , dimension(:)  , allocatable :: bfactor
  integer         , dimension(:)  , allocatable :: res_num_renumbered    ! This array of size natoms will contain the reordered residue number of each atom.
  integer                                       :: nres=0,natoms=0       ! The number of residues and atoms in the pdb file
  integer         , dimension(:)  , allocatable :: natomsSC,natomsAA     ! Arrays that will contain the number of atoms in each Amino Acid and each residue Side Chain.
  integer                                       :: res_num_old           ! A dummy variable that will contain the number of atoms in e
  character(len=3)                              :: res_nam_old           ! A dummy variable used to recognize duplicate Amino Acid entries in the pdb file (The AAs that two different conformations resolved for).
  character(len=1)                              :: res_cod_old           ! A dummy variable used to recognize duplicate Amino Acid entries in the pdb file (The AAs that two different conformations resolved for).
  character(len=6)                              :: pdb_name              ! The name of the pdb structure and the SINGLE chain contained in it.
  logical                                       :: side_chain_atom
  logical                                       :: not_same_residue
  logical         , dimension(:,:), allocatable :: no_residue_contact
  logical         , dimension(:)  , allocatable :: no_side_chain
! input files:
  character(len=300) :: pdb_in                                           ! input pdb file. ATTN: Any non ATOM records will be ignored!
  !character(len=300) :: elj_in                                          ! input file containing Eleisha's Seq.Ent for all pdb files. ATTN: The file should be in its original format written by Eleisha, that is, the columns ordering should have not been changes.
  character(len=300) :: output	                                         ! output file containing on each line, the Spearman correlations of Seq.Ent of a single pdb file with different Gaussian definitions of WCN (different STDEVs).
! input files variables:
  logical            :: output_exists, sum_out_exists                    ! if the output files already exist, these variables are true.
  integer            :: pdb_in_unit=11                                   ! input file units
  integer            :: output_unit=21                                   ! output file units
  integer            :: ios
! contact order variables:
  integer            :: ncontacts,ncontactsSC,ncontactsAA
  real*8             :: distance
  real*8, parameter  :: cmin = 0.d0; cmax = 30.0d0
  real*8             :: cutoff = 6.d0
  real*8             :: contact_order,contact_orderSC,contact_orderAA
! ELJ input file variables:
!  character(len=4)   :: elj_pdb_name
!  character(len=1)   :: elj_chain_id
!  integer            :: elj_res_num
!  real*8             :: elj_seqent, elj_ddgent  
!  real*8, dimension(:) , allocatable :: pdb_seqent, pdb_ddgent  
!  integer, dimension(:), allocatable :: pdb_resnum
! variables and parameters for wcn calculations and correlations:
  real*8 , parameter                      :: stdev_min = 0.0d0                                ! stdev stands for STandard DEViation.
  real*8 , parameter                      :: stdev_max = 50.d0                                ! stdev stands for STandard DEViation.
  real*8 , parameter                      :: stdev_stride = 0.2                               ! stdev stands for STandard DEViation.
  real*8 , parameter                      :: alpha_default = -2.d0                            ! The default value of alpha as in the original definition of wcn.
  integer, parameter                      :: nstdev = nint((stdev_max-stdev_min)/stdev_stride)
  real*8 , dimension(nstdev)              :: stdev,sp_cor,abs_sp_cor                          ! The array of contact numbers for each each CA atom. sp_cor stands for the spearman correlation of wcn and ELJ sequence entropies. abs stands for the absolute values of the spearman correlations.
  real*8 , dimension(:)     , allocatable :: wcn                                              ! The array of contact numbers for each each CA atom.
  real*8 , dimension(:,:)   , allocatable :: crdSC,crdAA                                      ! these vectors will contain the coordinates of the Center Of Mass of the residue Side Chains (SC) and the Amino Acid (AA).
  real*8 , dimension(:,:)   , allocatable :: crdCA                                            ! If an amino acid lacks a side chain, then the value of crdSC will be set to crd of CA atom of the same amino acid.
  real*8                                  :: spear                                            ! function
  real*8                                  :: sp_default                                       ! The value of Spearman correlation for the default case of alpha = -2 in WCN definition.
  character(len=6)                        :: dummy_char
  !character(len=100):: wcn_output_format,res_output_format,resnum_output_format
! other variables:
  integer            :: i,ires,j
! First get the command line arguments

  if (command_argument_count()/=2) then
    write(*,*)
    write(*,*) "Incorrect number of input arguments on the command line."
    write(*,*) "Correct use:"
    write(*,'(1A140)') "./a.out <input pdb file> <output Contact Order summary file>"
    write(*,*)
    stop
  end if
  call get_command_argument(1,pdb_in)
  call get_command_argument(2,output)
  !write (*,*) pdb_in
  !write (*,*) output
  !read(*,*)
  
! OPEN INPUT FILES:
  open(unit=pdb_in_unit,file=trim(adjustl(pdb_in)),status='old')
    
! FIRST DETERMINE THE FOUR-LETTER NAME OF THE PDB STRUCTURE:
  pdb_in = adjustr(pdb_in)
  i = len(adjustr(pdb_in))
  if (pdb_in(i-3:i) /= '.pdb') then ! .or. pdb_in(i-10:i-10) /= '_') then
    write(*,*) 'WARNING: pdb filename does not seem to be correct'
    write(*,*) 'pdb filename: ', pdb_in(i-9:i)
    write(*,*) 'The last 10 letters of pdb filename must be of the form: ***_*.pdb'
    STOP  ! read(*,*)
  else
    pdb_name = pdb_in(i-9:i-4)
    write(*,*); write(*,*) 'pdb name: ', pdb_name
  end if

! THEN DETERMINE THEW NUMBER OF ATOMS IN THE PDB FILE (ATTN:  Only ATOM records will be considered as pdb atoms)

  allocate (record(1))
  do
    read(pdb_in_unit,'(1A4)') record(1)
    if (trim(adjustl(record(1)))=='ATOM') exit
    !write(*,'(1A4)') record
    cycle
  end do
  
  deallocate (record)
  backspace(unit=pdb_in_unit)
  
  allocate(record(1),atom_num(1),atom_name(1),alt_loc_ind(1),res_name(1),chain_id(1),res_num(1),res_code(1),crd(1,3),occupancy(1),bfactor(1))
  
  do
    read(pdb_in_unit,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2)',iostat=ios) record(1),atom_num(1),atom_name(1),alt_loc_ind(1),res_name(1),chain_id(1),res_num(1),res_code(1),(crd(1,j),j=1,3),occupancy(1),bfactor(1)
    if (ios<0) then
      exit
    elseif (ios>0) then
      write(*,*) 'Sum Tin Wong! : PDB file broken.'; stop
    elseif (trim(adjustl(record(1)))=='TER') then
      write(*,*) 'Found TER record in pdb file: ', trim(adjustl(pdb_in))
    elseif (trim(adjustl(record(1)))=='END') then
      write(*,'(A3)') 'END'
      exit
    else
      if (trim(adjustl(record(1)))=='ATOM') natoms = natoms + 1
      if (trim(adjustl(record(1)))=='ATOM' .and. trim(adjustl(atom_name(1))) == 'CA') then
        if (nres == 0) then
          nres = nres + 1
        elseif (res_num_old == res_num(1) .and. res_nam_old == res_name(1) .and. res_cod_old == res_code(1)) then
            write(*,*) 'FATAL: Two CA atoms in the pdb file have the same res_num, res_name and res_code!', res_num(1), res_name(1), res_code(1)
        elseif (res_num_old == res_num(1) .and. res_nam_old == res_name(1) .and. res_cod_old /= res_code(1)) then
            write(*,*) 'WARNING: Duplicate Amino Acid found in structure!', res_num(1), res_name(1)
            res_num_old  = res_num(1)
            res_nam_old  = res_name(1)
            res_cod_old  = res_code(1)
            nres = nres + 1
            cycle
        else
          !write(*,*) 'Duplicate Amino Acid found in structure!', res_num(1), res_name(1)
          nres = nres + 1
        end if
        res_num_old = res_num(1)
        res_nam_old = res_name(1)
        res_cod_old  = res_code(1)
      end if
      cycle
    end if
  end do
  
  deallocate(record,atom_num,atom_name,alt_loc_ind,res_name,chain_id,res_num,res_code,crd,occupancy,bfactor)
  close(pdb_in_unit)
  
  write(*,'(1A9,1I8)') "natoms: ",natoms; write(*,'(1A9,1I8)') "nres: ",nres !; write(*,*)

! NOW RECORD THE CA COORDINATES OF THE PDB FILE:

  open(unit=pdb_in_unit,file=trim(adjustl(pdb_in)),status='old')
  allocate(crdSC(nres,3),crdAA(nres,3),crdCA(nres,3),natomsSC(nres),natomsAA(nres))
  allocate(record(natoms),atom_num(natoms),atom_name(natoms),alt_loc_ind(natoms),res_name(natoms),chain_id(natoms),res_num(natoms),res_code(natoms),crd(natoms,3),occupancy(natoms),bfactor(natoms),res_num_renumbered(natoms))
  allocate(no_side_chain(nres))
  do
    read(pdb_in_unit,'(1A4)') record(1)
    if (trim(adjustl(record(1)))=='ATOM') exit
    !write(*,'(1A4)') record
    cycle
  end do
  backspace(unit=pdb_in_unit)
  
  i = 1; ires = 1; natomsSC = 0; natomsAA = 0; crdSC = 0.d0; crdAA = 0.d0; res_num_renumbered(1) = 1; no_side_chain = .TRUE.
  do
    read(pdb_in_unit,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2)',iostat=ios) record(i),atom_num(i),atom_name(i),alt_loc_ind(i),res_name(i),chain_id(i),res_num(i),res_code(i),(crd(i,j),j=1,3),occupancy(i),bfactor(i)
    !write(*,*) 'Amir',i
    
    if (ios<0 .and. i<natoms) then
      write(*,'(1A100)') 'Sum Tin Wong! : reached the end of output PDB file, while expecting further records from the file.'; stop
    elseif (ios>0) then
      write(*,*) 'Sum Tin Wong! : input PDB file broken: ', pdb_name; stop  

    elseif (trim(adjustl(record(i)))=='ATOM') then
    
      side_chain_atom = trim(adjustl(atom_name(i))) /= 'N' .and. trim(adjustl(atom_name(i))) /= 'CA' .and. &
                        trim(adjustl(atom_name(i))) /= 'C' .and. trim(adjustl(atom_name(i))) /= 'O'  .and. &
                        trim(adjustl(atom_name(i))) /= 'OXT'
    
      if (i==natoms .and. ires/=nres) then

        write(*,*) 'i, ires: ', i,ires
        write(*,*) 'FATAL: Sum Tin Wong! input PDB file broken: ', pdb_name
        stop
      
      elseif (i == 1) then
      
        res_num_old  = res_num(i)
        res_nam_old  = res_name(i)
        res_cod_old  = res_code(i)
        res_num_renumbered(i) = ires
        !ires = ires + 1
        natomsAA(ires)  = natomsAA(ires) + 1
        crdAA(ires,1:3) = crdAA(ires,1:3) + crd(i,1:3)
        if (side_chain_atom) then
          no_side_chain(ires) = .FALSE.
          natomsSC(ires)  = natomsSC(ires) + 1
          crdSC(ires,1:3) = crdSC(ires,1:3) + crd(i,1:3)
        end if
        ! safety check:  First atom must always be 'N'
          if (trim(adjustl(atom_name(i))) /= 'N') then
            write(*,*) 'FATAL: The first atom of the PDB file is not "N" : ',trim(adjustl(atom_name(i))), pdb_name
            stop
          end if

      elseif (i > 1) then
        
        if (res_num_old == res_num(i) .and. res_nam_old == res_name(i) .and. res_cod_old == res_code(i)) then  ! The atom belongs to the same residue as the previous atom's.
          !write(*,*) 'atom name:', atom_name(i)
          res_num_renumbered(i) = ires
          natomsAA(ires)  = natomsAA(ires) + 1
          crdAA(ires,1:3) = crdAA(ires,1:3) + crd(i,1:3)
          if (trim(adjustl(atom_name(i))) == 'CA') crdCA(ires,1:3) = crd(i,1:3)
          if (side_chain_atom) then
            no_side_chain(ires) = .FALSE.
            natomsSC(ires)  = natomsSC(ires) + 1
            crdSC(ires,1:3) = crdSC(ires,1:3) + crd(i,1:3)
          end if
          if (i == natoms) then   ! Normalize crdSC and crdAA of the last residue to get center of mass:
            if (no_side_chain(ires)) then
              crdSC(ires,1:3) = crdCA(ires,1:3)
            else
              crdSC(ires,1:3) = crdSC(ires,1:3)/dble(natomsSC(ires))
            end if
            crdAA(ires,1:3) = crdAA(ires,1:3)/dble(natomsAA(ires))
            !write(*,*) res_name(i),res_code(i),natomsSC(ires),natomsAA(ires)
            !write(*,*) res_num_renumbered(i), ires
            !write(*,*) crdSC(ires,1:3)
            !write(*,*) crdAA(ires,1:3)
            !write(*,*) 'done',natoms
            exit
          end if
        else        !  It's an atom belonging to a new residue.
          ! First normalize crdSC and crdAA of the previous residue to get center of mass:
            if (no_side_chain(ires)) then
              crdSC(ires,1:3) = crdCA(ires,1:3)
            else
              crdSC(ires,1:3) = crdSC(ires,1:3)/dble(natomsSC(ires))
            end if
            crdAA(ires,1:3) = crdAA(ires,1:3)/dble(natomsAA(ires))
          res_num_old  = res_num(i)
          res_nam_old  = res_name(i)
          res_cod_old  = res_code(i)
          !write(*,*) res_name(i-1),res_code(i-1),natomsSC(ires),natomsAA(ires)
          !write(*,*) res_num_renumbered(i-1), ires
          !write(*,*) crdSC(ires,1:3)
          !write(*,*) crdAA(ires,1:3)
          !read(*,*)
          ires = ires + 1
          res_num_renumbered(i) = ires
          ! safety check: ires must be always <= nres
            if (nres < ires) then
              write(*,*) 'FATAL: the number of residues is more than expected in PDB file. Something is terribly wrong.', pdb_name
              stop
            end if
          natomsAA(ires)  = natomsAA(ires) + 1
          crdAA(ires,1:3) = crdAA(ires,1:3) + crd(i,1:3)
          if (trim(adjustl(atom_name(i))) == 'CA') crdCA(ires,1:3) = crd(i,1:3)
          if (side_chain_atom) then
            no_side_chain(ires) = .FALSE.
            natomsSC(ires)  = natomsSC(ires) + 1
            crdSC(ires,1:3) = crdSC(ires,1:3) + crd(i,1:3)
          end if
        end if
      end if
      !write(*,*) 'ordered_res_number: ', res_num_renumbered(i), ires
      i = i + 1
      cycle
    elseif (trim(adjustl(record(i)))=='TER') then
      cycle
    elseif (trim(adjustl(record(i)))=='END') then
      if (ires<=nres) then
        write(*,'(1A100)') 'Sum Tin Wong! : reached the end of output PDB file, while expecting further reading from the file.'
        write(*,*) 'number of residues read: ', i
        stop
      end if
      exit
    end if
  end do; close(pdb_in_unit)

! Now calculate the pairwise atomic distances, also pairwise AA and SC distances
  
  ! safety check point: 
    if (nres /= res_num_renumbered(natoms)) then
      write(*,*) 'disaster detected! Code has some serious flaws.'
      write(*,*) 'Total number of residues in the pdb file does not match the expected number.'
      stop
    end if
  allocate (no_residue_contact(nres,nres))
  no_residue_contact = .TRUE.     ! Initially assume no contact between any two atoms of a pair of residues.
  ncontacts          = 0
  do i = 1,natoms-1 
    do j = i+1,natoms
      not_same_residue = res_num(i) /= res_num(j) .or. res_name(i) /= res_name(j) .or. res_code(i) /= res_code(j)
      ! safety check point
        if ( not_same_residue .neqv. (res_num_renumbered(i) /= res_num_renumbered(j)) ) then
          write(*,*) 'disaster detected! Code has some serious flaws.'
          stop
        end if
      if ( not_same_residue .and. no_residue_contact(res_num_renumbered(i),res_num_renumbered(j)) ) then
        !write(*,*) 'Calculating CA contact order for residue number ', res_num(i)
        distance = sqrt( (crd(i,1)-crd(j,1))**2 + (crd(i,2)-crd(j,2))**2 + (crd(i,3)-crd(j,3))**2 )
        if (distance <= cutoff) then
          no_residue_contact(res_num_renumbered(i),res_num_renumbered(j)) = .FALSE.
          ncontacts = ncontacts + 1
          contact_order = contact_order + dble(abs(res_num_renumbered(j)-res_num_renumbered(i)))
          ! safety check point:
            if ( abs(res_num_renumbered(j)-res_num_renumbered(i)) /= res_num_renumbered(j)-res_num_renumbered(i) ) then
              write(*,*) 'Something is terribly wrong with res_num_renumbered'
              write(*,*) abs(res_num_renumbered(j)-res_num_renumbered(i)), res_num_renumbered(j)-res_num_renumbered(i)
              stop
            end if
        end if
      end if
    end do
  end do
  
  contact_order = contact_order/dble(nres*ncontacts)
  
  ncontactsSC   = 0
  ncontactsAA   = 0
  do i = 1,nres
    !write(*,*) 'Calculating AA/SC contact order for residue number ', i
    do j = i+1,nres
      distance = sqrt( (crdSC(i,1)-crdSC(j,1))**2 + (crdSC(i,2)-crdSC(j,2))**2 + (crdSC(i,3)-crdSC(j,3))**2 )
      if (distance <= 2.d0*cutoff) then
        ncontactsSC = ncontactsSC + 1
        contact_orderSC = contact_orderSC + dble(j-i)
      end if
      distance = sqrt( (crdAA(i,1)-crdAA(j,1))**2 + (crdAA(i,2)-crdAA(j,2))**2 + (crdAA(i,3)-crdAA(j,3))**2 )
      if (distance <= 2.d0*cutoff) then
        ncontactsAA = ncontactsAA + 1
        contact_orderAA = contact_orderAA + dble(j-i)
      end if
    end do
  end do
  contact_orderSC = contact_orderSC/dble(nres*ncontactsSC)
  contact_orderAA = contact_orderAA/dble(nres*ncontactsAA)
      
! NOW WRITE OUT THE FIRST FOUR LINES OF THE OUTPUT WCN FILE

  !write(*,*) contact_order,ncontacts
  !write(*,*) contact_orderSC,ncontactsSC
  !write(*,*) contact_orderAA,ncontactsAA
  ! open output files:
    inquire(file=trim(adjustl(output)),exist=output_exists)
    if (output_exists) then
      open(unit=output_unit,file=trim(adjustl(output)),status='old')
      ios=0
      do while(ios==0)
        read(output_unit,*,iostat=ios)
      end do
    else
      open(unit=output_unit,file=trim(adjustl(output)),status='new')
      write(output_unit,'(6A25)') 'pdb','natoms','nres','contact_order','contact_orderSC','contact_orderAA'
    end if

    write(output_unit,'(1A25,2I25,3F25.5)') pdb_name,natoms,nres,contact_order,contact_orderSC,contact_orderAA

    close(output_unit)

end program get_pdb_prob_CO
