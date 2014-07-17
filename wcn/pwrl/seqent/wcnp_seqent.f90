! This Fortran program attemots to find the exponent of the definition of the variable Weighted Contact Number (wcn) that results in the highest Spearman correlation between the sequence entropies of Eleisha Jackson (on Echave pdb data set) and wcn.
! GOALS:
!       -- To see if the best performing wcn exponent is pdb structure dependent, or is the same for all pdbs.
!       -- If it is the same for all pdbs, then is the average value of the exponent systematically different from the original definition of WCN (exponent of 2)?
!       -- and on the follow-up, could this possibly explain part of the diversity of the correlations seen among wcn-seq_ent?
! INPUT:
!       -- Path to the input pdb file
!       -- Path to Eleisha's qes. entropy file
!       -- Path and name of the output file containing Spearman correlations of wcn-seq.ent. for each pdb structure on each line.
!       -- Path and name of the output summary file containing the name of the pdbs on the first column, the best performing exponent, the corresponding Spearman corr., and the length of the protein, each on a separate column.
! METHOD:
!       The code searches for the best performing exponent for WCN in the stupidest and slowest, but simplest way possible, which does not matter, thanks for Fortran's miraculous runtime speed. It searches for the best exponent starting from exponent = 0 to 10, with 0.01 strides.
! Amir Shahmoradi, Wednesday 4:27 PM, June 25 2014, iCMB, UT Austin

program wcn_ent_cor
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
  integer                                       :: nres=0,natoms=0                                  ! The number of residues and atoms in the pdb file
  integer                                       :: res_num_old                                      ! A dummy variable used to recognize duplicate Amino Acid entries in the pdb file (The AAs that two different conformations resolved for).
  character(len=1)                              :: res_cod_old                                      ! A dummy variable used to recognize duplicate Amino Acid entries in the pdb file (The AAs that two different conformations resolved for).
  character(len=6)                              :: pdb_name                                         ! The name of the pdb structure and the SINGLE chain contained in it.
! input files:                                                                                      
  character(len=300)                            :: pdb_in                                           ! input pdb file. ATTN: Any non ATOM records will be ignored!
  character(len=300)                            :: elj_in                                           ! input file containing Eleisha's Seq.Ent for all pdb files. ATTN: The file should be in its original format written by Eleisha, that is, the columns ordering should have not been changes.
  character(len=300)                            :: exp_out                                          ! output file containing on each line, the Spearman correlations of Seq.Ent of a single pdb file with different definitions of WCN (different exponents).
  character(len=300)                            :: sum_out                                          ! output summary file containing the pdb names and the best performing exponents and the corresponding Spearman correlations and the length of the protein on the last column.
! input files variables:                                                                            
  logical                                       :: exp_out_exists, sum_out_exists                   ! if the output files already exist, these variables are true.
  integer                                       :: pdb_in_unit=11, elj_in_unit=12                   ! input file units
  integer                                       :: exp_out_unit=21, sum_out_unit=22                 ! output file units
  integer                                       :: ios
  character(len=100)                            :: exp_output_format
! ELJ input file variables:
  character(len=4)                              :: elj_pdb_name
  character(len=1)                              :: elj_chain_id
  integer                                       :: elj_res_num
  real*8                                        :: elj_seqent, elj_ddgent  
  real*8          , dimension(:)  , allocatable :: pdb_seqent, pdb_ddgent  
  integer         , dimension(:)  , allocatable :: pdb_resnum
! variables and parameters for wcn calculations and correlations:
  real*8          , parameter                   :: alpha_min = -30.d0                               ! alpha stands for the exponent of inter-residue distance in the definition of WCN. alpha = -2 corresponds to the original definition of WCN.
  real*8          , parameter                   :: alpha_max =  30.d0                               ! alpha stands for the exponent of inter-residue distance in the definition of WCN. alpha = -2 corresponds to the original definition of WCN.
  real*8          , parameter                   :: alpha_stride = 0.1                               ! alpha stands for the exponent of inter-residue distance in the definition of WCN. alpha = -2 corresponds to the original definition of WCN.
  real*8          , parameter                   :: alpha_default = -2.d0                            ! The default value of alpha as in the original definition of wcn.
  integer         , parameter                   :: nalpha = nint((alpha_max-alpha_min)/alpha_stride)
  real*8          , dimension(:)  , allocatable :: wcn                                              ! The array of contact numbers for each each CA atom.
  real*8          , dimension(:,:), allocatable :: CAcrd                                            ! this vector will contain the coordinates of the c-alpha carbons in the md crd file. dimension is (nres,3). 
  real*8          , dimension(nalpha)           :: alpha,sp_cor,abs_sp_cor                          ! The array of contact numbers for each each CA atom. sp_cor stands for the spearman correlation of wcn and ELJ sequence entropies. abs stands for the absolute values of the spearman correlations.
  real*8                                        :: spear                                            ! function
  real*8                                        :: sp_default                                       ! The value of Spearman correlation for the default case of alpha = -2 in WCN definition.
  character(len=6)                              :: dummy_char
! other variables:
  integer            :: i,ii,j
! First get the command line arguments

  if (command_argument_count()/=4) then
    write(*,*)
    write(*,*) "Incorrect number of input arguments on the command line."
    write(*,*) "Correct use:"
    write(*,'(1A140)') "./a.out <input: pdb file> <input: Eleisha's Seq.Ent file> <output: pdb exponent file> <output summary file: best performing exponents>"
    write(*,*)
    stop
  end if
  call get_command_argument(1,pdb_in)
  call get_command_argument(2,elj_in)
  call get_command_argument(3,exp_out)
  call get_command_argument(4,sum_out)
  
! OPEN INPUT FILES:
  open(unit=pdb_in_unit,file=trim(adjustl(pdb_in)),status='old')
  open(unit=elj_in_unit,file=trim(adjustl(elj_in)),status='old')
    
! FIRST DETERMINE THE FOUR-LETTER NAME OF THE PDB STRUCTURE:
  pdb_in = adjustr(pdb_in)
  i = len(adjustr(pdb_in))
  if (pdb_in(i-3:i) /= '.pdb' .or. pdb_in(i-10:i-10) /= '_') then
    write(*,*) 'WARNING: pdb filename does not seem to be correct'
    write(*,*) 'pdb filename: ', pdb_in(i-9:i)
    read(*,*)
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
      if (trim(adjustl(atom_name(1))) == 'CA') then
        if (nres == 0) then
          nres = nres + 1
        elseif (res_num_old == res_num(1) .and. res_cod_old == res_code(1)) then
            write(*,*) 'Duplicate Amino Acid found in structure!', res_num(1), res_name(1)
            res_num_old = res_num(1)
            res_cod_old = res_name(1)
            natoms = natoms + 1
            cycle
        else
          !write(*,*) 'Duplicate Amino Acid found in structure!', res_num(1), res_name(1)
          nres = nres + 1
        end if
        res_num_old = res_num(1)
        res_cod_old = res_name(1)
      end if
      natoms = natoms + 1
      cycle
    end if
  end do
  
  deallocate(record,atom_num,atom_name,alt_loc_ind,res_name,chain_id,res_num,res_code,crd,occupancy,bfactor)
  close(pdb_in_unit)
  
  write(*,'(1A9,1I8)') "natoms: ",natoms; write(*,'(1A9,1I8)') "nres: ",nres !; write(*,*)

! NOW READ IN THE SEQUENCE ENTROPY FROM ELEISHA'S INPUT FILE:
  
  allocate(pdb_seqent(nres),pdb_ddgent(nres),pdb_resnum(nres))
  read(elj_in_unit,*)
  ios = 0; ii = 0
  do while(ios==0)
    read(elj_in_unit,*,iostat=ios) i, elj_pdb_name, elj_chain_id, elj_res_num, elj_seqent, elj_ddgent
    if (ios<0 .and. (ii < nres)) then
      write(*,'(1A100)') 'Sum Tin Wong! : reached the end of ELJ input file, while expecting further records.'
      write(*,*) 'ELJ line: ', i
      write(*,*) 'ELJ pdb name: ', elj_pdb_name
      write(*,*) 'input pdb name: ', pdb_name
      stop
    elseif (ios<0) then
      exit
    elseif (ios>0) then
      write(*,*) 'IO stat = ', ios
      write(*,*) 'Sum Tin Wong! : input ELJ file broken at around line ~ ', i, ' for pdb ', elj_pdb_name; stop
    else
      if (pdb_name(1:4) == elj_pdb_name) then
        ii = ii + 1
        pdb_resnum(ii) = elj_res_num
        pdb_seqent(ii) = elj_seqent
        pdb_ddgent(ii) = elj_ddgent
        !write(*,'(1I10,2E15.5)') ii,elj_seqent,pdb_seqent(ii)
      end if
      if (ii>nres) then
        write(*,*) 'Sum Tin Wong in ELJ input file: '
        write(*,*) 'The number of residues in ELJ file is larger than that of the pdb file: '
        write(*,*) 'pdb nres: ', nres
        write(*,*) 'ELJ nres: ', ii
        stop
      end if
    end if
  end do
  
  
! NOW RECORD THE CA COORDINATES OF THE PDB FILE:

  open(unit=pdb_in_unit,file=trim(adjustl(pdb_in)),status='old')
  allocate(CAcrd(nres,3))
  allocate(record(natoms),atom_num(natoms),atom_name(natoms),alt_loc_ind(natoms),res_name(natoms),chain_id(natoms),res_num(natoms),res_code(natoms),crd(natoms,3),occupancy(natoms),bfactor(natoms))
  do
    read(pdb_in_unit,'(1A4)') record(1)
    if (trim(adjustl(record(1)))=='ATOM') exit
    !write(*,'(1A4)') record
    cycle
  end do
  backspace(unit=pdb_in_unit)
  
  i = 1; ii = 0
  do
    read(pdb_in_unit,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2)',iostat=ios) record(i),atom_num(i),atom_name(i),alt_loc_ind(i),res_name(i),chain_id(i),res_num(i),res_code(i),(crd(i,j),j=1,3),occupancy(i),bfactor(i)
    if (ios<0 .and. i<natoms) then
      write(*,'(1A100)') 'Sum Tin Wong! : reached the end of output PDB file, while expecting further records from the file.'; stop
    elseif (ios>0) then
      write(*,*) 'Sum Tin Wong! : input PDB file broken: ', pdb_name; stop
    elseif (trim(adjustl(record(i)))=='ATOM') then
      if (i>=natoms .and. ii/=nres) then
        write(*,*) 'i, ii: ', i,ii
        write(*,*) 'Sum Tin Wong! : input PDB file broken: ', pdb_name
        stop
      elseif (trim(adjustl(atom_name(i)))=='CA') then
        if (ii == 0) then
          ii = ii + 1
          CAcrd(ii,1:3) = crd(i,1:3)
          res_num_old = res_num(i)
          res_cod_old = res_name(i)
        elseif (res_num_old == res_num(i) .and. res_cod_old == res_code(i)) then
          write(*,*) 'WARNING: Duplicate Amino Acid found in structure!', res_num(i), res_name(i)
          i = i + 1
          cycle
        else
          ii = ii + 1
          CAcrd(ii,1:3) = crd(i,1:3)
          res_num_old = res_num(i)
          res_cod_old = res_name(i)
        end if
      end if
      i = i + 1
      if (ii == nres) exit
      cycle
    elseif (trim(adjustl(record(i)))=='TER') then
      cycle
    elseif (trim(adjustl(record(i)))=='END') then
      if (ii<=nres) then
        write(*,'(1A100)') 'Sum Tin Wong! : reached the end of output PDB file, while expecting further reading from the file.'
        write(*,*) 'number of residues read: ', i
        stop
      end if
      exit
    end if
  end do; close(pdb_in_unit)
  
  write(*,*) 'number of CA atoms in the pdb file: ', ii ; write(*,*)

! NOW CALCULATE WCN FOR EACH CA ATOM, FOR DIFFERENT EXPONENT (ALPHA) VALUES IN THE DEFINITION OF WCN AND MEASURE THE CORRESPONDING SPEARMAN CORRELAITONS OF WCN AND PDB_SEQENT.

  allocate(wcn(nres))
  call get_wcn(nres,CAcrd,alpha_default,wcn)
  sp_default = spear(nres,wcn(1:nres),pdb_seqent(1:nres))
  do i = 1,nalpha
    alpha(i) = alpha_min + dble(i)*alpha_stride
    call get_wcn(nres,CAcrd,alpha(i),wcn)
    !write(*,*) nres; read(*,*)
    !write(*,*) size(wcn),size(pdb_seqent),spear(wcn,pdb_seqent,nres)
    !write(*,*) (CAcrd(ii),ii=1,nres),'heeeyy'; read(*,*)
    sp_cor(i) = spear(nres,wcn(1:nres),pdb_seqent(1:nres))
    !write(*,*) wcn
    !write(*,*) spear(nres,wcn(1:nres),wcn(1:nres)); read(*,*)
  end do

! NOW WRITE OUT THE FIRST FOUR LINES OF THE OUTPUT WCN FILE
  ! open output files:
    write(dummy_char,'(1I6)') nalpha
    exp_output_format = trim(adjustl('(1A20,' // trim(adjustl(dummy_char)) // 'F20.6)'))
    inquire(file=trim(adjustl(exp_out)),exist=exp_out_exists)
    if (exp_out_exists) then
      open(unit=exp_out_unit,file=trim(adjustl(exp_out)),status='old')
      ios=0
      do while(ios==0)
        read(exp_out_unit,*,iostat=ios)
      end do
    else
      open(unit=exp_out_unit,file=trim(adjustl(exp_out)),status='new')
      write(exp_out_unit,exp_output_format) 'exponents',(alpha(i),i=1,nalpha)
    end if

    write(exp_out_unit,exp_output_format) pdb_name,(sp_cor(i),i=1,nalpha)
    
    inquire(file=trim(adjustl(sum_out)),exist=sum_out_exists)
    if (sum_out_exists) then
      open(unit=sum_out_unit,file=trim(adjustl(sum_out)),status='old')
      ios=0
      do while(ios==0)
        read(sum_out_unit,*,iostat=ios)
      end do
    else
      open(unit=sum_out_unit,file=trim(adjustl(sum_out)),status='new')
      write(sum_out_unit,'(6A20)') 'pdb_name','exp_best','sp_best','sp_default','sp_diff','nres'
    end if
  
  
  abs_sp_cor = abs(sp_cor)
  write(sum_out_unit,'(1A20,4F20.3,1I20)') pdb_name,alpha(maxloc(abs_sp_cor)),sp_cor(maxloc(abs_sp_cor)),sp_default,abs_sp_cor(maxloc(abs_sp_cor))-abs(sp_default),nres

end program wcn_ent_cor




! This Fortran subroutine takes in the 3D coordinates of a set of atoms (presumably C-alpha of pdb files) and also a exponent 'alpha' for the power-law definition of the weighted contact numbers of the atoms. On the output it gives the WCN of all the atoms that were given in the input.
! Amir Shahmoradi, Wednesday 3:38 PM, June 25 2014, iCMB, UT Austin

subroutine get_wcn(natoms,crd,alpha,wcn)
implicit none
integer, intent(in) :: natoms
!real*8 :: alpha_temp
real*8, intent(in) :: alpha,crd(natoms,3)
real*8, intent(out) :: wcn(natoms)
integer :: i,j
real*8 :: distance

wcn = 0.d0
!alpha_temp = -0.01d0
do i=1,natoms
  do j=1,natoms
    if (j/=i) then
      distance = sqrt((crd(i,1)-crd(j,1))**2 + (crd(i,2)-crd(j,2))**2 + (crd(i,3)-crd(j,3))**2)
      !wcn(i) = wcn(i) + distance**alpha_temp
      wcn(i) = wcn(i) + distance**alpha
    end if
  end do
end do
!write(*,*) wcn; read(*,*)
end subroutine get_wcn


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! The following codes are used for the calculation of the Spearman correlation coefficients.

function spear(n,data1,data2)
! USES betai,crank,erfcc,sort2
INTEGER, intent(in) :: n
real*8, intent(in), dimension(n)  :: data1,data2
real*8  :: d,probd,probrs,spear,zd,wksp1(n),wksp2(n)
real*8  :: aved,df,en,en3n,fac,sf,sg,t,vard,betai,erfcc
INTEGER :: j
wksp1 = data1
wksp2 = data2
call sort2(n,wksp1,wksp2)
call crank(n,wksp1,sf)
call sort2(n,wksp2,wksp1)
call crank(n,wksp2,sg)
d=0.d0
do j=1,n
  d=d+(wksp1(j)-wksp2(j))**2
end do
en=n
en3n=en**3-en
aved=en3n/6.-(sf+sg)/12.
fac=(1.-sf/en3n)*(1.-sg/en3n)
vard=((en-1.)*en**2*(en+1.)**2/36.)*fac
zd=(d-aved)/sqrt(vard)
!probd=erfcc(abs(zd)/1.4142136)
spear=(1.-(6./en3n)*(d+(sf+sg)/12.))/sqrt(fac)
fac=(1.+spear)*(1.-spear)
!if(fac.gt.0.)then
!  t=spear*sqrt((en-2.)/fac)
!  df=en-2.
!  probrs=betai(0.5*df,0.5,df/(df+t**2))
!else
!  probrs=0.
!endif
end function spear

SUBROUTINE sort2(n,arr,brr)
      INTEGER n,M,NSTACK
      real*8 arr(n),brr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      real*8 a,b,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do j=l+1,ir
          a=arr(j)
          b=brr(j)
          do i=j-1,l,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
          end do
          i=l-1
2         arr(i+1)=a
          brr(i+1)=b
        end do
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        temp=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=temp
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
          temp=brr(l)
          brr(l)=brr(ir)
          brr(ir)=temp
        endif
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
          temp=brr(l+1)
          brr(l+1)=brr(ir)
          brr(ir)=temp
        endif
        if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
          temp=brr(l)
          brr(l)=brr(l+1)
          brr(l+1)=temp
        endif
        i=l+1
        j=ir
        a=arr(l+1)
        b=brr(l+1)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        temp=brr(i)
        brr(i)=brr(j)
        brr(j)=temp
        goto 3
5       arr(l+1)=arr(j)
        arr(j)=a
        brr(l+1)=brr(j)
        brr(j)=b
        jstack=jstack+2
        if(jstack.gt.NSTACK) then
           write(*,*) 'NSTACK too small in sort2'
           stop
        end if
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
END SUBROUTINE sort2

SUBROUTINE crank(n,w,s)
      INTEGER n
      real*8 s,w(n)
      INTEGER j,ji,jt
      real*8 rank,t
      s=0.
      j=1
1     if(j.lt.n)then
        if(w(j+1).ne.w(j))then
          w(j)=j
          j=j+1
        else
          do jt=j+1,n
            if(w(jt).ne.w(j))goto 2
          end do
          jt=n+1
2         rank=0.5*(j+jt-1)
          do ji=j,jt-1
            w(ji)=rank
          end do
          t=jt-j
          s=s+t**3-t
          j=jt
        endif
      goto 1
      endif
      if(j.eq.n)w(n)=n
END SUBROUTINE crank