! This Fortran program attempts to find the best performing definition and free parameters of the definition of Contact Number as used in my previous work "structural prediction of ER" that results in the highest Spearman correlation between the B-factors of the representative residue atoms in pdb files (on Echave pdb data set) and CN.
! GOALS:
!       -- To see if all structures have approximately the best performing free parameters of the CN definition or not.
! INPUT:
!       -- CN definition to be used, with weighting options: h (Heaviside step function), p (power-law), e (exponential), g (Gaussian)
!       -- Path to the input file pdb_prop_CO.out containing pdb names and number of residues
!       -- Path to the input representative atomic coordinates and Bfactors (../../../properties/res_crd_bf)
!       -- Path to the input r4sJC Evolutionary Rates estimates (../../../jec_pdb_r4s.csv)
!       -- Path and name of the output file containing Spearman correlations of CN-bfac for each pdb structure on each line.
!       -- Path and name of the output summary file containing the name of the pdbs on the first column, the best performing value of free parameter of the model, the corresponding Spearman correlation.

! METHOD:
!       The code searches for the best performing cutoff_distance for CN in the stupidest and slowest, but simplest way possible, which does not matter, thanks for Fortran's miraculous runtime speed. It searches for the best value of the free parameter cutoff_distance.

! Amir Shahmoradi, Sunday 6:05 PM, March 29 2015, iCMB, UT Austin

program wcn_r4sJC_cor
implicit none
! pdb file variables:
  integer, parameter                            :: npdb=213                ! number of pdbs in pdb_prop_CO.out file
  integer                                       :: natoms, nres  
  character(len=3)                              :: res_name
  integer                                       :: res_num
  real*8          , dimension(:,:), allocatable :: crd
  real*8          , dimension(:)  , allocatable :: bfactor
  character(len=6)                              :: pdb,pdb_name            ! The name of the pdb structure and the SINGLE chain contained in it.
! input files and variables:
  character(len=1)   :: model                                              ! the wcn model and definition: h, p, e, g
  character(len=300) :: co_in                                              ! input pdb file. ATTN: Any non ATOM records will be ignored!
  character(len=300) :: crd_in                                             ! input pdb file. ATTN: Any non ATOM records will be ignored!
  character(len=300) :: jec_in                                             ! input file containing Echave's r4sJC estimates for all pdb files.
  character(len=300) :: exp_out	                                           ! output file containing on each line, the Spearman correlations of Seq.Ent of a single pdb file with different exponential definitions of WCN (different cutoff_distances).
  character(len=300) :: sum_out	                                           ! output summary file containing the pdb names and the best performing cutoff_distances and the corresponding Spearman correlations and the length of the protein on the last column.
! input files variables:
  integer            :: iostat, counter
  integer            :: co_in_unit=10,crd_in_unit=11                        ! input file units
  integer            :: jec_in_unit=12                                      ! input file units
  integer            :: exp_out_unit=21, sum_out_unit=22                    ! output file units
  integer            :: ios
  character(len=100) :: exp_output_format
! ELJ input file variables:
!  character(len=4)   :: elj_pdb_name
!  character(len=1)   :: elj_chain_id
!  integer            :: elj_res_num
!  real*8             :: elj_seqent, elj_ddgent  
!  real*8, dimension(:) , allocatable :: pdb_seqent, pdb_ddgent  
!  integer, dimension(:), allocatable :: pdb_resnum
! variables and parameters for wcn calculations and correlations:
  real*8                                  :: free_param_min, free_param_max,stride
  integer                                 :: nstride     ! number of strides
  real*8 , dimension(:)     , allocatable :: free_param,sp_cor,abs_sp_cor                       ! The array of the values of the free parameter of WCN definitions.
  real*8 , dimension(:)     , allocatable :: wcn                                                ! The array of contact numbers for each each CA atom.
  real*8                                  :: spear                                              ! function
  character(len=6)                        :: dummy_char
  !character(len=100):: wcn_output_format,res_output_format,resnum_output_format
! other variables:
  integer                                 :: i,ii,j
  real*8 :: tstart,tend

! First get the command line arguments
  if (command_argument_count()/=6) then
    write(*,*)
    write(*,*) "Incorrect number of input arguments on the command line."
    write(*,*) "Correct use:"
    write(*,'(1A180)') "./a.out <input: wcn weighting model: h, p, e, g> <input: pdb_prop_CO.out> <input: crd_bf file> <output: cor. vs. param. file> <output summary file: best performing params>"
    write(*,*)
    stop
  end if
  call get_command_argument(1,model)
  call get_command_argument(2,co_in)
  call get_command_argument(3,crd_in)
  call get_command_argument(4,jec_in)
  call get_command_argument(5,exp_out)
  call get_command_argument(6,sum_out)
  if ( model == 'h' ) then
    free_param_min = 0.d0
    free_param_max = 50.d0
    stride = 0.1d0
  elseif ( model == 'p' ) then
    free_param_min = -30.d0
    free_param_max = 30.d0
    stride = 0.1d0
  elseif ( model == 'e' ) then
    free_param_min = 0.d0
    free_param_max = 50.d0
    stride = 0.2d0
  elseif ( model == 'g' ) then
    free_param_min = 0.d0
    free_param_max = 50.d0
    stride = 0.2d0
  else
    write (*,*) 'invalid input model! : ', model
    write (*,*) 'program aborted'
    STOP
  end if

! allocate fre_param values:
  nstride = nint ( (free_param_max-free_param_min) / stride )
  allocate(free_param(nstride),sp_cor(nstride),abs_sp_cor(nstride))
  do i = 1,nstride
    free_param(i) = free_param_min + dble(i)*stride
  end do
  
! OPEN INPUT FILES:
  open(unit=co_in_unit,file=trim(adjustl(co_in)),status='old')
  open(unit=crd_in_unit,file=trim(adjustl(crd_in)),status='old')
  open(unit=jec_in_unit,file=trim(adjustl(jec_in)),status='old')
  read(co_in_unit,*)        ! read file header
  read(crd_in_unit,*)       ! read file header
! OPEN OUTPUT FILES:
  write(dummy_char,'(1I6)') nstride
  exp_output_format = trim(adjustl('(1A20,' // trim(adjustl(dummy_char)) // 'F20.6)'))
  open(unit=exp_out_unit,file=trim(adjustl(exp_out)),status='replace')
  write(exp_out_unit,exp_output_format) 'free_parameter',(free_param(i),i=1,nstride)
  open(unit=sum_out_unit,file=trim(adjustl(sum_out)),status='replace')
  write(sum_out_unit,'(4A20)') 'pdb','free_param_best','sp_cor_best','nres'
  
call cpu_time(tstart)
do ii = 1,npdb

  read(co_in_unit,*) pdb,natoms,nres
  write(*,*) "processing PDB number ", ii, " : ", pdb, " out of 213 pdbs"
  allocate(crd(3,nres),bfactor(nres),wcn(nres))
  ! Now read the crd file:
  do j = 1,nres
    read(crd_in_unit,*) pdb_name, res_name, res_num, crd(1:3,j), bfactor(j)
    if (pdb /= pdb_name) then
      write(*,*) 'Something is heavily fishy here!'
      write(*,*) 'pdb names do not match!', pdb, pdb_name; write(*,*) 'Program Aborted'
      STOP
    end if
  end do
  ! now calculate WCN for all values of parameters:
  do i = 1,nstride
    call wcn_finder(model,free_param(i),nres,crd,wcn)
    ! Now calculate the spearman correlation between wcn and the quantity of interest:
    sp_cor(i) = spear(nres,wcn,bfactor)
    if (isnan(sp_cor(i)) .or. abs(sp_cor(i)) > 1.d0) sp_cor(i) = 0.d0     ! replace NAN values with zero.
  end do
  deallocate(crd,bfactor,wcn)
  ! NOW WRITE OUT THE FIRST FOUR LINES OF THE OUTPUT WCN FILE
  write(exp_out_unit,exp_output_format) pdb,(sp_cor(i),i=1,nstride)
  abs_sp_cor = abs(sp_cor)
  write(sum_out_unit,'(1A20,2F20.3,1I20)') pdb,free_param(maxloc(abs_sp_cor)),sp_cor(maxloc(abs_sp_cor)),nres

end do
call cpu_time(tend)

write(*,*) new_line('a'), 'Overall it took ', tend-tstart, ' seconds for all 213 proteins.' 

close(co_in_unit)
close(crd_in_unit)
close(exp_out_unit)
close(sum_out_unit)

end program wcn_r4sJC_cor


! This Fortran subroutine takes in the 3D coordinates of a set of atoms. On the output it gives the WCN of all the atoms that were given in the input.
! Also in the input variables are the number of residues (nres) and the input model for wcn: h, p, e, g

subroutine wcn_finder(model,free_param,nres,crd,wcn)
  implicit none
  character(len=1), intent(in) :: model
  integer, intent(in)          :: nres
  real*8, intent(in)           :: free_param,crd(3,nres)
  real*8, intent(out)          :: wcn(nres)
  integer                      :: i,j
  real*8                       :: distance,free_param_sq
  wcn = 0.d0
  free_param_sq = free_param*free_param
  do i=1,nres
    do j=1,nres
      if (i /= j) then
        distance = sqrt( (crd(1,i)-crd(1,j))**2 + (crd(2,i)-crd(2,j))**2 + (crd(3,i)-crd(3,j))**2 )
        if (model == 'h') then
          if (distance <= free_param) wcn(i) = wcn(i) + 1.d0
        elseif (model == 'p') then
          wcn(i) = wcn(i) + distance**free_param
        elseif (model == 'e') then
          wcn(i) = wcn(i) + exp( -distance/free_param )
        elseif (model == 'g') then
          wcn(i) = wcn(i) + exp( -distance*distance/free_param_sq )
        end if
      end if
    end do
  end do
end subroutine wcn_finder

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
