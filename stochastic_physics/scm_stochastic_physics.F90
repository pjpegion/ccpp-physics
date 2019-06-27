module scm_stochastic_physics

      use machine
      use mersenne_twister, only: random_setseed,random_gauss,random_stat
      use GFS_typedefs,       only: GFS_control_type, GFS_init_type,GFS_coupling_type
      implicit none
      logical sppt_sfclimit
      real :: sppt_sigtop1,sppt_sigtop2,shum_sigefold
      real, dimension(5) :: sppt,sppt_tau
      real, dimension(5) :: shum,shum_tau
      real, dimension(5) :: sfc_tau
      real(kind=kind_dbl_prec), dimension(5) :: pertz0,pertshc,pertzt
      real(kind=kind_dbl_prec), dimension(5) :: pertlai,pertvegf,pertalb
      integer nsfcpert
      integer(8) ::iseed_sfc
      logical sppt_land
      logical do_sfcperts
      integer(8) ::iseed_sppt,iseed_shum
      logical sppt_logit
      logical do_shum,do_sppt,use_zmtnblck,pert_clds

 private
 public :: init_stochastic_physics,stochastic_physics_run

!derived types
 type random_pattern
    real, allocatable, dimension(:) :: stdev
    real, allocatable, dimension(:) :: tau
    real, allocatable, dimension(:) :: phi
    real, allocatable, dimension(:) :: var

    real, allocatable, dimension(:,:) :: rnoise   ! red noise holder

    integer :: npatterns
    type(random_stat) :: rstate
 end type random_pattern

 integer :: nsppt=0
 integer :: nshum=0
 integer :: nblks=0
 real*8,allocatable :: sl(:)

 real, allocatable :: vfact_sppt(:),vfact_shum(:)
 logical :: first_call=.true.
 type(random_pattern),     save :: rpattern_sppt(1),rpattern_shum(1)

contains

subroutine init_stochastic_physics(Model,errmsg,errflg)

implicit none
type(GFS_control_type),   intent(inout) :: Model
real :: PRSI(Model%levs),PRSL(Model%levs),dx,dtp
character(len=*),         intent(out)   :: errmsg
integer,                  intent(out)   :: errflg
integer :: k
character*2::proc
errmsg = ''
errflg = 0

! replace
dtp=Model%dtp
call init_stochdata(dtp,Model%input_nml_file,Model%fn_nml,Model%nlunit)
! check to see decomposition
Model%pert_clds=pert_clds
Model%do_sppt=do_sppt
Model%use_zmtnblck=use_zmtnblck
Model%do_shum=do_shum
Model%do_sfcperts=do_sfcperts             ! mg, sfc-perts
Model%nsfcpert=nsfcpert         ! mg, sfc-perts
Model%pertz0=pertz0         ! mg, sfc-perts
Model%pertzt=pertzt         ! mg, sfc-perts
Model%pertshc=pertshc         ! mg, sfc-perts
Model%pertlai=pertlai         ! mg, sfc-perts
Model%pertalb=pertalb         ! mg, sfc-perts
Model%pertvegf=pertvegf         ! mg, sfc-perts
Model%sppt_amp=sqrt(SUM(sppt(1:nsppt)**2))
print*,'sppt_amp=',Model%sppt_amp
if ( (.NOT. do_sppt) .AND. (.NOT. do_shum) .AND. (.NOT. do_sfcperts) ) return
allocate(sl(Model%levs))
do k=1,Model%levs
   sl(k)= 0.5*(Model%ak(k)/101300.+Model%bk(k)+Model%ak(k+1)/101300.0+Model%bk(k+1)) ! si are now sigmas
enddo
if (do_sppt) then
   allocate(vfact_sppt(Model%levs))
   do k=1,Model%levs
      if (sl(k) .lt. sppt_sigtop1 .and. sl(k) .gt. sppt_sigtop2) then
         vfact_sppt(k) = (sl(k)-sppt_sigtop2)/(sppt_sigtop1-sppt_sigtop2)
      else if (sl(k) .lt. sppt_sigtop2) then
          vfact_sppt(k) = 0.0
      else
          vfact_sppt(k) = 1.0
      endif
   enddo
   if (sppt_sfclimit) then
       vfact_sppt(2)=vfact_sppt(3)*0.5
       vfact_sppt(1)=0.0
   endif
   do k=1,MOdel%levs
      print *,'sppt vert profile',k,sl(k),vfact_sppt(k)
   enddo
endif

if (do_shum) then
   allocate(vfact_shum(Model%levs))
   do k=1,Model%levs
      vfact_shum(k) = exp((sl(k)-1.)/shum_sigefold)
      if (sl(k).LT. 2*shum_sigefold) then
         vfact_shum(k)=0.0
      endif
      print *,'shum vert profile',k,sl(k),vfact_shum(k)
   enddo
endif

end subroutine init_stochastic_physics

subroutine stochastic_physics_run(Model, Coupling, errmsg, errflg)
implicit none
type(GFS_control_type),   intent(in) :: Model
type(GFS_coupling_type),  intent(inout) :: Coupling
real(kind=kind_phys),allocatable :: wnoise(:)
character(len=*),         intent(out)   :: errmsg
integer,                  intent(out)   :: errflg
         
integer :: k
integer j,ierr,i
errmsg = ''
errflg = 0
print*,'In stochastic_physics_run'
if ( (.NOT. do_sppt) .AND. (.NOT. do_shum) ) return

allocate(wnoise(1))
if (do_sppt) then
   call random_gauss(wnoise,rpattern_sppt(1)%rstate)
   if (first_call) then
      rpattern_sppt(1)%rnoise(1,1)= wnoise(1)*rpattern_sppt(1)%stdev(1)
      do k=1,Model%levs
         Coupling%sppt_wts(1,k)= rpattern_sppt(1)%rnoise(1,1)*vfact_sppt(k)
      enddo
   else
      rpattern_sppt(1)%rnoise(1,1)=rpattern_sppt(1)%phi(1)*rpattern_sppt(1)%rnoise(1,1)+  & ! AR(1)
            sqrt(1-rpattern_sppt(1)%phi(1)**2)*rpattern_sppt(1)%stdev(1)*wnoise(1)
      do k=1,Model%levs
         Coupling%sppt_wts(1,k)= rpattern_sppt(1)%rnoise(1,1)*vfact_sppt(k)
      enddo
 
   endif
   if (sppt_logit) then
      Coupling%sppt_wts(1,:) = (2./(1.+exp(Coupling%sppt_wts(1,:))))
   else
       Coupling%sppt_wts(1,:)= Coupling%sppt_wts(1,:)+1.0
   endif
endif

if (do_shum) then
   call random_gauss(wnoise,rpattern_shum(1)%rstate)
   if (first_call) then
      rpattern_shum(1)%rnoise(1,1)= wnoise(1)*rpattern_shum(1)%stdev(1)
      do k=1,Model%levs
         Coupling%shum_wts(1,k)= rpattern_shum(1)%rnoise(1,1)*vfact_shum(k)
      enddo
   else
      rpattern_shum(1)%rnoise(1,1)=rpattern_shum(1)%phi(1)*rpattern_shum(1)%rnoise(1,1)+  & ! AR(1)
            sqrt(1-rpattern_shum(1)%phi(1)**2)*rpattern_shum(1)%stdev(1)*wnoise(1)
      do k=1,Model%levs
         Coupling%shum_wts(1,k)= rpattern_shum(1)%rnoise(1,1)*vfact_shum(k)
      enddo
   endif
endif
deallocate(wnoise)
first_call=.false.

end subroutine stochastic_physics_run

subroutine stochy_namelist(sz_nml,input_nml_file,fn_nml,nlunit,deltim,iret)
      
implicit none

 
      integer,              intent(out)   :: iret
      integer,              intent(in)    :: nlunit,sz_nml
      character(len=*),     intent(in)    :: input_nml_file(sz_nml)
      character(len=64),    intent(in)    :: fn_nml
      real,                 intent(in)    :: deltim
      real tol
      integer k,ios

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      namelist /nam_stochy/sppt,sppt_tau,sppt_logit, &
      iseed_shum,iseed_sppt,shum,shum_tau,& 
      sppt_sfclimit,pert_clds, &
      sppt_sigtop1,sppt_sigtop2,&
      shum_sigefold,use_zmtnblck, &
      nsfcpert,pertz0,pertshc,pertzt,pertlai, & ! mg, sfcperts
      pertvegf,pertalb,iseed_sfc,sfc_tau,sppt_land
      tol=0.01  ! tolerance for calculations
!     spectral resolution defintion
      ! can specify up to 5 values for the stochastic physics parameters 
      ! (each is an array of length 5)
      sppt             = -999.  ! stochastic physics tendency amplitude
      shum             = -999.  ! stochastic boundary layer spf hum amp   
! logicals
      do_sppt = .false.
      pert_clds = .false.
      use_zmtnblck = .false.
      do_shum = .false.
      ! mg, sfcperts
      do_sfcperts = .false.
      sppt_land = .false.
      nsfcpert = 0
! for sfcperts random patterns
      sfc_tau     = -999.       ! time scales
      iseed_sfc   = 0           ! random seeds (if 0 use system clock)
! for SKEB random patterns.
      sppt_tau         = -999.  ! time scales
      shum_tau         = -999.
      iseed_sppt       = 0      ! random seeds (if 0 use system clock)
      iseed_shum       = 0

! parameters to control vertical tapering of stochastic physics with
! height
      sppt_sigtop1 = 0.1
      sppt_sigtop2 = 0.025
      shum_sigefold = 0.2
! reduce amplitude of sppt near surface (lowest 2 levels)
      sppt_sfclimit = .false.
      sppt_logit        = .false. ! logit transform for sppt to bounded interval [-1,+1]

#ifdef INTERNAL_FILE_NML
      read(input_nml_file, nml=nam_stochy)
#else
      rewind (nlunit)
      open (unit=nlunit, file=fn_nml, READONLY, status='OLD', iostat=ios)
      read(nlunit,nam_stochy)
#endif

      print *,' in stochy_namelist'

! PJP stochastic physics additions
      IF (sppt(1) > 0 ) THEN
        do_sppt=.true.
      ENDIF
      IF (shum(1) > 0 ) THEN
        do_shum=.true.
!     shum parameter has units of 1/hour, to remove time step
!     dependence.
!     change shum parameter units from per hour to per timestep
         DO k=1,5
            IF (shum(k) .gt. 0.0) shum(k)=shum(k)*deltim/3600.0
         ENDDO
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     ! mg, sfcperts
      IF (pertz0(1) > 0 .OR. pertshc(1) > 0 .OR. pertzt(1) > 0 .OR. &
          pertlai(1) > 0 .OR. pertvegf(1) > 0 .OR. pertalb(1) > 0) THEN
        do_sfcperts=.true.
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  All checks are successful.
!
         print *, 'stochastic physics'
         print *, ' do_sppt : ', do_sppt
         print *, ' pert_clds : ', pert_clds
         print *, ' do_shum : ', do_shum
      iret = 0
!
      return
      end subroutine stochy_namelist 

! set up and initialize stochastic random patterns.


 subroutine init_stochdata(delt,input_nml_file,fn_nml,nlunit)

! initialize random patterns.  A spinup period of spinup_efolds times the
! temporal time scale is run for each pattern.  
   integer, intent(in) :: nlunit
   character(len=*),  intent(in) :: input_nml_file(:)
   character(len=64), intent(in) :: fn_nml
   real, intent(in) :: delt

   integer :: stochlun,iret,n
   print*,'in init stochdata'
   call stochy_namelist(size(input_nml_file,1),input_nml_file(:),fn_nml,nlunit,delt,iret)
   if (do_sppt.EQ. .false. .AND. do_shum.EQ. .false.) return
   stochlun=99
   iret=0
! determine number of random patterns to be used for each scheme.
   do n=1,size(sppt)
     if (sppt(n) > 0) then
        nsppt=nsppt+1
     else
        exit
     endif
   enddo
   print *,'nsppt = ',nsppt
   do n=1,size(shum)
     if (shum(n) > 0) then
        nshum=nshum+1
     else
        exit
     endif
   enddo
   print *,'nshum = ',nshum

   if (nsppt > 0) then
      print *, 'Initialize random pattern for SPPT',nsppt,sppt_tau(1:nsppt)
      call stochy_patterngenerator_init(delt,sppt_tau(1:nsppt),sppt(1:nsppt),iseed_sppt,rpattern_sppt(1),nsppt)
   endif
   if (nshum > 0) then
      print *, 'Initialize random pattern for SHUM',nshum,shum_tau(1:nshum)
      call stochy_patterngenerator_init(delt,shum_tau(1:nshum),shum(1:nshum),iseed_shum,rpattern_shum(1),nshum)
   endif
end subroutine init_stochdata


 subroutine stochy_patterngenerator_init(delt, tscale, stdev, iseed, rpattern, npatterns)
   real, intent(in) :: tscale(npatterns),stdev(npatterns)
   real, intent(in) :: delt
   integer, intent(in) :: npatterns
   integer(8), intent(inout) :: iseed
   type(random_pattern), intent(out) :: rpattern
! locals
   integer(8) count, count_rate, count_max, count_trunc
   integer(8) :: iscale = 10000000000
   integer count4,seed
    
!  propagate seed supplied from namelist to all patterns... (can remove PJP)
   ! seed computed on root, then bcast to all tasks and set.
   if (iseed == 0) then
     ! generate a random seed from system clock and ens member number
     call system_clock(count, count_rate, count_max)
     ! iseed is elapsed time since unix epoch began (secs)
     ! truncate to 4 byte integer
     count_trunc = iscale*(count/iscale)
     count4 = count - count_trunc !+ member_id
     print *,'using seed',count4
   else
     count4 = mod(iseed + 2147483648, 4294967296) - 2147483648
     print *,'using seed',count4,iseed!,member_id
   endif
   seed = count4
   rpattern%npatterns=npatterns 
   !print*,'allocated spec_coeffu',nlm,2,npatterns
   allocate(rpattern%tau(npatterns))
   allocate(rpattern%stdev(npatterns))
   allocate(rpattern%phi(npatterns))
   allocate(rpattern%rnoise(1,1))
   rpattern%tau(1:npatterns) = tscale(1:npatterns)
   rpattern%phi(1:npatterns) = exp(-delt/tscale(1:npatterns))
   rpattern%stdev(1:npatterns) = stdev(1:npatterns)
   call random_setseed(seed,rpattern%rstate)
 end subroutine stochy_patterngenerator_init
end module scm_stochastic_physics 
