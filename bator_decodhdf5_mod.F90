MODULE BATOR_DECODHDF5_MOD

!**************************************************************************************************
! Source Filename  | bator_decodhdf5_mod.F90
! Author           | F. GUILLAUME (CNRM/GMAP/OBS)
! Creation Date    | 31/07/2016
!--------------------------------------------------------------------------------------------------
! Description : decoding HDF5 files module
!--------------------------------------------------------------------------------------------------
! Subroutines :
!              InitHdf5               : initialize and check fortran library
!              PrefetchHdf5           : check HDF5 file type - evaluate nbobs & nbwagons
!              - ValidOdim            : add some checks for ODIM data
!              ExpandHdf5File         : call for decoding subroutines
!              PrintMessages          : 
!              GetDatasetDataspace    : get number of elements (up to 3D) stored in a dataset
!              GetData                : check if data are readable and get them
!              ReadData               : interface for the following subroutines :
!              - ReadDataInt1         : read scalar or 1D array integer data
!              - ReadDataInt2         : read 2D array integer data
!              - ReadDataInt3         : read 3D array integer data
!              - ReadDataReal1        : read scalar or 1D array real data
!              - ReadDataReal2        : read 2D array real data
!              - ReadDataReal3        : read 3D array real data
!              GetAttribute           : check whether attribute is readable and get it
!              ReadAttribute          : interface for the following subroutines :
!              - ReadAttributeStr1    : read scalar or 1D array string attribute
!              - ReadAttributeStr2    : read 2D array string attribute
!              - ReadAttributeInt1    : read scalar or 1D array integer attribute
!              - ReadAttributeInt2    : read 2D integer attribute
!              - ReadAttributeReal1   : read scalar or 1D array real attribute
!              - ReadAttributeReal2   : read 2D array real attribute
!              Hdf5Type               : select the Atomic needed to read an attribute or data
!              Odim                   : decode hdf5 odim files. contains :
!              - GetDAttributes       : get needed attributes (from head, dataset and data groups)
!              - GetFAttributes       : get needed attributes from 'Quality' groups
!              - BuildAttributes      : 
!              Mtvza                  : decode hdf5 mtvza files
!              maps                   : maps bins to get a profile
!**************************************************************************************************

USE PARKIND1, ONLY           : jpim, jpib, jprb
USE YOMHOOK, ONLY            : lhook, Dr_Hook
USE VARNO_MODULE, ONLY       : Varno
USE YOMCOCTP, ONLY           : NRADAR, NRADA1, NSATEM, NGTHRB
USE BATOR_MODULE, ONLY       : batout, rabsi,  zent,  zentsup, zwagon, nciotp, ncirfl, ncioch, clsid, &
                             & ncialt, ncilon, ncidat, ncilat, ncinlv, ncietm, ncistd, rabso,  nabso, &
                             & nabsi,  abor1, INbTypeHdf5, HdfConfig, tref_ficobs, HMtvza, HOdim,     &
                             & NbMaxMtvzaChan, LabelLength, TREF_RADAR, JPXRADS, ncilet, radians
USE BATOR_DATETIME_MOD, ONLY : TDate, NewDate, AddDate, ReturnDateArray, VerifDate, DiffDate
USE BATOR_UTIL_MOD, ONLY     : BATOR_FILTER_RADAR, BATOR_RADAR_WIND_CLEANER, UVCOM
USE YOMCST, ONLY             : rpi, ra
USE HDF5
USE ISO_C_BINDING

IMPLICIT NONE
PRIVATE

PUBLIC :: InitHdf5, PrefetchHdf5, ExpandHdf5File

  character(len=LabelLength),dimension(:),pointer   :: String1Buf    => NULL()
  character(len=LabelLength),dimension(:,:),pointer :: String2Buf    => NULL()
  real(kind=jprb),dimension(:),pointer              :: Real1Buf      => NULL()
  real(kind=jprb),dimension(:,:),pointer            :: Real2Buf      => NULL()
  real(kind=jprb),dimension(:,:,:),pointer          :: Real3Buf      => NULL()
  integer(kind=jpib),dimension(:),pointer           :: Inte1Buf      => NULL()
  integer(kind=jpib),dimension(:,:),pointer         :: Inte2Buf      => NULL()
  integer(kind=jpib),dimension(:,:,:),pointer       :: Inte3Buf      => NULL()
  integer(kind=HID_T)                               :: HdfIntKindType, HdfRealKindType
  
INTERFACE ReadAttribute
  MODULE PROCEDURE ReadAttributeStr1,ReadAttributeStr2,ReadAttributeInt1,ReadAttributeInt2,ReadAttributeReal1,ReadAttributeReal2
END INTERFACE
INTERFACE ReadData
  MODULE PROCEDURE ReadDataInt1,ReadDataInt2,ReadDataInt3,ReadDataReal1,ReadDataReal2,ReadDataReal3
END INTERFACE

CONTAINS

!**************************************************************************************************
! Subroutine/Function | InitHdf5
! Source Filename     | bator_decodhdf5_mod.F90
! Author              | F. Guillaume (CNRM/GMAP/OBS)
! Creation Date       | 31/07/2016
!--------------------------------------------------------------------------------------------------
! Description : initialize fortran interface for HDF5 library, and check its version is >=1.8.8.
!               or close fortran interface.
!--------------------------------------------------------------------------------------------------
! Global Variables : none
!--------------------------------------------------------------------------------------------------
! History :
!         - 03/2017 Guillaume F. : 1.0 release version
!         - 08/2017 Guillaume F. : add closing fortran interface fonctionality.
!**************************************************************************************************
SUBROUTINE InitHdf5(Fct)
  integer(kind=jpim),intent(in) :: Fct
  integer(kind=jpim)            :: MajNum, MinNum, RelNum, Error
  real(kind=jprb)               :: zhook_handle

  if (lhook) call dr_hook('InitHdf5',0,zhook_handle)

  select case (Fct)
  case(0)
    call h5open_f(Error)
    if (Error == 0) then
      call h5get_libversion_f(MajNum,MinNum,RelNum,Error)
      if (Error == 0 .and. MajNum >=1 .and. ((MinNum == 8 .and. RelNum > 7) .or. MinNum >8)) then
        write(batout,'("*** INFO - BATOR : HDF5 Fortran interface initialized")')
        write(batout, '("*** INFO - BATOR : HDF5 library V.",I1,".",I1,".",I2)') MajNum,MinNum,RelNum
        HdfIntKindType  = h5kind_to_type(jpib,H5_INTEGER_KIND)
        HdfRealKindType = h5kind_to_type(jprb,H5_REAL_KIND)
        call h5eset_auto_f(0,Error)                                    ! automatic error printing off 
      else
        call Abor1('** ERROR -- BATOR : HDF5 library version must be >= 1.8.8.')
      endif
    else
      call Abor1('** ERROR -- BATOR : cannot initialize HDF5 Fortran interface.')
    endif
  case(1)
    call h5close_f(Error)
    write(batout,'(/,"*** INFO - BATOR : HDF5 Fortran interface closed",/)')
  end select
 
  if (lhook) call dr_hook('InitHdf5',1,zhook_handle)
END SUBROUTINE InitHdf5

!**************************************************************************************************
! Subroutine/Function | PrefetchHdf5
! Source Filename     | bator_decodhdf5_mod.F90
! Author              | F. Guillaume (CNRM/GMAP/OBS)
! Creation Date       | 21/07/2016
!--------------------------------------------------------------------------------------------------
! Description : check HDF5 file type and structure (only first level elements).
!               Compute number of observations for allocating. 
!--------------------------------------------------------------------------------------------------
! Global Variables : 
!--------------------------------------------------------------------------------------------------
! History :
!         - 03/2017 Guillaume F. : 1.0 release version
!**************************************************************************************************
SUBROUTINE PrefetchHdf5(Kfic)
  integer(kind=jpim),intent(in)            :: Kfic
  integer(kind=HID_T)                      :: FileId, GroupId
  integer(kind=jpim)                       :: i, j, k
  integer(kind=jpim)                       :: Error, Iret
  integer(kind=jpim)                       :: NGrp, NDset, NAtt
  integer(kind=jpim)                       :: FoundAttributes, FoundGroups, FoundDatasets
  integer(kind=jpim)                       :: WaitedAtomic, WaitedRank, NMembers, ObjectType
  integer(kind=jpib)                       :: XValues, YValues, ZValues     
  logical                                  :: Available
  real(kind=jprb)                          :: zhook_handle
  character(len=LabelLength)               :: NomMembre, NomGrpData
  character(len=17)                        :: NomFic

! for ODIM only
  integer(kind=jpim)                       :: MaxRays, NGrpData
  real(kind=jprb),dimension(:),allocatable :: ElevationsList 
  real(kind=jprb)                          :: MaxLength

  if (lhook) call dr_hook('PrefetchHdf5',0,zhook_handle)

  write(NomFic,'(A,".",A)') trim(tref_ficobs(Kfic)%format),trim(tref_ficobs(Kfic)%nomfic)
  call h5fopen_f(trim(Nomfic),H5F_ACC_RDONLY_F,FileId,Error)
  if (Error == 0) then
    i         = 0
    Available = .TRUE.

! control all general attributes, first level groups and datasets (see param.cfg)
    do while (i < INbTypeHdf5 .and. .not.tref_ficobs(Kfic)%Valide)
      i                          = i + 1
      tref_ficobs(Kfic)%type     = hdfconfig(i)%ssensor
      tref_ficobs(Kfic)%Template = i
      FoundAttributes            = 0
      FoundDatasets              = 0
      FoundGroups                = 0
      NGrp                       = 0
      NDset                      = 0
      call h5aget_num_attrs_f(FileId,NAtt,Error)
      if (Error /= 0) then
        write(batout,'("** ERROR - BATOR : PrefetchHdf5() cannot get number of general attibutes")')
      else
        tref_ficobs(Kfic)%Valide = .TRUE.
        call h5gn_members_f(FileId,"/",NMembers,Error)
        if (Error /= 0) then
          write(batout,'("** ERROR - BATOR : PrefetchHdf5() cannot get number of first level objects")')
          tref_ficobs(Kfic)%Valide = .FALSE.
        else
          do j=0, NMembers - 1
            call h5gget_obj_info_idx_f(FileId,"/",j,NomMembre,ObjectType,Error)
            if (ObjectType == H5G_GROUP_F .and. Error == 0)   NGrp  = NGrp + 1
            if (ObjectType == H5G_DATASET_F .and. Error == 0) NDset = NDset + 1
          enddo
        endif
      endif
  
      if (tref_ficobs(Kfic)%Valide) then
        do k=1, size(hdfconfig(i)%GenAttrib(:))
          call h5aexists_f(FileId,hdfconfig(i)%GenAttrib(k)%Label,Available,Error)
          if (Available .and. Error == 0) then
            FoundAttributes = FoundAttributes + 1
          endif
        enddo
        if (FoundAttributes /= NAtt) tref_ficobs(Kfic)%Valide = .FALSE.
      endif

      if (tref_ficobs(Kfic)%Valide) then
        do k=1, size(hdfconfig(i)%Group(:))
          do j=0, NMembers - 1
            call h5gget_obj_info_idx_f(FileId,"/",j,NomMembre,ObjectType,Error)
            if (ObjectType == H5G_GROUP_F .and. Error == 0) then
              if (hdfconfig(i)%Group(k)%Num < 0) NomMembre = NomMembre(1:len(trim(hdfconfig(i)%Group(k)%Label)))
              if (trim(NomMembre) == trim(hdfconfig(i)%Group(k)%Label)) then
                FoundGroups = FoundGroups + 1
              endif
            endif
          enddo
        enddo
        if (FoundGroups /= NGrp) tref_ficobs(Kfic)%Valide = .FALSE.
      endif

      if (tref_ficobs(Kfic)%Valide) then
        do k=1, size(hdfconfig(i)%Dataset(:))
          do j=0, NMembers - 1
            call h5gget_obj_info_idx_f(FileId,"/",j,NomMembre,ObjectType,Error)
            if (ObjectType == H5G_DATASET_F .and. Error == 0) then
              if (hdfconfig(i)%Dataset(k)%Num < 0) NomMembre = NomMembre(1:len(trim(hdfconfig(i)%Dataset(k)%Label)))
              if (trim(NomMembre) == trim(hdfconfig(i)%Dataset(k)%Label)) then
                FoundDatasets = FoundDatasets + 1
              endif
            endif
          enddo
        enddo
        if (FoundDatasets /= NDset) tref_ficobs(Kfic)%Valide = .FALSE.
      endif
    enddo

    if (.not.tref_ficobs(Kfic)%Valide) then
      write(batout,'("** ERROR - BATOR : PrefetchHdf5() file ",A," skipped because unknown")') trim(NomFic)
      tref_ficobs(Kfic)%Valide = .FALSE.
      tref_ficobs(Kfic)%NbWag  = 0
    else
      select case (trim(hdfconfig(i)%ssensor))
      case ('odim')
        WaitedRank   = 0
        WaitedAtomic = Hdf5Type(0)
        call GetAttribute(FileId,HODIM%ConventionName,WaitedRank,WaitedAtomic,Iret)
        if (Iret /= 0) then 
          write(batout,'("** ERROR - BATOR : PrefetchHdf5() file ",A," Attribute Conventions missing.")') trim(NomFic)
          tref_ficobs(Kfic)%Valide = .FALSE.
        elseif (.not.any(HODIM%AllowedConventions == String1Buf(1))) then
          write(batout,'("** ERROR - BATOR : PrefetchHdf5() file ",A," unknown Conventions ",A)') trim(NomFic),trim(String1Buf(1))
          tref_ficobs(Kfic)%Valide = .FALSE.
        endif
        if (associated(String1Buf)) deallocate(String1Buf)
        MaxRays           =  0
        MaxLength         =  0._jprb
        if (tref_ficobs(Kfic)%Valide) then
          allocate(ElevationsList(NMembers),STAT=Error)
          if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate ElevationsList(:)")
          ElevationsList(:) = rabsi
          call ValidOdim(FileId," ",MaxRays,MaxLength,ElevationsList)
          do j=0, NMembers - 1
            call h5gget_obj_info_idx_f(FileId,"/",j,NomMembre,ObjectType,Error)
            if (ObjectType == H5G_GROUP_F .and. Error == 0 .and. &
              & NomMembre(1:len(trim(HODIM%GrpElevName))) == trim(HODIM%GrpElevName)) then
              call ValidOdim(FileId,NomMembre,MaxRays,MaxLength,ElevationsList)
              call h5gn_members_f(FileId,trim(NomMembre),NGrpData,Error)
              do k=0, NGrpData - 1
                call h5gget_obj_info_idx_f(FileId,trim(NomMembre),k,NomGrpData,ObjectType,Error)
                if (ObjectType == H5G_GROUP_F .and. Error == 0 .and. &
                  & NomGrpData(1:len(trim(HODIM%GrpParamName))) == trim(HODIM%GrpParamName)) then
                  call ValidOdim(FileId,trim(NomMembre)//"/"//NomGrpData,MaxRays,MaxLength,ElevationsList)
                endif
              enddo
            endif
          enddo
        endif
        TREF_FICOBS(kfic)%nbobs = MaxRays * nint(MaxLength/HODIM%Resolution) ! nint peut majorer le nombre
        TREF_FICOBS(kfic)%nbwag = TREF_FICOBS(kfic)%nbobs * TREF_FICOBS(kfic)%nbwag * count(ElevationsList /= rabsi)
        if (TREF_FICOBS(kfic)%nbwag == 0) tref_ficobs(Kfic)%Valide = .FALSE.
        if (allocated(ElevationsList)) deallocate(ElevationsList)

! calcul du nombre d'obs/wagon pour METEOR en prenant le dataset Latitude en référence.
      case ('mtvza')
        do j=0, NMembers - 1
          call h5gget_obj_info_idx_f(FileId,"/",j,NomMembre,ObjectType,Error)
          if (ObjectType == H5G_DATASET_F .and. Error == 0) then
            do k=lbound(HMTVZA,1), ubound(HMTVZA,1)
              if ( trim(HMTVZA(k)%DatasetNameRoot) == NomMembre(1:len(trim(HMTVZA(k)%DatasetNameRoot)))) exit
            enddo
            if (HMTVZA(k)%Lat == NomMembre(1:len(trim(HMTVZA(k)%Lat)))) then
              call GetDatasetDataspace(FileId,trim(NomMembre),XValues,YValues,ZValues,Iret)
              if (Iret == 0) then
                TREF_FICOBS(kfic)%nbobs = XValues * YValues * ZValues
                TREF_FICOBS(kfic)%nbwag = TREF_FICOBS(kfic)%nbwag * TREF_FICOBS(kfic)%nbobs
                if (TREF_FICOBS(kfic)%nbwag == 0) tref_ficobs(Kfic)%Valide = .FALSE.
               else
                write(batout,'("** ERROR - BATOR : PrefetchHdf5() file ",A," skipped because no OBS found")') trim(NomFic)
                tref_ficobs(Kfic)%Valide = .FALSE.
              endif
            endif
          endif
        enddo
      end select
    endif
    call h5fclose_f(FileId,Error)
  else
    write(batout,'("** ERROR - BATOR : PrefetchHdf5() file ",A," cannot be opened")') trim(NomFic)
    tref_ficobs(Kfic)%Valide = .FALSE.
  endif

  if (lhook) call dr_hook('PrefetchHdf5',1,zhook_handle)
END SUBROUTINE PrefetchHdf5

SUBROUTINE ValidOdim(IdIn,NomMembre,MaxRays,MaxLength,ElevationsList)
  integer(kind=HID_T),intent(in)             :: IdIn
  character(len=*),intent(in)                :: NomMembre
  integer(kind=jpim),intent(inout)           :: MaxRays 
  real(kind=jprb),intent(inout)              :: MaxLength
  real(kind=jprb),dimension(:),intent(inout) :: ElevationsList 
  integer(kind=HID_T)                        :: GroupId
  integer(kind=jpim)                         :: Error, Iret
  integer(kind=jpim)                         :: WaitedAtomic, WaitedRank
  integer(kind=jpim)                         :: NBins
  real(kind=jprb)                            :: RScale
  real(kind=jprb)                            :: zhook_handle

  if (lhook) call dr_hook('ValidOdim',0,zhook_handle)

  call h5gopen_f(IdIn,trim(NomMembre)//"/"//trim(HODIM%GrpWhereName),GroupId,Error)
  if (Error == 0) then
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(2)
    call GetAttribute(GroupId,HODIM%ElevName,WaitedRank,WaitedAtomic,Iret)
    if (Iret == 0) then
      if (.not.any(ElevationsList == Real1Buf(1))) then
        ElevationsList(count(ElevationsList /= rabsi)+1) = Real1Buf(1)
        deallocate(Real1Buf)
      endif
    endif
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(1)
    call GetAttribute(GroupId,HODIM%NRaysName,WaitedRank,WaitedAtomic,Iret)
    if (Iret == 0) then
      MaxRays = max(MaxRays,inte1buf(1))
      deallocate(Inte1Buf)
    endif
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(1)
    call GetAttribute(GroupId,HODIM%NBinsName,WaitedRank,WaitedAtomic,Iret)
    if (Iret == 0) then 
      NBins   = Inte1Buf(1)
      deallocate(Inte1Buf)
    endif
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(2)
    call GetAttribute(GroupId,HODIM%RScaleName,WaitedRank,WaitedAtomic,Iret)
    if (Iret == 0) then 
      RScale = Real1Buf(1)
      deallocate(real1buf)
    endif
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(2)
    call GetAttribute(groupId,HODIM%RStartName,WaitedRank,WaitedAtomic,Iret)
    if (Iret == 0) then
      if (Real1Buf(1) == 0._jprb) NBins = NBins + 1
      MaxLength = max(Real1Buf(1)*1000._jprb + RScale * NBins, MaxLength)
      deallocate(Real1Buf)
    endif
    call h5gclose_f(GroupId,Error)
  endif

  if (lhook) call dr_hook('ValidOdim',1,zhook_handle)
END SUBROUTINE ValidOdim


!**************************************************************************************************
! Subroutine/Function | ExpandHdf5File
! Source Filename     | bator_decodhdf5_mod.F90
! Author              | F. Guillaume (CNRM/GMAP/OBS)
! Creation Date       | 09/08/2016
!--------------------------------------------------------------------------------------------------
! Description : open HDF5 files and launch required decoding subroutine according to the sensor.
!--------------------------------------------------------------------------------------------------
! Global Variables : 
!--------------------------------------------------------------------------------------------------
! History :
!         - 04/2017 Guillaume F. : 1.0 release version
!**************************************************************************************************
SUBROUTINE ExpandHdf5File(Filename,Kobs,Kw,Kfic,TabSlots,AnalysisDate)
  type(TDate),dimension(:),intent(in) :: TabSlots
  type(TDate),intent(in)              :: AnalysisDate
  character(len=*),intent(in)         :: Filename
  integer(kind=jpim),intent(in)       :: Kfic
  integer(kind=jpim),intent(inout)    :: Kobs
  integer(kind=jpib),intent(inout)    :: Kw
  integer(kind=jpim)                  :: Error
  integer(kind=HID_T)                 :: FileId
  integer(kind=jpim), dimension(50)   :: iterr
  logical                             :: Conformity
  real(kind=jprb)                     :: zhook_handle

  if (lhook) call dr_hook('ExpandHdf5File',0,zhook_handle)

  iterr(:)     = 0

  call h5fopen_f(trim(Filename),H5F_ACC_RDONLY_F,FileId,Error)
  IF (Error /= 0) CALL Abor1('** ERROR - BATOR : ExpandHdf5File() Pb ouverture ')

  select case (trim(tref_ficobs(Kfic)%type))
  case ('odim')
    call odim(FileId,TabSlots,AnalysisDate,Kobs,Kw,Kfic,Conformity)
  case ('mtvza')
    call mtvza(FileId,TabSlots,Kobs,Kw,iterr,Conformity)
  end select
  call h5fclose_f(FileId,Error)
  If (.not. Conformity)  write(batout,'("WARNING - BATOR : file ",A," not conform.")') trim(FileName)
  call PrintMessages(HdfConfig(tref_ficobs(kfic)%Template)%ssensor,iterr)

  if (lhook) call dr_hook('ExpandHdf5File',1,zhook_handle)
END SUBROUTINE ExpandHdf5File

!**************************************************************************************************
! Subroutine/Function | PrintMessages
! Source Filename     | bator_decodehdf5_mod.F90
! Author              | F. Guillaume (CNRM/GMAP/OBS)
! Creation Date       | 06/09/2017
!--------------------------------------------------------------------------------------------------
! Description : print warning and information messages from decoding subroutines.
!--------------------------------------------------------------------------------------------------
! History :
!         - 09/2017 Guillaume F. : 1.0 release version
!**************************************************************************************************
SUBROUTINE PrintMessages(sensor,iterr)
  character(len=16),                 intent(in) :: sensor
  integer(kind=jpim), dimension(50), intent(in) :: iterr
  real(kind=jprb)                               :: zhook_handle

  if (lhook) call dr_hook('PrintMessages',0,zhook_handle)

  select case (sensor)
  case ('mtvza')
    if (iterr(1) > 0) write(batout,'(5X,"WARNING - Incorrect Lat/lon :",15X,I9)') iterr(1)
    if (iterr(2) > 0) write(batout,'(5X,"WARNING - Date outside interval :",11X,I9)') iterr(2)
  case default
  end select

  if (lhook) call dr_hook('PrintMessages',1,zhook_handle)
END SUBROUTINE PrintMessages

!**************************************************************************************************
! Subroutine/Function | GetData
! Source Filename     | bator_decodhdf5_mod.F90
! author              | F. Guillaume (CNRM/GMAP/OBS)
! creation date       | 31/07/2016
!--------------------------------------------------------------------------------------------------
! Description : initialization to get data from a dataset.
!--------------------------------------------------------------------------------------------------
! History :
!         - 04/2017 Guillaume F. : 1.0 release version
!**************************************************************************************************
SUBROUTINE GetData(IdIn,DataName,WaitedRank,WaitedAtomic,Iret)
  integer(kind=HID_T),intent(in)                 :: IdIn
  character(len=*),  intent(in)                  :: DataName
  integer(kind=jpim),intent(in)                  :: WaitedRank
  integer(kind=jpim),intent(in)                  :: WaitedAtomic
  integer(kind=jpim),intent(out)                 :: Iret
  integer(kind=jpim)                             :: Error
  integer(kind=jpim)                             :: CurrentRank
  integer(kind=jpim)                             :: CurrentAtomic
  integer(kind=HID_T)                            :: DatasetId
  integer(kind=HID_T)                            :: TypeId
  integer(kind=HID_T)                            :: SpaceId
  integer(kind=HSIZE_T),dimension(:),allocatable :: Dimensions
  integer(kind=HSIZE_T),dimension(:),allocatable :: MaxDimensions
  logical                                        :: IsSimple
  real(kind=jprb)                                :: zhook_handle
 
  if (lhook) call dr_hook('GetData',0,zhook_handle)
  Iret = 0
  call h5dopen_f(IdIn,DataName,DatasetId,Error)
  call h5dget_space_f(DatasetId,SpaceId,Error)
  call h5sis_simple_f(SpaceId,IsSimple,Error)
  if (Error /= 0 .or. .not.IsSimple) then
    Iret = 1
  else
    call h5sget_simple_extent_ndims_f(SpaceId,CurrentRank,Error)
    call h5dget_type_f(DatasetId,TypeId,Error)
    call h5tget_class_f(TypeId,CurrentAtomic,Error)
    if (WaitedRank /= CurrentRank .or. (WaitedAtomic /= NABSO .and. WaitedAtomic /= CurrentAtomic)) then
      Iret = 2
    else
      if (CurrentRank == 0) CurrentRank = 1       ! a scalar value <=> 1D array values
      allocate(Dimensions(CurrentRank),MaxDimensions(CurrentRank),STAT=Error)
      if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate Dimensions(:) or MaxDimensions(:)")
      select case (WaitedRank)
      case (0)
        Dimensions(CurrentRank)    = 1
        MaxDimensions(CurrentRank) = 1
        if (CurrentAtomic == H5T_INTEGER_F) call ReadData(DatasetId,Inte1Buf,Dimensions,MaxDimensions,WaitedAtomic,Iret)
        if (CurrentAtomic == H5T_FLOAT_F)   call ReadData(DatasetId,Real1Buf,Dimensions,MaxDimensions,Iret)
      case (1)
        call h5sget_simple_extent_dims_f(SpaceId,Dimensions,MaxDimensions,Error)
        if (CurrentAtomic == H5T_INTEGER_F) call ReadData(DatasetId,Inte1Buf,Dimensions,MaxDimensions,WaitedAtomic,Iret)
        if (CurrentAtomic == H5T_FLOAT_F)   call ReadData(DatasetId,Real1Buf,Dimensions,MaxDimensions,Iret)
      case (2)
        call h5sget_simple_extent_dims_f(SpaceId,Dimensions,MaxDimensions,Error)
        if (CurrentAtomic == H5T_INTEGER_F) call ReadData(DatasetId,Inte2Buf,Dimensions,MaxDimensions,WaitedAtomic,Iret)
        if (CurrentAtomic == H5T_FLOAT_F)   call ReadData(DatasetId,Real2Buf,Dimensions,MaxDimensions,Iret)
      case (3)
        call h5sget_simple_extent_dims_f(SpaceId,Dimensions,MaxDimensions,Error)
        if (CurrentAtomic == H5T_INTEGER_F) call ReadData(DatasetId,Inte3Buf,Dimensions,MaxDimensions,WaitedAtomic,Iret)
        if (CurrentAtomic == H5T_FLOAT_F)   call ReadData(DatasetId,Real3Buf,Dimensions,MaxDimensions,Iret)
      case default
        Iret = 4
      end select
      if (allocated(Dimensions)) deallocate(Dimensions)
      if (allocated(MaxDimensions)) deallocate(MaxDimensions)
    endif
    call h5tclose_f(TypeId,Error)
  endif
  call h5sclose_f(SpaceId,Error)
  call h5dclose_f(DatasetId,Error)
  if (lhook) call dr_hook('GetData',1,zhook_handle)
END SUBROUTINE GetData


!**************************************************************************************************
! Subroutine/Function | GetDatasetDataspace
! Source Filename     | bator_decodhdf5_mod.F90
! author              | F. Guillaume (CNRM/GMAP/OBS)
! creation date       | 21/04/2017
!--------------------------------------------------------------------------------------------------
! Description : get number of elements (scalar, 1D, 2D or 3D) stored in a dataset (x*y*z).
!--------------------------------------------------------------------------------------------------
! History :
!         - 03/2017 Guillaume F. : 1.0 release version
!**************************************************************************************************
SUBROUTINE GetDatasetDataspace(IdIn,DatasetName,XValues,YValues,ZValues,Iret)
  integer(kind=HID_T),intent(in)                 :: IdIn
  character(len=*),  intent(in)                  :: DatasetName
  integer(kind=jpim),intent(out)                 :: Iret
  integer(kind=jpib),intent(out)                 :: XValues, YValues, ZValues
  integer(kind=jpim)                             :: CurrentRank, VirtualRank
  integer(kind=jpim)                             :: Error
  integer(kind=HID_T)                            :: DatasetId
  integer(kind=HID_T)                            :: SpaceId
  integer(kind=HSIZE_T),dimension(:),allocatable :: Dimensions
  integer(kind=HSIZE_T),dimension(:),allocatable :: MaxDimensions
  logical                                        :: IsSimple
  real(kind=jprb)                                :: zhook_handle

  if (lhook) call dr_hook('GetDatasetDataspace',0,zhook_handle)
  XValues = 0
  YValues = 0
  ZValues = 0
  Iret    = 0
  call h5dopen_f(IdIn,DatasetName,DatasetId,Error)
  if (Error == 0) then
    call h5dget_space_f(DatasetId,SpaceId,Error)
    call h5sis_simple_f(SpaceId,IsSimple,Error)
    if (Error /= 0 .OR. .NOT.IsSimple) then
      Iret = 1
    else
      call h5sget_simple_extent_ndims_f(SpaceId,CurrentRank,Error)
      VirtualRank = CurrentRank
      if (CurrentRank == 0) VirtualRank = 1       ! a scalar value <=> 1D array values
      allocate(Dimensions(VirtualRank),MaxDimensions(VirtualRank),STAT=Error)
      if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate Dimensions(:) or MaxDimensions(:)")
      select case (CurrentRank)
      case (0)
        XValues = 1
        YValues = 1
        Zvalues = 1
      case (1)
        call h5sget_simple_extent_dims_f(SpaceId,Dimensions,MaxDimensions,Error)
        XValues = Dimensions(1)
        YValues = 1
        Zvalues = 1
      case (2)
        call h5sget_simple_extent_dims_f(SpaceId,Dimensions,MaxDimensions,Error)
        Xvalues = Dimensions(1)
        YValues = Dimensions(2)
        Zvalues = 1
      case (3)
        call h5sget_simple_extent_dims_f(SpaceId,Dimensions,MaxDimensions,Error)
        Xvalues = Dimensions(1)
        YValues = Dimensions(2)
        ZValues = Dimensions(3)
      case default
        Iret = 4
      end select
      if (allocated(Dimensions)) deallocate(Dimensions)
      if (allocated(MaxDimensions)) deallocate(MaxDimensions)
    endif
    call h5sclose_f(SpaceId,Error)
  else
    Iret = -1
  endif
  call h5dclose_f(DatasetId,Error)
  if (lhook) call dr_hook('GetDatasetDataspace',1,zhook_handle)
END SUBROUTINE GetDatasetDataspace

!**************************************************************************************************
! Subroutine/Function | ReadDataInt1, ReadDataInt2, ReadDataInt3
! Source Filename     | bator_decodhdf5_mod.F90
! author              | F. Guillaume (CNRM/GMAP/OBS)
! creation date       | 31/07/2016
!--------------------------------------------------------------------------------------------------
! Description : subroutines to read integer data (scalar, 1D, 2D or 3D).
!               WARNING - CAN RETURN REAL ARRAY if WaitedAtomic == NABSO
!--------------------------------------------------------------------------------------------------
! History :
!         - 04/2017 Guillaume F. : 1.0 release version
!**************************************************************************************************
SUBROUTINE ReadDataInt1(IdIn,Inte1Buf,Dimensions,MaxDimensions,WaitedAtomic,Iret)   ! scalar or 1D
  integer(kind=jpim),intent(inout)              :: Iret
  integer(HID_T),intent(in)                     :: IdIn
  integer(kind=jpim),intent(in)                 :: WaitedAtomic
  integer(kind=HSIZE_T),dimension(:),intent(in) :: Dimensions
  integer(kind=HSIZE_T),dimension(:),intent(in) :: MaxDimensions
  integer(kind=jpim)                            :: Error
  integer(kind=jpib),dimension(:),pointer       :: Inte1Buf
  real(kind=jprb)                               :: zhook_handle
  type(C_PTR)                                   :: FortPtr

  if (lhook) call dr_hook('ReadDataInt1',0,zhook_handle)
  allocate(Inte1Buf(Dimensions(1)),STAT=Error)
  if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate Inte1Buf(:)")
  FortPtr = C_LOC(Inte1Buf(1))
  call h5dread_f(IdIn,HdfIntKindType,FortPtr,Error)
  if (Error /= 0) Iret = 3
  if (WaitedAtomic == nabso) then
    allocate(Real1Buf(Dimensions(1)),STAT=Error)
    if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate Real1Buf(:)")
    Real1Buf = real(Inte1Buf,jprb)
    deallocate(Inte1Buf)
  endif
  if (lhook) call dr_hook('ReadDataInt1',1,zhook_handle)
END SUBROUTINE ReadDataInt1

SUBROUTINE ReadDataInt2(IdIn,Inte2Buf,Dimensions,MaxDimensions,WaitedAtomic,Iret)   ! 2D
  integer(kind=jpim),intent(inout)              :: Iret
  integer(HID_T),intent(in)                     :: IdIn
  integer(kind=jpim),intent(in)                 :: WaitedAtomic
  integer(kind=HSIZE_T),dimension(:),intent(in) :: Dimensions
  integer(kind=HSIZE_T),dimension(:),intent(in) :: MaxDimensions
  integer(kind=jpib),dimension(:,:),pointer     :: Inte2Buf
  integer(kind=jpim)                            :: Error
  real(kind=jprb)                               :: zhook_handle
  type(C_PTR)                                   :: FortPtr

  if (lhook) call dr_hook('ReadDataInt2',0,zhook_handle)
  allocate(Inte2Buf(Dimensions(1),Dimensions(2)),STAT=Error)
  if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate Inte2Buf(:)")
  FortPtr = C_LOC(Inte2Buf(1,1))
  call h5dread_f(IdIn,HdfIntKindType,FortPtr,Error)
  if (Error /= 0) Iret = 3
  if (WaitedAtomic == nabso) then
    allocate(Real2Buf(Dimensions(1),Dimensions(2)),STAT=Error)
    if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate Real2Buf(:)")
    Real2Buf = real(Inte2Buf,jprb)
    deallocate(Inte2Buf)
  endif
  if (lhook) call dr_hook('ReadDataInt2',1,zhook_handle)
END SUBROUTINE ReadDataInt2

SUBROUTINE ReadDataInt3(IdIn,Inte3Buf,Dimensions,MaxDimensions,WaitedAtomic,Iret)   ! 3D
  integer(kind=jpim),intent(inout)              :: Iret
  integer(HID_T),intent(in)                     :: IdIn
  integer(kind=jpim),intent(in)                 :: WaitedAtomic
  integer(kind=HSIZE_T),dimension(:),intent(in) :: Dimensions
  integer(kind=HSIZE_T),dimension(:),intent(in) :: MaxDimensions
  integer(kind=jpib),dimension(:,:,:),pointer   :: Inte3Buf
  integer(kind=jpim)                            :: Error
  real(kind=jprb)                               :: zhook_handle
  type(C_PTR)                                   :: FortPtr

  if (lhook) call dr_hook('ReadDataInt3',0,zhook_handle)
  allocate(Inte3Buf(Dimensions(1),Dimensions(2),Dimensions(3)),STAT=Error)
  if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate Inte3Buf(:)")
  FortPtr = C_LOC(Inte3Buf(1,1,1))
  call h5dread_f(IdIn,HdfIntKindType,FortPtr,Error)
  if (Error /= 0) Iret = 3
  if (WaitedAtomic == nabso) then
    allocate(Real3Buf(Dimensions(1),Dimensions(2),Dimensions(3)),STAT=Error)
    if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate Real3Buf(:)")
    Real3Buf = real(Inte3Buf,jprb)
    deallocate(Inte3Buf)
  endif
  if (lhook) call dr_hook('ReadDataInt3',1,zhook_handle)
END SUBROUTINE ReadDataInt3


!**************************************************************************************************
! Subroutine/Function | ReadDataReal1, ReadDataReal2, ReadDataReal3
! Source Filename     | bator_decodhdf5_mod.F90
! author              | F. Guillaume (CNRM/GMAP/OBS)
! creation date       | 31/07/2016
!--------------------------------------------------------------------------------------------------
! Description : subroutines to read real data (scalar, 1D, 2D or 3D).
!--------------------------------------------------------------------------------------------------
! History :
!         - 04/2017 Guillaume F. : 1.0 release version
!**************************************************************************************************
SUBROUTINE ReadDataReal1(IdIn,Real1Buf,Dimensions,MaxDimensions,Iret)   ! scalar or 1D
  integer(kind=jpim),intent(inout)              :: Iret
  integer(HID_T),intent(in)                     :: IdIn
  integer(kind=HSIZE_T),dimension(:),intent(in) :: Dimensions
  integer(kind=HSIZE_T),dimension(:),intent(in) :: MaxDimensions
  integer(kind=jpim)                            :: Error
  real(kind=jprb),dimension(:),pointer          :: Real1Buf
  real(kind=jprb)                               :: zhook_handle
  type(C_PTR)                                   :: FortPtr

  if (lhook) call dr_hook('ReadDataReal1',0,zhook_handle)
  allocate(Real1Buf(Dimensions(1)),STAT=Error)
  if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate Real1Buf(:)")
  FortPtr = C_LOC(Real1Buf(1))
  call h5dread_f(IdIn,HdfRealKindType,FortPtr,Error)
  if (Error /= 0) Iret = 3
  if (lhook) call dr_hook('ReadDataReal1',1,zhook_handle)
END SUBROUTINE ReadDataReal1

SUBROUTINE ReadDataReal2(IdIn,Real2Buf,Dimensions,MaxDimensions,Iret)   ! 2D
  integer(kind=jpim),intent(inout)              :: Iret
  integer(HID_T),intent(in)                     :: IdIn
  integer(kind=HSIZE_T),dimension(:),intent(in) :: Dimensions
  integer(kind=HSIZE_T),dimension(:),intent(in) :: MaxDimensions
  integer(kind=jpim)                            :: Error
  real(kind=jprb),dimension(:,:),pointer        :: Real2Buf
  real(kind=jprb)                               :: zhook_handle
  type(C_PTR)                                   :: FortPtr

  if (lhook) call dr_hook('ReadDataReal2',0,zhook_handle)
  allocate(Real2Buf(Dimensions(1),Dimensions(2)),STAT=Error)
  if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate Real2Buf(:)")
  FortPtr = C_LOC(Real2Buf(1,1))
  call h5dread_f(IdIn,HdfRealKindType,FortPtr,Error)
  if (Error /= 0) Iret = 3
  if (lhook) call dr_hook('ReadDataReal2',1,zhook_handle)
END SUBROUTINE ReadDataReal2

SUBROUTINE ReadDataReal3(IdIn,Real3Buf,Dimensions,MaxDimensions,Iret)   ! 3D
  integer(kind=jpim),intent(inout)              :: Iret
  integer(HID_T),intent(in)                     :: IdIn
  integer(kind=HSIZE_T),dimension(:),intent(in) :: Dimensions
  integer(kind=HSIZE_T),dimension(:),intent(in) :: MaxDimensions
  integer(kind=jpim)                            :: Error
  real(kind=jprb),dimension(:,:,:),pointer      :: Real3Buf
  real(kind=jprb)                               :: zhook_handle
  type(C_PTR)                                   :: FortPtr

  if (lhook) call dr_hook('ReadDataReal3',0,zhook_handle)
  allocate(Real3Buf(Dimensions(1),Dimensions(2),Dimensions(3)),STAT=Error)
  if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate Real3Buf(:)")
  FortPtr = C_LOC(Real3Buf(1,1,1))
  call h5dread_f(IdIn,HdfRealKindType,FortPtr,Error)
  if (Error /= 0) Iret = 3
  if (lhook) call dr_hook('ReadDataReal3',1,zhook_handle)
END SUBROUTINE ReadDataReal3

!**************************************************************************************************
! Subroutine/Function | GetAttribute
! Source Filename     | bator_decodhdf5_mod.F90
! author              | F. Guillaume (CNRM/GMAP/OBS)
! creation date       | 31/07/2016
!--------------------------------------------------------------------------------------------------
! Description : initialization to get an attribute.
!--------------------------------------------------------------------------------------------------
! History :
!         - 03/2017 Guillaume F. : 1.0 release version
!**************************************************************************************************
SUBROUTINE GetAttribute(IdIn,AttribName,WaitedRank,WaitedAtomic,Iret)
  integer(kind=HID_T),intent(in)                 :: IdIn
  character(len=*),intent(in)                    :: AttribName
  integer(kind=jpim),intent(in)                  :: WaitedRank
  integer(kind=jpim),intent(in)                  :: WaitedAtomic
  integer(kind=jpim),intent(out)                 :: Iret
  integer(kind=jpim)                             :: CurrentRank
  integer(kind=jpim)                             :: CurrentAtomic
  integer(kind=jpim)                             :: Error
  integer(kind=HID_T)                            :: AttribId
  integer(kind=HID_T)                            :: TypeId
  integer(kind=HID_T)                            :: SpaceId
  integer(kind=HSIZE_T),dimension(:),allocatable :: Dimensions
  integer(kind=HSIZE_T),dimension(:),allocatable :: MaxDimensions
  real(kind=jprb)                                :: zhook_handle
  logical                                        :: IsSimple

  if (lhook) call dr_hook('GetAttribute',0,zhook_handle)
  Iret = 0
  call h5aopen_f(IdIn,AttribName,AttribId,Error)
  if (Error == 0)  then                                 ! test d'existence de l'attribut
    call h5aget_space_f(AttribId,SpaceId,Error)
    call h5sis_simple_f(SpaceId,IsSimple,Error)
    if (Error /= 0 .OR. .NOT.IsSimple) then
      Iret = 1
    else
      call h5sget_simple_extent_ndims_f(SpaceId,CurrentRank,Error)
      call h5aget_type_f(AttribId,TypeId,Error)
      call h5tget_class_f(TypeId,CurrentAtomic,Error)
      if (WaitedRank /= CurrentRank .or. WaitedAtomic /= CurrentAtomic) then
        Iret = 2
      else
        if (CurrentRank == 0) CurrentRank = 1       ! a scalar value <=> 1D array values
        allocate(Dimensions(CurrentRank),MaxDimensions(CurrentRank),STAT=Error)
        if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate Dimensions(:) or MaxDimensions(:)")
        select case (WaitedRank)
        case (0)
          dimensions(CurrentRank)    = 1
          maxdimensions(CurrentRank) = 1
          if (CurrentAtomic == H5T_STRING_F)  call ReadAttribute(AttribId,TypeId,String1Buf,Dimensions,MaxDimensions,Iret)
          if (CurrentAtomic == H5T_INTEGER_F) call ReadAttribute(AttribId,Inte1Buf,Dimensions,MaxDimensions,Iret)
          if (CurrentAtomic == H5T_FLOAT_F)   call ReadAttribute(AttribId,Real1Buf,Dimensions,MaxDimensions,Iret)
        case (1)
          call h5sget_simple_extent_dims_f(SpaceId,Dimensions,MaxDimensions,Error)
          if (CurrentAtomic == H5T_STRING_F)  call ReadAttribute(AttribId,TypeId,String1Buf,Dimensions,MaxDimensions,Iret)
          if (CurrentAtomic == H5T_INTEGER_F) call ReadAttribute(AttribId,Inte1Buf,Dimensions,MaxDimensions,Iret)
          if (CurrentAtomic == H5T_FLOAT_F)   call ReadAttribute(AttribId,Real1Buf,Dimensions,MaxDimensions,Iret)
        case (2)
          call h5sget_simple_extent_dims_f(SpaceId,Dimensions,MaxDimensions,Error)
          if (CurrentAtomic == H5T_STRING_F)  call ReadAttribute(AttribId,TypeId,String2Buf,Dimensions,MaxDimensions,Iret)
          if (CurrentAtomic == H5T_INTEGER_F) call ReadAttribute(AttribId,Inte2Buf,Dimensions,MaxDimensions,Iret)
          if (CurrentAtomic == H5T_FLOAT_F)   call ReadAttribute(AttribId,Real2Buf,Dimensions,MaxDimensions,Iret)
        case default
          Iret = 4
        end select
        if (allocated(Dimensions)) deallocate(Dimensions)
        if (allocated(MaxDimensions)) deallocate(MaxDimensions)
      endif
      call h5tclose_f(TypeId,Error)
    endif
    call h5sclose_f(SpaceId,Error)
    call h5aclose_f(AttribId,Error)
  else
    Iret = -1
  endif
  if (lhook) call dr_hook('GetAttribute',1,zhook_handle)
END SUBROUTINE GetAttribute

!**************************************************************************************************
! Subroutine/Function | ReadAttributeStr1, ReadAttributeStr2
! Source Filename     | bator_decodhdf5_mod.F90
! author              | F. Guillaume (CNRM/GMAP/OBS)
! creation date       | 31/07/2016
!--------------------------------------------------------------------------------------------------
! Description : subroutines to read string attributes (scalar, 1D or 2D)
!--------------------------------------------------------------------------------------------------
! History :
!         - 03/2017 Guillaume F. : 1.0 release version
!**************************************************************************************************
SUBROUTINE ReadAttributeStr1(IdIn,TypeId,String1Buf,Dimensions,MaxDimensions,Iret)   ! scalar or 1D
  integer(kind=jpim),intent(inout)                :: Iret
  integer(HID_T),intent(in)                       :: IdIn
  integer(HID_T),intent(in)                       :: TypeId
  integer(kind=HSIZE_T),dimension(:),intent(in)   :: Dimensions
  integer(kind=HSIZE_T),dimension(:),intent(in)   :: MaxDimensions
  character(len=LabelLength),dimension(:),pointer :: String1Buf
  integer(SIZE_T)                                 :: Longueur
  integer(kind=jpim)                              :: TypeString
  integer(kind=jpim)                              :: Error
  real(kind=jprb)                                 :: zhook_handle

  if (lhook) call dr_hook('ReadAttributeStr1',0,zhook_handle)
  allocate(String1Buf(Dimensions(1)),STAT=Error)
  if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate String1Buf(:)")
  call h5tget_size_f(TypeId,Longueur,Error)
  call h5tget_strpad_f(TypeId,TypeString,Error)
  if (TypeString == 0) Longueur = Longueur - 1
  call h5aread_f(IdIn,TypeId,String1Buf,Dimensions,Error)
  if (Error /= 0 .or. Longueur > LabelLength) then
    Iret = 3
  else
    if (Longueur < LabelLength) String1Buf(:)(Longueur+1:LabelLength) = repeat(" ",(LabelLength-Longueur))
  endif
  if (lhook) call dr_hook('ReadAttributeStr1',1,zhook_handle)
END SUBROUTINE ReadAttributeStr1


SUBROUTINE ReadAttributeStr2(IdIn,TypeId,String2Buf,Dimensions,MaxDimensions,Iret)   ! 2D
  integer(kind=jpim),intent(inout)                  :: Iret
  integer(HID_T),intent(in)                         :: IdIn
  integer(SIZE_T)                                   :: Longueur
  integer(HID_T),intent(in)                         :: TypeId
  integer(kind=HSIZE_T),dimension(:),intent(in)     :: Dimensions
  integer(kind=HSIZE_T),dimension(:),intent(in)     :: MaxDimensions
  character(len=LabelLength),dimension(:,:),pointer :: String2Buf
  integer(kind=jpim)                                :: TypeString
  integer(kind=jpim)                                :: Error
  real(kind=jprb)                                   :: zhook_handle

  if (lhook) call dr_hook('ReadAttributeStr2',0,zhook_handle)
  allocate(String2Buf(Dimensions(1),Dimensions(2)),STAT=Error)
  if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate String2Buf(:)")
  call h5tget_size_f(TypeId,Longueur,Error)
  call h5tget_strpad_f(TypeId,TypeString,Error)
  if (TypeString == 0) Longueur = Longueur - 1
  call h5aread_f(IdIn,TypeId,String2Buf,Dimensions,Error)
  if (Error /= 0 .or. Longueur > LabelLength) then
    Iret = 3
  else
    if (Longueur < LabelLength)String2Buf(:,:)(Longueur+1:LabelLength) = repeat(" ",(LabelLength-Longueur))
  endif
  if (lhook) call dr_hook('ReadAttributeStr2',1,zhook_handle)
END SUBROUTINE ReadAttributeStr2

!**************************************************************************************************
! Subroutine/Function | ReadAttributeInt1, ReadAttributeInt2
! Source Filename     | bator_decodhdf5_mod.F90
! author              | F. Guillaume (CNRM/GMAP/OBS)
! creation date       | 31/07/2016
!--------------------------------------------------------------------------------------------------
! Description : subroutines to read integer attributes (scalar, 1D or 2D)
!--------------------------------------------------------------------------------------------------
! History :
!         - 03/2017 Guillaume F. : 1.0 release version
!**************************************************************************************************
SUBROUTINE ReadAttributeInt1(IdIn,Inte1Buf,Dimensions,MaxDimensions,Iret)   ! scalar or 1D
  integer(kind=jpim),intent(inout)              :: Iret
  integer(HID_T),intent(in)                     :: IdIn
  integer(kind=jpib),dimension(:),pointer       :: Inte1Buf
  integer(kind=HSIZE_T),dimension(:),intent(in) :: Dimensions
  integer(kind=HSIZE_T),dimension(:),intent(in) :: MaxDimensions
  integer(kind=jpim)                            :: Error
  real(kind=jprb)                               :: zhook_handle
  type(C_PTR)                                   :: FortPtr

  if (lhook) call dr_hook('ReadAttributeInt1',0,zhook_handle)
  allocate(Inte1Buf(Dimensions(1)),STAT=Error)
  if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate Inte1Buf(:)")
  FortPtr = C_LOC(Inte1Buf(1))
  call h5aread_f(IdIn,HdfIntKindType,FortPtr,Error)
  if (Error /= 0) Iret = 3
  if (lhook) call dr_hook('ReadAttributeInt1',1,zhook_handle)
END SUBROUTINE ReadAttributeInt1


SUBROUTINE ReadAttributeInt2(IdIn,Inte2Buf,Dimensions,MaxDimensions,Iret)   ! 2D
  integer(kind=jpim),intent(inout)              :: Iret
  integer(HID_T),intent(in)                     :: IdIn
  integer(kind=jpib),dimension(:,:),pointer     :: Inte2Buf
  integer(kind=HSIZE_T),dimension(:),intent(in) :: Dimensions
  integer(kind=HSIZE_T),dimension(:),intent(in) :: MaxDimensions
  integer(kind=jpim)                            :: Error
  real(kind=jprb)                               :: zhook_handle
  type(C_PTR)                                   :: FortPtr

  if (lhook) call dr_hook('ReadAttributeInt2',0,zhook_handle)
  allocate(Inte2Buf(Dimensions(1),Dimensions(2)),STAT=Error)
  if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate Inte2Buf(:)")
  FortPtr = C_LOC(Inte2Buf(1,1))
  call h5aread_f(IdIn,HdfIntKindType,FortPtr,Error)
  if (Error /= 0) Iret = 3
  if (lhook) call dr_hook('ReadAttributeInt2',1,zhook_handle)
END SUBROUTINE ReadAttributeInt2

!**************************************************************************************************
! Subroutine/Function | ReadAttributeReal1, ReadAttributeReal2
! Source Filename     | bator_decodhdf5_mod.F90
! author              | F. Guillaume (CNRM/GMAP/OBS)
! creation date       | 31/07/2016
!--------------------------------------------------------------------------------------------------
! Description : subroutines to read float attributes (scalar, 1D or 2D)
!--------------------------------------------------------------------------------------------------
! History :
!         - 03/2017 Guillaume F. : 1.0 release version
!**************************************************************************************************
SUBROUTINE ReadAttributeReal1(IdIn,Real1Buf,Dimensions,MaxDimensions,Iret)   ! scalar or 1D
  integer(kind=jpim),intent(inout)              :: Iret
  integer(HID_T),intent(in)                     :: IdIn
  integer(kind=HSIZE_T),dimension(:),intent(in) :: Dimensions
  integer(kind=HSIZE_T),dimension(:),intent(in) :: MaxDimensions
  integer(kind=jpim)                            :: Error
  real(kind=jprb),dimension(:),pointer          :: Real1Buf
  real(kind=jprb)                               :: zhook_handle
  type(C_PTR)                                   :: FortPtr

  if (lhook) call dr_hook('ReadAttributeReal1',0,zhook_handle)
  allocate(Real1Buf(Dimensions(1)),STAT=Error)
  if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate Real1Buf(:)")
  FortPtr = C_LOC(Real1Buf(1))
  call h5aread_f(IdIn,HdfRealKindType,FortPtr,Error)
  if (Error /= 0) Iret = 3
  if (lhook) call dr_hook('ReadAttributeReal1',1,zhook_handle)
END SUBROUTINE ReadAttributeReal1


SUBROUTINE ReadAttributeReal2(IdIn,Real2Buf,Dimensions,MaxDimensions,Iret)   ! 2D
  integer(kind=jpim),intent(inout)              :: Iret
  integer(HID_T),intent(in)                     :: IdIn
  integer(kind=HSIZE_T),dimension(:),intent(in) :: Dimensions
  integer(kind=HSIZE_T),dimension(:),intent(in) :: MaxDimensions
  integer(kind=jpim)                            :: Error
  real(kind=jprb),dimension(:,:),pointer        :: Real2Buf
  real(kind=jprb)                               :: zhook_handle
  type(C_PTR)                                   :: FortPtr

  if (lhook) call dr_hook('ReadAttributeReal2',0,zhook_handle)
  allocate(Real2Buf(Dimensions(1),Dimensions(2)),STAT=Error)
  if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate Real2Buf(:)")
  FortPtr = C_LOC(Real2Buf(1,1))
  call h5aread_f(IdIn,HdfRealKindType,FortPtr,Error)
  if (Error /= 0) Iret = 3
  if (lhook) call dr_hook('ReadAttributeReal2',1,zhook_handle)
END SUBROUTINE ReadAttributeReal2

!**************************************************************************************************
! Subroutine/Function | Hdf5Type
! Source Filename     | bator_decodhdf5_mod.F90
! author              | F. Guillaume (CNRM/GMAP/OBS)
! creation date       | 24/03/2017
!--------------------------------------------------------------------------------------------------
! Description : function to define which kind of Atomic is needed to read an attribute or data.
!--------------------------------------------------------------------------------------------------
! History :
!         - 03/2017 Guillaume F. : 1.0 release version
!**************************************************************************************************
function Hdf5Type(WaitedAtomic)
  integer(kind=jpim)            :: Hdf5Type
  integer(kind=jpim),intent(in) :: WaitedAtomic
  real(kind=jprb)               :: zhook_handle

  if (lhook) call dr_hook('Hdf5Type',0,zhook_handle)
  select case (WaitedAtomic)
  case (0)
    Hdf5Type = H5T_STRING_F
  case (1)
    Hdf5Type = H5T_INTEGER_F
  case (2)
    Hdf5Type = H5T_FLOAT_F
  case (-1)
    Hdf5Type = nabso
  case default
    Hdf5Type = nabsi
  end select
  if (lhook) call dr_hook('Hdf5Type',1,zhook_handle)
end function Hdf5Type

!**************************************************************************************************
! Subroutine/Function | Mtvza
! Source Filename     | bator_decodhdf5_mod.F90
! Author              | F. Guillaume (CNRM/GMAP/OBS)
! Creation Date       | 12/12/2016
!--------------------------------------------------------------------------------------------------
! Description : decod meteor-m2 data.
!--------------------------------------------------------------------------------------------------
! Global Variables : 
!--------------------------------------------------------------------------------------------------
! History :
!         - 09/2017 Guillaume F. : 1.0 release version
!**************************************************************************************************
SUBROUTINE Mtvza(FileId,TabSlots,Kobs,Kw,iterr,Conformity)
  TYPE(TDate),dimension(:),intent(in)            :: TabSlots
  integer(kind=jpim),intent(inout)               :: Kobs
  integer(kind=jpib),intent(inout)               :: Kw
  integer(HID_T),intent(in)                      :: FileId
  integer(kind=jpim),dimension(50),intent(inout) :: iterr
  logical,intent(out)                            :: Conformity
  integer(kind=jpim)                             :: Kobs0
  integer(kind=jpib)                             :: kw0
  integer(HID_T)                                 :: DsetId
  integer(kind=jpim)                             :: WaitedAtomic, WaitedRank, FirstDay, Scanline, Fov, RootLen
  integer(kind=jpim)                             :: Error, Iret
  integer(kind=jpim),DIMENSION(6)                :: DateRef
  integer(kind=jpim)                             :: NMembers,NDataset,ObjectType,CurrentChannel,i,j,k
  character(len=LabelLength)                     :: NomMembre
  integer(kind=jpim)                             :: SelectedSat = -1
  integer(kind=jpib),dimension(:,:),allocatable  :: TimeOfDay
  integer(kind=jpim),dimension(:,:),allocatable  :: Surface
  real(kind=jprb),dimension(:,:),allocatable     :: Lat
  real(kind=jprb),dimension(:,:),allocatable     :: Lon
  real(kind=jprb),dimension(:,:),allocatable     :: SunAzimuth
  real(kind=jprb),dimension(:,:),allocatable     :: SunZenith
  real(kind=jprb)                                :: TbMin,TbMax
  real(kind=jprb)                                :: zhook_handle
  type(TDate)                                    :: BeginDate
  type canal
    real(kind=jprb),dimension(:,:),pointer       :: Values      => NULL()
  end type canal
  type(canal),dimension(NbMaxMtvzaChan)          :: Channels

  if (lhook) call dr_hook('Mtvza',0,zhook_handle)

  Conformity = .TRUE.
  call h5gn_members_f(FileId,"/",NMembers,Error)
  if (Error /= 0) Conformity = .FALSE.


  i = lbound(HMTVZA,1)
  do while (i <= ubound(HMTVZA,1))
    j       = 0
    RootLen = len(trim(HMTVZA(i)%DatasetNameRoot))
    do while (j <= NMembers-1)
      call h5gget_obj_info_idx_f(FileId,"/",j,NomMembre,ObjectType,Error)
      if (ObjectType == H5G_DATASET_F) then
        if (NomMembre(1:RootLen) == trim(HMTVZA(i)%DatasetNameRoot)) then
          j           = NMembers
          SelectedSat = i
          Conformity  = .TRUE.
        else
          Conformity  = .FALSE.
          j           = j + 1
        endif
      else
        Conformity = .FALSE.
        j          = j + 1
      endif
    enddo
    i                 = i + 1
    if (Conformity) i = ubound(HMTVZA,1) + 1
  enddo

  if (Conformity) then
    do i=0, NMembers-1
      call h5gget_obj_info_idx_f(FileId,"/",i,NomMembre,ObjectType,Error)
      if (ObjectType == H5G_DATASET_F) then
        WaitedRank       = 2
        WaitedAtomic     = Hdf5Type(-1)
        call GetData(FileId,NomMembre,WaitedRank,WaitedAtomic,Iret)
        if (Iret == 0) then
          if (trim(HMTVZA(SelectedSat)%Time) == NomMembre(1:len(trim(HMTVZA(SelectedSat)%Time)))) then
            allocate(TimeOfDay(size(real2buf,1),size(real2buf,2)),STAT=Error)
            if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate TimeOfDay(:,:)")
            TimeOfDay = nint(real2buf/1000_jprb,jpib)
          elseif (trim(HMTVZA(SelectedSat)%Julien) == NomMembre(1:len(trim(HMTVZA(SelectedSat)%Julien)))) then
            DateRef(:)   = 0
            DateRef(2:3) = 1
            DateRef(1)   = 2000 + (nint(real2buf(1,1),jpim)/1000_jpim)
            BeginDate = AddDate(NewDate(DateRef),(mod(nint(real2buf(1,1),jpib),1000)-1)*86400)
            DateRef = ReturnDateArray(BeginDate)
          elseif (trim(HMTVZA(SelectedSat)%Lat) == NomMembre(1:len(trim(HMTVZA(SelectedSat)%Lat)))) then
            allocate(Lat(size(real2buf,1),size(real2buf,2)),STAT=Error)
            if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate Lat(:,:)")
            Lat      = real2buf
            Scanline = size(real2buf,2)
            Fov      = size(real2buf,1)
            where (abs(Lat) > 90._jprb) Lat(:,:) = rabsi
          elseif (trim(HMTVZA(SelectedSat)%Lon) == NomMembre(1:len(trim(HMTVZA(SelectedSat)%Lon)))) then
            allocate(Lon(size(real2buf,1),size(real2buf,2)),STAT=Error)
            if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate Lon(:,:)")
            Lon = real2buf
            where (abs(Lon) > 180._jprb) Lon(:,:) = rabsi
          elseif (trim(HMTVZA(SelectedSat)%Surf) == NomMembre(1:len(trim(HMTVZA(SelectedSat)%Surf)))) then
            allocate(Surface(size(real2buf,1),size(real2buf,2)),STAT=Error)
            if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate Surface(:,:)")
            Surface = nint(real2buf,jpim)
          elseif (trim(HMTVZA(SelectedSat)%SunAzimuth) == NomMembre(1:len(trim(HMTVZA(SelectedSat)%SunAzimuth)))) then
            allocate(SunAzimuth(size(real2buf,1),size(real2buf,2)),STAT=Error)
            if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate SunAzimuth(:,:)")
            SunAzimuth = real2buf
          elseif (trim(HMTVZA(SelectedSat)%SunZenith) == NomMembre(1:len(trim(HMTVZA(SelectedSat)%SunZenith)))) then
            allocate(SunZenith(size(real2buf,1),size(real2buf,2)),STAT=Error)
            if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate SunZenith(:,:)")
            SunZenith = real2buf
          else
            do j=1,HMTVZA(SelectedSat)%NbChannels
              if (trim(HMTVZA(SelectedSat)%NamChannels(j)) == NomMembre(1:len(trim(HMTVZA(SelectedSat)%NamChannels(j))))) then
                CurrentChannel = HMTVZA(SelectedSat)%Channels(j)
                allocate(Channels(CurrentChannel)%Values(size(real2buf,1),size(real2buf,2)),STAT=Error)
                if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate Channels(SelectedSat)%Values(:,:)")
                Channels(CurrentChannel)%Values = real2buf/100._jprb
                call h5dopen_f(FileId,NomMembre,DsetId,Error)
                WaitedRank   = 0
                WaitedAtomic = Hdf5Type(1)
                call GetAttribute(DsetId,HMTVZA(SelectedSat)%TbMaxAttrib,WaitedRank,WaitedAtomic,Iret)
                if (Iret == 0) then
                  TbMax        = Inte1Buf(1)
                  deallocate(inte1buf)
                  where (Channels(CurrentChannel)%Values > TbMax) Channels(CurrentChannel)%Values = rabsi
                endif
                WaitedRank   = 0
                WaitedAtomic = Hdf5Type(1)
                call GetAttribute(DsetId,HMTVZA(SelectedSat)%TbMinAttrib,WaitedRank,WaitedAtomic,Iret)
                if (Iret == 0) then
                  TbMin        = Inte1Buf(1)
                  deallocate(inte1buf)
                  where (Channels(CurrentChannel)%Values < TbMin) Channels(CurrentChannel)%Values = rabsi
                endif
                call h5dclose_f(DsetId,Error)
              endif
            enddo
          endif
          deallocate(real2buf)
        endif
      endif
    enddo
  
    kobs0 = kobs
    kw0   = kw

    do i=1, scanline
      DateRef = ReturnDateArray(AddDate(BeginDate,TimeOfDay(1,i)))
      if (.not.(VerifDate(DateRef,TabSlots))) then
        iterr(2) = iterr(2) + 1
        cycle
      endif
      do j=1, fov
        if (Lat(j,i) == rabsi .or. Lon(j,i) == rabsi) then 
          iterr(1) = iterr(1) + 1
          cycle
        endif

        kobs = kobs +  1

        zent(kobs,NCIOTP) = NSATEM
        zent(kobs,NCIRFL) = 11111
        zent(kobs,NCIOCH) = NGTHRB
        write(clsid(kobs),'(5x,i3.3)')  SelectedSat                    ! identite satellite
        zent(kobs,NCIDAT) = DateRef(1)*10000+DateRef(2)*100+DateRef(3)
        zent(kobs,NCIETM) = DateRef(4)*10000+DateRef(5)*100+DateRef(6)
        zent(kobs,NCISTD) = 0                                          ! numero d'amendement
        zent(kobs,NCINLV) = NbMaxMtvzaChan                             ! nombre de wagon
        zent(kobs,NCILAT) = Lat(j,i)
        zent(kobs,NCILON) = Lon(j,i)
        zentsup(kobs,1)   = HMTVZA(SelectedSat)%sensor                 ! type de capteur
        zentsup(kobs,2)   = i                                          ! scan line
        zentsup(kobs,5)   = 5                                          ! signification verticale
        zentsup(kobs,3)   = j                                          ! fov
        zentsup(kobs,4)   = surface(j,i)                               ! flag de surface
        zentsup(kobs,9)   = SunZenith(j,i)
        zentsup(kobs,10)  = SunAzimuth(j,i)

        do k=1,NbMaxMtvzaChan
          zwagon(kw+k,1) = Varno%rawbt
          zwagon(kw+k,5) = 2063
          zwagon(kw+k,2) = k
          if (associated(Channels(k)%Values)) then
            zwagon(kw+k,4) = Channels(k)%Values(j,i)
            deallocate(Channels(k)%Values)
          endif
        enddo
        kw = kw + NbMaxMtvzaChan
      enddo
    enddo

    deallocate(Lat)
    deallocate(Lon)
    deallocate(TimeOfDay)
    deallocate(SunAzimuth)
    deallocate(SunZenith)
    deallocate(Surface)
    write(batout,'(6x,"Selected Obs = ",I9,"  --> ",I12," datas.")')    Kobs-Kobs0, Kw-Kw0
    write(batout,'("Total selected Obs = ",I9,"  --> ",I12," datas.")') Kobs,       Kw
  endif

  if (lhook) call dr_hook('Mtvza',1,zhook_handle)
END SUBROUTINE Mtvza

!**************************************************************************************************
! Subroutine/Function | odim
! Source Filename     | bator_decodhdf5_mod.F90
! Author              | F. Guillaume (CNRM/GMAP/OBS)
! Creation Date       | 09/08/2016
! Modification        | 12/04/2021  
!                     | Idir DEHMOUS  : READ LOCAL RMI HDF5 RADAR FILES 
!--------------------------------------------------------------------------------------------------
! Description : decod odim radar data.
!--------------------------------------------------------------------------------------------------
! Global Variables : 
!--------------------------------------------------------------------------------------------------
! History : prototype - en cours de reecriture
!**************************************************************************************************
SUBROUTINE Odim(FileId,TabSlots,AnalysisDate,Kobs,Kw,Kfic,Conformity)
  TYPE(TDate),intent(in)                         :: AnalysisDate
  integer(HID_T),intent(in)                      :: FileId
  TYPE(TDate),dimension(:),intent(in)            :: TabSlots
  integer(kind=jpim),intent(in)                  :: Kfic
  integer(kind=jpim),intent(inout)               :: Kobs
  integer(kind=jpib),intent(inout)               :: Kw
  logical,intent(out)                            :: Conformity
  integer(HID_T)                                 :: GroupId
  integer(kind=jpim)                             :: Error, Iret, NumGDataset, NumGData, NumGQuality
  integer(kind=jpim)                             :: WaitedAtomic, WaitedRank
  integer(kind=jpim)                             :: NbConcernedData, CurrentMinElangle
  integer(kind=jpim)                             :: NMembers, NMembers1
  integer(kind=jpim)                             :: NDataset, NData, NQuality, ObjectType, NbNrays
  integer(kind=jpim)                             :: i, j, k, l
  integer(kind=jpib)                             :: DeltaTime
  real(kind=jprb)                                :: MaxRayLength, RayLength
  character(len=LabelLength)                     :: NomMembre, NomMembre1
  logical                                        :: LRMID 

  integer(kind=jpim)                             :: ispindir, kobs0, ilw, nb_obs, inlv, iobs
  integer(kind=jpim)                            ::  ifreq, icol, irow, iobs1
  integer(kind=jpim)                             :: i0, i1
  integer(kind=jpib)                             :: zaz_offset, nbw, iw, kw0
  real(kind=jpib)                                :: zconst,zsensib
  real(kind=jprb)                                :: zdimhpiy, RAE, zlat, zlon, zalt
  real(kind=jprb)                                :: zalpha, zdx, zdy, zphi
  real(kind=jprb)                                :: zthreshold
  real(kind=jprb),dimension(:,:),allocatable     :: ztent, ztwag
  real(kind=jprb),dimension(:),allocatable       :: zazim, zdist
  logical,dimension(:,:),allocatable             :: zthmask
  logical                                        :: empty  ! .FALSE. = doppler radar --> wind available
 

  TYPE GAttributes
    character(len=6)                             :: StartTime = '?'
    character(len=8)                             :: StartDate = '?'
    character(len=LabelLength)                   :: Quantity  = '?'
    real(kind=jprb)                              :: Gain      = rabsi
    real(kind=jprb)                              :: Offset    = rabsi
    real(kind=jprb)                              :: NoData    = rabsi
    real(kind=jprb)                              :: NoDetect  = rabsi
    real(kind=jprb)                              :: Elangle   = rabsi
    integer(kind=jpim)                           :: NBins     = nabsi
    integer(kind=jpim)                           :: NRays     = nabsi
    real(kind=jprb)                              :: RScale    = rabsi
    real(kind=jprb)                              :: Rstart    = rabsi
    character(len=LabelLength)                   :: Task      = '?'
  END TYPE GAttributes
  TYPE GrpData
    character(len=LabelLength*2+1)               :: Label
    TYPE(GAttributes)                            :: Attrib
    real(kind=jprb),dimension(:,:),pointer       :: Values => NULL()
  END TYPE GrpData
  TYPE GrpQual
    character(len=LabelLength*2+1)               :: Label
    TYPE(GAttributes)                            :: Attrib
    real(kind=jprb),dimension(:,:),pointer       :: Values => NULL()
  END TYPE GrpQual
  TYPE GrpDataset
    character(len=LabelLength)                   :: Label
    TYPE(GrpData),dimension(:),pointer           :: GData    => NULL()
    TYPE(GrpQual),dimension(:),pointer           :: GQuality => NULL()
    TYPE(GAttributes)                            :: Attrib
  END TYPE GrpDataset
  TYPE(GrpDataset),dimension(:),allocatable      :: FullDatasetList

  TYPE FilteredElangle
    real(kind=jprb)                              :: Elangle = rabsi
    TYPE(GrpData),pointer                        :: DBZH => NULL()
    integer(kind=jpib)                           :: DbZHDeltaTime = 86400
    TYPE(GrpData),pointer                        :: TH   => NULL()
    integer(kind=jpib)                           :: THDeltaTime = 86400
    TYPE(GrpData),pointer                        :: VRAD => NULL()
    integer(kind=jpib)                           :: VRadDeltaTime = 86400
    TYPE(GrpQual),pointer                        :: FLAG => NULL()
  END TYPE FilteredElangle
  TYPE(FilteredElangle),dimension(:),allocatable :: SelectedElangles

  TYPE Fichier
    character(len=5)                             :: Nod     = '?'
    integer(kind=jpim),dimension(6)              :: DateOpt = nabsi
    real(kind=jprb)                              :: Lat     = rabsi
    real(kind=jprb)                              :: Lon     = rabsi
    real(kind=jprb)                              :: Alt     = rabsi
    real(kind=jprb)                              :: BeamWidth = 1._jprb
    integer(kind=jpim)                           :: NPoints = nabsi
    integer(kind=jpim)                           :: NRayons = nabsi
    TYPE(GAttributes)                            :: Attrib
    TYPE(FilteredElangle),dimension(:),pointer   :: FinalElev
  END TYPE Fichier
  TYPE(Fichier)                                  :: Radar

  integer(kind=jpim),dimension(2)                :: SelectedNrays
  integer(kind=jpim),dimension(:,:),allocatable  :: NRaysPopulation
  integer(kind=jpim)                             :: NbSelectedElangles

  real(kind=jprb)                                :: zhook_handle

  if (lhook) call dr_hook('Odim',0,zhook_handle)

  Conformity = .TRUE.

! get top-level what group needed attributes
  call h5gopen_f(FileId,HODIM%GrpWhatName,GroupId,Error)
  if (Error /= 0) then
    Conformity = .FALSE.
  else
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(0)
    call GetAttribute(GroupId,HODIM%ObjectName,WaitedRank,WaitedAtomic,Iret)
    if (Iret /= 0) then
      Conformity = .FALSE.
    else
      if (trim(String1Buf(1)) == "PVOL" .or. trim(String1Buf(1)) == "SCAN" ) then
        write (Batout,'("*** INFO - BATOR : Type produit",9x,": ",A)') trim(String1Buf(1))
      else
        Conformity = .FALSE.
      endif
      deallocate(String1Buf)
    endif
  endif
  if (Conformity) then
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(0)
    call GetAttribute(GroupId,HODIM%SourceName,WaitedRank,WaitedAtomic,Iret)
    if (Iret /= 0) then
      Conformity = .FALSE.
    else
      if (index(String1Buf(1),'NOD') >0) then
        Radar%Nod = String1Buf(1)(index(String1Buf(1),'NOD')+4:index(String1Buf(1),'NOD')+8)
        if (any(HODIM%NodeNames == Radar%Nod)) then
          write (Batout,'("*** INFO - BATOR : Source NOD",11x,": ",A)') Radar%Nod
        else
          Conformity = .FALSE.
        endif
      else
        Conformity = .FALSE.
      endif
      deallocate(String1Buf)
    endif
  endif
  if (Conformity) then
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(0)
    call GetAttribute(GroupId,HODIM%DateName,WaitedRank,WaitedAtomic,Iret)
    if (Iret /= 0 ) then
      Conformity = .FALSE.
    else
      read (String1Buf(1),'(I4,I2,I2)') Radar%DateOpt(1),Radar%DateOpt(2),Radar%DateOpt(3)
      write (Batout,'("*** INFO - BATOR : Date optimale",8x,": ",I4.4,"-",I2.2,"-",I2.2)') &
           & Radar%DateOpt(1),Radar%DateOpt(2),Radar%DateOpt(3)
      deallocate(String1Buf)
      WaitedRank   = 0
      WaitedAtomic = Hdf5Type(0)
      call GetAttribute(GroupId,HODIM%TimeName,WaitedRank,WaitedAtomic,Iret)
      if (Iret /= 0) then
        Conformity = .FALSE.
      else
        read (String1Buf(1),'(I2,I2,I2)') Radar%DateOpt(4),Radar%DateOpt(5),Radar%DateOpt(6)
        write (Batout,'("*** INFO - BATOR : heure optimale",7x,": ",I2.2,":",I2.2,":",I2.2)') &
             & Radar%DateOpt(4),Radar%DateOpt(5),Radar%DateOpt(6)
        deallocate(String1Buf)
      endif
    endif
    if (.not.(VerifDate(Radar%DateOpt,TabSlots))) Conformity = .FALSE.
  endif
  call h5gclose_f(GroupId,Error)

! get top-level where group needed attributes
  if (Conformity) then
    call h5gopen_f(FileId,HODIM%GrpWhereName,GroupId,Error)
    if (Error /= 0) then
      Conformity = .FALSE.
    else
      WaitedRank   = 0
      WaitedAtomic = Hdf5Type(2)
      call GetAttribute(GroupId,HODIM%SiteHeightName,WaitedRank,WaitedAtomic,Iret)
      if (Iret /= 0) then
        Conformity = .FALSE.
      else
        Radar%Alt = real1buf(1)
        write (Batout,'("*** INFO - BATOR : Hauteur",14x,": ",F7.2)') Radar%Alt
        deallocate(real1buf)
      endif
    endif
  endif
  if (Conformity) then
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(2)
    call GetAttribute(GroupId,HODIM%SiteLatName,WaitedRank,WaitedAtomic,Iret)
    if (Iret /= 0) then
      Conformity = .FALSE.
    else
      Radar%Lat = real1buf(1)
      write (Batout,'("*** INFO - BATOR : lat",18x,": ",F11.6)') Radar%Lat
      deallocate(real1buf)
    endif
  endif
  if (Conformity) then
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(2)
    call GetAttribute(GroupId,HODIM%SiteLonName,WaitedRank,WaitedAtomic,Iret)
    if (Iret /= 0) then
      Conformity = .FALSE.
    else
      Radar%Lon = real1buf(1)
      write (Batout,'("*** INFO - BATOR : lon",18x,": ",F11.6)') Radar%Lon
      deallocate(real1buf)
    endif
  endif
  call h5gclose_f(GroupId,Error)

! get top-level how group needed attributes
  if (Conformity) then
    call h5gopen_f(FileId,HODIM%GrpHowName,GroupId,Error)
    if (Error /= 0) then
      Conformity = .FALSE.
    else
      WaitedRank   = 0
      WaitedAtomic = Hdf5Type(2)
      call GetAttribute(GroupId,HODIM%BeamWidthName,WaitedRank,WaitedAtomic,Iret)
      if (Iret /= 0) then 
        Conformity = .FALSE.
      else
        Radar%BeamWidth = Real1Buf(1)
        write (Batout,'("*** INFO - BATOR : BeamWidth",12x,": ",F7.2)') Radar%BeamWidth
        deallocate(Real1Buf)
      endif
    endif
  endif
  call h5gclose_f(GroupId,Error)

! get top-level standard attributes (same list as in dataset & data groups)
  if (Conformity) then
    call GetDAttributes(FileId,' ',Radar%Attrib)
  endif

! check whether there is at least one dataset (group)
  if (Conformity) then
    NDataset = 0
    call h5gn_members_f(FileId,"/",NMembers,Error)
    if (Error == 0) then
      do i=0, NMembers-1
        call h5gget_obj_info_idx_f(FileId,"/",i,NomMembre,ObjectType,Error)
        if (Error /= 0) then
          Conformity = .FALSE.
        else
          if (ObjectType == H5G_GROUP_F .and. NomMembre(1:len(trim(HODIM%GrpElevName))) == trim(HODIM%GrpElevName)) then
            NDataset = NDataset + 1
          endif
        endif
      enddo
    endif
    if (NDataset == 0) Conformity = .FALSE.
    write (Batout,'("*** INFO - BATOR : ",I3,"Dataset groups found.")') NDataset
  endif
if (Conformity) then
do i=0, NMembers - 1
call h5gget_obj_info_idx_f(FileId,"/",i,NomMembre,ObjectType,Error)
if (objectType == H5G_GROUP_F .and. NomMembre(1:len(trim(HODIM%GrpElevName))) == trim(HODIM%GrpElevName)) then
    call h5gn_members_f(FileId,NomMembre,NMembers1,Error)
       if ( NMembers1  .gt. 4  )  then
          LRMID=.FALSE.              ! NUMBER OF MEMBERS UNDER /dataset/ GREATER THAN 4  (OIFS FILE )
       else
          LRMID=.TRUE.               ! LOCAL RADAR DATA (RMI FILE )
       endif
endif  
enddo 
endif 



if ( LRMID ==  .TRUE.  ) then  
  write(*,*)   "LRMID "   , LRMID  
  if (Conformity) then
    allocate(FullDatasetList(NDataset),STAT=Error)
    if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate FullDatasetList(:)")
  
    do i=0, NMembers - 1
    call h5gget_obj_info_idx_f(FileId,"/",i,NomMembre,ObjectType,Error)

      if (objectType == H5G_GROUP_F .and. NomMembre(1:len(trim(HODIM%GrpElevName))) == trim(HODIM%GrpElevName)) then
        read (NomMembre(len(trim(HODIM%GrpElevName))+1:len(trim(NomMembre))),'(I3)') NumGDataset
        FullDatasetList(NumGDataset)%Label = trim(NomMembre)
        call h5gn_members_f(FileId,NomMembre,NMembers1,Error)
        NData    = 0
!        NQuality = 0

        do j=0, NMembers1 - 1
          call h5gget_obj_info_idx_f(FileId,NomMembre,j,NomMembre1,ObjectType,Error)
      
        if (ObjectType == H5G_GROUP_F .and. NomMembre1(1:len(trim(HODIM%GrpParamName))) == trim(HODIM%GrpParamName)) then
             NData = NData + 1   
        endif
! IDENTIFICATION OF QUALITY FALGS GROUPS (under /dataset[N]/data1)     
         NDset= NomMembre(len(trim(HODIM%GrpElevName))+1:len(trim(NomMembre))) 
         call h5gn_members_f(FileId, "/dataset"//trim(NDset)//"/data1" , NMembers2 ,Error)

          NQuality=0 
          ! ADD A SECOND LOOP FOR /data1/  SUBGROUP  
          do  ii =0,NMembers2-1
             call  h5gget_obj_info_idx_f(FileId, "/dataset"//trim(NDset)//"/data1"  , ii ,NomMembre2 , ObjectType,Error)

           if (ObjectType == H5G_GROUP_F .and. NomMembre2(1:len(trim(HODIM%GrpFlagName))) == trim(HODIM%GrpFlagName)) then
             NQuality = NQuality + 1
            endif     
          enddo  
         
        enddo
#ifdef DEBUG
        write (Batout,'(/,"Dataset ",A," contains",I3," data and",I3," quality.")') trim(NomMembre),NData,NQuality
#endif
        if (NData == 0) then
          Conformity = .FALSE.
        else
          allocate(FullDatasetList(NumGDataset)%GData(NData),STAT=Error)
          if (error /= 0)  call Abor1("** ERROR - BATOR : cannot allocate GData(:)")
          if (NQuality > 0) then
            allocate(FullDatasetList(NumGDataset)%GQuality(NQuality),STAT=Error)
            if (error /= 0)  call Abor1("** ERROR - BATOR : cannot allocate GQuality(:)")
          endif
          do j=0, NMembers1 - 1
            call h5gget_obj_info_idx_f(FileId,NomMembre,j,NomMembre1,ObjectType,Error)
            if (ObjectType == H5G_GROUP_F .and. NomMembre1(1:len(trim(HODIM%GrpParamName))) == trim(HODIM%GrpParamName)) then
              read (NomMembre1(len(trim(HODIM%GrpParamName))+1:len(trim(NomMembre1))),'(I3)') NumGData
              FullDatasetList(NumGDataset)%GData(NumGData)%Label = trim(NomMembre)//'/'//trim(NomMembre1)
            endif
          enddo  
          call h5gn_members_f(FileId, "/"//trim(NomMembre)//trim(NDset)//"/data1" , NMembers2 ,Error)
           do  ii =0,NMembers2-1
          call  h5gget_obj_info_idx_f(FileId,"/dataset"//trim(NDset)//"/data1"  , ii ,NomMembre2 , ObjectType,Error)
          ! MODIFICATION  (Looking for flag attributes and  the data1 subgroup under dataset group  )
          !                retrieving flags from  /dataset[X]/data1/quality[X]  rather
          !                than  /dataset[X]/quality[Y]
          if (ObjectType == H5G_GROUP_F .and. NomMembre2(1:len(trim(HODIM%GrpFlagName))) == trim(HODIM%GrpFlagName)) then
              read (NomMembre2(len(trim(HODIM%GrpFlagName))+1:len(trim(NomMembre2))),'(I3)') NumGQuality
              FullDatasetList(NumGDataset)%GQuality(NumGQuality)%Label =trim(NomMembre)//'/data1/'//trim(NomMembre2)
            endif
          enddo
          call GetDAttributes(FileId,NomMembre,FullDatasetList(NumGDataset)%Attrib)


#ifdef DEBUG
write (Batout,'(/,"*** INFO - BATOR : DatasetName",10x,": ",A)') trim(FullDatasetList(NumGDataset)%Label)
write (Batout,'("*** INFO - BATOR : startdate",12x,": ",A8)')    FullDatasetList(NumGDataset)%Attrib%Startdate
write (Batout,'("*** INFO - BATOR : starttime",12x,": ",A6)')    FullDatasetList(NumGDataset)%Attrib%Starttime
write (Batout,'("*** INFO - BATOR : quantity",13x,": ",A)')      trim(FullDatasetList(NumGDataset)%Attrib%Quantity)
write (Batout,'("*** INFO - BATOR : gain",17x,": ",F11.6)')      FullDatasetList(NumGDataset)%Attrib%gain
write (Batout,'("*** INFO - BATOR : offset",15x,": ",F11.6)')    FullDatasetList(NumGDataset)%Attrib%offset
write (Batout,'("*** INFO - BATOR : nodata",15x,": ",F11.6)')    FullDatasetList(NumGDataset)%Attrib%nodata
write (Batout,'("*** INFO - BATOR : undetect",13x,": ",F11.6)')  FullDatasetList(NumGDataset)%Attrib%nodetect
write (Batout,'("*** INFO - BATOR : elangle",14x,": ",F11.6)')   FullDatasetList(NumGDataset)%Attrib%elangle
write (Batout,'("*** INFO - BATOR : nrays",16x,": ",I4)')        FullDatasetList(NumGDataset)%Attrib%NRays
write (Batout,'("*** INFO - BATOR : nbins",16x,": ",I4)')        FullDatasetList(NumGDataset)%Attrib%nbins
write (Batout,'("*** INFO - BATOR : rscale",15x,": ",F11.6)')    FullDatasetList(NumGDataset)%Attrib%rscale
write (Batout,'("*** INFO - BATOR : rstart",15x,": ",F11.6)')    FullDatasetList(NumGDataset)%Attrib%rstart
write (Batout,'("*** INFO - BATOR : beamwidth",12x,": ",F11.6)') FullDatasetList(NumGDataset)%Attrib%beamwidth
#endif
          do j=1,size(FullDatasetList(NumGDataset)%GData)
            call GetDAttributes(FileId,FullDatasetList(NumGDataset)%GData(j)%Label,FullDatasetList(NumGDataset)%GData(j)%Attrib)
#ifdef DEBUG
write (Batout,'(/,"*** INFO - BATOR : DataName",13x,": ",A)')    trim(FullDatasetList(NumGDataset)%GData(j)%Label)
write (Batout,'("*** INFO - BATOR : startdate",12x,": ",A8)')    FullDatasetList(NumGDataset)%GData(j)%Attrib%Startdate
write (Batout,'("*** INFO - BATOR : starttime",12x,": ",A6)')    FullDatasetList(NumGDataset)%GData(j)%Attrib%Starttime
write (Batout,'("*** INFO - BATOR : quantity",13x,": ",A)')      trim(FullDatasetList(NumGDataset)%GData(j)%Attrib%Quantity)
write (Batout,'("*** INFO - BATOR : gain",17x,": ",F11.6)')      FullDatasetList(NumGDataset)%GData(j)%Attrib%gain
write (Batout,'("*** INFO - BATOR : offset",15x,": ",F11.6)')    FullDatasetList(NumGDataset)%GData(j)%Attrib%offset
write (Batout,'("*** INFO - BATOR : nodata",15x,": ",F11.6)')    FullDatasetList(NumGDataset)%GData(j)%Attrib%nodata
write (Batout,'("*** INFO - BATOR : undetect",13x,": ",F11.6)')  FullDatasetList(NumGDataset)%GData(j)%Attrib%nodetect
write (Batout,'("*** INFO - BATOR : elangle",14x,": ",F11.6)')   FullDatasetList(NumGDataset)%GData(j)%Attrib%elangle
write (Batout,'("*** INFO - BATOR : nrays",16x,": ",I4)')        FullDatasetList(NumGDataset)%GData(j)%Attrib%NRays
write (Batout,'("*** INFO - BATOR : nbins",16x,": ",I4)')        FullDatasetList(NumGDataset)%GData(j)%Attrib%nbins
write (Batout,'("*** INFO - BATOR : rscale",15x,": ",F11.6)')    FullDatasetList(NumGDataset)%GData(j)%Attrib%rscale
write (Batout,'("*** INFO - BATOR : rstart",15x,": ",F11.6)')    FullDatasetList(NumGDataset)%GData(j)%Attrib%rstart
write (Batout,'("*** INFO - BATOR : beamwidth",12x,": ",F11.6)') FullDatasetList(NumGDataset)%GData(j)%Attrib%beamwidth
#endif
            call BuildAttributesList(Radar%Attrib,FullDatasetList(NumGDataset)%Attrib,FullDatasetList(NumGDataset)%GData(j)%Attrib)
#ifdef DEBUG
write (Batout,'(/,"*** INFO - BATOR : DataName",13x,": ",A)')    trim(FullDatasetList(NumGDataset)%GData(j)%Label)
write (Batout,'("*** INFO - BATOR : startdate",12x,": ",A8)')    FullDatasetList(NumGDataset)%GData(j)%Attrib%Startdate
write (Batout,'("*** INFO - BATOR : starttime",12x,": ",A6)')    FullDatasetList(NumGDataset)%GData(j)%Attrib%Starttime
write (Batout,'("*** INFO - BATOR : quantity",13x,": ",A)')      trim(FullDatasetList(NumGDataset)%GData(j)%Attrib%Quantity)
write (Batout,'("*** INFO - BATOR : gain",17x,": ",F11.6)')      FullDatasetList(NumGDataset)%GData(j)%Attrib%gain
write (Batout,'("*** INFO - BATOR : offset",15x,": ",F11.6)')    FullDatasetList(NumGDataset)%GData(j)%Attrib%offset
write (Batout,'("*** INFO - BATOR : nodata",15x,": ",F11.6)')    FullDatasetList(NumGDataset)%GData(j)%Attrib%nodata
write (Batout,'("*** INFO - BATOR : undetect",13x,": ",F11.6)')  FullDatasetList(NumGDataset)%GData(j)%Attrib%nodetect
write (Batout,'("*** INFO - BATOR : elangle",14x,": ",F11.6)')   FullDatasetList(NumGDataset)%GData(j)%Attrib%elangle
write (Batout,'("*** INFO - BATOR : nrays",16x,": ",I4)')        FullDatasetList(NumGDataset)%GData(j)%Attrib%NRays
write (Batout,'("*** INFO - BATOR : nbins",16x,": ",I4)')        FullDatasetList(NumGDataset)%GData(j)%Attrib%nbins
write (Batout,'("*** INFO - BATOR : rscale",15x,": ",F11.6)')    FullDatasetList(NumGDataset)%GData(j)%Attrib%rscale
write (Batout,'("*** INFO - BATOR : rstart",15x,": ",F11.6)')    FullDatasetList(NumGDataset)%GData(j)%Attrib%rstart
write (Batout,'("*** INFO - BATOR : beamwidth",12x,": ",F11.6)') FullDatasetList(NumGDataset)%GData(j)%Attrib%beamwidth
#endif
          enddo
          do j=1,size(FullDatasetList(NumGDataset)%GQuality)
            call GetFAttributes(FileId,FullDatasetList(NumGDataset)%GQuality(j)%Label, &
                                   & FullDatasetList(NumGDataset)%GQuality(j)%Attrib)
            call BuildAttributesList(Radar%Attrib,FullDatasetList(NumGDataset)%Attrib, &
                                   & FullDatasetList(NumGDataset)%GQuality(j)%Attrib)
#ifdef DEBUG
write (Batout,'(/,"*** INFO - BATOR : gain",17x,": ",F11.6)') FullDatasetList(NumGDataset)%GQuality(j)%Attrib%gain
write (Batout,'("*** INFO - BATOR : offset",15x,": ",F11.6)') FullDatasetList(NumGDataset)%GQuality(j)%Attrib%offset
write (Batout,'("*** INFO - BATOR : task",17x,": ",A/)')      trim(FullDatasetList(NumGDataset)%GQuality(j)%Attrib%task)
#endif
          enddo
       endif
      endif
    enddo
  endif     ! END TEST ON conformity  

elseif ( LRMID== .FALSE. )then         ! OIFS FILE 
if (Conformity) then
    allocate(FullDatasetList(NDataset),STAT=Error)
    if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate FullDatasetList(:)")
    do i=0, NMembers - 1
      call h5gget_obj_info_idx_f(FileId,"/",i,NomMembre,ObjectType,Error)
      if (objectType == H5G_GROUP_F .and. NomMembre(1:len(trim(HODIM%GrpElevName))) == trim(HODIM%GrpElevName)) then
        read (NomMembre(len(trim(HODIM%GrpElevName))+1:len(trim(NomMembre))),'(I3)') NumGDataset
        FullDatasetList(NumGDataset)%Label = trim(NomMembre)
        call h5gn_members_f(FileId,NomMembre,NMembers1,Error)
        NData    = 0
        NQuality = 0
        do j=0, NMembers1 - 1
          call h5gget_obj_info_idx_f(FileId,NomMembre,j,NomMembre1,ObjectType,Error)
          if (ObjectType == H5G_GROUP_F .and. NomMembre1(1:len(trim(HODIM%GrpParamName))) == trim(HODIM%GrpParamName)) then
            NData = NData + 1
          endif
          if (ObjectType == H5G_GROUP_F .and. NomMembre1(1:len(trim(HODIM%GrpFlagName))) == trim(HODIM%GrpFlagName)) then
            NQuality = NQuality + 1
          endif
        enddo
#ifdef DEBUG
        write (Batout,'(/,"Dataset ",A," contains",I3," data and",I3," quality.")') trim(NomMembre),NData,NQuality
#endif
        if (NData == 0) then
          Conformity = .FALSE.
        else
          allocate(FullDatasetList(NumGDataset)%GData(NData),STAT=Error)
          if (error /= 0)  call Abor1("** ERROR - BATOR : cannot allocate GData(:)")
          if (NQuality > 0) then
            allocate(FullDatasetList(NumGDataset)%GQuality(NQuality),STAT=Error)
            if (error /= 0)  call Abor1("** ERROR - BATOR : cannot allocate GQuality(:)")
          endif
          do j=0, NMembers1 - 1
            call h5gget_obj_info_idx_f(FileId,NomMembre,j,NomMembre1,ObjectType,Error)
            if (ObjectType == H5G_GROUP_F .and. NomMembre1(1:len(trim(HODIM%GrpParamName))) == trim(HODIM%GrpParamName)) then
              read (NomMembre1(len(trim(HODIM%GrpParamName))+1:len(trim(NomMembre1))),'(I3)') NumGData
              FullDatasetList(NumGDataset)%GData(NumGData)%Label = trim(NomMembre)//'/'//trim(NomMembre1)
            endif
            if (ObjectType == H5G_GROUP_F .and. NomMembre1(1:len(trim(HODIM%GrpFlagName))) == trim(HODIM%GrpFlagName)) then
              read (NomMembre1(len(trim(HODIM%GrpFlagName))+1:len(trim(NomMembre1))),'(I3)') NumGQuality
              FullDatasetList(NumGDataset)%GQuality(NumGQuality)%Label = trim(NomMembre)//'/'//trim(NomMembre1)
            endif
          enddo
          call GetDAttributes(FileId,NomMembre,FullDatasetList(NumGDataset)%Attrib)
#ifdef DEBUG
write (Batout,'(/,"*** INFO - BATOR : DatasetName",10x,": ",A)') trim(FullDatasetList(NumGDataset)%Label)
write (Batout,'("*** INFO - BATOR : startdate",12x,": ",A8)')    FullDatasetList(NumGDataset)%Attrib%Startdate
write (Batout,'("*** INFO - BATOR : starttime",12x,": ",A6)')    FullDatasetList(NumGDataset)%Attrib%Starttime
write (Batout,'("*** INFO - BATOR : quantity",13x,": ",A)')      trim(FullDatasetList(NumGDataset)%Attrib%Quantity)
write (Batout,'("*** INFO - BATOR : gain",17x,": ",F11.6)')      FullDatasetList(NumGDataset)%Attrib%gain
write (Batout,'("*** INFO - BATOR : offset",15x,": ",F11.6)')    FullDatasetList(NumGDataset)%Attrib%offset
write (Batout,'("*** INFO - BATOR : nodata",15x,": ",F11.6)')    FullDatasetList(NumGDataset)%Attrib%nodata
write (Batout,'("*** INFO - BATOR : undetect",13x,": ",F11.6)')  FullDatasetList(NumGDataset)%Attrib%nodetect
write (Batout,'("*** INFO - BATOR : elangle",14x,": ",F11.6)')   FullDatasetList(NumGDataset)%Attrib%elangle
write (Batout,'("*** INFO - BATOR : nrays",16x,": ",I4)')        FullDatasetList(NumGDataset)%Attrib%NRays
write (Batout,'("*** INFO - BATOR : nbins",16x,": ",I4)')        FullDatasetList(NumGDataset)%Attrib%nbins
write (Batout,'("*** INFO - BATOR : rscale",15x,": ",F11.6)')    FullDatasetList(NumGDataset)%Attrib%rscale
write (Batout,'("*** INFO - BATOR : rstart",15x,": ",F11.6)')    FullDatasetList(NumGDataset)%Attrib%rstart
write (Batout,'("*** INFO - BATOR : beamwidth",12x,": ",F11.6)') FullDatasetList(NumGDataset)%Attrib%beamwidth
#endif
          do j=1,size(FullDatasetList(NumGDataset)%GData)
            call GetDAttributes(FileId,FullDatasetList(NumGDataset)%GData(j)%Label,FullDatasetList(NumGDataset)%GData(j)%Attrib)
#ifdef DEBUG
write (Batout,'(/,"*** INFO - BATOR : DataName",13x,": ",A)')    trim(FullDatasetList(NumGDataset)%GData(j)%Label)
write (Batout,'("*** INFO - BATOR : startdate",12x,": ",A8)')    FullDatasetList(NumGDataset)%GData(j)%Attrib%Startdate
write (Batout,'("*** INFO - BATOR : starttime",12x,": ",A6)')    FullDatasetList(NumGDataset)%GData(j)%Attrib%Starttime
write (Batout,'("*** INFO - BATOR : quantity",13x,": ",A)')      trim(FullDatasetList(NumGDataset)%GData(j)%Attrib%Quantity)
write (Batout,'("*** INFO - BATOR : gain",17x,": ",F11.6)')      FullDatasetList(NumGDataset)%GData(j)%Attrib%gain
write (Batout,'("*** INFO - BATOR : offset",15x,": ",F11.6)')    FullDatasetList(NumGDataset)%GData(j)%Attrib%offset
write (Batout,'("*** INFO - BATOR : nodata",15x,": ",F11.6)')    FullDatasetList(NumGDataset)%GData(j)%Attrib%nodata
write (Batout,'("*** INFO - BATOR : undetect",13x,": ",F11.6)')  FullDatasetList(NumGDataset)%GData(j)%Attrib%nodetect
write (Batout,'("*** INFO - BATOR : elangle",14x,": ",F11.6)')   FullDatasetList(NumGDataset)%GData(j)%Attrib%elangle
write (Batout,'("*** INFO - BATOR : nrays",16x,": ",I4)')        FullDatasetList(NumGDataset)%GData(j)%Attrib%NRays
write (Batout,'("*** INFO - BATOR : nbins",16x,": ",I4)')        FullDatasetList(NumGDataset)%GData(j)%Attrib%nbins
write (Batout,'("*** INFO - BATOR : rscale",15x,": ",F11.6)')    FullDatasetList(NumGDataset)%GData(j)%Attrib%rscale
write (Batout,'("*** INFO - BATOR : rstart",15x,": ",F11.6)')    FullDatasetList(NumGDataset)%GData(j)%Attrib%rstart
write (Batout,'("*** INFO - BATOR : beamwidth",12x,": ",F11.6)') FullDatasetList(NumGDataset)%GData(j)%Attrib%beamwidth
#endif
            call BuildAttributesList(Radar%Attrib,FullDatasetList(NumGDataset)%Attrib,FullDatasetList(NumGDataset)%GData(j)%Attrib)
#ifdef DEBUG
write (Batout,'(/,"*** INFO - BATOR : DataName",13x,": ",A)')    trim(FullDatasetList(NumGDataset)%GData(j)%Label)
write (Batout,'("*** INFO - BATOR : startdate",12x,": ",A8)')    FullDatasetList(NumGDataset)%GData(j)%Attrib%Startdate
write (Batout,'("*** INFO - BATOR : starttime",12x,": ",A6)')    FullDatasetList(NumGDataset)%GData(j)%Attrib%Starttime
write (Batout,'("*** INFO - BATOR : quantity",13x,": ",A)')      trim(FullDatasetList(NumGDataset)%GData(j)%Attrib%Quantity)
write (Batout,'("*** INFO - BATOR : gain",17x,": ",F11.6)')      FullDatasetList(NumGDataset)%GData(j)%Attrib%gain
write (Batout,'("*** INFO - BATOR : offset",15x,": ",F11.6)')    FullDatasetList(NumGDataset)%GData(j)%Attrib%offset
write (Batout,'("*** INFO - BATOR : nodata",15x,": ",F11.6)')    FullDatasetList(NumGDataset)%GData(j)%Attrib%nodata
write (Batout,'("*** INFO - BATOR : undetect",13x,": ",F11.6)')  FullDatasetList(NumGDataset)%GData(j)%Attrib%nodetect
write (Batout,'("*** INFO - BATOR : elangle",14x,": ",F11.6)')   FullDatasetList(NumGDataset)%GData(j)%Attrib%elangle
write (Batout,'("*** INFO - BATOR : nrays",16x,": ",I4)')        FullDatasetList(NumGDataset)%GData(j)%Attrib%NRays
write (Batout,'("*** INFO - BATOR : nbins",16x,": ",I4)')        FullDatasetList(NumGDataset)%GData(j)%Attrib%nbins
write (Batout,'("*** INFO - BATOR : rscale",15x,": ",F11.6)')    FullDatasetList(NumGDataset)%GData(j)%Attrib%rscale
write (Batout,'("*** INFO - BATOR : rstart",15x,": ",F11.6)')    FullDatasetList(NumGDataset)%GData(j)%Attrib%rstart
write (Batout,'("*** INFO - BATOR : beamwidth",12x,": ",F11.6)') FullDatasetList(NumGDataset)%GData(j)%Attrib%beamwidth
#endif
          enddo
          do j=1,size(FullDatasetList(NumGDataset)%GQuality)
            call GetFAttributes(FileId,FullDatasetList(NumGDataset)%GQuality(j)%Label, &
                                   & FullDatasetList(NumGDataset)%GQuality(j)%Attrib)
            call BuildAttributesList(Radar%Attrib,FullDatasetList(NumGDataset)%Attrib, &
                                   & FullDatasetList(NumGDataset)%GQuality(j)%Attrib)
#ifdef DEBUG
write (Batout,'(/,"*** INFO - BATOR : gain",17x,": ",F11.6)') FullDatasetList(NumGDataset)%GQuality(j)%Attrib%gain
write (Batout,'("*** INFO - BATOR : offset",15x,": ",F11.6)') FullDatasetList(NumGDataset)%GQuality(j)%Attrib%offset
write (Batout,'("*** INFO - BATOR : task",17x,": ",A/)')      trim(FullDatasetList(NumGDataset)%GQuality(j)%Attrib%task)
#endif
          enddo
       endif
      endif
    enddo
  endif   ! Conformity  
endif    ! END TEST ON FILE STRUCTURE  

 

! get the first two proportional most popular nrays
  if (Conformity) then
    allocate(NRaysPopulation(size(FullDatasetList)*3,2),STAT=Error)
    if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate NRaysPopulation(:,:)")
    NRaysPopulation(:,1) = nabsi
    NRaysPopulation(:,2) = 0
    SelectedNRays        = 0
    NbNrays              = 0
    NBConcernedData      = 0
    do i=1, size(FullDatasetlist)
      do j=1, size(FullDatasetList(i)%GData)
        if (any((/'TH   ','DBZH ','VRAD ','VRADH'/) == trim(FullDatasetList(i)%GData(j)%Attrib%Quantity))) then
          if (any(NraysPopulation(:,1) == FullDatasetList(i)%GData(j)%Attrib%NRays)) then
            do k=1, NbNrays 
              if (NraysPopulation(k,1) == FullDatasetList(i)%GData(j)%Attrib%NRays) NraysPopulation(k,2) = NraysPopulation(k,2)+1
            enddo
          else
            NbNrays = NbNrays + 1
            NraysPopulation(NbNrays,1) = FullDatasetList(i)%GData(j)%Attrib%NRays
            NraysPopulation(NbNrays,2) = NraysPopulation(NbNrays,2)+1
          endif
        else
          FullDatasetList(i)%GData(j)%Attrib%StartTime = '?'
          FullDatasetList(i)%GData(j)%Attrib%StartDate = '?'
          FullDatasetList(i)%GData(j)%Attrib%Quantity  = '?'
          FullDatasetList(i)%GData(j)%Attrib%Gain      = rabsi
          FullDatasetList(i)%GData(j)%Attrib%Offset    = rabsi
          FullDatasetList(i)%GData(j)%Attrib%NoData    = rabsi
          FullDatasetList(i)%GData(j)%Attrib%NoDetect  = rabsi
          FullDatasetList(i)%GData(j)%Attrib%Elangle   = rabsi
          FullDatasetList(i)%GData(j)%Attrib%NBins     = nabsi
          FullDatasetList(i)%GData(j)%Attrib%NRays     = nabsi
          FullDatasetList(i)%GData(j)%Attrib%RScale    = rabsi
          FullDatasetList(i)%GData(j)%Attrib%Rstart    = rabsi
          FullDatasetList(i)%GData(j)%Attrib%Task      = '?'
        endif
      enddo
    enddo
    SelectedNRays(1) = NraysPopulation(maxloc(NraysPopulation(:,2),DIM=1),1)
    NbConcernedData = maxval(NraysPopulation(:,2))
    NraysPopulation(maxloc(NraysPopulation(:,2),DIM=1),:) = nabsi
    do while (any(NraysPopulation(:,1) /= nabsi) .and. SelectedNRays(2) == 0)
      SelectedNRays(2) = NraysPopulation(maxloc(NraysPopulation(:,2),DIM=1),1)
      if (mod(maxval(SelectedNRays),minval(SelectedNRays)) /= 0) then
        SelectedNRays(2) = 0
        NraysPopulation(maxloc(NraysPopulation(:,2),DIM=1),2) = 0
        NraysPopulation(maxloc(NraysPopulation(:,2),DIM=1),1) = nabsi
      else
         NbConcernedData = NbConcernedData + maxval(NraysPopulation(:,2))
      endif
    enddo
    if (all(SelectedNRays <= 0)) Conformity = .FALSE.
    Radar%NRayons = maxval(SelectedNRays)
  endif

! select closest DBZH (+FLAG), TH, VRAD to analysis date for each elevations
  if (Conformity) then
    allocate(SelectedElangles(NbConcernedData),STAT=Error)
    if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate SelectedElangles(:")
    NbSelectedElangles = 0
    MaxRayLength       = 0
    do i=1, size(FullDatasetlist)
      do j=1, size(FullDatasetList(i)%GData)
        if (any(SelectedNRays == FullDatasetList(i)%GData(j)%Attrib%NRays)) then
          if (FullDatasetList(i)%GData(j)%Attrib%RStart == 0._jprb)  FullDatasetList(i)%GData(j)%Attrib%NBins = &
                                                                  & FullDatasetList(i)%GData(j)%Attrib%NBins + 1
          RayLength = FullDatasetList(i)%GData(j)%Attrib%RStart * 1000._jprb + &
                    & FullDatasetList(i)%GData(j)%Attrib%RScale * FullDatasetList(i)%GData(j)%Attrib%NBins
          DeltaTime = abs(DiffDate(NewDate(FullDatasetList(i)%GData(j)%Attrib%StartDate, &
                                 & FullDatasetList(i)%GData(j)%Attrib%StartTime),AnalysisDate))
          if (.not.any(SelectedElangles%Elangle == FullDatasetList(i)%GData(j)%Attrib%Elangle) .and. &
            & FullDatasetList(i)%GData(j)%Attrib%Elangle < 80._jprb ) then
            NbSelectedElangles = NbSelectedElangles + 1
            SelectedElangles(NbSelectedElangles)%Elangle = FullDatasetList(i)%GData(j)%Attrib%Elangle
          endif
          do k=1, NbSelectedElangles
            if (SelectedElangles(k)%Elangle == FullDatasetList(i)%GData(j)%Attrib%Elangle) then
              if (FullDatasetList(i)%GData(j)%Attrib%Quantity == "DBZH" .and. DeltaTime < SelectedElangles(k)%DbZHDeltaTime) then
                SelectedElangles(k)%DbZHDeltaTime = DeltaTime
                SelectedElangles(k)%DBZH => FullDatasetList(i)%GData(j)
                MaxRayLength  = max(RayLength,MaxRayLength)
                if (associated(FullDatasetList(i)%GQuality)) then
                  do l=1, size(FullDatasetList(i)%GQuality)
                    if (FullDatasetList(i)%GQuality(l)%Attrib%Task == HODIM%ChoosenTask) then
                      SelectedElangles(k)%FLAG => FullDatasetList(i)%GQuality(l)
                    endif
                  enddo
                else
                  nullify(SelectedElangles(k)%FLAG)
                endif
              endif
              if (FullDatasetList(i)%GData(j)%Attrib%Quantity == "TH" .and. DeltaTime < SelectedElangles(k)%THDeltaTime) then
                SelectedElangles(k)%THDeltaTime = DeltaTime
                SelectedElangles(k)%TH => FullDatasetList(i)%GData(j)
                MaxRayLength  = max(RayLength,MaxRayLength)
              endif
              if ((FullDatasetList(i)%GData(j)%Attrib%Quantity=="VRAD" .or. FullDatasetList(i)%GData(j)%Attrib%Quantity=="VRADH") &
                 & .and. DeltaTime < SelectedElangles(k)%VRadDeltaTime) then
                SelectedElangles(k)%VRadDeltaTime = DeltaTime
                SelectedElangles(k)%VRAD => FullDatasetList(i)%GData(j)
                MaxRayLength  = max(RayLength,MaxRayLength)
              endif
            endif
          enddo
        endif
      enddo
    enddo
    Radar%NPoints = nint(MaxRayLength/HODIM%Resolution)
    allocate(Radar%FinalElev(NbSelectedElangles),STAT=Error)
    if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate Radar%FinalElev(:)")
    write(Batout,'("*** INFO - BATOR : rays1= ",I4," rays2= ",I4)') SelectedNRays(1), SelectedNRays(2)
    write(Batout,'("*** INFO - BATOR : Radar%NRayons= ",I4," Radar%NPoints= ",I4)') Radar%NRayons, Radar%NPoints
! sort elevations, get data for each param
    do i=1, NbSelectedElangles
      CurrentMinElangle = minloc(abs(SelectedElangles(1:NbSelectedElangles)%Elangle),DIM=1)
      Radar%FinalElev(i) = SelectedElangles(CurrentMinElangle)
      SelectedElangles(CurrentMinElangle)%Elangle = rabsi
      write(Batout,'("*** INFO - BATOR : elevation = ",F11.6)') Radar%FinalElev(i)%Elangle
      if (associated(Radar%FinalElev(i)%DBZH)) then
        call h5gopen_f(FileId,Radar%FinalElev(i)%DBZH%Label,GroupId,Error)
        WaitedRank   = 2
        WaitedAtomic = Hdf5Type(-1)
        call GetData(GroupId,HODIM%GrpParamName,WaitedRank,WaitedAtomic,Iret)
        where (real2buf(:,:) == Radar%FinalElev(i)%DBZH%Attrib%NoData)   real2buf(:,:) = rabso
        where (real2buf(:,:) == Radar%FinalElev(i)%DBZH%Attrib%NoDetect) real2buf(:,:) = -rabsi
        where (real2buf(:,:) /= -rabsi .and. real2buf(:,:) /= rabso) &
             & real2buf(:,:) = real2buf(:,:) * Radar%FinalElev(i)%DBZH%Attrib%gain + Radar%FinalElev(i)%DBZH%Attrib%offset
        allocate(Radar%FinalElev(i)%DBZH%Values(Radar%NPoints,Radar%NRayons),STAT=Error)
        if (Error /= 0) call Abor1("* ERROR - BATOR : cannot allocate Radar%FinalElev(i)%DBZH%Values(:,:)")
        Radar%FinalElev(i)%DBZH%Values = rabsi
        call h5gclose_f(GroupId,Error)
        if (Radar%FinalElev(i)%DBZH%Attrib%RStart == 0._jprb) then
          Radar%FinalElev(i)%DBZH%Attrib%RStart = Radar%FinalElev(i)%DBZH%Attrib%RScale
        else
          Radar%FinalElev(i)%DBZH%Attrib%RStart = Radar%FinalElev(i)%DBZH%Attrib%RStart * 1000._jprb
        endif
        write(Batout,'("*** INFO - BATOR : selected dbzh = ",A)') trim(Radar%FinalElev(i)%DBZH%Label)
        call maps(Radar%FinalElev(i)%DBZH%Attrib%RScale,Radar%FinalElev(i)%DBZH%Attrib%RStart,&
                & Radar%FinalElev(i)%Elangle,real2buf,Radar%FinalElev(i)%DBZH%Values)
#ifdef DEBUG
        print *, radar%FinalElev(i)%DBZH%Values(1:10,1)
#endif
        if (associated(real2buf)) deallocate(real2buf)
      endif
      if (associated(Radar%FinalElev(i)%TH)) then 
        call h5gopen_f(FileId,Radar%FinalElev(i)%TH%Label,GroupId,Error)
        WaitedRank   = 2
        WaitedAtomic = Hdf5Type(-1)
        call GetData(GroupId,HODIM%GrpParamName,WaitedRank,WaitedAtomic,Iret)
        where (real2buf(:,:) == Radar%FinalElev(i)%TH%Attrib%NoData)   real2buf(:,:) = rabso
        where (real2buf(:,:) == Radar%FinalElev(i)%TH%Attrib%NoDetect) real2buf(:,:) = -rabsi
        where (real2buf(:,:) /= -rabsi .and. real2buf(:,:) /= rabso) &
             & real2buf(:,:) = real2buf(:,:) * Radar%FinalElev(i)%TH%Attrib%gain + Radar%FinalElev(i)%TH%Attrib%offset
        allocate(Radar%FinalElev(i)%TH%Values(Radar%NPoints,Radar%NRayons),STAT=Error)
        if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate Radar%FinalElev(i)%TH%Values(:,:)")
        Radar%FinalElev(i)%TH%Values = rabsi
        call h5gclose_f(GroupId,Error)
        if (Radar%FinalElev(i)%TH%Attrib%RStart == 0._jprb) then
          Radar%FinalElev(i)%TH%Attrib%RStart = Radar%FinalElev(i)%TH%Attrib%RScale
        else
          Radar%FinalElev(i)%TH%Attrib%RStart = Radar%FinalElev(i)%TH%Attrib%RStart * 1000._jprb
        endif
        write(Batout,'("*** INFO - BATOR : selected th = ",A)') trim(Radar%FinalElev(i)%TH%Label)
        call maps(Radar%FinalElev(i)%TH%Attrib%RScale,Radar%FinalElev(i)%TH%Attrib%RStart,&
                & Radar%FinalElev(i)%Elangle,real2buf,Radar%FinalElev(i)%TH%Values)
#ifdef DEBUG
        print *, radar%FinalElev(i)%TH%Values(1:10,1)
#endif
        if (associated(real2buf)) deallocate(real2buf)
      endif
      if (associated(Radar%FinalElev(i)%VRAD)) then
        call h5gopen_f(FileId,Radar%FinalElev(i)%VRAD%Label,GroupId,Error)
        WaitedRank   = 2
        WaitedAtomic = Hdf5Type(-1)
        call GetData(GroupId,HODIM%GrpParamName,WaitedRank,WaitedAtomic,Iret)
        where (real2buf(:,:) == Radar%FinalElev(i)%VRAD%Attrib%NoData)   real2buf(:,:) = rabso
        where (real2buf(:,:) == Radar%FinalElev(i)%VRAD%Attrib%NoDetect) real2buf(:,:) = -rabsi
        where (real2buf(:,:) /= -rabsi .and. real2buf(:,:) /= rabso) &
             & real2buf(:,:) = real2buf(:,:) * Radar%FinalElev(i)%VRAD%Attrib%gain + Radar%FinalElev(i)%VRAD%Attrib%offset
        allocate(Radar%FinalElev(i)%VRAD%Values(Radar%NPoints,Radar%NRayons),STAT=Error)
        if (Error /= 0) call Abor1("* ERROR - BATOR : cannot allocate Radar%FinalElev(i)%VRAD%Values(:,:)")
        Radar%FinalElev(i)%VRAD%Values = rabsi
        call h5gclose_f(GroupId,Error)
        if (Radar%FinalElev(i)%VRAD%Attrib%RStart == 0._jprb) then
          Radar%FinalElev(i)%VRAD%Attrib%RStart = Radar%FinalElev(i)%VRAD%Attrib%RScale
        else
          Radar%FinalElev(i)%VRAD%Attrib%RStart = Radar%FinalElev(i)%VRAD%Attrib%RStart * 1000._jprb
        endif
        write(Batout,'("*** INFO - BATOR : selected vrad = ",A)') trim(Radar%FinalElev(i)%VRAD%Label)
        call maps(Radar%FinalElev(i)%VRAD%Attrib%RScale,Radar%FinalElev(i)%VRAD%Attrib%RStart,&
                & Radar%FinalElev(i)%Elangle,real2buf,Radar%FinalElev(i)%VRAD%Values)
#ifdef DEBUG
        print *, radar%FinalElev(i)%VRAD%Values(1:10,1)
#endif
        if (associated(real2buf)) deallocate(real2buf)
      endif
      if (associated(Radar%FinalElev(i)%FLAG)) then
        call h5gopen_f(FileId,Radar%FinalElev(i)%FLAG%Label,GroupId,Error)
        WaitedRank   = 2
        WaitedAtomic = Hdf5Type(-1)
        call GetData(GroupId,HODIM%GrpParamName,WaitedRank,WaitedAtomic,Iret)
        real2buf(:,:) = real2buf(:,:) * Radar%FinalElev(i)%FLAG%Attrib%gain + Radar%FinalElev(i)%FLAG%Attrib%offset
        allocate(Radar%FinalElev(i)%FLAG%Values(Radar%NPoints,Radar%NRayons),STAT=Error)
        if (Error /= 0) call Abor1("* ERROR - BATOR : cannot allocate Radar%FinalElev(i)%FLAG%Values(:,:)")
        Radar%FinalElev(i)%FLAG%Values = rabso
        call h5gclose_f(GroupId,Error)
        write(Batout,'("*** INFO - BATOR : selected flag = ",A)') trim(Radar%FinalElev(i)%FLAG%Label)
        call maps(Radar%FinalElev(i)%DBZH%Attrib%RScale,Radar%FinalElev(i)%DBZH%Attrib%RStart,&
                & Radar%FinalElev(i)%Elangle,real2buf,Radar%FinalElev(i)%FLAG%Values)
#ifdef DEBUG
        print *, radar%FinalElev(i)%FLAG%Values(1:10,1)
#endif
        if (associated(real2buf)) deallocate(real2buf)
      endif
    enddo

    RAE        = RA*4./3. ! Effective radius of the Earth
    ispindir   = 1        ! Scanning direction : always clockwise in ODIM 2.2 norm
    zaz_offset = 0.       ! Azimuth offset (confondu avec rstart=r_binoffset dans CONRAD)
    zdimhpiy   = 360./Radar%NRayons
    iobs       = 0
    iw         = 0
    nb_obs     = Radar%NRayons * Radar%NPoints
    nbw        = 3 * NbSelectedElangles * nb_obs
    ilw        = TREF_FICOBS(kfic)%ilwag
    write(batout,'("NbElev : ",I2,"ilw : ",I2,"NbObs : ",I6)') NbSelectedElangles, ilw, nb_obs

! infos pour table RADAR_STATION
    JPXRADS = JPXRADS + 1
    TREF_RADAR(JPXRADS)%Ident = nabso
    WRITE(TREF_RADAR(JPXRADS)%Type,'(3x,A5)') Radar%Nod
    TREF_RADAR(JPXRADS)%lat       = Radar%Lat
    TREF_RADAR(JPXRADS)%lon       = Radar%Lon
    TREF_RADAR(JPXRADS)%stalt     = Radar%Alt       ! Height of the antenna above sea level (zhsta + zhant)
    TREF_RADAR(JPXRADS)%antenht   = 0._jprb         ! zhant
    TREF_RADAR(JPXRADS)%beamwidth = Radar%BeamWidth ! half-3 dB aperture

!    i_indic = TREF_RADAR(JPXRADS)%ident
!    write(batout,*),"i_indic vaut",i_indic

    if (.not.allocated(ztent)) allocate(ztent(nb_obs,ncilet),STAT=Error)
    if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate ztent(:,:)")
    if (.not.allocated(ztwag)) allocate(ztwag(nbw,ilw),STAT=Error)
    if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate ztwag(:,:)")
    if (.not.allocated(zdist)) allocate(zdist(nb_obs),STAT=Error)
    if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate zdist(:)")
    if (.not.allocated(zazim)) allocate(zazim(nb_obs),STAT=Error)
    if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate zazim(:)")
    ztent(:,:)      = rabsi
    ztent(:,ncinlv) = 0
    ztwag(:,:)      =rabso
    zazim(:)        =rabso

    do i= 1, Radar%NRayons
      do j= 1, Radar%NPoints
        iobs = iobs + 1
        inlv = 0
        do l=1, NbSelectedElangles
          if (l == 1) then
!           zdist(iobs) = (radar%Elevations(l)%rscale1 * (j - 0.5) + &
!                       & radar%Elevations(l)%rstart1)**2 
            zdist(iobs) = (1000. * (j - 0.5))**2
            zazim(iobs) = ispindir * zdimhpiy * (i - 0.5) + zaz_offset
            zalt        = RA + sin(Radar%FinalElev(l)%Elangle*RADIANS) * sqrt(zdist(iobs)) + zdist(iobs)/2./RAE
            if (abs((RA**2 + zalt**2 - zdist(iobs)) / (2. * RA * zalt)) <= 1.) then
              zalpha = acos((RA**2 + zalt**2 - zdist(iobs)) / (2. * RA * zalt))
              zlat = 180./RPI * asin( sin(Radar%Lat*RADIANS) * cos(zalpha) + &
                   & cos(Radar%Lat*RADIANS) * sin(zalpha) * cos(zazim(iobs)*RADIANS) )
              zlon = Radar%Lon  + 180./RPI * asin( sin(zalpha) * sin(zazim(iobs)*RADIANS) /cos(zlat*RADIANS) )
            else 
!               CALL Abor1('** ERROR Radar : cannot obtain lat, lon for an obs')
              write(Batout,*) "Attention: certaines obs sans lat,lon"
              exit
            endif
            i0              = 3 * (NbSelectedElangles * (iobs-1)) + 1
            i1              = 3 * NbSelectedElangles * iobs
            ztwag(i0:i1,5)  = zazim(iobs)                             ! azimuth in deg
            ztwag(i0:i1,9)  = sqrt(zdist(iobs))                       ! distance horizontale du point
            ztwag(i0:i1,10) = 0                                       ! polarisation
            ztwag(i0:i1,11) = 0                                       ! anaprop
          else
            zalt = RA + sin(Radar%FinalElev(l)%Elangle*RADIANS) * sqrt(zdist(iobs)) + zdist(iobs)/2./RAE  
          endif

! Detection threshold
          zsensib    = -110
          zconst     =  -71
          zthreshold = zsensib - zconst + 40 +20.0*log10(sqrt(zdist(iobs))/100000.)
          if (associated(Radar%FinalElev(l)%DBZH) .and. associated(Radar%FinalElev(l)%TH) .and. &
            & Radar%FinalElev(l)%DBZH%Values(j,i) /= rabso .and. Radar%FinalElev(l)%TH%Values(j,i) /= rabso) then ! rabsi ou rabso
            if ((sqrt(zdist(iobs))<160000.) .and. (Radar%FinalElev(l)%DBZH%Values(j,i) <= 100.)) then
              iw          = 3*((l-1) + NbSelectedElangles*(Radar%NPoints*(i-1)+j-1))+1
              ztwag(iw,1) = Varno%refl
              ztwag(iw,2) = Radar%FinalElev(l)%Elangle
              ztwag(iw,3) = zalt - RA + TREF_RADAR(JPXRADS)%stalt
              ztwag(iw,5) = zazim(iobs) 
              ztwag(iw,9) = sqrt(zdist(iobs))
              if (associated(Radar%FinalElev(l)%FLAG) .and. Radar%FinalElev(l)%FLAG%Values(j,i) /= rabso ) then
                inlv        = inlv + 1
                ztwag(iw,4) = Radar%FinalElev(l)%DBZH%Values(j,i)
                if (Radar%FinalElev(l)%FLAG%Values(j,i) >= 0.7) then 
                  if (Radar%FinalElev(l)%DBZH%Values(j,i) == -rabsi .and. Radar%FinalElev(l)%TH%Values(j,i) == -rabsi) then
!                  IF  (radar%Elevations(l)%dbzh1(j,i) == -rabsi) THEN 
                    ztwag(iw,8)   = 0.
                    ztwag(iw+1,8) = 0.
                    ztwag(iw,4)   = zthreshold
                  elseif (Radar%FinalElev(l)%DBZH%Values(j,i) > (Radar%FinalElev(l)%TH%Values(j,i) - 3.))  then 
                    ztwag(iw,8)   = 8.
                    ztwag(iw+1,8) = 8.
                  endif
                elseif (Radar%FinalElev(l)%FLAG%Values(j,i) < 0.7) then ! flag < 0.7
                  ztwag(iw,8)   = 1.
                  ztwag(iw+1,8) = 1.
                endif
!           IF  (radar%Elevations(l)%dbzh1(j,i) /= & 
!&              radar%Elevations(l)%th1(j,i)) THEN
!               zwagon(kw,8)=1.
!           ENDIF
        ! RH if Z valid   
                inlv        = inlv + 1
                iw          = 3*( (l-1) + NbSelectedElangles*(Radar%NPoints*(i-1)+j-1) )+2
                ztwag(iw,1) = Varno%rh
                ztwag(iw,2) = Radar%FinalElev(l)%Elangle
                ztwag(iw,3) = zalt - RA + TREF_RADAR(JPXRADS)%stalt
                ztwag(iw,4) = 0.
              endif
!           IF (radar%Elevations(l)%flag1(j,i) /= rabsi) THEN
!           zwagon(kw,8) = radar%Elevations(l)%flag1(j,i)
!           ELSE 
!           zwagon(kw,8) = 0.
!           ENDIF
            endif ! if dist < 160 km and dbzh1<=100

    ! associated dbzh and not th -- qu'est-ce que le i_indic ????
          elseif (associated(Radar%FinalElev(l)%DBZH) .and. Radar%FinalElev(l)%DBZH%Values(j,i) /= rabso) then
!&                  .not.(associated(radar%Elevations(l)%Th1)) .AND. &
!          write(*,'("- station indicateur            : ",I5.5)') i_indic
            if ((sqrt(zdist(iobs)) < 160000.) .and. (Radar%FinalElev(l)%DBZH%Values(j,i) <= 100.)) then
              iw          = 3*( (l-1) + NbSelectedElangles*(Radar%NPoints*(i-1)+j-1) )+1
              ztwag(iw,1) = Varno%refl
              ztwag(iw,2) = Radar%FinalElev(l)%Elangle
              ztwag(iw,3) = zalt - RA + TREF_RADAR(JPXRADS)%stalt
              ztwag(iw,5) = zazim(iobs)
              ztwag(iw,9) = sqrt(zdist(iobs))
              if  (associated(Radar%FinalElev(l)%FLAG) .and. Radar%FinalElev(l)%FLAG%Values(j,i) /= rabso ) then !flag exist
                inlv        = inlv + 1
                ztwag(iw,4) = Radar%FinalElev(l)%DBZH%Values(j,i)
                if (Radar%FinalElev(l)%FLAG%Values(j,i) >= 0.7) then
                  if (Radar%FinalElev(l)%DBZH%Values(j,i) == -rabsi ) then ! rain..
                    ztwag(iw,8) = 0.
                    ztwag(iw,4) = zthreshold
                  else
                    ztwag(iw,8)=8.
                    ztwag(iw+1,8)=8.
                  endif ! rain or no rain
                elseif (Radar%FinalElev(l)%FLAG%Values(j,i) < 0.7) then ! flag < 0.7
                  ztwag(iw,8)   = 3.
                  ztwag(iw+1,8) = 3.
                endif
       ! RH if Z valid
                inlv        = inlv + 1
                iw          = 3*( (l-1) + NbSelectedElangles*(Radar%NPoints*(i-1)+j-1) )+2
                ztwag(iw,1) = Varno%rh
                ztwag(iw,2) = Radar%FinalElev(l)%Elangle
                ztwag(iw,3) = zalt - RA + TREF_RADAR(JPXRADS)%stalt
                ztwag(iw,4) = 0.
              endif ! flag exist
            endif ! if dist < 160 km and dbzh1 <=100
          endif ! if associated dbzh1 and not th1

          if (associated(Radar%FinalElev(l)%VRAD)) then
!        IF ( associated(radar%Elevations(l)%vrad1) .and. &
!& radar%Elevations(l)%vrad1(j,i) /= rabso ) THEN
            iw          = 3*( (l-1) + NbSelectedElangles*(Radar%NPoints*(i-1)+j-1) )+3
            ztwag(iw,1) = Varno%dopp
            ztwag(iw,2) = Radar%FinalElev(l)%Elangle
            ztwag(iw,3) = zalt - RA + TREF_RADAR(JPXRADS)%stalt
            if ((sqrt(zdist(iobs)) < 150000.) .and. (abs(Radar%FinalElev(l)%VRAD%Values(j,i)) <= 100.)) then
!           iw = 3*( (l-1) + NbElangles*(radar%head%NPixels*(i-1)+j-1) )+3
!           ztwag(iw,1) = 195
!           ztwag(iw,2) = radar%Elevations(l)%elangle
              ztwag(iw,4) = -Radar%FinalElev(l)%VRAD%Values(j,i)
              inlv        = inlv + 1
              if (associated(Radar%FinalElev(l)%FLAG) .and. &
                & Radar%FinalElev(l)%FLAG%Values(j,i) /= rabso .and. Radar%FinalElev(l)%FLAG%Values(j,i) >= 0.7) then
                ztwag(iw,8) = 8
              else       
                ztwag(iw,8) = 2
              endif ! associated flag1 for radial wind
            endif ! condition on value and range from the radar
          else ! associated radial wind: nothing
          endif ! associated radial wind         
!         IF (associated(radar%Elevations(l)%dbzh1) .AND. &
!&                  associated(radar%Elevations(l)%Th1)==.false. .AND. &
!&                  radar%Elevations(l)%dbzh1(j,i)/=rabsi THEN 
        enddo 
!     IF (inlv /=0) THEN
        ztent(iobs,ncinlv) = inlv    ! all elevations in the same time: not
                                     ! NCINLV to add before
        ztent(iobs,NCILAT) = zlat    ! Pixel latitude
        ztent(iobs,NCILON) = zlon    ! Pixel longitude
        ztent(iobs,NCIOTP) = NRADAR
        ztent(iobs,NCIRFL) = 11111    ! Flag
        ztent(iobs,NCIOCH) = NRADA1        ! ODB codetype
        ztent(iobs,NCIDAT) = Radar%DateOpt(1)*10000 + Radar%DateOpt(2)*100 + Radar%DateOpt(3)   ! Observation date
        ztent(iobs,NCIETM) = Radar%DateOpt(4)*10000 + Radar%DateOpt(5)*100 + Radar%DateOpt(6)   ! Observation heure
!     ELSE       
!        iobs = iobs - 1   
!     ENDIF
      enddo
    enddo

! Spatial Filtering of DOW
!-------------------------
    do l=1,NbSelectedElangles ! fileltering for each elevation
      if (associated(Radar%FinalElev(l)%VRAD)) then 
        call bator_radar_wind_cleaner('DOPW',ztent,ztwag,nbw,ilw,l,NbSelectedElangles,Radar%NRayons,Radar%NPoints,.TRUE.,empty)
        if (.not.empty) then
          call bator_filter_radar('DOPW',ztent,ztwag,nbw,ilw,l,NbSelectedElangles,Radar%NRayons,Radar%NPoints)
          call bator_radar_wind_cleaner('DOPW',ztent,ztwag,nbw,ilw,l,NbSelectedElangles,Radar%NRayons,Radar%NPoints,.FALSE.,empty)
        else
          write(batout,*) 'Wind data not avalaible. ELEVATION REJECTED'
        endif
      endif
    enddo ! filtering all elevations

! sub-sampling ! all elevations already filled: only one sub-sampling to do 
! -------------------------------------------------------------------------
    iobs1 = 0
    ifreq = nint(HODIM%Sample/HODIM%Resolution)
    write(batout,*) "Thinning every ",HODIM%Sample," km (each ",ifreq," pixels)"

    if (.not.allocated(zthmask)) allocate(zthmask(Radar%NPoints,Radar%NRayons),STAT=Error)
    if (Error /= 0) call Abor1("** ERROR - BATOR : cannot allocate zthmask(:,:)")
    zthmask=.FALSE.
    do i=1, 2*Radar%NPoints
      do j=1, 2*Radar%NPoints
        if (mod(j,ifreq) /= 0 .or. mod(i,ifreq) /= 0 ) cycle
        zdx  = i-(Radar%NPoints+0.5)
        zdy  = j-(Radar%NPoints+0.5)
        icol = nint( sqrt( zdx**2 + zdy**2 ) + 0.5 )
        if (icol > radar%NPoints) cycle
        zphi = atan2(zdy,zdx) - RPI/2.
        if (zphi < 0.) zphi = zphi + 2.*RPI
        irow = nint(0.5 + radar%NRayons * zphi/(2.*RPI))
        if (irow > Radar%NRayons) cycle
        zthmask(icol,irow) = .TRUE.
      enddo
    enddo
    write(batout,*) "Number of obs left after sub-sampling ",count(zthmask)
    do irow=1, Radar%NRayons
      do icol=1, Radar%NPoints
        iobs1 = icol + (irow-1)*Radar%NPoints
        if (.not. zthmask(icol,irow)) ztent(iobs1,ncinlv) = 0
      enddo
    enddo
! all elevations filled 

! compactage
    kobs0 = kobs
    kw0   = kw
    do i= 1, nb_obs
      if (ztent(i,ncinlv) > 0) then
        kobs = kobs + 1
        zent(kobs,:) = ztent(i,:)
!     zent(kobs,NCINLV) = ztent(i,NCINLV)
!     zent(kobs,NCILAT) = ztent(i,NCILAT)
!     zent(kobs,NCILON) = ztent(i,NCILON)
!     zent(kobs,NCIOTP) = ztent(i,NCIOTP)
!     zent(kobs,NCIRFL) = ztent(i,NCIRFL)
!     zent(kobs,NCIOCH) = ztent(i,NCIOCH)
!     zent(kobs,NCIDAT) = ztent(i,NCIDAT)
!     zent(kobs,NCIETM) = ztent(i,NCIETM)
        write(clsid(kobs),'(3x,A5)') Radar%Nod ! indicatif NOD

        inlv = 0
        do j=1, NbSelectedElangles * 3
          iw = (j-1)+3*(NbSelectedElangles*(i-1))+1
          if (ztwag(iw,4) /= rabso) then
            kw           = kw+1
            zwagon(kw,:) = ztwag(iw,:)
!         zwagon(kw,1) = ztwag(iw,1)
!         zwagon(kw,2) = ztwag(iw,2)
!         zwagon(kw,4) = ztwag(iw,4)
!         zwagon(kw,5) = ztwag(iw,5)
!         zwagon(kw,9) = ztwag(iw,9)
            inlv         = inlv+1
          endif
        enddo
        if (ztent(i,ncinlv) /= inlv) write(batout,*) "!!!!!! iobs:",i," NCINLV:",ztent(i,ncinlv)," &
                                                   & inb:",inlv,"ztwag(iw,4)",ztwag(iw,4)
      endif
    enddo
! decompte final des obs et des wagons

!WRITE(batout,*) ,"la date de la derniere observation &
!&                 est:",zent(kobs,NCIDAT)
!WRITE(batout,*) ,"l heure de la derniere observation &
!&                 est:",zent(kobs,NCIETM)
!WRITE(batout,*) ,"inlv vaut:",inlv,"pour la derniere obs:",kobs,"et la &
!          &         valeur du dernier wagon vaut:",zwagon(kw,4),"correspondant &
!          &         au varno:",zwagon(kw,1)  
!IF (zwagon(kw,1) == 29) WRITE(batout,*),"pour le varno precedent:",zwagon(kw-1,1),zwagon(kw-1,4)

!print *,"la valeur de kobs-kobs0 vaut:",kobs-kobs0
    print *,"la valeur de iobs vaut:",iobs
!print *,"la valeur de kw-kw0 vaut:",kw-kw0
!print *,"la valeur de kobs0 vaut:",kobs0
  endif

  write(batout,'(6x,"Selected Obs = ",I9,"  --> ",I12," datas.")')    Kobs-Kobs0, Kw-Kw0
  write(batout,'("Total selected Obs = ",I9,"  --> ",I12," datas.")') Kobs,       Kw





!free memory
  do i=1, size(FullDatasetList)
    if (associated(FullDatasetList(i)%GData))    deallocate(FullDatasetList(i)%GData)
    if (associated(FullDatasetList(i)%GQuality)) deallocate(FullDatasetList(i)%GQuality)
  enddo
  do i=1, size(SelectedElangles)
    if (associated(SelectedElangles(i)%DBZH)) nullify(SelectedElangles(i)%DBZH)
    if (associated(SelectedElangles(i)%TH))   nullify(SelectedElangles(i)%TH)
    if (associated(SelectedElangles(i)%VRAD)) nullify(SelectedElangles(i)%VRAD)
    if (associated(SelectedElangles(i)%FLAG)) nullify(SelectedElangles(i)%FLAG)
  enddo
  do i=1,size(Radar%FinalElev)
    if (associated(Radar%FinalElev(i)%DBZH)) deallocate(Radar%FinalElev(i)%DBZH%Values)
    nullify(Radar%FinalElev(i)%DBZH)
    if (associated(Radar%FinalElev(i)%TH))   deallocate(Radar%FinalElev(i)%TH%Values)
    nullify(Radar%FinalElev(i)%TH)
    if (associated(Radar%FinalElev(i)%VRAD)) deallocate(Radar%FinalElev(i)%VRAD%Values)
    nullify(Radar%FinalElev(i)%VRAD)
    if (associated(Radar%FinalElev(i)%FLAG)) deallocate(Radar%FinalElev(i)%FLAG%Values)
    nullify(Radar%FinalElev(i)%FLAG)
  enddo
  if (allocated(FullDatasetList))  deallocate(FullDatasetList)
  if (associated(Radar%FinalElev)) deallocate(Radar%FinalElev)
  if (allocated(NRaysPopulation))  deallocate(NRaysPopulation)
  if (allocated(SelectedElangles)) deallocate(SelectedElangles)
  if (allocated(ztent))   deallocate (ztent)
  if (allocated(ztwag))   deallocate (ztwag)
  if (allocated(zazim))   deallocate (zazim)
  if (allocated(zdist))   deallocate(zdist)
  if (allocated(zthmask)) deallocate(zthmask)

  if (lhook) call dr_hook('Odim',1,zhook_handle)

CONTAINS

SUBROUTINE GetDAttributes(IdIn,NomMembre,Attributes)
  integer(HID_T),intent(in)             :: IdIn
  character(len=*),intent(in)           :: NomMembre
  type(GAttributes),intent(inout)       :: Attributes
  integer(kind=HID_T)                   :: GroupId
  integer(kind=jpim)                    :: Error, Iret
  integer(kind=jpim)                    :: WaitedAtomic, WaitedRank
  real(kind=jprb)                       :: zhook_handle

  if (lhook) call dr_hook('GetDAttributes',0,zhook_handle)

  call h5gopen_f(IdIn,trim(NomMembre)//"/"//trim(HODIM%GrpWhatName),GroupId,Error)
  if (Error == 0) then
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(0)
    call GetAttribute(GroupId,HODIM%StartDateName,WaitedRank,WaitedAtomic,Iret)
    if (Iret == 0) then
      Attributes%StartDate = String1Buf(1)
      deallocate(String1Buf)
    endif
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(0)
    call GetAttribute(GroupId,HODIM%StartTimeName,WaitedRank,WaitedAtomic,Iret)
    if (Iret == 0) then
      Attributes%StartTime = String1Buf(1)
      deallocate(String1Buf)
    endif
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(0)
    call GetAttribute(GroupId,HODIM%QuantityName,WaitedRank,WaitedAtomic,Iret)
    if (Iret == 0) then 
      Attributes%Quantity  = String1Buf(1)
      deallocate(String1Buf)
    endif
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(2)
    call GetAttribute(GroupId,HODIM%GainName,WaitedRank,WaitedAtomic,Iret)
    if (Iret == 0) then 
      Attributes%Gain = Real1Buf(1)
      deallocate(real1buf)
    endif
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(2)
    call GetAttribute(groupId,HODIM%OffsetName,WaitedRank,WaitedAtomic,Iret)
    if (Iret == 0) then
      Attributes%Offset = Real1Buf(1)
      deallocate(Real1Buf)
    endif
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(2)
    call GetAttribute(groupId,HODIM%NoDataName,WaitedRank,WaitedAtomic,Iret)
    if (Iret == 0) then
      Attributes%NoData = Real1Buf(1)
      deallocate(Real1Buf)
    endif
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(2)
    call GetAttribute(groupId,HODIM%NoDetectName,WaitedRank,WaitedAtomic,Iret)
    if (Iret == 0) then
      Attributes%NoDetect = Real1Buf(1)
      deallocate(Real1Buf)
    endif
    call h5gclose_f(GroupId,Error)
  endif

  call h5gopen_f(IdIn,trim(NomMembre)//"/"//trim(HODIM%GrpWhereName),GroupId,Error)
  if (Error == 0) then
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(1)
    call GetAttribute(GroupId,HODIM%NraysName,WaitedRank,WaitedAtomic,Iret)
    if (Iret == 0) then
      Attributes%NRays = Inte1Buf(1)
      deallocate(Inte1Buf)
    endif
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(1)
    call GetAttribute(GroupId,HODIM%NbinsName,WaitedRank,WaitedAtomic,Iret)
    if (Iret == 0) then
      Attributes%Nbins = Inte1Buf(1)
      deallocate(Inte1Buf)
    endif
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(2)
    call GetAttribute(GroupId,HODIM%ElevName,WaitedRank,WaitedAtomic,Iret)
    if (Iret == 0) then 
      Attributes%Elangle  = Real1Buf(1)
      deallocate(Real1Buf)
    endif
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(2)
    call GetAttribute(GroupId,HODIM%RscaleName,WaitedRank,WaitedAtomic,Iret)
    if (Iret == 0) then 
      Attributes%RScale = Real1Buf(1)
      deallocate(real1buf)
    endif
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(2)
    call GetAttribute(groupId,HODIM%RStartName,WaitedRank,WaitedAtomic,Iret)
    if (Iret == 0) then
      Attributes%RStart = Real1Buf(1)
      deallocate(Real1Buf)
    endif
    call h5gclose_f(GroupId,Error)
  endif

  if (lhook) call dr_hook('GetDAttributes',1,zhook_handle)
END SUBROUTINE GetDAttributes

SUBROUTINE GetFAttributes(IdIn,NomMembre,Attributes)
  integer(HID_T),intent(in)             :: IdIn
  character(len=*),intent(in)           :: NomMembre
  type(GAttributes),intent(inout)       :: Attributes
  integer(kind=HID_T)                   :: GroupId
  integer(kind=jpim)                    :: Error, Iret
  integer(kind=jpim)                    :: WaitedAtomic, WaitedRank
  real(kind=jprb)                       :: zhook_handle

  if (lhook) call dr_hook('GetFAttributes',0,zhook_handle)
  
  call h5gopen_f(IdIn,trim(NomMembre)//"/"//trim(HODIM%GrpWhatName),GroupId,Error)
  if (Error == 0) then
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(2)
    call GetAttribute(GroupId,HODIM%GainName,WaitedRank,WaitedAtomic,Iret)
    if (Iret == 0) then 
      Attributes%Gain = Real1Buf(1)
      deallocate(real1buf)
    endif
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(2)
    call GetAttribute(groupId,HODIM%OffsetName,WaitedRank,WaitedAtomic,Iret)
    if (Iret == 0) then
      Attributes%Offset = Real1Buf(1)
      deallocate(Real1Buf)
    endif
    call h5gclose_f(GroupId,Error)
  endif

  call h5gopen_f(IdIn,trim(NomMembre)//"/"//trim(HODIM%GrpHowName),GroupId,Error)
  if (Error == 0) then
    WaitedRank   = 0
    WaitedAtomic = Hdf5Type(0)
    call GetAttribute(GroupId,HODIM%TaskName,WaitedRank,WaitedAtomic,Iret)
    if (Iret == 0) then 
      Attributes%Task  = String1Buf(1)
      deallocate(String1Buf)
    endif
    call h5gclose_f(GroupId,Error)
  endif

  if (lhook) call dr_hook('GetFAttributes',1,zhook_handle)
END SUBROUTINE GetFAttributes

SUBROUTINE BuildAttributesList(RootAttributes,DatasetAttributes,DataAttributes)
  type(GAttributes),intent(in)          :: RootAttributes
  type(GAttributes),intent(inout)       :: DatasetAttributes
  type(GAttributes),intent(inout)       :: DataAttributes
  real(kind=jprb)                       :: zhook_handle

  if (lhook) call dr_hook('BuildAttributesList',0,zhook_handle)
  
  if (trim(DataAttributes%StartDate) == '?') then
    if (trim(DatasetAttributes%StartDate) == '?') then
      DataAttributes%StartDate = RootAttributes%Startdate
    else
      DataAttributes%StartDate = DatasetAttributes%StartDate
    endif
  endif
  if (trim(DataAttributes%StartTime) == '?') then
    if (trim(DatasetAttributes%StartTime) == '?') then
      DataAttributes%StartTime = RootAttributes%StartTime
    else
      DataAttributes%StartTime = DatasetAttributes%StartTime
    endif
  endif
  if (trim(DataAttributes%Quantity) == '?') then
    if (trim(DatasetAttributes%Quantity) == '?') then
      DataAttributes%Quantity = RootAttributes%Quantity
    else
      DataAttributes%Quantity = DatasetAttributes%Quantity
    endif
  endif
  if (DataAttributes%Gain == rabsi) then
    if (DatasetAttributes%Gain == rabsi) then
      DataAttributes%Gain = RootAttributes%Gain
      if (DataAttributes%Gain == rabsi) DataAttributes%Gain = 1
    else
      DataAttributes%Gain = DatasetAttributes%Gain
    endif
  endif
  if (DataAttributes%Offset == rabsi) then
    if (DatasetAttributes%Offset == rabsi) then
      DataAttributes%Offset = RootAttributes%Offset
      if (DataAttributes%Offset == rabsi) DataAttributes%Offset = 0
    else
      DataAttributes%Offset = DatasetAttributes%Offset
    endif
  endif
  if (DataAttributes%NoData == rabsi) then
    if (DatasetAttributes%NoData == rabsi) then
      DataAttributes%NoData = RootAttributes%NoData
    else
      DataAttributes%NoData = DatasetAttributes%NoData
    endif
  endif
  if (DataAttributes%NoDetect == rabsi) then
    if (DatasetAttributes%NoDetect == rabsi) then
      DataAttributes%NoDetect = RootAttributes%NoDetect
    else
      DataAttributes%NoDetect = DatasetAttributes%NoDetect
    endif
  endif
  if (DataAttributes%NRays == nabsi) then
    if (DatasetAttributes%NRays == nabsi) then
      DataAttributes%NRays = RootAttributes%NRays
    else
      DataAttributes%NRays = DatasetAttributes%NRays
    endif
  endif
  if (DataAttributes%Nbins == nabsi) then
    if (DatasetAttributes%Nbins == nabsi) then
      DataAttributes%Nbins = RootAttributes%Nbins
    else
      DataAttributes%Nbins = DatasetAttributes%Nbins
    endif
  endif
  if (DataAttributes%Elangle == rabsi) then
    if (DatasetAttributes%Elangle == rabsi) then
      DataAttributes%Elangle = RootAttributes%Elangle
    else
      DataAttributes%Elangle = DatasetAttributes%Elangle
    endif
  endif
  if (DataAttributes%RScale == rabsi) then
    if (DatasetAttributes%RScale == rabsi) then
      DataAttributes%RScale = RootAttributes%RScale
    else
      DataAttributes%RScale = DatasetAttributes%RScale
    endif
  endif
  if (DataAttributes%RStart == rabsi) then
    if (DatasetAttributes%RStart == rabsi) then
      DataAttributes%RStart = RootAttributes%RStart
    else
      DataAttributes%RStart = DatasetAttributes%RStart
    endif
  endif
  if (trim(DataAttributes%Task) == '?') then
    if (trim(DatasetAttributes%Task) == '?') then
      DataAttributes%Task = RootAttributes%Task
    else
      DataAttributes%Task = DatasetAttributes%Task
    endif
  endif

  if (lhook) call dr_hook('BuildAttributesList',1,zhook_handle)
END SUBROUTINE BuildAttributesList

END SUBROUTINE Odim

SUBROUTINE Maps(rscale,rstart,elangle,arrayin,elevation)
  real(kind=jprb),intent(in)                :: rscale, rstart, elangle
  real(kind=jprb),dimension(:,:),intent(in) :: arrayin
  real(kind=jprb),dimension(:,:),pointer    :: elevation
  integer(kind=jpim)                        :: nraysin, nraysout, nbinsin, ratio, i, j, x, y, indice
  real(kind=jprb)                           :: distance, CurrentDist, LesserDist, RadElev
  real(kind=jprb)                           :: zhook_handle

  if (lhook) call dr_hook('Maps',0,zhook_handle)

  nraysin     = size(arrayin,2)
  nraysout    = size(elevation,2)
  nbinsin     = size(arrayin,1)
  ratio       = nraysout / nraysin
  x           = 1
  y           = 0
  LesserDist  = 9999._jprb

#ifdef DEBUG
  write(Batout,'("rscale=",F5.0," rstart=",F5.0," NraysIn=",I4," NraysOut=",I4," NbinsIn=",I4)') rscale, rstart, &
                & nraysin, nraysout, nbinsin
#endif

  do i=1, nraysin
    do j=1, nbinsin
      distance = (j-1) * rscale + rstart
      RadElev = elangle*rpi/180._jprb
      if (RadElev > 0.35_jprb) distance = cos(RadElev)*distance
      indice = nint(distance/HODIM%Resolution)
      if (indice /= x) LesserDist = 9999._jprb
      CurrentDist = abs(HODIM%Resolution * indice - distance)
      x = indice
      y = 1+(i-1)*ratio
      if (CurrentDist <= (HODIM%Resolution/2._jprb) .and. CurrentDist < LesserDist .and. indice > 0 ) then
        LesserDist = CurrentDist
        elevation(x, y) = arrayin(j,i)
      endif
    enddo
  enddo

  if (lhook) call dr_hook('Maps',1,zhook_handle)
END SUBROUTINE Maps


END MODULE BATOR_DECODHDF5_MOD
