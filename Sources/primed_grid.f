!*******************************************************************************
!>  @file primed_grid.f
!>  @brief Contains module @ref primed_grid.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref primed_grid_class. This contains the
!>  state variables to define the currents and positions of the volumn integral.
!*******************************************************************************
      MODULE primed_grid
      USE stel_kinds, ONLY: rprec
      USE profiler
      USE bmw_parallel_context

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) primed grid base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a primed grid. This is grid the volume integral will
!>  be summed over.
!-------------------------------------------------------------------------------
      TYPE :: primed_grid_class
!>  Volume integration element.
         REAL (rprec)                            :: dvol
!>  Toroidal grid spacing.
         REAL (rprec)                            :: dv

!>  X position.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: x => null()
!>  Y position.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: y => null()
!>  Z position.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: z => null()

!>  Current density in the X direction.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: j_x => null()
!>  Current density in the Y direction.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: j_y => null()
!>  Current density in the Z direction.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: j_z => null()
      CONTAINS
         FINAL :: primed_grid_destruct
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the bmw_commandline_parser constructor.
!-------------------------------------------------------------------------------
      INTERFACE primed_grid_class
         MODULE PROCEDURE primed_grid_construct
      END INTERFACE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref primed_grid_class object.
!>
!>  Allocates memory and initializes a @ref primed_grid_class object depending
!>  on the option flags.
!>
!>  @param[in] num_v    Number of toroidal grid points.
!>  @param[in] flags    Number of toroidal grid points.
!>  @param[in] parallel @ref bmw_parallel_context_class object instance.
!>  @param[in] io_unit  Unit number to write messages to.
!>  @returns A pointer to a constructed @ref primed_grid_class object.
!-------------------------------------------------------------------------------
      FUNCTION primed_grid_construct(num_v, flags, siesta_file,                &
     &                               parallel, io_unit)
      USE bmw_state_flags

      IMPLICIT NONE

!  Declare Arguments
      TYPE (primed_grid_class), POINTER :: primed_grid_construct
      INTEGER, INTENT(in)                           :: num_v
      INTEGER, INTENT(in)                           :: flags
      CHARACTER (len=*), INTENT(in)                 :: siesta_file
      TYPE (bmw_parallel_context_class), INTENT(in) :: parallel
      INTEGER, INTENT(in)                           :: io_unit

!  local variables
      REAL (rprec)                                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      IF (BTEST(flags, bmw_state_flags_ju)) THEN
         primed_grid_construct =>                                              &
     &      primed_grid_construct_ju(num_v, parallel)
      ELSE IF (BTEST(flags, bmw_state_flags_jv)) THEN
         primed_grid_construct =>                                              &
     &      primed_grid_construct_jv(num_v, parallel)
      ELSE IF (BTEST(flags, bmw_state_flags_siesta)) THEN
         primed_grid_construct =>                                              &
     &      primed_grid_construct_siesta(num_v, siesta_file, parallel)
      ELSE
         primed_grid_construct =>                                              &
     &      primed_grid_construct_both(num_v, parallel)
      END IF

      IF (parallel%offset .eq. 0) THEN
         WRITE(io_unit,1000)
      END IF

      CALL profiler_set_stop_time('primed_grid_construct', start_time)

1000  FORMAT('Prime Grid Ready')

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref primed_grid_class object.
!>
!>  Allocates memory and initializes a @ref primed_grid_class object. This
!>  computes the currents and positions on the primed grid. Plasma currents are
!>  obtained from Curl(B).
!>
!>  @param[in] num_v Number of toroidal grid points.
!>  @param[in] parallel @ref bmw_parallel_context_class object instance.
!>  @returns A pointer to a constructed @ref primed_grid_class object.
!-------------------------------------------------------------------------------
      FUNCTION primed_grid_construct_both(num_v, parallel)
      USE stel_constants, ONLY: twopi, mu0
      USE read_wout_mod, ONLY: mnmax, mnmax_nyq, lasym, isigng, ns,            &
     &                         xm, xn, xm_nyq, xn_nyq,                         &
     &                         rmnc, zmns, currumnc, currvmnc,                 &
     &                         rmns, zmnc, currumns, currvmns!, gmnc   ! TEMP ,gmnc

      IMPLICIT NONE

!  Declare Arguments
      CLASS (primed_grid_class), POINTER :: primed_grid_construct_both
      INTEGER, INTENT(in)                            :: num_v
      CLASS (bmw_parallel_context_class), INTENT(in) :: parallel

!  local variables
      REAL (rprec)                                   :: start_time
      INTEGER                                        :: i
      REAL (rprec)                                   :: x
      REAL (rprec)                                   :: ds
      INTEGER                                        :: si
      INTEGER                                        :: ui
      INTEGER                                        :: vi
      REAL (rprec)                                   :: r
      REAL (rprec)                                   :: z
      REAL (rprec)                                   :: ru
      REAL (rprec)                                   :: zu
      REAL (rprec)                                   :: rv
      REAL (rprec)                                   :: zv
      REAL (rprec)                                   :: ju
      REAL (rprec)                                   :: jv
      REAL (rprec)                                   :: jr
      REAL (rprec)                                   :: jp
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: cosv
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: sinv
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: cosmu
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: sinmu
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: cosmu_nyq
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: sinmu_nyq
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: cosnv
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: sinnv
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: cosnv_nyq
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: sinnv_nyq
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: cosmn
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: sinmn
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: cosmn_nyq
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: sinmn_nyq
      INTEGER                                        :: num_v_use

!  local parameters
      INTEGER, PARAMETER                             :: num_u = 101
      REAL (rprec), PARAMETER                        :: du = twopi/num_u

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(primed_grid_construct_both)

      IF (num_v .eq. 1) THEN
         num_v_use = 360
      ELSE
         num_v_use = num_v
      END IF

      ALLOCATE(primed_grid_construct_both%x(ns,num_u,num_v_use))
      ALLOCATE(primed_grid_construct_both%y(ns,num_u,num_v_use))
      ALLOCATE(primed_grid_construct_both%z(ns,num_u,num_v_use))

      ALLOCATE(primed_grid_construct_both%j_x(ns,num_u,num_v_use))
      ALLOCATE(primed_grid_construct_both%j_y(ns,num_u,num_v_use))
      ALLOCATE(primed_grid_construct_both%j_z(ns,num_u,num_v_use))

      ALLOCATE(cosv(num_v_use))
      ALLOCATE(sinv(num_v_use))
      ALLOCATE(cosmu(mnmax,num_u))
      ALLOCATE(sinmu(mnmax,num_u))
      ALLOCATE(cosnv(mnmax,num_v_use))
      ALLOCATE(sinnv(mnmax,num_v_use))
      ALLOCATE(cosmu_nyq(mnmax_nyq,num_u))
      ALLOCATE(sinmu_nyq(mnmax_nyq,num_u))
      ALLOCATE(cosnv_nyq(mnmax_nyq,num_v_use))
      ALLOCATE(sinnv_nyq(mnmax_nyq,num_v_use))

      ds = 1.0/(ns - 1.0)
      primed_grid_construct_both%dv = twopi/num_v_use

      primed_grid_construct_both%dvol = isigng*ds*du                           &
     &   * primed_grid_construct_both%dv/(2.0*twopi)

!$OMP PARALLEL
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(i, x, si, ui, vi, r, z, ru, rv, zu, zv, ju, jv, jr, jp,         &
!$OMP&         cosmn, sinmn, cosmn_nyq, sinmn_nyq)

!  Multi process will do an all reduce so these arrays need to be initalized.
      IF (parallel%stride .gt. 1) THEN
!$OMP WORKSHARE
         primed_grid_construct_both%x = 0.0
         primed_grid_construct_both%y = 0.0
         primed_grid_construct_both%z = 0.0
         primed_grid_construct_both%j_x = 0.0
         primed_grid_construct_both%j_y = 0.0
         primed_grid_construct_both%j_z = 0.0
         cosv = 0.0
         sinv = 0.0
         cosmu = 0.0
         sinmu = 0.0
         cosnv = 0.0
         sinnv = 0.0
         cosmu_nyq = 0.0
         sinmu_nyq = 0.0
         cosnv_nyq = 0.0
         sinnv_nyq = 0.0
!$OMP END WORKSHARE
      END IF

!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO i = parallel%start(num_u), parallel%end(num_u)
         x = (i - 0.5)*du
         cosmu(:,i) = COS(xm*x)
         sinmu(:,i) = SIN(xm*x)
         cosmu_nyq(:,i) = COS(xm_nyq*x)
         sinmu_nyq(:,i) = SIN(xm_nyq*x)
      END DO
!$OMP END DO

!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO i = parallel%start(num_v_use), parallel%end(num_v_use)
         x = (i - 0.5)*primed_grid_construct_both%dv
         cosv(i) = COS(x)
         sinv(i) = SIN(x)
         cosnv(:,i) = COS(xn*x)
         sinnv(:,i) = SIN(xn*x)
         cosnv_nyq(:,i) = COS(xn_nyq*x)
         sinnv_nyq(:,i) = SIN(xn_nyq*x)
      END DO
!$OMP END DO

      IF (parallel%stride .gt. 1) THEN
!$OMP SINGLE
         CALL parallel%reduce(cosv)
         CALL parallel%reduce(sinv)
         CALL parallel%reduce(cosmu)
         CALL parallel%reduce(sinmu)
         CALL parallel%reduce(cosmu_nyq)
         CALL parallel%reduce(sinmu_nyq)
         CALL parallel%reduce(cosnv)
         CALL parallel%reduce(sinnv)
         CALL parallel%reduce(cosnv_nyq)
         CALL parallel%reduce(sinnv_nyq)
!$OMP END SINGLE
      END IF

      ALLOCATE(cosmn(mnmax))
      ALLOCATE(sinmn(mnmax))
      ALLOCATE(cosmn_nyq(mnmax_nyq))
      IF (lasym) THEN
         ALLOCATE(sinmn_nyq(mnmax_nyq))
      END IF

!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO i = parallel%start(ns*num_u*num_v_use),                               &
     &       parallel%end(ns*num_u*num_v_use)
         si = bmw_parallel_context_i(i, ns)
         ui = bmw_parallel_context_j(i, ns, num_u)
         vi = bmw_parallel_context_k(i, ns, num_u)

         cosmn = cosmu(:,ui)*cosnv(:,vi) + sinmu(:,ui)*sinnv(:,vi)
         sinmn = sinmu(:,ui)*cosnv(:,vi) - cosmu(:,ui)*sinnv(:,vi)
         cosmn_nyq = cosmu_nyq(:,ui)*cosnv_nyq(:,vi)                           &
     &             + sinmu_nyq(:,ui)*sinnv_nyq(:,vi)

         r = SUM(rmnc(:,si)*cosmn(:))
         z = SUM(zmns(:,si)*sinmn(:))

         ru = -SUM(xm*rmnc(:,si)*sinmn(:))
         rv =  SUM(xn*rmnc(:,si)*sinmn(:))
         zu =  SUM(xm*zmns(:,si)*cosmn(:))
         zv = -SUM(xn*zmns(:,si)*cosmn(:))

         ju = SUM(currumnc(:,si)*cosmn_nyq(:))*mu0
         jv = SUM(currvmnc(:,si)*cosmn_nyq(:))*mu0

         IF (lasym) THEN
            sinmn_nyq = sinmu_nyq(:,ui)*cosnv_nyq(:,vi)                        &
     &                - cosmu_nyq(:,ui)*sinnv_nyq(:,vi)

            r = r + SUM(rmns(:,si)*sinmn(:))
            z = z + SUM(zmnc(:,si)*cosmn(:))

            ru = ru + SUM(xm*rmns(:,si)*cosmn(:))
            rv = rv - SUM(xn*rmns(:,si)*cosmn(:))
            zu = zu - SUM(xm*zmnc(:,si)*sinmn(:))
            zv = zv + SUM(xn*zmnc(:,si)*sinmn(:))

            ju = ju + SUM(currumns(:,si)*sinmn_nyq(:))*mu0
            jv = jv + SUM(currvmns(:,si)*sinmn_nyq(:))*mu0
         END IF

         jr = ju*ru + jv*rv
         jp = jv*r
         primed_grid_construct_both%j_z(si,ui,vi) = ju*zu + jv*zv

         IF (si .eq. 1 .or. si .eq. ns) THEN
            jr = jr/2.0
            jp = jp/2.0
            primed_grid_construct_both%j_z(si,ui,vi) =                         &
     &         primed_grid_construct_both%j_z(si,ui,vi)/2.0
         END IF

         primed_grid_construct_both%j_x(si,ui,vi) = jr*cosv(vi)                &
     &                                            - jp*sinv(vi)
         primed_grid_construct_both%j_y(si,ui,vi) = jr*sinv(vi)                &
     &                                            + jp*cosv(vi)

         primed_grid_construct_both%x(si,ui,vi) = r*cosv(vi)
         primed_grid_construct_both%y(si,ui,vi) = r*sinv(vi)
         primed_grid_construct_both%z(si,ui,vi) = z
      END DO
!$OMP END DO

      DEALLOCATE(cosmn)
      DEALLOCATE(sinmn)
      DEALLOCATE(cosmn_nyq)
      IF (lasym) THEN
         DEALLOCATE(sinmn_nyq)
      END IF
!$OMP END PARALLEL

      DEALLOCATE(cosv)
      DEALLOCATE(sinv)
      DEALLOCATE(cosmu)
      DEALLOCATE(sinmu)
      DEALLOCATE(cosnv)
      DEALLOCATE(sinnv)
      DEALLOCATE(cosmu_nyq)
      DEALLOCATE(sinmu_nyq)
      DEALLOCATE(cosnv_nyq)
      DEALLOCATE(sinnv_nyq)

!  Multi process did not fill out the entire array. Get the missing pieces from
!  the other processes.
      IF (parallel%stride .gt. 1) THEN
         CALL parallel%reduce(primed_grid_construct_both%x)
         CALL parallel%reduce(primed_grid_construct_both%y)
         CALL parallel%reduce(primed_grid_construct_both%z)
         CALL parallel%reduce(primed_grid_construct_both%j_x)
         CALL parallel%reduce(primed_grid_construct_both%j_y)
         CALL parallel%reduce(primed_grid_construct_both%j_z)
      END IF

      CALL profiler_set_stop_time('primed_grid_construct_both',                &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref primed_grid_class object.
!>
!>  Allocates memory and initializes a @ref primed_grid_class object with an
!>  This computes the currents and positions on the primed grid. J^u plasma
!>  current is obtained from force balance.
!>
!>  J^u = (p' + J^v*B^u)/B^v
!>
!>  @param[in] num_v Number of toroidal grid points.
!>  @param[in] parallel @ref bmw_parallel_context_class object instance.
!>  @returns A pointer to a constructed @ref primed_grid_class object.
!-------------------------------------------------------------------------------
      FUNCTION primed_grid_construct_ju(num_v, parallel)
      USE stel_constants, ONLY: twopi, mu0
      USE read_wout_mod, ONLY: mnmax, mnmax_nyq, lasym, isigng, ns,            &
     &                         rmnc, rmns, currvmnc,                           &
     &                         zmnc, zmns, currvmns,                           &
     &                         bsupumnc, bsupvmnc, bsupumns, bsupvmns,         &
     &                         xm, xn, xm_nyq, xn_nyq, presf

      IMPLICIT NONE

!  Declare Arguments
      CLASS (primed_grid_class), POINTER :: primed_grid_construct_ju
      INTEGER, INTENT(in)                            :: num_v
      CLASS (bmw_parallel_context_class), INTENT(in) :: parallel

!  local variables
      REAL (rprec)                                   :: start_time
      INTEGER                                        :: i
      REAL (rprec)                                   :: x
      REAL (rprec)                                   :: ds
      INTEGER                                        :: si
      INTEGER                                        :: ui
      INTEGER                                        :: vi
      REAL (rprec)                                   :: r
      REAL (rprec)                                   :: z
      REAL (rprec)                                   :: ru
      REAL (rprec)                                   :: zu
      REAL (rprec)                                   :: rv
      REAL (rprec)                                   :: zv
      REAL (rprec)                                   :: ju
      REAL (rprec)                                   :: jv
      REAL (rprec)                                   :: bu
      REAL (rprec)                                   :: bv
      REAL (rprec)                                   :: jr
      REAL (rprec)                                   :: jp
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: cosv
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: sinv
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: cosmu
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: sinmu
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: cosmu_nyq
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: sinmu_nyq
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: cosnv
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: sinnv
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: cosnv_nyq
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: sinnv_nyq
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: cosmn
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: sinmn
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: cosmn_nyq
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: sinmn_nyq
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: rmnch
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: rmnsh
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: zmnch
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: zmnsh
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: currvmnch
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: currvmnsh
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: p_prime
      INTEGER                                        :: num_v_use

!  local parameters
      INTEGER, PARAMETER                             :: num_u = 101
      REAL (rprec), PARAMETER                        :: du = twopi/num_u

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(primed_grid_construct_ju)

      IF (num_v .eq. 1) THEN
         num_v_use = 360
      ELSE
         num_v_use = num_v
      END IF

      ALLOCATE(primed_grid_construct_ju%x(ns - 1,num_u,num_v_use))
      ALLOCATE(primed_grid_construct_ju%y(ns - 1,num_u,num_v_use))
      ALLOCATE(primed_grid_construct_ju%z(ns - 1,num_u,num_v_use))

      ALLOCATE(primed_grid_construct_ju%j_x(ns - 1,num_u,num_v_use))
      ALLOCATE(primed_grid_construct_ju%j_y(ns - 1,num_u,num_v_use))
      ALLOCATE(primed_grid_construct_ju%j_z(ns - 1,num_u,num_v_use))

      ALLOCATE(cosv(num_v_use))
      ALLOCATE(sinv(num_v_use))
      ALLOCATE(cosmu(mnmax,num_u))
      ALLOCATE(sinmu(mnmax,num_u))
      ALLOCATE(cosnv(mnmax,num_v_use))
      ALLOCATE(sinnv(mnmax,num_v_use))
      ALLOCATE(cosmu_nyq(mnmax_nyq,num_u))
      ALLOCATE(sinmu_nyq(mnmax_nyq,num_u))
      ALLOCATE(cosnv_nyq(mnmax_nyq,num_v_use))
      ALLOCATE(sinnv_nyq(mnmax_nyq,num_v_use))

      ALLOCATE(rmnch(mnmax,ns - 1))
      ALLOCATE(zmnsh(mnmax,ns - 1))
      ALLOCATE(currvmnch(mnmax_nyq,ns - 1))

      IF (lasym) THEN
         ALLOCATE(rmnsh(mnmax,ns - 1))
         ALLOCATE(zmnch(mnmax,ns - 1))
         ALLOCATE(currvmnsh(mnmax_nyq,ns - 1))
      END IF

      ALLOCATE(p_prime(ns - 1))

      ds = 1.0/(ns - 1.0)
      primed_grid_construct_ju%dv = twopi/num_v_use

      primed_grid_construct_ju%dvol = isigng*ds*du                             &
     &   * primed_grid_construct_ju%dv/(2.0*twopi)

!$OMP PARALLEL
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(i, x, si, ui, vi, r, z, ru, rv, zu, zv, ju, jv, bu, bv,         &
!$OMP&         cosmn, sinmn, cosmn_nyq, sinmn_nyq, jr, jp)

!  Multi process will do an all reduce so these arrays need to be initalized.
      IF (parallel%stride .gt. 1) THEN
!$OMP WORKSHARE
         primed_grid_construct_ju%x = 0.0
         primed_grid_construct_ju%y = 0.0
         primed_grid_construct_ju%z = 0.0
         primed_grid_construct_ju%j_x = 0.0
         primed_grid_construct_ju%j_y = 0.0
         primed_grid_construct_ju%j_z = 0.0
         cosv = 0.0
         sinv = 0.0
         cosmu = 0.0
         sinmu = 0.0
         cosnv = 0.0
         sinnv = 0.0
         cosmu_nyq = 0.0
         sinmu_nyq = 0.0
         cosnv_nyq = 0.0
         sinnv_nyq = 0.0
         rmnch = 0.0
         zmnsh = 0.0
         currvmnch = 0.0
         p_prime = 0.0
!$OMP END WORKSHARE
         IF (lasym) THEN
!$OMP WORKSHARE
            rmnsh = 0.0
            zmnch = 0.0
            currvmnsh = 0.0
!$OMP END WORKSHARE
         END IF
      END IF

!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO i = parallel%start(ns - 1), parallel%end(ns - 1)
         rmnch(:,i) = (rmnc(:,i + 1) + rmnc(:,i))/2.0
         zmnsh(:,i) = (zmns(:,i + 1) + zmns(:,i))/2.0

         currvmnch(:,i) = (currvmnc(:,i + 1) + currvmnc(:,i))/2.0

         p_prime(i) = (presf(i + 1) - presf(i))/ds

         IF (lasym) THEN
            rmnsh(:,i) = (rmns(:,i + 1) + rmns(:,i))/2.0
            zmnch(:,i) = (zmnc(:,i + 1) + zmnc(:,i))/2.0

            currvmnsh(:,i) = (currvmns(:,i + 1) + currvmns(:,i))/2.0
         END IF
      END DO
!$OMP END DO

!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO i = parallel%start(num_u), parallel%end(num_u)
         x = (i - 0.5)*du
         cosmu(:,i) = COS(xm*x)
         sinmu(:,i) = SIN(xm*x)
         cosmu_nyq(:,i) = COS(xm_nyq*x)
         sinmu_nyq(:,i) = SIN(xm_nyq*x)
      END DO
!$OMP END DO

!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO i = parallel%start(num_v_use), parallel%end(num_v_use)
         x = (i - 0.5)*primed_grid_construct_ju%dv
         cosv(i) = COS(x)
         sinv(i) = SIN(x)
         cosnv(:,i) = COS(xn*x)
         sinnv(:,i) = SIN(xn*x)
         cosnv_nyq(:,i) = COS(xn_nyq*x)
         sinnv_nyq(:,i) = SIN(xn_nyq*x)
      END DO
!$OMP END DO

      IF (parallel%stride .gt. 1) THEN
!$OMP SINGLE
         CALL parallel%reduce(cosv)
         CALL parallel%reduce(sinv)
         CALL parallel%reduce(cosmu)
         CALL parallel%reduce(sinmu)
         CALL parallel%reduce(cosmu_nyq)
         CALL parallel%reduce(sinmu_nyq)
         CALL parallel%reduce(cosnv)
         CALL parallel%reduce(sinnv)
         CALL parallel%reduce(cosnv_nyq)
         CALL parallel%reduce(sinnv_nyq)
         CALL parallel%reduce(rmnch)
         CALL parallel%reduce(zmnsh)
         CALL parallel%reduce(currvmnch)
         CALL parallel%reduce(p_prime)
         IF (lasym) THEN
            CALL parallel%reduce(rmnsh)
            CALL parallel%reduce(zmnch)
            CALL parallel%reduce(currvmnsh)
         END IF
!$OMP END SINGLE
      END IF

      ALLOCATE(cosmn(mnmax))
      ALLOCATE(sinmn(mnmax))
      ALLOCATE(cosmn_nyq(mnmax_nyq))
      IF (lasym) THEN
         ALLOCATE(sinmn_nyq(mnmax_nyq))
      END IF

!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO i = parallel%start(ns*num_u*num_v_use),                               &
     &       parallel%end(ns*num_u*num_v_use)
         si = bmw_parallel_context_i(i, ns)
         ui = bmw_parallel_context_j(i, ns, num_u)
         vi = bmw_parallel_context_k(i, ns, num_u)

         cosmn = cosmu(:,ui)*cosnv(:,vi) + sinmu(:,ui)*sinnv(:,vi)
         sinmn = sinmu(:,ui)*cosnv(:,vi) - cosmu(:,ui)*sinnv(:,vi)
         cosmn_nyq = cosmu_nyq(:,ui)*cosnv_nyq(:,vi)                           &
     &             + sinmu_nyq(:,ui)*sinnv_nyq(:,vi)

         r = SUM(rmnch(:,si)*cosmn(:))
         z = SUM(zmnsh(:,si)*sinmn(:))

         ru = -SUM(xm*rmnch(:,si)*sinmn(:))
         rv =  SUM(xn*rmnch(:,si)*sinmn(:))
         zu =  SUM(xm*zmnsh(:,si)*cosmn(:))
         zv = -SUM(xn*zmnsh(:,si)*cosmn(:))

         jv = SUM(currvmnch(:,si)*cosmn_nyq(:))

         bu = SUM(bsupumnc(:,si + 1)*cosmn_nyq(:))
         bv = SUM(bsupvmnc(:,si + 1)*cosmn_nyq(:))

         IF (lasym) THEN
            sinmn_nyq = sinmu_nyq(:,ui)*cosnv_nyq(:,vi)                        &
     &                - cosmu_nyq(:,ui)*sinnv_nyq(:,vi)

            r = r + SUM(rmnsh(:,si)*sinmn(:))
            z = z + SUM(zmnch(:,si)*cosmn(:))

            ru = ru + SUM(xm*rmnsh(:,si)*cosmn(:))
            rv = rv - SUM(xn*rmnsh(:,si)*cosmn(:))
            zu = zu - SUM(xm*zmnch(:,si)*sinmn(:))
            zv = zv + SUM(xn*zmnch(:,si)*sinmn(:))

            jv = jv + SUM(currvmnsh(:,si)*sinmn_nyq(:))

            bu = bu + SUM(bsupumns(:,si + 1)*sinmn_nyq(:))
            bv = bv + SUM(bsupvmns(:,si + 1)*sinmn_nyq(:))
         END IF

         ju = (p_prime(si) + jv*bu)/bv*mu0
         jv = jv*mu0

         jr = ju*ru + jv*rv
         jp = jv*r

         primed_grid_construct_ju%j_z(si,ui,vi) = ju*zu + jv*zv
         primed_grid_construct_ju%j_x(si,ui,vi) = jr*cosv(vi)                  &
     &                                          - jp*sinv(vi)
         primed_grid_construct_ju%j_y(si,ui,vi) = jr*sinv(vi)                  &
     &                                          + jp*cosv(vi)

         primed_grid_construct_ju%x(si,ui,vi) = r*cosv(vi)
         primed_grid_construct_ju%y(si,ui,vi) = r*sinv(vi)
         primed_grid_construct_ju%z(si,ui,vi) = z
      END DO
!$OMP END DO

      DEALLOCATE(cosmn)
      DEALLOCATE(sinmn)
      DEALLOCATE(cosmn_nyq)
      IF (lasym) THEN
         DEALLOCATE(sinmn_nyq)
      END IF
!$OMP END PARALLEL

      DEALLOCATE(cosv)
      DEALLOCATE(sinv)
      DEALLOCATE(cosmu)
      DEALLOCATE(sinmu)
      DEALLOCATE(cosnv)
      DEALLOCATE(sinnv)
      DEALLOCATE(cosmu_nyq)
      DEALLOCATE(sinmu_nyq)
      DEALLOCATE(cosnv_nyq)
      DEALLOCATE(sinnv_nyq)

      DEALLOCATE(rmnch)
      DEALLOCATE(zmnsh)
      DEALLOCATE(currvmnch)
      DEALLOCATE(p_prime)
      IF (lasym) THEN
         DEALLOCATE(rmnsh)
         DEALLOCATE(zmnch)
         DEALLOCATE(currvmnsh)
      END IF

!  Multi process did not fill out the entire array. Get the missing pieces from
!  the other processes.
      IF (parallel%stride .gt. 1) THEN
         CALL parallel%reduce(primed_grid_construct_ju%x)
         CALL parallel%reduce(primed_grid_construct_ju%y)
         CALL parallel%reduce(primed_grid_construct_ju%z)
         CALL parallel%reduce(primed_grid_construct_ju%j_x)
         CALL parallel%reduce(primed_grid_construct_ju%j_y)
         CALL parallel%reduce(primed_grid_construct_ju%j_z)
      END IF

      CALL profiler_set_stop_time('primed_grid_construct_ju',                  &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref primed_grid_class object.
!>
!>  Allocates memory and initializes a @ref primed_grid_class object with an
!>  This computes the currents and positions on the primed grid. J^v plasma
!>  current is obtained from force balance.
!>
!>  J^v = (J^u*B^v - p')/B^u
!>
!>  @param[in] num_v Number of toroidal grid points.
!>  @param[in] parallel @ref bmw_parallel_context_class object instance.
!>  @returns A pointer to a constructed @ref primed_grid_class object.
!-------------------------------------------------------------------------------
      FUNCTION primed_grid_construct_jv(num_v, parallel)
      USE stel_constants, ONLY: twopi, mu0
      USE read_wout_mod, ONLY: mnmax, mnmax_nyq, lasym, isigng, ns,            &
     &                         rmnc, rmns, currumnc,                           &
     &                         zmnc, zmns, currumns,                           &
     &                         bsupumnc, bsupvmnc, bsupumns, bsupvmns,         &
     &                         xm, xn, xm_nyq, xn_nyq, presf

      IMPLICIT NONE

!  Declare Arguments
      CLASS (primed_grid_class), POINTER :: primed_grid_construct_jv
      INTEGER, INTENT(in)                            :: num_v
      CLASS (bmw_parallel_context_class), INTENT(in) :: parallel

!  local variables
      REAL (rprec)                                   :: start_time
      INTEGER                                        :: i
      REAL (rprec)                                   :: x
      REAL (rprec)                                   :: ds
      INTEGER                                        :: si
      INTEGER                                        :: ui
      INTEGER                                        :: vi
      REAL (rprec)                                   :: r
      REAL (rprec)                                   :: z
      REAL (rprec)                                   :: ru
      REAL (rprec)                                   :: zu
      REAL (rprec)                                   :: rv
      REAL (rprec)                                   :: zv
      REAL (rprec)                                   :: ju
      REAL (rprec)                                   :: jv
      REAL (rprec)                                   :: bu
      REAL (rprec)                                   :: bv
      REAL (rprec)                                   :: jr
      REAL (rprec)                                   :: jp
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: cosv
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: sinv
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: cosmu
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: sinmu
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: cosmu_nyq
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: sinmu_nyq
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: cosnv
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: sinnv
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: cosnv_nyq
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: sinnv_nyq
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: cosmn
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: sinmn
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: cosmn_nyq
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: sinmn_nyq
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: rmnch
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: rmnsh
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: zmnch
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: zmnsh
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: currumnch
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: currumnsh
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: p_prime
      INTEGER                                        :: num_v_use

!  local parameters
      INTEGER, PARAMETER                             :: num_u = 101
      REAL (rprec), PARAMETER                        :: du = twopi/num_u

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(primed_grid_construct_jv)

      IF (num_v .eq. 1) THEN
         num_v_use = 360
      ELSE
         num_v_use = num_v
      END IF

      ALLOCATE(primed_grid_construct_jv%x(ns - 1,num_u,num_v_use))
      ALLOCATE(primed_grid_construct_jv%y(ns - 1,num_u,num_v_use))
      ALLOCATE(primed_grid_construct_jv%z(ns - 1,num_u,num_v_use))

      ALLOCATE(primed_grid_construct_jv%j_x(ns - 1,num_u,num_v_use))
      ALLOCATE(primed_grid_construct_jv%j_y(ns - 1,num_u,num_v_use))
      ALLOCATE(primed_grid_construct_jv%j_z(ns - 1,num_u,num_v_use))

      ALLOCATE(cosv(num_v_use))
      ALLOCATE(sinv(num_v_use))
      ALLOCATE(cosmu(mnmax,num_u))
      ALLOCATE(sinmu(mnmax,num_u))
      ALLOCATE(cosnv(mnmax,num_v_use))
      ALLOCATE(sinnv(mnmax,num_v_use))
      ALLOCATE(cosmu_nyq(mnmax_nyq,num_u))
      ALLOCATE(sinmu_nyq(mnmax_nyq,num_u))
      ALLOCATE(cosnv_nyq(mnmax_nyq,num_v_use))
      ALLOCATE(sinnv_nyq(mnmax_nyq,num_v_use))

      ALLOCATE(rmnch(mnmax,ns - 1))
      ALLOCATE(zmnsh(mnmax,ns - 1))
      ALLOCATE(currumnch(mnmax_nyq,ns - 1))

      IF (lasym) THEN
         ALLOCATE(rmnsh(mnmax,ns - 1))
         ALLOCATE(zmnch(mnmax,ns - 1))
         ALLOCATE(currumnsh(mnmax_nyq,ns - 1))
      END IF

      ALLOCATE(p_prime(ns - 1))

      ds = 1.0/(ns - 1.0)
      primed_grid_construct_jv%dv = twopi/num_v_use

      primed_grid_construct_jv%dvol = isigng*ds*du                             &
     &   * primed_grid_construct_jv%dv/(2.0*twopi)

!$OMP PARALLEL
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(i, x, si, ui, vi, r, z, ru, rv, zu, zv, ju, jv, bu, bv,         &
!$OMP&         cosmn, sinmn, cosmn_nyq, sinmn_nyq, jr, jp)

!  Multi process will do an all reduce so these arrays need to be initalized.
      IF (parallel%stride .gt. 1) THEN
!$OMP WORKSHARE
         primed_grid_construct_jv%x = 0.0
         primed_grid_construct_jv%y = 0.0
         primed_grid_construct_jv%z = 0.0
         primed_grid_construct_jv%j_x = 0.0
         primed_grid_construct_jv%j_y = 0.0
         primed_grid_construct_jv%j_z = 0.0
         cosv = 0.0
         sinv = 0.0
         cosmu = 0.0
         sinmu = 0.0
         cosnv = 0.0
         sinnv = 0.0
         cosmu_nyq = 0.0
         sinmu_nyq = 0.0
         cosnv_nyq = 0.0
         sinnv_nyq = 0.0
         rmnch = 0.0
         zmnsh = 0.0
         currumnch = 0.0
         p_prime = 0.0
!$OMP END WORKSHARE
         IF (lasym) THEN
            rmnsh = 0.0
            zmnch = 0.0
            currumnsh = 0.0
         END IF
      END IF

!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO i = parallel%start(ns - 1), parallel%end(ns - 1)
         si = bmw_parallel_context_i(i, ns)

         rmnch(:,si) = (rmnc(:,si + 1) + rmnc(:,si))/2.0
         zmnsh(:,si) = (zmns(:,si + 1) + zmns(:,si))/2.0

         currumnch(:,si) = (currumnc(:,si + 1) + currumnc(:,si))/2.0

         p_prime(si) = (presf(si + 1) - presf(si))/ds

         IF (lasym) THEN
            rmnsh(:,si) = (rmns(:,si + 1) + rmns(:,si))/2.0
            zmnch(:,si) = (zmnc(:,si + 1) + zmnc(:,si))/2.0

            currumnsh(:,si) = (currumns(:,si + 1) + currumns(:,si))/2.0
         END IF
      END DO
!$OMP END DO

!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO i = parallel%start(num_u), parallel%end(num_u)
         x = (i - 0.5)*du
         cosmu(:,i) = COS(xm*x)
         sinmu(:,i) = SIN(xm*x)
         cosmu_nyq(:,i) = COS(xm_nyq*x)
         sinmu_nyq(:,i) = SIN(xm_nyq*x)
      END DO
!$OMP END DO

!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO i = parallel%start(num_v_use), parallel%end(num_v_use)
         x = (i - 0.5)*primed_grid_construct_jv%dv
         cosv(i) = COS(x)
         sinv(i) = SIN(x)
         cosnv(:,i) = COS(xn*x)
         sinnv(:,i) = SIN(xn*x)
         cosnv_nyq(:,i) = COS(xn_nyq*x)
         sinnv_nyq(:,i) = SIN(xn_nyq*x)
      END DO
!$OMP END DO

      IF (parallel%stride .gt. 1) THEN
!$OMP SINGLE
         CALL parallel%reduce(cosv)
         CALL parallel%reduce(sinv)
         CALL parallel%reduce(cosmu)
         CALL parallel%reduce(sinmu)
         CALL parallel%reduce(cosmu_nyq)
         CALL parallel%reduce(sinmu_nyq)
         CALL parallel%reduce(cosnv)
         CALL parallel%reduce(sinnv)
         CALL parallel%reduce(cosnv_nyq)
         CALL parallel%reduce(sinnv_nyq)
         CALL parallel%reduce(rmnch)
         CALL parallel%reduce(zmnsh)
         CALL parallel%reduce(currumnch)
         CALL parallel%reduce(p_prime)
         IF (lasym) THEN
            CALL parallel%reduce(rmnsh)
            CALL parallel%reduce(zmnch)
            CALL parallel%reduce(currumnsh)
         END IF
!$OMP END SINGLE
      END IF

      ALLOCATE(cosmn(mnmax))
      ALLOCATE(sinmn(mnmax))
      ALLOCATE(cosmn_nyq(mnmax_nyq))
      IF (lasym) THEN
         ALLOCATE(sinmn_nyq(mnmax_nyq))
      END IF

!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO i = parallel%start(ns*num_u*num_v_use),                               &
     &       parallel%end(ns*num_u*num_v_use)
         si = bmw_parallel_context_i(i, ns)
         ui = bmw_parallel_context_j(i, ns, num_u)
         vi = bmw_parallel_context_k(i, ns, num_u)

         cosmn = cosmu(:,ui)*cosnv(:,vi) + sinmu(:,ui)*sinnv(:,vi)
         sinmn = sinmu(:,ui)*cosnv(:,vi) - cosmu(:,ui)*sinnv(:,vi)
         cosmn_nyq = cosmu_nyq(:,ui)*cosnv_nyq(:,vi)                           &
     &             + sinmu_nyq(:,ui)*sinnv_nyq(:,vi)

         r = SUM(rmnch(:,si)*cosmn(:))
         z = SUM(zmnsh(:,si)*sinmn(:))

         ru = -SUM(xm*rmnch(:,si)*sinmn(:))
         rv =  SUM(xn*rmnch(:,si)*sinmn(:))
         zu =  SUM(xm*zmnsh(:,si)*cosmn(:))
         zv = -SUM(xn*zmnsh(:,si)*cosmn(:))

         ju = SUM(currumnch(:,si)*cosmn_nyq(:))

         bu = SUM(bsupumnc(:,si + 1)*cosmn_nyq(:))
         bv = SUM(bsupvmnc(:,si + 1)*cosmn_nyq(:))

         IF (lasym) THEN
            sinmn_nyq = sinmu_nyq(:,ui)*cosnv_nyq(:,vi)                        &
     &                - cosmu_nyq(:,ui)*sinnv_nyq(:,vi)

            r = r + SUM(rmnsh(:,si)*sinmn(:))
            z = z + SUM(zmnch(:,si)*cosmn(:))

            ru = ru + SUM(xm*rmnsh(:,si)*cosmn(:))
            rv = rv - SUM(xn*rmnsh(:,si)*cosmn(:))
            zu = zu - SUM(xm*zmnch(:,si)*sinmn(:))
            zv = zv + SUM(xn*zmnch(:,si)*sinmn(:))

            ju = ju + SUM(currumnsh(:,si)*sinmn_nyq(:))

            bu = bu + SUM(bsupumns(:,si + 1)*sinmn_nyq(:))
            bv = bv + SUM(bsupvmns(:,si + 1)*sinmn_nyq(:))
         END IF

         jv = (ju*bv - p_prime(si))/bu*mu0
         ju = ju*mu0

         jr = ju*ru + jv*rv
         jp = jv*r

         primed_grid_construct_jv%j_z(si,ui,vi) = ju*zu + jv*zv
         primed_grid_construct_jv%j_x(si,ui,vi) = jr*cosv(vi)                  &
     &                                          - jp*sinv(vi)
         primed_grid_construct_jv%j_y(si,ui,vi) = jr*sinv(vi)                  &
     &                                          + jp*cosv(vi)

         primed_grid_construct_jv%x(si,ui,vi) = r*cosv(vi)
         primed_grid_construct_jv%y(si,ui,vi) = r*sinv(vi)
         primed_grid_construct_jv%z(si,ui,vi) = z
      END DO
!$OMP END DO

      DEALLOCATE(cosmn)
      DEALLOCATE(sinmn)
      DEALLOCATE(cosmn_nyq)
      IF (lasym) THEN
         DEALLOCATE(sinmn_nyq)
      END IF
!$OMP END PARALLEL

      DEALLOCATE(cosv)
      DEALLOCATE(sinv)
      DEALLOCATE(cosmu)
      DEALLOCATE(sinmu)
      DEALLOCATE(cosnv)
      DEALLOCATE(sinnv)
      DEALLOCATE(cosmu_nyq)
      DEALLOCATE(sinmu_nyq)
      DEALLOCATE(cosnv_nyq)
      DEALLOCATE(sinnv_nyq)

      DEALLOCATE(rmnch)
      DEALLOCATE(zmnsh)
      DEALLOCATE(currumnch)
      DEALLOCATE(p_prime)
      IF (lasym) THEN
         DEALLOCATE(rmnsh)
         DEALLOCATE(zmnch)
         DEALLOCATE(currumnsh)
      END IF

!  Multi process did not fill out the entire array. Get the missing pieces from
!  the other processes.
      IF (parallel%stride .gt. 1) THEN
         CALL parallel%reduce(primed_grid_construct_jv%x)
         CALL parallel%reduce(primed_grid_construct_jv%y)
         CALL parallel%reduce(primed_grid_construct_jv%z)
         CALL parallel%reduce(primed_grid_construct_jv%j_x)
         CALL parallel%reduce(primed_grid_construct_jv%j_y)
         CALL parallel%reduce(primed_grid_construct_jv%j_z)
      END IF

      CALL profiler_set_stop_time('primed_grid_construct_jv',                  &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref primed_grid_class object.
!>
!>  Allocates memory and initializes a @ref primed_grid_class object with an
!>  This computes the currents and positions on the primed grid. Plasma currents
!>  are obtained from Curl(B) of the siesta solution.
!>
!>  @param[in] num_v       Number of toroidal grid points.
!>  @param[in] siesta_file Name of the siesta restart file.
!>  @param[in] parallel    @ref bmw_parallel_context_class object instance.
!>  @returns A pointer to a constructed @ref primed_grid_class object.
!-------------------------------------------------------------------------------
      FUNCTION primed_grid_construct_siesta(num_v, siesta_file_name,           &
     &                                      parallel)
      USE stel_constants, ONLY: twopi, mu0
      USE read_wout_mod, ONLY: nfp, ns, rmnc, rmns, zmnc, zmns, isigng,        &
     &                         mnmax, xm, xn, lasym
      USE siesta_file

      IMPLICIT NONE

!  Declare Arguments
      CLASS (primed_grid_class), POINTER :: primed_grid_construct_siesta
      INTEGER, INTENT(in)                            :: num_v
      CHARACTER (len=*)                              :: siesta_file_name
      CLASS (bmw_parallel_context_class), INTENT(in) :: parallel

!  local variables
      REAL (rprec)                                   :: start_time
      INTEGER                                        :: i
      REAL (rprec)                                   :: x
      REAL (rprec)                                   :: ds
      INTEGER                                        :: si
      INTEGER                                        :: ui
      INTEGER                                        :: vi
      REAL (rprec)                                   :: r
      REAL (rprec)                                   :: z
      REAL (rprec)                                   :: rs
      REAL (rprec)                                   :: zs
      REAL (rprec)                                   :: ru
      REAL (rprec)                                   :: zu
      REAL (rprec)                                   :: rv
      REAL (rprec)                                   :: zv
      REAL (rprec)                                   :: js
      REAL (rprec)                                   :: ju
      REAL (rprec)                                   :: jv
      REAL (rprec)                                   :: jr
      REAL (rprec)                                   :: jp
      REAL (rprec)                                   :: jz
      INTEGER                                        :: m
      INTEGER                                        :: n
      INTEGER                                        :: ilow
      INTEGER                                        :: ihigh
      REAL (rprec)                                   :: wlow
      REAL (rprec)                                   :: whigh
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: cosv
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: sinv
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: cosmu
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: sinmu
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: cosnv
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: sinnv
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: cosmu_vmec
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: sinmu_vmec
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: cosnv_vmec
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: sinnv_vmec
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: cosmn
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE      :: sinmn
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: cosmn_vmec
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: sinmn_vmec
      REAL (rprec), DIMENSION(:), ALLOCATABLE        :: amnint
      INTEGER                                        :: num_v_use
      CLASS (siesta_file_class), POINTER             :: siesta

!  local parameters
      INTEGER, PARAMETER                             :: num_u = 101
      REAL (rprec), PARAMETER                        :: du = twopi/num_u

!  Start of executable code
      start_time = profiler_get_start_time()

      siesta => siesta_file_class(TRIM(siesta_file_name))

      ALLOCATE(primed_grid_construct_siesta)

      IF (num_v .eq. 1) THEN
         num_v_use = 360
      ELSE
         num_v_use = num_v
      END IF

      ALLOCATE(primed_grid_construct_siesta%x(siesta%nrad,num_u,               &
     &                                        num_v_use))
      ALLOCATE(primed_grid_construct_siesta%y(siesta%nrad,num_u,               &
     &                                        num_v_use))
      ALLOCATE(primed_grid_construct_siesta%z(siesta%nrad,num_u,               &
     &                                        num_v_use))

      ALLOCATE(primed_grid_construct_siesta%j_x(siesta%nrad,num_u,             &
     &                                          num_v_use))
      ALLOCATE(primed_grid_construct_siesta%j_y(siesta%nrad,num_u,             &
     &                                          num_v_use))
      ALLOCATE(primed_grid_construct_siesta%j_z(siesta%nrad,num_u,             &
     &                                          num_v_use))

      ALLOCATE(cosv(num_v_use))
      ALLOCATE(sinv(num_v_use))
      ALLOCATE(cosmu_vmec(mnmax,num_u))
      ALLOCATE(sinmu_vmec(mnmax,num_u))
      ALLOCATE(cosnv_vmec(mnmax,num_v_use))
      ALLOCATE(sinnv_vmec(mnmax,num_v_use))
      ALLOCATE(cosmu(0:siesta%mpol,num_u))
      ALLOCATE(sinmu(0:siesta%mpol,num_u))
      ALLOCATE(cosnv(-siesta%ntor:siesta%ntor,num_v_use))
      ALLOCATE(sinnv(-siesta%ntor:siesta%ntor,num_v_use))

      ds = 1.0/(ns - 1.0)
      primed_grid_construct_siesta%dv = twopi/num_v_use

      primed_grid_construct_siesta%dvol = isigng*ds*du                         &
     &   * primed_grid_construct_siesta%dv/(2.0*twopi)

!$OMP PARALLEL
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(i, x, si, ui, vi, n, m, r, ru, rv, z, zu, zv, rs, zs,           &
!$OMP&         js, ju, jv, jr, jp, jz, ilow, wlow, ihigh, whigh,               &
!$OMP&         cosmn, sinmn, cosmn_vmec, sinmn_vmec, amnint)

!  Multi process will do an all reduce so these arrays need to be initalized.
      IF (parallel%stride .gt. 1) THEN
!$OMP WORKSHARE
         primed_grid_construct_siesta%x = 0.0
         primed_grid_construct_siesta%y = 0.0
         primed_grid_construct_siesta%z = 0.0
         primed_grid_construct_siesta%j_x = 0.0
         primed_grid_construct_siesta%j_y = 0.0
         primed_grid_construct_siesta%j_z = 0.0
         cosmu_vmec = 0.0
         cosnv_vmec = 0.0
         sinmu_vmec = 0.0
         sinnv_vmec = 0.0
         cosmu = 0.0
         cosnv = 0.0
         sinmu = 0.0
         sinnv = 0.0
         cosv  = 0.0
         sinv  = 0.0
!$OMP END WORKSHARE
      END IF

!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO i = parallel%start(num_u), parallel%end(num_u)
         x = (i - 0.5)*du
         cosmu_vmec(:,i) = COS(xm*x)
         sinmu_vmec(:,i) = SIN(xm*x)

         DO m = 0, siesta%mpol
            cosmu(m,i) = COS(m*x)
            sinmu(m,i) = SIN(m*x)
         END DO
      END DO
!$OMP END DO

!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO i = parallel%start(num_v_use), parallel%end(num_v_use)
         x = (i - 0.5)*primed_grid_construct_siesta%dv
         cosv(i) = COS(x)
         sinv(i) = SIN(x)
         cosnv_vmec(:,i) = COS(xn*x)
         sinnv_vmec(:,i) = SIN(xn*x)
         DO n = -siesta%ntor, siesta%ntor
            cosnv(n,i) = COS(n*x)
            sinnv(n,i) = SIN(n*x)
         END DO
      END DO
!$OMP END DO

      IF (parallel%stride .gt. 1) THEN
!$OMP SINGLE
         CALL parallel%reduce(cosv)
         CALL parallel%reduce(sinv)
         CALL parallel%reduce(cosmu_vmec)
         CALL parallel%reduce(sinmu_vmec)
         CALL parallel%reduce(cosnv_vmec)
         CALL parallel%reduce(sinnv_vmec)
         CALL parallel%reduce(cosmu)
         CALL parallel%reduce(sinmu)
         CALL parallel%reduce(cosnv)
         CALL parallel%reduce(sinnv)
!$OMP END SINGLE
      END IF

      ALLOCATE(cosmn(0:siesta%mpol,-siesta%ntor:siesta%ntor))
      ALLOCATE(sinmn(0:siesta%mpol,-siesta%ntor:siesta%ntor))

      ALLOCATE(cosmn_vmec(mnmax))
      ALLOCATE(sinmn_vmec(mnmax))

      ALLOCATE(amnint(mnmax))

!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO i = parallel%start(siesta%nrad*num_u*num_v_use),                      &
     &       parallel%end(siesta%nrad*num_u*num_v_use)
         si = bmw_parallel_context_i(i, siesta%nrad)
         ui = bmw_parallel_context_j(i, siesta%nrad, num_u)
         vi = bmw_parallel_context_k(i, siesta%nrad, num_u)

         x = (si - 1.0)/(siesta%nrad - 1.0)

         CALL primed_grid_siesta_interpol(x, ilow, wlow, ihigh, whigh)

         cosmn_vmec = cosmu_vmec(:,ui)*cosnv_vmec(:,vi)                        &
     &              + sinmu_vmec(:,ui)*sinnv_vmec(:,vi)
         sinmn_vmec = sinmu_vmec(:,ui)*cosnv_vmec(:,vi)                        &
     &              - cosmu_vmec(:,ui)*sinnv_vmec(:,vi)

         DO n = -siesta%ntor, siesta%ntor
            DO m = 0, siesta%mpol
               cosmn(m,n) = cosmu(m, ui)*cosnv(n, vi)                          &
     &                    - sinmu(m, ui)*sinnv(n, vi)
               sinmn(m,n) = sinmu(m, ui)*cosnv(n, vi)                          &
     &                    + cosmu(m, ui)*sinnv(n, vi)
            END DO
         END DO

         amnint = wlow*rmnc(:,ilow) + whigh*rmnc(:,ihigh)
         r = SUM(amnint*cosmn_vmec)
         ru = -SUM(xm*amnint*sinmn_vmec)
         rv =  SUM(xn*amnint*sinmn_vmec)

         amnint = wlow*zmns(:,ilow) + whigh*zmns(:,ihigh)
         z = SUM(amnint*sinmn_vmec)
         zu =  SUM(xm*amnint*cosmn_vmec)
         zv = -SUM(xn*amnint*cosmn_vmec)

         amnint = 2.0*x*(rmnc(:,ihigh) - rmnc(:,ilow))*(ns - 1.0)
         rs = SUM(amnint*cosmn_vmec)

         amnint = 2.0*x*(zmns(:,ihigh) - zmns(:,ilow))*(ns - 1.0)
         zs = SUM(amnint*sinmn_vmec)

         js = SUM(siesta%jksupsmnsf(:,:,si)*sinmn)
         ju = SUM(siesta%jksupumncf(:,:,si)*cosmn)
         jv = SUM(siesta%jksupvmncf(:,:,si)*cosmn)
         IF (BTEST(siesta%flags, siesta_lasym_flag)) THEN
            IF (lasym) THEN
               amnint = wlow*rmns(:,ilow) + whigh*rmns(:,ihigh)
               r = r + SUM(amnint*sinmn_vmec)
               ru = ru + SUM(xm*amnint*cosmn_vmec)
               rv = rv - SUM(xn*amnint*cosmn_vmec)

               amnint = wlow*zmnc(:,ilow) + whigh*zmnc(:,ihigh)
               z = z + SUM(amnint*cosmn_vmec)
               zu = zu - SUM(xm*amnint*sinmn_vmec)
               zv = zv + SUM(xn*amnint*sinmn_vmec)

               amnint = 2.0*x*(rmns(:,ihigh) - rmns(:,ilow))*(ns - 1.0)
               rs = rs + SUM(amnint*sinmn_vmec)

               amnint = 2.0*x*(zmnc(:,ihigh) - zmnc(:,ilow))*(ns - 1.0)
               zs = zs + SUM(amnint*cosmn_vmec)
            END IF

            js = js + SUM(siesta%jksupsmncf(:,:,si)*cosmn)
            ju = ju + SUM(siesta%jksupumnsf(:,:,si)*sinmn)
            jv = jv + SUM(siesta%jksupvmnsf(:,:,si)*sinmn)

         END IF

         js = js/(siesta%b_factor*mu0)
         ju = ju/(siesta%b_factor*mu0)
         jv = jv/(siesta%b_factor*mu0)

         jr = js*rs + ju*ru + jv*rv
         jp = jv*r
         jz = js*zs + ju*zu + jv*zv

         IF (si .eq. 1 .or. si .eq. ns) THEN
            jr = jr/2.0
            jp = jp/2.0
            jz = jz/2.0
         END IF

         primed_grid_construct_siesta%j_x(si,ui,vi) = jr*cosv(vi)              &
     &                                              - jp*sinv(vi)
         primed_grid_construct_siesta%j_y(si,ui,vi) = jr*sinv(vi)              &
     &                                              + jp*cosv(vi)
         primed_grid_construct_siesta%j_z(si,ui,vi) = jz

         primed_grid_construct_siesta%x(si,ui,vi) = r*cosv(vi)
         primed_grid_construct_siesta%y(si,ui,vi) = r*sinv(vi)
         primed_grid_construct_siesta%z(si,ui,vi) = z
      END DO
!$OMP END DO

      DEALLOCATE(cosmn)
      DEALLOCATE(sinmn)

      DEALLOCATE(cosmn_vmec)
      DEALLOCATE(sinmn_vmec)

      DEALLOCATE(amnint)
!$OMP END PARALLEL

      DEALLOCATE(cosv)
      DEALLOCATE(sinv)
      DEALLOCATE(cosmu_vmec)
      DEALLOCATE(sinmu_vmec)
      DEALLOCATE(cosnv_vmec)
      DEALLOCATE(sinnv_vmec)
      DEALLOCATE(cosmu)
      DEALLOCATE(sinmu)
      DEALLOCATE(cosnv)
      DEALLOCATE(sinnv)

!  Multi process did not fill out the entire array. Get the missing pieces from
!  the other processes.
      IF (parallel%stride .gt. 1) THEN
         CALL parallel%reduce(primed_grid_construct_siesta%x)
         CALL parallel%reduce(primed_grid_construct_siesta%y)
         CALL parallel%reduce(primed_grid_construct_siesta%z)
         CALL parallel%reduce(primed_grid_construct_siesta%j_x)
         CALL parallel%reduce(primed_grid_construct_siesta%j_y)
         CALL parallel%reduce(primed_grid_construct_siesta%j_z)
      END IF

      DEALLOCATE(siesta)

      CALL profiler_set_stop_time('primed_grid_construct_siesta',              &
     &                            start_time)

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref primed_grid_class object.
!>
!>  Deallocates memory and uninitializes a @ref primed_grid_class object.
!>
!>  @param[inout] this A @ref primed_grid_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE primed_grid_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (primed_grid_class), INTENT(inout) :: this

!  Start of executable code
      IF (ASSOCIATED(this%x)) THEN
         DEALLOCATE(this%x)
         this%x => null()
      END IF

      IF (ASSOCIATED(this%y)) THEN
         DEALLOCATE(this%y)
         this%y => null()
      END IF

      IF (ASSOCIATED(this%z)) THEN
         DEALLOCATE(this%z)
         this%z => null()
      END IF

      IF (ASSOCIATED(this%j_x)) THEN
         DEALLOCATE(this%j_x)
         this%j_x => null()
      END IF

      IF (ASSOCIATED(this%j_y)) THEN
         DEALLOCATE(this%j_y)
         this%j_y => null()
      END IF

      IF (ASSOCIATED(this%j_z)) THEN
         DEALLOCATE(this%j_z)
         this%j_z => null()
      END IF

      END SUBROUTINE

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
      SUBROUTINE primed_grid_siesta_interpol(s, ilow, wlow,                    &
     &                                       ihigh, whigh)
      USE read_wout_mod, ONLY: ns

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), INTENT(in)  :: s
      INTEGER, INTENT(out)      :: ilow
      REAL (rprec), INTENT(out) :: wlow
      INTEGER, INTENT(out)      :: ihigh
      REAL (rprec), INTENT(out) :: whigh

!  local variables
      REAL (rprec)              :: start_time
      REAL (rprec)              :: wlow_r

!  Start of executable code
      start_time = profiler_get_start_time()

      wlow_r = s*s*(ns - 1.0) + 1.0
      ilow = FLOOR(wlow_r)
      IF (ilow .eq. ns) THEN
         ilow = ns - 1
      END IF
      ihigh = ilow + 1
      wlow = -wlow_r + 1.0 + ilow
      whigh = 1.0 - wlow

      END SUBROUTINE

      END MODULE
