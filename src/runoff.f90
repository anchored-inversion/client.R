! 1-D surface rainfall-runoff finite difference.
! Reference:
!   Huang and Lee, "Influences of spatially heterogeneous roughness on
!   flow hydrographs, Advances in Water Resources, 32 (2009): 1580--1587.
!
!   Vieux, "Distributed Hydrologic Modeling Using GIS", 2nd ed, 2004,
!   Kluwer. Chapter 6.
!
!   Engman, "Roughness coefficients for routing surface runoff",
!   J. Irrig. and Drain. Engrg, 112(1): 39--53, 1986.
!
! Uses the MacCormack method.
! Reference:
!   MacCormack, R. W., The effect of viscosity in hypervelocity impact
!     cratering, AIAA Paper, 69-354, 1969.
!   Anderson, J. D., Jr., Computational Fluid Dynamics: The Basics with
!     Applications, McGraw Hill, 1994.
!
! See function 'runoff.1d' in '../R'.
subroutine f_runoff_1d(dx, nx, dt, nt, h0, q0, qt_ub, rain, rough, s0, &
    discharge_space_idx, discharge_space_n, &
    discharge_time_idx, discharge_time_n, &
    flowdepth_space_idx, flowdepth_space_n, &
    flowdepth_time_idx, flowdepth_time_n, &
    discharge, flowdepth, flag)
double precision, intent(in) :: dx
    ! Size (m) of each cell, e.g. 0.5.
integer, intent(in) :: nx
    ! Number of cells.
    ! There are 'nx' cells and 'nx + 1' cell boundaries called
    ! 'sampling points'. The first sampling point is the upstream
    ! boundary.
double precision, intent(in) :: dt
    ! Time step (seconds), e.g. 0.5.
integer, intent(in) :: nt
    ! Number of time steps to calculate.
double precision, intent(inout) :: h0(nx + 1)
    ! Flow depth (m) initial condition at sampling points.
double precision, intent(inout) :: q0(nx + 1)
    ! Initial flow discharge per unit width (m^2 / sec)
    ! at sampling points.
double precision, intent(in) :: qt_ub(nt)
    ! Discharge input at upstream boundary for every subsequent
    ! time step.
double precision, intent(in) :: rain(nx * nt)
    ! Rainfall (m/sec) input in each cell during each time step.
    ! Example value: 10 mm/hr = 10 / 1000 / 3600 m/sec.
    ! The condition for cells 1 : nx in time step 1 is stored first,
    ! then time step 2, and so on.
double precision, intent(in) :: rough(nx)
    ! Average Manning's roughness coefficient (T / L^(1/3)) in each cell.
    ! A reasonable value is 0.02.
    ! See Vieux (2004).
double precision, intent(in) :: s0(nx)
    ! Average bed slope (L/L) in each cell.
    ! A reasonable value is 0.01.
integer, intent(in) :: discharge_space_n, discharge_time_n
integer, intent(in) :: discharge_space_idx(discharge_space_n)
integer, intent(in) :: discharge_time_idx(discharge_time_n)
    ! Space index and time index to be sampled and returned.
    ! Space index runs from 0 to nx.
integer, intent(in) :: flowdepth_space_n, flowdepth_time_n
integer, intent(in) :: flowdepth_space_idx(flowdepth_space_n)
integer, intent(in) :: flowdepth_time_idx(flowdepth_time_n)
double precision, intent(out) :: discharge(&
    discharge_space_n * discharge_time_n)
double precision, intent(out) :: flowdepth(&
    flowdepth_space_n * flowdepth_time_n)
! Return discharge (m^2 / sec) and flow depth (m) at requested
! sampling points and times.
! Flow velocity can be obtained from this as 'discharge / flowdepth'.
integer, intent(out) :: flag

    ! Index '1' is upstream; 'nx' is downstream.

    double precision :: rough_db(nx)
        ! Roughness at downstream boundaries, approximately.
    double precision :: dtx, g_sqr
    double precision :: h1(nx + 1), h2(nx + 1), q1(nx + 1)
    double precision :: sf(nx)
    integer :: it, ix, i, j
    logical :: ret_discharge, ret_flowdepth
    integer :: ret_discharge_n, ret_flowdepth_n
    integer :: discharge_time_now, flowdepth_time_now
        ! Point to the next element in 'discharge_time_idx'
        ! and 'flowdepth_time_idx'.

    ret_discharge = (discharge_space_n .gt. 0 &
        .and.  discharge_time_n .gt. 0)
    ret_flowdepth = (flowdepth_space_n .gt. 0 &
        .and.  flowdepth_time_n .gt. 0)
    discharge_time_now = 1
    flowdepth_time_now = 1

    rough_db = (/ (rough(1 : (nx-1)) + rough(2 : nx)) * 0.5, &
        rough(nx) /)
    dtx = dt / dx
    g_sqr = 3.13

    flag = 0

    do it = 1, nt
        h1(2 : (nx+1)) = rain(((it - 1) * nx + 1) : (it * nx)) * dt
        h2(1 : nx) = h1(2 : (nx+1))

        ! Prediction step.
        h1(2 : (nx+1)) = h1(2 : (nx+1)) &
            + h0(2 : (nx+1)) &
            - dtx * (q0(2 : (nx+1)) - q0(1 : nx))
            ! At downstream boundaries.
        h1(1) = h0(1)
        sf = s0 - (h1(2 : (nx+1)) - h1(1 : nx)) / dx
            ! Friction slope (L/L).
            ! At cell midpoints.
        sf = (/ (sf(1 : (nx-1)) + sf(2 : nx)) * .5, sf(nx) /)
            ! Friction slopt at downstream boundaries.

        q1(2 : (nx+1)) = h1(2 : (nx+1)) ** 1.66667 * sqrt(sf) / rough_db
        q1(1) = qt_ub(it)

        ! Correction step.
        h2(1 : nx) = h2(1 : nx) &
            + h0(1 : nx) &
            - dtx * (q1(2 : (nx+1)) - q1(1 : nx))
        h2(nx+1) = h0(nx+1)

        ! Summary step.
        h0 = (h1 + h2) / 2
        where (h0 < 0) h0 = 0
        sf = s0 - (h0(2 : (nx+1)) - h0(1 : nx)) / dx
        sf = (/ (sf(1 : (nx-1)) + sf(2 : nx)) * .5, sf(nx) /)
        q0(2 : (nx+1)) = h0(2 : (nx+1)) ** 1.66667 * sqrt(sf) / rough_db
        q0(1) = qt_ub(it)

        ! Courant-Friedrichs-Lewy (CFL) criterion of stability.
        if (any(dtx * (q0 + h0 * sqrt(h0) * g_sqr) > h0)) then
            ! 'numerical stability criterion (CFL) violdated'
            flag = 1
            exit
        end if

        if (ret_discharge .and. &
                (it .eq. discharge_time_idx(discharge_time_now))) then
            discharge(&
                ((discharge_time_now-1) * discharge_space_n + 1) :&
                (discharge_time_now * discharge_space_n) ) &
                = q0(discharge_space_idx + 1)
            discharge_time_now = discharge_time_now + 1
        endif

        if (ret_flowdepth .and. &
                (it .eq. flowdepth_time_idx(flowdepth_time_now))) then
            flowdepth(&
                ((flowdepth_time_now-1) * flowdepth_space_n + 1) :&
                (flowdepth_time_now * flowdepth_space_n) ) &
                = h0(flowdepth_space_idx + 1)
            flowdepth_time_now = flowdepth_time_now + 1
        endif
    end do

end subroutine f_runoff_1d

