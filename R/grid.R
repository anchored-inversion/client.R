# The model grid is described by a list with three members:
#   list(from = , by = , len = )
# where
#   'from' is coord of the _center_ of the first cell.
#   'by' is step size between cells (ie, cell centers).
#   'len' is number of cells in each dimension.
#
# The three elements are vectors of the same length;
# their length is the spatial dimensionality of the model.
#
# All dimensions should have the same length unit,
# so that 'x = a' and 'y = a' means the same physical length
# in the two directions.
#
# Usually the first dimension is 'X', west-east,
# the second is 'Y', south-north,
# the third is 'Z', bottom-top.



#' Get a coord matrix of all cell centers in the grid,
#' listing points with 'X' increasing the fastest, then 'Y',
#' then 'Z'.
#'
#' @export
grid.fullxyz <- function(grid)
{
    if (is.list(grid))
    {
        if (length(grid$from) == 1)
            matrix(
                seq(from = grid$from, by = grid$by, length = grid$len),
                ncol = 1)
        else
            as.matrix(
                expand.grid(lapply(
                    seq_along(grid$from),
                    function(i) do.call(seq, lapply(grid, '[', i))
                    )
                )
            )
    } else
        grid
}




## Find the (center) locations of grid cells.

#' @export
grid.ijk2xyz <- function(
    grid,
    ijk
        # A matrix of array subscripts if 'grid' is a list;
        # in this case, 'ijk' may be vector if 1-D or is a single point.
        # A vector of point indices if 'grid' is a matrix.
        # Fractions allowed, but may not be very useful.
    )
{
    if (is.list(grid))
    {
        ndim <- length(grid$from)

        if (is.matrix(ijk))
        {
            n <- nrow(ijk)
            stopifnot(
                ncol(ijk) == ndim,
                ijk > 0.5,
                    # 0.5 indicates the boundary, i.e.
                    # half cell before the cell with index 1.
                ijk < rep(grid$len + 0.5, each = n)
                )
            xyz <- rep(grid$from, each = n) +
                rep(grid$by, each = n) * (ijk - 1)
                # 'xyz' has the same 'dim' as 'ijk'.
        } else
        {
            stopifnot(
                ijk > 0.5,
                ijk < grid$len + 0.5,
                ndim == 1 || length(ijk) == ndim)
            xyz <- matrix(grid$from + grid$by * (ijk - 1), ncol = ndim)
        }
    } else
    {
        if (is.matrix(ijk))
        {
            stopifnot(ncol(ijk) == 1L)
            dim(ijk) <- NULL
        }
        xyz <- grid.ind2xyz(grid, ijk)
    }

    xyz
}




## Find the (center) locations of grid cells
## given cell indices.

#' @export
grid.ind2xyz <- function(
    grid,
    ind   # Vector of cell indices (column major).
    )
{
    if (is.list(grid))
        grid.ijk2xyz(grid, ind2sub(ind, grid$len))
    else
        grid[ind, , drop = FALSE]
}




## Find the grid indices of locations represented by
## (x,y,z) coordinates.
## Return an index vector (if 1D and 'xyz' is vector)
## or an index matrix (if > 1D, or 1D and 'xyz' is a 1-col matrix).

#' @export
grid.xyz2ijk <- function(
    grid,
        # as returned from 'flexgrid'.
    xyz,
        # Coord matrix. If space is 1D, can be a vector as well.
    rounded = TRUE
    )
{
    stopifnot(is.list(grid))
    ndim <- length(grid$from)

    if (is.matrix(xyz))
    {
        n <- nrow(xyz)
        stopifnot(
            ncol(xyz) == ndim,
            xyz > rep(grid$from - grid$by/2, each = n),
            xyz < rep(grid$from + grid$by * (grid$len - 0.5),
                each = n)
            )
        ijk <- (xyz - rep(grid$from, each = n)) /
            rep(grid$by, each = n) + 1
    } else
    {
        if (ndim == 1)
        {
            stopifnot(
                xyz > grid$from - grid$by/2,
                xyz < grid$from +
                    grid$by * (grid$len - 0.5)
                )
            ijk <- (xyz - grid$from) / grid$by + 1
        } else
        {
            stopifnot(
                length(xyz) == ndim,
                xyz > grid$from - grid$by/2,
                xyz < grid$from +
                    grid$by * (grid$len - 0.5)
                )
            ijk <- (xyz - grid$from) / grid$by + 1
        }
    }

    if (rounded)
        as.integer(round(ijk))
    else
        ijk
}

