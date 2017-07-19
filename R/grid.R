## The model grid is described by a list with three members:
##   list(from = , by = , len = )
## where
##   'from' is coord of the _center_ of the first cell.
##   'by' is step size between cells (ie, cell centers).
##   'len' is number of cells in each dimension.
##
## The three elements are vectors of the same length;
## their length is the spatial dimensionality of the model.
##
## All dimensions should have the same length unit,
## so that 'x = a' and 'y = a' means the same physical length
## in the two directions.
##
## Usually the first dimension is 'X', west-east,
## the second is 'Y', south-north,
## the third is 'Z', bottom-top.
##
## I chose these names for the elements such that, e.g.
## the sequence of the grid centers in the first dimension is
##   do.call(seq, lapply(grid, '[', 1))
##
## However, to accommodate some special needs,
## I also allow 'grid' to be a coord matrix for irregularly spaced
## points.
##
## This list or matrix object is the 'grid' argument
## to the following functions.


# #' @export
# is.grid <- function(x)
# {
# # A valid grid definition is either a triple (from, by, len) spec
# # or a coord matrix.
# 
#     (is.list(x) && length(x) == 3L &&
#         identical(sort(names(x)), c('by', 'from', 'len')) &&
#         length(unique(sapply(x, length))) == 1L
#         ) ||
#     (is.numeric(x) && is.matrix(x))
# }



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

# #' @export
# grid.ind2xyz <- function(
#     grid,
#     ind   # Vector of cell indices (column major).
#     )
# {
#     if (is.list(grid))
#         grid.ijk2xyz(grid, ind2sub(ind, grid$len))
#     else
#         grid[ind, , drop = FALSE]
# }




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




## Find the grid indices of locations represented by
## (x,y,z) coordinates.
## Return an index vector.

# #' @export
# grid.xyz2ind <- function(
#     grid,
#     xyz # Coord matrix. If space is 1D, can be a vector.
#     )
# {
#     if (is.list(grid))
#     {
#         sub2ind(grid.xyz2ijk(grid, xyz), grid$len)
#     } else
#     {
#         if (!is.matrix(xyz))
#             dim(xyz) <- c(1, length(xyz))
#         if (!is.matrix(grid))
#             dim(grid) <- c(1, length(grid))
# 
#         dist <- distmat(xyz, grid)
#         apply(dist, 1, which.min)
#     }
# }
# 



# #' Inverse of \code{sub2ind}.
# #'
# #' @export
# ind2sub <- function(ind, dims)
# {
#     ndim <- length(dims)
#     if (ndim == 1)
#         return(as.integer(c(ind)))
# 
#     z <- matrix(0, nrow = length(ind), ncol = ndim)
#     z[, 1] <- ind
#     n <- prod(dims)
#     for (i in ndim : 2)
#     {
#         n <- n / dims[i]
#         z[, i] <- (z[, 1] - 1) %/% n + 1
#         z[, 1] <- z[, 1] - (z[, i] - 1) * n
#     }
# 
#     if (nrow(z) == 1) as.integer(c(z)) else array(as.integer(z), dim=dim(z))
# }




# #' Convert array subscripts to single indices.
# #'
# #' Value: a vector.
# #'
# #' @export
# sub2ind <- function(
#     sub,
#         # Subscripts of elements.
#         # Three possibilities:
#         # (1) length(dims) > 1,
#         #     'sub' is matrix with length(dims) columns.
#         #     Each row identifies one point.
#         # (2) length(dims) > 1,
#         #     'sub' is vector of length length(dims).
#         #     Identifies one point.
#         # (3) length(dims) == 1,
#         #     'sub' is a vector or a column matrix.
#     dims
#         # 'dim' of the array.
#     )
# {
#     ndim <- length(dims)
#     if (ndim == 1)
#         return(c(sub))
# 
#     if (is.matrix(sub))
#     {
#         stopifnot(ncol(sub) == ndim)
#         z <- sub[, 1]
#         n <- 1
#         for (i in 2 : ndim)
#         {
#             n <- n * dims[i-1]
#             z <- z + (sub[,i] - 1) * n
#         }
#     } else
#     {
#         stopifnot(length(sub) == ndim)
#         z <- sub[1]
#         n <- 1
#         for (i in 2 : ndim)
#         {
#             n <- n * dims[i - 1]
#             z <- z + (sub[i] - 1) * n
#         }
#     }
# 
#     as.integer(z)
# }
# 

