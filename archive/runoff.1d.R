require(AnchoredInversionClient)
require(lattice)
options(error = utils::recover)

runoff.env <- AnchoredInversionClient::runoff.1d(seed = NULL)
attach(runoff.env)

# #--- plot the forward process and measurements ---
#
# tmp.forward.config <- forward.config
# for (i.config in seq_along(tmp.forward.config))
#     tmp.forward.config[[i.config]][[1]]$time.idx <-
#         1 : runoff.config[[i.config]]$nt
#
# z <- f.runoff(
#     f.field.transform(myfield, rev = TRUE),
#     grid = mygrid,
#     runoff.config = runoff.config,
#     return.config = tmp.forward.config)
# for (i.runoff in seq_along(z))
# {
#     for (i.forward in seq_along(tmp.forward.config[[i.runoff]]))
#     {
#         forward.type <- names(tmp.forward.config[[i.runoff]])[i.forward]
#
#         y <- z[[i.runoff]][[forward.type]]
#             # if 'flowdepth, 'm'; # if 'discharge', m^2 / s
#         x <- seq_along(y) * runoff.config[[i.runoff]]$dt / 60
#             # minutes
#         if (forward.type == 'discharge')
#         {
#             y <- y * 60  # m^2 / min
#             ylab <- expression(paste('Discharge (', m^2, '/min)', sep = ''))
#         } else
#         {
#             ylab <- 'Flow depth (m)'
#         }
#
#         x11(width = 4, height = 4)
#         print(xyplot(y ~ x,
#             type = 'l',
#             xlab = 'Time (min)',
#             ylab = ylab,
#             main = forward.type
#             ))
#     }
# }


field_value_range = range(myfield) +
    c(-1, 1) * runif(2, 2, 10) * diff(range(myfield))
    #   c(-1, 1) * diff(range(myfield)),
    # A guessed range of the field values.
    # Use a wide range to make the problem more difficult.
# However, the field is defined on the entire real line.

source('anchor.common.R')
