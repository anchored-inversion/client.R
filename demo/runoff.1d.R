require(AnchoredInversionClient)
require(lattice)
options(error = utils::recover)

runoff.env <- AnchoredInversionClient::runoff.1d(seed = NULL)
attach(runoff.env)

#--- plot the actual field ---

# trellis.device(color = FALSE, width = 5, height = 4)
# print(xyplot(
#     f.field.transform(myfield, rev = TRUE) ~ AnchoredInversionClient::grid.fullxyz(mygrid),
#     type = 'l',
#     xlab = 'X',
#     ylab = 'Roughness')
#     )
# 
# trellis.device(color = FALSE, width = 5, height = 4)
# print(xyplot(
#     myfield ~ anchoredInversionClient::grid.fullxyz(mygrid),
#     type = 'l',
#     xlab = 'X',
#     ylab = 'Log roughness')
#     )
# 
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


#------------------------------------
# Model inference, summary, plotting.
#------------------------------------

# save.file.prefix <- 'runoff'

source('anchor.common.R')

# forward.data.x <- do.call(c,
#     lapply(forward.config,
#         lapply,
#         function(x) x$time.idx / 60)
#             # Time in 'minutes'
#     )
# 
# z <- plot.anchorit.forward.sample(ai.out,
#     data.forward.groups = forward.data.groups,
#     x.at = forward.data.x,
#     xlab = 'Time (minutes)')
# for (zz in unpack.lattice.plots(z)) {x11(); print(zz)}

# rm(z, zz)
