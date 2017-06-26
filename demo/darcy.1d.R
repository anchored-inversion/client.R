require(AnchoredInversionExamples)
require(lattice)
require(AnchoredInversionUtils)

darcy.env <- AnchoredInversionExamples::darcy.1d(seed = NULL)
attach(darcy.env)

#--- plot the actual field ---

trellis.device(color = FALSE, width = 5, height = 4)
print(xyplot(
    f.field.transform(myfield, rev = TRUE)
        ~ AnchoredInversionUtils::grid.fullxyz(mygrid),
    type = 'l',
    xlab = 'X',
    ylab = 'Conductivity')
    )

trellis.device(color = FALSE, width = 5, height = 4)
print(xyplot(
    myfield ~ AnchoredInversionUtils::grid.fullxyz(mygrid),
    type = 'l',
    xlab = 'X',
    ylab = 'Log conductivity')
    )

#--- plot the forward process and measurements ---

forward.curve <- list(
    x = AnchoredInversionUtils::grid.fullxyz(mygrid),
    y = f.darcy(f.field.transform(myfield, rev = TRUE))
    )

trellis.device(color = FALSE, width = 4, height = 4)
print(xyplot(
    forward.curve$y ~ forward.curve$x,
    panel = function(...) {
            panel.xyplot(...)
            lpoints(forward.data.x, forward.data.y, pch = 20)
        },
    type = 'l',
    lty = 'dashed',
    xlab = 'X (spatial coord)',
    ylab = 'Z (head)'
    ))


#------------------------------------
# Model inference, summary, plotting.
#------------------------------------

source('anchor.common.R')

trellis.device(color = TRUE, width = 6, height = 7)
z <- AnchoredInversionUtils::plot.forward.reproduction(
        forward.data = forward.data,
        forward.sample = lapply(
            lapply(forward.samples, function(v) do.call(rbind, v)),
            apply, 2, boxplot.stats),
        x.at = forward.data.x,
        xlab = 'X coord')
print(z)


