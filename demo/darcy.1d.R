require(AnchoredInversionClient)
options(error = utils::recover)

#=====================================================
# Make up synthetic field and data; prep for modeling.

seed <- sample(1000, 1)
#Tough seed cases: 721, 115
set.seed(seed)
cat('set.seed called with parameter seed=', seed)

n.linear.data <- 1
n.forward <- 12

myfield <- AnchoredInversionClient::make.darcy.field()

field_value_range = range(myfield) +
    c(-1, 1) * runif(2, 2, 10) * diff(range(myfield))
    #   c(-1, 1) * diff(range(myfield)),
    # A guessed range of the field values.
    # Use a wide range to make the problem more difficult.
# However, the field is defined on the entire real line.

mygrid <- list(
    from = .5,
    by = 1,
    len = length(myfield)
    )

stopifnot(n.forward < mygrid$len / 3)
forward.data.idx <- seq(from = 1, to = mygrid$len, len = n.forward + 2)
forward.data.idx <- round(forward.data.idx[2 : (n.forward + 1)])
    # Uniform sampling
#forward.data.idx <- sample(3 : (mygrid$len - 2), n.forward, replace = FALSE)
    # Random sampling
f.forward <- AnchoredInversionClient::make.darcy.forward.function(mygrid, forward.data.idx)

# forward.data.x <- grid.ijk2xyz(mygrid, forward.data.idx)
# forward.data.y <- darcy.forward.transform(forward.data, rev = TRUE)
# These two are useful for plotting, otherwise not needed
# in the inverse algorithm.

linear.data <- AnchoredInversionClient::make.darcy.linear.data(mygrid, myfield, n.linear.data)

forward.data <- f.forward(myfield)


#=========================
# Now the actual modeling.

N <- 6 #12
n.start <- 1000 #2000
n.finish <- 1000 #2000
r <- (n.finish/n.start) ^ (1/(N-1))
n.samples <- round(n.start * r^(0 : (N-1L)))
n_simulations <- 10

source('anchor.common.R')


