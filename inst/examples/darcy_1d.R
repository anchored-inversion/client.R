require(AnchoredInversionClient)
options(error = utils::recover)

#=====================================================
# Make up synthetic field and data; prep for modeling.

seed <- sample(1000, 1)
#Tough seed cases: 721, 115
set.seed(seed)
cat('set.seed called with parameter seed=', seed)

n_linear_data <- 1
n_forward <- 12

myfield <- AnchoredInversionClient::make_darcy_field()

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

stopifnot(n_forward < mygrid$len / 3)
forward_data_idx <- seq(from = 1, to = mygrid$len, len = n_forward + 2)
forward_data_idx <- round(forward_data_idx[2 : (n_forward + 1)])
    # Uniform sampling
#forward_data_idx <- sample(3 : (mygrid$len - 2), n_forward, replace = FALSE)
    # Random sampling
f_forward <- AnchoredInversionClient::make_darcy_forward_function(mygrid, forward_data_idx)

# forward_data_x <- AnchoredInversionClient::grid_ijk2xyz(mygrid, forward_data_idx)
# forward_data_y <- AnchoredInversionClient::darcy_forward_transform(forward_data, rev = TRUE)
# These two are useful for plotting, otherwise not needed
# in the inverse algorithm.

linear_data <- AnchoredInversionClient::make_darcy_linear_data(mygrid, myfield, n_linear_data)

forward_data <- f_forward(myfield)


#=========================
# Now the actual modeling.

N <- 4 #12
n_start <- 300 #2000
n_finish <- 300 #2000
r <- (n_finish/n_start) ^ (1/(N-1))
n_samples <- round(n_start * r^(0 : (N-1L)))
n_simulations <- 10

source('anchor_common.R')


