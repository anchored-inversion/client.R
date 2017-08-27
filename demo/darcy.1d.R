require(AnchoredInversionClient)
options(error = utils::recover)

darcy.env <- AnchoredInversionClient::darcy.1d(seed = NULL)
attach(darcy.env)

field_value_range = range(myfield) +
    c(-1, 1) * runif(2, 2, 10) * diff(range(myfield))
    #   c(-1, 1) * diff(range(myfield)),
    # A guessed range of the field values.
    # Use a wide range to make the problem more difficult.
# However, the field is defined on the entire real line.

N <- 4 #12
n.start <- 1000 #2000
n.finish <- 1000 #2000
r <- (n.finish/n.start) ^ (1/(N-1))
n.samples <- round(n.start * r^(0 : (N-1L)))

n_simulations <- 1000

source('anchor.common.R')


