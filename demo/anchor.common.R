require(AnchoredInversion)

####### THE REAL DEAL ##########

N <- 4 #12
n.start <- 2000
n.finish <- 2000
r <- (n.finish/n.start) ^ (1/(N-1))
n.samples <- round(n.start * r^(0 : (N-1L)))


cat('\n')
cat('field grid:', paste(mygrid$len, sep = ' x '), '\n')
cat('# of forward data:', length(forward.data), '\n')
cat('sample sizes:', n.samples, '\n')
cat('       total:', sum(n.samples), '\n')
cat('\n')

v <- AnchoredInversion::init_model(
        grid = mygrid,
        data_linear = linear.data,
        field_value_range = range(myfield) +
            c(-1, 1) * runif(2, 2, 10) * diff(range(myfield)),
            #   c(-1, 1) * diff(range(myfield)),
            # A guessed range of the field values.
            # Use a wide range to make the problem more difficult.
            # However, the field is defined on the entire real line.
        data_forward = forward.data
        )
user_id <- v$user_id
task_id <- v$task_id
stamp <- v$stamp

forward.samples <- c()

for (iter in seq_along(n.samples))
{
    cat('\n=== iteration', iter, '===\n')

    forward.sample <- list()
    while (length(forward.sample) < n.samples[iter])
    {
        n_sim <- trunc((n.samples[iter] - length(forward.sample)) * 1.2)
        cat('\nRequesting', n_sim, 'field realizations... ...\n')
        v <- AnchoredInversion::request_fields(
                user_id = user_id,
                task_id = task_id,
                n_sim = n_sim,
                stamp = stamp)

        cat('\nRunning forward model on', length(v$fields), 'field realizations... ...\n')
        forwards <- lapply(v$fields, f.forward)
        forwards_good <- Filter(function(v) !all(is.na(v)), forwards)
        cat('\n   ', length(forwards) - length(forwards_good), 'forward results are invalid\n')

        cat('\nSubmitting', length(forwards), 'forward results (including invalid ones if any)... ...\n')
        stamp <- AnchoredInversion::submit_forwards(
            user_id = user_id,
            task_id = task_id,
            forward_values = forwards,
            stamp = v$stamp)

        forward.sample <- c(forward.sample,
            Filter(function(v) !all(is.na(v)), forwards))
    }

    forward.samples <- c(forward.samples, list(forward.sample))

    cat('\nUpdating approx to posterior... ...\n')
    stamp <- AnchoredInversion::update_model(
        user_id = user_id,
        task_id = task_id)
}


# Summaries

cat('\n')
AnchoredInversion::summarize_model(user_id, task_id)

# TODO: plotting support
# z <- AnchoredInversion::plot_model(user_id, task_id)
# for (x in AnchoredInversionUtils::unpack.lattice.plots(z)) { x11(); print(x)}


### Field simulations ###

cat('\n')
n_sim <= 1000
cat('Requesting', n_sim, 'field simulations...\n')
myfields <- AnchoredInversion::simulate_fields(user_id, task_id, n_sim = n_sim)
cat('Summarizing simulated fields...\n')
myfields_summary <- summary(AnchoredInversionUtils::field.ensemble(myfields, mygrid), field.ref = myfield)
# TODO: plotting support
#z <- plot(myfields_summary)
#for (x in AnchoredInversionUtils::unpack.lattice.plots(z)) { x11(); print(x)}

# TODO: plotting support
# Plot a few simulations.
#x11()
#print(plot(field.ensemble(myfields[1:3], mygrid), field.ref = myfield))


