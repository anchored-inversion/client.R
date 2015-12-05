
####### THE REAL DEAL ##########

N <- 12
n.start <- 2000
n.finish <- 2000
r <- (n.finish/n.start) ^ (1/(N-1))
n.samples <- round(n.start * r^(0 : (N-1L)))


cat('\n')
cat('field grid:', paste(mygrid$len, sep = ' x '), '\n')
cat('# of forward data:', length(forward.data), '\n')
cat('sample sizes:', n.samples, '\n')
cat('       total:', sum(n.samples), '\n')

v <- AnchoredInversion::init_anchorit(
    grid = mygrid,
    data.linear.grid = linear.data,
    field.value.range = range(myfield) +
        c(-1, 1) * runif(2, 2, 10) * diff(range(myfield)),
#       c(-1, 1) * diff(range(myfield)),
        # A guessed range of the field values.
        # Use a wide range to make the problem more difficult.
        # However, the field is defined on the entire real line.
    forward.data = forward.data
    )
task_id <- v$task_id
stamp <- v$stamp


for (iter in seq_along(n.samples))
{
    cat('\n=== iteration', iter, '===\n')

    cat('\nRequesting field realizations... ...\n')
    v <- AnchoredInversion::request_anchorit_fields(
        task_id = task_id,
        n_sample = n.samples[iter],
        stamp = stamp,
        verbose = 2L)
    fields <- v$fields
    stamp <- v$stamp

    cat('\nRunning forward model on field realizations... ...\n')
    forwards <- lapply(fields, f.forward)

    cat('\nSubmitting forward results... ...\n')
    stamp <- AnchoredInversion::submit_anchorit_forwards(
        task_id = task_id,
        forward_values = forwards,
        stamp = stamp)

    cat('\nUpdating approx to posterior... ...\n')
    stamp <- AnchoredInversion::update_anchorit(
        task_id = task_id,
        verbose = 2L)
}


# Summaries

cat('\n')
AnchoredInversion::summarize_anchorit(task_id)

z <- AnchoredInversion::plot_anchorit(task_id, field_ref = myfield,
    forward_data_groups = forward.data.groups)
for (x in AnchoredInversionUtils::unpack.lattice.plots(z)) { x11(); print(x)}


### Field simulations ###

cat('\n')
cat('simulating fields...\n')
myfields <- AnchoredInversion::simulate_anchorit(task_id, n = 1000)
z <- plot(summary(myfields, field.ref = myfield))
for (x in AnchoredInversionUtils::unpack.lattice.plots(z)) { x11(); print(x)}

# Plot a few simulations.
x11()
print(plot(AnchoredInversion::simulate_anchorit(task_id, n = 3), field.ref = myfield))

