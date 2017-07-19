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

cat('Opening session...\n')
cookies <- AnchoredInversionClient::open_session()

cat('Creating a new task...\n')
task_id <- AnchoredInversionClient::create_task(cookies)
#cat('    task_id:', task_id, '\n')

field_value_range = range(myfield) +
    c(-1, 1) * runif(2, 2, 10) * diff(range(myfield))
    #   c(-1, 1) * diff(range(myfield)),
    # A guessed range of the field values.
    # Use a wide range to make the problem more difficult.
# However, the field is defined on the entire real line.

cat('Initializing model...\n')
stamp <- AnchoredInversionClient::init_model(
                                             task_id=task_id,
                                             mygrid=mygrid,
                                             field_value_range=field_value_range,
                                             forward.data=forward.data,
                                             linear.data=linear.data,
                                             cookies=cookies)
#cat('    stamp:', stamp, '\n')

forward.samples <- c()

for (iter in seq_along(n.samples))
{
    cat('\n=== iteration', iter, '===\n')
    stamp <- AnchoredInversionClient::update_model(
                                                   n.samples=n.samples[iter],
                                                   task_id=task_id,
                                                   f.forward=f.forward,
                                                   cookies=cookies,
                                                   stamp=stamp)
}


summ <- AnchoredInversionClient::summarize_task(task_id=task_id, cookies=cookies)
AnchoredInversionClient::print.summary(summ)


cat('\n')
n_sim <- 1000
cat('Requesting', n_sim, 'field simulations...\n')
simulations <- AnchoredInversionClient::simulate_fields(n_sim, task_id=task_id, cookies=cookies)

#cat('Summarizing simulated fields...\n')
#myfields_summary <- summary(AnchoredInversionUtils::field.ensemble(myfields, mygrid), field.ref = myfield)

