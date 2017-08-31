cat('\n')
cat('field grid:', paste(mygrid$len, sep = ' x '), '\n')
cat('# of forward data:', length(forward.data), '\n')
cat('sample sizes:', n.samples, '\n')
cat('       total:', sum(n.samples), '\n')
cat('\n')

cat('Opening a session...\n')
cookies <- AnchoredInversionClient::open_session()

cat('Creating a new task...\n')
task_id <- AnchoredInversionClient::create_task(cookies)
cat('    task_id:', task_id, '\n')

cat('Initializing model...\n')
stamp <- AnchoredInversionClient::init_model(
             task_id=task_id,
             mygrid=mygrid,
             field_value_range=field_value_range,
             forward.data=forward.data,
             linear.data=linear.data,
             cookies=cookies)

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
AnchoredInversionClient::print_summary(summ)

cat('\n')
cat('Requesting', n_simulations, 'field simulations...\n')
simulations <- AnchoredInversionClient::simulate_fields(n_simulations, task_id=task_id, cookies=cookies)

