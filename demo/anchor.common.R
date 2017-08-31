cat('\n')
cat('field grid:', paste(mygrid$len, sep = ' x '), '\n')
cat('# of forward data:', length(forward.data), '\n')
cat('sample sizes:', n.samples, '\n')
cat('       total:', sum(n.samples), '\n')
cat('\n')

cat('Opening a session...\n')
cookies <- AnchoredInversionClient::open_session()

cat('Creating a new project...\n')
project_id <- AnchoredInversionClient::create_project(cookies)
cat('    project_id:', project_id, '\n')

cat('Initializing model...\n')
stamp <- AnchoredInversionClient::init_model(
             project_id=project_id,
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
                   project_id=project_id,
                   f.forward=f.forward,
                   cookies=cookies,
                   stamp=stamp)
}

summ <- AnchoredInversionClient::summarize_project(project_id=project_id, cookies=cookies)
AnchoredInversionClient::print_summary(summ)

cat('\n')
cat('Requesting', n_simulations, 'field simulations...\n')
simulations <- AnchoredInversionClient::simulate_fields(n_simulations, project_id=project_id, cookies=cookies)

