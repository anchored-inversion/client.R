cat('\n')
cat('field grid:', paste(mygrid$len, sep = ' x '), '\n')
cat('# of forward data:', length(forward.data), '\n')
cat('sample sizes:', n.samples, '\n')
cat('       total:', sum(n.samples), '\n')
cat('\n')

cat('Logging in...\n')
AnchoredInversionClient::login_demo()

cat('Getting a project to work on...\n')
project_ids <- AnchoredInversionClient::get_project_ids()
cat('projects:')
print(project_ids)
project_id <- project_ids[1]
cat('    project_id:', project_id, '\n')

cat('Setting project to work on\n')
AnchoredInversionClient::set_project(project_id)

cat('Clearing existing models...\n')
AnchoredInversionClient::clear_models()

cat('Initializing model...\n')
AnchoredInversionClient::init_model(
     mygrid=mygrid,
     field_value_range=field_value_range,
     forward.data=forward.data,
     linear.data=linear.data
     )

for (iter in seq_along(n.samples))
{
    cat('\n=== iteration', iter, '===\n')
    AnchoredInversionClient::update_model(n.samples=n.samples[iter], f.forward=f.forward)
}

summ <- AnchoredInversionClient::summarize_project()
AnchoredInversionClient::print_summary(summ)

cat('\n')
cat('Requesting', n_simulations, 'field simulations...\n')
simulations <- AnchoredInversionClient::simulate_fields(n_simulations)

cat('Logging out...')
AnchoredInversionClient::logout()


