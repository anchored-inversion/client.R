cat('\n')
cat('field grid:', paste(mygrid$len, sep = ' x '), '\n')
cat('# of forward data:', length(forward_data), '\n')
cat('sample sizes:', n_samples, '\n')
cat('       total:', sum(n_samples), '\n')
cat('\n')

# 'domain' and 'port' are for developers of this package doing local tests.
# End user should omit these parameters or set them to NULL, hence use internal default.
domain <- NULL # 'http://localhost'
port <- NULL # 8000
sess <- AnchoredInversionClient::Session$new(domain = domain, port = port)


cat('Logging in...\n')
# You need to use your own account ID and email here.
user_id <- Sys.getenv('DEMO_USER_ID')
user_email <- Sys.getenv('DEMO_USER_EMAIL')
sess$login(user_id = user_id, user_email = user_email)

cat('Setting project to work on\n')
# You need to use your own project ID here.
project_id <- Sys.getenv('DEMO_PROJECT_ID')
sess$set_project(project_id)

cat('Clearing existing models in the project...\n')
sess$clear_models()

cat('Initializing model...\n')
sess$init_models(
     grid=mygrid,
     field_value_range=field_value_range,
     data_forward=forward_data,
     data_linear=linear_data
     )

for (iter in seq_along(n_samples))
{
    cat('\n=== iteration', iter, '===\n')
    sess$update_models(n_samples=n_samples[iter], f_forward=f_forward)
}

summ <- sess$summarize_project()
AnchoredInversionClient::print_summary(summ)

cat('\n')
cat('Requesting', n_simulations, 'field simulations...\n')
simulations <- sess$simulate_fields(n_simulations)
cat('  simulated', length(simulations), 'fields\n')

cat('Logging out...')
sess$logout()


