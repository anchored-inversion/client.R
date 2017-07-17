require(jsonlite)
require(httr)

API_URL <- 'http://localhost:8000/api'
USER_ID <- 'abc'
USER_EMAIL <- 'demo@anchored-inversion.com'


json_dumps <- function(x) {
    jsonlite::toJSON(x)
}


json_loads <- function(x) {
    jsonlite::fromJSON(x)
}


set_cookie <- function(cookies) {
    cookie <- setNames(as.list(cookies$value), cookies$name)
    do.call(httr::set_cookies, cookie)
}


do_post <- function(url, cookies=NULL, ...) {
    url = paste(API_URL, url, sep='')
    if (is.null(cookies)) {
        z = httr::POST(url, body=list(...))
    } else {
        z = httr::POST(url, set_cookie(cookies), body=list(...))
    }
    if (is.null(cookies)) {
        list(value=httr::content(z), cookies=httr::cookies(z))
    } else {
        httr::content(z)
    }
}


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

cat('Opening session...\n')
z <- do_post('/open_session', user_id=USER_ID, user_email=USER_EMAIL)
cookies <- z$cookies

cat('Creating a new task...\n')
task_id <- do_post('/create_task', cookies=cookies)
#cat('    task_id:', task_id, '\n')

cat('Initializing model...\n')
field_value_range = range(myfield) +
    c(-1, 1) * runif(2, 2, 10) * diff(range(myfield))
    #   c(-1, 1) * diff(range(myfield)),
    # A guessed range of the field values.
    # Use a wide range to make the problem more difficult.
    # However, the field is defined on the entire real line.
stamp <- do_post('/init_model',
                 cookies = cookies,
                 task_id = task_id,
                 grid = json_dumps(mygrid),
                 data_linear=json_dumps(linear.data),
                 field_value_range = json_dumps(field_value_range),
                 data_forward = json_dumps(forward.data)
                 )
#cat('    stamp:', stamp, '\n')

forward.samples <- c()

for (iter in seq_along(n.samples))
{
    cat('\n=== iteration', iter, '===\n')

    forward.sample <- list()
    while (length(forward.sample) < n.samples[iter])
    {
        n_sim <- trunc((n.samples[iter] - length(forward.sample)) * 1.2)
        cat('\nRequesting', n_sim, 'field realizations... ...\n')
        z <- do_post('/request_fields',
                     cookies = cookies,
                     task_id = task_id,
                     n = n_sim,
                     stamp = stamp)
        z <- json_loads(z)
        stamp <- z$stamp
        fields <- z$fields
        fields <- split(fields, row(fields)) # from matrix to list
        # cat('    stamp:', stamp, '\n')
        cat('    fields:', length(fields), 'x', length(fields[[1]]), '\n')

        cat('\nRunning forward model on',
            length(fields), 'field realizations... ...\n')
        forwards <- lapply(fields, f.forward)
        forwards_good <- Filter(function(v) !all(is.na(v)), forwards)
        cat('\n   ', length(forwards) - length(forwards_good),
            'forward results are invalid\n')

        cat('\nSubmitting', length(forwards),
            'forward results (including invalid ones if any)... ...\n')
        stamp <- do_post('/submit_forwards',
                         cookies = cookies,
                         task_id = task_id,
                         forward_values = json_dumps(do.call(rbind, forwards)),
                         stamp = stamp)
        # JSON converts R matrix to list of lists ([[...], [...],...]),
        # each row being a member list.
        # `NA` in numerical arrays become `null` in JSON.
        # `null` in JSON becomes 'NA' (string) after transported to server.

        forward.sample <- c(forward.sample, forwards_good)
    }

    forward.samples <- c(forward.samples, list(forward.sample))

    cat('\nUpdating approx to posterior... ...\n')
    stamp <- do_post('/update_model',
                     cookies = cookies,
                     task_id = task_id)
}


# Summaries
# TODO: get and print summary


### Field simulations ###

cat('\n')
n_sim <= 1000
cat('Requesting', n_sim, 'field simulations...\n')
z <- do_post('/simulate', cookies=cookies, task_id=task_id, n=n_sim)
simulations <- json_loads(z)
simulations <- split(simulations, row(simulations))

#cat('Summarizing simulated fields...\n')
#myfields_summary <- summary(AnchoredInversionUtils::field.ensemble(myfields, mygrid), field.ref = myfield)

