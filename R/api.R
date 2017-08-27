
API_URL <- 'http://localhost:8000/api'


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


#' @export
open_session <- function(user_id='abc', user_email='demo@anchored-inversion.com')
{
    z <- do_post('/open_session', user_id=user_id, user_email=user_email)
    z$cookies
}


#' @export
create_task <- function(cookies)
{
    task_id <- do_post('/create_task', cookies=cookies)
    task_id
}


#' @export
init_model <- function(task_id, mygrid, field_value_range, forward.data, linear.data, cookies)
{
    if (is.null(linear.data)) {
        stamp <- do_post(
                 '/init_model',
                 cookies = cookies,
                 task_id = task_id,
                 grid = json_dumps(mygrid),
                 field_value_range = json_dumps(field_value_range),
                 data_forward = json_dumps(forward.data)
                 )
    } else {
        stamp <- do_post(
                 '/init_model',
                 cookies = cookies,
                 task_id = task_id,
                 grid = json_dumps(mygrid),
                 field_value_range = json_dumps(field_value_range),
                 data_linear = json_dumps(linear.data),
                 data_forward = json_dumps(forward.data)
                 )
    }
    stamp
}


#' @export
update_model <- function(n.samples, task_id, f.forward, cookies, stamp)
{
    n.samp <- 0
    while (n.samp < n.samples)
    {
        n_sim <- trunc((n.samples - n.samp) * 1.2)
        flog.info('Requesting %s field realizations... ...', n_sim)
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
        flog.debug('    fields: %s x %s', length(fields), length(fields[[1]]))

        flog.info('Running forward model on %s field realizations... ...', length(fields))
        forwards <- lapply(fields, f.forward)
        n.good <- sum(sapply(forwards, function(v) if (all(is.na(v))) 0 else 1))
        flog.info('   %s forward results are invalid', length(forwards) - n.good)

        flog.info('Submitting %s forward results (including invalid ones if any)... ...',
                  length(forwards))
        stamp <- do_post('/submit_forwards',
                         cookies = cookies,
                         task_id = task_id,
                         forward_values = json_dumps(do.call(rbind, forwards)),
                         stamp = stamp)
        # JSON converts R matrix to list of lists ([[...], [...],...]),
        # each row being a member list.
        # `NA` in numerical arrays become `null` in JSON.
        # `null` in JSON becomes 'NA' (string) after transported to server.
        n.samp <- n.samp + n.good
    }

    flog.info('Updating approx to posterior... ...')
    stamp <- do_post('/update_model',
                     cookies = cookies,
                     task_id = task_id)
    stamp
}


#' @export
summarize_task <- function(task_id, cookies)
{
    summ <- do_post('/summarize_task', cookies=cookies, task_id=task_id)
    json_loads(summ)
}


#' @export
visualize_task <- function(task_id, cookies)
{
    msg <- do_post('/visualize_task', cookies=cookies, task_id=task_id)
    msg
}


#' @export
showcase_task <- function(task_id, cookies)
{
    msg <- do_post('/showcase_task', cookies=cookies, task_id=task_id)
    msg
}


#' @export
simulate_fields <- function(n, task_id, cookies)
{
    z <- do_post('/simulate', cookies=cookies, task_id=task_id, n=n)
    simulations <- json_loads(z)
    simulations <- split(simulations, row(simulations))
    simulations
}

