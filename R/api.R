API_URL <- 'http://localhost:8000'

make_url <- function(url) {
    paste(API_URL, url, sep='')
}


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


global_cookies <- new.env(parent = emptyenv())


get_cookies <- function() {
    z <- mget('cookies', global_cookies, inherits = FALSE, ifnotfound = list(NULL))
    z[[1]]
}


set_cookies <- function(value) {
    assign('cookies', value, pos = global_cookies, inherits = FALSE)
}


http_call <- function(method, url, ...) {
    cookies <- get_cookies()
    if (is.null(cookies)) {
        z = method(url, ...)
    } else {
        z = method(url, set_cookie(cookies), ...)
    }
    set_cookies(httr::cookies(z))
    zz <- httr::content(z)
    if (is.null(zz)) invisible() else json_loads(zz)
}


#' @export
http_get <- function(url, ...) {
    http_call(httr::GET, url, query = list(...))
}


#' @export
http_post <- function(url, ...) {
    http_call(httr::POST, url, body = list(...))
}


#' @export
login_demo <- function()
{
    http_post(make_url('/login_demo'))
}


#' @export
logout <- function()
{
    http_post(make_url('/logout'))
}


#' @export
get_project_ids <- function()
{
    http_get(make_url('/user/projects'))
}


#' @export
set_project <- function(project_id)
{
    http_post(make_url('/user/set_project'), project_id = project_id)
}


#' @export
clear_models <- function()
{
    http_post(make_url('/user/project/models/clear'))
}


#' @export
init_model <- function(mygrid, field_value_range, forward.data, linear.data)
{
    if (is.null(linear.data)) {
        http_post(
                 make_url('/user/project/models/init'),
                 grid = json_dumps(mygrid),
                 field_value_range = json_dumps(field_value_range),
                 data_forward = json_dumps(forward.data)
                 )
    } else {
        http_post(
                 make_url('/user/project/models/init'),
                 grid = json_dumps(mygrid),
                 field_value_range = json_dumps(field_value_range),
                 data_linear = json_dumps(linear.data),
                 data_forward = json_dumps(forward.data)
                 )
    }
}


#' @export
update_model <- function(n.samples, f.forward)
{
    n.samp <- 0
    while (n.samp < n.samples)
    {
        n_sim <- trunc((n.samples - n.samp) * 1.2)
        flog.info('Requesting %s field realizations... ...', n_sim)
        fields <- http_post(make_url('/user/project/models/request_fields'), n = n_sim)
        fields <- split(fields, row(fields)) # from matrix to list
        # cat('    stamp:', stamp, '\n')
        flog.debug('    fields: %s x %s', length(fields), length(fields[[1]]))

        flog.info('Running forward model on %s field realizations... ...', length(fields))
        forwards <- lapply(fields, f.forward)
        n.good <- sum(sapply(forwards, function(v) if (all(is.na(v))) 0 else 1))
        flog.info('   %s forward results are invalid', length(forwards) - n.good)

        flog.info('Submitting %s forward results (including invalid ones if any)... ...',
                  length(forwards))
        http_post(make_url('/user/project/models/submit_forwards'),
                  forward_values = json_dumps(do.call(rbind, forwards)))
        # JSON converts R matrix to list of lists ([[...], [...],...]),
        # each row being a member list.
        # `NA` in numerical arrays become `null` in JSON.
        # `null` in JSON becomes 'NA' (string) after transported to server.
        n.samp <- n.samp + n.good
    }

    flog.info('Updating approx to posterior... ...')
    http_post(make_url('/user/project/models/update'))
}


#' @export
summarize_project <- function()
{
    http_get(make_url('/user/project/summary'))
}


#' @export
simulate_fields <- function(n)
{
    simulations <- http_get(make_url('/user/project/request_fields'), n=n)
    simulations <- split(simulations, row(simulations))
    simulations
}

