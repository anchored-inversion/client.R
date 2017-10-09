#' @importFrom futile.logger flog.debug flog.info flog.warn flog.error
NULL

#' HTTP client API
#'
#' @docType topic
#' @name HTTP API
#' @seealso \code{\link{Session}}
NULL


json_dumps <- function(x) {
    jsonlite::toJSON(x)
}


json_loads <- function(x) {
    jsonlite::fromJSON(x)
}


http_call <- function(method, url, cookies, ...) {
    if (is.null(cookies)) {
        z <- method(url, ...)
    } else {
        cookies <- do.call(httr::set_cookies, setNames(as.list(cookies$value), cookies$name))
        z <- method(url, cookies, ...)
    }
    if (httr::status_code(z) != 200) {
        cat('Something went wrong!\n')
        cat('status_code is', httr::status_code(z), '\n')
        cat(httr::content(z), '\n')
        stop('halting the program... ...')
            # Of course we'd like to do something more elegant. Maybe later.
    }
    cookies <- httr::cookies(z)
    value <- httr::content(z)
    if (!is.null(value)) value <- json_loads(value)
    list(cookies = cookies, value = value)
}


# Documentation for this class is hand written in file 'man/Session.Rd'
# because roxygen2 has difficulties in documenting R6 classes.
# Do not add roxygen2 comments here except for the '@export' line.

#' @export
Session <- R6::R6Class("Session",
    public = list(
        initialize = function(base_url = NULL) {
            if (!is.null(base_url)) {
                while (substr(base_url, nchar(base_url), nchar(base_url)) == '/') {
                    base_url <- substr(base_url, 1, nchar(base_url) - 1)
                }
                private$base_url <- base_url
            }
        },

        login_demo = function() {
            if (!is.null(private$cookies)) {
                stop('You are already logged in. Please log out first before attempting another log-in.')
            }
            private$do_post('/login_demo')
        },

        logout = function() {
            if (!is.null(private$cookies)) {
                private$do_post('/logout')
                # TODO: the '/logout' API should clear the cookie.
                private$cookies <- NULL
            }
            private$project_id_ <- NULL
        },

        set_project = function(project_id) {
            private$ensure_logged_in()
            private$do_post('/user/set_project', project_id = project_id)
            private$project_id_ <- project_id
        },

        clear_models = function() {
            private$ensure_in_project()
            private$do_post('/user/project/models/clear')
        },

        init_models = function(grid, field_value_range, data_forward, data_linear) {
            private$ensure_in_project()
            args <- list(
                url = '/user/project/models/init',
                grid = json_dumps(grid),
                field_value_range = json_dumps(field_value_range),
                data_forward = json_dumps(data_forward))
            if (!is.null(data_linear)) args$data_linear <- json_dumps(data_linear)
            do.call(private$do_post, args)
        },

        update_models = function(n_samples, f_forward) {
            private$ensure_in_project()
            n_samp <- 0
            while (n_samp < n_samples)
            {
                n_sim <- trunc((n_samples - n_samp) * 1.2)
                flog.info('Requesting %s field realizations... ...', n_sim)
                fields <- private$do_post('/user/project/models/request_fields', n = n_sim)
                fields <- split(fields, row(fields)) # from matrix to list
                # cat('    stamp:', stamp, '\n')
                flog.debug('    fields: %s x %s', length(fields), length(fields[[1]]))

                flog.info('Running forward model on %s field realizations... ...', length(fields))
                forwards <- lapply(fields, f_forward)
                n_good <- sum(sapply(forwards, function(v) if (all(is.na(v))) 0 else 1))
                flog.info('   %s forward results are invalid', length(forwards) - n_good)

                flog.info('Submitting %s forward results (including invalid ones if any)... ...',
                          length(forwards))
                private$do_post('/user/project/models/submit_forwards',
                          forward_values = json_dumps(do.call(rbind, forwards)))
                # JSON converts R matrix to list of lists ([[...], [...],...]),
                # each row being a member list.
                # `NA` in numerical arrays become `null` in JSON.
                # `null` in JSON becomes 'NA' (string) after transported to server.
                n_samp <- n_samp + n_good
            }

            flog.info('Updating approx to posterior... ...')
            private$do_post('/user/project/models/update')
        },

        summarize_project = function() {
            private$ensure_in_project()
            private$do_get('/user/project/summary')
        },

        simulate_fields = function(n) {
            private$ensure_in_project()
            simulations <- private$do_get('/user/project/request_fields', n=n)
            split(simulations, row(simulations))
        },

        print = function(...) {
            cat('<Session>\n')
            cat('  base_url:', private$base_url, '\n')
            cat('  cookies:\n')
            if (!is.null(private$cookies))
                print(private$cookies)
            cat('  project_id:', private$project_id_, '\n')
        },

        finalize = function() {
            self$logout()
        }
    ),

    active = list(
        # Obtain a vector of the IDs of all the projects in the account
        # that is currently logged in.
        projects = function() {
            private$ensure_logged_in()
            private$do_get('/user/projects')
        },

        # Obtain the `project_id` that is set in the latest call
        # to `set_project`.
        project_id = function() {
            private$project_id_
        }
    ),

    private = list(
        cookies = NULL,
        base_url = 'http://api.anchored-inversion.com',
        project_id_ = NULL,

        do_get = function(url, ...) {
            z <- http_call(httr::GET,
                      url = paste0(private$base_url, url),
                      cookies = private$cookies,
                      query = list(...))
            private$cookies <- z$cookies
            if (is.null(z$value)) invisible() else z$value
        },

        do_post = function(url, ...) {
            z <- http_call(httr::POST,
                      url = paste0(private$base_url, url),
                      cookies = private$cookies,
                      body = list(...))
            private$cookies <- z$cookies
            if (is.null(z$value)) invisible() else z$value
        },

        ensure_logged_in = function() {
            if (is.null(private$cookies)) {
                stop('Please log in first by calling "login_demo()".')
            }
        },

        ensure_in_project = function() {
            if (is.null(private$project_id_)) {
                stop('Please first select a project by calling "set_project(project_id)".')
            }
        }
    ),

    lock_class = TRUE,
    cloneable = FALSE
)

