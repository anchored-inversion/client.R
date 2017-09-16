#' @importFrom futile.logger flog.debug flog.info flog.warn flog.error
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
    cookies <- httr::cookies(z)
    value <- httr::content(z)
    if (!is.null(value)) value <- json_loads(value)
    list(cookies = cookies, value = value)
}


#' Class `Session`.
#'
#' This is the main entry point that the user (the "client")
#' uses to communicate with the Anchored Inversion server.
#' The user would create a new `Session` object by `Session$new()`,
#' then call the object's methods to communicate with the server.
#'
#' This uses the OOP style of programming provided by the package
#' `R6`. This is an OOP system that is very different from
#' R's built-in `S3` and `S4` systems. It is close to the OOP
#' style of mainstream languages such as `Python`, `Java`, and `C++`.
#' This style is convenient for the object to persist connection
#' states as the client-server communication progresses.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords data
#' @export
#'
#' @section Methods:
#' \describe{
#'     \item{\code{new(base_url)}}{%
#'         Create a new object of this class.
#'         @param base_url Base URL portion of the API server.
#'            The general user should ignore this.}
#'     \item{\code{login_demo()}}{Log into demo account.}
#'     \item{\code{set_project(project_id)}}{%
#'         Set the `project_id` of the project to work on.
#'         The project becomes the context of subsequent communications,
#'         and does not need to be specified again until changed by
#'         another call to this method.}
#'     \item{\code{clear_models()}}{%
#'         Clear all model iterations in this project.
#'         The only thing that remains is the problem definition
#'         and certain settings. Afterwards the user would
#'         restart model iterations from ground zero.}
#'     \item{\code{init_models(...)}}{%
#'         Initiate the model iteration.
#'         Call this method in a new project or after calling `clear_models`
#'         on an existing project.
#'
#'         @param grid Grid definition. See `grid.R`.
#'         @param field_value_range An estimated range for the field values.
#'         @param data_forward Vector of forward data.
#'         @param data_linear Linear data.}
#'     \item{\code{update_models(...)}}{%
#'         Iterate the model.
#'
#'         @param n_samples Number of field realizations to request and use
#'           for this iteration.
#'         @param f_forward Function of the forward process. This function
#'           must take one field realization and return a vector of the forward values.}
#'     \item{\code{summarize_project(...)}}{%
#'         Obtain a summary of the model iterations of the project.}
#'     \item{\code{simulate_fields(n)}}{%
#'         Request field realizations (i.e. simulations) using
#'         the latest iteration in the project.}
#'     \item{\code{logout()}}{}
#' }
#'
#' @section Properties:
#' \describe{
#'     \item{\code{projects()}}{%
#'         Obtain a vector of the IDs of all the projects in the account
#'         that is currently logged in.}
#'     \item{\code{project_id()}}{%
#'         Obtain the `project_id` that is set in the latest call
#'         to `set_project`.}
#' }
Session <- R6::R6Class("Session",
    public = list(
        initialize = function(base_url = NULL) {
            # TODO: remove trailing '/' in parameter `base_url`.
            if (!is.null(base_url)) private$base_url <- base_url
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
        #' Obtain a vector of the IDs of all the projects in the account
        #' that is currently logged in.
        projects = function() {
            private$ensure_logged_in()
            private$do_get('/user/projects')
        },

        #' Obtain the `project_id` that is set in the latest call
        #' to `set_project`.
        project_id = function() {
            private$project_id_
        }
    ),

    private = list(
        cookies = NULL,
        base_url = 'http://localhost:8000',
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

