\name{Session}
\alias{Session}
\title{Class `Session`.}

\description{
This is the main entry point that the user (the "client")
uses to communicate with the Anchored Inversion server.
The user would create a new `Session` object by `Session$new()`,
then call the object's methods to communicate with the server.
}

\details{
This uses the OOP style of programming provided by the package
`R6`. This is an OOP system that is very different from
R's built-in `S3` and `S4` systems. It is close to the OOP
style of mainstream languages such as `Python`, `Java`, and `C++`.
This style is convenient for the object to persist connection
states as the client-server communication progresses.

Some methods of this class are simple wrappers of corresponding functions
in the public API. Please see the API reference for more details.
}

\section{Methods}{

\subsection{\code{new(base_url)}}{
    Instantiate the class.

    \emph{Arguments:}

    \describe{
        \item{\code{base_url}}{Base URL portion of the API server.
          The general user should ignore this.}
    }

    \emph{Value:}

    A new instance of the \code{Session} class.
}

\subsection{\code{login_demo}}{
    Log into demo account.

    See API reference for details.
}

\subsection{\code{set_project(project_id)}}{
    Set the \code{project_id} of the project to work on.

    See API reference for details.
}

\subsection{\code{clear_models()}}{
    Clear all model iterations in this project.

    See API reference for details.
}

\subsection{\code{init_models(grid, field_value_range, data_forward, data_linear)}}{
    Initiate the model iteration.

    See API reference for details.
}

\subsection{\code{update_models(n_samples, f_forward)}}{
    Iterate the model.

    \describe{
        \item{\code{n_samples}}{Number of field realizations to request and use
          for this iteration.}
        \item{\code{f_forward}}{Function of the forward process. This function
          must take one field realization and return a vector of the forward values.}
    }

    See `request_fields`, `submit_forwards`, and `update` in API reference
    for more details.
}

\subsection{\code{summarize_project()}}{
    Obtain a summary of the model iterations of the project.
}

\subsection{\code{simulate_fields(n)}}{
    Request field realizations (i.e. simulations).

    See API reference for details.
}

\subsection{\code{logout()}}{
    Log out of the current user account.

    See API reference for details.
}
}


\section{Properties}{

\subsection{\code{projects()}}{
    Obtain a vector of the IDs of all the projects in the account
    that is currently logged in.
}

\subsection{\code{project_id()}}{
    Obtain the \code{project_id} that is set in the latest call
    to \code{set_project}.
}
}

