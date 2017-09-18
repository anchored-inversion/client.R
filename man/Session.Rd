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

    All subsequent method calls on this object
    are in the context of this account,
    until \code{logout} is called.
}

\subsection{\code{set_project(project_id)}}{
    Set the \code{project_id} of the project to work on.

    The project becomes the context of subsequent communications to the server,
    until another call to this method resets the project.

    \describe{
        \item{\code{project_id}}{The ID of the project}
    }
}

\subsection{\code{clear_models()}}{
    Clear all model iterations in this project.
    After this method is called,
    the only thing that remains is the problem definition
    and some settings.

    Typically, this method is called when the user is ready
    to discard all previous work on this project and restart
    from scratch.
}

\subsection{\code{init_models(grid, field_value_range, data_forward, data_linear)}}{
    Initiate the model iteration.

    This method is called in a new project or after \code{clear_models}
    is called on an existing project.

    \describe{
        \item{\code{grid}}{Grid definition. See `grid.R`.}
        \item{\code{field_value_range}}{An estimated range for the field values.}
        \item{\code{data_forward}}{Vector of forward data.}
        \item{\code{data_linear}}{Linear data.}
    }
}

\subsection{\code{update_models(n_samples, f_forward)}}{
    Iterate the model.

    \describe{
        \item{\code{n_samples}}{Number of field realizations to request and use
          for this iteration.}
        \item{\code{f_forward}}{Function of the forward process. This function
          must take one field realization and return a vector of the forward values.}
    }
}

\subsection{\code{summarize_project()}}{
    Obtain a summary of the model iterations of the project.
}

\subsection{\code{simulate_fields(n)}}{
    Request field realizations (i.e. simulations) using
    the latest iteration in the project.

    \describe{
        \item{\code{n}}{Number of field realizations requested.}
    }
}

\subsection{\code{logout()}}{
    Log out of the current user account.
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

