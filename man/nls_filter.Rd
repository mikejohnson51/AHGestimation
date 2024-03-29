% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filters.R
\name{nls_filter}
\alias{nls_filter}
\title{Implements NLS filtering}
\usage{
nls_filter(df, allowance = 0.5)
}
\arguments{
\item{df}{a data.frame with at least a Q and one other AHG field (Y. TW, V).}

\item{allowance}{how much deviation from observed should 
be allowed (default = .5)}
}
\value{
data.frame
}
\description{
An NLS fit provides the best relation by relation fit. 
For each provided relationship, an NLS fit is computed and used to estimate 
the predicted {V,TW,Y} for a given Q. If the actual value is outside 
the specified allowance it is removed.
}
\seealso{
Other filters: 
\code{\link{date_filter}()},
\code{\link{mad_filter}()},
\code{\link{qva_filter}()},
\code{\link{significance_check}()}
}
\concept{filters}
