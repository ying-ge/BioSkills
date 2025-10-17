\name{Helix Type Filters}
\alias{isConflictingHelix}
\alias{isDuplicatingHelix}
\alias{isOverlappingHelix}
\title{
    Logical filters of helix by type
}
\description{
    Given a helix data frame, checks if helices are conflicting, duplicating, or
    overlapping, and returns an array of numeric values, where 0 is FALSE and 1 is
    TRUE.  Values in between 0 and 1 occur when a single helix has multiple
    basepairs with different values, the number observed in this case is the mean
    of the basepair values within the helix. See details for exact definition of
    the three types of events.
}
\usage{
    isConflictingHelix(helix)
    isDuplicatingHelix(helix)
    isOverlappingHelix(helix, query)
}
\arguments{
    \item{helix}{
        A helix data frame
    }
    \item{query}{
        For \code{isOverlappingHelix}, a helix data structure against which
        \code{helix} will be checked for overlap against.
    }
}
\details{
    Helices of length greater than 1 are internally expanded into basepairs of
    length 1, after which the following conditions are evaluated:

    A \bold{conflicting} basepair is one where at least one of its two positions
    is used by either end of another basepair.

    A \bold{duplicating} basepair is one where both of its positions are used
    by both ends of another basepair.

    An \bold{overlapping} basepair is one in \code{helix} where both of its
    positions are used by both ends of another basepair in the \code{query}
    structure.

    In the case of \emph{conflicting} and \emph{duplicating} basepairs, for a
    set of basepairs that satisfies this condition, the basepair situation
    highest on the data frame will be exempt from the condition. i.e.  Say 5
    basepairs are all duplicates of each other, the top 1 will return FALSE,
    while the bottom 4 will return TRUE.  This assumes some significant meaning
    to the ordering of rows prior to using this function.  This is to be used
    with \code{which} to filter out basepairs that satisfy these conditions,
    leaving a set of basepairs free of these events.

    If the original input had helices greater than length 1, then after applying
    all of the above, TRUE is treated as 1, FALSE as 0, and the average of values
    from each basepair is taken as the value for the helix in question.
}
\value{
    Returns an array of numerics corresponding to each row of \code{helix}, giving
    the average conditional status of the helix, where 0 signifying all basepairs
    are FALSE, and 1 where all basepairs are TRUE.
}
\author{
    Daniel Lai
}

\examples{
    data(helix)

    conflicting <- isConflictingHelix(helix)
    duplicating <- isDuplicatingHelix(helix)

    # Nonsensical covariation plot
    plotCovariance(fasta, helix)

    # Plot nonconflicting helices
    plotCovariance(fasta, helix[(!conflicting & !duplicating), ])

    # Similar result
    plotCovariance(fasta, helix, conflict.col = "lightgrey")
}
\keyword{logic}

