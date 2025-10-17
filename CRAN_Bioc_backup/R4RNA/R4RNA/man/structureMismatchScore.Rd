\name{Structure Mismatch Score}
\alias{structureMismatchScore}
\title{
    Scores how a basepair structure fits a sequence
}
\description{
    Calculates a score that indicates how badly a set of basepairs (i.e. a
    secondary structure) fits with a sequence.  A perfect fit is a structure
    where all basepairs form valid basepairs (A:U, G:C, G:U, and equivalents)
    and has a score of 0.  Each basepair that forms a non-canonical pairing
    or pairs to gaps increases the score by 1, and each base-pair with a
    single-sided gap increases the score by 2.
}
\usage{
    structureMismatchScore(msa, helix, one.gap.penalty = 2, two.gap.penalty = 2,
                 invalid.penalty = 1)
}
\arguments{
    \item{msa}{
        A multiple sequence alignment.  Can be either a \code{Biostrings}
        \code{XStringSet} object or a named array of strings like ones obtained
        from converting XStringSet with \code{as.character}.
    }
    \item{helix}{
        A helix data.frame
    }
    \item{one.gap.penalty}{
        Penalty score for basepairs with one of the bases being a gap
    }
    \item{two.gap.penalty}{
        Penalty score for basepairs with both bases being a gaps
    }
    \item{invalid.penalty}{
        Penalty score for non-canonical basepairs
    }
}
\value{
    Returns an array of mismatch scores.
}
\author{
    Jeff Proctor, Daniel Lai
}
\examples{
    data(helix)
    mismatch <- structureMismatchScore(fasta, known)
    
    # Sort by increasing mismatch
    sorted_fasta <- fasta[order(mismatch)]
}
\keyword{math}
