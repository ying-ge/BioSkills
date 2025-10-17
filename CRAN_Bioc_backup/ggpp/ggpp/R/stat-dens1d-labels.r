#' @title Replace labels in data based on 1D density
#'
#' @description \code{stat_dens1d_labels()} Sets values mapped to the
#'   \code{label} aesthetic to \code{""} or a user provided character string
#'   based on the local density in regions of a plot panel. Its main use is
#'   together with repulsive geoms from package \code{\link[ggrepel]{ggrepel}}
#'   to restrict labeling to the low density tails of a distribution. By default
#'   the data are handled all together, but it is also possible to control
#'   labeling separately in each tail.
#'
#'   If there is no mapping to \code{label} in \code{data}, the mapping is set
#'   to \code{rownames(data)}, with a message.
#'
#' @details \code{stat_dens1d_labels()} is designed to work together with
#'   geometries from package 'ggrepel'. To avoid text labels being plotted over
#'   unlabelled points the corresponding rows in data need to be retained but
#'   labels replaced with the empty character string, \code{""}. Function
#'   \code{\link{stat_dens1d_filter}} cannot be used with the repulsive geoms
#'   from 'ggrepel' because it drops the observations.
#'
#'   \code{stat_dens1d_labels()} can be useful also in other situations, as the
#'   substitution character string can be set by the user by passing an argument
#'   to \code{label.fill}. If this argument is \code{NULL} the unselected rows
#'   are filtered out.
#'
#'   The local density of observations along \emph{x} or \emph{y} is computed
#'   with function \code{\link[stats]{density}} and used to select observations,
#'   passing to the geom all the rows in its \code{data} input but with with the
#'   text of labels replaced in those "not kept". The default is to select
#'   observations in sparse regions of the plot, but the selection can be
#'   inverted so that only observations in the densest regions are returned.
#'   Specific observations can be protected from having the label replaced by
#'   passing a suitable argument to \code{keep.these}. Logical and integer
#'   vectors function as indexes to rows in \code{data}, while a character
#'   vector is compared to values in the variable mapped to the \code{label}
#'   aesthetic. A function passed as argument to keep.these will receive as
#'   argument the values in the variable mapped to \code{label} and should
#'   return a character, logical or numeric vector as described above.
#'
#'   How many labels are retained intact in addition to those in
#'   \code{keep.these} is controlled with arguments passed to \code{keep.number}
#'   and \code{keep.fraction}. \code{keep.number} sets the maximum number of
#'   observations selected, whenever \code{keep.fraction} results in fewer
#'   observations selected, it is obeyed. If \code{xintercept} is a finite value
#'   within the \emph{x} range of the data and \code{pool.along} is passed
#'   \code{"none"} the data are split into two groups and \code{keep.number} and
#'   \code{keep.fraction} are applied separately to each tail with density still
#'   computed jointly from all observations. If the length of \code{keep.number}
#'   and \code{keep.fraction} is one, half this value is used each tail, if
#'   their length is two, the first value is use for the left tail and the
#'   second value for the right tail (or if using \code{orientation = "y"} the
#'   lower and upper tails, respectively).
#'
#'   Computation of density and of the default bandwidth require at least
#'   two observations with different values. If data do not fulfill this
#'   condition, they are kept only if \code{keep.fraction = 1}. This is correct
#'   behavior for a single observation, but can be surprising in the case of
#'   multiple observations.
#'
#'   Parameters \code{keep.these} and \code{exclude.these} make it possible to
#'   force inclusion or exclusion of labels after the density is computed.
#'   In case of conflict, \code{exclude.these} overrides \code{keep.these}.
#'
#' @note Which points are kept and which not depends on how dense and flexible
#'   is the density curve estimate. This depends on the values passed as
#'   arguments to parameters \code{n}, \code{bw} and \code{kernel}. It is
#'   also important to be aware that both \code{geom_text()} and
#'   \code{geom_text_repel()} can avoid overplotting by discarding labels at
#'   the plot rendering stage, i.e., what is plotted may differ from what is
#'   returned by this statistic.
#'
#' @param mapping The aesthetic mapping, usually constructed with
#'   \code{\link[ggplot2]{aes}} or \code{\link[ggplot2]{aes_}}. Only needs to be
#'   set at the layer level if you are overriding the plot defaults.
#' @param data A layer specific dataset - only needed if you want to override
#'   the plot defaults.
#' @param geom The geometric object to use display the data.
#' @param keep.fraction numeric vector of length 1 or 2 [0..1]. The fraction of
#'   the observations (or rows) in \code{data} to be retained.
#' @param keep.number integer vector of length 1 or 2. Set the maximum number of
#'   observations to retain, effective only if obeying \code{keep.fraction}
#'   would result in a larger number.
#' @param keep.sparse logical If \code{TRUE}, the default, observations from the
#'   more sparse regions are retained, if \code{FALSE} those from the densest
#'   regions.
#' @param keep.these,exclude.these character vector, integer vector, logical
#'   vector or function that takes one or more variables in data selected by
#'   \code{these.target}. Negative integers behave as in R's extraction methods.
#'   The rows from \code{data} indicated by \code{keep.these} and
#'   \code{exclude.these} are kept or excluded irrespective of the local
#'   density.
#' @param these.target character, numeric or logical selecting one or more
#'   column(s) of \code{data}. If \code{TRUE} the whole \code{data} object is
#'   passed.
#' @param pool.along character, one of \code{"none"} or \code{"x"},
#'   indicating if selection should be done pooling the observations along the
#'   \emph{x} aesthetic, or separately on either side of \code{xintercept}.
#' @param xintercept numeric The split point for the data filtering.
#' @param invert.selection logical If \code{TRUE}, the complement of the
#'   selected rows are returned.
#' @param bw numeric or character The smoothing bandwidth to be used. If
#'   numeric, the standard deviation of the smoothing kernel. If character, a
#'   rule to choose the bandwidth, as listed in \code{\link[stats]{bw.nrd}}.
#' @param adjust numeric A multiplicative bandwidth adjustment. This makes it
#'   possible to adjust the bandwidth while still using the a bandwidth
#'   estimator through an argument passed to \code{bw}. The larger the value
#'   passed to \code{adjust} the stronger the smoothing, hence decreasing
#'   sensitivity to local changes in density.
#' @param kernel character See \code{\link{density}} for details.
#' @param n numeric Number of equally spaced points at which the density is to
#'   be estimated for applying the cut point. See \code{\link{density}} for
#'   details.
#' @param return.density logical vector of lenght 1. If \code{TRUE} add columns
#'   \code{"density"} and \code{"keep.obs"} to the returned data frame.
#' @param orientation	character The aesthetic along which density is computed.
#'   Given explicitly by setting orientation to either "x" or "y".
#' @param label.fill character vector of length 1 or a function.
#' @param position The position adjustment to use for overlapping points on this
#'   layer
#' @param show.legend logical. Should this layer be included in the legends?
#'   \code{NA}, the default, includes if any aesthetics are mapped. \code{FALSE}
#'   never includes, and \code{TRUE} always includes.
#' @param inherit.aes If \code{FALSE}, overrides the default aesthetics, rather
#'   than combining with them. This is most useful for helper functions that
#'   define both data and aesthetics and shouldn't inherit behaviour from the
#'   default plot specification, e.g. \code{\link[ggplot2]{borders}}.
#' @param ... other arguments passed on to \code{\link[ggplot2]{layer}}. This
#'   can include aesthetics whose values you want to set, not map. See
#'   \code{\link[ggplot2]{layer}} for more details.
#' @param na.rm	a logical value indicating whether NA values should be stripped
#'   before the computation proceeds.
#'
#' @seealso \code{\link[stats]{density}} used internally.
#'
#' @family statistics returning a subset of data
#'
#' @return A plot layer instance. Using as output \code{data} the input
#'   \code{data} after value substitution based on a 1D the filtering criterion.
#'
#' @examples
#'
#' random_string <-
#'   function(len = 6) {
#'     paste(sample(letters, len, replace = TRUE), collapse = "")
#'   }
#'
#' # Make random data.
#' set.seed(1005)
#' d <- tibble::tibble(
#'   x = rnorm(100),
#'   y = rnorm(100),
#'   group = rep(c("A", "B"), c(50, 50)),
#'   lab = replicate(100, { random_string() })
#' )
#'
#' # using defaults
#' ggplot(data = d, aes(x, y, label = lab)) +
#'   geom_point() +
#'   stat_dens1d_labels()
#'
#' ggrepel.installed <- requireNamespace("ggrepel", quietly = TRUE)
#' if (ggrepel.installed) {
#'   library(ggrepel)
#'
#' # using defaults
#'   ggplot(data = d, aes(x, y, label = lab)) +
#'     geom_point() +
#'     stat_dens1d_labels(geom = "text_repel")
#'
#' # if no mapping to label is found, it is set row names
#'   ggplot(data = d, aes(x, y)) +
#'     geom_point() +
#'     stat_dens1d_labels(geom = "text_repel")
#'
#'   ggplot(data = d, aes(x, y)) +
#'     geom_point() +
#'     stat_dens1d_labels(geom = "text_repel", pool.along = "none")
#'
#'   ggplot(data = d, aes(x, y)) +
#'     geom_point() +
#'     stat_dens1d_labels(geom = "text_repel",
#'                        keep.number = c(0, 10), pool.along = "none")
#'
#'   ggplot(data = d, aes(x, y)) +
#'     geom_point() +
#'     stat_dens1d_labels(geom = "text_repel",
#'                        keep.fraction = c(0, 0.2), pool.along = "none")
#'
#' # using defaults, along y-axis
#'   ggplot(data = d, aes(x, y, label = lab)) +
#'     geom_point() +
#'     stat_dens1d_labels(orientation = "y", geom = "text_repel")
#'
#' # example labelling with coordiantes
#'   ggplot(data = d, aes(x, y, label = sprintf("x = %.2f\ny = %.2f", x, y))) +
#'     geom_point() +
#'     stat_dens1d_filter(colour = "red") +
#'     stat_dens1d_labels(geom = "text_repel", colour = "red", size = 3)
#'
#'   ggplot(data = d, aes(x, y, label = lab, colour = group)) +
#'     geom_point() +
#'     stat_dens1d_labels(geom = "text_repel")
#'
#'   ggplot(data = d, aes(x, y, label = lab, colour = group)) +
#'     geom_point() +
#'     stat_dens1d_labels(geom = "text_repel", label.fill = NA)
#'
#' # we keep labels starting with "a" across the whole plot, but all in sparse
#' # regions. To achieve this we pass as argument to label.fill a fucntion
#' # instead of a character string.
#'   label.fun <- function(x) {ifelse(grepl("^a", x), x, "")}
#'   ggplot(data = d, aes(x, y, label = lab, colour = group)) +
#'     geom_point() +
#'     stat_dens1d_labels(geom = "text_repel", label.fill = label.fun)
#' }
#'
#' # Using geom_debug() we can see that all 100 rows in \code{d} are
#' # returned. But only those labelled in the previous example still contain
#' # the original labels.
#'
#' gginnards.installed <- requireNamespace("gginnards", quietly = TRUE)
#' if (gginnards.installed) {
#'   library(gginnards)
#'
#'   ggplot(data = d, aes(x, y, label = lab)) +
#'     geom_point() +
#'     stat_dens1d_labels(geom = "debug")
#'
#'   ggplot(data = d, aes(x, y, label = lab)) +
#'     geom_point() +
#'     stat_dens1d_labels(geom = "debug", return.density = TRUE)
#'
#'   ggplot(data = d, aes(x, y, label = lab)) +
#'     geom_point() +
#'     stat_dens1d_labels(geom = "debug", label.fill = NULL, return.density = TRUE)
#'
#'   ggplot(data = d, aes(x, y, label = lab)) +
#'     geom_point() +
#'     stat_dens1d_labels(geom = "debug", label.fill = NA, return.density = TRUE)
#'
#'   ggplot(data = d, aes(x, y, label = lab)) +
#'     geom_point() +
#'     stat_dens1d_labels(geom = "debug", label.fill = FALSE, return.density = TRUE)
#' }
#'
#' @export
#'
stat_dens1d_labels <-
  function(mapping = NULL,
           data = NULL,
           geom = "text",
           position = "identity",
           ...,
           keep.fraction = 0.10,
           keep.number = Inf,
           keep.sparse = TRUE,
           keep.these = FALSE,
           exclude.these = FALSE,
           these.target = "label",
           pool.along = c("x", "none"),
           xintercept = 0,
           invert.selection = FALSE,
           bw = "SJ",
           kernel = "gaussian",
           adjust = 1,
           n = 512,
           orientation = c("x", "y"),
           label.fill = "",
           return.density = FALSE,
           na.rm = TRUE,
           show.legend = FALSE,
           inherit.aes = TRUE) {

    pool.along <- rlang::arg_match(pool.along)
    orientation <- rlang::arg_match(orientation)

    if (length(label.fill) > 1L) {
      stop("Length for 'label.fill' is not 0 or 1: ", label.fill)
    }
    if (is.numeric(label.fill)) {
      stop("'label.fill' should not be a 'numeric' value: ", label.fill)
    }
    if (any(is.na(keep.fraction) | keep.fraction < 0 | keep.fraction > 1)) {
      stop("Out of range or missing value for 'keep.fraction': ", keep.fraction)
    }
    if (any(is.na(keep.number) | keep.number < 0)) {
      stop("Out of range or missing value for 'keep.number': ", keep.number)
    }

    ggplot2::layer(
      stat = StatDens1dLabels, data = data, mapping = mapping, geom = geom,
      position = position, show.legend = show.legend, inherit.aes = inherit.aes,
      params = list(na.rm = na.rm,
                    keep.fraction = keep.fraction,
                    keep.number = keep.number,
                    keep.sparse = keep.sparse,
                    keep.these = keep.these,
                    exclude.these = exclude.these,
                    these.target = these.target,
                    pool.along = pool.along,
                    xintercept = xintercept,
                    invert.selection = invert.selection,
                    bw = bw,
                    adjust = adjust,
                    kernel = kernel,
                    n = n,
                    orientation = orientation,
                    label.fill = label.fill,
                    return.density = return.density,
                    ...)
    )
  }

#' @rdname ggpp-ggproto
#' @format NULL
#' @usage NULL
#' @export
StatDens1dLabels <-
  ggplot2::ggproto(
    "StatDens1dLabels",
    ggplot2::Stat,
    compute_panel =
      function(data,
               scales,
               keep.fraction,
               keep.number,
               keep.sparse,
               keep.these,
               exclude.these,
               these.target = "label",
               pool.along,
               xintercept,
               invert.selection,
               bw,
               kernel,
               adjust,
               n,
               orientation,
               label.fill,
               return.density) {

        force(data)
        if (!exists("label", data) && !is.null(label.fill)) {
          data[["label"]] <- rownames(data)
        }

        keep.these <- these2logical(these = keep.these,
                                    data = data,
                                    these.target = these.target)

        exclude.these <- these2logical(these = exclude.these,
                                       data = data,
                                       these.target = these.target)

        # discard redundant splits and make list of logical vectors
        if (pool.along != "x" &&
            xintercept < max(data[[orientation]]) &&
            xintercept > min(data[[orientation]])) {

          selectors <-list(low.tail = data[[orientation]] <= xintercept,
                           high.tail = data[[orientation]] > xintercept)
          if (length(keep.fraction) != 2L) {
            keep.fraction <- rep_len(keep.fraction, length.out = 2)
          }
          if (length(keep.number) != 2L) {
            if (length(keep.number) == 1L) {
              keep.number <- keep.number %/% 2
            }
            keep.number <- rep_len(keep.number, length.out = 2)
          }
          num.rows <- sapply(selectors, sum) # selectors are logical
        } else {
          keep.fraction <- keep.fraction[[1]] # can be a vector or a list
          keep.number <- keep.number[[1]]
          num.rows <- nrow(data)
          selectors <- list(all = rep.int(TRUE, times = num.rows))
        }

        # vectorized
        too.large.frac <- num.rows * keep.fraction > keep.number
        keep.fraction[too.large.frac] <-
          keep.number[too.large.frac] / num.rows[too.large.frac]

        # density on a grid
        # data with fewer than 2 rows is as a special case as density() fails
        if (length(unique(data[[orientation]])) >= 2L) {
          dens <-
            stats::density(data[[orientation]],
                           bw = bw, kernel = kernel, adjust = adjust, n = n,
                           from = scales[[orientation]]$dimension()[1],
                           to = scales[[orientation]]$dimension()[2])

          # estimate density at each observations coordinates
          fdens <- stats::splinefun(dens$x, dens$y) # y contains estimate of density
          dens <- fdens(data[[orientation]])
        } else {
          if (nrow(data) > 1L) {
            message("Density not computed, too few distinct values in '", orientation, "'")
          }
          dens <- rep_len(1, nrow(data))
        }
        # we construct one logical vector by adding observations/label to be kept
        # we may have a list of 1 or 2 logical vectors
        keep <- logical(nrow(data))
        for (i in seq_along(selectors)) {
          if (keep.fraction[i] == 1) {
            keep[ selectors[[i]] ] <- TRUE
          } else if (keep.fraction[i] != 0 && length(selectors[[i]]) >= 2L) {
            if (keep.sparse) {
              keep[ selectors[[i]] ] <-
                dens[ selectors[[i]] ] < stats::quantile(dens[ selectors[[i]] ],
                                                         keep.fraction[i],
                                                         names = FALSE,
                                                         type = 8)
            } else {
              keep[ selectors[[i]] ] <-
                dens[ selectors[[i]] ] >= stats::quantile(dens[ selectors[[i]] ],
                                                          1 - keep.fraction[i],
                                                          names = FALSE,
                                                          type = 8)
            }
          }
        }
        keep <- (keep | keep.these) & !exclude.these

        if (invert.selection){
          keep <- !keep
        }

        if (return.density) {
          data[["keep.obs"]] <- keep
          data[["density"]] <- dens
        }

        if (is.null(label.fill)) {
          data <- data[keep, ]
        } else if (is.function(label.fill)) {
          data[["label"]][!keep] <- label.fill(data[["label"]][!keep])
        } else if (is.na(label.fill)) {
          # NA_logical_, the default NA, cannot always be assigned to character
          label.fill <- NA_character_
          data[["label"]][!keep] <- label.fill
        } else if (is.character(label.fill)) {
          data[["label"]][!keep] <- label.fill
        } else if (is.logical(label.fill)) {
          if (label.fill) {
            data[["label"]][!keep] <- ""
          } # if FALSE data is not modified
        } else {
          stop("'label.fill' is : ", mode(label.fill),
               " instead of 'character' or 'function'.")
        }
        data
      },

    required_aes = "x|y"
  )
