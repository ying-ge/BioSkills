#' Reference points on the margins
#'
#' Margin points can supplement a 2d display with annotations. Margin points
#' can highlight individual values along a margin. The geometries
#' \code{geom_x_margin_point()} and \code{geom_y_margin_point()} behave
#' similarly \code{geom_vline()} and \code{geom_hline()} and share their "double
#' personality" as both annotations and geometries.
#'
#' @param mapping The aesthetic mapping, usually constructed with
#'   \code{\link[ggplot2]{aes}}. Only needs to be set at the layer level if you
#'   are overriding the plot defaults.
#' @param data A layer specific dataset - only needed if you want to override
#'   the plot defaults.
#' @param stat The statistical transformation to use on the data for this layer,
#'   as a string.
#' @param na.rm If \code{FALSE} (the default), removes missing values with a
#'   warning.  If \code{TRUE} silently removes missing values.
#' @param position Position adjustment, either as a string, or the result of a
#'   call to a position adjustment function.
#' @param ... other arguments passed on to \code{\link[ggplot2]{layer}}. This
#'   can include aesthetics whose values you want to set, not map. See
#'   \code{\link[ggplot2]{layer}} for more details.
#' @param show.legend logical. Should this layer be included in the legends?
#'   \code{NA}, the default, includes if any aesthetics are mapped. \code{FALSE}
#'   never includes, and \code{TRUE} always includes.
#' @param inherit.aes If \code{FALSE}, overrides the default aesthetics, rather
#'   than combining with them. This is most useful for helper functions that
#'   define both data and aesthetics and shouldn't inherit behaviour from the
#'   default plot specification, e.g. \code{\link[ggplot2]{borders}}.
#' @param sides A string that controls which sides of the plot the rugs appear
#'   on. It can be set to a string containing any combination of \code{"trbl"},
#'   for top, right, bottom, and left.
#' @param point.shift numeric value expressed in npc units for the shift of the
#'   rug points inwards from the edge of the plotting area.
#' @param xintercept,yintercept numeric Parameters that control the position of
#'   the marginal points. If these are set, data, mapping and show.legend are
#'   overridden.
#'
#' @family Geometries for marginal annotations in ggplots
#'
#' @return A plot layer instance.
#'
#' @export
#'
#' @examples
#' p <- ggplot(mtcars, aes(wt, mpg)) +
#'        geom_point()
#' p
#' p + geom_x_margin_point(xintercept = 3.5)
#' p + geom_y_margin_point(yintercept = c(18, 28, 15))
#' p + geom_x_margin_point(data = data.frame(x = c(2.5, 4.5)),
#'                         mapping = aes(xintercept = x))
#' p + geom_x_margin_point(data = data.frame(x = c(2.5, 4.5)),
#'                         mapping = aes(xintercept = x),
#'                         sides = "tb")
#'
geom_x_margin_point <- function(mapping = NULL, data = NULL,
                     stat = "identity", position = "identity",
                     ...,
                     xintercept,
                     sides = "b",
                     point.shift = 0.017,
                     na.rm = FALSE,
                     show.legend = FALSE,
                     inherit.aes = FALSE) {

  # Act like an annotation
  if (!missing(xintercept)) {
    data <- as.data.frame(list(xintercept = xintercept))
    mapping <- aes(xintercept = xintercept)
    show.legend <- FALSE
  }

  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomXMarginPoint,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      sides = sides,
      point.shift = point.shift,
      na.rm = na.rm,
      ...
    )
  )
}

#' @rdname ggpp-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomXMarginPoint <-
  ggproto("GeomXMarginPoint", Geom,
          required_aes = c("xintercept"),
          non_missing_aes = c("size", "shape", "colour"),
          default_aes = aes(shape = 99, colour = "black", size = 2.5,
                            fill = "red", alpha = NA, stroke = 0.5),

  draw_panel = function(data, panel_params, coord, sides = "b",
                        point.shift = 0.01, na.rm = FALSE) {
    if (is.character(data$shape)) {
      data$shape <- translate_shape_string(data$shape)
    }

    rugpoints <- list()
    data <- coord$transform(data, panel_params)

    # For coord_flip, coord$tranform does not flip the sides where to
    # draw the rugs. We have to flip them.
    flipped <- inherits(coord, 'CoordFlip')
    if (flipped) {
      sides <- chartr('tblr', 'rlbt', sides)
    }

    gp <- gpar(
      col = alpha(data$colour, data$alpha),
      fill = alpha(data$fill, data$alpha),
      # Stroke is added around the outside of the point
      fontsize = data$size * .pt + data$stroke * ggplot2::.stroke / 2,
      lwd = data$stroke * ggplot2::.stroke / 2
    )
    if (!flipped && !is.null(data$xintercept) ) {
      if (grepl("b", sides)) {
        rugpoints$x_b <- pointsGrob(
          x = unit(data$xintercept, "native"),
          y = unit(rep(point.shift, nrow(data)), "npc"),
          pch = ifelse(data$shape == 99, 25, data$shape),
          gp = gp
        )
      }

      if (grepl("t", sides)) {
        rugpoints$x_t <- pointsGrob(
          x = unit(data$xintercept, "native"),
          y = unit(rep(1 - point.shift, nrow(data)), "npc"),
          pch = ifelse(data$shape == 99, 24, data$shape),
          gp = gp
        )
      }
    }

    # needed for handling flipped coords
    if (flipped && !is.null(data$yintercept)) {
      if (grepl("l", sides)) {
        rugpoints$y_l <- pointsGrob(
          x = unit(rep(point.shift, nrow(data)), "npc"),
          y = unit(data$yintercept, "native"),
          pch = ifelse(data$shape == 99, 23, data$shape),
          gp = gp
        )
      }
      if (grepl("r", sides)) {
        rugpoints$y_r <- pointsGrob(
          x = unit(rep(1 - point.shift, nrow(data)), "npc"),
          y = unit(data$yintercept, "native"),
          pch = ifelse(data$shape == 99, 23, data$shape),
          gp = gp
        )
      }
    }

    gTree(children = do.call("gList", rugpoints))
  },


  draw_key = ggplot2:::draw_key_point
)

## y axis marging
#' @rdname geom_x_margin_point
#' @export
#'
geom_y_margin_point <- function(mapping = NULL, data = NULL,
                                 stat = "identity", position = "identity",
                                 ...,
                                 yintercept,
                                 sides = "l",
                                 point.shift = 0.017,
                                 na.rm = FALSE,
                                 show.legend = FALSE,
                                 inherit.aes = FALSE) {

  # Act like an annotation
  if (!missing(yintercept)) {
    data <- as.data.frame(list(yintercept = yintercept))
    mapping <- aes(yintercept = yintercept)
    show.legend <- FALSE
  }

  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomYMarginPoint,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      sides = sides,
      point.shift = point.shift,
      na.rm = na.rm,
      ...
    )
  )
}

#' @rdname ggpp-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomYMarginPoint <-
  ggproto("GeomYMarginPoint", Geom,
          required_aes = c("yintercept"),
          non_missing_aes = c("size", "shape", "colour"),
          default_aes = aes(shape = 99, colour = "black", size = 2.5,
                            fill = "red", alpha = NA, stroke = 0.5),

          draw_panel = function(data, panel_params, coord, sides = "l",
                                point.shift = 0.01, na.rm = FALSE) {
            if (is.character(data$shape)) {
              data$shape <- translate_shape_string(data$shape)
            }

            rugpoints <- list()
            data <- coord$transform(data, panel_params)

            # For coord_flip, coord$tranform does not flip the sides where to
            # draw the rugs. We have to flip them.
            flipped <- inherits(coord, 'CoordFlip')
            if (flipped) {
              sides <- chartr('tblr', 'rlbt', sides)
            }

            gp <- gpar(
              col = alpha(data$colour, data$alpha),
              fill = alpha(data$fill, data$alpha),
              # Stroke is added around the outside of the point
              fontsize = data$size * .pt + data$stroke * ggplot2::.stroke / 2,
              lwd = data$stroke * ggplot2::.stroke / 2
            )
            if (!flipped && !is.null(data$yintercept)) {
              if (grepl("l", sides)) {
                rugpoints$y_l <- pointsGrob(
                  x = unit(rep(point.shift, nrow(data)), "npc"),
                  y = unit(data$yintercept, "native"),
                  pch = ifelse(data$shape == 99, 23, data$shape),
                  gp = gp
                )
              }
              if (grepl("r", sides)) {
                rugpoints$y_r <- pointsGrob(
                  x = unit(rep(1 - point.shift, nrow(data)), "npc"),
                  y = unit(data$yintercept, "native"),
                  pch = ifelse(data$shape == 99, 23, data$shape),
                  gp = gp
                )
              }
            }

            # needed for handling flipped coords
            if (flipped && !is.null(data$xintercept)) {
              if (grepl("b", sides)) {
                rugpoints$x_b <- pointsGrob(
                  x = unit(data$xintercept, "native"),
                  y = unit(rep(point.shift, nrow(data)), "npc"),
                  pch = ifelse(data$shape == 99, 25, data$shape),
                  gp = gp
                )
              }

              if (grepl("t", sides)) {
                rugpoints$x_t <- pointsGrob(
                  x = unit(data$xintercept, "native"),
                  y = unit(rep(1 - point.shift, nrow(data)), "npc"),
                  pch = ifelse(data$shape == 99, 24, data$shape),
                  gp = gp
                )
              }
            }

            gTree(children = do.call("gList", rugpoints))
          },


          draw_key = ggplot2:::draw_key_point
  )
