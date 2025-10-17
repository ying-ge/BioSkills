#' @title Linked Text
#' @include ggp2-margins.R utilities.R ggpp-legend-draw.R
#'
#' @description Linked text geometries are most useful for adding data labels to
#'   plots. `geom_text_s()` and `geom_label_s()` add text to the plot and for
#'   nudged positions link the original location to the nudged text with a
#'   segment or arrow.
#'
#' @details Geometries \code{geom_text_s()} and \code{geom_label_s()} have an
#'   interface similar to that of \code{\link[ggplot2]{geom_text}} and
#'   \code{\link[ggplot2]{geom_label}}, but support additional features.
#'   Similarly to \code{geom_text_repel()} and \code{geom_label_repel()} when
#'   used together with position functions defined in package 'ggpp' they draw a
#'   segment linking the label at a displaced position to the original position,
#'   usually a point corresponding to an observation to which the label refers.
#'   Another difference is that they allow control of to which graphical
#'   elements the mappings to colour and alpha aesthetics are applied.
#'   Differently to \code{geom_label()}, \code{geom_label_s()} obeys aesthetic
#'   mappings to \code{linewidth} and \code{linetype} applied to the line at the
#'   edge of the label box. These features are reflected in the plot key, except
#'   for the segment, assumed not to be used to display information
#'   only in coordination with other graphic elements.
#'
#'   In \code{geom_label_s()} the default \code{fill} is similar to
#'   \code{"white"} but with its \code{alpha} component set to 0.75. This
#'   differs from \code{"white"} used in \code{geom_label()}: the default fill
#'   is semitransparent with the intention that accidental occlusion of
#'   observations is obvious irrespective of the order in which layers are added
#'   to the plot.
#'
#'   Layer functions \code{geom_text_s()} and \code{geom_label_s()} use by
#'   default \code{\link{position_nudge_keep}} which is backwards compatible
#'   with \code{\link[ggplot2]{position_nudge}}. In contrast to
#'   \code{\link[ggplot2]{position_nudge}}, \code{\link{position_nudge_keep}}
#'   and all other position functions defined in packages 'ggpp' and 'ggrepel'
#'   keep the original coordinates, thus allowing the plotting of connecting
#'   segments and arrows.
#'
#'   Differently to \code{geom_text_repel()} and \code{geom_label_repel()},
#'   \code{geom_text_s()} and \code{geom_label_s()} do not make use of
#'   additional aesthetics for the segments or boxes, but instead allow the
#'   choice of which elements are targeted by the aesthetics and which are
#'   rendered in a default colour. In the grammar of graphics using the same
#'   aesthetic with multiple meanings is not allowed, thus, the approach used in
#'   the geometry layer functions from package 'ggpp' attempts to enforce this.
#'
#' @section Plot boundaries and clipping: Note that when you change the scale
#'   limits for \emph{x} and/or \emph{y} of a plot, text labels stay the same
#'   size, as determined by the \code{size} aesthetic, given in millimetres. The
#'   actual size as seen in the plotted output is decided during the rendering
#'   of the plot to a graphics device. Limits are expanded only to include the
#'   anchor point of the labels because the "width" and "height" of a text
#'   element are 0 (as seen by ggplot2). Text labels do have height and width,
#'   but in grid units, not data units.
#'
#' @section Alignment: You can modify text alignment with the \code{vjust} and
#'   \code{hjust} aesthetics. These can either be a number between 0
#'   (right/bottom) and 1 (top/left) or a character (\code{"left"},
#'   \code{"middle"}, \code{"right"}, \code{"bottom"}, \code{"center"},
#'   \code{"top"}). In addition, you can use special alignments for
#'   justification including \code{"position"}, \code{"inward"} and
#'   \code{"outward"}. Inward always aligns text towards the center of the
#'   plotting area, and outward aligns it away from the center of the plotting
#'   area. If tagged with \code{_mean} or \code{_median} (e.g.,
#'   \code{"outward_mean"}) the mean or median of the data in the panel along
#'   the corresponding axis is used as center. If the characters following the
#'   underscore represent a number (e.g., \code{"outward_10.5"}) the reference
#'   point will be this value in data units. Position justification is computed
#'   based on the direction of the displacement of the position of the label so
#'   that each individual text or label is justified outwards from its original
#'   position. The default justification is \code{"position"}.
#'
#'   If no position displacement is applied, or a position function defined in
#'   'ggplot2' is used, these geometries behave similarly to the corresponding
#'   ones from package 'ggplot2' with a default justification of \code{0.5} and
#'   no segment drawn.
#'
#' @section Differences from earlier versions:
#'   The user interface is for the most part stable starting from 'ggpp' (==
#'   0.5.7). In 'ggpp' (== 0.5.0) support for aesthetics related to segments was
#'   removed, and replaced by parameters and a new mechanism for targeting the
#'   usual \code{colour} and \code{alpha} aesthetics to text, border, and
#'   segment.
#'
#' @param mapping Set of aesthetic mappings created by
#'   \code{\link[ggplot2]{aes}}. If specified and with \code{inherit.aes = TRUE}
#'   (the default), it is combined with the default mapping at the top level of
#'   the plot. You only need to supply \code{mapping} if there isn't a mapping
#'   defined for the plot.
#' @param data A data frame. If specified, overrides the default data frame
#'   defined at the top level of the plot.
#' @param stat The statistical transformation to use on the data for this layer,
#'   as a string.
#' @param position Position adjustment, either as a string, or the result of a
#'   call to a position adjustment function.
#' @param parse If \code{TRUE}, the labels will be parsed into expressions and
#'   displayed as described in \code{?plotmath}.
#' @param na.rm If \code{FALSE} (the default), removes missing values with a
#'   warning.  If \code{TRUE} silently removes missing values.
#' @param show.legend logical. Should this layer be included in the legends?
#'   \code{NA}, the default, includes a legend if any aesthetics are mapped.
#'   \code{FALSE} never includes it, and \code{TRUE} always includes it.
#' @param inherit.aes If \code{FALSE}, overrides the default aesthetics, rather
#'   than combining them. This is most useful for helper functions that define
#'   both data and aesthetics and shouldn't inherit behaviour from the default
#'   plot specification, e.g., \code{\link[ggplot2]{borders}}.
#' @param check_overlap If \code{TRUE}, text that overlaps previous text in the
#'   same layer will not be plotted. \code{check_overlap} takes place at draw
#'   time and in the order of the data, thus its action depends of the size at
#'   which the plot is drawn.
#' @param size.unit How the `size` aesthetic is interpreted: as millimetres
#'   (`"mm"`, default), points (`"pt"`), centimetres (`"cm"`), inches (`"in"`),
#'   or picas (`"pc"`).
#' @param ... other arguments passed on to \code{\link[ggplot2]{layer}}. There
#'   are three types of arguments you can use here:
#'
#'   \itemize{ \item Aesthetics: to set an aesthetic to a fixed value, like
#'   \code{colour = "red"} or \code{size = 3}. \item Other arguments to the
#'   layer, for example you override the default \code{stat} associated with the
#'   layer. \item Other arguments passed on to the stat. }
#' @param nudge_x,nudge_y Horizontal and vertical adjustments to nudge the
#'   starting position of each text label. The units for \code{nudge_x} and
#'   \code{nudge_y} are the same as for the data units on the x-axis and y-axis.
#' @param default.colour,default.color A colour definition to use for elements
#'   not targeted by the colour aesthetic.
#' @param colour.target,color.target A vector of character strings; \code{"all"},
#'   \code{"text"}, \code{"segment"}, \code{"box"}, \code{"box.line"}, and
#'   \code{"box.fill"} or \code{"none"}.
#' @param default.alpha numeric in [0..1] A transparency value to use for
#'   elements not targeted by the alpha aesthetic.
#' @param alpha.target A vector of character strings; \code{"all"},
#'   \code{"text"}, \code{"segment"}, \code{"box"}, \code{"box.line"}, and
#'   \code{"box.fill"} or \code{"none"}.
#' @param add.segments logical Display connecting segments or arrows between
#'   original positions and displaced ones if both are available.
#' @param box.padding,point.padding numeric By how much each end of the segments
#'   should shortened in mm.
#' @param segment.linewidth numeric Width of the segments or arrows in mm.
#' @param min.segment.length numeric Segments shorter that the minimum length
#'   are not rendered, in mm.
#' @param arrow specification for arrow heads, as created by
#'   \code{\link[grid]{arrow}}
#'
#' @section Aesthetics: Layer functions \code{geom_text_s()} and
#'   \code{geom_label_s()} require aesthetics \code{x}, \code{y} and
#'   \code{label} and support aesthetics: \code{alpha}, \code{colour},
#'   \code{group}, \code{size} (of text), \code{family}, \code{fontface},
#'   \code{lineheight}, \code{hjust} and \code{vjust}. In addition,
#'   \code{geom_text_s} supports \code{angle} and \code{geom_label_s} supports
#'   \code{fill}, \code{linewidth} and \code{linetype}. See
#'   \code{\link[ggplot2]{aes_colour_fill_alpha}},
#'   \code{\link[ggplot2]{aes_linetype_size_shape}},
#'   \code{\link[ggplot2]{aes_position}}, and
#'   \code{\link[ggplot2]{aes_group_order}}.
#'
#'   In 'ggplot2' \code{linewidth} when applied to the border of the box drawn
#'   by \code{geom_label()} is given in points rather than in mm because of a
#'   historical error in the code. In other geometries such as
#'   \code{geom_segment()} \code{linewidth} is given in mm. As in
#'   \code{geom_label_s()} it is important to remain consistent among
#'   different \code{linewidth} specifications, mm are used both for the box
#'   border and linking segment. To imitate the behaviour of \code{geom_label()}
#'   a correction factor of 0.75 (more exactly 1 pt = 0.7528 mm) can be used for
#'   the line width of the border of the box.
#'
#' @section Position functions: Many layer functions from package 'ggpp' are
#'   designed to work seamlessly with position functions that keep, rather than
#'   discard, the original \code{x} and \code{y} positions in \code{data} when
#'   computing a new displaced position. See \code{\link{position_nudge_keep}},
#'   \code{\link{position_dodge_keep}}, \code{\link{position_jitter_keep}},
#'   \code{\link{position_nudge_center}}, \code{\link{position_nudge_line}},
#'   \code{\link{position_nudge_to}}, \code{\link{position_dodgenudge}},
#'   \code{\link{position_jitternudge}}, and \code{\link{position_stacknudge}}
#'   for examples and details of their use.
#'
#' @seealso \code{\link[ggplot2]{geom_text}}, \code{\link[ggplot2]{geom_label}}
#'   and other documentation of package 'ggplot2'.
#'
#' @return A plot layer instance.
#'
#' @export
#'
#' @examples
#'
#' my.cars <- mtcars[c(TRUE, FALSE, FALSE, FALSE), ]
#' my.cars$name <- rownames(my.cars)
#'
#' # no nudging
#' ggplot(my.cars, aes(wt, mpg, label = name)) +
#'   geom_text_s() +
#'   expand_limits(x = c(2, 6))
#'
#' # base plot
#' p <- ggplot(my.cars, aes(wt, mpg, label = name)) +
#'        geom_point()
#'
#' # Using nudging
#' p +
#'   geom_text_s(nudge_x = 0.12) +
#'   expand_limits(x = 6.2)
#' p +
#'   geom_text_s(nudge_x = -0.12) +
#'   expand_limits(x = 1.5)
#' p +
#'   geom_text_s(nudge_x = 0.12,
#'               arrow = arrow(length = grid::unit(1.5, "mm")),
#'               point.padding = 0.4) +
#'   expand_limits(x = 6.2)
#' p +
#'   geom_text_s(nudge_y = 0.1, nudge_x = 0.07) +
#'   expand_limits(x = 6.2)
#' p +
#'   geom_text_s(nudge_y = 1, angle = 90) +
#'   expand_limits(y = 30)
#' p +
#'   geom_text_s(angle = 90, nudge_y = 1,
#'               arrow = arrow(length = grid::unit(1.5, "mm")),
#'               colour.target = "segment", colour = "red") +
#'   expand_limits(y = 30)
#' p +
#'   geom_text_s(aes(colour = factor(cyl)),
#'               angle = 90, nudge_y = 1,
#'               arrow = arrow(length = grid::unit(1.5, "mm")),
#'               alpha.target = "segment", alpha = 0.3) +
#'   expand_limits(y = 30)
#'
#' p +
#'   geom_label_s(nudge_x = 0.12) +
#'   expand_limits(x = 6.2)
#' p +
#'   geom_label_s(nudge_x = 0.12, linetype = "dotted", linewidth = 0.3) +
#'   expand_limits(x = 6.2)
#' p +
#'   geom_label_s(aes(colour = factor(cyl)),
#'                nudge_x = 0.12,
#'                colour.target = "box",
#'                linewidth = 0.5,
#'                label.r = unit(0, "lines")) +
#'   expand_limits(x = 6.2)
#' p +
#'   geom_label_s(nudge_x = 0.12, linewidth = 0) +
#'   expand_limits(x = 6.2)
#'
#' # No segments
#' p +
#'   geom_label_s(nudge_x = 0.05, segment.linewidth = 0) +
#'   expand_limits(x = 6.2)
#'
#' # Nudging away from arbitrary point
#' p +
#'   geom_label_s(hjust = "outward_1", nudge_x = 0.12) +
#'   expand_limits(x = 6.2)
#' p +
#'   geom_label_s(hjust = "inward_3", nudge_y = 0.4)
#'
#' p +
#'   geom_label_s(nudge_y = 1, angle = 90) +
#'   expand_limits(y = 30)
#'
#' # Add aesthetic mappings and adjust arrows
#' p +
#'   geom_text_s(aes(colour = factor(cyl)),
#'               angle = 90,
#'               nudge_y = 1,
#'               arrow = arrow(angle = 20,
#'                             length = grid::unit(1.5, "mm"),
#'                             ends = "first",
#'                             type = "closed")) +
#'   scale_colour_discrete(l = 40) + # luminance, make colours darker
#'   expand_limits(y = 27)
#'
#' p +
#'   geom_text_s(aes(colour = factor(cyl)),
#'               angle = 90,
#'               nudge_y = 1,
#'               arrow = arrow(angle = 20,
#'                             length = grid::unit(1.5, "mm"),
#'                             ends = "first",
#'                             type = "closed")) +
#'   scale_colour_discrete(l = 40) + # luminance, make colours darker
#'   expand_limits(y = 27)
#'
#' p +
#'   geom_label_s(aes(colour = factor(cyl)),
#'               colour.target = c("box", "text"),
#'               nudge_x = 0.3,
#'               arrow = arrow(angle = 20,
#'                             length = grid::unit(1/3, "lines"))) +
#'   scale_colour_discrete(l = 40) + # luminance, make colours darker
#'   expand_limits(x = 7)
#'
#' p +
#'   geom_label_s(aes(colour = factor(cyl)),
#'               nudge_x = 0.3,
#'               colour.target = c("box", "segment"),
#'               linewidth = 0.5,
#'               arrow = arrow(angle = 20,
#'                             length = grid::unit(1/3, "lines"))) +
#'   scale_colour_discrete(l = 40) + # luminance, make colours darker
#'   expand_limits(x = 7)
#'
#' p +
#'   geom_label_s(aes(colour = factor(cyl), fill = factor(cyl)),
#'               nudge_x = 0.3,
#'               alpha.target = "box",
#'               alpha = 0.1,
#'               linewidth = 0.5,
#'               arrow = arrow(angle = 20,
#'                             length = grid::unit(1/3, "lines"))) +
#'   scale_colour_discrete(l = 40) + # luminance, make colours darker
#'   expand_limits(x = 7)#' # Scale height of text, rather than sqrt(height)
#'
#' p +
#'   geom_text_s(aes(size = wt), nudge_x = -0.1) +
#'   scale_radius(range = c(3,6)) + # override scale_area()
#'     expand_limits(x = c(1.8, 5.5))
#'
geom_text_s <- function(mapping = NULL,
                        data = NULL,
                        stat = "identity",
                        position = "identity",
                        ...,
                        parse = FALSE,
                        nudge_x = 0,
                        nudge_y = 0,
                        default.colour = "black",
                        default.color = default.colour,
                        colour.target = "text",
                        color.target = colour.target,
                        default.alpha = NA,
                        alpha.target = "all",
                        add.segments = TRUE,
                        box.padding = 0.25,
                        point.padding = 1e-06,
                        segment.linewidth = 0.5,
                        min.segment.length = 0,
                        arrow = NULL,
                        check_overlap = FALSE,
                        size.unit = "mm",
                        na.rm = FALSE,
                        show.legend = NA,
                        inherit.aes = TRUE) {

  colour.target <-
    rlang::arg_match(color.target,
                     values = c("all", "text", "segment", "none"),
                     multiple = TRUE)
  alpha.target <-
    rlang::arg_match(alpha.target,
                     values = c("all", "text", "segment", "none"),
                     multiple = TRUE)


  if (!missing(nudge_x) || !missing(nudge_y)) {
    if (!missing(position) && position != "identity") {
      rlang::abort("You must specify either `position` or `nudge_x`/`nudge_y`.")
    }
    # original position needed for "position" justification
    position <-
      position_nudge_center(nudge_x, nudge_y, kept.origin = "original")
  }

  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomTextS,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      parse = parse,
      default.colour = default.color,
      colour.target = colour.target,
      default.alpha = default.alpha,
      alpha.target = alpha.target,
      add.segments = add.segments,
      box.padding = box.padding,
      point.padding = point.padding,
      segment.linewidth = segment.linewidth,
      min.segment.length = min.segment.length,
      arrow = arrow,
      check_overlap = check_overlap,
      size.unit = size.unit,
      na.rm = na.rm,
      ...
    )
  )
}

#' @rdname ggpp-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomTextS <-
  ggplot2::ggproto("GeomTextS", ggplot2::Geom,
                   required_aes = c("x", "y", "label"),

                   non_missing_aes = "angle",

                   default_aes = ggplot2::aes(
                     colour = "black",
                     size = 3.88,
                     angle = 0,
                     hjust = "position",
                     vjust = "position",
                     alpha = NA,
                     family = "",
                     fontface = 1,
                     lineheight = 1.2
                   ),

                   draw_panel = function(data,
                                         panel_params,
                                         coord, #panel_scales,
                                         parse = FALSE,
                                         size.unit = "mm",
                                         default.colour = "black",
                                         colour.target = "all",
                                         default.alpha = NA,
                                         alpha.target = "all",
                                         na.rm = FALSE,
                                         check_overlap = FALSE,
                                         add.segments = TRUE,
                                         box.padding = 0.25,
                                         point.padding = 1e-06,
                                         segment.linewidth = 0.5,
                                         min.segment.length = 0,
                                         arrow = NULL) {

                     add.segments <- add.segments && any(c("x_orig", "y_orig") %in% colnames(data))

                     data$label <- as.character(data$label)
                     data <- subset(data, !is.na(label) & label != "")
                     if (nrow(data) == 0L) {
                       return(nullGrob())
                     }

                     lab <- data$label
                     if (parse) {
                       lab <- parse_safe(lab)
                     }

                     data <- coord$transform(data, panel_params)

                     if (all(c("x_orig", "y_orig") %in% colnames(data))) {
                       data_orig <- data.frame(x = data$x_orig, y = data$y_orig)
                       data_orig <- coord$transform(data_orig, panel_params)
                       data$x_orig <- data_orig$x
                       data$y_orig <- data_orig$y
                     }

                     if (is.character(data$vjust)) {
                       data$vjust <-
                         compute_just2d(data = data,
                                        coord = coord,
                                        panel_params = panel_params,
                                        just = data$vjust,
                                        a = "y", b = "x")
                     }
                     if (is.character(data$hjust)) {
                       data$hjust <-
                         compute_just2d(data = data,
                                        coord = coord,
                                        panel_params = panel_params,
                                        just = data$hjust,
                                        a = "x", b = "y")
                     }

                     size.unit <- resolve_text_unit(size.unit)

                     if (add.segments) {
                       segments.data <-
                         shrink_segments(data,
                                         point.padding = point.padding,
                                         box.padding = box.padding,
                                         min.segment.length = min.segment.length)
                     }
                     # loop needed as gpar is not vectorized
                     all.grobs <- grid::gList()

                     for (row.idx in 1:nrow(data)) {
                       row <- data[row.idx, , drop = FALSE]
                       text.alpha <-
                         ifelse(any(alpha.target %in% c("all", "text")),
                                row$alpha, default.alpha)
                       segment.alpha <-
                         ifelse(any(alpha.target %in% c("all", "segment")),
                                row$alpha, default.alpha)
                       user.grob <- grid::textGrob(
                         lab[row.idx],
                         row$x, row$y, default.units = "native",
                         hjust = row$hjust, vjust = row$vjust,
                         rot = row$angle,
                         gp = gpar(
                           col = ifelse(any(colour.target %in% c("all", "text")),
                                        ggplot2::alpha(row$colour, text.alpha),
                                        ggplot2::alpha(default.colour, text.alpha)),
                           fontsize = row$size * size.unit,
                           fontfamily = row$family,
                           fontface = row$fontface,
                           lineheight = row$lineheight
                         ),
                         check.overlap = check_overlap
                       )

                       # give unique name to each grob
                       user.grob$name <- paste("text.s.grob", row$group, row.idx, sep = ".")

                       if (add.segments) {
                         segment.row <- segments.data[row.idx, , drop = FALSE]
                         if (segment.row$too.short) {
                           segment.grob <- grid::nullGrob()
                         } else {
                           segment.grob <-
                             grid::segmentsGrob(x0 = segment.row$x,
                                                y0 = segment.row$y,
                                                x1 = segment.row$x_orig,
                                                y1 = segment.row$y_orig,
                                                arrow = arrow,
                                                gp = grid::gpar(
                                                  col = if (segment.linewidth == 0) NA else # lwd = 0 is invalid in 'grid'
                                                    ifelse(any(colour.target %in% c("all", "segment")),
                                                           ggplot2::alpha(row$colour, segment.alpha),
                                                           ggplot2::alpha(default.colour, segment.alpha)),
                                                  lwd = (if (segment.linewidth == 0) 0.5 else segment.linewidth) * ggplot2::.stroke),
                                                name = paste("text.s.segment", row$group, row.idx, sep = "."))
                         }
                         all.grobs <- grid::gList(all.grobs, segment.grob, user.grob)
                       } else {
                         all.grobs <- grid::gList(all.grobs, user.grob)
                       }
                     }

                     # name needs to be unique within plot, so we would need to know other layers
#                     grid::grobTree(children = all.grobs, name = "geom.text.s.panel")
                     grid::grobTree(children = all.grobs)

                   },

                   draw_key = draw_key_text_s
  )

# heavily modified from geom-text.r from 'ggplot2' 3.1.0
# when just is "outward" or "inward" or one of its variants and the geom
# supports the angle aesthetics we need to take into account that justification
# is relative to the text, rather than the plot axes. By using compute_split()
# we add support for definitions of "inward" and "outward" relative to
# arbitrary positions along the axis.
#
# We support "position" (could be called "away") when nudging or other
# displacement has been applied and the original position saved.
#
# This function can handle either hjust or vjust, but only one at a time.
compute_just2d <- function(data,
                           coord,
                           panel_params,
                           just,
                           a = "x",
                           b = a) {
  just <- ifelse(just == "away", "position", just)
  if (a != b) {
    angle <- data$angle
  } else {
    angle <- 0
  }
  if (any(grepl("outward|inward|position", just))) {
    # ensure all angles are in -360...+360
    angle <- angle %% 360
    # ensure correct behaviour for angles in -360...+360
    angle <- ifelse(angle > 180, angle - 360, angle)
    angle <- ifelse(angle < -180, angle + 360, angle)
    rotated_forward <-
      grepl("outward|inward|position", just) & (angle > 45 & angle < 135)
    rotated_backwards <-
      grepl("outward|inward|position", just) & (angle < -45 & angle > -135)

    ab <- ifelse(rotated_forward | rotated_backwards, b, a)
    swap_ab <- rotated_backwards | abs(angle) > 135

    if (any(just == "position")) {
      # compute justification based on saved original position
      # works only with special position functions from 'ggpp' and 'ggrepel'
      ab_orig <- paste(ab, "_orig", sep = "")
      position <- just == "position"
      if (!all(unique(ab_orig) %in% colnames(data))) {
        # this needs to be silent as it is normal behaviour
        just[position] <- "middle"
      } else {
        just[position] <-
          c("left", "middle", "right")[2L + 1L *
                                         sign(data[[ab_orig[1L]]][position] -
                                                data[[ab[1L]]][position])]
      }
    }

    if (any(grepl("outward|inward", just))) {
      # we allow tags for setting center/middle position for just
      just_used <- unique(just)
      just_special <- grep("_mean$|_median$|.*[0-9].*", just_used, value = TRUE)
      middle <- rep(0.5, length(just))
      for (j in just_special) { # skipped if just_special has length zero
        j_selector <- just == j
        if (j %in% c("outward_mean", "inward_mean")) {
          middle[j_selector & !swap_ab] <- mean(data[[a]])
          middle[j_selector & swap_ab] <- mean(data[[b]])
        } else if (j %in% c("outward_median", "inward_median")) {
          middle[j_selector & !swap_ab] <- stats::median(data[[a]])
          middle[j_selector & swap_ab] <- stats::median(data[[b]])
        } else {
          middle[j_selector & swap_ab] <- stats::median(data[[b]])
          middle_a <- as.numeric(gsub("outward_|inward_", "", unique(just)))
          if (a == "x") {
            tmp_data <- tibble::tibble(x = middle_a, y = data[[b]])
            middle[j_selector & !swap_ab] <- coord$transform(tmp_data, panel_params)$x
          } else {
            tmp_data <- tibble::tibble(y = middle_a, x = data[[b]])
            middle[j_selector & !swap_ab] <- coord$transform(tmp_data, panel_params)$y
          }
        }
      }

      just <- gsub("_.*$", "", just)

      # what follows is like in ggplot2 except for split_at
      obs <- data[[a]]
      obs[swap_ab] <- data[[b]][swap_ab]

      inward <- just == "inward"
      just[inward] <- c("left", "middle", "right")[just_dir(obs[inward],
                                                            split_at = middle[inward])]
      outward <- just == "outward"
      just[outward] <- c("right", "middle", "left")[just_dir(obs[outward],
                                                             split_at = middle[outward])]
    }
  }

  unname(c(left = 0, center = 0.5, right = 1,
           bottom = 0, middle = 0.5, top = 1)[just])
}

# modified from geom-text.r from 'ggplot2' 3.1.0 to support arbitrary split
# for "inward" and "outward" justification
just_dir <- function(x, tol = 0.001, split_at = 0.5) {
  out <- rep(2L, length(x))
  out[x < split_at - tol] <- 1L
  out[x > split_at + tol] <- 3L
  out
}

# as in pull request to ggplot2
# can be used by other geoms
compute_just <- function(just, a, b = a, angle = 0) {
  #  As justification direction is relative to the text, not the plotting area
  #  we need to swap x and y if text direction is rotated so that hjust is
  #  applied along y and vjust along x.
  if (any(grepl("outward|inward", just))) {
    # ensure all angles are in -360...+360
    angle <- angle %% 360
    # ensure correct behaviour for angles in -360...+360
    angle <- ifelse(angle > 180, angle - 360, angle)
    angle <- ifelse(angle < -180, angle + 360, angle)
    rotated_forward <-
      grepl("outward|inward", just) & (angle > 45 & angle < 135)
    rotated_backwards <-
      grepl("outward|inward", just) & (angle < -45 & angle > -135)

    ab <- ifelse(rotated_forward | rotated_backwards, b, a)
    just_swap <- rotated_backwards | abs(angle) > 135
    inward <-
      (just == "inward" & !just_swap | just == "outward" & just_swap)
    just[inward] <- c("left", "middle", "right")[just_dir(ab[inward])]
    outward <-
      (just == "outward" & !just_swap) | (just == "inward" & just_swap)
    just[outward] <- c("right", "middle", "left")[just_dir(ab[outward])]

  }

  unname(c(left = 0, center = 0.5, right = 1,
           bottom = 0, middle = 0.5, top = 1)[just])
}

# shorten segments to add padding
# code based on https://stackoverflow.com/questions/22649781/
#
# box.padding and point.pading givern in mm
#
shrink_segments <- function(data,
                            box.padding = 0,
                            point.padding = 0,
                            min.segment.length = 0.5) {
  stopifnot("'box.padding' must be >= 0" = box.padding >= 0,
            "'point.padding' must be >= 0" =  point.padding >= 0)
  segments.data <- data[ , c("x_orig", "y_orig", "x", "y")]
  starting.length <- apply(segments.data, 1,
                           function(x) stats::dist(rbind(x[1:2], x[3:4])))

  zero.len <- starting.length < max(1e-5, (point.padding + box.padding) / 100)

  # padding origin
  if (point.padding != 0) {
    p <- point.padding / starting.length[!zero.len] / 100
    p <- 1 - p
    segments.data[!zero.len, "x_orig"] <- data[!zero.len, "x"] + p * (data[!zero.len, "x_orig"] - data[!zero.len, "x"])
    segments.data[!zero.len, "y_orig"] <- data[!zero.len, "y"] + p * (data[!zero.len, "y_orig"] - data[!zero.len, "y"])
  }
  # padding position
  if (box.padding != 0) {
    p <- box.padding / starting.length[!zero.len] / 100
    p <- 1 - p
    segments.data[!zero.len, "x"] <- data[!zero.len, "x_orig"] + p * (data[!zero.len, "x"] - data[!zero.len, "x_orig"])
    segments.data[!zero.len, "y"] <- data[!zero.len, "y_orig"] + p * (data[!zero.len, "y"] - data[!zero.len, "y_orig"])
  }
  final.length <- apply(segments.data, 1,
                        function(x) stats::dist(rbind(x[1:2], x[3:4])))
  segments.data$too.short <- final.length < min.segment.length
  segments.data
}

# from ggplot2
resolve_text_unit <- function(unit) {
  unit <- rlang::arg_match0(unit, c("mm", "pt", "cm", "in", "pc"))
  switch(
    unit,
    "mm" = .pt,
    "cm" = .pt * 10,
    "in" = 72.27,
    "pc" = 12,
    1
  )
}

# from ggplot2 utilities-grid.R
#
# Name ggplot grid object
# Convenience function to name grid objects
#
# @keyword internal
ggname <- function(prefix, grob) {
  grob$name <- grobName(grob, prefix)
  grob
}
