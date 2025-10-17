#' Position scales for continuous data (npcx & npcy)
#'
#' \code{scale_npcx_continuous()} and \code{scale_npcy_continuous()} are
#' scales for continuous npcx and npcy aesthetics expressed in "npc" units.
#' There are no variants. Obviously limits are always the full range of "npc"
#' units and transformations meaningless. These scales are used by the newly
#' defined aesthetics \code{npcx} and \code{npcy}.
#'
#' @param ... Other arguments passed on to \code{continuous_scale()}
#' @name scale_continuous_npc
#'
#' @return A \code{"Scale"} object.
#'
#' @export
#'
scale_npcx_continuous <- function(...) {
  ggplot2::continuous_scale(
    aesthetics = "npcx",
    palette = identity, name = NULL, breaks = NULL,
    minor_breaks = NULL, labels = NULL, limits = c(NA_real_, NA_real_),
    expand = c(0, 0, 0, 0), oob = scales::censor, na.value = NA_real_, transform = "identity",
    guide = "none", position = "bottom", super = ggplot2::ScaleContinuousPosition,
    ...
  )
}

#' @rdname scale_continuous_npc
#' @export
scale_npcy_continuous <- function(...) {
  ggplot2::continuous_scale(
    aesthetics = "npcy",
    palette = identity, name = NULL, breaks = NULL,
    minor_breaks = NULL, labels = NULL, limits = c(NA_real_, NA_real_),
    expand = c(0, 0, 0, 0), oob = scales::censor, na.value = NA_real_, transform = "identity",
    guide = "none", position = "bottom", super = ggplot2::ScaleContinuousPosition,
    ...
  )
}
