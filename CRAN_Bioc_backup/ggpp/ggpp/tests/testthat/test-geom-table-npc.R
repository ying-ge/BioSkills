context("geom_tb_npc")

library(tibble)
library(dplyr)

test_that("character label gives error in geom_table_npc", {
  my.df <- data.frame(x = (0:10) / 10, y = (0:10) / 10, tb = letters[1:11])
  expect_error(print(ggplot(my.df, aes(npcx = x, npcy = y, label = tb)) +
                       geom_table_npc()))
})

test_that("numeric label gives error in geom_table_npc", {
  my.df <- data.frame(x = (0:10) / 10, y = (0:10) / 10, tb = 1:11)
  expect_error(print(ggplot(my.df, aes(npcx = x, npcy = y, label = tb)) +
                       geom_table_npc()))
})

get_tb <- function() {
  mtcars %>%
    group_by(cyl) %>%
    summarize(wt = mean(wt), mpg = mean(mpg)) %>%
    ungroup() %>%
    mutate(
      wt = sprintf("%.2f", wt),
      mpg = sprintf("%.1f", mpg)
    )
}

test_that("geom_table_npc works as expected", {
  p <- ggplot(mtcars, aes(wt, mpg, colour = factor(cyl))) +
    geom_point()
  tb <- get_tb()

  tbbl <- tibble(x = 0.95, y = 0.95, tb = list(tb))
  result <- expect_silent(
    p + geom_table_npc(data = tbbl, aes(npcx = x, npcy = y, label = tb))
  )
  expect_s3_class(result, "ggplot")

  df <- data.frame(x = 0.95, y = 0.95, tb = I(list(tb)))
  result <- expect_silent(
    p + geom_table_npc(data = df, aes(npcx = x, npcy = y, label = tb))
  )
  expect_s3_class(result, "ggplot")
})


# test_that("data.frame", {
#   my.df <- data.frame(x = (0:10) / 10, y = (0:10) / 10, tb = letters[1:10])
#   expect_warning(ggplot(my.df, aes(npcx = x, npcy = y, label = tb)) +
#                    geom_table_npc())
#   vdiffr::expect_doppelganger("geom_table_npc_text_label",
#                               suppressWarnings(
#                                 ggplot(my.df, aes(npcx = x, npcy = y, label = tb)) +
#                                   geom_table_npc()
#                               )
#   )
# })

test_that("multiple_rows_tb_npc", {
  tb <- tibble(Z1 = LETTERS[2:4], z1 = letters[4:2])
  tbb <- tibble(Z2 = LETTERS[2:4], z2 = letters[4:2])
  tbbb <- tibble(Z3 = LETTERS[2:4], z3 = letters[4:2])
  my.tb <- tibble(x = c(0.1, 0.5, 0.9), y = c(0.1, 0.5, 0.9), tb = list(t1 = tb, t2 = tbb, t3 = tbbb))
  vdiffr::expect_doppelganger("geom_table_npc_multi_row",
                              ggplot() +
                                geom_table_npc(data = my.tb,
                                           mapping = aes(npcx = x, npcy = y, label = tb)) +
                                lims(x = c(0, 1), y = c(0, 1))
                              )
})

test_that("numbers_tb works with geom_table_npc", {
  my_data.tb <- tibble(x = -5:5, y = -5:5)
  tb <- tibble(a = 2:4, b = 4:2)
  my.tb <- tibble(x = 0.5, y = 0.5, tb = list(tb))
  vdiffr::expect_doppelganger("geom_table_npc_geom_data",
                              ggplot() +
                                geom_table_npc(data = my.tb,
                                           mapping = aes(npcx = x, npcy = y, label = tb)) +
                                lims(x = c(0, 1), y = c(0, 1))
                              )
  vdiffr::expect_doppelganger("geom_table_npc_plot_data",
                              ggplot(data = my.tb) +
                                geom_table_npc(mapping = aes(npcx = x, npcy = y, label = tb)) +
                                lims(x = c(0, 1), y = c(0, 1))
  )
  vdiffr::expect_doppelganger("geom_table_npc_vjust",
                              ggplot() +
                                geom_table_npc(data = my.tb,
                                               vjust = 1,
                                               mapping = aes(npcx = x, npcy = y, label = tb)) +
                                lims(x = c(0, 1), y = c(0, 1))
                              )
  vdiffr::expect_doppelganger("geom_table_npc_vjust_hjust",
                              ggplot() +
                                geom_table_npc(data = my.tb,
                                               vjust = 1, hjust = 1,
                                               mapping = aes(npcx = x, npcy = y, label = tb)) +
                                lims(x = c(0, 1), y = c(0, 1))
                              )
  vdiffr::expect_doppelganger("geom_table_npc_with_points",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           mapping = aes(npcx = x, npcy = y, label = tb))
                              )
})

test_that("alpha targets work with geom_table_npc", {
  my_data.tb <- tibble(x = -5:5, y = -5:5)
  tb <- tibble(a = 2:4, b = 4:2)
  my.tb <- tibble(x = 0.5, y = 0.5, tb = list(tb))
  vdiffr::expect_doppelganger("geom_table_npc_alpha_aes",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           alpha = 0.25,
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_alpha_all",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                            alpha = 0.25,
                                           alpha.target = "all",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_alpha_table",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           alpha = 0.25,
                                           alpha.target = "table",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_alpha_canvas",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           alpha = 0.25,
                                           alpha.target = "table.canvas",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_alpha_rules",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           alpha = 0.25,
                                           alpha.target = "table.rules",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_alpha_base",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           alpha = 0.25,
                                           alpha.target = "table.text",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_alpha_box",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           alpha = 0.25,
                                           alpha.target = "box",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_alpha_none",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           alpha = 0.25,
                                           alpha.target = "none",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

})

test_that("colour targets work in geom_table_npc", {
  my_data.tb <- tibble(x = -5:5, y = -5:5)
  tb <- tibble(a = 2:4, b = 4:2)
  my.tb <- tibble(x = 0.5, y = 0.5, tb = list(tb))

  vdiffr::expect_doppelganger("geom_table_npc_colour_na",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                            colour = NA,
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_colour_aes",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                            colour = "red",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_default_colour",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           default.colour = "blue",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_default_colour_none",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           default.colour = "blue",
                                           colour.target = "none",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_colour_all",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           colour = "red",
                                           colour.target = "all",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_colour_table",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           colour = "red",
                                           colour.target = "table",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_colour_rules",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           colour = "red",
                                           colour.target = "table.rules",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_colour_base",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           colour = "red",
                                           colour.target = "table.text",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_colour_box",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           colour = "red",
                                           colour.target = "box",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_colour_none",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           colour = "red",
                                           colour.target = "none",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

})

test_that("colour targets and alpha targets work together in geom_table_npc", {
  my_data.tb <- tibble(x = -5:5, y = -5:5)
  tb <- tibble(a = 2:4, b = 4:2)
  my.tb <- tibble(x = 0.5, y = 0.5, tb = list(tb))
  vdiffr::expect_doppelganger("geom_table_npc_colour_alpha_aes",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           colour = "red",
                                           alpha = 0.25,
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_colour_all_alpha",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           alpha = 0.25,
                                           alpha.target = "all",
                                           colour = "red",
                                           colour.target = "all",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_colour_table_alpha_table",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           alpha = 0.25,
                                           alpha.target = "table",
                                           colour = "red",
                                           colour.target = "table",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_colour_rules_alpha_table",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           alpha = 0.25,
                                           alpha.target = "table",
                                           colour = "red",
                                           colour.target = "table.rules",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_colour_rules_alpha_canvas",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           alpha = 0.25,
                                           alpha.target = "table.canvas",
                                           colour = "red",
                                           colour.target = "table.rules",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_default_colour_rules_alpha_canvas",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           alpha = 0.25,
                                           alpha.target = "table.canvas",
                                           default.colour = "blue",
                                           colour = "red",
                                           colour.target = "table.rules",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_colour_base_alpha_multiple",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           alpha = 0.25,
                                           alpha.target = c("table.rules",
                                                            "table.canvas"),
                                           colour = "red",
                                           colour.target = "table.text",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_colour_box_alpha_none",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                            alpha = 0.1,
                                           alpha.target = "none",
                                           colour = "red",
                                           colour.target = "box",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

})

test_that("fill and alpha work together in geom_table_npc", {
  my_data.tb <- tibble(x = -5:5, y = -5:5)
  tb <- tibble(a = 2:4, b = 4:2)
  my.tb <- tibble(x = 0.5, y = 0.5, tb = list(tb))
  vdiffr::expect_doppelganger("geom_table_npc_fill_aes",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           fill = "yellow",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_fill_alpha_all",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           alpha = 0.25,
                                           alpha.target = "all",
                                           fill = "yellow",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_fill_alpha_table",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           alpha = 0.25,
                                           alpha.target = "table",
                                           fill = "yellow",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_fill_alpha_segment",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           alpha = 0.25,
                                           alpha.target = "segment",
                                           fill = "yellow",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_fill_alpha_canvas",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           alpha = 0.25,
                                           alpha.target = "table.canvas",
                                           fill = "yellow",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_fill_alpha_multiple",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           alpha = 0.25,
                                           alpha.target = c("table.rules",
                                                            "table.canvas"),
                                           fill = "yellow",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_fill_alpha_none",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           alpha = 0.1,
                                           alpha.target = "none",
                                           fill = "yellow",
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

  vdiffr::expect_doppelganger("geom_table_npc_fill_alpha_na",
                              ggplot(my_data.tb, aes(x, y)) +
                                geom_point() +
                                geom_table_npc(data = my.tb,
                                           alpha = 0.1,
                                           alpha.target = "none",
                                           fill = NA,
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )

})

test_that("letters_tb works in geom_table_npc", {
  tb <- tibble(a = LETTERS[2:4], b = letters[4:2])
  my.tb <- tibble(x = 0.5, y = 0.5, tb = list(tb))
  vdiffr::expect_doppelganger("geom_table_npc_letters",
                              ggplot() +
                                geom_table_npc(data = my.tb,
                                           mapping = aes(npcx = x, npcy = y, label = tb))
                              )
})

test_that("parsed_tb works with geom_table_npc", {
  tb <- tibble("alpha" = c("x[2]~\"=\"~a^2", "sqrt(y)"),
               "beta" = c("x[2]~\"=\"~b^2", "sqrt(1/y)"))
  my.tb <- tibble(x = 0.5, y = 0.5, tb = list(tb))
  vdiffr::expect_doppelganger("geom_table_npc_parsed_all",
                              ggplot() +
                                geom_table_npc(my.tb, mapping = aes(npcx = x, npcy = y, label = tb), parse = TRUE)
                              )

  tb <- tibble("alpha" = c("x[2]~\"=\"~a^2", "text"),
               "beta" = c("x[2]~\"=\"~b^2", "1200"))
  my.tb <- tibble(x = 0.5, y = 0.5, tb = list(tb))
  vdiffr::expect_doppelganger("geom_table_npc_parsed_partial",
                              ggplot() +
                                geom_table_npc(my.tb, mapping = aes(npcx = x, npcy = y, label = tb), parse = TRUE)
                              )
  vdiffr::expect_doppelganger("geom_table_npc_parsed_partial_dark",
                              ggplot() +
                                geom_table_npc(my.tb,
                                           table.theme = ttheme_gtdark,
                                           mapping = aes(npcx = x, npcy = y, label = tb),
                                           parse = TRUE)
  )
})

test_that("geom_table_npc_string_left_hjust", {
  tb <- tibble(a = 2:4, b = 4:2)
  my.tb <- tibble(x = 0.5, y = 0.5, tb = list(tb))
x <- geom_table_npc(data = my.tb,
                mapping = aes(npcx = x, npcy = y, label = tb),
                table.hjust = "left")
expect_equal(x$geom_params$table.hjust, 0)
})

test_that("letters_tb works with geom_table_npc", {
  tb <- tibble(a = LETTERS[2:4], b = letters[4:2])
  my.tb <- tibble(x = 0.5, y = 0.5, tb = list(tb))
  vdiffr::expect_doppelganger("geom_table_npc_letters",
                              ggplot() +
                                geom_table_npc(data = my.tb,
                                           mapping = aes(npcx = x, npcy = y, label = tb))
  )
})

# this does not render the plot so it tests only the stat function, not the
# functions in StatTableNpc object.
test_that("string_center_hjust works with geom_table_npc", {
  tb <- tibble(a = 2:4, b = 4:2)
  my.tb <- tibble(x = 0, y = 0, tb = list(tb))
  p <- ggplot() +
    geom_table_npc(data = my.tb,
                   table.hjust = "center",
                   mapping = aes(npcx = x, npcy = y, label = tb))
  expect_equal(p$layers[[1]]$geom_params$table.hjust, 0.5)

  result <- layer_data(p)[, c("npcx", "npcy", "label", "hjust", "vjust")]
  expect_identical(result$label[[1]], tb)
})
