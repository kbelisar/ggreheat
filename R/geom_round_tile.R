### The following is adapted code to create rounded geom_tile corners from both the ggplot2 geom_tile function source code and from
### boB Rudis' package statebins (adapted from geomRrect: https://github.com/hrbrmstr/statebins/blob/master/R/geom-rrect.r)

#' @name `%||%`
#' @noRd

`%||%` <- ggplot2:::`%||%`


#' @name geom_round_tile
#' @noRd

geom_round_tile <- function(mapping = NULL, data = NULL,
                            stat = "identity", position = "identity",
                            ...,
                            linewidth = NULL,
                            na.rm = FALSE,
                            show.legend = NA,
                            inherit.aes = TRUE) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomRoundTile,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      ...
    )
  )
}


GeomRoundTile <- ggplot2::ggproto("GeomRoundTile", ggplot2::GeomRect,

                                  setup_data = function(data, params) {

                                    data$width <- data$width %||% params$width %||% stats::ave(data$x, data$PANEL, FUN = function(x) ggplot2::resolution(x, FALSE, TRUE))
                                    data$height <- data$height %||% params$height %||% stats::ave(data$y, data$PANEL, FUN = function(y) ggplot2::resolution(y, FALSE, TRUE))

                                    transform(data,
                                              xmin = x - width / 2,  xmax = x + width / 2,  width = NULL,
                                              ymin = y - height / 2, ymax = y + height / 2, height = NULL
                                    )
                                  },

                                  default_aes = ggplot2::aes(colour = NA, fill = NA, alpha = NA, width = NA, height = NA),

                                  required_aes = c("x", "y"),

                                  draw_panel = function(self, data, panel_params, coord) {

                                    coords <- coord$transform(data, panel_params)

                                    lapply(1:length(coords$xmin), function(i) {

                                      grid::roundrectGrob(
                                        coords$xmin[i], coords$ymax[i],
                                        width = (coords$xmax[i] - coords$xmin[i]),
                                        height = (coords$ymax[i] - coords$ymin)[i],
                                        r = grid::unit(15,"pt"),
                                        default.units = "native",
                                        just = c("left", "top"),
                                        gp = grid::gpar(
                                          col = coords$colour[i],
                                          fill = ggplot2::alpha(coords$fill[i], coords$alpha[i]),
                                          lwd = coords$size[i],
                                          lty = coords$linetype[i]
                                        )
                                      )

                                    }) -> gl

                                    grobs <- do.call(grid::gList, gl)

                                    ggplot2:::ggname("geom_round_tile", grid::grobTree(children = grobs))

                                  },

                                  non_missing_aes = c("xmin", "xmax", "ymin", "ymax"),

                                  draw_key = ggplot2::draw_key_polygon
)

