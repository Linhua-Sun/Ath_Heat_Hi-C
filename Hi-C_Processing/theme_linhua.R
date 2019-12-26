theme_linhua <-
    function (base_size = 20,
              base_family = "",
              border = T,
              margin = TRUE,
              legend = c("top", "bottom", "left", "right", "none"),
              x.text.angle = 0)
    {
        half_line <- base_size / 2
        if (!is.numeric(legend))
            legend <- match.arg(legend)
        if (x.text.angle > 5)
            xhjust <- 1
        else
            xhjust <- NULL
        if (border) {
            panel.border <- element_rect(fill = NA,
                                         colour = "black",
                                         size = 1)
            axis.line <- element_blank()
        }
        else {
            panel.border <- element_blank()
            axis.line = element_line(colour = "green", size = 0.5)
        }
        if (margin)
            plot.margin <- margin(half_line, half_line, half_line,
                                  half_line)
        else
            plot.margin <- unit(c(0.5, 0.3, 0.3, 0.3), "mm")
        .theme <-
            theme_bw(base_size = base_size, base_family = base_family) %+replace%
            theme(
                panel.border = panel.border,
                #panel.grid.major = element_line(colour = "grey",linetype = 1,size = 0.1),
                #panel.grid.minor = element_line(colour = "grey",linetype = 1,size = 0.01),
                panel.grid.major = element_line(colour = "grey92"),
                panel.grid.minor = element_line(colour = "grey92", size = 0.25),
                axis.line = axis.line,
                axis.text = element_text(
                    color = "black",
                    family = base_family,
                    size = base_size - 4
                ),
                axis.title = element_text(
                    family = base_family,
                    size = base_size,
                    colour = "black"
                ),
                axis.title.y = element_text(
                    angle = 90,
                    margin = margin(r = half_line * 1.5),
                    vjust = 1
                ),
                axis.title.x = element_text(margin = margin(t = half_line *
                                                                1.5), vjust = 1),
                legend.key = element_blank(),
                strip.background = element_rect(
                    fill = "#F2F2F2",
                    colour = "black",
                    size = 1
                ),
                strip.text = element_text(family = base_family, size = base_size-2),
                plot.margin = plot.margin,
                legend.position = legend,
                complete = TRUE,
                aspect.ratio = 1
            )
        if (x.text.angle != 0)
            .theme <-
            .theme + theme(axis.text.x = element_text(angle = x.text.angle,
                                                      hjust = xhjust))
        .theme
    }
