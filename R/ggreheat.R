#' This function is used to produce a heatmap using a single correlation coefficient (Pearson's r or Spearman's rho). It can also create side-by-side heatmaps for
#' a single binary grouping variable. When calculating Pearson's r, confidence intervals produced from Fisher Z transformations can also be displayed alongside the
#' correlation coefficient. If creating side-by-side heatmaps, the confidence intervals which do not overlap between the two groups will be highlighted.

#' @name ggreheat
#' @param data A data frame consisting only of the variables to be included in the heatmap.
#' @param vars The list of variables to be included in the heatmap.
#' @param type The type of correlation coefficient to calculate. The default is "pearson". One of c("pearson","spearman").
#' @param by A binary variable to create side-by-side heatmaps. If defined, side-by-side heatmaps will be created by the binary variable.
#' @param alpha The alpha value to determine significance, in which significant p-values will have the corresponding correlation coefficient in the heatmap denoted in bold text. The default alpha value is 0.05.
#' @param include_ci Whether confidence intervals calculated using Fisher Z transformations should be displayed in the heatmap. Only available when 'type' is "pearson". The default is TRUE.
#' @param coef_size The size of the correlation coefficient in the heatmap. The default is 7.
#' @param colour_ci The colour of the outline for highlighting Fisher Z transformations CIs that do not overlap. Available when type is "pearson" and a 'by' variable is defined. The default is a dark grey.
#' @param tile_type The type of tile corner to use in the heatmap. One of c("square","round"), with the default being "square".
#' @param triangle Whether the heatmap takes on a square or triangle form. The default is FALSE.
#' @param vars_rename Optional names to reassign to the variables listed in the vars argument.
#' @export

globalVariables(c("Var1","Var2","coefficient","pval_sign","CI","no_overlap"))

ggreheat <- function(data, vars = NULL, type = "pearson", by = NULL, alpha = 0.05, include_ci = FALSE, coef_size = 7, colour_ci = NULL, tile_type = "square", triangle = FALSE, vars_rename = NULL){

  data <- as.data.frame(data)

  if(type=="spearman" & include_ci==TRUE) stop(rlang::format_error_bullets(c( x = c("CIs are only available when 'type' is 'pearson'."))), call. = FALSE)
  if(!type %in% c("pearson","spearman")) stop(rlang::format_error_bullets(c( x = c("Currently only Pearson's r or Spearman's rho are supported. Please select one of c('pearson','spearman')"))), call. = FALSE)
  if(!tile_type %in% c("square","round")) stop(rlang::format_error_bullets(c( x = c("Not a valid 'tile_type' argument. Please select one of c('square','round')"))), call. = FALSE)

  if(!is.null(by)){
    data$by_var <- data[,by]
  }

  if(!is.null(by) & length(unique(data$by_var))!=2) stop(rlang::format_error_bullets(c( x = c("The 'by' variable is not a binary variable."))), call. = FALSE)

  if(!is.null(vars_rename)){
    if(is.null(by)){
      data <- data[c(vars)]
      names(data) <- c(vars_rename)
    }
    if(!is.null(by)){
      data <- data[c("by_var",vars)]
      names(data) <- c("by_var",vars_rename)
    }
    vars <- vars_rename
  }

  corr_type <- ifelse(type=="pearson","Pearson Correlation","Spearman Correlation")
  colour_ci <- ifelse(is.null(colour_ci),"#333333",colour_ci) ## default grey colour for non-overlapping CI

  ### ---CORRELATIONS (NO BY VARIABLE)

  if(is.null(by)){

    corr_dat <- {}
    corr <- {}
    pval <- {}
    ci_corr <- {}

    suppressWarnings({
      for(var in vars){
        corr_dat[[var]] <- apply(data[c(vars)], 2, stats::cor.test, data[[var]], method = type, use = "pairwise.complete.obs")
      }
    })

    for(i in 1:length(corr_dat)){
      corr[[i]] <- sapply(corr_dat[[i]][1:length(corr_dat)], "[[","estimate")
      pval[[i]] <- sapply(corr_dat[[i]][1:length(corr_dat)], "[[", "p.value")
    }

    corr_final <- data.frame("corr" = do.call(rbind.data.frame,corr))
    colnames(corr_final) <- vars
    rownames(corr_final) <- vars

    pval_final <- data.frame("pval" = do.call(rbind.data.frame,pval))
    colnames(pval_final) <- vars
    rownames(pval_final) <- vars

    if(include_ci==TRUE){
      for(i in 1:length(corr_dat)){
        ci_corr[[i]] <- sapply(corr_dat[[i]][1:length(corr_dat)], "[[","conf.int", simplify = FALSE)
      }

      ci_final <- data.frame("X" = 1:length(ci_corr))
      for(i in 1:length(ci_corr)){
        ci_dat <- as.data.frame(t(as.data.frame(ci_corr[[i]])))
        ci_dat_frame <- data.frame("CI" = paste0(format(round(ci_dat$V1,2),nsmall = 2)," \u2014 ",format(round(ci_dat$V2,2),nsmall = 2)))
        ci_final <- cbind(ci_final,as.data.frame(ci_dat_frame))
      }

      ci_final <- ci_final[c(-1)] ## remove X placeholder columns
      colnames(ci_final) <- vars
      rownames(ci_final) <- vars
    }

    if(triangle==TRUE){
      corr_final <- upper_tri(corr_final)
      pval_final <- upper_tri(pval_final)
    }
    ### LONG format of matrix
    corr_dat <- base_melt(corr_final)
    pval_dat <- base_melt(pval_final)

    pval_dat$pval_sign <- ifelse(pval_dat$value<alpha,paste0("p < ",alpha),"N.S.")
    pval_dat$pval_sign[pval_dat$Var1==pval_dat$Var2] <- NA

    corr_dat$coefficient <- corr_dat$value

    ### MERGE p-values with the correlation coefficients

    corr_dat <- merge(corr_dat[c("Var1","Var2","coefficient")],
                      pval_dat[c("Var1","Var2","pval_sign")], by = c("Var1","Var2"), all = TRUE)

    if(include_ci==TRUE){

      if(triangle==TRUE){
        ci_final <- upper_tri(ci_final)
      }
      ci_dat <- base_melt(ci_final)

      ci_dat$CI <- ci_dat$value
      ci_dat$CI[ci_dat$Var1==ci_dat$Var2] <- NA
      corr_dat <- merge(corr_dat,ci_dat[c("Var1","Var2","CI")], by = c("Var1","Var2"), all = TRUE)
    }
  }

  ### ---CORRELATIONS (BY GROUPING VARIABLE)

  if(!is.null(by)){

    group_corr <- data.frame("Var1" = NULL, "Var2" = NULL, "by_var" = NULL, "coefficient" = NULL, "pval_sign" = NULL)

    for(group_u in unique(data$by_var)){

      g_data <- data[(data$by_var==group_u),]
      corr_dat_g <- {}
      corr_g <- {}
      pval_g <- {}
      ci_corr_g <- {}

      suppressWarnings({
        for(var in vars){
          corr_dat_g[[var]] <- apply(g_data[c(vars)], 2, stats::cor.test, g_data[[var]], method = type, use = "pairwise.complete.obs")
        }
      })

      for(i in 1:length(corr_dat_g)){
        corr_g[[i]] <- sapply(corr_dat_g[[i]][1:length(corr_dat_g)], "[[","estimate")
        pval_g[[i]] <- sapply(corr_dat_g[[i]][1:length(corr_dat_g)], "[[", "p.value")
      }

      corr_final <- data.frame("corr" = do.call(rbind.data.frame,corr_g))
      colnames(corr_final) <- vars
      rownames(corr_final) <- vars

      pval_final <- data.frame("pval" = do.call(rbind.data.frame,pval_g))
      colnames(pval_final) <- vars
      rownames(pval_final) <- vars

      if(include_ci==TRUE){
        for(i in 1:length(corr_dat_g)){
          ci_corr_g[[i]] <- sapply(corr_dat_g[[i]][1:length(corr_dat_g)], "[[","conf.int", simplify = FALSE)
        }

        ci_final <- data.frame("X" = 1:length(ci_corr_g))
        for(i in 1:length(ci_corr_g)){
          ci_dat <- as.data.frame(t(as.data.frame(ci_corr_g[[i]])))
          ci_dat$lower <- apply(ci_dat[c("V1","V2")], 1, function(x) min(x))
          ci_dat$upper <- apply(ci_dat[c("V1","V2")], 1, function(x) max(x))
          ci_dat_frame <- data.frame("CI" = paste0(format(round(ci_dat$lower,2),nsmall = 2)," \u2014 ",format(round(ci_dat$upper,2),nsmall = 2)))
          ci_final <- cbind(ci_final,as.data.frame(ci_dat_frame))
        }

        ci_final <- ci_final[c(-1)] ## remove X placeholder columns
        colnames(ci_final) <- vars
        rownames(ci_final) <- vars
      }

      ### LONG format of matrix
      corr_dat <- base_melt(corr_final)
      pval_dat <- base_melt(pval_final)

      pval_dat$pval_sign <- ifelse(pval_dat$value<alpha,paste0("p < ",alpha),"N.S.")
      pval_dat$pval_sign[pval_dat$Var1==pval_dat$Var2] <- NA

      corr_dat$coefficient <- corr_dat$value

      corr_dat$by_var <- group_u

      ### MERGE p-values with the correlation coefficients

      corr_dat <- merge(corr_dat[c("by_var","Var1","Var2","coefficient")],
                        pval_dat[c("Var1","Var2","pval_sign")], by = c("Var1","Var2"), all = TRUE)

      if(include_ci==TRUE){
        ci_dat <- base_melt(ci_final)

        ci_dat$CI <- ci_dat$value
        ci_dat$CI[ci_dat$Var1==ci_dat$Var2] <- NA
        corr_dat <- merge(corr_dat,ci_dat[c("Var1","Var2","CI")], by = c("Var1","Var2"), all = TRUE)
      }

      group_corr <- rbind(group_corr,corr_dat)
    }

    if(include_ci==TRUE){

      group_corr$lower_ci <- as.numeric(gsub(" ","",(sub("(.*)\u2014.*", "\\1", group_corr$CI))))
      group_corr$upper_ci <- as.numeric(gsub(" ","",(sub(".*\u2014(.*)", "\\1", group_corr$CI))))

      group_corr$by_num <- as.numeric(as.factor(group_corr$by_var))

      corr_wide <- stats::reshape(as.data.frame(group_corr[c("Var1","Var2","lower_ci","upper_ci","by_num")]), idvar = c("Var1","Var2"), timevar = "by_num",
                           v.names = c("lower_ci","upper_ci"), sep = "_", direction = "wide")

      corr_wide_ci <- data.frame("Var1" = NULL, "VAR2" = NULL, "no_overlap" = NULL)

      for(r in 1:NROW(corr_wide)){
        corr_dat <- corr_wide[(r),]
        if(corr_dat$Var1==corr_dat$Var2){
          corr_dat$no_overlap <- FALSE
        }
        if(corr_dat$Var1!=corr_dat$Var2){
          corr_dat$no_overlap <- ifelse(corr_dat$lower_ci_1 <= corr_dat$upper_ci_2 && corr_dat$lower_ci_2 <= corr_dat$upper_ci_1,FALSE,TRUE)
        }
        corr_wide_ci <- rbind(corr_wide_ci,corr_dat[c("Var1","Var2","no_overlap")])
      }

      group_corr <- merge(group_corr,corr_wide_ci, by = c("Var1","Var2"), all.x = TRUE)

    }

    }

  if(is.null(by)){ ## SINGLE HEATMAPS

    if(tile_type=="square"){
      ## linewidth units == 0.01 width/height unit
      hm_plot <- ggplot2::ggplot(data = corr_dat, ggplot2::aes(Var2, Var1)) +
        ggplot2::geom_tile(corr_dat, mapping = ggplot2::aes(Var2, Var1, fill = coefficient), width = 0.9, height = 0.9) +
        ggplot2::scale_fill_gradient2(low = "#5f6488", high = "#5f8883", mid = "#c8c8c9",
                                      midpoint = 0, limit = c(-1,1), space = "Lab",
                                      name = corr_type) +
        ggplot2::scale_y_discrete(position = "right") + ggplot2::theme_minimal() + ggplot2::coord_fixed(clip = "off") +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       axis.title = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.border = ggplot2::element_blank(),
                       legend.text = ggplot2::element_text(size = 15, face = "bold", hjust = "0.5"),
                       legend.title = ggplot2::element_text(face = "bold",size = 15, hjust = 0.5),
                       legend.direction = "horizontal",
                       legend.position = "bottom",
                       axis.text.y = ggplot2::element_text(face = "bold", size = 20, colour = "black"),
                       axis.text.x = ggplot2::element_text(face = "bold", size = 20, angle = 30, vjust = 1, hjust = 1, colour = "black"),
                       axis.line = ggplot2::element_line(linewidth = 2, lineend = "square"),
                       axis.ticks.length = ggplot2::unit(3,"mm"),
                       axis.ticks = ggplot2::element_line(linewidth = 2)) +
        ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = 30, barheight = 1, title.position = "top", ticks.linewidth = 0))
    }

    if(tile_type=="round"){
      ## linewidth units == 0.01 width/height unit
      hm_plot <- ggplot2::ggplot(data = corr_dat, ggplot2::aes(Var2, Var1)) +
        geom_round_tile(corr_dat, mapping = ggplot2::aes(Var2, Var1, fill = coefficient), width = 0.9, height = 0.9) +
        ggplot2::scale_fill_gradient2(low = "#5f6488", high = "#5f8883", mid = "#c8c8c9",
                                      midpoint = 0, limit = c(-1,1), space = "Lab",
                                      name = corr_type) +
        ggplot2::scale_y_discrete(position = "right") + ggplot2::theme_minimal() + ggplot2::coord_fixed(clip = "off") +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       axis.title = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.border = ggplot2::element_blank(),
                       legend.text = ggplot2::element_text(size = 15, face = "bold", hjust = "0.5"),
                       legend.title = ggplot2::element_text(face = "bold",size = 15, hjust = 0.5),
                       legend.direction = "horizontal",
                       legend.position = "bottom",
                       axis.text.y = ggplot2::element_text(face = "bold", size = 20, colour = "black"),
                       axis.text.x = ggplot2::element_text(face = "bold", size = 20, angle = 30, vjust = 1, hjust = 1, colour = "black"),
                       axis.line = ggplot2::element_line(linewidth = 2, lineend = "square"),
                       axis.ticks.length = ggplot2::unit(3,"mm"),
                       axis.ticks = ggplot2::element_line(linewidth = 2)) +
        ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = 30, barheight = 1, title.position = "top", ticks.linewidth = 0))
    }

  if(include_ci==FALSE){
    hm_plot <- hm_plot +
      ggplot2::geom_text(corr_dat[(corr_dat$pval_sign=="N.S." & !is.na(corr_dat$pval_sign)),],
                         mapping = ggplot2::aes(Var2, Var1, label = format(round(coefficient,2),nsmall = 2), colour = pval_sign),
                         size = coef_size, colour = "#333333", show.legend = FALSE) +
      ggplot2::geom_text(corr_dat[(corr_dat$pval_sign!="N.S." & !is.na(corr_dat$pval_sign)),],
                         mapping = ggplot2::aes(Var2, Var1, label = format(round(coefficient,2),nsmall = 2), colour = pval_sign),
                         size = coef_size, fontface = "bold", colour = "#000000", show.legend = FALSE)
  }

  if(include_ci==TRUE){
    hm_plot <- hm_plot +
      ggplot2::geom_text(corr_dat[(corr_dat$pval_sign=="N.S." & !is.na(corr_dat$pval_sign)),],
                         mapping = ggplot2::aes(Var2, Var1, label = format(round(coefficient,2),nsmall = 2), colour = pval_sign),
                         size = coef_size, colour = "#333333", show.legend = FALSE, vjust = 0) +
      ggplot2::geom_text(corr_dat[(corr_dat$pval_sign!="N.S." & !is.na(corr_dat$pval_sign)),],
                         mapping = ggplot2::aes(Var2, Var1, label = format(round(coefficient,2),nsmall = 2), colour = pval_sign),
                         size = coef_size, fontface = "bold", colour = "#000000", show.legend = FALSE, vjust = 0) +
      ggplot2::geom_text(corr_dat,
                         mapping = ggplot2::aes(Var2, Var1, label = CI),
                         size = coef_size*0.7, fontface = "italic", colour = "#000000", show.legend = FALSE, vjust = 2)
  }
  }

  if(!is.null(by)){ ## SIDE-BY-SIDE HEATMAPS

    if(tile_type=="square"){
      ## linewidth units == 0.01 width/height unit
      hm_plot <- ggplot2::ggplot(data = group_corr, ggplot2::aes(Var2, Var1)) + ggplot2::facet_wrap(~by_var) +
        ggplot2::geom_tile(group_corr, mapping = ggplot2::aes(Var2, Var1,fill = coefficient), width = 0.9, height = 0.9) +
        ggplot2::scale_fill_gradient2(low = "#5f6488", high = "#5f8883", mid = "#c8c8c9",
                                      midpoint = 0, limit = c(-1,1), space = "Lab",
                                      name = corr_type) +
        ggplot2::scale_y_discrete(position = "right") + ggplot2::theme_minimal() + ggplot2::coord_fixed(clip = "off") +
        ggplot2::theme(strip.text = ggplot2::element_text(size = 25, face = "bold"),
                       panel.spacing = ggplot2::unit("20","pt"),
                       strip.background = ggplot2::element_blank(),
                       axis.title.x = ggplot2::element_blank(),
                       axis.title = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.border = ggplot2::element_blank(),
                       legend.text = ggplot2::element_text(size = 15, face = "bold", hjust = "0.5"),
                       legend.title = ggplot2::element_text(face = "bold",size = 15, hjust = 0.5),
                       legend.direction = "horizontal",
                       legend.position = "bottom",
                       axis.text.y = ggplot2::element_text(face = "bold", size = 20, colour = "black"),
                       axis.text.x = ggplot2::element_text(face = "bold", size = 20, angle = 30, vjust = 1, hjust = 1, colour = "black"),
                       axis.line = ggplot2::element_line(linewidth = 2, lineend = "square"),
                       axis.ticks.length = ggplot2::unit(3,"mm"),
                       axis.ticks = ggplot2::element_line(linewidth = 2)) +
        ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = 60, barheight = 1, title.position = "top", ticks.linewidth = 0))
    }

    if(tile_type=="round"){
      ## linewidth units == 0.01 width/height unit
      hm_plot <- ggplot2::ggplot(data = group_corr, ggplot2::aes(Var2, Var1)) + ggplot2::facet_wrap(~by_var) +
        geom_round_tile(group_corr, mapping = ggplot2::aes(Var2, Var1,fill = coefficient), width = 0.9, height = 0.9) +
        ggplot2::scale_fill_gradient2(low = "#5f6488", high = "#5f8883", mid = "#c8c8c9",
                                      midpoint = 0, limit = c(-1,1), space = "Lab",
                                      name = corr_type) +
        ggplot2::scale_y_discrete(position = "right") + ggplot2::theme_minimal() + ggplot2::coord_fixed(clip = "off") +
        ggplot2::theme(strip.text = ggplot2::element_text(size = 25, face = "bold"),
                       panel.spacing = ggplot2::unit("20","pt"),
                       strip.background = ggplot2::element_blank(),
                       axis.title.x = ggplot2::element_blank(),
                       axis.title = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.border = ggplot2::element_blank(),
                       legend.text = ggplot2::element_text(size = 15, face = "bold", hjust = "0.5"),
                       legend.title = ggplot2::element_text(face = "bold",size = 15, hjust = 0.5),
                       legend.direction = "horizontal",
                       legend.position = "bottom",
                       axis.text.y = ggplot2::element_text(face = "bold", size = 20, colour = "black"),
                       axis.text.x = ggplot2::element_text(face = "bold", size = 20, angle = 30, vjust = 1, hjust = 1, colour = "black"),
                       axis.line = ggplot2::element_line(linewidth = 2, lineend = "square"),
                       axis.ticks.length = ggplot2::unit(3,"mm"),
                       axis.ticks = ggplot2::element_line(linewidth = 2)) +
        ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = 60, barheight = 1, title.position = "top", ticks.linewidth = 0))

    }

    if(include_ci==FALSE){
      hm_plot <- hm_plot +
        ggplot2::geom_text(group_corr[(group_corr$pval_sign=="N.S." & !is.na(group_corr$pval_sign)),],
                           mapping = ggplot2::aes(Var2, Var1, label = format(round(coefficient,2),nsmall = 2), colour = pval_sign),
                           size = coef_size, colour = "#333333", show.legend = FALSE) +
        ggplot2::geom_text(group_corr[(group_corr$pval_sign!="N.S." & !is.na(group_corr$pval_sign)),],
                           mapping = ggplot2::aes(Var2, Var1, label = format(round(coefficient,2),nsmall = 2), colour = pval_sign),
                           size = coef_size, fontface = "bold", colour = "#000000", show.legend = FALSE)
    }

    if(include_ci==TRUE){

      if(tile_type=="square"){
        hm_plot <- hm_plot +
          ggplot2::geom_tile(group_corr, mapping = ggplot2::aes(Var2, Var1, colour = no_overlap), alpha = 0, width = 0.88, height = 0.88, linewidth = 2, show.legend = FALSE) +
          ggplot2::scale_colour_manual("",values = c(ggplot2::alpha("#333333",0),colour_ci)) +
          ggplot2::geom_text(group_corr[(group_corr$pval_sign=="N.S." & !is.na(group_corr$pval_sign)),],
                             mapping = ggplot2::aes(Var2, Var1, label = format(round(coefficient,2),nsmall = 2), colour = pval_sign),
                             size = coef_size, colour = "#333333", show.legend = FALSE, vjust = 0) +
          ggplot2::geom_text(group_corr[(group_corr$pval_sign!="N.S." & !is.na(group_corr$pval_sign)),],
                             mapping = ggplot2::aes(Var2, Var1, label = format(round(coefficient,2),nsmall = 2), colour = pval_sign),
                             size = coef_size, fontface = "bold", colour = "#000000", show.legend = FALSE, vjust = 0) +
          ggplot2::geom_text(group_corr,
                             mapping = ggplot2::aes(Var2, Var1, label = CI),
                             size = coef_size*0.7, fontface = "italic", colour = "#000000", show.legend = FALSE, vjust = 2)
      }
      if(tile_type=="round"){
        hm_plot <- hm_plot +
          geom_round_tile(group_corr, mapping = ggplot2::aes(Var2, Var1, colour = no_overlap), alpha = 0, width = 0.88, height = 0.88, linewidth = 2, show.legend = FALSE) +
          ggplot2::scale_colour_manual("",values = c(ggplot2::alpha("#333333",0),colour_ci)) +
          ggplot2::geom_text(group_corr[(group_corr$pval_sign=="N.S." & !is.na(group_corr$pval_sign)),],
                             mapping = ggplot2::aes(Var2, Var1, label = format(round(coefficient,2),nsmall = 2), colour = pval_sign),
                             size = coef_size, colour = "#333333", show.legend = FALSE, vjust = 0) +
          ggplot2::geom_text(group_corr[(group_corr$pval_sign!="N.S." & !is.na(group_corr$pval_sign)),],
                             mapping = ggplot2::aes(Var2, Var1, label = format(round(coefficient,2),nsmall = 2), colour = pval_sign),
                             size = coef_size, fontface = "bold", colour = "#000000", show.legend = FALSE, vjust = 0) +
          ggplot2::geom_text(group_corr,
                             mapping = ggplot2::aes(Var2, Var1, label = CI),
                             size = coef_size*0.7, fontface = "italic", colour = "#000000", show.legend = FALSE, vjust = 2)

      }
     }
  }

    return(hm_plot)

}


