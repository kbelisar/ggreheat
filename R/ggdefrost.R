#' This function is used to produce a single heatmap which combines two different correlation coefficients (Pearson's R and Spearman's Rho). Variables will be ordered from left to right and bottom to top as identified in 'vars1' and' vars2'.

#' @name ggdefrost
#' @param data A data frame consisting of the variables to be included in the heatmap.
#' @param vars1 The first of two possible list of variables to be included in the heatmap.
#' @param vars2 The second of two possible list of variables to be included in the heatmap.
#' @param type1 The type of correlation coefficient to calculate for the variables listed in the 'vars1' argument. One of c("pearson","spearman"), with the default being "pearson".
#' @param type2 The type of correlation coefficient to calculate for the variables listed in the 'vars2' argument. One of c("pearson","spearman"), with the default being "spearman".
#' @param diff The characteristic used to differentiate the variables listed in the 'vars2' argument. One of c("outline","round"), with the default being "outline".
#' @param alpha The alpha value to determine significance, in which significant p-values will have the corresponding correlation coefficient in the heatmap denoted in bold text. The default alpha value is 0.05.
#' @param coef_size The size of the correlation coefficient in the heatmap. The default is 7.
#' @param colour The colour of the outline for variables identified in 'vars2' when 'diff' is set to "outline". The default is a dark grey.
#' @param triangle Whether the heatmap takes on a square or triangle form. The default is FALSE.
#' @param vars_rename Optional names to reassign to the variables (in order of the variables identified in the vars1 and vars2 arguments).
#' @examples
#'
#' ### --- Load Simulated Data
#'
#' data("reheat_data")
#'
#' ### --- Heatmap with Pearson and Spearman correlations
#'
#' ### The default 'type1' argument is Pearson correlations for 'vars1' variables
#' ### The default 'type2' argument is Spearman correlations for 'vars2' variables
#' ### Variables in 'vars2' argument will be denoted with a round tile
#'
#' ggdefrost(reheat_data,
#' vars1 = c("audit_score","sip_score","drinks_per_week","cudit_score","ftnd_score"),
#' vars2 = c("cannabis_assist","cigarette_assist","ecigarette_assist"),
#' diff = "round")
#'
#' @return A ggplot graphical object.
#' @export

globalVariables(c("Var1","Var2","coefficient","pval_sign","differential"))

ggdefrost <- function(data, vars1 = NULL, vars2 = NULL, type1 = "pearson", type2 = "spearman", diff = "outline", alpha = 0.05, coef_size = 7, colour = NULL, triangle = FALSE, vars_rename = NULL){

  data <- as.data.frame(data)
  if(!type1 %in% c("pearson","spearman") & !type2 %in% c("pearson","spearman")) stop(rlang::format_error_bullets(c( x = c("Currently only Pearson's r or Spearman's rho are supported. Please select one of c('pearson','spearman')"))), call. = FALSE)
  if(!diff %in% c("outline","round")) stop(rlang::format_error_bullets(c( x = c("Not a valid 'diff' argument. Please select one of c('outline','round')"))), call. = FALSE)

  vars <- c(vars1,vars2) ### combined variable name list

  if(is.null(vars_rename)){
    var_order <- vars
  }

  if(!is.null(vars_rename)){

    data <- data[c(vars)]
    names(data) <- c(vars_rename)

    var_order <- vars_rename
    vars1 <- var_order[1:length(vars1)]
    vars2 <- var_order[length(vars1):length(vars_rename)]
  }


  colour2 <- ifelse(is.null(colour),"#333333",colour) ## default grey colour for outline of second variables (vars2)

  corr_dat1 <- {}
  corr_dat2 <- {}
  corr1 <- {}
  pval1 <- {}
  corr2 <- {}
  pval2 <- {}

  ### --- PEARSON CORRELATIONS
  suppressWarnings({
  if(!is.null(vars1)){

    for(var in var_order){
      corr_dat1[[var]] <- apply(data[c(var_order)], 2, stats::cor.test, data[[var]], method = type1, use = "pairwise.complete.obs")
    }
  for(i in 1:length(corr_dat1)){
    corr1[[i]] <- sapply(corr_dat1[[i]][1:length(corr_dat1)], "[[","estimate")
    pval1[[i]] <- sapply(corr_dat1[[i]][1:length(corr_dat1)], "[[", "p.value")
  }
    corr1 <- data.frame("corr" = do.call(rbind.data.frame,corr1))
    colnames(corr1) <- var_order
    rownames(corr1) <- var_order

    pval1 <- data.frame("pval" = do.call(rbind.data.frame,pval1))
    colnames(pval1) <- var_order
    rownames(pval1) <- var_order
  }
  })

  ### --- SPEARMAN CORRELATIONS
  suppressWarnings({
  if(!is.null(vars2)){
    for(var in var_order){
      corr_dat2[[var]] <- apply(data[c(var_order)], 2, stats::cor.test, data[[var]], method = type2, use = "pairwise.complete.obs")
    }
  for(i in 1:length(corr_dat2)){
    corr2[[i]] <- sapply(corr_dat2[[i]][1:length(corr_dat2)], "[[","estimate")
    pval2[[i]] <- sapply(corr_dat2[[i]][1:length(corr_dat2)], "[[", "p.value")
  }
    corr2 <- data.frame("corr" = do.call(rbind.data.frame,corr2))
    colnames(corr2) <- var_order
    rownames(corr2) <- var_order

    pval2 <- data.frame("pval" = do.call(rbind.data.frame,pval2))
    colnames(pval2) <- var_order
    rownames(pval2) <- var_order
  }
  })

  corr_final <- rbind(corr1[(rownames(corr1) %in% vars1),],corr2[(rownames(corr2) %in% vars2),])
  pval_final <- rbind(pval1[(rownames(pval1) %in% vars1),],pval2[(rownames(pval2) %in% vars2),])

  corr_final <- corr_final[c(var_order)]
  corr_final <- corr_final[match(var_order,rownames(corr_final)),]
  pval_final <- pval_final[c(var_order)]
  pval_final <- pval_final[match(var_order,rownames(pval_final)),]

  if(triangle==TRUE){
    corr_final <- upper_tri(corr_final)
    pval_final <- upper_tri(pval_final)
  }
  corr_dat <- base_melt(corr_final)
  pval_dat <- base_melt(pval_final)

  pval_dat$pval_sign <- ifelse(pval_dat$value<alpha,paste0("p < ",alpha),"N.S.")
  pval_dat$pval_sign[pval_dat$Var1==pval_dat$Var2] <- NA

  corr_dat$coefficient <- corr_dat$value

    ### MERGE p-values with the correlation coefficients

    corr_dat <- merge(corr_dat[c("Var1","Var2","coefficient")],
                      pval_dat[c("Var1","Var2","pval_sign")], by = c("Var1","Var2"), all = TRUE)

    corr_dat$differential <- ifelse(corr_dat$Var1 %in% vars2 | corr_dat$Var2 %in% vars2,TRUE,FALSE)

    if(diff=="outline"){
      hm_plot <- ggplot2::ggplot(data = corr_dat, ggplot2::aes(Var2, Var1)) +
        ggplot2::geom_tile(corr_dat, mapping = ggplot2::aes(Var2, Var1,fill = coefficient), width = 0.9, height = 0.9) +
        ggplot2::geom_tile(corr_dat, mapping = ggplot2::aes(Var2, Var1, colour = differential), alpha = 0, width = 0.88, height = 0.88, linewidth = 2, show.legend = FALSE) +
        ggplot2::scale_fill_gradient2(low = "#5f6488", high = "#5f8883", mid = "#c8c8c9",
                                      midpoint = 0, limit = c(-1,1), space = "Lab",
                                      name = "Correlation") +
        ggplot2::scale_colour_manual("",values = c(ggplot2::alpha("#c8c8c9",0),colour2)) +
        ggplot2::scale_y_discrete(position = "right") + ggplot2::theme_minimal() + ggplot2::coord_fixed(clip = "off") +
        ggplot2::geom_text(corr_dat[(corr_dat$pval_sign=="N.S." & !is.na(corr_dat$pval_sign)),],
                           mapping = ggplot2::aes(Var2, Var1, label = format(round(coefficient,2),nsmall = 2), colour = pval_sign),
                           size = coef_size, colour = "#333333", show.legend = FALSE) +
        ggplot2::geom_text(corr_dat[(corr_dat$pval_sign!="N.S." & !is.na(corr_dat$pval_sign)),],
                           mapping = ggplot2::aes(Var2, Var1, label = format(round(coefficient,2),nsmall = 2), colour = pval_sign),
                           size = coef_size, fontface = "bold", colour = "#000000", show.legend = FALSE) +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       axis.title = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.border = ggplot2::element_blank(),
                       legend.text = ggplot2::element_text(size = 15, face = "bold", hjust = "0.5"),
                       plot.title = ggplot2::element_text(size = 30, face = "bold", hjust = 0.5),
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

    if(diff=="round"){
      hm_plot <- ggplot2::ggplot(data = corr_dat, ggplot2::aes(Var2, Var1, fill = coefficient, colour = coefficient)) +
        ggplot2::geom_tile(corr_dat[(corr_dat$differential==FALSE),], mapping = ggplot2::aes(Var2, Var1), width = 0.9, height = 0.9) +
        geom_round_tile(corr_dat[(corr_dat$differential==TRUE),], mapping = ggplot2::aes(Var2, Var1), width = 0.9, height = 0.9, show.legend = FALSE) +
        ggplot2::scale_fill_gradient2(low = "#5f6488", high = "#5f8883", mid = "#c8c8c9",
                                      midpoint = 0, limit = c(-1,1), space = "Lab",
                                      name = "Correlation") +
        ggplot2::scale_colour_gradient2(low = "#5f6488", high = "#5f8883", mid = "#c8c8c9",
                                        midpoint = 0, limit = c(-1,1), space = "Lab",
                                        name = "Correlation") +
        ggplot2::scale_y_discrete(position = "right") + ggplot2::theme_minimal() + ggplot2::coord_fixed(clip = "off") +
        ggplot2::geom_text(corr_dat[(corr_dat$pval_sign=="N.S." & !is.na(corr_dat$pval_sign)),],
                           mapping = ggplot2::aes(Var2, Var1, label = format(round(coefficient,2),nsmall = 2), colour = pval_sign),
                           size = coef_size, colour = "#333333", show.legend = FALSE) +
        ggplot2::geom_text(corr_dat[(corr_dat$pval_sign!="N.S." & !is.na(corr_dat$pval_sign)),],
                           mapping = ggplot2::aes(Var2, Var1, label = format(round(coefficient,2),nsmall = 2), colour = pval_sign),
                           size = coef_size, fontface = "bold", colour = "#000000", show.legend = FALSE) +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       axis.title = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.border = ggplot2::element_blank(),
                       legend.text = ggplot2::element_text(size = 15, face = "bold", hjust = "0.5"),
                       plot.title = ggplot2::element_text(size = 30, face = "bold", hjust = 0.5),
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

    return(hm_plot)
}


