### function for multinomial regression with binary exposure


analyse_binary_multinomial <- function(datasetwide, datasetlong, varname, traitname, missingnessplot = FALSE) {

  # datasetwide, datasetlong, varname, traitname

  ### initilise empty lists for saving and returning plots and results

  plots_ls <- list()
  multinomial_res_ls <- list()

  datasetlong[["predclass"]] <- as.factor(as.character(datasetlong[["predclass"]]))
  datasetwide[["predclass"]] <- as.factor(as.character(datasetwide[["predclass"]]))

  ### create formula for multinomial regression

  formula_i <- formula(paste0("predclass ~ ", varname))

  for (reference_i in c("main", "relevel")) {

    ### if loop is 'relevel' then change the reference class to "3"

    if (reference_i == "relevel") {
      datasetwide[["predclass"]] <-
        relevel(as.factor(datasetwide[["predclass"]]), ref = 3)

      datasetlong[["predclass"]] <-
        relevel(as.factor(datasetlong[["predclass"]]), ref = 3)
      palette_choice4 <- palette_choice_ma
      
    } else {
      print(paste0("This is the main analysis (reference level ", levels(datasetlong[["predclass"]])[1], ")"))
      palette_choice4 <- palette_choice_ir
      
    }

    ### multinomial logistic regression naive to probabilities

    mres <- multinom(formula = formula_i, data = datasetwide)

    tidy_mres <-
      tidy(mres,
        conf.int = TRUE,
        conf.level = 0.95,
        exponentiate = TRUE
      )

    tidy_mres$model <- "naive"
    tidy_mres$exposure <- varname

    print("Completed naive model estimation")
    print(tidy_mres)

    ### multinomial logistic regression weighted for probabilities
    mres_weighted <-
      multinom(formula_i, data = datasetlong, weights = probabilities)

    # extract and exponentiate coefficients
    tidy_mres_weighted <-
      tidy(mres_weighted,
        conf.int = TRUE,
        conf.level = 0.95,
        exponentiate = TRUE
      )

    tidy_mres_weighted$model <- "weighted"
    tidy_mres_weighted$exposure <- varname

    print("Completed weighted model estimation")
    print(tidy_mres_weighted)

    ### bind and drop intercept from results

    tidy_mres_df <- tidy_mres[which(tidy_mres$term != "(Intercept)"), ]

    tidy_mres_weighted_df <- tidy_mres_weighted[which(tidy_mres_weighted$term != "(Intercept)"), ]

    # join names
    tidy_mres_df <- left_join(tidy_mres_df, class_names, by=c("y.level" = "class"))
    tidy_mres_weighted_df <- left_join(tidy_mres_weighted_df, class_names, by=c("y.level" = "class"))
    
    refclass <- as.character(class_names$class_abbreviation)[which(class_names$class==levels(datasetwide[["predclass"]])[1])]
    refclasslong <- as.character(class_names$class_abbreviation)[which(class_names$class==levels(datasetwide[["predclass"]])[1])]
    
    ### create a 'comparison' column of statistical test for plot
    tidy_mres_df$comparison <- paste0(refclass, " v ", tidy_mres_df$class_abbreviation)
    tidy_mres_weighted_df$comparison <- paste0(refclasslong, " v ", tidy_mres_weighted_df$class_abbreviation)
    
    ### remove periods if any exist from names

    tidy_mres_df$exposure <- as.factor(gsub("\\.", " ", as.character(tidy_mres_df$exposure)))

    tidy_mres_weighted_df$exposure <- as.factor(gsub("\\.", " ", as.character(tidy_mres_weighted_df$exposure)))

    write.csv(tidy_mres_weighted_df, paste0(output_dir, "/", "multinomial_", traitname, "_weighted_res_", reference_i, ".csv"), col.names = T, row.names = F, quote = F)

    title_lb <- paste0(gsub("Self reported ", "", as.character(tidy_mres_weighted_df$exposure[1])))
    print(title_lb)

    ### plot multinomial results
    gg01 <- ggplot(tidy_mres_df, aes(comparison, estimate, colour = comparison)) +
      geom_point(alpha = 0.8) +
      geom_hline(yintercept = 1, color = "coral", alpha = 0.4) +
      geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
      xlab("Class comparisons") +
      scale_colour_manual(values=palette_choice4) +
      ylab("Estimate (RR)") +
      theme_classic() +
      labs(colour = "", title = title_lb) +
      theme(legend.position = "none")

    # theme(
    #  panel.background = element_rect(fill = "white", colour = "black"),
    #  strip.background = element_rect(fill = "white", colour = "white")
    # ) +

    png(
      filename = paste0(output_dir, "/", "multinomial_naive_", traitname, "_", reference_i, "_", optimum_classes, "class_scatterplot.png"),
      width = 15, height = 10, units = "cm", res = 300
    )

    print(gg01)

    dev.off()

    ###

    gg02 <- ggplot(tidy_mres_weighted_df, aes(comparison, estimate, colour = comparison)) +
      geom_point(alpha = 0.8) +
      geom_hline(yintercept = 1, color = "coral", alpha = 0.4) +
      geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
      xlab("Class comparisons") +
      scale_colour_manual(values=palette_choice4) +
      ylab("Estimate (RR)") +
      theme_classic() +
      labs(colour = "", title = title_lb) +
      theme(legend.position = "none")

    png(
      filename = paste0(output_dir, "/", "multinomial_weighted_", traitname, "_", reference_i, "_", optimum_classes, "class_scatterplot.png"),
      width = 15, height = 10, units = "cm", res = 300
    )

    print(gg02)

    dev.off()

    ###

    ###

    tidy_mres_weighted_df$estimatelog <- log(tidy_mres_weighted_df$estimate)
    tidy_mres_weighted_df$conf.lowlog <- log(tidy_mres_weighted_df$conf.low)
    tidy_mres_weighted_df$conf.highlog <- log(tidy_mres_weighted_df$conf.high)

    print(tidy_mres_weighted_df)
    
    gg02log <- ggplot(tidy_mres_weighted_df, aes(comparison, estimatelog, colour = comparison)) +
      geom_point(alpha = 0.8) +
      geom_hline(yintercept = 0, color = "coral", alpha = 0.4) +
      geom_pointrange(aes(ymin = conf.lowlog, ymax = conf.highlog)) +
      xlab("Class comparisons") +
      scale_colour_manual(values=palette_choice4) +
      ylab("Estimate log(RR)") +
      theme_classic() +
      ylim(-1, 6) +
      labs(colour = "", title = title_lb) +
      theme(legend.position = "none")

    png(
      filename = paste0(output_dir, "/", "multinomial_weighted_log_", traitname, "_", reference_i, "_", optimum_classes, "class_scatterplot.png"),
      width = 15, height = 10, units = "cm", res = 300
    )

    print(gg02log)

    dev.off()

    ### bind results and plots in to lists to return

    plots_ls[[paste0("naive_", reference_i)]] <- gg01
    plots_ls[[paste0("weighted_", reference_i)]] <- gg02
    plots_ls[[paste0("weightedlog_", reference_i)]] <- gg02log

    multinomial_res_ls[[paste0("naive_", reference_i)]] <- tidy_mres_df
    multinomial_res_ls[[paste0("weighted_", reference_i)]] <- tidy_mres_weighted_df

    ### tabulate abosulte number and percentages in each class also

    tabulations <- list()

    for (i in 1:optimum_classes) {
      tb1 <- tabyl(datasetwide[[varname]][datasetwide[["predclass"]] == i], show.na = FALSE)

      names(tb1)[1] <- paste0("class_", i)

      tb1$categories <- tb1[, 1]

      tabulations[[i]] <- tb1
    }

    tabulations2 <- bind_rows(tabulations, .id = c("class"))

    tabulations2 <- left_join(tabulations2, class_names[,c("class","class_abbreviation")])
    
    tabulations2$class <- tabulations2$class_abbreviation
    
    tabulations2$class_abbreviation <- NULL
    
    #tabulations2 <- tabulations2[c("categories", "n", "percent", "valid_percent", "class")]

    tabulations2 <- tabulations2 %>% select(any_of(c("categories", "n", "percent", "valid_percent", "class")))
    
    tabulations2$categories[is.na(tabulations2$categories)] <- "Missing - NA"
    
    # if valid_percent missing because there are no NAs in var then percent is the valid_percent
    if (is.null(tabulations2$valid_percent)) {
      tabulations2$valid_percent <- tabulations2$percent
    } else {}
    
    tabulations2$`N (percent)` <-
      paste0(tabulations2$n, " (", round(tabulations2$valid_percent * 100, digits = 2), ")")

    ### remove periods if any exist

    tabulations2[[traitname]] <- as.factor(gsub("\\.", " ", as.character(tabulations2$categories)))

    write.table(tabulations2[c("categories", "n", "percent", "class", "N (percent)")],
      file = paste0(output_dir, "/", traitname, "_", optimum_classes, "class_tabulations.csv"),
      row.names = FALSE, quote = FALSE, sep = ","
    )
    ### tabulate numbers and subset any rows of NA
    
    missing_duration <- NROW(datasetwide[[varname]]) - NROW(tabulations2[which(datasetwide[[varname]] != "NA"), ])
 
    tabulations2_valid <- tabulations2[which(tabulations2$categories != "Missing - NA"), ]
    
    tabulations2_valid[[traitname]] <- droplevels(tabulations2_valid[[traitname]])
    
    print(tabulations2_valid)
    print(levels(tabulations2_valid[[traitname]]))
    
    ### plot percent comparisons in each class
    
    gg1 <- ggplot(tabulations2_valid, aes(get(traitname), valid_percent)) +
      geom_bar(aes(fill = class), position = "dodge", stat = "identity", alpha = 0.8) +
      scale_fill_brewer(palette = palette_choice) +
      ylab("Proportion of class") +
      theme_classic() +
      xlab(paste0("Response (N= ", NROW(datasetwide[which(datasetwide[[varname]] != "NA"), ]), ")")) +
      geom_text(aes(get(traitname),
        y = valid_percent + 0.03, group = class,
        label = format(valid_percent, nsmall = 0, digits = 2, scientific = FALSE, size = 4)
      ),
      color = "cornflowerblue", position = position_dodge(.9), hjust = .5
      ) +
      labs(fill = "Class", colour = "", title = paste0(tidy_mres_weighted_df$exposure[1]))


    png(
      filename = paste0(output_dir, "/", "percent_", traitname, "_", optimum_classes, "class_barplot.png"),
      width = 20, height = 13, units = "cm", res = 300
    )

    print(gg1)

    dev.off()

    ### plot of absolute numbers in each class

    gg2 <- ggplot(tabulations2_valid, aes(get(traitname), n)) +
      geom_bar(aes(fill = class), position = "dodge", stat = "identity", alpha = 0.8) +
      scale_fill_brewer(palette = palette_choice) +
      ylab("N in class") +
      theme_classic() +
      xlab(paste0("Response (N= ", NROW(datasetwide[which(datasetwide[[varname]] != "NA"), ]), ")")) +
      geom_text(aes(get(traitname),
        y = n + 500, group = class,
        label = format(n, nsmall = 0, digits = 2, scientific = FALSE, size = 4)
      ),
      color = "cornflowerblue", position = position_dodge(.9), hjust = (0.5)
      ) +
      labs(fill = "Class", colour = "", title = paste0(tidy_mres_weighted_df$exposure[1]))

    png(
      filename = paste0(output_dir, "/", "number_", traitname, "_", optimum_classes, "class_barplot.png"),
      width = 20, height = 13, units = "cm", res = 300
    )

    print(gg2)

    dev.off()

    if (missingnessplot == TRUE) {

      ### plot percent comparisons in each class

      gg1 <- ggplot(tabulations2, aes(get(traitname), percent)) +
        geom_bar(aes(fill = class), position = "dodge", stat = "identity", alpha = 0.8) +
        scale_fill_brewer(palette = palette_choice) +
        ylab("Proportion of class") +
        theme_classic() +
        xlab(paste0("Class (N= ", NROW(datasetwide[which(datasetwide[[varname]] != "NA"), ]), ")")) +
        geom_text(aes(get(traitname),
          y = percent + 0.03, group = class,
          label = format(percent, nsmall = 0, digits = 2, scientific = FALSE, size = 4)
        ),
        color = "cornflowerblue", position = position_dodge(.9), hjust = .5
        ) +
        labs(fill = "Class", colour = "", title = paste0(tidy_mres_weighted_df$exposure[1]))


      png(
        filename = paste0(output_dir, "/", "percent_NAs_", traitname, "_", optimum_classes, "class_barplot.png"),
        width = 20, height = 13, units = "cm", res = 300
      )

      print(gg1)

      dev.off()

      ### plot of absolute numbers in each class

      gg2 <- ggplot(tabulations2, aes(get(traitname), n)) +
        geom_bar(aes(fill =class ), position = "dodge", stat = "identity", alpha = 0.8) +
        scale_fill_brewer(palette = palette_choice) +
        ylab("N in class") +
        theme_classic() +
        xlab(paste0("Class (N= ", NROW(datasetwide[which(datasetwide[[varname]] != "NA"), ]), ")")) +
        geom_text(aes(get(traitname),
          y = n + 500, group = class,
          label = format(n, nsmall = 0, digits = 2, scientific = FALSE, size = 4)
        ),
        color = "cornflowerblue", position = position_dodge(.9), hjust = (0.5)
        ) +
        labs(fill = "Class", colour = "", title = paste0(tidy_mres_weighted_df$exposure[1]))

      png(
        filename = paste0(output_dir, "/", "number_NAs_", traitname, "_", optimum_classes, "class_barplot.png"),
        width = 20, height = 13, units = "cm", res = 300
      )

      print(gg2)

      dev.off()
      
    } else {}
  }

  res_ls <- list(plots_ls, multinomial_res_ls)

  return(res_ls)
}
