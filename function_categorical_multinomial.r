analyse_categorical_multinomial <- function(datasetwide,
                                            datasetlong,
                                            varname,
                                            traitname) {

  ### initilise empty lists for saving and returning plots and results

  plots_ls <- list()
  multinomial_res_ls <- list()

  ### dummy code (create comparions) of categorical var - split in to binary vars

  datasetwide_spl <- datasetwide[which(!is.na(datasetwide[[varname]])), ]

  res_dummy <- as.data.frame(model.matrix(formula(paste0("~", varname, " -1")), data = datasetwide_spl))

  names(res_dummy) <- unlist(lapply(strsplit(names(res_dummy), varname), "[[", 2))

  datasetwide <- cbind(datasetwide_spl, res_dummy)

  ### repeat for long dataset

  datasetlong_spl <- datasetlong[which(!is.na(datasetlong[[varname]])), ]

  res_dummy <- as.data.frame(model.matrix(formula(paste0("~", varname, " -1")), data = datasetlong_spl))

  names(res_dummy) <- unlist(lapply(strsplit(names(res_dummy), varname), "[[", 2))

  datasetlong <- cbind(datasetlong_spl, res_dummy)

  print(names(res_dummy))

  ### loop over either catgegory 1 or cat 3 as the reference category

  for (reference_i in c("main", "relevel")) {
    naive_res <- list()
    weighted_res <- list()

    datasetlong[["predclass"]] <- as.factor(as.character(datasetlong[["predclass"]]))

    datasetwide[["predclass"]] <- as.factor(as.character(datasetwide[["predclass"]]))

    ### if loop is 'relevel' then change the reference class to "3"

    if (reference_i == "relevel") {
      datasetwide[["predclass"]] <-
        relevel(as.factor(datasetwide[["predclass"]]), ref = 3)

      datasetlong[["predclass"]] <-
        relevel(as.factor(datasetlong[["predclass"]]), ref = 3)
    } else {
      print("This is the main analysis (reference level '1')")
    }

    ### create formula for multinomial regression

    for (varname_levels in names(res_dummy)) {
      formula_i <- formula(paste0("predclass ~ ", varname_levels))

      print(formula_i)

      ### multinomial logistic regression naive to probabilities
      mres <- multinom(formula_i, data = datasetwide)

      tidy_mres <- tidy(mres, conf.int = TRUE, conf.level = 0.95, exponentiate = TRUE)

      tidy_mres$model <- "naive"
      tidy_mres$exposure <- varname_levels
      naive_res[[varname_levels]] <- tidy_mres

      ### multinomial logistic regression weighted for probabilities
      mres_weighted <- multinom(formula_i, data = datasetlong, weights = probabilities)

      # extract and exponentiate coefficients
      tidy_mres_weighted <- tidy(mres_weighted, conf.int = TRUE, conf.level = 0.95, exponentiate = TRUE)
      tidy_mres_weighted$model <- "weighted"
      tidy_mres_weighted$exposure <- varname_levels

      weighted_res[[varname_levels]] <- tidy_mres_weighted
    }

    ### bind and drop intercept from results
    naive_res_df <- bind_rows(naive_res)
    naive_res_df <- naive_res_df[which(naive_res_df$term != "(Intercept)"), ]

    weighted_res_df <- bind_rows(weighted_res)
    weighted_res_df <- weighted_res_df[which(weighted_res_df$term != "(Intercept)"), ]

    ### create a 'comparison' column of statistical test for plot
    
    naive_res_df$comparison <-
      paste0(levels(datasetwide[["predclass"]])[1], "v", naive_res_df$`y.level`)
    
    weighted_res_df$comparison <-
      paste0(levels(datasetlong[["predclass"]])[1], "v", weighted_res_df$`y.level`)
    
    ### remove periods if any exist from names
    
    naive_res_df$exposure <- as.factor(gsub("\\.", " ", as.character(naive_res_df$exposure)))
    weighted_res_df$exposure <- as.factor(gsub("\\.", " ", as.character(weighted_res_df$exposure)))

    ### bind results in to lists to return

    multinomial_res_ls[[paste0("naive_", reference_i)]] <- naive_res_df
    multinomial_res_ls[[paste0("weighted_", reference_i)]] <- weighted_res_df

    print(weighted_res_df)

    write.csv(weighted_res_df, paste0(output_dir,"/","multinomial_",traitname,"_weighted_res_",reference_i,".csv"), 
              col.names=T, row.names=F, quote = F)
    
    ### plot multinomial results

    gg01 <- ggplot(naive_res_df, aes(comparison, estimate, colour = comparison)) +
      geom_point(alpha = 0.8) +
      geom_hline(yintercept = 1, color = "coral", alpha = 0.4) +
      geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
      xlab("Class comparisons") +
      ylab("Estimate (RR)") +
      facet_wrap(~exposure) +
      scale_fill_brewer(palette = palette_choice) +
      theme_classic() +
      theme(
        panel.background = element_rect(fill = "white", colour = "black"),
        strip.background = element_rect(fill = "white", colour = "white")
      ) +
      labs(colour = "")

    png(
      filename = paste0(
        output_dir, "/", "multinomial_naive_", traitname,
        "_", reference_i, "_", optimum_classes, "class_scatterplot.png"
      ),
      width = 22, height = 11, units = "cm", res = 300
    )

    print(gg01)

    dev.off()

    ###
    gg02 <- ggplot(weighted_res_df, aes(comparison, estimate, colour = comparison)) +
      geom_point(alpha = 0.8) +
      geom_hline(yintercept = 1, color = "coral", alpha = 0.4) +
      geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
      xlab("Class comparisons") +
      ylab("Estimate (RR)") +
      facet_wrap(~exposure) +
      scale_fill_brewer(palette = palette_choice) +
      theme_classic() +
      theme(
        panel.background = element_rect(fill = "white", colour = "black"),
        strip.background = element_rect(fill = "white", colour = "white")
      ) +
      labs(colour = "")

    png(
      filename = paste0(
        output_dir, "/", "multinomial_weighted_", traitname,
        "_", reference_i, "_", optimum_classes, "class_scatterplot.png"
      ),
      width = 22, height = 11, units = "cm", res = 300
    )

    print(gg02)

    dev.off()

    ### bind plots in to lists to return

    plots_ls[[paste0("naive_", reference_i)]] <- gg01
    plots_ls[[paste0("weighted_", reference_i)]] <- gg02

    print("Plotted mutinomial results!")

    ### summarise by most probable class

    tabulations <- list()

    for (i in 1:optimum_classes) {
      tb1 <- tabyl(datasetwide[[varname]][datasetwide$predclass == i])

      names(tb1)[1] <- paste0("class_", i)

      tb1$categories <- tb1[, 1]

      tabulations[[i]] <- tb1
    }

    tabulations2 <- bind_rows(tabulations, .id = c("class"))[c("categories", "n", "percent", "class")]

    tabulations2$`N (percent)` <-
      paste0(tabulations2$n, " (", round(tabulations2$percent*100, digits = 2), ")")

    tabulations2$exposure <- as.factor(gsub("\\.", " ", as.character(tabulations2$categories))) # remove periods if any exist
    
    tabulations2$categories <- as.factor(gsub("\\.", " ", as.character(tabulations2$categories))) # remove periods if any exist
    
    ### tabulate number of missing rows (hence individuals with NA for exposure) and subset
    missing <- NROW(datasetwide[[varname]]) - NROW(tabulations2[which(datasetwide[[varname]] != "NA"), ])
    tabulations2 <- tabulations2[which(tabulations2$categories != "NA"), ]

    print(tabulations2)

    write.table(tabulations2[c("categories", "n", "percent", "class", "N (percent)")],
                file = paste0(output_dir, "/", traitname, "_", optimum_classes, "class_tabulations.csv"),
                row.names = FALSE, quote = FALSE, sep = ","
    )
    
    ### plot raw comparisons by most probable class
    png(
      filename = paste0(
        output_dir, "/", "percent_", traitname,
        "_", optimum_classes, "class_barplot.png"
      ),
      width = 20, height = 13, units = "cm", res = 300
    )

    print(ggplot(tabulations2, aes(class, percent)) +
      geom_bar(aes(fill = exposure), position = "dodge", stat = "identity", alpha = 0.8) +
      scale_fill_brewer(palette = palette_choice) +
      ylab("Proportion of class") +
      xlab(paste0("Class (N= ", NROW(datasetwide[which(datasetwide[[varname_levels]] != "NA"), ]), ")")) +
      labs(fill = gsub("\\.", " ", paste0(varname))) +
      theme_classic() +
      geom_text(aes(class,
        y = percent + 0.03,
        group = exposure,
        label = format(
          percent,
          nsmall = 0,
          digits = 1,
          scientific = FALSE, size = 4
        )
      ),
      color = "cornflowerblue",
      position = position_dodge(.9),
      hjust = .5
      ))

    dev.off()
  }

  res_ls <- list(plots_ls, multinomial_res_ls)

  return(res_ls)
}
