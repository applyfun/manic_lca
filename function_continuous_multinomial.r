### function for multinomial logistic regression with continuous predictor
### loop over either having class 1 or class 3 as reference class

analyse_continuous_multinomial <-
  function(datasetwide, datasetlong, varname, traitname) {

    ### initilise empty lists for saving and returning plots and results

    plots_ls <- list()
    multinomial_res_ls <- list()

    datasetlong[["predclass"]] <- as.factor(as.character(datasetlong[["predclass"]]))
    datasetwide[["predclass"]] <- as.factor(as.character(datasetwide[["predclass"]]))

    for (reference_i in c("main", "relevel")) {

      ### if loop is 'relevel' then change the reference class to "3"

      if (reference_i == "relevel") {
        datasetwide[["predclass"]] <-
          relevel(as.factor(datasetwide[["predclass"]]), ref = 3)

        datasetlong[["predclass"]] <-
          relevel(as.factor(datasetlong[["predclass"]]), ref = 3)
      } else {
        print("This is the main analysis (reference level '1')")
      }


      formula_i <-
        formula(paste0("predclass ~ ", varname))

      ### multinomial logistic regression naive to probabilities
      print(formula_i)

      mres <- multinom(formula_i, data = datasetwide)
      tidy_mres <-
        tidy(mres,
          conf.int = TRUE,
          conf.level = 0.95,
          exponentiate = TRUE
        )

      tidy_mres$model <- "naive"

      tidy_mres$exposure <- varname

      naive_res_df <- tidy_mres

      print("Finished estimating naive model")

      ### multinomial logistic regression weighted for probabilities (bias adjusted)

      mres_weighted <-
        multinom(formula_i, data = datasetlong, weights = probabilities)

      ### extract and exponentiate coefficients, add columns for labelling purposes

      tidy_mres_weighted <-
        tidy(mres_weighted,
          conf.int = TRUE,
          conf.level = 0.95,
          exponentiate = TRUE
        )

      tidy_mres_weighted$model <- "weighted"

      tidy_mres_weighted$exposure <- varname

      weighted_res_df <- tidy_mres_weighted

      print(weighted_res_df)

      print("Finished estimating weighted model")

      ### bind and drop intercept term from results

      naive_res_df <-
        naive_res_df[which(naive_res_df$term != "(Intercept)"), ]

      weighted_res_df <-
        weighted_res_df[which(weighted_res_df$term != "(Intercept)"), ]

      ### create a 'comparison' column of statistical test for plot labels
      naive_res_df$comparison <-
        paste0(
          levels(datasetlong$predclass)[1],
          "v",
          naive_res_df$`y.level`
        )

      weighted_res_df$comparison <-
        paste0(
          levels(datasetlong$predclass)[1],
          "v",
          naive_res_df$`y.level`
        )

      ### remove periods if any exist from names
      naive_res_df$exposure <-
        as.factor(gsub("\\.", " ", as.character(naive_res_df$exposure)))

      weighted_res_df$exposure <-
        as.factor(gsub("\\.", " ", as.character(weighted_res_df$exposure)))

      write.csv(weighted_res_df, paste0(output_dir,"/","multinomial_",traitname,"_weighted_res_",reference_i,".csv"), 
                col.names=T, row.names=F, quote = F)
      
      ### plot multinomial results

      gg01 <- ggplot(naive_res_df, aes(comparison, estimate, colour = comparison)) +
        geom_point(alpha = 0.8) +
        geom_hline(
          yintercept = 1,
          color = "coral",
          alpha = 0.4
        ) +
        geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
        xlab("Class comparisons") +
        scale_fill_brewer(palette = palette_choice) +
        ylab(paste0("Risk ratio (per SD increase in ", traitname, ")")) +
        theme_classic() +
        labs(title = paste0(naive_res_df$exposure[1])) +
        theme(
          panel.background = element_rect(fill = "white", colour = "black"),
          strip.background = element_rect(
            fill = "white", colour =
              "white"
          )
        )

      png(
        filename = paste0(
          output_dir,
          "/",
          "multinomial_naive_",
          reference_i,
          "_", traitname, "_",
          optimum_classes,
          "class_scatterplot.png"
        ),
        width = 20,
        height = 10,
        units = "cm",
        res = 300
      )

      print(gg01)

      dev.off()

      ###

      gg02 <- ggplot(
        weighted_res_df,
        aes(comparison, estimate, colour = comparison)
      ) +
        geom_point(alpha = 0.8) +
        geom_hline(
          yintercept = 1,
          color = "coral",
          alpha = 0.4
        ) +
        geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
        xlab("Class comparisons") +
        scale_fill_brewer(palette = palette_choice) +
        ylab(paste0("Risk ratio (per SD increase in ", traitname, ")")) +
        theme_classic() +
        labs(title = paste0(weighted_res_df$exposure[1])) +
        theme(
          panel.background = element_rect(fill = "white", colour = "black"),
          strip.background = element_rect(
            fill = "white", colour =
              "white"
          )
        )

      png(
        filename = paste0(
          output_dir,
          "/",
          "multinomial_weighted_",
          reference_i,
          "_", traitname, "_",
          optimum_classes,
          "class_scatterplot.png"
        ),
        width = 20,
        height = 10,
        units = "cm",
        res = 300
      )

      print(gg02)

      dev.off()

      ### bind results and plots in to lists to return

      plots_ls[[paste0("naive_", reference_i)]] <- gg01
      plots_ls[[paste0("weighted_", reference_i)]] <- gg02

      multinomial_res_ls[[paste0("naive_", reference_i)]] <- naive_res_df
      multinomial_res_ls[[paste0("weighted_", reference_i)]] <- weighted_res_df

      gc()
    }

    ### plot density plots of continuous variable by most likely class

    datasetwide[["predclass"]] <-
      relevel(as.factor(datasetwide[["predclass"]]), ref = 3)
    datasetwide[["predclass"]] <-
      relevel(as.factor(datasetwide[["predclass"]]), ref = 2)
    datasetwide[["predclass"]] <-
      relevel(as.factor(datasetwide[["predclass"]]), ref = 1)
    
    gg1 <-
      ggplot(datasetwide, aes(x = get(varname), fill = predclass)) +
      geom_density(alpha = 0.3) +
      xlab(paste0(gsub("\\.", " ", as.character(varname)), " N=", NROW(datasetwide))) +
      ggtitle(paste0(gsub("\\.", " ", as.character(varname)))) +
      labs(fill = "Class") +
      theme_classic()

    png(
      filename = paste0(
        output_dir,
        "/",
        "",
        traitname,
        "_",
        optimum_classes,
        "class_densityplot.png"
      ),
      width = 20,
      height = 15,
      units = "cm",
      res = 300
    )

    print(gg1)

    dev.off()


    res_ls <- list(plots_ls, multinomial_res_ls)

    return(res_ls)
  }
