# ==========================================================================
#         DO WHAT THE F**K YOU WANT TO PUBLIC LICENSE 
#              Version 2, December 2004 
#
#     Copyright (C)2026 Émile Caron <emile.caron@umontreal.ca> &
#                  2026 Sébastien Béland <sebastien.beland@umontreal.ca>
#
#     Everyone is permitted to copy and distribute verbatim or modified 
#     copies of this license document, and changing it is allowed as long 
#     as the name is changed. 
#
#     DO WHAT THE F**K YOU WANT TO PUBLIC LICENSE 
#     TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION 
#
#     0. You just DO WHAT THE F**K YOU WANT TO.
# =========================================================================

library(shiny)
library(difR)
library(markdown)
library(shinythemes)

# 1. Génération des données

# Paramètres "vrais" des items pour les deux groupes (dichotomique, 2PL)
make_item_params <- function(n_items, seed, dif_items, type_dif, delta_b, delta_a) {
  set.seed(seed)
  a_ref <- round(runif(n_items, 0.6, 1.6), 2)          # discriminations
  b_ref <- round(seq(-1.5, 1.5, length.out = n_items), 2) # difficultés étagées
  a_foc <- a_ref
  b_foc <- b_ref
  for (i in dif_items) {
    if (type_dif %in% c("uniforme", "mixte")) b_foc[i] <- b_ref[i] + delta_b
    if (type_dif %in% c("non_uniforme", "mixte")) a_foc[i] <- a_foc[i] + delta_a
  }
  data.frame(
    item = paste0("Item", 1:n_items),
    a_ref = a_ref, b_ref = b_ref,
    a_foc = a_foc, b_foc = b_foc,
    DIF = 1:n_items %in% dif_items,
    type = ifelse(1:n_items %in% dif_items, type_dif, "")
  )
}

# Simulations des réponses dichotomiques : P(X=1) = 1 / (1 + exp(-a (theta - b)))
gen_responses <- function(params, n_ref, n_foc, mu_foc, seed) {
  set.seed(seed + 1)
  theta_ref <- rnorm(n_ref, 0, 1)
  theta_foc <- rnorm(n_foc, mu_foc, 1)
  n_items <- nrow(params)
  p2 <- function(th, a, b) plogis(outer(th, b, "-") * matrix(a, length(th), n_items, byrow = TRUE))
  X_ref <- matrix(runif(n_ref * n_items), n_ref) < p2(theta_ref, params$a_ref, params$b_ref)
  X_foc <- matrix(runif(n_foc * n_items), n_foc) < p2(theta_foc, params$a_foc, params$b_foc)
  Data <- as.data.frame(rbind(
    if (n_ref > 0) X_ref * 1 else NULL,
    if (n_foc > 0) X_foc * 1 else NULL
  ))
  names(Data) <- paste0("Item", 1:n_items)
  Data$group <- c(rep(0, n_ref), rep(1, n_foc))  # 0 = référence, 1 = focal
  Data
}

# Génération polytomique "maison" (modèle GRM, réponse graduée)
# P(X >= k | thêta) = plogis(a * (thêta - b - seuil_k))  pour k = 1 ... ncat-1
# Uniform DIF = décalage de b, non-uniforme = variation de a (comme en 2PL)
gen_poly_data <- function(n_items, seed, dif_items, type_dif, delta_b, delta_a,
                          n_ref, n_foc, mu_foc, ncat) {
  set.seed(seed)
  a <- round(runif(n_items, 0.6, 1.6), 2)            # pentes
  b <- round(seq(-1.5, 1.5, length.out = n_items), 2) # difficultés
  # Seuils de catégories : répartis autour de b, croissants (convention GRM :
  # pour que P(X >= 1) >= P(X >= 2) >= ..., il faut s_1 < s_2 < ...)
  steps <- round(seq(-1.2, 1.2, length.out = ncat - 1), 2)
  a_foc <- a; b_foc <- b
  for (i in dif_items) {
    if (type_dif %in% c("uniforme", "mixte")) b_foc[i] <- b[i] + delta_b
    if (type_dif %in% c("non_uniforme", "mixte")) a_foc[i] <- a_foc[i] + delta_a
  }
  # Simulation d'un score GRM : catégorie = nb de seuils franchis
  sim_group <- function(n, mu, a_g, b_g) {
    th <- rnorm(n, mu, 1)
    scores <- matrix(0L, n, n_items)
    for (k in seq_along(steps)) {
      pk <- t(sapply(th, function(x) plogis(a_g * (x - b_g - steps[k]))))
      scores <- scores + (matrix(runif(n * n_items), n) < pk)
    }
    scores
  }
  X_ref <- sim_group(n_ref, 0, a, b)
  X_foc <- sim_group(n_foc, mu_foc, a_foc, b_foc)
  Data <- as.data.frame(rbind(X_ref, X_foc))
  names(Data) <- paste0("Item", 1:n_items)
  Data$group <- c(rep(0, n_ref), rep(1, n_foc))
  params <- data.frame(
    item = paste0("Item", 1:n_items),
    a_ref = a, b_ref = b, a_foc = a_foc, b_foc = b_foc,
    DIF = 1:n_items %in% dif_items,
    type = ifelse(1:n_items %in% dif_items, type_dif, "")
  )
  list(params = params, Data = Data, dif_items = dif_items,
       steps = steps, ncat = ncat)
}

# Extraction robuste des items signalés par une méthode difR
get_flags <- function(res) {
  d <- res$DIFitems
  if (is.null(d) || is.character(d)) integer(0) else as.integer(d)
}

# Appel robuste : essaie avec les arguments complets, retombe sur le minimum
try_call <- function(f, full, minimal) {
  tryCatch(do.call(f, full),
           error = function(e) tryCatch(do.call(f, minimal),
                                        error = function(e2) e))
}

# Lancement d'une méthode difR selon le type d'items
run_difr_method <- function(meth, Data, alpha, purify, item_type) {
  common <- list(Data = Data, group = "group", focal.name = 1,
                 alpha = alpha, purify = purify)
  minimal <- list(Data = Data, group = "group", focal.name = 1, alpha = alpha)
  f <- switch(meth,
              # Méthodes dichotomiques
              MH = function() do.call(difMH, common),
              ILHTE = function() try_call(ILHTEdif, common, minimal),
              Logistic = function() do.call(difLogistic, c(common, list(type = "both"))),
              Lord = function() do.call(difLord, c(common, list(model = "2PL", engine = "ltm"))),
              Raju = function() do.call(difRaju, c(common, list(model = "2PL", engine = "ltm"))),
              SIBTEST = function() do.call(difSIBTEST, common),
              Std = function() do.call(difStd, common),
              # Méthodes polytomiques
              MantelPoly = function() try_call(difMantelPoly, common, minimal),
              # Quade : conçue pour le polytomique, mais offerte aussi en
              # dichotomique à titre didactique (avertissement côté server)
              Quade = function() try_call(difQuade, common, minimal),
              PolyLogistic = function() try_call(difPolyLogistic,
                                                 c(common, list(type = "both")),
                                                 c(minimal, list(type = "both"))),
              NULL)
  if (is.null(f)) stop("Unknown method: ", meth)
  if (item_type == "poly" && meth %in% c("MH", "ILHTE", "Logistic", "Lord", "Raju", "SIBTEST", "Std"))
    stop(meth, " is not available for polytomous items.")
  if (item_type == "dicho" && meth %in% c("MantelPoly", "PolyLogistic"))
    stop(meth, " is not available for dichotomous items.")
  # NB : difQuade est offerte dans les deux modes (outil didactique) ; le
  # server affiche un avertissement si elle est utilisée sur du dichotomique.
  tryCatch(f(), error = function(e) e)
}

# 2. UI

ui <- fluidPage(
  theme = shinytheme("spacelab"),
  # Logo dans le titre : www/logo.png s'il existe, sinon l'icône Font Awesome
  titlePanel(
    if (file.exists("www/logo.png")) {
      tagList(tags$img(src = "logo.png", height = "40px",
                       alt = "logo", style = "vertical-align:middle; margin-right:10px;"),
              "DIF Simulator with difR")
    } else {
      tagList(icon("chart-line"), "DIF Simulator with difR")
    }
  ),
  sidebarLayout(
    sidebarPanel(
      h4(tagList(icon("database"), "Data source")),
      radioButtons("data_source", "Source",
                   choices = c("Dichotomous (simulated)" = "sim",
                               "Polytomous (simulated)" = "simpoly",
                               "Upload your own data (CSV)" = "upload")),
      conditionalPanel(
        condition = "input.data_source == 'sim'",
        h4(tagList(icon("dice"), "Simulated data")),
        numericInput("seed", "Random seed", 0, min = 0),
        sliderInput("n_items", "Number of items", 5, 40, 20, step = 1),
        numericInput("n_ref", "Respondents — reference group", 500, min = 50),
        numericInput("n_foc", "Respondents — focal group", 500, min = 50),
        sliderInput("mu_foc", "Mean ability of the focal group (theta)",
                    -1.5, 1.5, 0, step = 0.1),
        hr(),
        h4(tagList(icon("syringe"), "Injected DIF (ground truth)")),
        sliderInput("n_dif", "Number of items with DIF", 0, 10, 3),
        selectInput("type_dif", "Type of DIF",
                    c("uniform" = "uniforme",
                      "non-uniform" = "non_uniforme",
                      "mixed" = "mixte")),
        sliderInput("delta_b", "Uniform magnitude (b shift)",
                    -2, 2, 0.6, step = 0.1),
        sliderInput("delta_a", "Non-uniform magnitude (a change)",
                    -1, 1, 0.5, step = 0.1)
      ),
      conditionalPanel(
        condition = "input.data_source == 'simpoly'",
        h4(tagList(icon("dice"), "Simulated data (polytomous)")),
        tags$p(style = "color:#666; font-size:0.85em;",
               "Generated under a graded response model (GRM)."),
        numericInput("seed_poly", "Random seed", 0, min = 0),
        sliderInput("n_items_poly", "Number of items", 5, 40, 20, step = 1),
        sliderInput("ncat", "Number of response categories", 2, 6, 3, step = 1),
        numericInput("n_ref_poly", "Respondents — reference group", 500, min = 50),
        numericInput("n_foc_poly", "Respondents — focal group", 500, min = 50),
        sliderInput("mu_foc_poly", "Mean ability of the focal group (theta)",
                    -1.5, 1.5, 0, step = 0.1),
        hr(),
        h4(tagList(icon("syringe"), "Injected DIF (ground truth)")),
        sliderInput("n_dif_poly", "Number of items with DIF", 0, 10, 3),
        selectInput("type_dif_poly", "Type of DIF",
                    c("uniform" = "uniforme",
                      "non-uniform" = "non_uniforme",
                      "mixed" = "mixte")),
        sliderInput("delta_b_poly", "Uniform magnitude (b shift)",
                    -2, 2, 0.6, step = 0.1),
        sliderInput("delta_a_poly", "Non-uniform magnitude (a change)",
                    -1, 1, 0.5, step = 0.1)
      ),
      conditionalPanel(
        condition = "input.data_source == 'upload'",
        h4(tagList(icon("upload"), "Uploaded data")),
        fileInput("file", "Upload a CSV file",
                  accept = c(".csv", ".txt")),
        tags$p(style = "color:#666; font-size:0.85em;",
               "One row per respondent, one column per item (0/1 for dichotomous, ",
               "0 to K-1 for polytomous), plus one column coding the group ",
               "(e.g. 0 = reference, 1 = focal). ",
               "Rows with missing values (NA) are removed automatically."),
        selectInput("item_type_upload", "Item type",
                    c("Dichotomous" = "dicho",
                      "Polytomous (ordinal)" = "poly")),
        checkboxInput("header", "File has a header row", TRUE),
        selectInput("sep", "Field separator",
                    c("Comma" = ",", "Semicolon" = ";", "Tab" = "\t")),
        uiOutput("group_col_ui"),
        selectInput("focal_name", "Value of the focal group",
                    choices = c("1")),
        tags$p(style = "color:#666; font-size:0.85em;",
               "No group variable in the file? Two options:"),
        checkboxInput("split_random", "Randomly split respondents into two groups", FALSE),
        sliderInput("prop_foc", "Proportion in the focal group",
                    0.1, 0.9, 0.5, step = 0.05)
      ),
      hr(),
      h4(tagList(icon("magnifying-glass-chart"), "Detection")),
      # Liste des méthodes adaptée au type d'items
      uiOutput("methods_ui"),
      sliderInput("alpha", "Alpha level", 0.01, 0.10, 0.05, step = 0.01),
      checkboxInput("purify", "Iterative purification", FALSE),
      actionButton("go", tagList(icon("play"), "Generate and analyse"),
                   class = "btn-primary"),
      tags$p(style = "color:#999; font-size:0.85em;",
             "Lord and Raju require the 'ltm' package.")
    ),
    mainPanel(
      tabsetPanel(
        id = "tabs",
        tabPanel(tagList(icon("book"), "Documentation"), value = "docs",
                 uiOutput("docs")
        ),
        tabPanel(tagList(icon("bullseye"), "Ground truth"), value = "groundtruth",
                 h4("Simulated item parameters"),
                 tableOutput("param_table"),
                 plotOutput("icc_plot"),
                 sliderInput("item_curves", "Item to display (ICC by group)",
                             1, 40, 1, step = 1),
                 plotOutput("score_plot"),
                 plotOutput("scatter_plot")
        ),
        tabPanel(tagList(icon("search"), "Detection"), value = "detection",
                 h4("Method verdicts vs ground truth"),
                 tableOutput("verdict_table"),
                 verbatimTextOutput("confusion"),
                 plotOutput("method_plot")
        ),
        tabPanel(tagList(icon("terminal"), "Detailed outputs"), value = "details",
                 verbatimTextOutput("raw_output")
        )
      )
    )
  )
)

# 3. Server

server <- function(input, output, session) {
  
  # Type d'items courant (dicho / poly) selon la source
  item_type <- reactive({
    switch(input$data_source,
           sim = "dicho", simpoly = "poly",
           upload = input$item_type_upload)
  })
  
  # Liste des méthodes disponibles selon le type d'items
  output$methods_ui <- renderUI({
    if (item_type() == "dicho") {
      checkboxGroupInput("methods", "difR methods (dichotomous)",
                         choices = c("MH", "Logistic", "ILHTE",
                                     "Lord", "Raju", "SIBTEST", "Std",
                                     "Quade"),
                         selected = c("MH", "Logistic", "Lord"))
    } else {
      checkboxGroupInput("methods", "difR methods (polytomous)",
                         choices = c("MantelPoly", "Quade", "PolyLogistic"),
                         selected = c("MantelPoly", "PolyLogistic"))
    }
  })
  
  # Lecture du fichier téléversé
  uploaded <- reactive({
    req(input$file)
    df <- tryCatch(
      read.csv(input$file$datapath, header = input$header, sep = input$sep,
               stringsAsFactors = FALSE, na.strings = c("NA", "")),
      error = function(e) NULL
    )
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df
  })
  
  # Choix dynamique de la colonne de groupe
  output$group_col_ui <- renderUI({
    df <- uploaded()
    if (is.null(df)) return(NULL)
    cols <- names(df)
    selectInput("group_col", "Column coding the group (optional)",
                choices = c("(none — no group variable)" = "", cols))
  })
  
  # Détection automatique des valeurs possibles pour le groupe focal
  observe({
    df <- uploaded()
    if (is.null(df)) return
    req(input$group_col)
    vals <- sort(unique(na.omit(df[[input$group_col]])))
    updateSelectInput(session, "focal_name",
                      choices = as.character(vals),
                      selected = as.character(vals[length(vals)]))
  })
  
  # Jeu de données sans NA + vecteur de groupe (relancé par le bouton)
  clean_data <- eventReactive(input$go, {
    if (input$data_source == "sim") {
      n_items <- input$n_items
      dif_items <- if (input$n_dif > 0) {
        sort(sample(seq_len(n_items), min(input$n_dif, n_items)))
      } else integer(0)
      params <- make_item_params(n_items, input$seed, dif_items,
                                 input$type_dif, input$delta_b, input$delta_a)
      Data <- gen_responses(params, input$n_ref, input$n_foc,
                            input$mu_foc, input$seed)
      # Étendue du slider d'items selon le nombre d'items
      updateSliderInput(session, "item_curves", max = n_items,
                        value = min(input$item_curves, n_items))
      list(params = params, Data = Data, dif_items = dif_items,
           steps = NULL, n_removed = 0,
           simulated = TRUE, item_type = "dicho")
    } else if (input$data_source == "simpoly") {
      n_items <- input$n_items_poly
      dif_items <- if (input$n_dif_poly > 0) {
        sort(sample(seq_len(n_items), min(input$n_dif_poly, n_items)))
      } else integer(0)
      out <- gen_poly_data(n_items, input$seed_poly, dif_items,
                           input$type_dif_poly, input$delta_b_poly,
                           input$delta_a_poly,
                           input$n_ref_poly, input$n_foc_poly,
                           input$mu_foc_poly, input$ncat)
      updateSliderInput(session, "item_curves", max = n_items,
                        value = min(input$item_curves, n_items))
      list(params = out$params, Data = out$Data, dif_items = out$dif_items,
           steps = out$steps, ncat = out$ncat, n_removed = 0,
           simulated = TRUE, item_type = "poly")
    } else {
      df <- uploaded()
      if (is.null(df)) stop("Please upload a valid CSV file first.")
      # Suppression des lignes avec NA
      df <- df[stats::complete.cases(df), , drop = FALSE]
      n_removed <- nrow(uploaded()) - nrow(df)
      if (nrow(df) < 10) stop("Not enough complete rows left after removing NA.")
      # Groupe : colonne du fichier ou répartition aléatoire
      if (isTRUE(input$split_random) || input$group_col == "") {
        n <- nrow(df)
        grp <- rep(0, n)
        grp[sample(n, round(input$prop_foc * n))] <- 1
      } else {
        gcol <- input$group_col
        grp <- ifelse(as.character(df[[gcol]]) == input$focal_name, 1, 0)
        if (length(unique(grp)) < 2) stop("The selected group column does not define two groups.")
        df[[gcol]] <- NULL
      }
      # Binarisation des items en 0/1 (uniquement en mode dichotomique)
      items <- names(df)
      if (input$item_type_upload == "dicho") {
        bad <- vapply(items, function(j) !all(df[[j]] %in% c(0, 1)), logical(1))
        if (any(bad)) {
          for (j in items[bad]) df[[j]] <- ifelse(as.numeric(df[[j]]) >= 1, 1, 0)
        }
      } else {
        # En polytomique : on force les scores au format numérique entier
        for (j in items) df[[j]] <- suppressWarnings(as.integer(df[[j]]))
        if (any(vapply(items, function(j) any(is.na(df[[j]])), logical(1))))
          stop("Some item columns are not numeric.")
      }
      Data <- df
      Data$group <- grp
      updateSliderInput(session, "item_curves", max = length(items),
                        value = min(input$item_curves, length(items)))
      list(params = NULL, Data = Data, dif_items = integer(0),
           steps = NULL, n_removed = n_removed,
           simulated = FALSE, item_type = input$item_type_upload)
    }
  })
  
  # Au clic sur "Generate and analyse" : bascule sur l'onglet Ground truth
  observeEvent(input$go, {
    updateTabsetPanel(session, "tabs", selected = "groundtruth")
  })
  
  results <- eventReactive(input$go, {
    d <- clean_data()
    # Avertissement didactique : Quade sur données dichotomiques
    if (d$item_type == "dicho" && "Quade" %in% input$methods) {
      showNotification(
        "Note: difQuade was designed for polytomous (ordinal) items. The
         results on dichotomous data are provided for experimentation
         purposes only — interpret with caution.",
        type = "warning", duration = NA)  # duration = NA : reste affiché
    }
    lapply(input$methods, function(m)
      list(method = m, res = run_difr_method(m, d$Data, input$alpha,
                                             input$purify, d$item_type)))
  })
  
  # --- Vérité terrain -------------------------------------------------
  output$param_table <- renderTable({
    d <- clean_data()
    if (d$simulated) {
      if (d$item_type == "poly") {
        # Version compacte pour le polytomique : a, b par groupe + type de DIF
        p <- d$params
        data.frame(item = p$item, a_ref = p$a_ref, b_ref = p$b_ref,
                   a_foc = p$a_foc, b_foc = p$b_foc,
                   DIF = ifelse(p$DIF, "DIF", "—"),
                   type = p$type)
      } else d$params
    } else {
      data.frame(Message = "Parameters are unknown for uploaded data.")
    }
  }, rownames = FALSE)
  
  output$icc_plot <- renderPlot({
    d <- clean_data()
    if (!d$simulated) return(invisible(NULL))
    i <- input$item_curves
    p <- d$params[i, ]
    th <- seq(-4, 4, length.out = 200)
    
    if (d$item_type == "dicho") {
      # 2PL : courbe unique par groupe
      curve_ref <- plogis(p$a_ref * (th - p$b_ref))
      curve_foc <- plogis(p$a_foc * (th - p$b_foc))
      matplot(th, cbind(curve_ref, curve_foc), type = "l", lty = 1, lwd = 2,
              col = c("steelblue", "firebrick"), ylim = c(0, 1),
              xlab = "Theta (ability)", ylab = "P(correct response)",
              main = sprintf("%s — %s", p$item,
                             if (p$DIF) paste("item with DIF (", p$type, ")") else "item without DIF"))
      abline(v = 0, lty = 3, col = "grey")
      legend("topleft", c("Reference group", "Focal group"),
             col = c("steelblue", "firebrick"), lwd = 2, bty = "n")
    } else {
      # GRM : toutes les courbes dans un seul graphique
      # Couleur = catégorie (comme plot(grm) de ltm), style = groupe
      # (plein = référence, pointillé = focal)
      steps <- d$steps   # seuils k = 1 ... ncat-1
      thr_ref <- p$b_ref + steps
      thr_foc <- p$b_foc + steps
      # Courbes cumulées par groupe (matrices 200 x (ncat-1))
      cum_ref <- sapply(thr_ref, function(bk) plogis(p$a_ref * (th - bk)))
      cum_foc <- sapply(thr_foc, function(bk) plogis(p$a_foc * (th - bk)))
      # Probabilités par catégorie : P(X=k) = cumulée_k - cumulée_{k+1}
      K <- length(steps) + 1
      cat_ref <- cbind(1 - cum_ref[, 1], if (K > 2) cum_ref[, 1:(K-2)] - cum_ref[, 2:(K-1)] else NULL, cum_ref[, K-1])
      cat_foc <- cbind(1 - cum_foc[, 1], if (K > 2) cum_foc[, 1:(K-2)] - cum_foc[, 2:(K-1)] else NULL, cum_foc[, K-1])
      # Palette : une couleur par catégorie (max 6 catégories)
      cat_cols <- c("#000000", "#E41A1C", "#4DAF4A", "#377EB8", "#FF7F00", "#984EA3")
      cat_cols <- cat_cols[1:K]
      matplot(th, cat_ref, type = "l", lty = 1, lwd = 2, col = cat_cols,
              ylim = c(0, 1), xlab = "Theta (ability)",
              ylab = "Probability",
              main = sprintf("%s (GRM) — %s", p$item,
                             if (p$DIF) paste("item with DIF (", p$type, ")") else "item without DIF"))
      matlines(th, cat_foc, lty = 2, lwd = 2, col = cat_cols)
      # Légende : catégories en couleur + groupes en style de trait
      legend("topright",
             c(paste0("Category ", 0:(K-1)), "",
               "Reference (solid)", "Focal (dashed)"),
             col = c(cat_cols, "white", "grey40", "grey40"),
             lty = c(rep(1, K), 0, 1, 2), lwd = 2, bty = "n", cex = 0.85,
             seg.len = 2.5)
    }
  })
  
  output$score_plot <- renderPlot({
    d <- clean_data()
    dd <- d$Data
    sc_ref <- rowSums(dd[dd$group == 0, -ncol(dd), drop = FALSE])
    sc_foc <- rowSums(dd[dd$group == 1, -ncol(dd), drop = FALSE])
    if (length(sc_ref) == 0 || length(sc_foc) == 0) return(invisible(NULL))
    ymax <- max(c(density(sc_ref)$y, density(sc_foc)$y)) * 1.1
    plot(density(sc_ref), col = "steelblue", lwd = 2, ylim = c(0, ymax),
         main = paste0("Observed score distributions",
                       if (!d$simulated) paste0("  (", d$n_removed, " rows with NA removed)") else ""),
         xlab = "Total score")
    lines(density(sc_foc), col = "firebrick", lwd = 2)
    legend("topright",
           c(sprintf("Reference (mean %.1f)", mean(sc_ref)),
             sprintf("Focal (mean %.1f)", mean(sc_foc))),
           col = c("steelblue", "firebrick"), lwd = 2, bty = "n")
  })
  
  # Nuage de points : score moyen par item, groupe focal vs référence
  output$scatter_plot <- renderPlot({
    d <- clean_data()
    dd <- d$Data
    items <- names(dd)[names(dd) != "group"]
    ref <- dd[dd$group == 0, items, drop = FALSE]
    foc <- dd[dd$group == 1, items, drop = FALSE]
    m_ref <- colMeans(ref)
    m_foc <- colMeans(foc)
    # Couleur : items signalés par les méthodes cochées
    rs <- results()
    flagged <- integer(0)
    for (r in rs) if (!inherits(r$res, "error")) flagged <- union(flagged, get_flags(r$res))
    col <- ifelse(seq_along(items) %in% flagged, "firebrick", "steelblue")
    plot(m_ref, m_foc, col = col, pch = 19, cex = 1.4,
         xlim = range(c(m_ref, m_foc)), ylim = range(c(m_ref, m_foc)),
         xlab = "Mean item score — reference group",
         ylab = "Mean item score — focal group",
         main = "Scatter plot of mean item scores (off-diagonal = DIF suspects)")
    abline(0, 1, lty = 2, col = "grey60")
    text(m_ref, m_foc, labels = seq_along(items), pos = 4, cex = 0.7, col = col)
    legend("bottomright",
           c("Not flagged", "Flagged as DIF"),
           col = c("steelblue", "firebrick"), pch = 19, bty = "n")
  })
  
  # --- Détection ------------------------------------------------------
  output$verdict_table <- renderTable({
    d <- clean_data(); rs <- results()
    n_items <- if (d$simulated) nrow(d$params) else sum(names(d$Data) != "group")
    verdicts <- lapply(rs, function(r) {
      flags <- if (inherits(r$res, "error")) NA else get_flags(r$res)
      ifelse(seq_len(n_items) %in% flags, "DIF", "—")
    })
    truth <- if (d$simulated) ifelse(d$params$DIF, "DIF", "—") else rep("?", n_items)
    out <- data.frame(Item = if (d$simulated) d$params$item else paste0("Item", 1:n_items),
                      Truth = truth, stringsAsFactors = FALSE)
    for (k in seq_along(rs)) out[[rs[[k]]$method]] <- verdicts[[k]]
    out
  }, rownames = FALSE)
  
  output$confusion <- renderText({
    d <- clean_data(); rs <- results()
    if (!d$simulated) {
      # Données réelles : pas de vérité terrain, on liste juste les items signalés
      txt <- vapply(rs, function(r) {
        if (inherits(r$res, "error")) {
          sprintf("%-12s : ERROR — %s", r$method, conditionMessage(r$res))
        } else {
          sprintf("%-12s : %d item(s) flagged as DIF", r$method, length(get_flags(r$res)))
        }
      }, character(1))
      return(paste(c("Summary per method:", txt), collapse = "\n"))
    }
    true <- d$dif_items
    lines <- vapply(rs, function(r) {
      if (inherits(r$res, "error")) {
        sprintf("%-12s : ERROR — %s", r$method, conditionMessage(r$res))
      } else {
        flags <- get_flags(r$res)
        tp <- sum(true %in% flags)
        fp <- sum(!(flags %in% true))
        fn <- sum(!(true %in% flags))
        sprintf("%-12s : %d DIF detected | %d true positives | %d false positives | %d missed",
                r$method, length(flags), tp, fp, fn)
      }
    }, character(1))
    paste(c("Summary per method (vs ground truth):", lines), collapse = "\n")
  })
  
  output$method_plot <- renderPlot({
    rs <- results()
    ok <- Filter(function(r) !inherits(r$res, "error"), rs)
    if (length(ok) == 0) return(invisible(NULL))
    # Graphique de la première méthode disponible (statistiques de DIF)
    plot(ok[[1]]$res, main = paste("DIF statistics —", ok[[1]]$method))
  })
  
  # --- Sorties brutes --------------------------------------------------
  output$raw_output <- renderPrint({
    rs <- results()
    txt <- lapply(rs, function(r) {
      if (inherits(r$res, "error")) {
        paste0("== ", r$method, " ==\nERROR: ", conditionMessage(r$res))
      } else {
        paste0(capture.output(print(r$res)), collapse = "\n")
      }
    })
    cat(paste(txt, collapse = "\n\n----------------------------------------\n\n"))
  })
  
  # --- Documentation ----------------------------------------------------
  # Le contenu vient d'un fichier externe "README.md", placé à côté
  # du app.R. On le relit à chaque ouverture de l'onglet : il suffit de
  # modifier le fichier pour que l'app soit à jour (aucun redémarrage requis).
  # NB : la docu utilise des symboles Unicode (θ, ≥, sₖ...) au lieu de LaTeX
  # -> pas besoin de MathJax, un simple rendu markdown suffit, et ça règle
  # tous les problèmes d'affichage (chevauchements, backslashes perdus...).
  docs_path <- "README.md"
  output$docs <- renderUI({
    if (!file.exists(docs_path)) {
      return(tags$p(style = "color:#B22222;",
                    paste0("File '", docs_path, "' was not found next to app.R ",
                           "— the documentation cannot be displayed.")))
    }
    html <- markdown::markdownToHTML(docs_path, fragment.only = TRUE)
    div(id = "docs-content", HTML(html))
  })
}

shinyApp(ui = ui, server = server)