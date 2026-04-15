clusterModel_new = function (data, cluster.info, admin.info, X = NULL, X.unit = NULL, 
          X.pixel = NULL, admin, CI = 0.95, model = c("bym2", "iid"), 
          stratification = FALSE, aggregation = FALSE, nested = FALSE, 
          interact = FALSE, overdisp.mean = 0, overdisp.prec = 0.4, 
          pc.u = 1, pc.alpha = 0.01, pc.u.phi = 0.5, pc.alpha.phi = 2/3, 
          ...) 
{
  cvf <- function(p, sd) {
    pmax(sd/p, sd/(1 - p))
  }
  if ("age" %in% colnames(data)) {
    fit <- clusterModel_u5mr(data = data, cluster.info = cluster.info, 
                             admin.info = admin.info, X = X, X.unit = X.unit, 
                             X.pixel = X.pixel, admin = admin, CI = CI, model = model, 
                             stratification = stratification, aggregation = aggregation, 
                             nested = nested, interact = interact, overdisp.mean = overdisp.mean, 
                             overdisp.prec = overdisp.prec, pc.u = pc.u, pc.alpha = pc.alpha, 
                             pc.u.phi = pc.u.phi, pc.alpha.phi = pc.alpha.phi, 
                             ...)
    return(fit)
  }
  if (sum(is.na(data$value)) > 0) {
    data <- data[rowSums(is.na(data)) == 0, ]
    message("Removing NAs in indicator response")
  }
  if (!is.null(admin.info)) {
    admin.info.output = admin.info
  }
  else {
    admin.info.output = NULL
  }
  if (is.null(X) == FALSE && dim(X)[1] != dim(admin.info$data)[1]) {
    message("Not valid covariates format. No covariates model is fitted instead")
    X = NULL
  }
  if (is.null(X.unit) == FALSE && is.null(X.pixel) == FALSE) {
    checkset = setequal(colnames(X.pixel)[!colnames(X.pixel) %in% 
                                            c("admin2.name.full", "admin1.name", "Population", 
                                              "strata")], colnames(X.unit)[colnames(X.unit) != 
                                                                             "cluster"])
    if (checkset) {
      X.unit.model = TRUE
    }
    else {
      message("covariates names are not matched in X.unit and X.pixel, non covariate model is fitted instead")
      X.unit.model = FALSE
      X.pixel = NULL
      X.unit = NULL
    }
  }
  else {
    X.unit.model = FALSE
  }
  if (unique(is.na(admin.info$data$urban)) && stratification == 
      T) {
    message("No urban/rural proportion found. Unstratified model is fitted instead")
    stratification = FALSE
  }
  if (unique(is.na(admin.info$data$population)) && aggregation == 
      T) {
    message("No population found")
    aggregation = FALSE
  }
  if (model == "bym2") {
    Amat <- admin.info$mat
  }
  else {
    Amat <- NULL
  }
  admin.mat <- Amat
  admin.info <- admin.info$data
  if (!is.null(cluster.info)) {
    modt <- dplyr::left_join(data, cluster.info$data, by = "cluster")
  } else {
    modt <- data
  }
  if (model == "bym2") {
    if (!("LONGNUM" %in% names(modt))) {
      stop("model = 'bym2' requires LONGNUM in data or cluster.info.")
    }
    modt <- modt[!is.na(modt$LONGNUM), ]
  }
  modt$strata.full <- paste(modt$admin1.name, modt$strata)
  if (nested & admin > 1) {
  }
  else if (nested & admin == 1) {
    message("Nested model is designed for Admin 2 or finer level. An Admin 1 model will be fitted")
    nested = FALSE
  }
  c.dat.tmp <- modt %>% group_by(cluster) %>% mutate(n = length(cluster)) %>% 
    mutate(value = sum(value, na.rm = T)) %>% ungroup() %>% 
    distinct(cluster, .keep_all = TRUE)
  admin_name_table <- admin.info
  nregion <- dim(admin_name_table)[1]
  admin_name_table$sID <- 1:dim(admin_name_table)[1]
  if (admin == 1) {
    c.dat.tmp[(dim(c.dat.tmp)[1] + 1):(dim(c.dat.tmp)[1] + 
                                         nregion), paste0("admin", admin, ".name")] <- admin_name_table[, 
                                                                                                        which(colnames(admin_name_table) == paste0("admin", 
                                                                                                                                                   admin, ".name"))]
    c.dat.tmp$ID <- 1:dim(c.dat.tmp)[1]
    c.dat.tmp$sID <- admin_name_table$sID[match(as.data.frame(c.dat.tmp)[, 
                                                                         which(colnames(c.dat.tmp) == paste0("admin", admin, 
                                                                                                             ".name"))], admin_name_table[, which(colnames(admin_name_table) == 
                                                                                                                                                    paste0("admin", admin, ".name"))])]
    if (is.null(X) == FALSE) {
      c.dat.tmp <- left_join(c.dat.tmp, X, by = "admin1.name")
    }
    if (is.null(X.unit) == FALSE) {
      c.dat.tmp <- left_join(c.dat.tmp, X.unit, by = "cluster")
    }
  }
  else {
    c.dat.tmp[(dim(c.dat.tmp)[1] + 1):(dim(c.dat.tmp)[1] + 
                                         nregion), "admin2.name.full"] <- admin_name_table[, 
                                                                                           which(colnames(admin_name_table) == "admin2.name.full")]
    c.dat.tmp$ID <- 1:dim(c.dat.tmp)[1]
    c.dat.tmp$sID <- admin_name_table$sID[match(as.data.frame(c.dat.tmp)[, 
                                                                         which(colnames(c.dat.tmp) == "admin2.name.full")], 
                                                admin_name_table[, which(colnames(admin_name_table) == 
                                                                           "admin2.name.full")])]
    if (nested) {
      c.dat.tmp[(dim(c.dat.tmp)[1] - (nregion - 1)):(dim(c.dat.tmp)[1]), 
                "admin1.name"] <- admin_name_table[, which(colnames(admin_name_table) == 
                                                             "admin1.name")]
    }
    if (is.null(X) == FALSE) {
      c.dat.tmp <- left_join(c.dat.tmp, X, by = "admin2.name.full")
    }
    if (is.null(X.unit) == FALSE) {
      c.dat.tmp <- left_join(c.dat.tmp, X.unit, by = "cluster")
    }
  }
  if (stratification == F) {
    if (model == "iid") {
      if (is.null(X) == FALSE) {
        cvrt <- paste(" + ", colnames(X)[2:dim(X)[2]], 
                      collapse = "")
      }
      else if (is.null(X.unit) == FALSE) {
        cvrt <- paste(" + ", colnames(X.unit)[2:dim(X.unit)[2]], 
                      collapse = "")
      }
      else {
        cvrt = NULL
      }
      if (nested) {
        sptl <- "+ f(sID, model = model, graph = admin.mat, hyper = list(prec = list(prior = \"pc.prec\", param = c(pc.u , pc.alpha))))"
        formula_string <- paste("value ~  -1+ admin1.name", 
                                cvrt, sptl)
      }
      else {
        sptl <- "+ f(sID, model = model, graph = admin.mat, hyper = list(prec = list(prior = \"pc.prec\", param = c(pc.u , pc.alpha))))"
        formula_string <- paste("value ~ 1", cvrt, sptl)
      }
      formula <- as.formula(formula_string)
    }
    else if (model == "bym2") {
      if (is.null(X) == FALSE) {
        cvrt <- paste(" + ", colnames(X)[2:dim(X)[2]], 
                      collapse = "")
      }
      else if (is.null(X.unit) == FALSE) {
        cvrt <- paste(" + ", colnames(X.unit)[2:dim(X.unit)[2]], 
                      collapse = "")
      }
      else {
        cvrt = NULL
      }
      if (nested) {
        sptl <- "+ f(sID, model = model, graph = admin.mat,hyper = list( prec = list(prior = \"pc.prec\", param = c(pc.u , pc.alpha)), phi = list(prior = 'pc', param = c(pc.u.phi , pc.alpha.phi))))"
        formula_string <- paste("value ~ -1+ admin1.name", 
                                cvrt, sptl)
      }
      else {
        sptl <- "+ f(sID, model = model, graph = admin.mat,  hyper = list( prec = list(prior = \"pc.prec\", param = c(pc.u , pc.alpha)), phi = list(prior = 'pc', param = c(pc.u.phi , pc.alpha.phi))))"
        formula_string <- paste("value ~ 1", cvrt, sptl)
      }
      formula <- as.formula(formula_string)
    }
  }
  else if (stratification) {
    c.dat.tmp$strata = as.factor(c.dat.tmp$strata)
    if (model == "iid") {
      if (is.null(X) == FALSE) {
        cvrt <- paste(" + ", colnames(X)[2:dim(X)[2]], 
                      collapse = "")
      }
      else if (is.null(X.unit) == FALSE) {
        cvrt <- paste(" + ", colnames(X.unit)[2:dim(X.unit)[2]], 
                      collapse = "")
      }
      else {
        cvrt = NULL
      }
      sptl <- "+ f(sID, model = model, graph = admin.mat, hyper = list(prec = list(prior = \"pc.prec\", param = c(pc.u , pc.alpha))))"
      formula_string <- paste("value ~ 1 + strata", cvrt, 
                              sptl)
      formula <- as.formula(formula_string)
      if (nested) {
        if (interact) {
          sptl <- "+ f(sID, model = model, graph = admin.mat, hyper = list(prec = list(prior = \"pc.prec\", param = c(pc.u , pc.alpha))))"
          formula_string <- paste("value ~ -1+admin1.name:strata", 
                                  cvrt, sptl)
        }
        else {
          sptl <- "+ f(sID, model = model, graph = admin.mat, hyper = list(prec = list(prior = \"pc.prec\", param = c(pc.u , pc.alpha))))"
          formula_string <- paste("value ~ -1+ admin1.name + strata", 
                                  cvrt, sptl)
        }
      }
      else {
        sptl <- "+ f(sID, model = model, graph = admin.mat, hyper = list(prec = list(prior = \"pc.prec\", param = c(pc.u , pc.alpha))))"
        formula_string <- paste("value ~ 1 + strata", 
                                cvrt, sptl)
      }
      formula <- as.formula(formula_string)
    }
    else if (model == "bym2") {
      if (is.null(X) == FALSE) {
        cvrt <- paste(" + ", colnames(X)[2:dim(X)[2]], 
                      collapse = "")
      }
      else if (is.null(X.unit) == FALSE) {
        cvrt <- paste(" + ", colnames(X.unit)[2:dim(X.unit)[2]], 
                      collapse = "")
      }
      else {
        cvrt = NULL
      }
      if (nested) {
        if (interact) {
          sptl <- "+ f(sID, model = model, graph = admin.mat, hyper = list( prec = list(prior = \"pc.prec\", param = c(pc.u , pc.alpha)), phi = list(prior = 'pc', param = c(pc.u.phi , pc.alpha.phi))))"
          formula_string <- paste("value ~ -1+admin1.name:strata", 
                                  cvrt, sptl)
        }
        else {
          sptl <- "+ f(sID, model = model, graph = admin.mat, hyper = list( prec = list(prior = \"pc.prec\", param = c(pc.u , pc.alpha)), phi = list(prior = 'pc', param = c(pc.u.phi , pc.alpha.phi))))"
          formula_string <- paste("value ~ -1+ admin1.name + strata", 
                                  cvrt, sptl)
        }
      }
      else {
        sptl <- "+ f(sID, model = model, graph = admin.mat,  hyper = list( prec = list(prior = \"pc.prec\", param = c(pc.u , pc.alpha)), phi = list(prior = 'pc', param = c(pc.u.phi , pc.alpha.phi))))"
        formula_string <- paste("value ~ 1+ strata", 
                                cvrt, sptl)
      }
      formula <- as.formula(formula_string)
    }
  }
  if (stratification && sum(!unique(c.dat.tmp$strata) %in% 
                            c("urban", "rural", NA))) {
    stop("The variable strata in the input data can only be 'urban' or 'rural' currently.")
  }
  overdisp.mean = overdisp.mean
  overdisp.prec = overdisp.prec
  control.family <- list(hyper = list(rho = list(param = c(overdisp.mean, 
                                                           overdisp.prec), initial = overdisp.mean)))
  imod <- INLA::inla(formula, family = "betabinomial", data = c.dat.tmp, 
                     Ntrials = n, control.predictor = list(compute = TRUE, 
                                                           link = 1), control.compute = list(config = TRUE, 
                                                                                             waic = TRUE), control.family = control.family)
  nsamp <- 1000
  samp <- INLA::inla.posterior.sample(n = nsamp, result = imod, 
                                      intern = TRUE)
  if (X.unit.model == FALSE) {
    draw.u <- matrix(NA, nsamp, nregion)
    draw.r <- matrix(NA, nsamp, nregion)
    draw.all <- matrix(NA, nsamp, nregion)
    for (i in 1:length(samp)) {
      tmp <- samp[[i]]$latent
      s.effect <- data.frame(s.effect = tmp[paste0("sID:", 
                                                   1:nregion), 1])
      s.effect$sID <- as.integer(sub("sID:(.*)", "\\1", 
                                     rownames(s.effect)))
      l.com = left_join(admin_name_table, s.effect, by = "sID")
      if (nested) {
        NN <- admin.info %>% group_by(admin1.name) %>% 
          mutate(N = length(admin2.name.full)) %>% select(c("admin1.name", 
                                                            "N")) %>% arrange(admin1.name) %>% distinct(admin1.name, 
                                                                                                        .keep_all = TRUE)
        n <- dim(NN)[1]
        if (interact) {
          rows_to_extract <- paste0("admin1.name", NN$admin1.name, 
                                    ":stratarural:1")
          intercept.fix.rural <- as.data.frame(tmp[rownames(tmp) %in% 
                                                     rows_to_extract, , drop = FALSE])
          intercept.fix.rural$admin1.name <- sub("admin1.name(.*):stratarural:1", 
                                                 "\\1", rownames(intercept.fix.rural))
          colnames(intercept.fix.rural)[colnames(intercept.fix.rural) == 
                                          "V1"] <- "intercept.rural"
          rows_to_extract <- paste0("admin1.name", NN$admin1.name, 
                                    ":strataurban:1")
          intercept.fix.urban <- as.data.frame(tmp[rownames(tmp) %in% 
                                                     rows_to_extract, , drop = FALSE])
          intercept.fix.urban$admin1.name <- sub("admin1.name(.*):strataurban:1", 
                                                 "\\1", rownames(intercept.fix.urban))
          colnames(intercept.fix.urban)[colnames(intercept.fix.urban) == 
                                          "V1"] <- "intercept.urban"
          l.com = left_join(l.com, intercept.fix.rural, 
                            by = "admin1.name")
          l.com = left_join(l.com, intercept.fix.urban, 
                            by = "admin1.name")
        }
        else {
          rows_to_extract <- paste0("admin1.name", NN$admin1.name, 
                                    ":1")
          intercept.fix <- as.data.frame(tmp[rows_to_extract, 
                                             , drop = FALSE])
          intercept.fix$admin1.name <- sub("admin1.name(.*):1", 
                                           "\\1", rownames(intercept.fix))
          colnames(intercept.fix)[colnames(intercept.fix) == 
                                    "V1"] <- "intercept"
          l.com = left_join(l.com, intercept.fix, by = "admin1.name")
        }
      }
      else {
        l.com$intercept <- tmp["(Intercept):1", 1]
      }
      if (is.null(X) == FALSE) {
        l.com$covariates = colSums(tmp[paste0(colnames(X)[!colnames(X) %in% 
                                                            c("admin2.name.full", "admin1.name")], ":1"), 
                                       1] %*% t(as.matrix(X[, colnames(X)[!colnames(X) %in% 
                                                                            c("admin2.name.full", "admin1.name")]])))
      }
      else {
        l.com$covariates = 0
      }
      if (stratification) {
        if (interact) {
        }
        else {
          if ("stratarural:1" %in% rownames(tmp)) {
            str.effect <- tmp["stratarural:1", 1]
            str.effect.u <- 0
          }
          else {
            str.effect <- 0
            str.effect.u <- tmp["strataurban:1", 1]
          }
        }
      }
      if (stratification == FALSE) {
        draw.all[i, ] <- SUMMER::expit(l.com$s.effect + 
                                         l.com$intercept + l.com$covariates)
      }
      else if (stratification) {
        if (interact) {
          draw.u[i, ] <- SUMMER::expit(l.com$s.effect + 
                                         l.com$intercept.urban + l.com$covariates)
          draw.r[i, ] <- SUMMER::expit(l.com$s.effect + 
                                         l.com$intercept.rural + l.com$covariates)
          draw.all[i, ] <- draw.u[i, ] * admin.info$urban + 
            draw.r[i, ] * (1 - admin.info$urban)
        }
        else {
          draw.u[i, ] <- SUMMER::expit(l.com$s.effect + 
                                         l.com$intercept + str.effect.u + l.com$covariates)
          draw.r[i, ] <- SUMMER::expit(l.com$s.effect + 
                                         l.com$intercept + str.effect + l.com$covariates)
          draw.all[i, ] <- draw.u[i, ] * admin.info$urban + 
            draw.r[i, ] * (1 - admin.info$urban)
        }
      }
    }
    if (admin == 1) {
      if (stratification == FALSE) {
        if (is.null(X.pixel) == FALSE) {
          post.all <- apply(draw.all, 2, mean)
          post.all.sd <- apply(draw.all, 2, sd)
          post.all.median <- apply(draw.all, 2, median)
          post.all.ci <- apply(draw.all, 2, quantile, 
                               probs = c((1 - CI)/2, 1 - (1 - CI)/2))
          admin1.bb.res <- data.frame(admin1.name = admin.info$admin1.name, 
                                      mean = post.all, median = post.all.median, 
                                      sd = post.all.sd, var = post.all.sd^2, lower = post.all.ci[1, 
                                      ], upper = post.all.ci[2, ], cv = cvf(post.all.median, 
                                                                            post.all.sd))
        }
        else {
          admin1.bb.res <- data.frame(cbind(admin1.name = admin.info$admin1.name[tail(c.dat.tmp$sID, 
                                                                                      n = nregion)], mean = tail(imod$summary.fitted.values, 
                                                                                                                 n = nregion)[, 1], median = colMedians(draw.all), 
                                            sd = tail(imod$summary.fitted.values, n = nregion)[, 
                                                                                               2], var = tail(imod$summary.fitted.values, 
                                                                                                              n = nregion)[, 2]^2, lower = tail(imod$summary.fitted.values, 
                                                                                                                                                n = nregion)[, 3], upper = tail(imod$summary.fitted.values, 
                                                                                                                                                                                n = nregion)[, 5]))
          admin1.bb.res$mean <- as.numeric(admin1.bb.res$mean)
          admin1.bb.res$median <- as.numeric(admin1.bb.res$median)
          admin1.bb.res$sd <- as.numeric(admin1.bb.res$sd)
          admin1.bb.res$var <- admin1.bb.res$sd^2
          admin1.bb.res$lower <- as.numeric(admin1.bb.res$lower)
          admin1.bb.res$upper <- as.numeric(admin1.bb.res$upper)
          admin1.bb.res$cv = cvf(admin1.bb.res$median, 
                                 admin1.bb.res$sd)
        }
      }
      else if (stratification) {
        post.u <- apply(draw.u, 2, mean)
        post.r <- apply(draw.r, 2, mean)
        post.all <- apply(draw.all, 2, mean)
        post.u.sd <- apply(draw.u, 2, sd)
        post.r.sd <- apply(draw.r, 2, sd)
        post.all.sd <- apply(draw.all, 2, sd)
        post.u.median <- apply(draw.u, 2, median)
        post.r.median <- apply(draw.r, 2, median)
        post.all.median <- apply(draw.all, 2, median)
        post.u.ci <- apply(draw.u, 2, quantile, probs = c((1 - 
                                                             CI)/2, 1 - (1 - CI)/2))
        post.r.ci <- apply(draw.r, 2, quantile, probs = c((1 - 
                                                             CI)/2, 1 - (1 - CI)/2))
        post.all.ci <- apply(draw.all, 2, quantile, 
                             probs = c((1 - CI)/2, 1 - (1 - CI)/2))
        admin1.bb.res <- data.frame(admin1.name = rep(admin.info$admin1.name, 
                                                      3), mean = c(post.u, post.r, post.all), median = c(post.u.median, 
                                                                                                         post.r.median, post.all.median), sd = c(post.u.sd, 
                                                                                                                                                 post.r.sd, post.all.sd), var = c(post.u.sd^2, 
                                                                                                                                                                                  post.r.sd^2, post.all.sd^2), lower = c(post.u.ci[1, 
                                                                                                                                                                                  ], post.r.ci[1, ], post.all.ci[1, ]), upper = c(post.u.ci[2, 
                                                                                                                                                                                  ], post.r.ci[2, ], post.all.ci[2, ]), cv = c(cvf(post.u.median, 
                                                                                                                                                                                                                                   post.u.sd), cvf(post.r.median, post.r.sd), 
                                                                                                                                                                                                                               cvf(post.all.median, post.all.sd)), type = c(rep("urban", 
                                                                                                                                                                                                                                                                                nregion), rep("rural", nregion), rep("full", 
                                                                                                                                                                                                                                                                                                                     nregion)))
      }
      if (aggregation == T) {
        post.all <- draw.all %*% admin.info$population/sum(admin.info$population)
        agg.natl <- data.frame(mean = mean(post.all), 
                               median = median(post.all), sd = sd(post.all), 
                               var = var(post.all), lower = quantile(post.all, 
                                                                     probs = c((1 - CI)/2, 1 - (1 - CI)/2))[1], 
                               upper = quantile(post.all, probs = c((1 - 
                                                                       CI)/2, 1 - (1 - CI)/2))[2])
        agg.natl$cv = cvf(agg.natl$median, agg.natl$sd)
        rownames(agg.natl) = NULL
      }
      if (aggregation == FALSE) {
        cm = list(res.admin1 = admin1.bb.res, inla = imod, 
                  admin1_post = draw.all, urban_post = draw.u, 
                  rural_post = draw.r, admin.info = admin.info.output, 
                  admin = admin)
        attr(cm, "class") = "clusterModel"
        attr(cm, "domain.names") <- admin.info$admin1.name
        return(cm)
      }
      else {
        cm = list(res.admin1 = admin1.bb.res, agg.natl = agg.natl, 
                  inla = imod, admin1_post = draw.all, nation_post = post.all, 
                  urban_post = draw.u, rural_post = draw.r, 
                  admin.info = admin.info, admin = admin)
        attr(cm, "class") = "clusterModel"
        attr(cm, "domain.names") <- admin.info$admin1.name
        return(cm)
      }
    }
    else if (admin == 2) {
      if (stratification == FALSE) {
        admin2.bb.res <- data.frame(cbind(admin2.name.full = admin.info$admin2.name.full[tail(c.dat.tmp$sID, 
                                                                                              n = nregion)], mean = tail(imod$summary.fitted.values, 
                                                                                                                         n = nregion)[, 1], median = colMedians(draw.all), 
                                          sd = tail(imod$summary.fitted.values, n = nregion)[, 
                                                                                             2], var = tail(imod$summary.fitted.values, 
                                                                                                            n = nregion)[, 2]^2, lower = tail(imod$summary.fitted.values, 
                                                                                                                                              n = nregion)[, 3], upper = tail(imod$summary.fitted.values, 
                                                                                                                                                                              n = nregion)[, 5]))
        admin2.bb.res$mean <- as.numeric(admin2.bb.res$mean)
        admin2.bb.res$median <- as.numeric(admin2.bb.res$median)
        admin2.bb.res$sd <- as.numeric(admin2.bb.res$sd)
        admin2.bb.res$var <- admin2.bb.res$sd^2
        admin2.bb.res$lower <- as.numeric(admin2.bb.res$lower)
        admin2.bb.res$upper <- as.numeric(admin2.bb.res$upper)
        admin2.bb.res$cv = cvf(admin2.bb.res$median, 
                               admin2.bb.res$sd)
        admin2.bb.res <- left_join(admin2.bb.res, distinct(admin.info), 
                                   by = "admin2.name.full")
      }
      else if (stratification) {
        post.u <- apply(draw.u, 2, mean)
        post.r <- apply(draw.r, 2, mean)
        post.all <- apply(draw.all, 2, mean)
        post.u.sd <- apply(draw.u, 2, sd)
        post.r.sd <- apply(draw.r, 2, sd)
        post.all.sd <- apply(draw.all, 2, sd)
        post.u.median <- apply(draw.u, 2, median)
        post.r.median <- apply(draw.r, 2, median)
        post.all.median <- apply(draw.all, 2, median)
        post.u.ci <- apply(draw.u, 2, quantile, probs = c((1 - 
                                                             CI)/2, 1 - (1 - CI)/2))
        post.r.ci <- apply(draw.r, 2, quantile, probs = c((1 - 
                                                             CI)/2, 1 - (1 - CI)/2))
        post.all.ci <- apply(draw.all, 2, quantile, 
                             probs = c((1 - CI)/2, 1 - (1 - CI)/2))
        admin2.bb.res <- data.frame(admin2.name.full = rep(admin.info$admin2.name.full, 
                                                           3), mean = c(post.u, post.r, post.all), median = c(post.u.median, 
                                                                                                              post.r.median, post.all.median), sd = c(post.u.sd, 
                                                                                                                                                      post.r.sd, post.all.sd), var = c(post.u.sd^2, 
                                                                                                                                                                                       post.r.sd^2, post.all.sd^2), lower = c(post.u.ci[1, 
                                                                                                                                                                                       ], post.r.ci[1, ], post.all.ci[1, ]), upper = c(post.u.ci[2, 
                                                                                                                                                                                       ], post.r.ci[2, ], post.all.ci[2, ]), cv = c(cvf(post.u.median, 
                                                                                                                                                                                                                                        post.u.sd), cvf(post.r.median, post.r.sd), 
                                                                                                                                                                                                                                    cvf(post.all.median, post.all.sd)), type = c(rep("urban", 
                                                                                                                                                                                                                                                                                     nregion), rep("rural", nregion), rep("full", 
                                                                                                                                                                                                                                                                                                                          nregion)))
        admin2.bb.res <- left_join(admin2.bb.res, distinct(admin.info), 
                                   by = "admin2.name.full")
      }
      if (aggregation) {
        weight = admin.info$population/admin.info$population.admin1
        post.all <- data.table::data.table(t(weight * 
                                               t(draw.all)))
        colnames(post.all) <- admin.info$admin1.name
        subgroups <- split.default(post.all, names(post.all))
        sums_list <- lapply(subgroups, function(subgroup) {
          rowSums(subgroup)
        })
        admin1.samp <- do.call(cbind, sums_list)
        agg.admin1 <- data.frame(mean = colMeans(admin1.samp), 
                                 median = colMedians(admin1.samp), sd = apply(admin1.samp, 
                                                                              2, sd), var = apply(admin1.samp, 2, var), 
                                 lower = apply(admin1.samp, 2, quantile, probs = c((1 - 
                                                                                      CI)/2, 1 - (1 - CI)/2))[1, ], upper = apply(admin1.samp, 
                                                                                                                                  2, quantile, probs = c((1 - CI)/2, 1 - (1 - 
                                                                                                                                                                            CI)/2))[2, ], cv = cvf(colMedians(admin1.samp), 
                                                                                                                                                                                                   apply(admin1.samp, 2, sd)))
        agg.admin1$admin1.name = rownames(agg.admin1)
        rownames(agg.admin1) = NULL
        agg.admin1 <- agg.admin1 %>% select("admin1.name", 
                                            "mean", "median", "sd", "var", "lower", "upper", 
                                            "cv")
        unique(admin.info$population.admin1)/sum(unique(admin.info$population.admin1))
        post.all <- admin1.samp %*% unique(admin.info$population.admin1)/sum(unique(admin.info$population.admin1))
        agg.natl <- data.frame(mean = mean(post.all), 
                               median = median(post.all), sd = sd(post.all), 
                               var = var(post.all), lower = quantile(post.all, 
                                                                     probs = c((1 - CI)/2, 1 - (1 - CI)/2))[1], 
                               upper = quantile(post.all, probs = c((1 - 
                                                                       CI)/2, 1 - (1 - CI)/2))[2])
        agg.natl$cv = cvf(agg.natl$median, agg.natl$sd)
        rownames(agg.natl) = NULL
      }
      if (aggregation == FALSE) {
        cm = list(res.admin2 = admin2.bb.res, inla = imod, 
                  admin2_post = draw.all, urban_post = draw.u, 
                  rural_post = draw.r, admin.info = admin.info.output, 
                  admin = admin)
        attr(cm, "class") = "clusterModel"
        attr(cm, "domain.names") <- admin.info$admin2.name.full
        return(cm)
      }
      else {
        cm = list(res.admin2 = admin2.bb.res, agg.admin1 = agg.admin1, 
                  agg.natl = agg.natl, inla = imod, admin2_post = draw.all, 
                  admin1_post = admin1.samp, nation_post = post.all, 
                  urban_post = draw.u, rural_post = draw.r, 
                  admin.info = admin.info.output, admin = admin)
        attr(cm, "class") = "clusterModel"
        attr(cm, "domain.names") <- admin.info$admin2.name.full
        return(cm)
      }
    }
  }
  else if (X.unit.model) {
    draw.u <- matrix(NA, nsamp, length(unique(X.pixel$admin2.name.full)))
    draw.r <- matrix(NA, nsamp, length(unique(X.pixel$admin2.name.full)))
    draw.all <- matrix(NA, nsamp, length(unique(X.pixel$admin2.name.full)))
    draw.u1 <- matrix(NA, nsamp, length(unique(admin.info$admin1.name)))
    draw.r1 <- matrix(NA, nsamp, length(unique(admin.info$admin1.name)))
    draw.all1 <- matrix(NA, nsamp, length(unique(admin.info$admin1.name)))
    draw.u0 <- matrix(NA, nsamp, 1)
    draw.r0 <- matrix(NA, nsamp, 1)
    draw.all0 <- matrix(NA, nsamp, 1)
    for (i in 1:length(samp)) {
      covpixel = X.pixel
      tmp <- samp[[i]]$latent
      s.effect <- data.frame(s.effect = tmp[paste0("sID:", 
                                                   1:nregion), 1])
      s.effect$sID <- as.integer(sub("sID:(.*)", "\\1", 
                                     rownames(s.effect)))
      l.com = left_join(admin_name_table, s.effect, by = "sID")
      if (nested) {
        NN <- admin.info %>% group_by(admin1.name) %>% 
          mutate(N = length(admin2.name.full)) %>% select(c("admin1.name", 
                                                            "N")) %>% arrange(admin1.name) %>% distinct(admin1.name, 
                                                                                                        .keep_all = TRUE)
        n <- dim(NN)[1]
        if (interact) {
          rows_to_extract <- paste0("admin1.name", NN$admin1.name, 
                                    ":stratarural:1")
          intercept.fix.rural <- as.data.frame(tmp[rownames(tmp) %in% 
                                                     rows_to_extract, , drop = FALSE])
          intercept.fix.rural$admin1.name <- sub("admin1.name(.*):stratarural:1", 
                                                 "\\1", rownames(intercept.fix.rural))
          colnames(intercept.fix.rural)[colnames(intercept.fix.rural) == 
                                          "V1"] <- "intercept.rural"
          rows_to_extract <- paste0("admin1.name", NN$admin1.name, 
                                    ":strataurban:1")
          intercept.fix.urban <- as.data.frame(tmp[rownames(tmp) %in% 
                                                     rows_to_extract, , drop = FALSE])
          intercept.fix.urban$admin1.name <- sub("admin1.name(.*):strataurban:1", 
                                                 "\\1", rownames(intercept.fix.urban))
          colnames(intercept.fix.urban)[colnames(intercept.fix.urban) == 
                                          "V1"] <- "intercept.urban"
          l.com = left_join(l.com, intercept.fix.rural, 
                            by = "admin1.name")
          l.com = left_join(l.com, intercept.fix.urban, 
                            by = "admin1.name")
          l.com$intercept.rural = ifelse(is.na(l.com$intercept.rural), 
                                         0, l.com$intercept.rural)
          l.com$intercept.urban = ifelse(is.na(l.com$intercept.urban), 
                                         0, l.com$intercept.urban)
        }
        else {
          rows_to_extract <- paste0("admin1.name", NN$admin1.name, 
                                    ":1")
          intercept.fix <- as.data.frame(tmp[rows_to_extract, 
                                             , drop = FALSE])
          intercept.fix$admin1.name <- sub("admin1.name(.*):1", 
                                           "\\1", rownames(intercept.fix))
          colnames(intercept.fix)[colnames(intercept.fix) == 
                                    "V1"] <- "intercept"
          l.com = left_join(l.com, intercept.fix, by = "admin1.name")
        }
      }
      else {
        l.com$intercept <- tmp["(Intercept):1", 1]
      }
      if (stratification) {
        if (interact) {
          l.com$str.effect = l.com$intercept.rural
          l.com$str.effect.u = l.com$intercept.urban
          l.com$intercept = 0
        }
        else {
          if ("stratarural:1" %in% rownames(tmp)) {
            l.com$str.effect <- tmp["stratarural:1", 
                                    1]
            l.com$str.effect.u <- 0
          }
          else {
            l.com$str.effect <- 0
            l.com$str.effect.u <- tmp["strataurban:1", 
                                      1]
          }
        }
      }
      else {
        l.com$str.effect <- 0
        l.com$str.effect.u <- 0
      }
      if (admin == 2) {
        covpixel = left_join(covpixel, l.com[, c("admin2.name.full", 
                                                 "intercept", "s.effect", "str.effect.u", "str.effect")], 
                             by = "admin2.name.full")
        covpixel$covariates = colSums(tmp[paste0(colnames(X.unit)[colnames(X.unit) != 
                                                                    "cluster"], ":1"), 1] %*% t(as.matrix(covpixel[, 
                                                                                                                   colnames(X.unit)[colnames(X.unit) != "cluster"]])))
      }
      else {
        covpixel = left_join(covpixel, l.com[, c("admin1.name", 
                                                 "intercept", "s.effect", "str.effect.u", "str.effect")], 
                             by = "admin1.name")
        covpixel$covariates = colSums(tmp[paste0(colnames(X.unit)[colnames(X.unit) != 
                                                                    "cluster"], ":1"), 1] %*% t(as.matrix(covpixel[, 
                                                                                                                   colnames(X.unit)[colnames(X.unit) != "cluster"]])))
      }
      if (stratification) {
        covpixel$p.u = SUMMER::expit(covpixel$s.effect + 
                                       covpixel$intercept + covpixel$covariates + 
                                       covpixel$str.effect.u * covpixel$strata)
        covpixel$p.r = SUMMER::expit(covpixel$s.effect + 
                                       covpixel$intercept + covpixel$covariates + 
                                       covpixel$str.effect * (1 - covpixel$strata))
        covpixel$p = SUMMER::expit(covpixel$s.effect + 
                                     covpixel$intercept + covpixel$covariates + 
                                     covpixel$str.effect * (1 - covpixel$strata) + 
                                     covpixel$str.effect.u * covpixel$strata)
        unique_areas <- sort(unique(covpixel$admin2.name.full))
        draw.u[i, ] = covpixel %>% filter(strata == 
                                            1) %>% group_by(admin2.name.full) %>% summarise(weighted_value = sum(p.u * 
                                                                                                                   Population)/sum(Population), .groups = "drop") %>% 
          right_join(data.frame(admin2.name.full = unique_areas), 
                     by = "admin2.name.full") %>% arrange(admin2.name.full) %>% 
          pull(weighted_value)
        draw.r[i, ] = covpixel %>% filter(strata == 
                                            0) %>% group_by(admin2.name.full) %>% summarise(weighted_value = sum(p.r * 
                                                                                                                   Population)/sum(Population), .groups = "drop") %>% 
          right_join(data.frame(admin2.name.full = unique_areas), 
                     by = "admin2.name.full") %>% arrange(admin2.name.full) %>% 
          pull(weighted_value)
        draw.all[i, ] = covpixel %>% group_by(admin2.name.full) %>% 
          summarise(weighted_value = sum(p * Population)/sum(Population), 
                    .groups = "drop") %>% right_join(data.frame(admin2.name.full = unique_areas), 
                                                     by = "admin2.name.full") %>% pull(weighted_value)
        unique_areas <- sort(unique(covpixel$admin1.name))
        draw.u1[i, ] = covpixel %>% filter(strata == 
                                             1) %>% group_by(admin1.name) %>% summarise(weighted_value = sum(p.u * 
                                                                                                               Population)/sum(Population), .groups = "drop") %>% 
          right_join(data.frame(admin1.name = unique_areas), 
                     by = "admin1.name") %>% arrange(admin1.name) %>% 
          mutate(weighted_value = ifelse(is.na(weighted_value), 
                                         0, weighted_value))
        pull(weighted_value)
        draw.r1[i, ] = covpixel %>% filter(strata == 
                                             0) %>% group_by(admin1.name) %>% summarise(weighted_value = sum(p.r * 
                                                                                                               Population)/sum(Population), .groups = "drop") %>% 
          right_join(data.frame(admin1.name = unique_areas), 
                     by = "admin1.name") %>% arrange(admin1.name) %>% 
          mutate(weighted_value = ifelse(is.na(weighted_value), 
                                         0, weighted_value))
        pull(weighted_value)
        draw.all1[i, ] = covpixel %>% group_by(admin1.name) %>% 
          summarise(weighted_value = sum(p * Population)/sum(Population), 
                    .groups = "drop") %>% right_join(data.frame(admin1.name = unique_areas), 
                                                     by = "admin1.name") %>% pull(weighted_value)
        draw.u0[i, ] = covpixel %>% filter(strata == 
                                             1) %>% summarise(weighted_value = sum(p.u * 
                                                                                     Population)/sum(Population), .groups = "drop") %>% 
          pull(weighted_value)
        draw.r0[i, ] = covpixel %>% filter(strata == 
                                             0) %>% summarise(weighted_value = sum(p.r * 
                                                                                     Population)/sum(Population), .groups = "drop") %>% 
          pull(weighted_value)
        draw.all0[i, ] = covpixel %>% summarise(weighted_value = sum(p * 
                                                                       Population)/sum(Population), .groups = "drop") %>% 
          pull(weighted_value)
      }
      else {
        covpixel$p = SUMMER::expit(covpixel$s.effect + 
                                     covpixel$intercept + covpixel$covariates)
        unique_areas <- sort(unique(covpixel$admin2.name.full))
        draw.all[i, ] = covpixel %>% group_by(admin2.name.full) %>% 
          summarise(weighted_value = sum(p * Population)/sum(Population), 
                    .groups = "drop") %>% right_join(data.frame(admin2.name.full = unique_areas), 
                                                     by = "admin2.name.full") %>% pull(weighted_value)
        unique_areas <- sort(unique(covpixel$admin1.name))
        draw.all1[i, ] = covpixel %>% group_by(admin1.name) %>% 
          summarise(weighted_value = sum(p * Population)/sum(Population), 
                    .groups = "drop") %>% right_join(data.frame(admin1.name = unique_areas), 
                                                     by = "admin1.name") %>% pull(weighted_value)
        draw.all0[i, ] = covpixel %>% summarise(weighted_value = sum(p * 
                                                                       Population)/sum(Population), .groups = "drop") %>% 
          pull(weighted_value)
      }
    }
    if (stratification) {
      post.u <- apply(draw.u, 2, mean, na.rm = TRUE)
      post.r <- apply(draw.r, 2, mean, na.rm = TRUE)
      post.all <- apply(draw.all, 2, mean)
      post.u.sd <- apply(draw.u, 2, sd, na.rm = TRUE)
      post.r.sd <- apply(draw.r, 2, sd, na.rm = TRUE)
      post.all.sd <- apply(draw.all, 2, sd)
      post.u.median <- apply(draw.u, 2, median, na.rm = TRUE)
      post.r.median <- apply(draw.r, 2, median, na.rm = TRUE)
      post.all.median <- apply(draw.all, 2, median)
      post.u.ci <- apply(draw.u, 2, quantile, probs = c((1 - 
                                                           CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
      post.r.ci <- apply(draw.r, 2, quantile, probs = c((1 - 
                                                           CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
      post.all.ci <- apply(draw.all, 2, quantile, probs = c((1 - 
                                                               CI)/2, 1 - (1 - CI)/2))
      admin2.bb.res <- data.frame(admin2.name.full = rep(sort(unique(covpixel$admin2.name.full)), 
                                                         3), mean = c(post.u, post.r, post.all), median = c(post.u.median, 
                                                                                                            post.r.median, post.all.median), sd = c(post.u.sd, 
                                                                                                                                                    post.r.sd, post.all.sd), var = c(post.u.sd^2, 
                                                                                                                                                                                     post.r.sd^2, post.all.sd^2), lower = c(post.u.ci[1, 
                                                                                                                                                                                     ], post.r.ci[1, ], post.all.ci[1, ]), upper = c(post.u.ci[2, 
                                                                                                                                                                                     ], post.r.ci[2, ], post.all.ci[2, ]), cv = c(cvf(post.u.median, 
                                                                                                                                                                                                                                      post.u.sd), cvf(post.r.median, post.r.sd), cvf(post.all.median, 
                                                                                                                                                                                                                                                                                     post.all.sd)), type = c(rep("urban", length(unique(X.pixel$admin2.name.full))), 
                                                                                                                                                                                                                                                                                                             rep("rural", length(unique(X.pixel$admin2.name.full))), 
                                                                                                                                                                                                                                                                                                             rep("full", length(unique(X.pixel$admin2.name.full)))))
      post.u <- apply(draw.u1, 2, mean, na.rm = TRUE)
      post.r <- apply(draw.r1, 2, mean, na.rm = TRUE)
      post.all <- apply(draw.all1, 2, mean)
      post.u.sd <- apply(draw.u1, 2, sd, na.rm = TRUE)
      post.r.sd <- apply(draw.r1, 2, sd, na.rm = TRUE)
      post.all.sd <- apply(draw.all1, 2, sd)
      post.u.median <- apply(draw.u1, 2, median, na.rm = TRUE)
      post.r.median <- apply(draw.r1, 2, median, na.rm = TRUE)
      post.all.median <- apply(draw.all1, 2, median)
      post.u.ci <- apply(draw.u1, 2, quantile, probs = c((1 - 
                                                            CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
      post.r.ci <- apply(draw.r1, 2, quantile, probs = c((1 - 
                                                            CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
      post.all.ci <- apply(draw.all1, 2, quantile, probs = c((1 - 
                                                                CI)/2, 1 - (1 - CI)/2))
      admin1.bb.res <- data.frame(admin1.name = rep(unique(admin.info$admin1.name), 
                                                    3), mean = c(post.u, post.r, post.all), median = c(post.u.median, 
                                                                                                       post.r.median, post.all.median), sd = c(post.u.sd, 
                                                                                                                                               post.r.sd, post.all.sd), var = c(post.u.sd^2, 
                                                                                                                                                                                post.r.sd^2, post.all.sd^2), lower = c(post.u.ci[1, 
                                                                                                                                                                                ], post.r.ci[1, ], post.all.ci[1, ]), upper = c(post.u.ci[2, 
                                                                                                                                                                                ], post.r.ci[2, ], post.all.ci[2, ]), cv = c(cvf(post.u.median, 
                                                                                                                                                                                                                                 post.u.sd), cvf(post.r.median, post.r.sd), cvf(post.all.median, 
                                                                                                                                                                                                                                                                                post.all.sd)), type = c(rep("urban", length(unique(X.pixel$admin1.name))), 
                                                                                                                                                                                                                                                                                                        rep("rural", length(unique(X.pixel$admin1.name))), 
                                                                                                                                                                                                                                                                                                        rep("full", length(unique(X.pixel$admin1.name)))))
      post.u <- apply(draw.u0, 2, mean, na.rm = TRUE)
      post.r <- apply(draw.r0, 2, mean, na.rm = TRUE)
      post.all <- apply(draw.all0, 2, mean)
      post.u.sd <- apply(draw.u0, 2, sd, na.rm = TRUE)
      post.r.sd <- apply(draw.r0, 2, sd, na.rm = TRUE)
      post.all.sd <- apply(draw.all0, 2, sd)
      post.u.median <- apply(draw.u0, 2, median, na.rm = TRUE)
      post.r.median <- apply(draw.r0, 2, median, na.rm = TRUE)
      post.all.median <- apply(draw.all0, 2, median)
      post.u.ci <- apply(draw.u0, 2, quantile, probs = c((1 - 
                                                            CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
      post.r.ci <- apply(draw.r0, 2, quantile, probs = c((1 - 
                                                            CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
      post.all.ci <- apply(draw.all0, 2, quantile, probs = c((1 - 
                                                                CI)/2, 1 - (1 - CI)/2))
      admin0.bb.res <- data.frame(mean = c(post.u, post.r, 
                                           post.all), median = c(post.u.median, post.r.median, 
                                                                 post.all.median), sd = c(post.u.sd, post.r.sd, 
                                                                                          post.all.sd), var = c(post.u.sd^2, post.r.sd^2, 
                                                                                                                post.all.sd^2), lower = c(post.u.ci[1, ], post.r.ci[1, 
                                                                                                                ], post.all.ci[1, ]), upper = c(post.u.ci[2, 
                                                                                                                ], post.r.ci[2, ], post.all.ci[2, ]), cv = c(cvf(post.u.median, 
                                                                                                                                                                 post.u.sd), cvf(post.r.median, post.r.sd), cvf(post.all.median, 
                                                                                                                                                                                                                post.all.sd)), type = c(rep("urban", 1), rep("rural", 
                                                                                                                                                                                                                                                             1), rep("full", 1)))
      cm = list(agg.admin2 = admin2.bb.res, agg.admin1 = admin1.bb.res, 
                agg.admin0 = admin0.bb.res, inla = imod, admin2_post = draw.all, 
                urban2_post = draw.u, rural2_post = draw.r, 
                admin1_post = draw.all1, urban1_post = draw.u1, 
                rural1_post = draw.r1, admin0_post = draw.all0, 
                urban0_post = draw.u0, rural0_post = draw.r0, 
                admin.info = admin.info.output, admin = admin)
      attr(cm, "class") = "clusterModel"
      if (admin == 1) {
        attr(cm, "domain.names") <- admin.info$admin1.name
      }
      else {
        attr(cm, "domain.names") <- admin.info$admin2.name.full
      }
      return(cm)
    }
    else {
      post.all <- apply(draw.all, 2, mean)
      post.all.sd <- apply(draw.all, 2, sd)
      post.all.median <- apply(draw.all, 2, median)
      post.all.ci <- apply(draw.all, 2, quantile, probs = c((1 - 
                                                               CI)/2, 1 - (1 - CI)/2))
      if (admin == 2) {
        admin2.bb.res <- data.frame(admin2.name.full = sort(unique(X.pixel$admin2.name.full)), 
                                    mean = c(post.all), median = c(post.all.median), 
                                    sd = c(post.all.sd), var = c(post.all.sd^2), 
                                    lower = c(post.all.ci[1, ]), upper = c(post.all.ci[2, 
                                    ]), cv = cvf(post.all.median, post.all.sd))
      }
      else {
        admin2.bb.res <- data.frame(admin2.name.full = sort(unique(X.pixel$admin2.name.full)), 
                                    mean = c(post.all), median = c(post.all.median), 
                                    sd = c(post.all.sd), var = c(post.all.sd^2), 
                                    lower = c(post.all.ci[1, ]), upper = c(post.all.ci[2, 
                                    ]), cv = cvf(post.all.median, post.all.sd))
      }
      post.all <- apply(draw.all1, 2, mean)
      post.all.sd <- apply(draw.all1, 2, sd)
      post.all.median <- apply(draw.all1, 2, median)
      post.all.ci <- apply(draw.all1, 2, quantile, probs = c((1 - 
                                                                CI)/2, 1 - (1 - CI)/2))
      admin1.bb.res <- data.frame(admin1.name = unique(admin.info$admin1.name), 
                                  mean = c(post.all), median = c(post.all.median), 
                                  sd = c(post.all.sd), var = c(post.all.sd^2), 
                                  lower = c(post.all.ci[1, ]), upper = c(post.all.ci[2, 
                                  ]), cv = cvf(post.all.median, post.all.sd))
      post.all <- apply(draw.all0, 2, mean)
      post.all.sd <- apply(draw.all0, 2, sd)
      post.all.median <- apply(draw.all0, 2, median)
      post.all.ci <- apply(draw.all0, 2, quantile, probs = c((1 - 
                                                                CI)/2, 1 - (1 - CI)/2))
      admin0.bb.res <- data.frame(mean = c(post.all), 
                                  median = c(post.all.median), sd = c(post.all.sd), 
                                  var = c(post.all.sd^2), lower = c(post.all.ci[1, 
                                  ]), upper = c(post.all.ci[2, ]), cv = cvf(post.all.median, 
                                                                            post.all.sd))
      cm = list(res.admin2 = admin2.bb.res, res.admin1 = admin1.bb.res, 
                res.admin0 = admin0.bb.res, inla = imod, admin2_post = draw.all, 
                urban2_post = draw.u, rural2_post = draw.r, 
                admin1_post = draw.all1, urban1_post = draw.u1, 
                rural1_post = draw.r1, admin0_post = draw.all0, 
                urban0_post = draw.u0, rural0_post = draw.r0, 
                admin.info = admin.info.output, admin = admin)
      attr(cm, "class") = "clusterModel"
      if (admin == 1) {
        attr(cm, "domain.names") <- admin.info$admin1.name
      }
      else {
        attr(cm, "domain.names") <- admin.info$admin2.name.full
      }
      return(cm)
    }
  }
}
