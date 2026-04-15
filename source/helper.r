# inverse logit ---------
expit <- function(x) 1 / (1 + exp(-x))

# logit ---------
logit <- function(p) log(p / (1 - p))

# Map pixel coordinates to administrative regions ---------
get_pixel_adm_grid <- function(pixel_grid,
                               admin.level,
                               admin.poly){
  
  # Match stratification admin
  pixel_grid_sf <- st_as_sf(pixel_grid, coords = c("x", "y"),
                            crs = st_crs(admin.poly))
  
  adm_match <- st_join(pixel_grid_sf,admin.poly, join = st_intersects)
  
  pixel_grid_df <- data.frame(
    x = pixel_grid$x,
    y = pixel_grid$y
  )
  
  if(admin.level==1){
    pixel_grid_df$adm.name =adm_match[[paste0('NAME_1')]]
    pixel_grid_df$adm.name.full = pixel_grid_df$adm.name
  }else{
    
    pixel_grid_df$adm.name = adm_match[[paste0('NAME_',admin.level)]]
    pixel_grid_df$upper.adm.name = adm_match[[paste0('NAME_',admin.level-1)]]
    pixel_grid_df$adm.name.full = paste0(adm_match[[paste0('NAME_',admin.level-1)]],
                                         '_',
                                         adm_match[[paste0('NAME_',admin.level)]])
  }
  
  
  
  return(pixel_grid_df[complete.cases(pixel_grid_df),])
}


# simulate data --------
run_sampling_stratum_level <- function(
    admin1.name,
    cluster_frame,
    n_iterations = n_iterations,
    n_clusters = NULL,       # number of sampled cluster
    n_indiv = 30,            # numbered of sampled individual in each cluster
    sigma0 = 0.5,            # Admin2 sigma
    sigma1 = 0.2,            # Cluster sigma
    sigma2 = 0.05,           # invidual sigma
    res_ad1 = NULL,          # Admin1 direct.est
    m = prevalence_value,
    coverage_level = 0.8,
    seed = 2024
) {
  set.seed(seed)
  
  
  cluster_sub <- cluster_frame[cluster_frame$admin1.name == admin1.name, ]
  cluster_ids <- cluster_sub$cluster_global_id
  M_c_total <- cluster_sub$M_c_total
  admin2_names_full <- cluster_sub$admin2.name.full
  admin2_names <- cluster_sub$admin2.name
  N <- nrow(cluster_sub)
  strata_vec <- cluster_sub$strata
  
  
  if (is.null(n_clusters)) {
    base_n <- sum(direct_merge$sampled_admin2_cluster[direct_merge$admin1.name == admin1.name], na.rm = TRUE)
    
    max_n <- admin2_cluster_pop %>%
      dplyr::filter(admin1.name == !!admin1.name) %>%
      dplyr::summarise(N = dplyr::n()) %>%
      dplyr::pull(N)
    
    n_clusters <- max(1L, min(as.integer(round(base_n * scale_cluster_sample)), max_n))
  }
  
  
  # baseline m
  if (is.null(m)) {
    m <- res_ad1$res.admin1$direct.est[res_ad1$res.admin1$admin1.name == admin1.name]
  }
  
  # alpha (admin2 level random effect) + e_c (cluster errir)
  alpha_admin2_map <- rnorm(length(unique(admin2_names_full)), 0, sigma0)
  names(alpha_admin2_map) <- unique(admin2_names_full)
  e_c <- rnorm(N, 0, sigma1)
  
  
  individual_df <- purrr::pmap_dfr(
    list(
      cluster_id = cluster_ids,
      M_c        = M_c_total,
      e_cluster  = e_c,
      admin2_full= admin2_names_full,
      admin2     = admin2_names,
      strata_c   = strata_vec          
    ),
    function(cluster_id, M_c, e_cluster, admin2_full, admin2, strata_c) { 
      alpha_admin2 <- alpha_admin2_map[admin2_full]
      e_ck <- rnorm(M_c, 0, sigma2)
      eta <- logit(m) + alpha_admin2 + e_cluster + e_ck   
      p_ck <- expit(eta)
      y <- rbinom(M_c, 1, p_ck)
      tibble(
        cluster_id = cluster_id,
        admin2.name.full = admin2_full,
        admin2.name = admin2,
        y = y,
        strata = strata_c            
      )
    }
  ) %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::mutate(
      householdID = dplyr::row_number(),
      value = y,
      cluster = cluster_id,
      weight = 1,
      admin1.name = admin1.name
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(value, admin2.name.full, admin2.name, cluster, householdID, weight, admin1.name, strata)
  
  
  
  sampled_modt_list <- vector("list", n_iterations)
  
  for (i in 1:n_iterations) {
    # first stage PPS sampling cluster
    total_pop_admin1 <- unique(direct_merge$total_pop_admin1[direct_merge$admin1.name == admin1.name])
    sampled_clusters <- sample(
      cluster_ids,
      size = n_clusters,
      prob = n_clusters * M_c_total / total_pop_admin1,
      replace = FALSE
    )
    
    sampled_df <- individual_df[individual_df$cluster %in% sampled_clusters, ]
    
    # second stage: sample n_indiv in each cluster and gives weight
    sampled_indiv <- sampled_df %>%
      group_by(cluster) %>%
      slice_sample(n = n_indiv) %>%
      ungroup() %>%
      mutate(weight = round(total_pop_admin1 / (n_clusters * n_indiv)))
    
    sampled_modt_list[[i]] <- sampled_indiv
  }
  
  return(list(
    full_population = individual_df,
    sampled_indiv_list = sampled_modt_list
  ))
}



# models analysis based on one sim dataset--------
process_one_sim <- function(sampled_data, truth_df, admin.info2) {
  

  ## Direct estimation
  res_ad2 <- directEST_national(
    data = sampled_data,
    cluster.info = NULL,
    admin = 2,
    aggregation = FALSE,
    var.fix = FALSE,
    all.fix = FALSE
  )$res.admin2
  

  ## Identify bad areas for FH
  res_ad2_fix <- directEST_national(
    data = sampled_data,
    cluster.info = NULL,
    admin = 2,
    aggregation = FALSE,
    var.fix = TRUE
  )
  
  bad_admin2 <- res_ad2_fix$fixed_areas
  bad_clusters <- unique(
    subset(sampled_data, admin2.name.full %in% bad_admin2)$cluster
  )
  
  ## only keep areas that actually appear in the sample
  areas_in_sample <- unique(sampled_data$admin2.name.full)
  
  admin.info2_sub <- list(
    data = admin.info2$data %>%
      dplyr::filter(admin2.name.full %in% areas_in_sample),
    mat = NULL
  )
  

  ## FH model
  smth_res_ad2 <- fhModel_new(
    data = subset(sampled_data, !cluster %in% bad_clusters),
    cluster.info = NULL,
    admin.info = admin.info2_sub,
    admin = 2,
    model = "iid",
    aggregation = FALSE,
    var.fix = FALSE,
    nested = FALSE
  )
  

  ## Unit-level model

  cl_res_ad2 <- clusterModel_new(
    data = sampled_data,
    cluster.info = NULL,
    admin.info = admin.info2_sub,
    model = "iid",
    stratification = FALSE,
    admin = 2,
    aggregation = FALSE,
    CI = 0.95
  )
  

  ## Tidy outputs
  direct_df <- res_ad2 %>%
    transmute(
      admin2.name.full,
      method = "Direct",
      mean = direct.est,
      lower = direct.lower,
      upper = direct.upper
    )
  
  fh_df <- smth_res_ad2$res.admin2 %>%
    transmute(
      admin2.name.full,
      method = "FH",
      mean = mean,
      lower = lower,
      upper = upper
    )
  
  unit_df <- cl_res_ad2$res.admin2 %>%
    transmute(
      admin2.name.full,
      method = "Unit-level",
      mean = mean,
      lower = lower,
      upper = upper
    )
  
  ## Keep only areas where direct exists, so all three methods are comparable
  areas_in_direct <- unique(direct_df$admin2.name.full)
  
  fh_df <- fh_df %>%
    filter(admin2.name.full %in% areas_in_direct)
  
  unit_df <- unit_df %>%
    filter(admin2.name.full %in% areas_in_direct)
  
  ## Merge with truth and compute metrics
  bind_rows(direct_df, fh_df, unit_df) %>%
    left_join(truth_df, by = "admin2.name.full") %>%
    mutate(
      covered = (true_mean >= lower & true_mean <= upper),
      width   = upper - lower,
      bias    = mean - true_mean
    )
}


# run simulations for one admin1  ----------

run_one_admin1_analysis <- function(ad1, admin2_cluster_pop, res_ad1,
                                    prevalence_value,
                                    n_iterations = n_iterations,
                                    n_indiv = 30,
                                    seed = 2024) {
  
  cat("Start admin1:", ad1, "...\n")
  
  res <- run_sampling_stratum_level(
    admin1.name   = ad1,
    cluster_frame = admin2_cluster_pop,
    n_iterations  = n_iterations,
    n_clusters    = NULL,
    n_indiv       = n_indiv,
    res_ad1       = res_ad1,
    m             = prevalence_value,
    seed          = seed
  )
  
  ## truth from full population
  truth_df <- res$full_population %>%
    group_by(admin2.name.full) %>%
    summarise(
      true_mean = mean(value),
      .groups = "drop"
    )
  
  ## simulation-level detailed results
  sim_results <- purrr::map_dfr(
    res$sampled_indiv_list,
    ~process_one_sim(.x, truth_df, admin.info2),
    .id = "simulation"
  ) %>%
    mutate(
      simulation = as.integer(simulation),
      admin1.name = ad1
    )
  
  cat("Finish admin1:", ad1, "\n")
  cat("---\n")
  
  list(
    res = res,
    truth_df = truth_df,
    sim_results = sim_results
  )
}



# run simulations for all admin1 ----------

run_all_admin1_analysis <- function(admin1_list, admin2_cluster_pop, res_ad1,
                                    prevalence_value,
                                    n_iterations = n_iterations,
                                    n_indiv = 30,
                                    seed = 2024) {
  
  out_list <- purrr::map(
    admin1_list,
    ~run_one_admin1_analysis(
      ad1 = .x,
      admin2_cluster_pop = admin2_cluster_pop,
      res_ad1 = res_ad1,
      prevalence_value = prevalence_value,
      n_iterations = n_iterations,
      n_indiv = n_indiv,
      seed = seed
    )
  )
  
  names(out_list) <- admin1_list
  
  sim_results_all <- purrr::map_dfr(out_list, "sim_results")
  
  list(
    by_admin1 = out_list,
    sim_results_all = sim_results_all
  )
}


