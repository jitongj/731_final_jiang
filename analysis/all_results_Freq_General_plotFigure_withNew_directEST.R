# *******************
# task: This script is for simulation data of Zambia, 
# author: Jitong Jiang
# ********************
library(purrr)
library(ggrepel)
library(devtools)
library(surveyPrev)
library(INLA)
library(geodata)
library(ggpattern)
library(SUMMER)
library(rdhs)
library(ggplot2)
library(patchwork)
library(tidyr)
library(kableExtra)
library(png)
library(grid) 
library(sf)
library(viridis)
library(gridExtra)
library(here)
library(dplyr)
library(purrr)
library(ggplot2)
library(matrixStats)
library(forcats)


## ======= sim_f1_s1_p1_r0.042_nationalWeight.RData ======= ## You can just load this data and run the plots!
## load(here::here("data/sim_f1_s1_p1_r0.042_nationalWeight.RData"))
scale_cluster_frame  <- 1.0    # eg 0.5, 1, 2
scale_cluster_sample <- 1.0    # eg 0.5, 1, 2
scale_population     <- 1.0    # eg 0.5, 1, 2
prevalence_value     <- 0.042  # prevalence rate
## ================================= ##


# **********************
# 1. Basic setting   -------
# *********************
source(here::here("source", "helper.R"))
source(here::here("source", "directEST_national.R"))
source(here::here("source", "clusterModel_new.R"))
source(here::here("source", "fhModel_new.R"))

indicator <- "CN_NUTS_C_WH2" 

year <- 2018
frame_year <-2010
country <- "Zambia"
country.abbrev = "ZMB"


dhsData <- getDHSdata(country = country, indicator = indicator, year = year)
data <- getDHSindicator(dhsData, indicator = indicator)
data0 <- data %>%
  filter(!is.na(value))


geo <- getDHSgeo(country = country, year = year)
poly.adm1 <- geodata::gadm(country=country.abbrev, level=1, path=tempdir())
poly.adm1 <- sf::st_as_sf(poly.adm1)
poly.adm2 <- geodata::gadm(country=country.abbrev, level=2, path=tempdir())
poly.adm2 <- sf::st_as_sf(poly.adm2) %>%
  mutate(admin2.name.full = paste0(NAME_1, "_", NAME_2))

cluster.info <- clusterInfo(geo=geo, poly.adm1=poly.adm1, poly.adm2=poly.adm2, 
                            by.adm1 = "NAME_1",by.adm2 = "NAME_2")
admin.info1 <- adminInfo(poly.adm = poly.adm1, admin = 1, by.adm = "NAME_1")
admin.info2 <- adminInfo(poly.adm = poly.adm2, admin = 2,
                         by.adm = "NAME_2", by.adm.upper = "NAME_1")


# **********************
# 2. Prepare pop den for admin1 and admin2   -------
# **********************
setwd(here::here("data", country))

country_shp_analysis <- readRDS('country_shp_analysis.rds')

pop.abbrev <- tolower(country.abbrev)
pop_file <- paste0(pop.abbrev,'_ppp_',frame_year,'_1km_Aggregated_UNadj.tif')

if(!file.exists(pop_file)){
  if (frame_year < 2021){
    url <- paste0("https://data.worldpop.org/GIS/Population/Global_2000_2020_1km_UNadj/", 
                  frame_year, "/", toupper(pop.abbrev),"/",      
                  pop.abbrev,'_ppp_',frame_year,'_1km_Aggregated_UNadj.tif')
    
    download.file(url, pop_file, method = "libcurl", mode = "wb")
    
  } else{
    # download female and male tif
    url_female <- paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2021_2022_1km_UNadj/unconstrained/",
                         frame_year, "/total_female_male/", toupper(pop.abbrev), "/",      
                         pop.abbrev,'_f_total_', frame_year,'_1km_UNadj.tif')
    
    url_male <- paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2021_2022_1km_UNadj/unconstrained/",
                       frame_year, "/total_female_male/", toupper(pop.abbrev), "/",      
                       pop.abbrev,'_m_total_', frame_year,'_1km_UNadj.tif')
    
    file_female <- paste0("female_", pop.abbrev, "_", frame_year, ".tif")
    file_male <- paste0("male_", pop.abbrev, "_", frame_year, ".tif")
    
    download.file(url_female, file_female, method = "libcurl", mode = "wb")
    download.file(url_male, file_male, method = "libcurl", mode = "wb")
    
    # combine
    rast_female <- terra::rast(file_female)
    rast_male <- terra::rast(file_male)
    
    total_pop <- rast_female + rast_male
    terra::writeRaster(total_pop, pop_file, overwrite = TRUE)
  }
}

if (!exists("total_pop")) {
  total_pop <- terra::rast(pop_file)
}

#frame_pop <- terra::rast('nga_ppp_2006_1km_Aggregated_UNadj.tif') # pop at frame year
k0_1_pop <- total_pop #terra::rast('nga_k0_1_2018_1km.tif') # kid from 0 to 1 at survey year

pixel_grid <- as.data.frame(terra::xyFromCell(total_pop, 1:ncell(total_pop)))

pixel_adm_grid <- get_pixel_adm_grid(pixel_grid = pixel_grid,
                                     admin.level = 2,
                                     admin.poly = country_shp_analysis[[2+1]])


# pixel pop
pixel_adm_grid$pop_den <- terra::extract(total_pop, pixel_adm_grid[,c('x','y')])[,2]
pixel_adm_grid$k0_1_pop <- terra::extract(k0_1_pop, pixel_adm_grid[,c('x','y')])[,2]
pixel_adm_grid <- pixel_adm_grid %>%
  filter(!is.na(pop_den))  


colnames(pixel_adm_grid) <- c("x","y","admin2.name","admin1.name", "admin2.name.full", "pop_den", "k0_1_pop")

# calculate admin1 total pop, k0-1
admin1_pop <- pixel_adm_grid %>%
  group_by(admin1.name) %>%
  summarise(
    total_pop_admin1 = sum(pop_den, na.rm = TRUE),
    k0_1_pop_admin1 = sum(k0_1_pop, na.rm = TRUE)
  )


# calculate admin2 total pop, k0-1
admin2_pop <- pixel_adm_grid %>%
  group_by(admin2.name.full) %>%
  summarise(
    total_pop_admin2 = sum(pop_den, na.rm = TRUE),
    k0_1_pop_admin2 = sum(k0_1_pop, na.rm = TRUE)
  )



# **********************
# 3. direct estimate   -------
# **********************

res_ad1 <- directEST_national(data = data,
                     cluster.info = cluster.info,
                     admin = 1)
options(survey.adjust.domain.lonely = TRUE)
options(survey.lonely.psu = "adjust")
res_ad2 <- directEST_national(data = data,
                     cluster.info = cluster.info,
                     admin = 2,
                     aggregation = FALSE,
                     var.fix = FALSE, all.fix=FALSE)
res_ad2_fix <- directEST_national(data = data,
                                cluster.info = cluster.info,
                                admin = 2,
                                aggregation = FALSE,
                                var.fix = TRUE, all.fix=FALSE)

## individual level data
myData_tmp <- data0 %>%
  group_by(cluster) %>%   
  dplyr::summarise(      
    strata = unique(strata),  
    Ntrials = n(),      
    value = sum(value),   
    households_number = n_distinct(householdID)
  )
# combine individual level data with cluster info
myData <-merge(myData_tmp, cluster.info[["data"]], by="cluster")

# get admin2's direct estimate
direct <- res_ad2$res.admin2 %>%  
  filter(!is.na(direct.est)) %>% 
  distinct()


# **********************
# 4. calculate some admin2 data   -------
# **********************

# Calculate the number of admin2 instances in each admin1 instance in dhs.
count_admin2 <- myData %>%
  group_by(admin1.name) %>%
  summarise(count_admin2 = length(unique(admin2.name.full)), .groups = "drop")

# Calculate the number of admin2 var << 1e-30 in each admin1 in dhs.
count_admin2_direct_var_0 <- direct %>%
  group_by(admin1.name) %>%
  summarise(count_admin2_direct_var_0 = sum(direct.var < 1e-30 , na.rm = TRUE), .groups = "drop")



# Calculate the proportion of admin2 with var = 0 in each admin1 in dhs.
count_admin2_direct_var_0 <- count_admin2_direct_var_0 %>%
  left_join(count_admin2, by = "admin1.name") %>%
  mutate(ratio = count_admin2_direct_var_0 / count_admin2)

# Calculate the number of clusters in each admin2 in dhs.
sampled_admin2_cluster <- myData %>%
  group_by(admin1.name, admin2.name, admin2.name.full) %>%
  summarise(sampled_admin2_cluster = n(), .groups = "drop")



# **********************
# 5. total number of clusters (i.e., in census)  -------
# **********************

#### Load admin1 total cluster:
frame_ea <- readRDS(here::here("data", country, paste0(pop.abbrev, "_frame_ea.rds")))


#### Estimate admin2 total cluster:

# Merge admin2's direct estimate
## Merge admin2's direct estimate
direct_merge <- merge(
  direct[, c("admin1.name", "admin2.name", "admin2.name.full",
             "direct.est", "direct.var","direct.se")], 
  count_admin2_direct_var_0, by = c("admin1.name")
)
direct_merge <- merge(
  direct_merge, sampled_admin2_cluster,
  by = c("admin1.name", "admin2.name", "admin2.name.full")
)
direct_merge <- merge(direct_merge, admin1_pop, by = "admin1.name")
direct_merge <- merge(direct_merge, admin2_pop, by = "admin2.name.full")

direct_merge <- direct_merge %>%
  dplyr::mutate(
    total_pop_admin1 = round(total_pop_admin1 * scale_population),
    total_pop_admin2 = round(total_pop_admin2 * scale_population)
  )


direct_merge <- direct_merge %>%
  dplyr::left_join(
    frame_ea %>% dplyr::select(admin1.name = strata, total_cluster_admin1 = total),
    by = "admin1.name"
  ) %>%
  dplyr::mutate(
    total_cluster_admin1_scaled = round(total_cluster_admin1 * scale_cluster_frame),
    
    ratio_admin2 = total_pop_admin2 / total_pop_admin1,
    estimated_cluster_admin2 = round(ratio_admin2 * total_cluster_admin1_scaled)
  )


poly.adm2 <- poly.adm2 %>%
  dplyr::mutate(admin2.name.full = paste(NAME_1, NAME_2, sep = "_"))


# **********************
# 6. number of individual per clusters (i.e., in census)   -------
# **********************
set.seed(2024)

admin2_cluster_pop <- purrr::pmap_dfr(
  list(
    admin1_name = direct_merge$admin1.name,
    admin2_name = direct_merge$admin2.name.full,
    total_pop = direct_merge$total_pop_admin2,
    N_clusters = direct_merge$estimated_cluster_admin2
  ),
  function(admin1_name, admin2_name, total_pop, N_clusters) {
    if (is.na(N_clusters) || N_clusters == 0) return(NULL)
    
    raw_weights <- rgamma(N_clusters, 1)
    cluster_weights <- raw_weights / sum(raw_weights)
    M_c_total <- pmax(round(cluster_weights * total_pop), 30)
    M_c_total <- round(M_c_total / sum(M_c_total) * total_pop)
    
    tibble(
      admin1.name = admin1_name,
      admin2.name.full = admin2_name,
      cluster_id_within_admin2 = 1:N_clusters,
      M_c_total = M_c_total
    )
  }
) %>%
  mutate(cluster_global_id = row_number())  



# **********************
# 7. direct merge   -------
# **********************

direct_merge_org <- direct_merge

direct_merge <- direct_merge %>%
  group_by(admin1.name) %>%
  mutate(
    mean_est_cluster = round(mean(estimated_cluster_admin2[estimated_cluster_admin2 >= 100], na.rm = TRUE)),
    estimated_cluster_admin2 = ifelse(estimated_cluster_admin2 < 100, mean_est_cluster, estimated_cluster_admin2),
    total_cluster_admin1 = sum(estimated_cluster_admin2, na.rm = TRUE),
    ratio_admin2 = estimated_cluster_admin2 / total_cluster_admin1,
    total_pop_admin1 = unique(total_pop_admin1),
    total_pop_admin2 = round(total_pop_admin1 * ratio_admin2),
    estimated_cluster_admin2 = round(estimated_cluster_admin2),
    total_cluster_admin1 = round(total_cluster_admin1)
  ) %>%
  ungroup() %>%
  select(-mean_est_cluster)



q_admin2 <- myData %>%
  dplyr::group_by(admin2.name.full) %>%
  dplyr::summarise(
    n = dplyr::n(),
    n_u = sum(strata == "urban", na.rm = TRUE),
    q_hat = ifelse(n > 0, n_u / n, NA_real_)
  ) %>%
  dplyr::ungroup()

admin2_cluster_pop <- purrr::pmap_dfr(
  list(
    admin1_name      = direct_merge$admin1.name,
    admin2_name_full = direct_merge$admin2.name.full,
    admin2_name      = direct_merge$admin2.name,   
    total_pop        = direct_merge$total_pop_admin2,
    N_clusters       = direct_merge$estimated_cluster_admin2
  ),
  function(admin1_name, admin2_name_full, admin2_name, total_pop, N_clusters) {
    if (is.na(N_clusters) || N_clusters == 0) return(NULL)
    
    # 1) Allocate cluster size
    raw_weights     <- rgamma(N_clusters, 1)
    cluster_weights <- raw_weights / sum(raw_weights)
    M_c_total <- pmax(round(cluster_weights * total_pop), 30L)
    M_c_total <- round(M_c_total / sum(M_c_total) * total_pop)
    
    # 2) City target percentage (if none, add it back to admin1 or the national logic).
    q2 <- q_admin2$q_hat[match(admin2_name_full, q_admin2$admin2.name.full)]
    q2 <- min(max(q2, 0), 1)
    
    # 3) Number of target city clusters
    N_u <- min(N_clusters, max(0, round(q2 * N_clusters)))
    
    # 4) Select city clusters by population weighting
    if (N_u == 0) {
      strata_vec <- rep("rural", N_clusters)
    } else if (N_u == N_clusters) {
      strata_vec <- rep("urban", N_clusters)
    } else {
      urban_idx  <- sample(seq_len(N_clusters), size = N_u, prob = M_c_total, replace = FALSE)
      strata_vec <- ifelse(seq_len(N_clusters) %in% urban_idx, "urban", "rural")
    }
    
    tibble::tibble(
      admin1.name             = admin1_name,
      admin2.name.full        = admin2_name_full,
      admin2.name             = admin2_name,  
      cluster_id_within_admin2= seq_len(N_clusters),
      M_c_total               = as.integer(M_c_total),
      strata                  = strata_vec
    )
  }
) %>%
  dplyr::mutate(cluster_global_id = dplyr::row_number())


admin1_list <- unique(direct_merge$admin1.name)

# **********************
# 8. Run all admin1 and save raw results ----
# **********************
n_iterations = 100

all_admin1_out <- run_all_admin1_analysis(
  admin1_list = admin1_list,
  admin2_cluster_pop = admin2_cluster_pop,
  res_ad1 = res_ad1,
  prevalence_value = prevalence_value,
  n_iterations = n_iterations,
  n_indiv = 30,
  seed = 2024
)

save(
  all_admin1_out,
  file = file.path("data", "all_admin1_out.RData")
)


# ***************************
# 9. Ploting -----
# ***************************
load("~/Desktop/731/731_final/data/sim_data.RData")

sim_results_all_raw <- all_admin1_out$sim_results_all
by_admin1_results <- all_admin1_out$by_admin1

## Coverage

coverage_summary_all <- sim_results_all_raw %>%
  group_by(admin1.name, method) %>%
  summarise(
    mean_coverage = mean(covered, na.rm = TRUE),
    .groups = "drop"
  )

coverage_by_sim_all <- sim_results_all_raw %>%
  group_by(admin1.name, simulation, method) %>%
  summarise(
    coverage = mean(covered, na.rm = TRUE),
    .groups = "drop"
  )

p_cov_box_all <- ggplot(
  coverage_by_sim_all,
  aes(x = admin1.name, y = coverage, fill = method)
) +
  geom_boxplot(
    width = 0.5,   
    alpha = 0.8,
    position = position_dodge(width = 0.7) 
  ) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(
    x = "Admin1 Region",
    y = "Coverage",
    fill = "Method",
    title = "Coverage"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(p_cov_box_all)


## Bias

bias_summary_all <- sim_results_all_raw %>%
  group_by(admin1.name, method) %>%
  summarise(
    mean_bias = mean(bias, na.rm = TRUE),
    .groups = "drop"
  )

bias_by_sim_all <- sim_results_all_raw %>%
  group_by(admin1.name, simulation, method) %>%
  summarise(
    bias = mean(bias, na.rm = TRUE),
    .groups = "drop"
  )

p_bias_box_all <- ggplot(
  bias_by_sim_all,
  aes(x = admin1.name, y = bias, fill = method)
) +
  geom_boxplot(
    width = 0.5,
    alpha = 0.8,
    position = position_dodge(width = 0.7)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "Admin1 Region",
    y = "Bias",
    fill = "Method",
    title = "Bias"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(p_bias_box_all)

## Interval Width

width_summary_all <- sim_results_all_raw %>%
  group_by(admin1.name, method) %>%
  summarise(
    mean_width = mean(width, na.rm = TRUE),
    .groups = "drop"
  )

width_by_sim_all <- sim_results_all_raw %>%
  group_by(admin1.name, simulation, method) %>%
  summarise(
    width = mean(width, na.rm = TRUE),
    .groups = "drop"
  )

p_width_box_all <- ggplot(
  width_by_sim_all,
  aes(x = admin1.name, y = width, fill = method)
) +
  geom_boxplot(
    width = 0.5,
    alpha = 0.8,
    position = position_dodge(width = 0.7)
  ) +
  labs(
    x = "Admin1 Region",
    y = "Interval Width",
    fill = "Method",
    title = "Interval Width"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(p_width_box_all)


## Save summaries

summary_tables_all <- list(
  sim_results_all_raw = sim_results_all_raw,
  coverage_summary_all = coverage_summary_all,
  bias_summary_all = bias_summary_all,
  width_summary_all = width_summary_all
)

summary_tables_all