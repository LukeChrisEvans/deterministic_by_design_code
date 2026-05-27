# =========================================================
# Simulations 
# Four scenarios:
# 1. Process effects 
# 2. Confounding
# 3. Compression
# 4. Detection bias
#
# Three analysis types:
#  1. OLS with log species richness ~ time at each site, extract slopes and regress on Z
#  2. Construct multi-level trend model, extract slopes and regress on Z
#  3. Construct multi-level interaction model with Z*year
#
#
#
# Outputs:
#   - coefficient table
#   - raw trajectory panel plot
#   - slope-vs-Z panel plot
#
# LCE final change 27/05/26
# =========================================================

library(tidyverse)
library(patchwork)
library(lme4)
library(lmerTest)
library(flextable)
library(broom.mixed)
library(broom)
library(interplot)
library(officer)
set.seed(1234)

# =========================================================
# 1. Global parameters
# =========================================================

N_sites <- 200    # sites
T_years <- 20     # annual observations
z_sd <- 1         # SD of gradient Z
sigma_proc <- 0.8 # process noise

# Common site and time structure
site_df <- data.frame(
  site = factor(seq_len(N_sites)),
  Z = rnorm(N_sites, 0, z_sd)
)

time_df <- tibble(
  year = 1:(T_years)
)

base_df <- tidyr::crossing(site_df, time_df)
View(base_df)

# =========================================================
# 2. Generate data sets
# =========================================================


# =========================================================
### Process effect 
# =========================================================
process_df <- base_df
c0 <- 0.02 ;cZ <- -0.025 ; sigma_c <- 0.01

site_pars_process = process_df %>%
  distinct(site, Z) %>%
  mutate(Y_true = rpois(n(),10))

grow_vec_process <- c()
for(i in seq_along(site_pars_process$Y_true)){
  y_out <- c()
  log_lambda <- c()
  
  y_out[1] <- site_pars_process$Y_true[i]
  Z_here <- site_pars_process$Z[i]
  
  for(j in 2:20){
    eps <-  rnorm(1, 0, sigma_c)
    log_lambda[j] <- log(y_out[j-1] + 1e-6) + c0 + cZ * Z_here + eps
    y_out[j] <- rpois(1, exp(log_lambda[j]))
  }
  
  grow_vec_process <- c(grow_vec_process,y_out)
}


process_df$Y_true <- grow_vec_process
process_df$Y_obs <- rbinom(n=length(process_df$Y_true),
                           size = process_df$Y_true,
                           prob =1)


process_df %>% ggplot(aes(year,Y_obs,group=as.factor(site)))+
  geom_line(alpha = 0.35, linewidth = 0.3)

process_df %>% ggplot(aes(year,log(Y_obs+1),group=as.factor(site)))+
  geom_line(alpha = 0.35, linewidth = 0.3)


# =========================================================
### Confounding 
# =========================================================
confound_df = base_df
aU = -0.5 ;y0_intercept = 10; y0_sd = 0.25; sigma_z_c = 0.1; sigma_c = 0.1;c0 = 0.02;cZ = -0.05; sigma_proc = 0.08; bZU=0.08; Ua =0

site_pars_confound = confound_df %>%
  distinct(site, Z) %>%
  rename(U=Z)     %>%                    # # we swap Z for U here for convenience 
  mutate(log_lambda_y0 = log(y0_intercept) + aU * U + rnorm(n(), 0, y0_sd),
         Y_true = rpois(n(), lambda = exp(log_lambda_y0))
          )

grow_vec_confound <- c()
for(i in seq_along(site_pars_confound$Y_true)){
  y_out <- c()
  log_lambda <- c()
  
  y_out[1] <- site_pars_confound$Y_true[i]
  U_here <- site_pars_confound$U[i] 
  
  for(j in 2:20){
    Z2 <- Ua + bZU * U_here + rnorm(1, 0, sigma_z_c) # U influencing Z
    C  <- c0 + cZ * Z2  + rnorm(1, 0, sigma_c) # 
    
    eps <-  rnorm(1, 0, sigma_proc)
    log_lambda[j] <- log(y_out[j-1] + 1e-6) + C + eps
    y_out[j] <- rpois(1, exp(log_lambda[j]))
  }
  grow_vec_confound <- c(grow_vec_confound,y_out)
}


confound_df$Y_true <- grow_vec_confound
confound_df$Y_obs <- rbinom(n=length(confound_df$Y_true),
                            size = confound_df$Y_true,
                            prob =1)


confound_df %>% ggplot(aes(year,Y_obs,group=as.factor(site)))+
  geom_line(alpha = 0.35, linewidth = 0.3)

confound_df %>% ggplot(aes(year,log(Y_obs+1),group=as.factor(site)))+
  geom_line(alpha = 0.35, linewidth = 0.3)


# =========================================================
### Compression 
# =========================================================
compression_df <- base_df
aZ = 10 ;y0_intercept = 30; y0_sd = 1.5; step_sd = 1; K = 40 ; global_growth = 1.5

site_pars_compression <- compression_df %>%
  distinct(site, Z) %>%
  mutate(mu_y0 = y0_intercept + aZ * Z + rnorm(n(), 0, y0_sd),
         mu_y0 = pmax(1e-6, pmin(K, mu_y0)),
            
         log_lambda_y0 = log(mu_y0),
         
         Y_true = rpois(n(), exp(log_lambda_y0)),
         Y_true = pmax(1, pmin(K, Y_true))
  )

grow_vec_compression <- c()

for(i in seq_along(site_pars_compression$Y_true)){
  y_out <- c()
  y_out[1] <- site_pars_compression$Y_true[i]
  
  for(j in 2:20){
    mu_j <- y_out[j - 1] + global_growth + rnorm(1, 0, step_sd)
    mu_j <- pmax(1e-6, pmin(K, mu_j))
    log_lambda_j <- log(mu_j)
    y_out[j] <- rpois(1, exp(log_lambda_j))
    y_out[j] <- pmax(0, pmin(K, y_out[j]))
  }
  
  grow_vec_compression <- c(grow_vec_compression,y_out)
}

compression_df$Y_true <- grow_vec_compression
compression_df$Y_obs <- rbinom(n=length(compression_df$Y_true),
                           size = compression_df$Y_true,
                           prob =1)

#
compression_df %>% ggplot(aes(year,Y_obs,group=as.factor(site)))+
  geom_line(alpha = 0.35, linewidth = 0.3)

#
compression_df %>% ggplot(aes(year,log(Y_obs+1),group=as.factor(site)))+
  geom_line(alpha = 0.35, linewidth = 0.3)


# =========================================================
### Detection Bias
# =========================================================
detection_df <- base_df
detection_df$Y_true <- rpois(nrow(base_df), 10)

a0 = 0.35; alpha_g = 0.3;betaZ = 1.2     

max_year <- max(detection_df$year)

detection_df <- detection_df %>%
  mutate(
    weight = 1 - year / max_year,
    p_it = plogis(
      qlogis(a0) +
        alpha_g * year +
        betaZ * Z * weight
    )
  )

detection_df$Y_obs <- rbinom(nrow(detection_df),size = pmax(1,round(detection_df$Y_true)),prob =detection_df$p_it)

detection_df %>% ggplot(aes(year,Y_obs,group=as.factor(site)))+
  geom_line(alpha = 0.35, linewidth = 0.3)

detection_df %>% ggplot(aes(year,Y_true,group=as.factor(site)))+
  geom_line(alpha = 0.35, linewidth = 0.3)


# =========================================================
# 3. Conduct OLS with log species richness ~ time, then extract trend and regress on z
# =========================================================
# add logs
process_df$lg_Y_obs <- log(process_df$Y_obs +1)
confound_df$lg_Y_obs <- log(confound_df$Y_obs +1)
compression_df$lg_Y_obs <- log(compression_df$Y_obs +1)
detection_df$lg_Y_obs <- log(detection_df$Y_obs + 1)

run_site_slopes <- function(df){
  coefs <- c()
  Zs <- c()
  for(i in unique(df$site)){
    t_site <- filter(df,site==i)
    t_mod <- lm(lg_Y_obs ~ year ,data=t_site)
    coefs <- c(coefs,coef(t_mod)["year"])
    Zs <- c(Zs,t_site$Z[1])
  }
  
  cbind.data.frame(site=unique(df$site),coef=coefs,Z=Zs)
}

# extract obs ~ year
trends_process <- run_site_slopes(process_df)
trends_confound <- run_site_slopes(confound_df)
trends_compression <- run_site_slopes(compression_df)
trends_detection <- run_site_slopes(detection_df)

# coef ~ Z
process_model <-lm(coef~Z,data= trends_process)
confound_model <-lm(coef~Z,data= trends_confound)
compression_model <-lm(coef~Z,data= trends_compression)
detection_model <-lm(coef~Z,data= trends_detection)


summary(process_model)
summary(confound_model)
summary(compression_model)
summary(detection_model)

# TABLE & PLOTS
mod_list=list(
  "Process \n effect (log)" = process_model,
  "Confounding (log)" = confound_model,
  "Compression (log)" = compression_model,
  "Detection \n bias (log)" = detection_model
)


table_1 <- bind_rows(lapply(names(mod_list), function(nm) {
  tidy(mod_list[[nm]]) |>
    transmute(
      term,
      model = nm,
      cell = sprintf("%.3f (%.3f)", estimate, std.error)
    )
})) |>
  pivot_wider(names_from = model, values_from = cell)

table_1[is.na(table_1)] <- ""

# add the actul direct effects
direct_effect_row <- tibble(
  term = "True direct effect of Z",
  `Process \n effect (log)` = "-0.025", #  cZ
  `Confounding (log)` = "-0.05", # cZ
  `Compression (log)` = "0",
  `Detection \n bias (log)` = "0"
)


table_1 <- bind_rows(table_1, direct_effect_row)
flextable(table_1) |> autofit() # view table

# Create Word document and save Table
doc <- read_docx()
# Add flextable
doc <- body_add_flextable(doc, flextable(table_1) |> autofit())
# Save
print(doc, target = "C:/Users/lukeva/OneDrive - UKCEH/Documents/Becks paper/Tables/Table_1.docx")
print(doc, target = "~/Tables/Table_1.docx")


## plots
p1 <- process_df %>% ggplot(aes(year,Y_obs,group=site,color=Z))+
  geom_line(alpha = 0.35, linewidth = 0.5)+
  theme_classic()+
  scale_color_viridis_b()+
  xlab("Year")+
  ylab("Richness")+
  theme(legend.position = "none") + ggtitle( 'a) Process')

p2 <- confound_df %>% ggplot(aes(year,Y_obs,group=site,color=Z))+
  geom_line(alpha = 0.35, linewidth = 0.5)+
  theme_classic()+
  scale_color_viridis_b()+
  xlab("Year")+
  ylab("")+
  theme(legend.position = "none") + ggtitle('b) Confounding')

p3 <- compression_df %>% ggplot(aes(year,Y_obs,group=site,color=Z))+
  geom_line(alpha = 0.35, linewidth = 0.5)+
  theme_classic()+
  scale_color_viridis_b()+
  xlab("Year")+
  ylab("")+
  theme(legend.position = "none")+ ggtitle( 'c) Compression')

p4 <- detection_df %>% ggplot(aes(year,Y_obs,group=site,color=Z))+
  geom_line(alpha = 0.35, linewidth = 0.5)+
  theme_classic()+
  scale_color_viridis_b()+
  xlab("Year")+
  ylab("")+
  theme(legend.position = "none")+ ggtitle( 'd) Detection \nbias')

#
lg_p1 <- process_df %>% ggplot(aes(year,lg_Y_obs,group=site,color=Z))+
  geom_line(alpha = 0.35, linewidth = 0.5)+
  theme_classic()+
  scale_color_viridis_b()+
  xlab("Year")+
  ylab("log (richness+1)")+
  theme(legend.position = "none")

lg_p2 <- confound_df %>% ggplot(aes(year,lg_Y_obs,group=site,color=Z))+
  geom_line(alpha = 0.35, linewidth = 0.5)+
  theme_classic()+
  scale_color_viridis_b()+
  xlab("Year")+
  ylab("")+
  theme(legend.position = "none") 

lg_p3 <- compression_df %>% ggplot(aes(year,lg_Y_obs,group=site,color=Z))+
  geom_line(alpha = 0.35, linewidth = 0.5)+
  theme_classic()+
  scale_color_viridis_b()+
  xlab("Year")+
  ylab("")+
  theme(legend.position = "none")

lg_p4 <- detection_df %>% ggplot(aes(year,lg_Y_obs,group=site,color=Z))+
  geom_line(alpha = 0.35, linewidth = 0.5)+
  theme_classic()+
  scale_color_viridis_b()+
  xlab("Year")+
  ylab("")+
  theme(legend.position = "none")


# slope plots
slopes_p1 <- trends_process %>% ggplot(aes(Z,coef))+
  geom_point(color="blue",size=2,alpha=0.5)+
  theme_classic()+
  ylab("Biodiversity change")+
  xlab("Spatial gradient (Z)")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_smooth(aes(x = Z, y = coef), method = "lm", formula = y ~ x,color="black",se=F)

slopes_p2 <- trends_confound %>% ggplot(aes(Z,coef))+
  geom_point(color="blue",size=2,alpha=0.5)+
  theme_classic()+
  ylab("")+
  xlab("Spatial gradient (Z)")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_smooth(aes(x = Z, y = coef), method = "lm", formula = y ~ x,color="black",se=F)

slopes_p3 <- trends_compression %>% ggplot(aes(Z,coef))+
  geom_point(color="blue",size=2,alpha=0.5)+
  theme_classic()+
  ylab("")+
  xlab("Spatial gradient (Z)")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_smooth(aes(x = Z, y = coef), method = "lm", formula = y ~ x,color="black",se=F)

slopes_p4 <- trends_detection %>% ggplot(aes(Z,coef))+
  geom_point(color="blue",size=2,alpha=0.5)+
  theme_classic()+
  ylab("")+
  xlab("Spatial gradient (Z)")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_smooth(aes(x = Z, y = coef), method = "lm", formula = y ~ x,color="black",se=F)



panelplt <- (p1 | p2 | p3 | p4)/ (lg_p1 | lg_p2 | lg_p3| lg_p4)/ (slopes_p1 | slopes_p2 | slopes_p3 | slopes_p4) + plot_layout(tag_level = 'new')


ggsave("~/figures/panel_1.png",
       panelplt, bg="white")

# =========================================================
# 4. Construct multi-level Poisson model and extract slopes, 
#    next, regress slopes on Z using OLS
# =========================================================

run_multi_level_slopes <- function(df){
  
  multi_level <- lmer(lg_Y_obs ~ year + (1 + year | site), control = lmerControl(optimizer = "bobyqa"), data = df)
  
  data.frame(
    site  = rownames(coef(multi_level)$site),
    trend = coef(multi_level)$site[,"year"],
    Z     = distinct(df,Z)$Z
  )
}

# extract year trends from multilevel
trends_process_multi_level <- run_multi_level_slopes(process_df)
trends_confound_multi_level <- run_multi_level_slopes(confound_df)
trends_compression_multi_level <- run_multi_level_slopes(compression_df)
trends_detection_multi_level <- run_multi_level_slopes(detection_df)

# model with ols against Z
process_model_multi_level <-lm(trend~Z,data= trends_process_multi_level)
confound_model_multi_level <-lm(trend~Z,data= trends_confound_multi_level)
compression_model_multi_level <-lm(trend~Z,data= trends_compression_multi_level)
detection_model_multi_level <-lm(trend~Z,data= trends_detection_multi_level)


# =========================================================
# 5. Construct multi-level model with time-by-z interaction 
# =========================================================
run_interaction_models  <- function(df){

  scaled_df <- df %>% # scale for stability 
    mutate(
      year_sc = as.numeric(scale(year, center = TRUE, scale = TRUE)),
      Z_sc = as.numeric(scale(Z, center = TRUE, scale = TRUE))
    )
  
  pois_mod <- glmer(
    Y_obs ~ year_sc * Z_sc + (1 | site),
    family = poisson(link = "log"),
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
    data = scaled_df
  )
  
  
  pois_mod
}


interaction_process <- run_interaction_models(process_df)
interaction_confound <- run_interaction_models(confound_df)
interaction_compression <- run_interaction_models(compression_df)
interaction_detection <- run_interaction_models(detection_df)

# =========================================================
# 6. Results tables
# =========================================================
summary(process_model_multi_level)
summary(confound_model_multi_level)
summary(compression_model_multi_level)
summary(detection_model_multi_level)

# TABLE & PLOTS
multilevel_list=list(
  "Process \n effect (log)" = process_model_multi_level,
  "Confounding (log)" = confound_model_multi_level,
  "Compression (log)" = compression_model_multi_level,
  "Detection \n bias (log)" = detection_model_multi_level
)


table_s1 <- bind_rows(lapply(names(multilevel_list), function(nm) {
  tidy(multilevel_list[[nm]]) |>
    transmute(
      term,
      model = nm,
      cell = sprintf("%.3f (%.3f)", estimate, std.error)
    )
})) |>
  pivot_wider(names_from = model, values_from = cell)

table_s1[is.na(table_s1)] <- ""

# add the actul direct effects
direct_effect_row <- tibble(
  term = "True direct effect of Z",
  `Process \n effect (log)` = "-0.025", #  cZ
  `Confounding (log)` = "-0.05", # cZ
  `Compression (log)` = "0",
  `Detection \n bias (log)` = "0"
)


table_s1 <- bind_rows(table_s1, direct_effect_row)

# Create Word document and save Table
doc <- read_docx()
# Add flextable
doc <- body_add_flextable(doc, flextable(table_s1) |> autofit())
# Save
print(doc, target = "~/Tables/Table_S1.docx")


##################################
summary(interaction_process)
summary(interaction_confound)
summary(interaction_compression)
summary(interaction_detection)

interaction_list<-list(
  "Process \n effect" = interaction_process,
  "Confounding" = interaction_confound,
  "Compression" = interaction_compression,
  "Detection \n bias" = interaction_detection
  )


table_s2 <- bind_rows(lapply(names(interaction_list), function(nm) {
  tidy(interaction_list[[nm]]) |>
    mutate(
      term = recode(term,
                    `(Intercept)` = "Intercept",
                    `year_sc` = "Year (scaled)",
                    `Z_sc` = "Z (scaled)",
                    `year_sc:Z_sc` = "Year × Z (scaled)"
      )
    ) |>
    transmute(
      term,
      model = nm,
      cell = sprintf("%.3f (%.3f)", estimate, std.error)
    )
})) |>
  pivot_wider(names_from = model, values_from = cell)

table_s2[is.na(table_s2)] <- ""
table_s2<-table_s2[1:4,] # remove sd


# Create Word document
doc <- read_docx()
# Add flextable
doc <- body_add_flextable(doc, flextable(table_s2) |> autofit())

# Save
print(doc, target = "~/Tables/Table_S2.docx")



# =========================================================
# 7. Plots
# =========================================================
multi_level_p1 <- trends_process_multi_level %>% ggplot(aes(Z,trend))+
  geom_point(color="blue",size=2,alpha=0.5)+
  theme_classic()+
  ylab("Biodiversity change")+
  xlab("Spatial gradient (Z)")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_smooth(aes(x = Z, y = trend), method = "lm", formula = y ~ x,color="black",se=F)

multi_level_p2 <- trends_confound_multi_level %>% ggplot(aes(Z,trend))+
  geom_point(color="blue",size=2,alpha=0.5)+
  theme_classic()+
  ylab("")+
  xlab("Spatial gradient (Z)")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_smooth(aes(x = Z, y = trend), method = "lm", formula = y ~ x,color="black",se=F)

multi_level_p3 <- trends_compression_multi_level %>% ggplot(aes(Z,trend))+
  geom_point(color="blue",size=2,alpha=0.5)+
  theme_classic()+
  ylab("")+
  xlab("Spatial gradient (Z)")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_smooth(aes(x = Z, y = trend), method = "lm", formula = y ~ x,color="black",se=F)

multi_level_p4 <- trends_detection_multi_level %>% ggplot(aes(Z,trend))+
  geom_point(color="blue",size=2,alpha=0.5)+
  theme_classic()+
  ylab("")+
  xlab("Spatial gradient (Z)")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_smooth(aes(x = Z, y = trend), method = "lm", formula = y ~ x,color="black",se=F)



panelplt_s1 <- (lg_p1 + ggtitle( 'a) Process') | lg_p2 + ggtitle( 'b) Confounding')  | lg_p3 + ggtitle( 'c) Compression')| lg_p4 + ggtitle( 'd) Detection\n bias'))/ (multi_level_p1 | multi_level_p2 | multi_level_p3 | multi_level_p4) + plot_layout(tag_level = 'new')

ggsave("~/figures/panel_s1.png",
       panelplt_s1, bg="white")


############
# inter plots
i1 <- interplot(interaction_process, var1 = "year_sc", var2 = "Z_sc") +
  theme_classic()+
  xlab("Z (scaled)") +
  ylab("Year (scaled)") +
  geom_hline(yintercept = 0,lty="dashed")+
  ylim(-0.3,0.56)

i2 <- interplot(interaction_confound, var1 = "year_sc", var2 = "Z_sc") +
  theme_classic()+
  xlab("Z (scaled)") +
  ylab("") +
  geom_hline(yintercept = 0,lty="dashed")+
  ylim(-0.3,0.56)

i3 <- interplot(interaction_compression, var1 = "year_sc", var2 = "Z_sc") +
  theme_classic()+
  xlab("Z (scaled)") +
  ylab("") +
  geom_hline(yintercept = 0,lty="dashed")+
  ylim(-0.3,0.56)

i4 <- interplot(interaction_detection, var1 = "year_sc", var2 = "Z_sc") +
  theme_classic()+
  xlab("Z (scaled)") +
  ylab("") +
  geom_hline(yintercept = 0,lty="dashed")+
  ylim(-0.3,0.56)


panelplt_s2 <- (lg_p1 + ggtitle( 'a) Process') | lg_p2 + ggtitle( 'b) Confounding')  | lg_p3 + ggtitle( 'c) Compression')| lg_p4 + ggtitle( 'd) Detection\n bias')) / (i1 | i2 | i3| i4)

ggsave("~/figures/panel_s2.png",
       panelplt_s2, bg="white")

