library(lme4)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)


#### Simulating two groups of people ####
n_in_group = 20e3
group2_eff = -1

parity0_eff = -0.5
parity1_eff = 0
parity2_eff = 0.3
parity3_eff = 0.5




#### Mothers' random effects/averages ####
mu1 = rnorm(n_in_group, 0, 1)
mu2 = rnorm(n_in_group, 0, 1) + group2_eff

df = bind_rows(data.frame("mu"=mu1, "g"="1"), data.frame("mu"=mu2, "g"="2"))
df$maternal_id = 1:nrow(df)




#### Outcomes of each pregnancy ####
# i.e. mother's mean + beta_parity + gaussian noise
df$y0 = df$mu + parity0_eff + rnorm(nrow(df))
df$y1 = df$mu + parity1_eff + rnorm(nrow(df))
df$y2 = df$mu + parity2_eff + rnorm(nrow(df))
df$y3 = df$mu + parity3_eff + rnorm(nrow(df))
df = gather(df, key="pa", value="y", y0:y3)

# Store the no-group-effect value 
df$y_nogreff = df$y + 1*(df$g=="2")

# Same data, but assuming mothers in group1 stop after the first kid, group2 after 2 kids:
df2 = filter(df, g!="1" | pa %in% c("y0", "y1"))
group_by(df2, g, pa) %>%
  summarize(mean(y), mean(y_nogreff), n())




#### Run the models ####
# If everyone has equal number of kids:
coefs_eqkids_lin = coef(lm(y ~ pa, data=group_by(df, maternal_id) %>% sample_n(1)))
coefs_eqkids_mix = fixef(lmer(y ~ pa + (1 | maternal_id), data=df))
# approx equal

# If number of kids differs by groups, but no group effect:
coefs_diffkids_lin = coef(lm(y_nogreff ~ pa, data=group_by(df2, maternal_id) %>% sample_n(1)))
coefs_diffkids_mix = fixef(lmer(y_nogreff ~ pa + (1 | maternal_id), data=df2))
# approx equal

# If number of kids differs by groups, and group effect:
coefs_diffkids_greff_lin = coef(lm(y ~ pa, data=group_by(df2, maternal_id) %>% sample_n(1)))
coefs_diffkids_greff_mix = fixef(lmer(y ~ pa + (1 | maternal_id), data=df2))
# NOT equal: LIN bad, MIX better, although still underestimated

results = rbind(coefs_eqkids_lin,coefs_eqkids_mix, coefs_diffkids_lin,
coefs_diffkids_mix, coefs_diffkids_greff_lin, coefs_diffkids_greff_mix)




#### Saving ####
fwrite(results, snakemake@output[[1]],sep=",")

cat("Results from simulation: \n") 
print(results)

