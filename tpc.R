library(tidyverse)
library(slider)
library(janitor)
library(nls.multstart)
library(broom)
library(minpack.lm)
library(ggridges)
library(gt)
library(gridExtra)
library(grid)
library(flextable)
library(officer)

#RStudio Version 2026.01.2+418 (2026.01.2+418)
#to use rTPC with all equations we need to install from github. first install pak (agree to install the "remotes" dependency), then install rTPC via github.
install.packages("pak")
library("pak")
pak::pak("padpadpadpad/rTPC")
library(rTPC)

############# FORMAT INPUT DATA ################################# #####

data<-readRDS("data/settlement")%>%
  select(treatment,spat)%>%
  rename(temp=treatment,rate=spat)%>%
  mutate(rate=rate+0.001)

get_model_names() 
output_overall<-data.frame(model=character(),aic=numeric())

for (i in c(1:49)){
  mod<-get_model_names()[i]
  params<-  setdiff(names(formals(get(mod, asNamespace("rTPC")))), "temp");params
  start_vals <- get_start_vals(data$temp, data$rate, model_name = mod)
  low_lims <- get_lower_lims(data$temp, data$rate, model_name = mod)
  upper_lims <- get_upper_lims(data$temp, data$rate, model_name = mod)
  formula <- as.formula(paste("rate ~",sprintf("%s(temp = temp, %s)",mod,paste(params, collapse = ", "))));formula
  fit <- try(nls_multstart(formula,
                           data = data,
                           iter = 500,
                           start_lower = start_vals - 10,
                           start_upper = start_vals + 10,
                           lower = low_lims,
                           upper = upper_lims,
                           supp_errors = 'Y'),    silent = TRUE)
  if (inherits(fit, "try-error") || is.null(fit)) {
    next
  }

  output_overall <- output_overall %>%
    bind_rows(tibble(model = mod, aic = AIC(fit))) %>%
    distinct(model,aic,.keep_all=TRUE)%>%
    arrange(aic)
    
}

output_overall
#best fit is taylorsexton_1972 based on minimum AIC

###model selection AIC table###

output_overall_table<-output_overall%>%
  rename(`Model`=model,
         `AIC`=aic)

output_overall_table_word<-flextable(output_overall_table)%>%
  theme_booktabs()%>%
  colformat_double(j="AIC",digits=2)%>%
  flextable::font(fontname="Times New Roman",part="all")%>%
  fontsize(size=12,part="all")%>%
  set_caption("Table S10. Thermal Performance Model selection using AIC comparisons.")%>%
  bold(part="header")%>%
  align(align="center",part="all")%>%
  align(j="Model",align="left",part="all")%>%
  autofit()

output_overall_table_word

save_as_docx(output_overall_table_word,
             path="output/TPC_model_comparison_table.docx")

############################ APPLY ########################## #####
get_model_names() #taylorsexton_1972
mod<-get_model_names()[44]
params<-  setdiff(names(formals(get(mod, asNamespace("rTPC")))), "temp");params

start_vals <- get_start_vals(data$temp, data$rate, model_name = mod)
low_lims <- get_lower_lims(data$temp, data$rate, model_name = mod)
upper_lims <- get_upper_lims(data$temp, data$rate, model_name = mod)
formula <- as.formula(paste("rate ~",sprintf("%s(temp = temp, %s)",mod,paste(params, collapse = ", "))));formula

global_fit <- nls_multstart(formula,
                            data = data,
                            iter = 500,
                            start_lower = start_vals - 10,
                            start_upper = start_vals + 10,
                            lower = low_lims,
                            upper = upper_lims,
                            supp_errors = 'Y')

new_data <- data.frame(temp = seq(min(data$temp), max(data$temp), 0.1))
predicted <- augment(global_fit, newdata = new_data)%>%rename(rate=2)

calc_params(global_fit) #thermal optimum at 33.778 degrees

###global fit table###
global_fit_params<-calc_params(global_fit)%>%
  pivot_longer(cols=everything(),
               names_to="parameter",
               values_to="estimate")%>%
  mutate(parameter=recode(parameter,
                          rmax="Maximum Rate",
                          topt="Topt (°C)",
                          ctmin="CTmin (°C)",
                          ctmax="CTmax (°C)",
                          e="Activation Energy",
                          eh="Deactivation Energy",
                          q10="Q10",
                          thermal_safety_margin="Thermal Safety Margin (°C)",
                          thermal_tolerance="Thermal Tolerance (°C)",
                          breadth="Thermal Breadth (°C)",
                          skewness="Skewness"))%>%
  rename(`Parameter`=parameter,
         `Estimate`=estimate)

global_fit_params_word<-flextable(global_fit_params)%>%
  theme_booktabs()%>%
  colformat_double(j="Estimate",digits=3)%>%
  flextable::font(fontname="Times New Roman",part="all")%>%
  fontsize(size=12,part="all")%>%
  set_caption("Table S11. Estimated Thermal Performance Curve parameters.")%>%
  bold(part="header")%>%
  align(align="center",part="all")%>%
  align(j="Parameter",align="left",part="all")%>%
  autofit()

global_fit_params_word

save_as_docx(global_fit_params_word, path="output/global_fit_parameters_table.docx")

###settlement and tcp figure###
ggplot(data)+
  geom_point(aes(temp,rate))+
  geom_smooth(aes(temp,rate),data=predicted, se = TRUE)

topt <- 33.778

topt_y <- predicted %>%
  filter(abs(temp - topt) == min(abs(temp - topt))) %>%
  pull(rate)

x_positions_1 <- c(32)
x_positions_2 <- c(33,35,36)
x_positions_3 <- c(34,38)

quartz(w=5.5,h=3.5)

ggplot(data) +
  geom_jitter(aes(temp, rate, color = as.factor(temp)), height = 0, width = 0.2) +
  geom_line(aes(temp, rate), data = predicted, linewidth = 1) +
  geom_hline(aes(yintercept = 1.041667), linetype="dashed", size=0.5)+
  scale_x_continuous(breaks=seq(27,38,1)) +
  xlab("Temperature (°C)") +
  ylab("Settled Juveniles (1 Week)") +
  theme_classic(base_size = 8) +
  theme(legend.position = "none") +
  annotate("text",27,15,label="Treatment p<0.001",hjust=0, size=2.98) +
  annotate("text",27,14,label="Thermal Optimum 33.778°C",hjust=0, size=2.98) +
  annotate("segment", x = topt, y = topt_y + 3, xend = topt, yend = 3.8, linewidth = 0.4,
           arrow = arrow(length = unit(0.18, "cm"))) +
  annotate("text", x = topt + - 0.2, y = topt_y + 3.5, label = "italic(T[opt])", parse = TRUE,
           hjust = 0, size = 3.2) +
  map(x_positions_1, ~annotate("text", .x, -1, label="***", size=5))+
  map(x_positions_2, ~annotate("text", .x, -1, label="**", size=5))+
  map(x_positions_3, ~annotate("text", .x, -1, label="*", size=5))
