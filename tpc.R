library(tidyverse)
library(slider)
library(janitor)
library(nls.multstart)
library(broom)
library(minpack.lm)
library(ggridges)


install.packages("devtools")
library(devtools)
devtools::install_github("padpadpadpad/rTPC")
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
    arrange(aic) %>%
    distinct()
}
#best fit is taylorsexton_1972 based on minimum AIC

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
