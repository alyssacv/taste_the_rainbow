library(tidyverse)
library(janitor)
library(rio)
#library(readr)
#library(mgcv)
library(readxl)
library(survival)
library(survminer)
library(lubridate)
library(emmeans)
library(broom)
library(cowplot)
library(gridExtra)
library(patchwork)
library(rlang)
library(betareg)
library(statmod)



###############################temperature profiles################################# #####
profiles<-read_csv("data/temp_logs/apex_rainbow_profiles.csv")%>%
  mutate(profile=as.factor(profile))

x_limits <- c(0, 385)
y_limits <- c(27, 38.1)

plot_1a<-ggplot(profiles)+
  geom_line(aes(minute,setpoint,color=profile))+
  xlab("Time (Minutes)")+
  ylab("Temperature (°C)")+
  theme_classic(base_size=8)+
  scale_x_continuous(breaks=seq(0,385,60), limits = x_limits)+
  scale_y_continuous(breaks=seq(27,38,1), limits = y_limits)+
  scale_color_discrete(name="Treatment (°C)")+
  theme(legend.position = "right")+
  guides(color=guide_legend(reverse=TRUE));plot_1a


filenames <- list.files(path = "data/temp_logs", pattern = "\\d+\\.csv$",full.names=TRUE)  

actual_temperature <- lapply(filenames,function(i) read.csv(i)%>% mutate(file=i))%>% bind_rows()%>%
  clean_names()%>%
  mutate(treatment=parse_number(file))%>%
  select(-file)%>%
  mutate(date=force_tz(mdy_hm(date),"US/Hawaii"))%>%
  #filter(treatment==35)%>%
  filter(date>="2021-08-12 10:00",date<="2021-08-12 16:25")%>%
  mutate(treatment=as.factor(treatment))%>%
  mutate(start=force_tz(ymd_hm("2021-08-12 10:00"),"US/Hawaii"))%>%
  mutate(duration=as.numeric(difftime(date,start,units="mins")))


quartz(width = 7.2, height = 3.5)

plot_1b<-ggplot(actual_temperature)+
  geom_line(aes(duration,temp,color=treatment))+
  xlab("Time (Minutes)")+
  ylab("Temperature (°C)")+
  scale_color_discrete(name="Treatment (°C)")+
  theme_classic(base_size = 8)+
  theme(axis.ticks.y = element_line(), plot.margin = margin(0, 0, 0, 0, "cm"), legend.position = "right")+
  scale_x_continuous(breaks=seq(0,385,60),limits = x_limits)+
  scale_y_continuous(breaks=seq(27,38,1),limits = y_limits);plot_1b

legend<- get_legend(plot_1a + theme(legend.box.margin = margin(0, 0, 0, 0, "cm"), 
                                    legend.text = element_text(size = 7),  
                                    legend.title = element_text(size = 8),
                                    legend.key.size = unit(0.5, "cm")))

compare<-plot_grid(plot_1a+theme(legend.position="none", plot.margin = margin(0,0,0,0.5,"cm")), 
                   plot_1b+theme(legend.position = "none",axis.ticks.y=element_line(),
                                 plot.margin = margin(0, 0, 0, 0.5, "cm")),
                   nrow = 1, ncol=2, labels = c("A", "B"), hjust = c(-0.5,0.5), align='hv', scale=0.9)
plot_1c<-plot_grid(compare, legend, rel_widths = c(0.6, 0.1))
plot_1c


###############################larval survivorship########################### #####
larval_survivorship_24hrs<-read_xlsx("data/larvae_counts.xlsx", sheet = 2)%>%
  clean_names()%>%
  select(-c(x5,x6,x7))%>%
  pivot_longer(!treatment,names_to = "replicate",values_to = "count")%>%
  mutate(time=24)

larval_survivorship_0hrs<-read_xlsx("data/larvae_counts.xlsx", sheet = 1)%>%
  clean_names()%>%
  pivot_longer(!treatment,names_to = "replicate",values_to = "count")%>%
  mutate(time=0)

larval_survivorship<-bind_rows(larval_survivorship_0hrs, larval_survivorship_24hrs)%>%
  mutate(count=case_when(count>30~30,TRUE~count))%>%
  mutate(proportion = count/30)

table(larval_survivorship$count)

max(larval_survivorship$percentage)

quartz(width=5.5, height=3.5) 
ggplot(larval_survivorship)+
  geom_jitter(aes(treatment, proportion, color=as.factor(time)),height = 0, width = 0.4)+
  geom_smooth(aes(treatment, proportion, color=as.factor(time),fill=as.factor(time)),alpha=0.2)+ 
  scale_x_continuous(breaks=seq(27,38,1))+
  xlab("Treatment (°C)")+
  ylab("Porportion of Intact Larvae")+
  scale_color_discrete("Time (Hours)")+
  scale_fill_discrete("Time (Hours)")+
  theme_classic(base_size = 8)+
  theme(legend.position = c(0.141,0.3),legend.title = element_text (size=8))+
  annotate("text",27,0.6,label="Treatment p<0.001\nTime p<0.001\nTreatment*Time p<0.001",hjust=0, size=2.89)

larval_survivorship_model<-glm(proportion~as.factor(treatment)*as.factor(time),family = "quasibinomial", data = larval_survivorship)

summary(larval_survivorship_model)

car::Anova(larval_survivorship_model)

emm_larval_survivorship<- emmeans(larval_survivorship_model, ~ treatment | time)
contrast(emm_larval_survivorship,method="trt.vs.ctrl", ref=1)

### residuals ###

plot(larval_survivorship_model$residuals)

plot(larval_survivorship_model$fitted.values, larval_survivorship_model$residuals)
abline(h = 0, col = "red")

dispersion <- sum(residuals(larval_survivorship_model, type = "deviance")^2) / larval_survivorship_model$df.residual
dispersion


###############################larval settlement############################# #####
initial_settlement <- lapply(1:12, function(i) read_excel("data/settlement_counts.xlsx", sheet = i))%>%
  bind_rows()%>%
  clean_names()%>%
  rename(treatment=temp)%>%
  group_by(treatment)%>%
  mutate(average= mean(spat, na.rm = TRUE))%>%
  filter(plug!=156)

table(initial_settlement$treatment,initial_settlement$replicate)

levels(as.factor(initial_settlement$treatment))

str(initial_settlement)

x_positions_1 <- c(32)
x_positions_2 <- c(33,35,36)
x_positions_3 <- c(34,38)

quartz(w=5.5,h=3.5)

ggplot(initial_settlement)+
  geom_jitter(aes(treatment, spat, color=as.factor(treatment)),height = 0, width = 0.2)+
  geom_line(aes(treatment,average))+
  geom_hline(aes(yintercept = 1.041667), linetype="dashed",size=0.5)+
  geom_point(aes(treatment, average, color=as.factor(treatment)),size=4,pch=21,fill="white", stroke=1)+
  scale_x_continuous(breaks=seq(27,38,1))+
  xlab("Temperature (°C)")+
  ylab("Settled Juveniles (1 Week)")+
  scale_color_discrete("Treatment (°C)")+
  theme_classic(base_size = 8)+
  theme(legend.position = "none")+
  annotate("text",27,15,label="Treatment p<0.001",hjust=0, size=2.98)+
  map(x_positions_1, ~annotate("text", .x, -1, label="***", size=5))+
  map(x_positions_2, ~annotate("text", .x, -1, label="**", size=5))+
  map(x_positions_3, ~annotate("text", .x, -1, label="*", size=5))
 
  
initial_settlement_model<-MASS::glm.nb(spat~as.factor(treatment),data=initial_settlement)
summary(initial_settlement_model)

car::Anova(initial_settlement_model)

emm_treatment_ambient <- emmeans(initial_settlement_model, ~ treatment)
contrast(emm_treatment_ambient,"trt.vs.ctrl", ref=1)

### residuals ###

plot(initial_settlement_model$residuals)

plot(initial_settlement_model$fitted.values, initial_settlement_model$residuals)
abline(h = 0, col = "red")

dispersion <- sum(residuals(initial_settlement_model, type = "deviance")^2) / initial_settlement_model$df.residual
dispersion


###############################spat survivorship############################# #####
week_1_survivorship<- lapply(1:12, function(i) read_excel("data/week_1_survival_counts.xlsx", sheet = i)) %>% bind_rows() %>% clean_names()

week_2_survivorship<- lapply(1:12, function(i) read_excel("data/week_2_survival_counts.xlsx", sheet = i))%>% bind_rows()%>% clean_names()

week_3_survivorship<- lapply(1:12, function(i) read_excel("data/week_3_survival_counts.xlsx", sheet = i))%>% bind_rows()%>% clean_names()

week_5_survivorship<- lapply(1:12, function(i) read_excel("data/week_5_survival_counts.xlsx", sheet = i))%>% bind_rows()%>% clean_names()

week_24_survivorship<- lapply(1:12, function(i) read_excel("data/week_24_survival_counts.xlsx", sheet = i))%>% bind_rows()%>% clean_names()

cox_data<-bind_rows(week_1_survivorship,week_2_survivorship,week_3_survivorship,week_5_survivorship,week_24_survivorship)%>%
  ungroup()%>%
  select(-aggregate,-individual,-replicate)%>%
  rename(alive=dce)%>%
  rename(treatment=temp)%>%
  group_by(plug)%>%
  arrange(plug,week)%>%
  mutate(initial=case_when(week==1~alive))%>%
  fill(initial,.direction="down")%>%
  mutate(dead = initial-alive)%>%
  select(plug,treatment,week,alive,dead)%>%
  gather(status, count, -plug,-treatment,-week)%>%
  arrange(week)%>%
  mutate(count = map(count, ~ rep(1, .x))) %>%
  unnest()%>%
  select(-count)%>%
  mutate(status=case_when(status=="alive"~0,status=="dead"~1))%>%
  arrange(treatment)

final_cox<-coxph(Surv(week,status)~treatment,data = cox_data);final_cox

treatment_values<- unique(cox_data$treatment)

fit <- survfit(final_cox, newdata = data.frame(treatment = treatment_values))

summary(final_cox)
 
quartz(width=5.5, height=3.5)
plot<-ggsurvplot(fit, 
           conf.int = TRUE, 
           conf.int.size =0.3,
           legend.labs = paste("", treatment_values), 
           ggtheme = theme_minimal(), 
           data = cox_data,
           size = 0.5)

(plot$plot)+
  scale_fill_discrete(name = "Treatment (°C)")+
  scale_color_discrete(name="Treatment (°C)")+
  xlab("Weeks Post-Settlement")+
  ylab("Survival Probability (%)")+
  theme_classic(base_size = 8)+
  theme(legend.title = element_text(size = 8),legend.text=element_text(size=7),legend.position = c(0.3,0.2),
        legend.key.size =  unit(1, "cm"),legend.key.width = unit(0.3, "cm"), 
        legend.key.height = unit(0.3,"cm"), legend.box.spacing = unit(0.3, "cm"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  scale_y_continuous(breaks=seq(0,1.00,0.20))+
  annotate("text",20,0.7,label="Treatment p=0.103", size=2.89)
