# example of BLUPS

# getting blups (best linear predictors from the model)
blups <- ranef(RE0.reml) t.cultivar <- blups$cultivar
t.cultivar <- rownames_to_column(t.cultivar, var = "cultivar")
# this is deviation from the meta-analytic mean -ve values
# mean more pollinator dependence and +ve means less than the
# average cultivar
colnames(t.cultivar) <- c("Cultivar", "Deviation", "SE", "Lower_bound", "Upper_bound")
# getting cutliver simple
cult <- data1 %>% distinct(cultivar, .keep_all = TRUE) %>% select(cultivar, cultivar.simple) %>% arrange(cultivar)
# knitr::kable(t.cultivar) col.names = c(’Cultivar’,
# ’Deviation’, ’SE’, ’Lower bound’, ’Upper bound’))
cult.mean <- RE0.reml$b + t.cultivar$Deviation 
cult.se <- sqrt(RE0.reml$se^2 + t.cultivar$SE^2) 
cult.lb <- cult.mean - cult.se * qnorm(0.975) 
cult.ub <- cult.mean + cult.se * qnorm(0.975)
t.cultivar2 <- tibble(Cultivar = t.cultivar$Cultivar, Type = cult$cultivar.simple, Mean = cult.mean, SE = cult.se, Lower_bound = cult.lb, Upper_bound = cult.ub) %>% 
  arrange(Type)


png(file="figs/figS3 cultivars.png",
    width=180,height=280,res=600,units="mm")

cul.plot <- ggplot(data = t.cultivar2, aes(x = Mean, y = Cultivar)) +
  geom_errorbarh(aes(xmin = Lower_bound, xmax = Upper_bound),  
                 height = 0, show.legend = FALSE, size = 0.5, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.5) +
  geom_vline(xintercept = RE0.reml$b, linetype = 1, colour = "red") +
  geom_point(aes(fill = Type), size = 3, shape = 21) + 
  xlim(-1.6, 1) +
  theme_bw() +
  labs(x = "lnRR (effect size): cultivar-specific means", y = "") + 
  theme(legend.position= c(0.99, 0.01), legend.justification = c(1, 0)) +
  theme(legend.title = element_text(size = 9)) +
  theme(legend.direction="horizontal") +
  theme(axis.text.y = element_blank()) +
  theme(axis.text.y = element_text(size = 10, colour ="black",
                                   hjust = 0.5))
cul.plot
dev.off()