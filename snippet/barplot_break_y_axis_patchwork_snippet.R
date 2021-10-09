
d1 <- d.cognitive_performance[, -1]
ht1 <- ht.cognitive_performance
group <- factor(d.cognitive_performance[, 1])


## calculate group mean, lower and upper
## d1s is array of 3 dimensions

d.summary <- lapply(colnames(d1), function(x, g = group) {
  out <- vaggregate(d1[, x], unname(g), mean_cl_normal, mult = 1)
  
  out <- matrix(unlist(out), nrow = nrow(out), ncol = ncol(out)) %>% t
  data.frame(
    Group = levels(g),
    variable = x,
    y = out[,1],
    ymin = out[, 2], 
    ymax = out[, 3]
  )
}) %>% do.call(rbind, .)

## multcompLetters work on lower.tri
d.sig_letters <- lapply(1:nrow(ht1), function(i, g = group) {
  lipid.name <- rownames(ht1)[i]
  
  p_matrix <- matrix(
    1, nrow = nlevels(g), ncol = nlevels(g),
    dimnames = list(levels(g), levels(g))    
  )
  
  p_matrix[lower.tri(p_matrix)] <- ht1[i, stringr::str_detect(colnames(ht1), 'Dunn')]
  
  sig_letters <- multcompView::multcompLetters(
    p_matrix, threshold = 0.05
  )$Letters
  
  out <- data.frame(Group = names(sig_letters), variable = lipid.name, sig_letters = sig_letters)
  rownames(out) <- NULL
  out
}) %>% do.call(rbind, .)

testthat::expect_equal(
  d.summary$Group, d.sig_letters$Group
)

testthat::expect_equal(
  d.summary$variable, d.sig_letters$variable
)

df0 <- merge(d.summary, d.sig_letters, by = c('Group', 'variable'))

p0 <- ggplot(df0, aes(x = variable, y = y, fill = Group)) +
  geom_bar(stat = 'identity', position = position_dodge2()) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax),
                width = 0.2,
                position = position_dodge(width = 0.9)) +
  geom_text(aes(y = 1.1 * ymax + 2, label = sig_letters),
            size = 3,
            position = position_dodge(width = 0.9)) +
  labs(x = '', y = '') +
  theme_classic()

p1 <- p0 + coord_cartesian(ylim = c(0, 20))

p2 <- p0 + coord_cartesian(ylim = c(20, 50)) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0), add = c(0, 0))) +
  theme(
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p3 <- p0 + coord_cartesian(ylim = c(50, 280)) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0), add = c(0, 0))) +
  theme(
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p3 / p2 / p1 + plot_layout(guides = 'collect')
