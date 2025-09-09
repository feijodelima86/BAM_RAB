library("vegan")
library("RVAideMemoire")

dis <- vegdist(my.data, method="euclidean")

# Permanova using adonis

results <- adonis2(formula = dis ~ groups,  permutations = 999)

pairwise.permanova <- pairwise.perm.manova(dis, groups, nperm = 999, p.method = "fdr", F = FALSE, R2 = FALSE)

# Betadisper, different approach that touches on beta diversity.

mod <- betadisper(dis, groups, sqrt.dist=T, bias.adjust=T, type = "centroid")

pairwise.variance <- permutest(mod, pairwise = TRUE, permutations = 999)
