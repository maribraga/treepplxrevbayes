### Compare outputs ###

#install_github("treeppl/treepplr")

library(tidyverse)
library(ape)
library(evolnets)
library(treepplr)
library(MCMCpack)
library(coda)
library(kdensity)
library(igraph)
library(patchwork)
library(jsonlite)
library(rjson)


###### REVBAYES #######

# read files
log1 <- read.table("output/out.1.compare.log", header = TRUE)
log2 <- read.table("output/out.2.compare.log", header = TRUE)

# take only columns of interest
chain1 <- log1[,c(1,5:6,8:11)]
chain2 <- log2[,c(1,5:6,8:11)]

# give them better column names
colnames(chain1) <- colnames(chain2) <- c("iteration","clock","beta", "gain_01", "loss_10", "gain_12", "loss_21")

its <- seq(10000,100000,100)

chain1 <- filter(chain1, iteration %in% its)
chain2 <- filter(chain2, iteration %in% its)

gelman.diag(mcmc.list(as.mcmc(chain1), as.mcmc(chain2)))

effectiveSize(chain1)
effectiveSize(chain2)


#### Parameter estimates ####

parameters <- chain1 %>%
  pivot_longer(cols = 2:7, names_to = "parameter", values_to = "value")

parameters %>%
  group_by(parameter) %>%
  summarise(mean = mean(value),
            sd = sd(value))

# plot_param <- ggplot(parameters, aes(parameter, value)) +
#   geom_violin(col = "grey40") +
#   stat_summary(fun.data = "median_hilow", color = "#E76F51", size = 0.5) +
#   theme_bw()
# plot_param

params <- colnames(trace[[1]])[c(5:6,8:11)]

paramdens_list <- list()
for(i in seq_along(params)){
  paramdens_list[[i]] <- plotTrace(trace[[i]], vars = params)
}

# Bayes factor

d_prior <- dexp(x=0, rate=1)

kd_beta <- kdensity(x = chain1$beta,
                    kernel='gamma',
                    support=c(0,Inf),
                    bw = 0.02)
max = kd_beta(0)

(BF <- d_prior/max)  # beta not different from zero


#### Character history ####

tree <- read_tree_from_revbayes("data/tree_Rev.tre")
host_tree <- read.tree("data/host_tree.tre")

matrix <- read.csv("data/matrix.csv", row.names = 1) %>% as.matrix()

history1 <- read_history("output/out.1.compare.history.txt", burnin = 0) %>%
  filter(iteration %in% its)

history2 <- read_history("output/out.2.compare.history.txt", burnin = 0) %>%
  filter(iteration %in% its)

# bug in revbayes - history contains all iterations

effective_rate(history1, tree)
effective_rate(history2, tree)

rate_gl(history1, tree)
rate_gl(history2, tree)
## evolnets TO DO add option for 3 states

mod <- mycomputeModules(matrix)

at_nodes1 <- posterior_at_nodes(history1, tree, host_tree, state = c(1,2))
at_nodes2 <- posterior_at_nodes(history1, tree, host_tree, state = c(1,2))

p_asr12 <- plot_matrix_phylo(matrix, at_nodes1, tree, host_tree, modules = mod, threshold = 0.9)
p_asr22 <- plot_matrix_phylo(matrix, at_nodes2, tree, host_tree, modules = mod, threshold = 0.9)
p_asr11 <- plot_matrix_phylo(matrix, at_nodes1, tree, host_tree, modules = mod, threshold = 0.9, type = "repertoires", repertoire = "fundamental")
p_asr21 <- plot_matrix_phylo(matrix, at_nodes2, tree, host_tree, modules = mod, threshold = 0.9, type = "repertoires", repertoire = "fundamental")


plot(tree, show.node.label = TRUE)
axisPhylo()

pp11 <- at_nodes1$post_states[,,1] %>%
  as.data.frame() %>%
  rownames_to_column(var = "node") %>%
  pivot_longer(2:4, names_to = "host", values_to = "pp11")

pp12 <- at_nodes1$post_states[,,2] %>%
  as.data.frame() %>%
  rownames_to_column(var = "node") %>%
  pivot_longer(2:4, names_to = "host", values_to = "pp12")

pp21 <- at_nodes2$post_states[,,1] %>%
  as.data.frame() %>%
  rownames_to_column(var = "node") %>%
  pivot_longer(2:4, names_to = "host", values_to = "pp21")

pp22 <- at_nodes2$post_states[,,2] %>%
  as.data.frame() %>%
  rownames_to_column(var = "node") %>%
  pivot_longer(2:4, names_to = "host", values_to = "pp22")

pps <- pp11 %>% left_join(pp12) %>% left_join(pp21) %>% left_join(pp22) %>%
  unite("pair", 1:2, remove = TRUE) %>%
  pivot_longer(2:5, names_to = "chain", values_to = "pp")

ggplot(pps, aes(chain, pair)) +
  geom_raster(aes(alpha = pp))  +
  theme_bw() +
  theme(axis.text.y = element_blank())




#### Likelihoods ####

liks1 <- log1[,c(1:4)]
liks2 <- log2[,c(1:4)]

mean(liks1$Likelihood)
mean(norm_weights(liks1$Likelihood))

ggplot() +
  geom_density(aes(Posterior), data = liks1, fill = "red", col = "red", alpha = 0.6) +
  geom_density(aes(Posterior), data = liks2, fill = "darkred", col = "darkred", alpha = 0.6) +
  geom_density(aes(Likelihood), data = liks1, fill = "blue", col = "blue", alpha = 0.6) +
  geom_density(aes(Likelihood), data = liks2, fill = "darkblue", col = "darkblue", alpha = 0.6) +
  theme_bw()






#### TREEPPL ####

## 10k particles, 2 sweeps, twice

out10k_2s_1 <- separate_outs("output/out10k2.json") %>%
  lapply(tidy_samples)

out10k_2s_2 <- separate_outs("output/out10k2_2.json") %>%
  lapply(tidy_samples)

out20k_1 <- bind_rows(out10k_2s_1) #, .id = "20k_1"
out20k_2 <- bind_rows(out10k_2s_2) # , .id = "20k_2"

out40k <- bind_rows(out20k_1, out20k_2)

## 1k particles, 20 sweeps

out1k_20s <- separate_outs("output/out1k20.json") %>%
  lapply(tidy_samples)

out20k_3 <- bind_rows(out1k_20s)

out60k <- bind_rows(out40k, out20k_3)

## 2k particles, 20 sweeps

out2k_20s <- separate_outs("output/out2k20.json") %>%
  lapply(tidy_samples)

out2k_20s <- read_treeppl_output("output/out2k20.json")

out40k_1 <- bind_rows(out2k_20s)

out100k <- bind_rows(out40k_1, out60k)

## Renormalize weights ##

out100k <- out100k %>%
  mutate(re_norm_w = nweights/sum(nweights))

# to test in Tracer
write.table(out100k[1:2000,], "output/treeppl_logs.log",
            row.names = FALSE, quote = FALSE, sep = "\t")

## Posterior ##

ggplot(out100k) +
  geom_histogram(aes(weights)) +
  theme_bw()

ggplot(out100k) +
  geom_histogram(aes(re_norm_w)) +
  theme_bw()


#### _subsampling ####

out2k_20s[[1]] %>%
  arrange(desc(norm_weights)) %>%
  as_tibble()

out2k_20s[[2]] %>%
  arrange(desc(norm_weights)) %>%
  as_tibble()

out2k_20s[[3]] %>%
  arrange(desc(norm_weights)) %>%
  as_tibble()
















# revbayes x treeppl

g1 <- ggplot() +
  geom_density(aes(clock, after_stat(scaled)), fill = "steelblue", col = NA, alpha = 0.7, data = chain1) +
  geom_col(aes(clock, nweights), width = 0.03, fill = "darkorange", data = out100k) +
  theme_bw()

g2 <- ggplot(out100k) +
  geom_density(aes(beta, after_stat(scaled)), fill = "steelblue", col = NA, alpha = 0.7, data = chain1) +
  geom_col(aes(beta, nweights), width = 0.03, fill = "darkorange", data = out100k) +
  theme_bw()

g3 <- ggplot(out100k) +
  geom_density(aes(gain_01, after_stat(scaled)), fill = "steelblue", col = NA, alpha = 0.7, data = chain1) +
  geom_col(aes(gain_01, nweights), width = 0.003, fill = "darkorange", data = out100k) +
  theme_bw()

g4 <- ggplot(out100k) +
  geom_density(aes(loss_10, after_stat(scaled)), fill = "steelblue", col = NA, alpha = 0.7, data = chain1) +
  geom_col(aes(loss_10, nweights), width = 0.003, fill = "darkorange", data = out100k) +
  theme_bw()

g5 <- ggplot(out100k) +
  geom_density(aes(gain_12, after_stat(scaled)), fill = "steelblue", col = NA, alpha = 0.7, data = chain1) +
  geom_col(aes(gain_12, nweights), width = 0.003, fill = "darkorange", data = out100k) +
  theme_bw()

g6 <- ggplot(out100k) +
  geom_density(aes(loss_21, after_stat(scaled)), fill = "steelblue", col = NA, alpha = 0.5, data = chain1) +
  geom_col(aes(loss_21, nweights), width = 0.003, fill = "darkorange", data = out100k) +
  theme_bw()

g1 + g2 + g3 + g4 + g5 + g6 + plot_layout(ncol = 2, nrow = 3, byrow = TRUE)


wmeans <- apply(out100k[,2:7], 2, function(x, a1) weighted.mean(x, a1), a1 = out100k$nweights)
wmeans

wmeans_pruned <- apply(out100k_pruned[,2:7], 2, function(x, a1) weighted.mean(x, a1), a1 = out100k_pruned$nweights)
wmeans_pruned

apply(chain1[,2:7], 2, mean)
apply(chain2[,2:7], 2, mean)



#### Likelihoods ####

# remove the samples with too high nweight

highw <- out100k %>%
  filter(nweights > 0.2) %>%
  arrange(desc(nweights))

out100k_pruned <- out100k %>%
  filter(nweights < 0.5)

ggplot(out100k_pruned) +
  geom_col(aes(clock, nweights), width = 0.03, fill = "red") +
  theme_bw()


out2 <- read_treeppl_output("/Users/mari/repos/treeppl_playground/treeppl_revbayes/output/out1k2.json")







