########################################## Setup ###################################################

kenBlue = "#3A8EE4"
A = LIPS::highway$A
y = LIPS::highway$y
lambda = LIPS::highway$lambda
M = LIPS::highway$MarkovBasis
L = LIPS::highway$LatticeBasis

a = c()
for(i in 8:28){
  a[i-7] = lpSolve::lp(direction = "max",
                       objective.in = diag(28)[i,],
                       const.mat = A,
                       const.dir = "==",
                       const.rhs = y,
                       all.int = TRUE)$solution[i]
}

b = sum(lpSolve::lp(direction = "max",
                    objective.in = c(rep(0,7), rep(1,21)),
                    const.mat = A,
                    const.dir = "==",
                    const.rhs = y,
                    all.int = TRUE)$solution[8:28])

########################################## Samplers ################################################
##########################################          Uniform Model ##################################

# FMB
set.seed(2026)
x_FMB = LIPS::xSampler(A = A,
                       y = y,
                       B = M,
                       nSample = 1e+06,
                       nBurnin = 1,
                       nChains = 1,
                       thinning = 10)

set.seed(2026)
x_1MB = LIPS::extensionSampler(A = A,
                               y = y,
                               B = L,
                               p = 1,
                               loggamma = 1,
                               nSample = 1e+06,
                               nBurnin = 1,
                               nChains = 1,
                               thinning = 10)

# 1MB
set.seed(2026)
x_1MB_w = LIPS::extensionSampler(A = A,
                                 y = y,
                                 B = L,
                                 p = 1,
                                 loggamma = 4,
                                 w = 3,
                                 nSample = 1e+06,
                                 nBurnin = 1,
                                 nChains = 1,
                                 thinning = 10)

# HR
set.seed(2026)
x_hyperrectangle = LIPS::extensionSampler(A = A,
                                          y = y,
                                          B = L,
                                          extension = "hyperrectangle",
                                          A2Idx = 8:ncol(A),
                                          a = a,
                                          loggamma = -13,
                                          nSample = 1e+06,
                                          nBurnin = 1,
                                          nChains = 1,
                                          thinning = 10)

set.seed(2026)
x_hyperrectangle_w = LIPS::extensionSampler(A = A,
                                            y = y,
                                            B = L,
                                            extension = "hyperrectangle",
                                            A2Idx = 8:ncol(A),
                                            a = a,
                                            w = 0.5,
                                            loggamma = 2,
                                            nSample = 1e+06,
                                            nBurnin = 1,
                                            nChains = 1,
                                            thinning = 10)

# KS
set.seed(2026)
x_knapsack = LIPS::extensionSampler(A = A,
                                    y = y,
                                    B = L,
                                    extension = "knapsack",
                                    A2Idx = 8:ncol(A),
                                    a = rep(1,21),
                                    b = b,
                                    loggamma = -13.5,
                                    nSample = 1e+06,
                                    nBurnin = 1,
                                    nChains = 1,
                                    thinning = 10)

set.seed(2026)
x_knapsack_w = LIPS::extensionSampler(A = A,
                                      y = y,
                                      B = L,
                                      extension = "knapsack",
                                      A2Idx = 8:ncol(A),
                                      a = rep(1,21),
                                      b = b,
                                      loggamma = 3,
                                      w = 1,
                                      nSample = 1e+06,
                                      nBurnin = 1,
                                      nChains = 1,
                                      thinning = 10)

any(sapply(x_1MB, function(chain) any(apply(chain, 1, function(row) any(row < 0)))))
any(sapply(x_1MB, function(chain) any(apply(chain, 1, function(row) all(row >= 0)))))

any(sapply(x_1MB_w, function(chain) any(apply(chain, 1, function(row) any(row < 0)))))
any(sapply(x_1MB_w, function(chain) any(apply(chain, 1, function(row) all(row >= 0)))))

any(sapply(x_hyperrectangle, function(chain) any(apply(chain, 1, function(row) any(row < 0)))))
any(sapply(x_hyperrectangle, function(chain) any(apply(chain, 1, function(row) all(row >= 0)))))

any(sapply(x_hyperrectangle_w, function(chain) any(apply(chain, 1, function(row) any(row < 0)))))
any(sapply(x_hyperrectangle_w, function(chain) any(apply(chain, 1, function(row) all(row >= 0)))))

any(sapply(x_knapsack, function(chain) any(apply(chain, 1, function(row) any(row < 0)))))
any(sapply(x_knapsack, function(chain) any(apply(chain, 1, function(row) all(row >= 0)))))

any(sapply(x_knapsack_w, function(chain) any(apply(chain, 1, function(row) any(row < 0)))))
any(sapply(x_knapsack_w, function(chain) any(apply(chain, 1, function(row) all(row >= 0)))))

x_1MB_pos = LIPS::only_positive(x_1MB)
x_1MB_w_pos = LIPS::only_positive(x_1MB_w)
x_hyperrectangle_pos = LIPS::only_positive(x_hyperrectangle)
x_hyperrectangle_w_pos = LIPS::only_positive(x_hyperrectangle_w)
x_knapsack_pos = LIPS::only_positive(x_knapsack)
x_knapsack_w_pos = LIPS::only_positive(x_knapsack_w)

##########################################          Poisson Model ##################################

set.seed(2026)
x_FMB_pois = LIPS::xSampler(A = A,
                            y = y,
                            B = M,
                            Model = "Pois",
                            lambda = lambda,
                            nSample = 1e+06,
                            nBurnin = 1,
                            nChains = 1,
                            thinning = 10)

set.seed(2026)
x_1MB_pois = LIPS::extensionSampler(A = A,
                                    y = y,
                                    B = L,
                                    Model = "poisson",
                                    lambda = lambda,
                                    p = 1,
                                    loggamma = -25.5,
                                    nSample = 1e+06,
                                    nBurnin = 1,
                                    nChains = 1,
                                    thinning = 10)

set.seed(2026)
x_1MB_w_pois = LIPS::extensionSampler(A = A,
                                      y = y,
                                      B = L,
                                      Model = "poisson",
                                      lambda = lambda,
                                      p = 1,
                                      loggamma = 0,
                                      w = 3,
                                      nSample = 1e+06,
                                      nBurnin = 1,
                                      nChains = 1,
                                      thinning = 10)

set.seed(2026)
x_hyperrectangle_pois = LIPS::extensionSampler(A = A,
                                               y = y,
                                               B = L,
                                               Model = "poisson",
                                               lambda = lambda,
                                               extension = "hyperrectangle",
                                               a = a,
                                               A2Idx = 8:28,
                                               loggamma = -3,
                                               nSample = 1e+06,
                                               nBurnin = 1,
                                               nChains = 1,
                                               thinning = 10)

set.seed(2026)
x_hyperrectangle_w_pois = LIPS::extensionSampler(A = A,
                                                 y = y,
                                                 B = L,
                                                 Model = "poisson",
                                                 lambda = lambda,
                                                 extension = "hyperrectangle",
                                                 a = a,
                                                 A2Idx = 8:28,
                                                 loggamma = -1.5,
                                                 w = 0.5,
                                                 nSample = 1e+06,
                                                 nBurnin = 1,
                                                 nChains = 1,
                                                 thinning = 10)

set.seed(2026)
x_knapsack_pois = LIPS::extensionSampler(A = A,
                                         y = y,
                                         B = L,
                                         Model = "poisson",
                                         lambda = lambda,
                                         extension = "knapsack",
                                         a = rep(1, 21),
                                         b = b,
                                         A2Idx = 8:28,
                                         loggamma = -3.5,
                                         nSample = 1e+06,
                                         nBurnin = 1,
                                         nChains = 1,
                                         thinning = 10)

set.seed(2026)
x_knapsack_w_pois = LIPS::extensionSampler(A = A,
                                           y = y,
                                           B = L,
                                           Model = "poisson",
                                           lambda = lambda,
                                           extension = "knapsack",
                                           a = rep(1, 21),
                                           b = b,
                                           A2Idx = 8:28,
                                           loggamma = -1.5,
                                           w = 0.5,
                                           nSample = 1e+06,
                                           nBurnin = 1,
                                           nChains = 1,
                                           thinning = 10)

any(sapply(x_1MB_pois, function(chain) any(apply(chain, 1, function(row) any(row < 0)))))
any(sapply(x_1MB_pois, function(chain) any(apply(chain, 1, function(row) all(row >= 0)))))

any(sapply(x_1MB_w_pois, function(chain) any(apply(chain, 1, function(row) any(row < 0)))))
any(sapply(x_1MB_w_pois, function(chain) any(apply(chain, 1, function(row) all(row >= 0)))))

any(sapply(x_hyperrectangle_pois, function(chain) any(apply(chain, 1, function(row) any(row < 0)))))
any(sapply(x_hyperrectangle_pois, function(chain) any(apply(chain, 1, function(row) all(row >= 0)))))

any(sapply(x_hyperrectangle_w_pois, function(chain) any(apply(chain, 1, function(row) any(row < 0)))))
any(sapply(x_hyperrectangle_w_pois, function(chain) any(apply(chain, 1, function(row) all(row >= 0)))))

any(sapply(x_knapsack_pois, function(chain) any(apply(chain, 1, function(row) any(row < 0)))))
any(sapply(x_knapsack_pois, function(chain) any(apply(chain, 1, function(row) all(row >= 0)))))

any(sapply(x_knapsack_w_pois, function(chain) any(apply(chain, 1, function(row) any(row < 0)))))
any(sapply(x_knapsack_w_pois, function(chain) any(apply(chain, 1, function(row) all(row >= 0)))))

x_1MB_pois_pos = LIPS::only_positive(x_1MB_pois)
x_1MB_w_pois_pos = LIPS::only_positive(x_1MB_w_pois)

x_hyperrectangle_pois_pos = LIPS::only_positive(x_hyperrectangle_pois)
x_hyperrectangle_w_pois_pos = LIPS::only_positive(x_hyperrectangle_w_pois)

x_knapsack_pois_pos = LIPS::only_positive(x_knapsack_pois)
x_knapsack_w_pois_pos = LIPS::only_positive(x_knapsack_w_pois)

########################################## Plots ###################################################
##########################################       Uniform model #####################################
##########################################                     Positive AND Negative ###############

pdf("highway_unif.pdf", width = 15, height = 13)
layout_matrix <- rbind(
  c(0,  0,  1,  2),
  c(0,  3,  7,  8),
  c(6, 4, 9, 10),
  c(0, 5, 11, 12)
)

layout(
  layout_matrix,
  widths  = c(1, 0.25, 1, 1),
  heights = c(0.30, 1, 1, 1)
)

par(mar = c(3.5, 3.5, 2, 1), mgp = c(1.6, 0.5, 0))

# Column title boxes

plot.new()
text(0.5, 0.5, "w=0", cex = 1.5, font = 2)

plot.new()
text(0.5, 0.5, "w \u2260 0", cex = 1.5, font = 2)

# Row title boxes
plot.new()
text(0.5, 0.5, "1MB-Extension", srt = 90, cex = 1.5, font = 2)

plot.new()
text(0.5, 0.5, "Hyperrectangle", srt = 90, cex = 1.5, font = 2)

plot.new()
text(0.5, 0.5, "Knapsack", srt = 90, cex = 1.5, font = 2)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(x_FMB[,6]),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5,
     main = "FMB")
lines(x_FMB[1:1e+05,6], col = kenBlue)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(do.call(rbind,x_1MB)[,6]),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5)
lines(x_1MB[[1]][,6], col = kenBlue)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(do.call(rbind,x_1MB_w)[,6]),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5)
lines(x_1MB_w[[1]][,6], col = kenBlue)

plot(NA, xlim = c(0,1e+05),
     ylim = range(do.call(rbind,x_hyperrectangle)[,6]),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5)
lines(x_hyperrectangle[[1]][,6], col = kenBlue)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(do.call(rbind,x_hyperrectangle_w)[,6]),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5)
lines(x_hyperrectangle_w[[1]][,6], col = kenBlue)

plot(NA, xlim = c(0,1e+05),
     ylim = range(do.call(rbind,x_knapsack)[,6]),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5)
lines(x_knapsack[[1]][,6], col = kenBlue)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(do.call(rbind,x_knapsack_w)[,6]),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5)
lines(x_knapsack_w[[1]][,6], col = kenBlue)
dev.off()

##########################################                     Only positive #######################

pdf("highway_unif_pos.pdf", width = 15, height = 13)
layout_matrix <- rbind(
  c(0,  0,  1,  2),
  c(0,  3,  7,  8),
  c(6, 4, 9, 10),
  c(0, 5, 11, 12)
)

layout(
  layout_matrix,
  widths  = c(1, 0.25, 1, 1),
  heights = c(0.30, 1, 1, 1)
)

par(mar = c(3.5, 3.5, 2, 1), mgp = c(1.6, 0.5, 0))

# Column title boxes

plot.new()
text(0.5, 0.5, "w=0", cex = 1.5, font = 2)

plot.new()
text(0.5, 0.5, "w \u2260 0", cex = 1.5, font = 2)

# Row title boxes
plot.new()
text(0.5, 0.5, "1MB-Extension", srt = 90, cex = 1.5, font = 2)

plot.new()
text(0.5, 0.5, "Hyperrectangle", srt = 90, cex = 1.5, font = 2)

plot.new()
text(0.5, 0.5, "Knapsack", srt = 90, cex = 1.5, font = 2)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(x_FMB[,6]),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5,
     main = "FMB")
lines(x_FMB[1:1e+05,6], col = kenBlue)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(do.call(rbind,x_1MB)[,6]),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5)
lines(x_1MB_pos[[1]][,6], col = kenBlue)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(do.call(rbind,x_1MB_w)[,6]),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5)
lines(x_1MB_w_pos[[1]][,6], col = kenBlue)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(do.call(rbind,x_hyperrectangle_pos)[,6], na.rm = TRUE),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5)
lines(x_hyperrectangle_pos[[1]][,6], col = kenBlue)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(do.call(rbind,x_hyperrectangle_w_pos)[,6], na.rm = TRUE),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5)
lines(x_hyperrectangle_w_pos[[1]][,6], col = kenBlue)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(do.call(rbind,x_knapsack_pos)[,6], na.rm = TRUE),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5)
lines(x_knapsack_pos[[1]][,6], col = kenBlue)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(do.call(rbind,x_knapsack_w_pos)[,6], na.rm = TRUE),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5)
lines(x_knapsack_w_pos[[1]][,6], col = kenBlue)
dev.off()

##########################################                     Comparative ESS #####################

x_unif_FMB = LIPS::to_draws_df(x_FMB, matrix = TRUE)
x_unif_1MB = LIPS::to_draws_df(x_1MB, matrix = FALSE)
x_unif_1MB_w = LIPS::to_draws_df(x_1MB_w, matrix = FALSE)
x_unif_hyperrectangle = LIPS::to_draws_df(x_hyperrectangle, matrix = FALSE)
x_unif_hyperrectangle_w = LIPS::to_draws_df(x_hyperrectangle_w, matrix = FALSE)
x_unif_knapsack = LIPS::to_draws_df(x_knapsack, matrix = FALSE)
x_unif_knapsack_w = LIPS::to_draws_df(x_knapsack_w, matrix = FALSE)

pos_prop_1MB = 1-mean(sapply(x_1MB_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))/sapply(x_1MB_pos, nrow))
pos_prop_1MB_w = 1-mean(sapply(x_1MB_w_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))/sapply(x_1MB_pos, nrow))
pos_prop_hyperrectangle = 1-mean(sapply(x_hyperrectangle_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))/sapply(x_1MB_pos, nrow))
pos_prop_hyperrectangle_w = 1-mean(sapply(x_hyperrectangle_w_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))/sapply(x_1MB_pos, nrow))
pos_prop_knapsack = 1-mean(sapply(x_knapsack_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))/sapply(x_knapsack_pos, nrow))
pos_prop_knapsack_w = 1-mean(sapply(x_knapsack_w_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))/sapply(x_knapsack_pos, nrow))


ess_unif_FMB = posterior::summarise_draws(x_unif_FMB, "ess_basic")$ess_basic
ess_unif_1MB = posterior::summarise_draws(x_unif_1MB, "ess_basic")$ess_basic*pos_prop_1MB
ess_unif_1MB_w = posterior::summarise_draws(x_unif_1MB_w, "ess_basic")$ess_basic*pos_prop_1MB_w
ess_unif_hyperrectangle = posterior::summarise_draws(x_unif_hyperrectangle, "ess_basic")$ess_basic*pos_prop_hyperrectangle
ess_unif_hyperrectangle_w = posterior::summarise_draws(x_unif_hyperrectangle_w, "ess_basic")$ess_basic*pos_prop_hyperrectangle_w
ess_unif_knapsack = posterior::summarise_draws(x_unif_knapsack, "ess_basic")$ess_basic*pos_prop_knapsack
ess_unif_knapsack_w = posterior::summarise_draws(x_unif_knapsack_w, "ess_basic")$ess_basic*pos_prop_knapsack_w

par(mfrow = c(1,1))
barplot(log(rbind(ess_unif_FMB,
              ess_unif_1MB,
              ess_unif_1MB_w,
              ess_unif_hyperrectangle,
              ess_unif_hyperrectangle_w,
              ess_unif_knapsack,
              ess_unif_knapsack_w))/log(10),
        beside = TRUE,
        legend.text = c("FMB", "1MB", "1MB (w = 3)", "HR", "HR (w = 0.5)", "KS", "KS (w = 1)"),
        col = 1:7,
        names.arg = paste0("x", 1:ncol(A)),
        ylab = expression(log[10](ESS)),
        xlab = "Variable")
box()

boxplot(log(cbind(ess_unif_FMB,
                  ess_unif_1MB,
                  ess_unif_1MB_w,
                  ess_unif_hyperrectangle,
                  ess_unif_hyperrectangle_w,
                  ess_unif_knapsack,
                  ess_unif_knapsack_w))/log(10),
        names = c("FMB", "1MB", "1MB (w = 3)", "HR", "HR (w = 0.5)", "KS", "KS (w = 1)"),
        ylab = expression(log[10](ESS)),
        col = kenBlue,
        outcol = kenBlue,
        outpch = 19)

par(mfrow = c(1,1))
boxplot(log(cbind(ess_unif_1MB,
                  ess_unif_1MB_w,
                  ess_unif_hyperrectangle,
                  ess_unif_hyperrectangle_w,
                  ess_unif_knapsack,
                  ess_unif_knapsack_w)/ess_unif_FMB)/log(10),
        names = c("1MB", "1MB (w = 3)", "HR", "HR (w = 0.5)", "KS", "KS (w = 1)"),
        ylab = expression(log[10](ESS/ESS_FMB)),
        col = kenBlue,
        outcol = kenBlue,
        outpch = 19)

##########################################       Poisson model #####################################
##########################################                     Positive AND Negative ###############

pdf("highway_pois.pdf", width = 15, height = 13)
layout_matrix <- rbind(
  c(0,  0,  1,  2),
  c(0,  3,  7,  8),
  c(6, 4, 9, 10),
  c(0, 5, 11, 12)
)

layout(
  layout_matrix,
  widths  = c(1, 0.25, 1, 1),
  heights = c(0.30, 1, 1, 1)
)

par(mar = c(3.5, 3.5, 2, 1), mgp = c(1.6, 0.5, 0))

# Column title boxes
plot.new()
text(0.5, 0.5, "FMB", cex = 1.5, font = 2)

plot.new()
text(0.5, 0.5, "w=0", cex = 1.5, font = 2)

plot.new()
text(0.5, 0.5, "w \u2260 0", cex = 1.5, font = 2)

# Row title boxes
plot.new()
text(0.5, 0.5, "1MB-Extension", srt = 90, cex = 1.5, font = 2)

plot.new()
text(0.5, 0.5, "Hyperrectangle", srt = 90, cex = 1.5, font = 2)

plot.new()
text(0.5, 0.5, "Knapsack", srt = 90, cex = 1.5, font = 2)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(x_FMB_pois[,6]),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5,
     main = "FMB")
lines(x_FMB_pois[1:1e+05,6], col = kenBlue)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(do.call(rbind,x_1MB_pois)[,6]),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5)
lines(x_1MB_pois[[1]][,6], col = kenBlue)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(do.call(rbind,x_1MB_w_pois)[,6]),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5)
lines(x_1MB_w_pois[[1]][,6], col = kenBlue)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(do.call(rbind,x_hyperrectangle_pois)[,6]),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5)
lines(x_hyperrectangle_pois[[1]][,6], col = kenBlue)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(do.call(rbind,x_hyperrectangle_w_pois)[,6]),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5)
lines(x_hyperrectangle_w_pois[[1]][,6], col = kenBlue)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(do.call(rbind,x_knapsack_pois)[,6]),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5)
lines(x_knapsack_pois[[1]][,6], col = kenBlue)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(do.call(rbind,x_knapsack_w_pois)[,6]),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5)
lines(x_knapsack_w_pois[[1]][,6], col = kenBlue)
dev.off()

##########################################                     Only positive #######################

pdf("highway_pois_pos.pdf", width = 15, height = 13)
layout_matrix <- rbind(
  c(0,  0,  1,  2),
  c(0,  3,  7,  8),
  c(6, 4, 9, 10),
  c(0, 5, 11, 12)
)

layout(
  layout_matrix,
  widths  = c(1, 0.25, 1, 1),
  heights = c(0.30, 1, 1, 1)
)

par(mar = c(3.5, 3.5, 2, 1), mgp = c(1.6, 0.5, 0))

# Column title boxes

plot.new()
text(0.5, 0.5, "w=0", cex = 1.5, font = 2)

plot.new()
text(0.5, 0.5, "w \u2260 0", cex = 1.5, font = 2)

# Row title boxes
plot.new()
text(0.5, 0.5, "1MB-Extension", srt = 90, cex = 1.5, font = 2)

plot.new()
text(0.5, 0.5, "Hyperrectangle", srt = 90, cex = 1.5, font = 2)

plot.new()
text(0.5, 0.5, "Knapsack", srt = 90, cex = 1.5, font = 2)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(x_FMB_pois[,6]),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5,
     main = "FMB")
lines(x_FMB_pois[1:1e+05,6], col = kenBlue)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(do.call(rbind,x_1MB_pois_pos)[,6], na.rm = TRUE),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5)
lines(x_1MB_pois_pos[[1]][,6], col = kenBlue)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(do.call(rbind,x_1MB_w_pois_pos)[,6], na.rm = TRUE),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5)
lines(x_1MB_w_pois_pos[[1]][,6], col = kenBlue)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(do.call(rbind,x_hyperrectangle_pois_pos)[,6], na.rm = TRUE),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5)
lines(x_hyperrectangle_pois_pos[[1]][,6], col = kenBlue)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(do.call(rbind,x_hyperrectangle_w_pois_pos)[,6], na.rm = TRUE),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5)
lines(x_hyperrectangle_w_pois_pos[[1]][,6], col = kenBlue)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(do.call(rbind,x_knapsack_pois_pos)[,6], na.rm = TRUE),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5)
lines(x_knapsack_pois_pos[[1]][,6], col = kenBlue)

plot(NA,
     xlim = c(0,1e+05),
     ylim = range(do.call(rbind,x_knapsack_w_pois_pos)[,6], na.rm = TRUE),
     xlab = "Iteration",
     ylab = expression("x"[6]),
     cex.lab = 1.5)
lines(x_knapsack_w_pois_pos[[1]][,6], col = kenBlue)
dev.off()

##########################################                     Comparative ESS #####################

x_pois_FMB = LIPS::to_draws_df(x_FMB_pois, matrix = TRUE)
x_pois_1MB = LIPS::to_draws_df(x_1MB_pois, matrix = FALSE)
x_pois_1MB_w = LIPS::to_draws_df(x_1MB_w_pois, matrix = FALSE)
x_pois_hyperrectangle = LIPS::to_draws_df(x_hyperrectangle_pois, matrix = FALSE)
x_pois_hyperrectangle_w = LIPS::to_draws_df(x_hyperrectangle_w_pois, matrix = FALSE)
x_pois_knapsack = LIPS::to_draws_df(x_knapsack_pois, matrix = FALSE)
x_pois_knapsack_w = LIPS::to_draws_df(x_knapsack_w_pois, matrix = FALSE)

pos_prop_1MB_pois = 1-mean(sapply(x_1MB_pois_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))/sapply(x_1MB_pois_pos, nrow))
pos_prop_1MB_w_pois = 1-mean(sapply(x_1MB_w_pois_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))/sapply(x_1MB_w_pois_pos, nrow))
pos_prop_hyperrectangle_pois = 1-mean(sapply(x_hyperrectangle_pois_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))/sapply(x_hyperrectangle_pois_pos, nrow))
pos_prop_hyperrectangle_w_pois = 1-mean(sapply(x_hyperrectangle_w_pois_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))/sapply(x_hyperrectangle_w_pois_pos, nrow))
pos_prop_knapsack_pois = 1-mean(sapply(x_knapsack_pois_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))/sapply(x_knapsack_pois_pos, nrow))
pos_prop_knapsack_w_pois = 1-mean(sapply(x_knapsack_w_pois_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))/sapply(x_knapsack_w_pois_pos, nrow))

ess_pois_FMB = posterior::summarise_draws(x_pois_FMB, "ess_basic")$ess_basic
ess_pois_1MB = posterior::summarise_draws(x_pois_1MB, "ess_basic")$ess_basic*pos_prop_1MB_pois
ess_pois_1MB_w = posterior::summarise_draws(x_pois_1MB_w, "ess_basic")$ess_basic*pos_prop_1MB_w_pois
ess_pois_hyperrectangle = posterior::summarise_draws(x_pois_hyperrectangle, "ess_basic")$ess_basic*pos_prop_hyperrectangle_pois
ess_pois_hyperrectangle_w = posterior::summarise_draws(x_pois_hyperrectangle_w, "ess_basic")$ess_basic*pos_prop_hyperrectangle_w_pois
ess_pois_knapsack = posterior::summarise_draws(x_pois_knapsack, "ess_basic")$ess_basic*pos_prop_knapsack_pois
ess_pois_knapsack_w = posterior::summarise_draws(x_pois_knapsack_w, "ess_basic")$ess_basic*pos_prop_knapsack_w_pois

par(mfrow = c(1,1))
barplot(log(rbind(ess_pois_FMB,
                  ess_pois_1MB,
                  ess_pois_1MB_w,
                  ess_pois_hyperrectangle,
                  ess_pois_hyperrectangle_w,
                  ess_pois_knapsack,
                  ess_pois_knapsack_w))/log(10),
        beside = TRUE,
        legend.text = c("FMB", "1MB", "1MB (w = 3)", "HR", "HR (w = 0.5)", "KS", "KS (w = 0.5)"),
        col = 1:7,
        names.arg = paste0("x", 1:ncol(A)),
        xlab = "Variable",
        ylab = expression(log[10](ESS)))
box()

par(mfrow = c(1,1))
boxplot(log(cbind(ess_pois_FMB,
                  ess_pois_1MB,
                  ess_pois_1MB_w,
                  ess_pois_hyperrectangle,
                  ess_pois_hyperrectangle_w,
                  ess_pois_knapsack,
                  ess_pois_knapsack_w))/log(10),
        names = c("FMB", "1MB", "1MB (w = 3)", "HR", "HR (w = 0.5)", "KS", "KS (w = 0.5)"),
        ylab = expression(log[10](ESS)),
        col = kenBlue,
        outcol = kenBlue,
        pch = 19)

par(mfrow = c(1,1))
boxplot(log(cbind(ess_unif_1MB,
                  ess_unif_1MB_w,
                  ess_unif_hyperrectangle,
                  ess_unif_hyperrectangle_w,
                  ess_unif_knapsack,
                  ess_unif_knapsack_w)/ess_pois_FMB)/log(10),
        names = c("1MB", "1MB (w = 3)", "HR", "HR (w = 0.5)", "KS", "KS (w = 0.5)"),
        ylab = expression(log[10](ESS/ESS_FMB)),
        col = kenBlue,
        outcol = kenBlue,
        outpch = 19)
