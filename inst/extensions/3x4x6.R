########################################### Setup ###################################################

kenBlue = "#3A8EE4"

A = LIPS::contingencyTable$A
M = LIPS::contingencyTable$MarkovBasis
y = LIPS::contingencyTable$y
lambda = LIPS::contingencyTable$lambda
B = LIPS::contingencyTable$BasicMoves
CPLB = LIPS::contingencyTable$LatticeBasis

# CPLB
A1 = A[,1]
A1_idx = 1
for(i in 1:ncol(A)){
  if(qr(A1)$rank < qr(cbind(A1, A[,i]))$rank){
    A1 = cbind(A1, A[,i])
    A1_idx = c(A1_idx, i)
  }
}
A2_idx = setdiff(1:ncol(A), A1_idx)

a = c()
j = 1
for(i in A2_idx){
  a[j] = lpSolve::lp(direction = "max",
                     objective.in = diag(72)[i,],
                     const.mat = A,
                     const.rhs = y,
                     const.dir = "=",
                     all.int = TRUE)$solution[i]
  j = j+1
}

objective = rep(0, ncol(A))
objective[A2_idx] = 1

b = sum(lpSolve::lp(direction = "max",
                    objective.in = objective,
                    const.mat = A,
                    const.rhs = y,
                    const.dir = "==",
                    all.int = TRUE)$solution[A2_idx])

#Check if CPLB and 1MB are in FMB (Answer: Yes!) => Fair comparison

isInFMB = apply(B, 2, function(col) any(apply((col == M) | (-col == M), 2, all)))
isInFMB = apply(CPLB, 2, function(col) any(apply((col == M) | (-col == M), 2, all)))

########################################### Samplers ###############################################
###########################################          Uniform Model ################################

# FMB
set.seed(2026)
x_FMB = LIPS::xSampler(A = A,
                       y = y,
                       B = M,
                       Model = "Unif",
                       proposal = "Random",
                       nSample = 1e+06,
                       nBurnin = 1,
                       thinning = 10,
                       nChains = 1)

set.seed(2026)
x_FMB_systematic = LIPS::xSampler(A = A,
                                   y = y,
                                   B = M,
                                   Model = "Unif",
                                   proposal = "NonRandom",
                                   nSample = 1e+06,
                                   nBurnin = 1e+05,
                                   thinning = 10,
                                  nChains = 1)

# 1-extension
set.seed(2026)
x_1MB = LIPS::extensionSampler(A = A,
                               y = y,
                               B = B,
                               extension = "p",
                               p = rep(1, ncol(A)),
                               loggamma = -13,
                               nSample = 1e+06,
                               nBurnin = 1,
                               thinning = 10,
                               nChains = 1)

set.seed(2026)
x_1MB_w = LIPS::extensionSampler(A = A,
                                 y = y,
                                 B = B,
                                 extension = "p",
                                 p = rep(1, ncol(A)),
                                 loggamma = -1.5,
                                 w = 3,
                                 nSample = 1e+06,
                                 nBurnin = 1e+05,
                                 thinning = 10,
                                 nChains = 1)

x_1MB_pos = LIPS::only_positive(x_1MB)
sapply(x_1MB_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))
x_1MB_w_pos = LIPS::only_positive(x_1MB_w)
sapply(x_1MB_w_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))


# Hyperrectangle
set.seed(2026)
x_hyperrectangle = LIPS::extensionSampler(A = A,
                                          y = y,
                                          extension = "hyperrectangle",
                                          B = CPLB,
                                          A2Idx = A2_idx,
                                          a = a,
                                          loggamma = -13,
                                          nSample = 1e+06,
                                          nBurnin = 1e+05,
                                          thinning = 10,
                                          nChains = 1)

set.seed(2026)
x_hyperrectangle_w = LIPS::extensionSampler(A = A,
                                            y = y,
                                            extension = "hyperrectangle",
                                            B = CPLB,
                                            A2Idx = A2_idx,
                                            a = a,
                                            w = 2,
                                            loggamma = -3,
                                            nSample = 1e+06,
                                            nBurnin = 1e+05,
                                            thinning = 10,
                                            nChains = 1)

x_hyperrectangle_pos = LIPS::only_positive(x_hyperrectangle)
sapply(x_hyperrectangle_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))
x_hyperrectangle_w_pos = LIPS::only_positive(x_hyperrectangle_w)
sapply(x_hyperrectangle_w_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))

# Knapsack
set.seed(2026)
x_knapsack = LIPS::extensionSampler(A = A,
                                    y = y,
                                    extension = "knapsack",
                                    B = CPLB,
                                    A2Idx = A2_idx,
                                    a = rep(1, ncol(A)-nrow(A)),
                                    b = b,
                                    loggamma = -13,
                                    nSample = 1e+06,
                                    nBurnin = 1,
                                    thinning = 10,
                                    nChains = 1)

set.seed(2026)
x_knapsack_w = LIPS::extensionSampler(A = A,
                                      y = y,
                                      extension = "knapsack",
                                      B = CPLB,
                                      A2Idx = A2_idx,
                                      a = rep(1, ncol(A)-nrow(A)),
                                      b = b,
                                      loggamma = -3.5,
                                      w = 2,
                                      nSample = 1e+06,
                                      nBurnin = 1,
                                      thinning = 10,
                                      nChains = 1)

x_knapsack_pos = LIPS::only_positive(x_knapsack)
sapply(x_knapsack_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))
x_knapsack_w_pos = LIPS::only_positive(x_knapsack_w)
sapply(x_knapsack_w_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))

###########################################          Poisson Model #############################

# FMB
set.seed(2026)
x_pois_FMB = LIPS::xSampler(A = A,
                            y = y,
                            B = M,
                            Model = "Pois",
                            lambda = lambda,
                            proposal = "Random",
                            nSample = 1e+06,
                            nBurnin = 1,
                            thinning = 10,
                            nChains = 1)

set.seed(2026)
x_pois_FMB_systematic = LIPS::xSampler(A = A,
                                        y = y,
                                        B = M,
                                        Model = "Unif",
                                        proposal = "NonRandom",
                                        nSample = 1e+06,
                                        nBurnin = 1e+05,
                                        thinning = 10,
                                        nChains = 1)

# 1MB
set.seed(2026)
x_pois_1MB = LIPS::extensionSampler(A = A,
                                    y = y,
                                    B = B,
                                    Model = "poisson",
                                    lambda = lambda,
                                    extension = "p",
                                    p = 1,
                                    loggamma = -17,
                                    nSample = 1e+06,
                                    nBurnin = 1,
                                    thinning = 10,
                                    nChains = 1)

set.seed(2026)
x_pois_1MB_w = LIPS::extensionSampler(A = A,
                                      y = y,
                                      B = B,
                                      Model = "poisson",
                                      lambda = lambda,
                                      extension = "p",
                                      p = 1,
                                      loggamma = -4.5,
                                      w = 3,
                                      nSample = 1e+06,
                                      nBurnin = 1,
                                      thinning = 10,
                                      nChains = 1)

x_pois_1MB_pos = LIPS::only_positive(x_pois_1MB)
sapply(x_pois_1MB_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))
x_pois_1MB_w_pos = LIPS::only_positive(x_pois_1MB_w)
sapply(x_pois_1MB_w_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))

set.seed(2026)
x_pois_hyperrectangle = LIPS::extensionSampler(A = A,
                                               y = y,
                                               B = CPLB,
                                               Model = "poisson",
                                               lambda = lambda,
                                               extension = "hyperrectangle",
                                               a = a,
                                               A2Idx = A2_idx,
                                               loggamma = -16,
                                               nSample = 1e+06,
                                               nBurnin = 1,
                                               thinning = 10,
                                               nChains = 1)

set.seed(2026)
x_pois_hyperrectangle_w = LIPS::extensionSampler(A = A,
                                                 y = y,
                                                 B = CPLB,
                                                 Model = "poisson",
                                                 lambda = lambda,
                                                 extension = "hyperrectangle",
                                                 a = a,
                                                 A2Idx = A2_idx,
                                                 loggamma = -10,
                                                 w = 1,
                                                 nSample = 1e+06,
                                                 nBurnin = 1,
                                                 thinning = 10,
                                                 nChains = 1)

x_pois_hyperrectangle_pos = LIPS::only_positive(x_pois_hyperrectangle)
sapply(x_pois_hyperrectangle_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))
x_pois_hyperrectangle_w_pos = LIPS::only_positive(x_pois_hyperrectangle_w)
sapply(x_pois_hyperrectangle_w_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))

set.seed(2026)
x_pois_knapsack = LIPS::extensionSampler(A = A,
                                         y = y,
                                         B = CPLB,
                                         Model = "poisson",
                                         lambda = lambda,
                                         extension = "knapsack",
                                         a = rep(1, ncol(A) - nrow(A)),
                                         b = b,
                                         A2Idx = A2_idx,
                                         loggamma = -15,
                                         nSample = 1e+06,
                                         nBurnin = 1,
                                         thinning = 10,
                                         nChains = 1)

set.seed(2026)
x_pois_knapsack_w = LIPS::extensionSampler(A = A,
                                           y = y,
                                           B = CPLB,
                                           Model = "poisson",
                                           lambda = lambda,
                                           extension = "knapsack",
                                           a = rep(1, ncol(A) - nrow(A)),
                                           b = b,
                                           A2Idx = A2_idx,
                                           loggamma = -10,
                                           w = 1,
                                           nSample = 1e+06,
                                           nBurnin = 1,
                                           thinning = 10,
                                           nChains = 1)

x_pois_knapsack_pos = LIPS::only_positive(x_pois_knapsack)
sapply(x_pois_knapsack_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))
x_pois_knapsack_w_pos = LIPS::only_positive(x_pois_knapsack_w)
sapply(x_pois_knapsack_w_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))


########################################### Plots    ###############################################
###########################################          Uniform Model #################################
###########################################                        Systematic vs Random ###############

# par(mfrow = c(1,2))
# plot(NA,
#      xlim = c(1,1e+05),
#      ylim = c(min(x_FMB[,59]), max(x_FMB[,59])),
#      xlab = "Iteration",
#      ylab = expression("x"[59]))
# lines(x_FMB[1:1e+05,59], type = "l", col = 1)
# lines(x_FMB[1e+05+(1:1e+05),59], type = "l", col = 2)
# lines(x_FMB[2e+05+(1:1e+05),59], type = "l", col = 3)
# lines(x_FMB[3e+05+(1:1e+05),59], type = "l", col = 4)
#
# plot(NA,
#      xlim = c(1,1e+05),
#      ylim = c(min(x_FMB_systematic[,59]), max(x_FMB_systematic[,59])),
#      xlab = "Iteration",
#      ylab = expression("x"[59]))
# lines(x_FMB_systematic[1:1e+05,59], type = "l", col = 1)
# lines(x_FMB_systematic[1e+05+(1:1e+05),59], type = "l", col = 2)
# lines(x_FMB_systematic[2e+05+(1:1e+05),59], type = "l", col = 3)
# lines(x_FMB_systematic[3e+05+(1:1e+05),59], type = "l", col = 4)

###########################################                        Positive AND negative ##############

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
     xlim = c(0, 1e+05),
     ylim = c(min(x_FMB[,59]),max(x_FMB[,59])),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5,
     main = "FMB")
lines(x_FMB[1:1e+05,59], col = kenBlue)

plot(NA,
     xlim = c(0, 1e+05),
     ylim = range(do.call(rbind, x_1MB)[,59]),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5)
lines(x_1MB[[1]][,59], col = kenBlue)

plot(NA,
     xlim = c(0, 1e+05),
     ylim = range(do.call(rbind, x_1MB_w)[,59]),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5)
lines(x_1MB_w[[1]][,59], col = kenBlue)

plot(NA,
     xlim = c(0, 1e+05),
     ylim = range(do.call(rbind, x_hyperrectangle)[,59]),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5)
lines(x_hyperrectangle[[1]][,59], col = kenBlue)

plot(NA,
     xlim = c(0, 1e+05),
     ylim = range(do.call(rbind, x_hyperrectangle_w)[,59]),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5)
lines(x_hyperrectangle_w[[1]][,59], col = kenBlue)

plot(NA,
     xlim = c(0, 1e+05),
     ylim = range(do.call(rbind, x_knapsack)[,59]),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5)
lines(x_knapsack[[1]][,59], col = kenBlue)

plot(NA,
     xlim = c(0, 1e+05),
     ylim = range(do.call(rbind, x_knapsack_w)[,59]),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5)
lines(x_knapsack_w[[1]][,59], col = kenBlue)

###########################################                        Only positive ######################

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
     xlim = c(0, 1e+05),
     ylim = c(min(x_FMB[,59]),max(x_FMB[,59])),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5,
     main = "FMB")
lines(x_FMB[1:1e+05,59], col = kenBlue)

plot(NA,
     xlim = c(0, 1e+05),
     ylim = range(do.call(rbind, x_1MB)[,59]),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5)
lines(x_1MB_pos[[1]][,59], col = kenBlue)

plot(NA,
     xlim = c(0, 1e+05),
     ylim = range(do.call(rbind, x_1MB_w)[,59]),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5)
lines(x_1MB_w_pos[[1]][,59], col = kenBlue)

plot(NA,
     xlim = c(0, 1e+05),
     ylim = range(do.call(rbind, x_hyperrectangle_pos)[,59], na.rm = TRUE),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5)
lines(x_hyperrectangle_pos[[1]][,59], col = kenBlue)

plot(NA,
     xlim = c(0, 1e+05),
     ylim = range(do.call(rbind, x_hyperrectangle_w_pos)[,59], na.rm = TRUE),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5)
lines(x_hyperrectangle_w_pos[[1]][,59], col = kenBlue)

plot(NA,
     xlim = c(0, 1e+05),
     ylim = range(do.call(rbind, x_knapsack_pos)[,59], na.rm = TRUE),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5)
lines(x_knapsack_pos[[1]][,59], col = kenBlue)

plot(NA,
     xlim = c(0, 1e+05),
     ylim = range(do.call(rbind, x_knapsack_w_pos)[,59], na.rm = TRUE),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5)
lines(x_knapsack_w_pos[[1]][,59], col = kenBlue)

###########################################                        Comparative ESS #################

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
pos_prop_knapsack = 1-mean(sapply(x_knapsack_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))/sapply(x_1MB_pos, nrow))
pos_prop_knapsack_w = 1-mean(sapply(x_knapsack_w_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))/sapply(x_1MB_pos, nrow))


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
        legend.text = c("FMB", "1MB", "1MB (w = 3)", "HR", "HR (w = 2)", "KS", "KS (w = 2)"),
        col = 1:7,
        names.arg = paste0("x", 1:ncol(A)),
        xlab = "Variable",
        ylab = expression(log[10](ESS)),
        cex.names = 0.5)

boxplot(log(cbind(ess_unif_FMB,
                  ess_unif_1MB,
                  ess_unif_1MB_w,
                  ess_unif_hyperrectangle,
                  ess_unif_hyperrectangle_w,
                  ess_unif_knapsack,
                  ess_unif_knapsack_w))/log(10),
        names = c("FMB", "1MB", "1MB (w = 3)", "HR", "HR (w = 2)", "KS", "KS (w = 2)"),
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
        names = c("1MB", "1MB (w = 3)", "HR", "HR (w = 1)", "KS", "KS (w = 1)"),
        ylab = expression(log[10](ESS/ESS_FMB)),
        col = kenBlue,
        outcol = kenBlue,
        outpch = 19)

###########################################          Poisson Model #################################
###########################################                        Systematic vs Random ###############

# par(mfrow = c(1,2))
# plot(NA,
#      xlim = c(1,1e+05),
#      ylim = c(min(x_pois_FMB[,59]), max(x_pois_FMB[,59])),
#      xlab = "Iteration",
#      ylab = expression("x"[59]))
# lines(x_pois_FMB[1:1e+05,59], type = "l", col = 1)
# lines(x_pois_FMB[1e+05+(1:1e+05),59], type = "l", col = 2)
# lines(x_pois_FMB[2e+05+(1:1e+05),59], type = "l", col = 3)
# lines(x_pois_FMB[3e+05+(1:1e+05),59], type = "l", col = 4)
#
# plot(NA,
#      xlim = c(1,1e+05),
#      ylim = c(min(x_pois_FMB_systematic[,59]), max(x_pois_FMB_systematic[,59])),
#      xlab = "Iteration",
#      ylab = expression("x"[59]))
# lines(x_pois_FMB_systematic[1:1e+05,59], type = "l", col = 1)
# lines(x_pois_FMB_systematic[1e+05+(1:1e+05),59], type = "l", col = 2)
# lines(x_pois_FMB_systematic[2e+05+(1:1e+05),59], type = "l", col = 3)
# lines(x_pois_FMB_systematic[3e+05+(1:1e+05),59], type = "l", col = 4)

###########################################                        Positive AND negative ##############

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
     xlim = c(0, 1e+05),
     ylim = c(min(x_pois_FMB[,59]),max(x_pois_FMB[,59])),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5,
     main = "FMB")
lines(x_pois_FMB[1:1e+05,59], col = kenBlue)

plot(NA,
     xlim = c(0, 1e+05),
     ylim = range(x_pois_1MB[[1]][,59]),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5)
lines(x_pois_1MB[[1]][,59], col = kenBlue)

plot(NA,
     xlim = c(0, 1e+05),
     ylim = range(x_pois_1MB_w[[1]][,59]),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5)
lines(x_pois_1MB_w[[1]][,59], col = kenBlue)

plot(NA,
     xlim = c(0, 1e+05),
     ylim = range(x_pois_hyperrectangle[[1]][,59]),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5)
lines(x_pois_hyperrectangle[[1]][,59], col = kenBlue)

plot(NA,
     xlim = c(0, 1e+05),
     ylim = range(x_pois_hyperrectangle_w[[1]][,59]),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5)
lines(x_pois_hyperrectangle_w[[1]][,59], col = kenBlue)

plot(NA,
     xlim = c(0, 1e+05),
     ylim = range(x_pois_knapsack[[1]][,59]),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5)
lines(x_pois_knapsack[[1]][,59], col = kenBlue)

plot(NA,
     xlim = c(0, 1e+05),
     ylim = range(x_pois_knapsack_w[[1]][,59]),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5)
lines(x_pois_knapsack_w[[1]][,59], col = kenBlue)

###########################################                        Only positive ######################
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
     xlim = c(0, 1e+05),
     ylim = range(x_pois_FMB[,59]),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5,
     main = "FMB")
lines(x_pois_FMB[1:1e+05,59], col = kenBlue)

plot(NA,
     xlim = c(0, 1e+05),
     ylim = range(x_pois_1MB_pos[[1]][,59], na.rm = TRUE),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5)
lines(x_pois_1MB_pos[[1]][,59], col = kenBlue)

plot(NA,
     xlim = c(0, 1e+05),
     ylim = range(x_pois_1MB_w_pos[[1]][,59], na.rm = TRUE),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5)
lines(x_pois_1MB_w_pos[[1]][,59], col = kenBlue)

plot(NA,
     xlim = c(0, 1e+05),
     ylim = range(x_pois_hyperrectangle_pos[[1]][,59], na.rm = TRUE),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5)
lines(x_pois_hyperrectangle_pos[[1]][,59], col = kenBlue)

plot(NA,
     xlim = c(0, 1e+05),
     ylim = range(x_pois_hyperrectangle_w_pos[[1]][,59], na.rm = TRUE),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5)
lines(x_pois_hyperrectangle_w_pos[[1]][,59], col = kenBlue)

plot(NA,
     xlim = c(0, 1e+05),
     ylim = range(x_pois_knapsack_pos[[1]][,59], na.rm = TRUE),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5)
lines(x_pois_knapsack_pos[[1]][,59], col = kenBlue)

plot(NA,
     xlim = c(0, 1e+05),
     ylim = range(x_pois_knapsack_w_pos[[1]][,59], na.rm = TRUE),
     xlab = "Iteration",
     ylab = expression(x[59]),
     cex.lab = 1.5)
lines(x_pois_knapsack_w_pos[[1]][,59], col = kenBlue)

###########################################                        Comparative ESS #################

x_FMB_pois = LIPS::to_draws_df(x_pois_FMB, matrix = TRUE)
x_1MB_pois = LIPS::to_draws_df(x_pois_1MB, matrix = FALSE)
x_1MB_w_pois = LIPS::to_draws_df(x_pois_1MB_w, matrix = FALSE)
x_hyperrectangle_pois = LIPS::to_draws_df(x_pois_hyperrectangle, matrix = FALSE)
x_hyperrectangle_w_pois = LIPS::to_draws_df(x_pois_hyperrectangle_w, matrix = FALSE)
x_knapsack_pois = LIPS::to_draws_df(x_pois_knapsack, matrix = FALSE)
x_knapsack_w_pois = LIPS::to_draws_df(x_pois_knapsack_w, matrix = FALSE)

pos_prop_1MB_pois = 1-mean(sapply(x_pois_1MB_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))/sapply(x_pois_1MB_pos, nrow))
pos_prop_1MB_w_pois = 1-mean(sapply(x_pois_1MB_w_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))/sapply(x_pois_1MB_w_pos, nrow))
pos_prop_hyperrectangle_pois = 1-mean(sapply(x_pois_hyperrectangle_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))/sapply(x_pois_hyperrectangle_pos, nrow))
pos_prop_hyperrectangle_w_pois = 1-mean(sapply(x_pois_hyperrectangle_w_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))/sapply(x_pois_hyperrectangle_w_pos, nrow))
pos_prop_knapsack_pois = 1-mean(sapply(x_pois_knapsack_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))/sapply(x_pois_knapsack_pos, nrow))
pos_prop_knapsack_w_pois = 1-mean(sapply(x_pois_knapsack_w_pos, function(chain) sum(apply(chain, 1, function(row) any(is.na(row)))))/sapply(x_pois_knapsack_w_pos, nrow))

ess_pois_FMB = posterior::summarise_draws(x_FMB_pois, "ess_basic")$ess_basic
ess_pois_1MB = posterior::summarise_draws(x_1MB_pois, "ess_basic")$ess_basic*pos_prop_1MB_pois
ess_pois_1MB_w = posterior::summarise_draws(x_1MB_w_pois, "ess_basic")$ess_basic*pos_prop_1MB_w_pois
ess_pois_hyperrectangle = posterior::summarise_draws(x_hyperrectangle_pois, "ess_basic")$ess_basic*pos_prop_hyperrectangle_pois
ess_pois_hyperrectangle_w = posterior::summarise_draws(x_hyperrectangle_w_pois, "ess_basic")$ess_basic*pos_prop_hyperrectangle_w_pois
ess_pois_knapsack = posterior::summarise_draws(x_knapsack_pois, "ess_basic")$ess_basic*pos_prop_knapsack_pois
ess_pois_knapsack_w = posterior::summarise_draws(x_knapsack_w_pois, "ess_basic")$ess_basic*pos_prop_knapsack_w_pois

par(mfrow = c(1,1))
barplot(log(rbind(ess_pois_FMB,
                  ess_pois_1MB,
                  ess_pois_1MB_w,
                  ess_pois_hyperrectangle,
                  ess_pois_hyperrectangle_w,
                  ess_pois_knapsack,
                  ess_pois_knapsack_w))/log(10),
        beside = TRUE,
        legend.text = c("FMB", "1MB", "1MB (w = 3)", "HR", "HR (w = 1)", "KS", "KS (w = 1)"),
        col = 1:7,
        names.arg = paste0("x", 1:ncol(A)),
        xlab = "Variable",
        ylab = expression(x[59]))
box()

par(mfrow = c(1,1))
boxplot(log(cbind(ess_pois_FMB,
                  ess_pois_1MB,
                  ess_pois_1MB_w,
                  ess_pois_hyperrectangle,
                  ess_pois_hyperrectangle_w,
                  ess_pois_knapsack,
                  ess_pois_knapsack_w))/log(10),
        names = c("FMB", "1MB", "1MB (w = 3)", "HR", "HR (w = 1)", "KS", "KS (w = 1)"),
        ylab = expression(log[10](ESS)),
        col = kenBlue,
        outcol = kenBlue,
        outpch = 19)

par(mfrow = c(1,1))
boxplot(log(cbind(ess_pois_1MB,
                  ess_pois_1MB_w,
                  ess_pois_hyperrectangle,
                  ess_pois_hyperrectangle_w,
                  ess_pois_knapsack,
                  ess_pois_knapsack_w)/ess_pois_FMB)/log(10),
        names = c("1MB", "1MB (w = 3)", "HR", "HR (w = 1)", "KS", "KS (w = 1)"),
        ylab = expression(log[10](ESS/ESS_FMB)),
        col = kenBlue,
        outcol = kenBlue,
        outpch = 19)
