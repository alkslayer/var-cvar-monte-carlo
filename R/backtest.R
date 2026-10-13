# Backtest VaR: kroczące prognozy bez look-ahead bias i testy statystyczne.

# xlogy: x * log(y) z konwencją 0 * log(0) = 0 (potrzebne w testach LR).
xlogy <- function(x, y) ifelse(x == 0, 0, x * log(y))

# test_kupca: test POF (proportion of failures) Kupca.
# H0: odsetek przekroczeń = p = 1 - alfa. LR ~ chi2(1).
test_kupca <- function(I, alfa) {
  p <- 1 - alfa; Tn <- length(I); x <- sum(I); ph <- x / Tn
  lr <- -2 * (xlogy(Tn - x, 1 - p) + xlogy(x, p) - xlogy(Tn - x, 1 - ph) - xlogy(x, ph))
  c(LR_uc = lr, p_uc = pchisq(lr, 1, lower.tail = FALSE))
}

# test_christoffersena: test niezależności przekroczeń (łańcuch Markowa 1. rzędu)
# oraz test łączny (conditional coverage) = Kupiec + niezależność, LR_cc ~ chi2(2).
test_christoffersena <- function(I, alfa) {
  I <- as.integer(I); a <- I[-length(I)]; b <- I[-1]
  n00 <- sum(a == 0 & b == 0); n01 <- sum(a == 0 & b == 1)
  n10 <- sum(a == 1 & b == 0); n11 <- sum(a == 1 & b == 1)
  p01 <- n01 / (n00 + n01); p11 <- if (n10 + n11 > 0) n11 / (n10 + n11) else 0
  p   <- (n01 + n11) / (n00 + n01 + n10 + n11)
  lr_ind <- -2 * (xlogy(n00 + n10, 1 - p) + xlogy(n01 + n11, p) -
                    xlogy(n00, 1 - p01) - xlogy(n01, p01) - xlogy(n10, 1 - p11) - xlogy(n11, p11))
  lr_uc <- test_kupca(I, alfa)[["LR_uc"]]
  lr_cc <- lr_uc + lr_ind
  c(n11 = n11, LR_ind = lr_ind, p_ind = pchisq(lr_ind, 1, lower.tail = FALSE),
    LR_cc = lr_cc, p_cc = pchisq(lr_cc, 2, lower.tail = FALSE))
}

# prognozy_kroczace: dla każdego dnia t > okno liczy kwantyl (1 - alfa) zwrotu log portfela
# funkcją f, która dostaje WYŁĄCZNIE okno R[(t-okno):(t-1), ] (dane sprzed dnia prognozy).
# f(okno_R, stan) zwraca wektor kwantyli (po jednym na alfę); stan pozwala co REFIT_CO dni
# przeliczać kosztowne parametry (np. nu) i zapamiętać je do kolejnych dni.
prognozy_kroczace <- function(R, f, okno = OKNO_BT, alfy = ALFY, co_ile = REFIT_CO) {
  R <- as.matrix(R); Tn <- nrow(R); dni <- (okno + 1):Tn
  wynik <- matrix(NA_real_, length(dni), length(alfy), dimnames = list(NULL, paste0("q", alfy)))
  stan <- NULL
  for (j in seq_along(dni)) {
    t <- dni[j]
    okno_R <- R[(t - okno):(t - 1), , drop = FALSE]
    odswiez <- (j - 1) %% co_ile == 0
    out <- f(okno_R, stan, odswiez)
    wynik[j, ] <- out$q; stan <- out$stan
  }
  wynik
}

# Funkcje prognozy dla poszczególnych metod (argumenty: okno, stan, odswiez).
f_historyczna <- function(okno_R, stan, odswiez) {
  rp <- as.numeric(okno_R %*% WAGI)
  list(q = quantile(rp, 1 - ALFY, names = FALSE, type = 7), stan = NULL)
}
f_normalna <- function(okno_R, stan, odswiez) {
  rp <- as.numeric(okno_R %*% WAGI)
  list(q = mean(rp) + qnorm(1 - ALFY) * sd(rp), stan = NULL)
}
f_mc_t <- function(okno_R, stan, odswiez) {
  nu <- if (odswiez || is.null(stan)) dopasuj_nu(okno_R) else stan
  sim <- sim_mc_t(N_SIM_BT, colMeans(okno_R), cov(okno_R), nu)
  list(q = quantile(sim, 1 - ALFY, names = FALSE, type = 7), stan = nu)
}

# podsumuj_backtest: liczba przekroczeń, oczekiwana liczba i testy dla jednej metody i alfy.
podsumuj_backtest <- function(r, q, alfa) {
  I <- as.integer(r < q)
  c(dni = length(I), przekroczenia = sum(I), oczekiwane = length(I) * (1 - alfa),
    odsetek = mean(I), test_kupca(I, alfa), test_christoffersena(I, alfa))
}
