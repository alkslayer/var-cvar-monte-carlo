# Symulacje Monte Carlo zwrotów portfela (wektoryzowane, bez pętli po symulacjach).
# Każda funkcja zwraca wektor n symulowanych h-dniowych zwrotów log portfela.
# Zwroty aktywów są symulowane dzień po dniu i sumowane w horyzoncie h
# (dla t-Studenta i bootstrapu suma h dni NIE ma już rozkładu wyjściowego).

# sumuj_horyzont: wektor n*h zwrotów dziennych -> n zwrotów h-dniowych.
sumuj_horyzont <- function(x, n, h) {
  if (h == 1) return(x)
  rowSums(matrix(x, nrow = n, ncol = h))
}

# sim_mc_normal: wielowymiarowy rozkład normalny N(mu, Sigma) przez rozkład Choleskiego.
# chol() zwraca macierz górnotrójkątną U (Sigma = U'U), więc X = Z %*% U ma kowariancję Sigma.
sim_mc_normal <- function(n, mu, Sigma, w = WAGI, h = 1) {
  U <- chol(Sigma)
  Z <- matrix(rnorm(n * h * length(mu)), ncol = length(mu))
  X <- sweep(Z %*% U, 2, mu, "+")
  sumuj_horyzont(as.numeric(X %*% w), n, h)
}

# sim_mc_t: wielowymiarowy t-Student o nu stopniach swobody i kowariancji Sigma.
# X = mu + sqrt((nu-2)/W) * Z U, W ~ chi2(nu): wspólny mianownik dla wszystkich aktywów
# danego dnia daje grube ogony i zależność w ogonach (krachy dotykają wszystkich naraz).
sim_mc_t <- function(n, mu, Sigma, nu, w = WAGI, h = 1) {
  U <- chol(Sigma)
  m <- n * h
  Z <- matrix(rnorm(m * length(mu)), ncol = length(mu))
  skala <- sqrt((nu - 2) / rchisq(m, df = nu))
  X <- sweep((Z %*% U) * skala, 2, mu, "+")
  sumuj_horyzont(as.numeric(X %*% w), n, h)
}

# sim_bootstrap: losowanie ze zwracaniem całych dni historycznych (wierszy macierzy zwrotów),
# co zachowuje korelacje i grube ogony z danych, ale niszczy zależność w czasie.
sim_bootstrap <- function(n, R, w = WAGI, h = 1) {
  rp <- as.numeric(as.matrix(R) %*% w)
  sumuj_horyzont(rp[sample.int(length(rp), n * h, replace = TRUE)], n, h)
}

# loglik_mvt: log-wiarygodność wielowymiarowego t (średnia mu, kowariancja Sigma, nu stopni swobody).
# Macierz skali S = Sigma * (nu-2)/nu, żeby kowariancja rozkładu była równa Sigma.
loglik_mvt <- function(nu, R, mu, Sigma) {
  R <- as.matrix(R); k <- ncol(R); S <- Sigma * (nu - 2) / nu
  d2 <- mahalanobis(R, mu, S)
  sum(lgamma((nu + k) / 2) - lgamma(nu / 2) - k / 2 * log(nu * pi) -
        0.5 * as.numeric(determinant(S)$modulus) - (nu + k) / 2 * log1p(d2 / nu))
}

# dopasuj_nu: estymacja nu metodą największej wiarygodności (profil po nu przy mu, Sigma z próby).
dopasuj_nu <- function(R, przedzial = c(2.05, 100)) {
  R <- as.matrix(R)
  optimize(loglik_mvt, przedzial, R = R, mu = colMeans(R), Sigma = cov(R),
           maximum = TRUE)$maximum
}
