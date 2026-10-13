# Miary ryzyka: VaR i CVaR (Expected Shortfall) w PLN.
# Konwencja: VaR i CVaR to DODATNIE kwoty straty. Strata = V * (1 - exp(r)),
# gdzie r to logarytmiczny zwrot portfela w danym horyzoncie.
# Zwrot portfela przybliżamy jako r_p = w' r (średnia ważona zwrotów log aktywów).

# strata_pln: zamiana zwrotu logarytmicznego portfela na stratę w PLN.
strata_pln <- function(r, V = WARTOSC) V * (1 - exp(r))

# zwrot_portfela: dzienne zwroty log portfela (wektor lub xts) przy wagach w.
zwrot_portfela <- function(R, w = WAGI) {
  if (is.xts(R)) return(xts(as.numeric(coredata(R) %*% w), index(R)))
  as.numeric(R %*% w)
}

# sumy_kroczace: h-dniowe zwroty log jako sumy zachodzących na siebie okien.
sumy_kroczace <- function(r, h) {
  r <- as.numeric(r)
  if (h == 1) return(r)
  cs <- cumsum(c(0, r))
  cs[(h + 1):length(cs)] - cs[1:(length(cs) - h)]
}

# var_cvar_probka: VaR i CVaR z próbki zwrotów (historycznej lub symulowanej).
# VaR = strata odpowiadająca kwantylowi (1 - alfa); CVaR = średnia strata w ogonie poza VaR.
var_cvar_probka <- function(r, alfa, V = WARTOSC) {
  q <- quantile(r, 1 - alfa, names = FALSE, type = 7)
  c(VaR = strata_pln(q, V), CVaR = mean(strata_pln(r[r <= q], V)))
}

# var_cvar_normalny: VaR i CVaR przy r ~ N(m, s^2), m = h*mu, s = sigma*sqrt(h).
# CVaR liczony dokładnie dla straty V(1 - exp(r)) (warunkowa średnia rozkładu lognormalnego).
var_cvar_normalny <- function(mu, sigma, alfa, h = 1, V = WARTOSC) {
  m <- mu * h; s <- sigma * sqrt(h); p <- 1 - alfa
  q <- m + qnorm(p) * s
  e_exp <- exp(m + s^2 / 2) * pnorm((q - m - s^2) / s) / p   # E[exp(r) | r <= q]
  c(VaR = strata_pln(q, V), CVaR = V * (1 - e_exp))
}

# var_cvar_t: VaR i CVaR dla r = m + s * sqrt((nu-2)/nu) * T_nu (t-Student o wariancji s^2).
# Używane do kontroli MC-t: przy wielowymiarowym t ze wspólnym nu portfel ma dokładnie ten rozkład.
var_cvar_t <- function(m, s, nu, alfa, V = WARTOSC) {
  sk <- s * sqrt((nu - 2) / nu); p <- 1 - alfa
  q  <- m + sk * qt(p, nu)
  gest <- function(x) dt((x - m) / sk, nu) / sk
  e_exp <- integrate(function(x) exp(x) * gest(x), -Inf, q)$value / p
  c(VaR = strata_pln(q, V), CVaR = V * (1 - e_exp))
}

# statystyki_opisowe: średnia, zmienność, skośność, kurtoza nadwyżkowa, test Jarque-Bera.
statystyki_opisowe <- function(x) {
  x <- as.numeric(x); n <- length(x); z <- (x - mean(x)) / sd(x)
  sk <- mean(z^3); ku <- mean(z^4) - 3
  jb <- n / 6 * (sk^2 + ku^2 / 4)
  c(n = n, srednia_dzienna = mean(x), sd_dzienne = sd(x),
    zmiennosc_roczna = sd(x) * sqrt(252), min = min(x), max = max(x),
    skosnosc = sk, kurtoza_nadwyzkowa = ku, JB_p = pchisq(jb, 2, lower.tail = FALSE))
}

# sprawdz_spojnosc: kontrola sensowności tabeli wyników (kolumny metoda, alfa, horyzont, VaR, CVaR).
# Zatrzymuje skrypt, jeśli VaR99 < VaR95 albo CVaR < VaR (z tolerancją numeryczną 1 PLN).
sprawdz_spojnosc <- function(tab, tol = 1) {
  stopifnot(all(tab$CVaR >= tab$VaR - tol))
  kl <- split(tab, list(tab$metoda, tab$horyzont), drop = TRUE)
  for (k in kl) {
    k <- k[order(k$alfa), ]
    stopifnot(all(diff(k$VaR) >= -tol), all(diff(k$CVaR) >= -tol))
  }
  invisible(TRUE)
}
