# 04. Monte Carlo: skorelowane zwroty (Cholesky), rozkład normalny, t-Student, bootstrap.
# Zbieżność VaR/CVaR względem liczby symulacji.

source("R/config.R")
suppressPackageStartupMessages(library(MASS))

d  <- readRDS(SCIEZKA_ZWR)
R  <- coredata(d$R)
mu <- colMeans(R); Sigma <- cov(R)
rp <- as.numeric(R %*% WAGI)
mu_p <- sum(WAGI * mu); sd_p <- sqrt(as.numeric(t(WAGI) %*% Sigma %*% WAGI))

# --- Stopnie swobody t-Studenta z danych -----------------------------------
nu <- dopasuj_nu(R)
# Dla porównania: nu z jednowymiarowego dopasowania t (MLE, MASS::fitdistr) dla każdej spółki i portfela.
nu_jedn <- sapply(c(asplit(R, 2), list(Portfel = rp)), function(x)
  suppressWarnings(fitdistr(x / sd(x), "t"))$estimate[["df"]])
nu_tab <- data.frame(seria = c("Wielowymiarowy t (5 spółek, MLE)", NAZWY[colnames(R)], "Portfel"),
                     nu = round(c(nu, nu_jedn), 2))
print(nu_tab)
write.csv(nu_tab, file.path(DIR_TAB, "04_stopnie_swobody.csv"), row.names = FALSE)

# Kontrola dopasowania: profil log-wiarygodności po nu
siatka_nu <- seq(2.5, 15, by = 0.25)
ll <- sapply(siatka_nu, loglik_mvt, R = R, mu = mu, Sigma = Sigma)
p_nu <- ggplot(data.frame(nu = siatka_nu, ll = ll - max(ll)), aes(nu, ll)) +
  geom_line(colour = PALETA[1], linewidth = 0.7) +
  geom_vline(xintercept = nu, linetype = "dashed", colour = SZARY) +
  labs(title = "Profil log-wiarygodności wielowymiarowego t względem nu",
       subtitle = sprintf("Maksimum przy nu = %.2f", nu), x = "nu", y = "logL - max") +
  theme_projekt()
zapisz_wykres(p_nu, "04_profil_nu", w = 7, h = 4)

# --- Główne wyniki MC (N_SIM symulacji, set.seed) --------------------------
symuluj <- function(metoda, n, h) switch(metoda,
  "MC normalny"  = sim_mc_normal(n, mu, Sigma, h = h),
  "MC t-Student" = sim_mc_t(n, mu, Sigma, nu, h = h),
  "MC bootstrap" = sim_bootstrap(n, R, h = h))
metody_mc <- c("MC normalny", "MC t-Student", "MC bootstrap")

set.seed(SEED)
wyniki <- do.call(rbind, lapply(metody_mc, function(m) do.call(rbind, lapply(HORYZONTY, function(h) {
  sim <- symuluj(m, N_SIM, h)
  do.call(rbind, lapply(ALFY, function(a)
    data.frame(metoda = m, alfa = a, horyzont = h, t(var_cvar_probka(sim, a)))))
}))))
sprawdz_spojnosc(wyniki)
wyniki$VaR_proc <- wyniki$VaR / WARTOSC; wyniki$CVaR_proc <- wyniki$CVaR / WARTOSC
print(transform(wyniki, VaR = round(VaR), CVaR = round(CVaR)))

# --- Kontrole poprawności: MC vs wzory analityczne -------------------------
# 1) MC normalny przy 10^6 symulacji musi trafić w VaR/CVaR parametryczny (h = 1 i 10).
# 2) MC t przy h = 1: portfel z wielowymiarowego t ze wspólnym nu ma dokładnie rozkład t,
#    więc znamy dokładną wartość VaR/CVaR.
set.seed(SEED + 1)
kontrola <- do.call(rbind, lapply(HORYZONTY, function(h) {
  sim <- sim_mc_normal(1e6, mu, Sigma, h = h)
  do.call(rbind, lapply(ALFY, function(a) {
    mc <- var_cvar_probka(sim, a); an <- var_cvar_normalny(mu_p, sd_p, a, h)
    data.frame(model = "normalny", alfa = a, horyzont = h, VaR_MC = mc[["VaR"]], VaR_wzor = an[["VaR"]],
               CVaR_MC = mc[["CVaR"]], CVaR_wzor = an[["CVaR"]])
  }))
}))
sim_t <- sim_mc_t(1e6, mu, Sigma, nu, h = 1)
kontrola <- rbind(kontrola, do.call(rbind, lapply(ALFY, function(a) {
  mc <- var_cvar_probka(sim_t, a); an <- var_cvar_t(mu_p, sd_p, nu, a)
  data.frame(model = "t-Student", alfa = a, horyzont = 1, VaR_MC = mc[["VaR"]], VaR_wzor = an[["VaR"]],
             CVaR_MC = mc[["CVaR"]], CVaR_wzor = an[["CVaR"]])
})))
kontrola$blad_VaR_proc  <- 100 * (kontrola$VaR_MC / kontrola$VaR_wzor - 1)
kontrola$blad_CVaR_proc <- 100 * (kontrola$CVaR_MC / kontrola$CVaR_wzor - 1)
print(kontrola)
stopifnot(abs(kontrola$blad_VaR_proc) < 1, abs(kontrola$blad_CVaR_proc) < 1.5)
write.csv(kontrola, file.path(DIR_TAB, "04_kontrola_mc.csv"), row.names = FALSE)
rm(sim_t)

# --- Zbieżność względem liczby symulacji -----------------------------------
# Dla każdej z 30 powtórek generujemy 10^6 scenariuszy i liczymy VaR/CVaR 99% (1 dzień)
# na pierwszych N z nich. Pętla jest tylko po powtórkach; same symulacje są wektorowe.
N_siatka <- round(10^seq(2, 6, by = 0.25))
powtorki <- 30
set.seed(SEED + 2)
zb <- do.call(rbind, lapply(metody_mc, function(m) do.call(rbind, lapply(seq_len(powtorki), function(k) {
  sim <- symuluj(m, max(N_siatka), 1)
  do.call(rbind, lapply(N_siatka, function(n) {
    v <- var_cvar_probka(sim[1:n], 0.99)
    data.frame(metoda = m, powtorka = k, N = n, VaR = v[["VaR"]], CVaR = v[["CVaR"]])
  }))
}))))
zb_sum <- do.call(rbind, lapply(split(zb, list(zb$metoda, zb$N)), function(s) data.frame(
  metoda = s$metoda[1], N = s$N[1],
  VaR_med = median(s$VaR), VaR_lo = quantile(s$VaR, 0.05), VaR_hi = quantile(s$VaR, 0.95),
  CVaR_med = median(s$CVaR), CVaR_lo = quantile(s$CVaR, 0.05), CVaR_hi = quantile(s$CVaR, 0.95),
  VaR_bladwzgl = sd(s$VaR) / mean(s$VaR))))
write.csv(zb_sum, file.path(DIR_TAB, "04_zbieznosc.csv"), row.names = FALSE)

# Wartości referencyjne: wzór (normalny, t) i wartość z pełnych danych (bootstrap = historyczny).
ref <- data.frame(metoda = metody_mc,
                  VaR = c(var_cvar_normalny(mu_p, sd_p, 0.99)[["VaR"]], var_cvar_t(mu_p, sd_p, nu, 0.99)[["VaR"]],
                          var_cvar_probka(rp, 0.99)[["VaR"]]),
                  CVaR = c(var_cvar_normalny(mu_p, sd_p, 0.99)[["CVaR"]], var_cvar_t(mu_p, sd_p, nu, 0.99)[["CVaR"]],
                           var_cvar_probka(rp, 0.99)[["CVaR"]]))
zb_long <- rbind(
  data.frame(miara = "VaR 99%", zb_sum[, c("metoda", "N")], med = zb_sum$VaR_med, lo = zb_sum$VaR_lo, hi = zb_sum$VaR_hi),
  data.frame(miara = "CVaR 99%", zb_sum[, c("metoda", "N")], med = zb_sum$CVaR_med, lo = zb_sum$CVaR_lo, hi = zb_sum$CVaR_hi))
ref_long <- rbind(data.frame(miara = "VaR 99%", metoda = ref$metoda, y = ref$VaR),
                  data.frame(miara = "CVaR 99%", metoda = ref$metoda, y = ref$CVaR))
zb_long$miara <- factor(zb_long$miara, c("VaR 99%", "CVaR 99%"))
ref_long$miara <- factor(ref_long$miara, c("VaR 99%", "CVaR 99%"))
zb_long$metoda <- factor(zb_long$metoda, metody_mc); ref_long$metoda <- factor(ref_long$metoda, metody_mc)

p_zb <- ggplot(zb_long, aes(N, med)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = metoda), alpha = 0.25) +
  geom_line(aes(colour = metoda), linewidth = 0.7) +
  geom_hline(data = ref_long, aes(yintercept = y), linetype = "dashed", colour = "#0b0b0b") +
  facet_grid(miara ~ metoda, scales = "free_y") +
  scale_x_log10(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = function(x) fmt_pln(x)) +
  scale_colour_manual(values = KOLORY_METOD, guide = "none") +
  scale_fill_manual(values = KOLORY_METOD, guide = "none") +
  labs(title = "Zbieżność Monte Carlo: VaR i CVaR 99% (1 dzień) a liczba symulacji",
       subtitle = "Linia: mediana z 30 powtórek; pas: 5-95%. Przerywana: wartość dokładna (wzór lub pełna próba)",
       x = "Liczba symulacji N (skala log)", y = "PLN") +
  theme_projekt()
zapisz_wykres(p_zb, "04_zbieznosc", w = 10, h = 6)

# Błąd względny ~ 1/sqrt(N): wykres log-log
p_blad <- ggplot(zb_sum, aes(N, VaR_bladwzgl, colour = metoda)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.5) +
  geom_line(data = data.frame(N = N_siatka, VaR_bladwzgl = 0.6 / sqrt(N_siatka)),
            aes(N, VaR_bladwzgl), inherit.aes = FALSE, linetype = "dashed", colour = SZARY) +
  scale_x_log10(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_log10(labels = percent) +
  scale_colour_manual(values = KOLORY_METOD) +
  labs(title = "Błąd Monte Carlo VaR 99% maleje jak 1/sqrt(N)",
       subtitle = "Odchylenie standardowe z 30 powtórek / średnia; przerywana: nachylenie 1/sqrt(N)",
       x = "Liczba symulacji N", y = "Względny błąd standardowy") +
  theme_projekt()
zapisz_wykres(p_blad, "04_blad_mc", w = 8, h = 4.5)

saveRDS(list(wyniki = wyniki, nu = nu, nu_tab = nu_tab, kontrola = kontrola, zb_sum = zb_sum),
        "output/04_monte_carlo.rds")
write.csv(wyniki, file.path(DIR_TAB, "04_var_mc.csv"), row.names = FALSE)
message("OK: Monte Carlo, nu = ", round(nu, 2))
