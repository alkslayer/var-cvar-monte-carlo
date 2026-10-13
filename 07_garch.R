# 07. Rozszerzenie: GARCH(1,1) (rugarch) ze zmienną zmiennością vs modele ze stałą zmiennością.
# Backtest na tych samych dniach co w 06, okno kroczące 500 dni, reestymacja co 20 dni.

source("R/config.R")
suppressPackageStartupMessages(library(rugarch))

d  <- readRDS(SCIEZKA_ZWR)
rp <- zwrot_portfela(d$R)
# GARCH estymujemy na zwrotach w procentach (x100): przy zwrotach rzędu 0,01 optymalizator
# rugarch często nie odwraca hesjanu. Wyniki (kwantyle, sigma) dzielimy potem przez 100.
rp100 <- rp * 100
bt <- readRDS("output/06_backtest_prognozy.rds")

spec_garch <- function(rozklad, model = "sGARCH", srednia = TRUE) ugarchspec(
  variance.model = list(model = model, garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = srednia),
  distribution.model = rozklad)

# --- Backtest kroczący (ugarchroll) ----------------------------------------
# n.start = 500: pierwsza prognoza dotyczy obserwacji 501, tak jak w 06.
# refit.window = "moving": estymacja tylko na ostatnich 500 dniach PRZED dniem prognozy.
kroczacy <- function(rozklad, model = "sGARCH", srednia = TRUE) {
  # suppressWarnings: rugarch ostrzega o nieodwróconym hesjanie (dotyczy błędów standardowych,
  # nie samych oszacowań); zbieżność sprawdzamy jawnie niżej.
  roll <- suppressWarnings(ugarchroll(spec_garch(rozklad, model, srednia), data = rp100, n.start = OKNO_BT,
                     refit.every = REFIT_CO, refit.window = "moving", window.size = OKNO_BT,
                     solver = "hybrid", calculate.VaR = TRUE, VaR.alpha = 1 - ALFY, keep.coef = TRUE))
  # Okna, w których estymacja się nie zbiegła, liczymy ponownie (resume).
  proby <- 0
  while (!is.null(convergence(roll)) && convergence(roll) != 0 && proby < 3) {
    roll <- resume(roll, solver = "hybrid"); proby <- proby + 1
  }
  stopifnot(convergence(roll) == 0)
  roll
}
set.seed(SEED + 30)
t0 <- Sys.time()
roll_n <- kroczacy("norm")
roll_t <- kroczacy("std")
message("GARCH kroczący policzony w ", round(difftime(Sys.time(), t0, units = "secs")), " s")

v_n <- as.data.frame(roll_n, which = "VaR"); v_t <- as.data.frame(roll_t, which = "VaR")
daty_g <- as.Date(rownames(v_n))
# Zgodność dat i zwrotów z backtestem 06 (te same dni, te same realizacje).
stopifnot(identical(daty_g, bt$data), max(abs(v_n$realized / 100 - bt$r)) < 1e-12)
bt$garchn_95 <- v_n[, 1] / 100; bt$garchn_99 <- v_n[, 2] / 100
bt$garcht_95 <- v_t[, 1] / 100; bt$garcht_99 <- v_t[, 2] / 100
bt$sigma_garch_t <- as.data.frame(roll_t)$Sigma / 100
saveRDS(bt, "output/07_backtest_prognozy.rds")

# Kontrola braku look-ahead: współczynniki dla prognozy na dzień t estymowane
# na oknie kończącym się przed t (rugarch zapisuje zakres estymacji każdego refitu).
koef <- roll_t@model$coef
okno_ok <- vapply(seq_along(koef), function(k) {
  koniec <- as.Date(koef[[k]]$index)          # data ostatniej obserwacji w oknie estymacji
  pierwsza_prognoza <- bt$data[(k - 1) * REFIT_CO + 1]
  koniec < pierwsza_prognoza
}, logical(1))
stopifnot(all(okno_ok))
message("Kontrola look-ahead GARCH: OK (", length(koef), " reestymacji)")

# --- Testy: wszystkie metody -----------------------------------------------
metody_bt <- c(hist = "Historyczna", norm = "Normalna", mct = "MC t-Student",
               garchn = "GARCH-norm", garcht = "GARCH-t")
testy <- do.call(rbind, lapply(names(metody_bt), function(m) do.call(rbind, lapply(ALFY, function(a) {
  q <- bt[[paste0(m, "_", a * 100)]]
  data.frame(metoda = metody_bt[[m]], alfa = a, t(podsumuj_backtest(bt$r, q, a)),
             sredni_VaR_PLN = mean(strata_pln(q)))
}))))
for (i in seq_len(nrow(testy))) {   # kontrola z rugarch::VaRTest
  m <- names(metody_bt)[metody_bt == testy$metoda[i]]
  vt <- VaRTest(1 - testy$alfa[i], bt$r, bt[[paste0(m, "_", testy$alfa[i] * 100)]])
  stopifnot(vt$actual.exceed == testy$przekroczenia[i], abs(vt$cc.LRstat - testy$LR_cc[i]) < 1e-6)
}
print(testy, digits = 3)
write.csv(testy, file.path(DIR_TAB, "07_backtest_wszystkie.csv"), row.names = FALSE)

# Marzec 2020 dla wszystkich metod
w_stresie <- bt$data >= STRES_OD & bt$data <= STRES_DO
stres <- do.call(rbind, lapply(names(metody_bt), function(m) do.call(rbind, lapply(ALFY, function(a) {
  q <- bt[[paste0(m, "_", a * 100)]][w_stresie]; r <- bt$r[w_stresie]
  data.frame(metoda = metody_bt[[m]], alfa = a, dni = length(r), przekroczenia = sum(r < q),
             oczekiwane = round(length(r) * (1 - a), 2),
             p_dwumian = pbinom(sum(r < q) - 1, length(r), 1 - a, lower.tail = FALSE),
             VaR_srednio_PLN = mean(strata_pln(q)))
}))))
print(stres, digits = 3)
write.csv(stres, file.path(DIR_TAB, "07_marzec2020_wszystkie.csv"), row.names = FALSE)

# Przekroczenia w podziale na lata (VaR 99%)
bt$rok <- format(bt$data, "%Y")
lata <- do.call(rbind, lapply(names(metody_bt), function(m)
  data.frame(metoda = metody_bt[[m]], aggregate(list(przekroczenia = bt$r < bt[[paste0(m, "_99")]]),
                                                list(rok = bt$rok), sum))))
lata_sz <- reshape(lata, idvar = "rok", timevar = "metoda", direction = "wide")
names(lata_sz) <- sub("przekroczenia.", "", names(lata_sz))
lata_sz$dni <- as.integer(table(bt$rok)[lata_sz$rok])
print(lata_sz)
write.csv(lata_sz, file.path(DIR_TAB, "07_przekroczenia_lata.csv"), row.names = FALSE)

# --- Odporność: inne specyfikacje GARCH-t ----------------------------------
# Model główny (sGARCH ze stałą w średniej) wybrany z góry; tu sprawdzamy, czy wnioski
# zależą od asymetrii (GJR: spadki zwiększają zmienność bardziej niż wzrosty) i od średniej.
warianty <- expand.grid(model = c("sGARCH", "gjrGARCH"), srednia = c(TRUE, FALSE), stringsAsFactors = FALSE)
odpornosc <- do.call(rbind, lapply(seq_len(nrow(warianty)), function(i) {
  v <- as.data.frame(kroczacy("std", warianty$model[i], warianty$srednia[i]), which = "VaR")
  do.call(rbind, lapply(seq_along(ALFY), function(k) {
    s <- podsumuj_backtest(v$realized / 100, v[, k] / 100, ALFY[k])
    data.frame(model = warianty$model[i], srednia_w_modelu = warianty$srednia[i], alfa = ALFY[k],
               przekroczenia = s[["przekroczenia"]], oczekiwane = s[["oczekiwane"]],
               p_kupiec = s[["p_uc"]], p_christoffersen_cc = s[["p_cc"]])
  }))
}))
print(odpornosc, digits = 3)
write.csv(odpornosc, file.path(DIR_TAB, "07_garch_odpornosc.csv"), row.names = FALSE)

# --- Prognoza GARCH-t na najbliższy dzień (ostatnie 500 dni) ---------------
fit <- suppressWarnings(ugarchfit(spec_garch("std"), data = tail(rp100, OKNO_BT), solver = "hybrid"))
wsp <- coef(fit)
print(round(wsp, 5))
fc <- ugarchforecast(fit, n.ahead = 1)
m1 <- as.numeric(fitted(fc)) / 100; s1 <- as.numeric(sigma(fc)) / 100; nu_g <- wsp[["shape"]]
garch_1d <- do.call(rbind, lapply(ALFY, function(a)
  data.frame(metoda = "GARCH-t", alfa = a, horyzont = 1, t(var_cvar_t(m1, s1, nu_g, a)))))
# Horyzont 10 dni: symulacja ścieżek GARCH z bieżącego stanu (zmienność ewoluuje w czasie).
set.seed(SEED + 31)
sim10 <- ugarchsim(fit, n.sim = 10, m.sim = 50000, startMethod = "sample")
r10 <- colSums(fitted(sim10)) / 100
garch_10d <- do.call(rbind, lapply(ALFY, function(a)
  data.frame(metoda = "GARCH-t", alfa = a, horyzont = 10, t(var_cvar_probka(r10, a)))))
garch_now <- rbind(garch_1d, garch_10d)
sprawdz_spojnosc(garch_now)
garch_now$sigma_dzienna_teraz <- s1
garch_now$sigma_dzienna_dlugookresowa <- sqrt(wsp[["omega"]] / (1 - wsp[["alpha1"]] - wsp[["beta1"]])) / 100
print(garch_now)
write.csv(garch_now, file.path(DIR_TAB, "07_garch_prognoza.csv"), row.names = FALSE)
write.csv(data.frame(parametr = names(wsp), wartosc = wsp), file.path(DIR_TAB, "07_garch_parametry.csv"),
          row.names = FALSE)

# --- Wykresy ---------------------------------------------------------------
dl <- rbind(data.frame(data = bt$data, metoda = "Historyczna", q = bt$hist_99),
            data.frame(data = bt$data, metoda = "GARCH-t", q = bt$garcht_99))
prz <- rbind(data.frame(data = bt$data, r = bt$r, metoda = "Historyczna")[bt$r < bt$hist_99, ],
             data.frame(data = bt$data, r = bt$r, metoda = "GARCH-t")[bt$r < bt$garcht_99, ])
dl$metoda <- factor(dl$metoda, c("Historyczna", "GARCH-t"))
prz$metoda <- factor(prz$metoda, c("Historyczna", "GARCH-t"))
p_g <- ggplot() +
  geom_col(data = bt, aes(data, -strata_pln(r)), fill = "#c9cbc4", width = 1) +
  geom_line(data = dl, aes(data, -strata_pln(q), colour = metoda), linewidth = 0.5) +
  geom_point(data = prz, aes(data, -strata_pln(r), colour = metoda), size = 1.4) +
  facet_wrap(~metoda, ncol = 1) +
  scale_colour_manual(values = KOLORY_METOD, guide = "none") +
  scale_y_continuous(labels = function(x) fmt_pln(x)) +
  labs(title = "VaR 99%: stała zmienność (historyczny) vs GARCH-t",
       subtitle = "GARCH podnosi VaR zaraz po pierwszych dużych spadkach, więc przekroczenia nie tworzą skupisk",
       x = NULL, y = "PLN") +
  theme_projekt()
zapisz_wykres(p_g, "07_garch_vs_hist", w = 10, h = 6.5)

okno_2020 <- bt$data >= as.Date("2020-02-01") & bt$data <= as.Date("2020-05-31")
dl5 <- do.call(rbind, lapply(names(metody_bt), function(m)
  data.frame(data = bt$data, metoda = metody_bt[[m]], q = bt[[paste0(m, "_99")]])))
dl5$metoda <- factor(dl5$metoda, metody_bt)
p_20 <- ggplot() +
  geom_col(data = bt[okno_2020, ], aes(data, -strata_pln(r)), fill = "#c9cbc4") +
  geom_line(data = dl5[dl5$data %in% bt$data[okno_2020], ], aes(data, -strata_pln(q), colour = metoda),
            linewidth = 0.8) +
  scale_colour_manual(values = KOLORY_METOD) +
  scale_y_continuous(labels = function(x) fmt_pln(x)) +
  labs(title = "Marzec 2020: VaR 99% wszystkich metod",
       subtitle = "GARCH reaguje w kilka dni, modele ze stałą zmiennością w oknie 500 dni prawie wcale",
       x = NULL, y = "PLN") +
  theme_projekt()
zapisz_wykres(p_20, "07_marzec2020", w = 9, h = 5)

saveRDS(list(testy = testy, stres = stres, lata = lata_sz, garch_now = garch_now, wsp = wsp,
             odpornosc = odpornosc),
        "output/07_garch.rds")
message("OK: GARCH")
