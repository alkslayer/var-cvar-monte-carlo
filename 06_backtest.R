# 06. Backtest VaR 1-dniowego: okno kroczące 500 dni, testy Kupca i Christoffersena,
# osobna ocena okresu stresu (marzec 2020).

source("R/config.R")
suppressPackageStartupMessages(library(rugarch))  # tylko do niezależnej kontroli testów (VaRTest)

d  <- readRDS(SCIEZKA_ZWR)
R  <- coredata(d$R)
daty <- index(d$R)
rp <- as.numeric(R %*% WAGI)
dni  <- (OKNO_BT + 1):nrow(R)          # dni prognozy (dzień t korzysta z dni t-500..t-1)

# --- Prognozy kroczące -----------------------------------------------------
set.seed(SEED + 20)
t0 <- Sys.time()
q_hist <- prognozy_kroczace(R, f_historyczna)
q_norm <- prognozy_kroczace(R, f_normalna)
q_mct  <- prognozy_kroczace(R, f_mc_t)
message("Prognozy policzone w ", round(difftime(Sys.time(), t0, units = "secs")), " s")

# --- Kontrola braku look-ahead bias ----------------------------------------
# Psujemy wszystkie dane od dnia t_test (mnożymy przez 10). Prognozy na dni <= t_test
# muszą zostać identyczne, bo korzystają tylko z danych sprzed dnia prognozy.
t_test <- 1500
R_zepsute <- R; R_zepsute[t_test:nrow(R), ] <- R_zepsute[t_test:nrow(R), ] * 10
j_ok <- which(dni <= t_test)
stopifnot(identical(prognozy_kroczace(R_zepsute, f_historyczna)[j_ok, ], q_hist[j_ok, ]),
          identical(prognozy_kroczace(R_zepsute, f_normalna)[j_ok, ], q_norm[j_ok, ]))
# Dzień t_test + 1 musi się już zmienić (okno zawiera zepsuty dzień t_test).
stopifnot(!identical(prognozy_kroczace(R_zepsute, f_normalna)[j_ok[length(j_ok)] + 1, ],
                     q_norm[j_ok[length(j_ok)] + 1, ]))
message("Kontrola look-ahead: OK")

bt <- data.frame(data = daty[dni], r = rp[dni],
                 hist_95 = q_hist[, 1], hist_99 = q_hist[, 2],
                 norm_95 = q_norm[, 1], norm_99 = q_norm[, 2],
                 mct_95  = q_mct[, 1],  mct_99  = q_mct[, 2])
saveRDS(bt, "output/06_backtest_prognozy.rds")

# --- Testy -----------------------------------------------------------------
metody_bt <- c(hist = "Historyczna", norm = "Normalna", mct = "MC t-Student")
testy <- do.call(rbind, lapply(names(metody_bt), function(m) do.call(rbind, lapply(ALFY, function(a) {
  q <- bt[[paste0(m, "_", a * 100)]]
  data.frame(metoda = metody_bt[[m]], alfa = a, t(podsumuj_backtest(bt$r, q, a)))
}))))
print(testy, digits = 3)

# Niezależna kontrola liczby przekroczeń i statystyk LR funkcją rugarch::VaRTest.
for (i in seq_len(nrow(testy))) {
  m <- names(metody_bt)[metody_bt == testy$metoda[i]]
  q <- bt[[paste0(m, "_", testy$alfa[i] * 100)]]
  vt <- VaRTest(1 - testy$alfa[i], bt$r, q)
  stopifnot(vt$actual.exceed == testy$przekroczenia[i],
            abs(vt$uc.LRstat - testy$LR_uc[i]) < 1e-6,
            abs(vt$cc.LRstat - testy$LR_cc[i]) < 1e-6)
}
message("Kontrola testów z rugarch::VaRTest: OK")
write.csv(testy, file.path(DIR_TAB, "06_backtest_testy.csv"), row.names = FALSE)

# --- Okres stresu: marzec 2020 ---------------------------------------------
w_stresie <- bt$data >= STRES_OD & bt$data <= STRES_DO
stres <- do.call(rbind, lapply(names(metody_bt), function(m) do.call(rbind, lapply(ALFY, function(a) {
  q <- bt[[paste0(m, "_", a * 100)]][w_stresie]; r <- bt$r[w_stresie]
  data.frame(metoda = metody_bt[[m]], alfa = a, dni = length(r), przekroczenia = sum(r < q),
             oczekiwane = round(length(r) * (1 - a), 2),
             p_dwumian = pbinom(sum(r < q) - 1, length(r), 1 - a, lower.tail = FALSE),
             VaR_srednio_PLN = mean(strata_pln(q)), max_strata_PLN = max(strata_pln(r)))
}))))
print(stres, digits = 3)
write.csv(stres, file.path(DIR_TAB, "06_backtest_marzec2020.csv"), row.names = FALSE)

# --- Wykresy ---------------------------------------------------------------
dl <- rbind(data.frame(data = bt$data, metoda = "Historyczna", q = bt$hist_99),
            data.frame(data = bt$data, metoda = "Normalna", q = bt$norm_99),
            data.frame(data = bt$data, metoda = "MC t-Student", q = bt$mct_99))
dl$metoda <- factor(dl$metoda, metody_bt)
prz <- bt[bt$r < bt$hist_99, ]
p_bt <- ggplot() +
  geom_col(data = bt, aes(data, -strata_pln(r)), fill = "#c9cbc4", width = 1) +
  geom_line(data = dl, aes(data, -strata_pln(q), colour = metoda), linewidth = 0.5) +
  geom_point(data = prz, aes(data, -strata_pln(r)), colour = "#0b0b0b", size = 1.2) +
  scale_colour_manual(values = KOLORY_METOD) +
  scale_y_continuous(labels = function(x) fmt_pln(x)) +
  labs(title = "Backtest VaR 99% (okno 500 dni): wynik dzienny portfela i prognozowany VaR",
       subtitle = "Szare słupki: zysk/strata dnia; linie: -VaR; czarne kropki: przekroczenia VaR historycznego",
       x = NULL, y = "PLN") +
  theme_projekt()
zapisz_wykres(p_bt, "06_backtest", w = 10, h = 5)

okno_2020 <- bt$data >= as.Date("2020-02-01") & bt$data <= as.Date("2020-05-31")
p_20 <- ggplot() +
  geom_col(data = bt[okno_2020, ], aes(data, -strata_pln(r)), fill = "#c9cbc4") +
  geom_line(data = dl[dl$data %in% bt$data[okno_2020], ], aes(data, -strata_pln(q), colour = metoda),
            linewidth = 0.8) +
  scale_colour_manual(values = KOLORY_METOD) +
  scale_y_continuous(labels = function(x) fmt_pln(x)) +
  labs(title = "Stres COVID: luty-maj 2020, VaR 99%",
       subtitle = "Modele ze stałą zmiennością w oknie 500 dni reagują z opóźnieniem tygodni",
       x = NULL, y = "PLN") +
  theme_projekt()
zapisz_wykres(p_20, "06_backtest_2020", w = 9, h = 4.8)
message("OK: backtest")
