# 03. VaR i CVaR: metoda historyczna i parametryczna (normalna), 95%/99%, 1 i 10 dni.

source("R/config.R")

d  <- readRDS(SCIEZKA_ZWR)
rp <- as.numeric(zwrot_portfela(d$R))
mu <- mean(rp); sigma <- sd(rp)

siatka <- expand.grid(alfa = ALFY, horyzont = HORYZONTY)

# Historyczna: kwantyl empiryczny. Dla 10 dni zachodzące na siebie 10-dniowe sumy
# (ok. 2900 obserwacji zamiast 293 rozłącznych; kosztem zależności między obserwacjami).
hist <- do.call(rbind, lapply(seq_len(nrow(siatka)), function(i) {
  a <- siatka$alfa[i]; h <- siatka$horyzont[i]
  data.frame(metoda = "Historyczna", alfa = a, horyzont = h,
             t(var_cvar_probka(sumy_kroczace(rp, h), a)))
}))

# Parametryczna normalna: r_h ~ N(h*mu, h*sigma^2) (zwroty dzienne niezależne).
norm <- do.call(rbind, lapply(seq_len(nrow(siatka)), function(i) {
  a <- siatka$alfa[i]; h <- siatka$horyzont[i]
  data.frame(metoda = "Normalna", alfa = a, horyzont = h,
             t(var_cvar_normalny(mu, sigma, a, h)))
}))

# Wariant kontrolny: historyczny 1-dniowy przeskalowany pierwiastkiem z czasu (reguła sqrt(10)).
hist_sqrt <- subset(hist, horyzont == 1)
hist_sqrt$metoda <- "Historyczna x sqrt(10)"; hist_sqrt$horyzont <- 10
hist_sqrt$VaR <- hist_sqrt$VaR * sqrt(10); hist_sqrt$CVaR <- hist_sqrt$CVaR * sqrt(10)

wyniki <- rbind(hist, norm)
sprawdz_spojnosc(wyniki)
wyniki$VaR_proc  <- wyniki$VaR / WARTOSC
wyniki$CVaR_proc <- wyniki$CVaR / WARTOSC
print(transform(wyniki, VaR = round(VaR), CVaR = round(CVaR)))

# Ile kosztuje przybliżenie r_p = w'r? Dokładny dzienny zwrot log portfela o stałych wagach:
# log(sum w_i * exp(r_i)). Porównujemy historyczny VaR 99% z obu definicji.
rp_dokl <- log(as.numeric(exp(coredata(d$R)) %*% WAGI))
roznica_aproksymacji <- var_cvar_probka(rp, 0.99)[["VaR"]] - var_cvar_probka(rp_dokl, 0.99)[["VaR"]]
message("Różnica VaR 99% (w'r vs dokładny zwrot portfela): ", round(roznica_aproksymacji, 1), " PLN")

saveRDS(list(wyniki = wyniki, hist_sqrt = hist_sqrt, mu = mu, sigma = sigma,
             roznica_aproksymacji = roznica_aproksymacji),
        "output/03_var_klasyczny.rds")
write.csv(rbind(wyniki[, 1:5], hist_sqrt), file.path(DIR_TAB, "03_var_klasyczny.csv"),
          row.names = FALSE)

# --- Wykres: rozkład strat 1-dniowych z zaznaczonym VaR i CVaR 99% ----------
straty <- strata_pln(rp)
v99 <- subset(wyniki, horyzont == 1 & alfa == 0.99)
linie <- data.frame(metoda = rep(v99$metoda, 2), miara = rep(c("VaR 99%", "CVaR 99%"), each = 2),
                    x = c(v99$VaR, v99$CVaR))
p <- ggplot(data.frame(s = straty), aes(s)) +
  geom_histogram(binwidth = 250, fill = "#c9dcf3", colour = "white", linewidth = 0.1) +
  geom_vline(data = linie, aes(xintercept = x, colour = metoda, linetype = miara), linewidth = 0.7) +
  scale_colour_manual(values = KOLORY_METOD) +
  scale_linetype_manual(values = c("VaR 99%" = "solid", "CVaR 99%" = "dashed")) +
  scale_x_continuous(labels = function(x) fmt_pln(x)) +
  coord_cartesian(xlim = c(-6000, 9000)) +
  labs(title = "Dzienne straty portfela 100 000 PLN (dodatnie = strata)",
       subtitle = "Normalny VaR 99% jest blisko historycznego, ale CVaR normalny wyraźnie niedoszacowuje ogona",
       x = "Strata dzienna [PLN]", y = "Liczba dni") +
  theme_projekt()
zapisz_wykres(p, "03_rozklad_strat", w = 9, h = 5)
message("OK: VaR klasyczny")
