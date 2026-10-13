# 02. Eksploracja: rozkłady, QQ-plot, grube ogony, autokorelacja, korelacje.

source("R/config.R")
suppressPackageStartupMessages(library(patchwork))

d  <- readRDS(SCIEZKA_ZWR)
R  <- d$R
rp <- zwrot_portfela(R)
dane <- cbind(coredata(R), Portfel = as.numeric(rp))
colnames(dane) <- c(NAZWY[colnames(R)], "Portfel")

# --- Statystyki opisowe ----------------------------------------------------
stat <- t(apply(dane, 2, statystyki_opisowe))
stat <- data.frame(seria = rownames(stat), round(stat, 5), row.names = NULL)
print(stat)
write.csv(stat, file.path(DIR_TAB, "02_statystyki.csv"), row.names = FALSE)

# --- Histogramy z gęstością normalną ---------------------------------------
dl <- data.frame(seria = factor(rep(colnames(dane), each = nrow(dane)), levels = colnames(dane)),
                 r = as.numeric(dane))
gest_n <- do.call(rbind, lapply(split(dl, dl$seria), function(s) {
  x <- seq(-0.12, 0.12, length.out = 400)
  data.frame(seria = s$seria[1], x = x, y = dnorm(x, mean(s$r), sd(s$r)))
}))
p_hist <- ggplot(dl, aes(r)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.004,
                 fill = PALETA[1], colour = "white", linewidth = 0.2) +
  geom_line(data = gest_n, aes(x, y), colour = PALETA[2], linewidth = 0.6) +
  facet_wrap(~seria, scales = "free_y") +
  coord_cartesian(xlim = c(-0.1, 0.1)) +
  scale_x_continuous(labels = percent) +
  labs(title = "Rozkład dziennych zwrotów log a rozkład normalny",
       subtitle = "Słupki: dane; linia: N(średnia, sd) z próby. Wyższy szczyt i grubsze ogony niż w normalnym",
       x = "Zwrot dzienny", y = "Gęstość") +
  theme_projekt()
zapisz_wykres(p_hist, "02_histogramy", w = 9, h = 6)

# --- QQ-ploty --------------------------------------------------------------
qq <- do.call(rbind, lapply(split(dl, dl$seria), function(s) {
  z <- sort((s$r - mean(s$r)) / sd(s$r))
  data.frame(seria = s$seria[1], teor = qnorm(ppoints(length(z))), emp = z)
}))
p_qq <- ggplot(qq, aes(teor, emp)) +
  geom_abline(colour = SZARY, linewidth = 0.4) +
  geom_point(size = 0.6, colour = PALETA[1], alpha = 0.6) +
  facet_wrap(~seria) +
  labs(title = "QQ-plot standaryzowanych zwrotów względem N(0,1)",
       subtitle = "Punkty poza przekątną w ogonach = grube ogony (więcej skrajnych dni niż w normalnym)",
       x = "Kwantyle teoretyczne N(0,1)", y = "Kwantyle empiryczne") +
  theme_projekt()
zapisz_wykres(p_qq, "02_qqplot", w = 9, h = 6)

# Liczba dni poniżej -3 sd: dane vs rozkład normalny
ogon <- data.frame(seria = colnames(dane),
                   dni_ponizej_3sd = colSums(scale(dane) < -3),
                   oczekiwane_normalny = round(nrow(dane) * pnorm(-3), 1))
print(ogon)
write.csv(ogon, file.path(DIR_TAB, "02_ogony.csv"), row.names = FALSE)

# --- Autokorelacja zwrotów i kwadratów zwrotów portfela --------------------
acf_df <- function(x, nazwa, lag.max = 30) {
  a <- acf(x, lag.max = lag.max, plot = FALSE)$acf[-1]
  data.frame(lag = seq_along(a), acf = a, seria = nazwa)
}
rpn <- as.numeric(rp)
acfs <- rbind(acf_df(rpn, "Zwroty portfela"), acf_df(rpn^2, "Kwadraty zwrotów portfela"))
acfs$seria <- factor(acfs$seria, levels = c("Zwroty portfela", "Kwadraty zwrotów portfela"))
pas <- 1.96 / sqrt(length(rpn))
p_acf <- ggplot(acfs, aes(lag, acf)) +
  geom_hline(yintercept = c(-pas, pas), linetype = "dashed", colour = SZARY) +
  geom_hline(yintercept = 0, colour = SZARY) +
  geom_segment(aes(xend = lag, yend = 0), colour = PALETA[1], linewidth = 0.8) +
  facet_wrap(~seria) +
  labs(title = "Autokorelacja (ACF) portfela",
       subtitle = "Zwroty: prawie brak pamięci. Kwadraty: silna, trwała autokorelacja = skupiska zmienności",
       x = "Opóźnienie (dni)", y = "ACF") +
  theme_projekt()
zapisz_wykres(p_acf, "02_acf", w = 9, h = 4)

lb <- data.frame(
  seria = c("zwroty", "kwadraty zwrotów"),
  LjungBox_lag10_p = c(Box.test(rpn, 10, "Ljung-Box")$p.value,
                       Box.test(rpn^2, 10, "Ljung-Box")$p.value),
  ACF_lag1 = c(acfs$acf[1], acfs$acf[31])
)
print(lb)
write.csv(lb, file.path(DIR_TAB, "02_ljungbox.csv"), row.names = FALSE)

# --- Macierz korelacji -----------------------------------------------------
kor <- cor(coredata(R)); dimnames(kor) <- list(NAZWY[colnames(R)], NAZWY[colnames(R)])
write.csv(round(kor, 3), file.path(DIR_TAB, "02_korelacje.csv"))
kd <- as.data.frame(as.table(kor))
p_kor <- ggplot(kd, aes(Var1, Var2, fill = Freq)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", Freq)), size = 3.5, colour = "#0b0b0b") +
  scale_fill_gradient(low = "#f1f6fd", high = "#2a78d6", limits = c(0, 1)) +
  labs(title = "Korelacje dziennych zwrotów (2015-2026)", x = NULL, y = NULL) +
  theme_projekt() + theme(legend.position = "right", panel.grid = element_blank())
zapisz_wykres(p_kor, "02_korelacje", w = 6.5, h = 5)

# Korelacja w okresach wysokiej i niskiej zmienności. Reżim wyznacza zmienność z POPRZEDNICH
# 21 dni (górny kwintyl = stres). Nie dzielimy wg dzisiejszego zwrotu portfela, bo selekcja
# po sumie zwrotów sztucznie zaniża korelację w wybranej podpróbie.
vol_pop <- as.numeric(stats::lag(zoo::rollapply(rp, 21, sd, align = "right", fill = NA), 1))
stres <- !is.na(vol_pop) & vol_pop > quantile(vol_pop, 0.8, na.rm = TRUE)
sr_kor <- function(m) { k <- cor(m); mean(k[upper.tri(k)]) }
kor_rezim <- data.frame(dni = c("wszystkie", "pozostałe 80% dni", "20% dni o najwyższej zmienności z ostatnich 21 sesji"),
                        srednia_korelacja = c(sr_kor(coredata(R)), sr_kor(coredata(R)[!stres, ]),
                                              sr_kor(coredata(R)[stres, ])))
print(kor_rezim)
write.csv(kor_rezim, file.path(DIR_TAB, "02_korelacje_rezim.csv"), row.names = FALSE)

# --- Zmienność krocząca portfela (kontekst dla backtestu) ------------------
vol <- na.omit(zoo::rollapply(rp, 21, sd, align = "right") * sqrt(252))
p_vol <- ggplot(data.frame(data = index(vol), vol = as.numeric(vol)), aes(data, vol)) +
  geom_line(colour = PALETA[1], linewidth = 0.5) +
  scale_y_continuous(labels = percent) +
  labs(title = "Zmienność portfela (21-dniowe okno, w skali roku)",
       subtitle = "Zmienność nie jest stała: skoki w 03.2020 i 02.2022", x = NULL, y = NULL) +
  theme_projekt()
zapisz_wykres(p_vol, "02_zmiennosc", w = 9, h = 4)
message("OK: eksploracja")
