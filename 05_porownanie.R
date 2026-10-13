# 05. Porównanie metod: tabela zbiorcza, wrażliwość na poziom ufności i horyzont.

source("R/config.R")

d  <- readRDS(SCIEZKA_ZWR)
R  <- coredata(d$R)
rp <- as.numeric(R %*% WAGI)
mu <- colMeans(R); Sigma <- cov(R)
mu_p <- mean(rp); sd_p <- sd(rp)
w03 <- readRDS("output/03_var_klasyczny.rds")
w04 <- readRDS("output/04_monte_carlo.rds")
nu  <- w04$nu
METODY <- c("Historyczna", "Normalna", "MC normalny", "MC t-Student", "MC bootstrap")

# --- Tabela zbiorcza -------------------------------------------------------
tab <- rbind(w03$wyniki, w04$wyniki)
tab$metoda <- factor(tab$metoda, METODY)
sprawdz_spojnosc(tab)
tab <- tab[order(tab$horyzont, tab$alfa, tab$metoda), ]
# Tabela szeroka: wiersz = metoda, kolumny = VaR/CVaR dla każdej pary (alfa, horyzont).
szeroka <- do.call(cbind, unname(lapply(split(tab, list(tab$alfa, tab$horyzont)), function(s) {
  s <- s[order(s$metoda), ]
  x <- data.frame(round(s$VaR), round(s$CVaR))
  names(x) <- paste0(c("VaR_", "CVaR_"), s$alfa[1] * 100, "_", s$horyzont[1], "d")
  x
})))
szeroka <- data.frame(metoda = METODY, szeroka)
print(szeroka)
write.csv(szeroka, file.path(DIR_TAB, "05_porownanie.csv"), row.names = FALSE)
write.csv(tab, file.path(DIR_TAB, "05_porownanie_dlugie.csv"), row.names = FALSE)

# Wykres słupkowy: VaR i CVaR wszystkich metod w 4 wariantach (alfa x horyzont)
tl <- rbind(data.frame(tab[, c("metoda", "alfa", "horyzont")], miara = "VaR", wartosc = tab$VaR),
            data.frame(tab[, c("metoda", "alfa", "horyzont")], miara = "CVaR", wartosc = tab$CVaR))
tl$panel <- factor(sprintf("%g%%, %d %s", tl$alfa * 100, tl$horyzont, ifelse(tl$horyzont == 1, "dzień", "dni")),
                   levels = c("95%, 1 dzień", "99%, 1 dzień", "95%, 10 dni", "99%, 10 dni"))
tl$miara <- factor(tl$miara, c("VaR", "CVaR"))
p_bar <- ggplot(tl, aes(metoda, wartosc, fill = metoda, alpha = miara)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  facet_wrap(~panel, scales = "free_y") +
  scale_fill_manual(values = KOLORY_METOD, guide = "none") +
  scale_alpha_manual(values = c(VaR = 0.55, CVaR = 1)) +
  scale_y_continuous(labels = function(x) fmt_pln(x), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "VaR i CVaR portfela 100 000 PLN: porównanie metod",
       subtitle = "Jaśniejszy słupek: VaR; ciemniejszy: CVaR. Różnice metod rosną w ogonie (99%, CVaR)",
       x = NULL, y = "PLN") +
  theme_projekt() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
zapisz_wykres(p_bar, "05_porownanie", w = 10, h = 6.5)

# --- Wrażliwość na poziom ufności (1 dzień) --------------------------------
alfy_s <- seq(0.90, 0.995, by = 0.005)
set.seed(SEED + 10)
sym <- list("MC normalny" = sim_mc_normal(N_SIM, mu, Sigma),
            "MC t-Student" = sim_mc_t(N_SIM, mu, Sigma, nu),
            "MC bootstrap" = sim_bootstrap(N_SIM, R))
wa <- do.call(rbind, lapply(alfy_s, function(a) rbind(
  data.frame(metoda = "Historyczna", alfa = a, t(var_cvar_probka(rp, a))),
  data.frame(metoda = "Normalna", alfa = a, t(var_cvar_normalny(mu_p, sd_p, a))),
  do.call(rbind, lapply(names(sym), function(m) data.frame(metoda = m, alfa = a, t(var_cvar_probka(sym[[m]], a)))))
)))
wa$metoda <- factor(wa$metoda, METODY)
wal <- rbind(data.frame(wa[, 1:2], miara = "VaR", y = wa$VaR), data.frame(wa[, 1:2], miara = "CVaR", y = wa$CVaR))
wal$miara <- factor(wal$miara, c("VaR", "CVaR"))
p_alfa <- ggplot(wal, aes(alfa, y, colour = metoda)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~miara) +
  scale_colour_manual(values = KOLORY_METOD) +
  scale_x_continuous(labels = percent) +
  scale_y_continuous(labels = function(x) fmt_pln(x)) +
  labs(title = "Wrażliwość na poziom ufności (horyzont 1 dzień)",
       subtitle = "Przy 90-95% metody są zgodne; powyżej ~97,5% rozkład normalny odstaje w dół",
       x = "Poziom ufności", y = "PLN") +
  theme_projekt()
zapisz_wykres(p_alfa, "05_wrazliwosc_alfa", w = 10, h = 4.8)
write.csv(wa, file.path(DIR_TAB, "05_wrazliwosc_alfa.csv"), row.names = FALSE)

# --- Wrażliwość na horyzont (VaR 99%, 1..20 dni) ---------------------------
H <- 1:20
n_h <- 50000
set.seed(SEED + 11)
# Dzienne scenariusze portfela (n_h x 20), skumulowane po kolumnach -> zwroty 1..20-dniowe.
dzienne <- list(
  "MC normalny"  = matrix(sim_mc_normal(n_h * max(H), mu, Sigma), n_h),
  "MC t-Student" = matrix(sim_mc_t(n_h * max(H), mu, Sigma, nu), n_h),
  "MC bootstrap" = matrix(sim_bootstrap(n_h * max(H), R), n_h))
skum <- lapply(dzienne, function(m) { for (j in 2:ncol(m)) m[, j] <- m[, j - 1] + m[, j]; m })
wh <- do.call(rbind, lapply(H, function(h) rbind(
  data.frame(metoda = "Historyczna", horyzont = h, t(var_cvar_probka(sumy_kroczace(rp, h), 0.99))),
  data.frame(metoda = "Normalna", horyzont = h, t(var_cvar_normalny(mu_p, sd_p, 0.99, h))),
  do.call(rbind, lapply(names(skum), function(m) data.frame(metoda = m, horyzont = h,
                                                             t(var_cvar_probka(skum[[m]][, h], 0.99)))))
)))
wh$metoda <- factor(wh$metoda, METODY)
reg_sqrt <- data.frame(horyzont = H, y = var_cvar_probka(rp, 0.99)[["VaR"]] * sqrt(H))
p_h <- ggplot(wh, aes(horyzont, VaR, colour = metoda)) +
  geom_line(data = reg_sqrt, aes(horyzont, y), inherit.aes = FALSE, linetype = "dotted", colour = SZARY) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.2) +
  scale_colour_manual(values = KOLORY_METOD) +
  scale_y_continuous(labels = function(x) fmt_pln(x)) +
  labs(title = "Wrażliwość VaR 99% na horyzont",
       subtitle = "Kropkowana: historyczny 1-dniowy x sqrt(h). Historyczny wielodniowy rośnie szybciej niż metody iid",
       x = "Horyzont (dni sesyjne)", y = "VaR 99% [PLN]") +
  theme_projekt()
zapisz_wykres(p_h, "05_wrazliwosc_horyzont", w = 9, h = 5)
write.csv(wh, file.path(DIR_TAB, "05_wrazliwosc_horyzont.csv"), row.names = FALSE)

saveRDS(list(tab = tab, szeroka = szeroka, wa = wa, wh = wh), "output/05_porownanie.rds")
message("OK: porównanie")
