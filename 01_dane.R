# 01. Dane: pobranie z Yahoo Finance, kontrola jakości, zwroty logarytmiczne.
# Surowe dane pobierane są tylko raz (data/raw/). Ponowne pobranie: REFRESH=1 Rscript 01_dane.R

source("R/config.R")

# --- Pobranie (lub odczyt lokalnej kopii) ----------------------------------
if (!file.exists(SCIEZKA_RAW) || Sys.getenv("REFRESH") == "1") {
  raw <- pobierz_yahoo(TICKERY, DATA_OD)
  attr(raw, "pobrano") <- Sys.time()
  saveRDS(raw, SCIEZKA_RAW)
  message("Pobrano dane z Yahoo: ", format(Sys.time()))
} else {
  raw <- readRDS(SCIEZKA_RAW)
  message("Używam lokalnej kopii z ", format(attr(raw, "pobrano")))
}

# --- Kontrola jakości surowych danych --------------------------------------
qc <- kontrola_jakosci(raw)
print(qc)
write.csv(qc, file.path(DIR_TAB, "qc_surowe.csv"), row.names = FALSE)

# Wspólne dni z zerowym wolumenem: święta GPW i dziury w danych (cena przepisana).
zera <- dni_wolumen_zero(raw)
# Czy w te dni cena się zmienia? (powinna być płaska, jeśli to sztuczne wiersze)
zmiana_ceny <- sapply(raw, function(x) {
  rc <- diff(log(x$Close))
  sum(abs(rc[index(rc) %in% zera]) > 1e-10, na.rm = TRUE)
})
stopifnot(all(zmiana_ceny == 0))
write.csv(data.frame(data = as.character(zera)), file.path(DIR_TAB, "qc_dni_wolumen_zero.csv"),
          row.names = FALSE)
message("Dni z zerowym wolumenem dla wszystkich spółek: ", length(zera))

# Splity zgłaszane przez Yahoo (Close jest już skorygowany o splity).
splity <- do.call(rbind, lapply(TICKERY, function(t) {
  s <- tryCatch(getSplits(t, from = DATA_OD), error = function(e) NULL)
  if (is.null(s) || !length(s) || all(is.na(s))) return(NULL)
  data.frame(ticker = t, data = as.character(index(s)), wspolczynnik = as.numeric(s))
}))
print(splity)
write.csv(if (is.null(splity)) data.frame() else splity,
          file.path(DIR_TAB, "qc_splity.csv"), row.names = FALSE)

# Korekty dywidendowe: różnica zwrotu z Adjusted i Close.
dyw <- dywidendy_i_splity(raw, DATA_OD, DATA_DO)
write.csv(dyw, file.path(DIR_TAB, "qc_korekty_adjusted.csv"), row.names = FALSE)
message("Dni z korektą Adjusted > 1 pp: ", nrow(dyw))

# --- Czyszczenie i zwroty --------------------------------------------------
ceny <- czyste_ceny(raw, DATA_OD, DATA_DO)
R    <- zwroty_log(ceny)
stopifnot(!anyNA(R), all(is.finite(R)), end(R) == DATA_DO)

luki <- luki_sesyjne(index(ceny))
print(luki)
write.csv(luki, file.path(DIR_TAB, "qc_luki.csv"), row.names = FALSE)

odst <- odstajace(R, k = 6)
print(odst)
write.csv(odst, file.path(DIR_TAB, "qc_odstajace.csv"), row.names = FALSE)

podsumowanie <- data.frame(
  od = as.character(start(R)), do = as.character(end(R)), dni_ze_zwrotem = nrow(R),
  usuniete_dni_wolumen_zero = sum(index(raw[[1]]) %in% zera &
                                    index(raw[[1]]) <= DATA_DO),
  max_luka_dni = max(as.integer(diff(index(ceny))))
)
print(podsumowanie)
write.csv(podsumowanie, file.path(DIR_TAB, "qc_podsumowanie.csv"), row.names = FALSE)

saveRDS(list(ceny = ceny, R = R), SCIEZKA_ZWR)
write.csv(data.frame(data = index(R), coredata(R)), "data/processed/zwroty.csv", row.names = FALSE)

# --- Wykres cen (indeks = 100 na starcie) ----------------------------------
idx <- sweep(coredata(ceny), 2, coredata(ceny)[1, ], "/") * 100
df <- data.frame(data = rep(index(ceny), ncol(idx)),
                 spolka = factor(rep(NAZWY[colnames(idx)], each = nrow(idx)), levels = NAZWY),
                 indeks = as.numeric(idx))
p <- ggplot(df, aes(data, indeks, colour = spolka)) +
  geom_line(linewidth = 0.5) +
  scale_y_log10() +
  scale_colour_manual(values = setNames(PALETA[1:5], NAZWY)) +
  labs(title = "Ceny skorygowane (dywidendy, splity), 2015-01-02 = 100",
       subtitle = "Skala logarytmiczna", x = NULL, y = "Indeks",
       caption = "Źródło: Yahoo Finance") +
  theme_projekt()
zapisz_wykres(p, "01_ceny")
message("OK: ", nrow(R), " dni zwrotów, ", start(R), " .. ", end(R))
