# Pobieranie notowań z Yahoo Finance i kontrola jakości danych.

suppressPackageStartupMessages({
  library(quantmod)
  library(xts)
})

# pobierz_yahoo: pobiera OHLCV + Adjusted dla tickerów, zwraca listę obiektów xts.
pobierz_yahoo <- function(tickery, od) {
  out <- lapply(tickery, function(t) {
    x <- getSymbols(t, src = "yahoo", from = od, auto.assign = FALSE)
    colnames(x) <- c("Open", "High", "Low", "Close", "Volume", "Adjusted")
    x
  })
  setNames(out, tickery)
}

# kontrola_jakosci: tabela zakresu dat, braków i dni z zerowym wolumenem dla każdego tickera.
kontrola_jakosci <- function(lista) {
  do.call(rbind, lapply(names(lista), function(t) {
    x <- lista[[t]]
    data.frame(
      ticker        = t,
      pierwszy_dzien = as.character(start(x)),
      ostatni_dzien  = as.character(end(x)),
      wiersze       = nrow(x),
      braki_adjusted = sum(is.na(x$Adjusted)),
      wolumen_zero  = sum(x$Volume == 0, na.rm = TRUE),
      ujemne_lub_zero_ceny = sum(x$Adjusted <= 0, na.rm = TRUE)
    )
  }))
}

# dni_wolumen_zero: daty, w których WSZYSTKIE tickery mają zerowy wolumen
# (sztuczne wiersze Yahoo w święta / dziury w danych, cena przepisana z dnia poprzedniego).
dni_wolumen_zero <- function(lista) {
  v <- do.call(merge, lapply(lista, function(x) x$Volume))
  index(v)[rowSums(v == 0, na.rm = TRUE) == ncol(v)]
}

# czyste_ceny: macierz cen skorygowanych (xts) po usunięciu braków,
# dni z zerowym wolumenem i obserwacji spoza [od, do].
czyste_ceny <- function(lista, od, do) {
  ceny <- do.call(merge, lapply(lista, function(x) x$Adjusted))
  colnames(ceny) <- names(lista)
  zera <- dni_wolumen_zero(lista)
  ceny <- ceny[!(index(ceny) %in% zera)]
  ceny <- ceny[paste0(od, "/", do)]
  ceny[complete.cases(ceny)]
}

# zwroty_log: dzienne zwroty logarytmiczne (pierwszy wiersz NA usunięty).
zwroty_log <- function(ceny) {
  r <- diff(log(ceny))
  r[-1]
}

# dywidendy_i_splity: dni, w których zwrot z Adjusted różni się od zwrotu z Close
# o więcej niż prog (korekta dywidendowa) oraz duże skoki Close (możliwy split).
dywidendy_i_splity <- function(lista, od, do, prog = 0.01) {
  do.call(rbind, lapply(names(lista), function(t) {
    x  <- lista[[t]][paste0(od, "/", do)]
    x  <- x[x$Volume > 0 & complete.cases(x)]
    rc <- diff(log(x$Close)); ra <- diff(log(x$Adjusted))
    d  <- as.numeric(ra - rc)
    i  <- which(abs(d) > prog)
    if (!length(i)) return(NULL)
    data.frame(ticker = t, data = as.character(index(rc)[i]),
               zwrot_close = round(as.numeric(rc[i]), 4),
               zwrot_adjusted = round(as.numeric(ra[i]), 4),
               roznica = round(d[i], 4))
  }))
}

# odstajace: zwroty większe co do modułu niż k odchyleń standardowych danego aktywa.
odstajace <- function(r, k = 6) {
  s <- apply(r, 2, sd)
  do.call(rbind, lapply(colnames(r), function(t) {
    i <- which(abs(r[, t]) > k * s[t])
    if (!length(i)) return(NULL)
    data.frame(ticker = t, data = as.character(index(r)[i]),
               zwrot = round(as.numeric(r[i, t]), 4),
               ile_sd = round(as.numeric(r[i, t]) / s[t], 1))
  }))
}

# luki_sesyjne: przerwy między kolejnymi notowaniami dłuższe niż max_dni dni kalendarzowych.
luki_sesyjne <- function(daty, max_dni = 4) {
  d <- diff(daty)
  i <- which(d > max_dni)
  data.frame(od = as.character(daty[i]), do = as.character(daty[i + 1]),
             dni_kalendarzowych = as.integer(d[i]))
}
