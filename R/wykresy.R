# Wspólny wygląd wykresów: paleta, motyw ggplot2, zapis do PNG.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

# Paleta kategorii w stałej kolejności (sprawdzona pod kątem daltonizmu).
PALETA <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
            "#008300", "#4a3aa7", "#e34948")
SZARY  <- "#8a8984"

# Stałe kolory metod, żeby ta sama metoda miała ten sam kolor na każdym wykresie.
KOLORY_METOD <- c(
  "Historyczna"    = PALETA[1],
  "Normalna"       = PALETA[2],
  "MC normalny"    = PALETA[3],
  "MC t-Student"   = PALETA[4],
  "MC bootstrap"   = PALETA[5],
  "GARCH-norm"     = PALETA[7],
  "GARCH-t"        = PALETA[8]
)

# theme_projekt: minimalistyczny motyw z dyskretną siatką.
theme_projekt <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3),
      plot.title       = element_text(face = "bold", colour = "#0b0b0b"),
      plot.subtitle    = element_text(colour = "#52514e"),
      plot.caption     = element_text(colour = "#52514e", size = rel(0.8)),
      axis.text        = element_text(colour = "#52514e"),
      legend.position  = "top",
      legend.title     = element_blank(),
      strip.text       = element_text(face = "bold", hjust = 0),
      plot.background  = element_rect(fill = "white", colour = NA)
    )
}

# zapisz_wykres: zapis wykresu do output/figures/<nazwa>.png
zapisz_wykres <- function(p, nazwa, w = 8, h = 5) {
  ggsave(file.path(DIR_FIG, paste0(nazwa, ".png")), p,
         width = w, height = h, dpi = 150, bg = "white")
  invisible(p)
}

# fmt_pln: format kwoty w PLN do tabel i podpisów.
fmt_pln <- function(x) format(round(x), big.mark = " ", scientific = FALSE, trim = TRUE)
