# Ustawienia całego projektu: portfel, okres, parametry metod.
# Każdy skrypt 0X_*.R zaczyna od source("R/config.R").

# --- Portfel ---------------------------------------------------------------
TICKERY <- c("PKO.WA", "PZU.WA", "CDR.WA", "KGH.WA", "PKN.WA")
NAZWY   <- c(PKO.WA = "PKO BP", PZU.WA = "PZU", CDR.WA = "CD Projekt",
             KGH.WA = "KGHM", PKN.WA = "Orlen")
WAGI    <- setNames(rep(1 / length(TICKERY), length(TICKERY)), TICKERY)
WARTOSC <- 100000  # wartość portfela w PLN

# --- Okres -----------------------------------------------------------------
DATA_OD <- as.Date("2015-01-01")
# Ostatni wiarygodny dzień: od 2026-10-05 do 2026-10-09 Yahoo podaje zerowy
# wolumen i płaskie ceny (dziura w danych), a 2026-10-12 to niepełna sesja (NA).
DATA_DO <- as.Date("2026-10-02")

# --- Metody ----------------------------------------------------------------
ALFY      <- c(0.95, 0.99)   # poziomy ufności
HORYZONTY <- c(1, 10)        # horyzonty w dniach sesyjnych
N_SIM     <- 100000          # liczba symulacji MC w wynikach głównych
SEED      <- 20261013
OKNO_BT   <- 500             # długość okna kroczącego w backteście
REFIT_CO  <- 20              # co ile dni ponowna estymacja nu (MC-t) i GARCH
N_SIM_BT  <- 20000           # liczba symulacji MC-t na każdy dzień backtestu
STRES_OD  <- as.Date("2020-03-01")
STRES_DO  <- as.Date("2020-03-31")

# --- Ścieżki ---------------------------------------------------------------
SCIEZKA_RAW  <- "data/raw/yahoo_ceny.rds"
SCIEZKA_ZWR  <- "data/processed/zwroty.rds"
DIR_FIG      <- "output/figures"
DIR_TAB      <- "output/tables"

# Funkcje projektu
for (f in c("R/wykresy.R", "R/dane.R", "R/ryzyko.R", "R/mc.R", "R/backtest.R")) source(f)
