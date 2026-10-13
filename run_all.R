# Uruchamia cały projekt od początku: skrypty 01-07, każdy w osobnej czystej sesji Rscript,
# a potem renderuje raport.qmd do output/raport.html.
# Użycie (w katalogu projektu): Rscript run_all.R

skrypty <- sprintf("%02d_%s.R", 1:7, c("dane", "eksploracja", "var_klasyczny", "monte_carlo",
                                        "porownanie", "backtest", "garch"))
rscript <- file.path(R.home("bin"), "Rscript")
for (s in skrypty) {
  message("\n==== ", s, " ====")
  t0 <- Sys.time()
  kod <- system2(rscript, c("--vanilla", s))
  if (kod != 0) stop("Błąd w skrypcie ", s)
  message("---- ", s, " OK (", round(difftime(Sys.time(), t0, units = "secs")), " s)")
}

# --- Raport ----------------------------------------------------------------
# Quarto: najpierw z PATH, potem kopia dołączona do RStudio (na Macu nie ma jej w PATH).
quarto <- Sys.which("quarto")
if (!nzchar(quarto)) {
  kandydat <- "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto"
  if (file.exists(kandydat)) quarto <- kandydat
}
if (nzchar(quarto)) {
  message("\n==== raport.qmd (Quarto: ", quarto, ") ====")
  kod <- system2(quarto, c("render", "raport.qmd", "--to", "html"))
  if (kod != 0) stop("Błąd renderowania raportu")
  # Quarto przy pojedynczym pliku zapisuje HTML obok .qmd; przenosimy do output/.
  file.rename("raport.html", "output/raport.html")
} else {
  # Plan awaryjny bez Quarto: rmarkdown czyta ten sam plik (nagłówek format: html jest pomijany,
  # więc wynik wygląda prościej). Instalacja Quarto: patrz README.md.
  message("\nNie znaleziono Quarto. Renderuję przez rmarkdown::render (wersja uproszczona).")
  rmarkdown::render("raport.qmd", output_format = rmarkdown::html_document(toc = TRUE, self_contained = TRUE),
                    output_file = "raport.html", output_dir = "output")
}
message("\nGotowe: output/raport.html")
