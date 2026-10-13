# Ryzyko portfela: VaR i CVaR metodą Monte Carlo (R)

> **Projekt edukacyjny. Nie jest rekomendacją inwestycyjną.**

**Raport online:** https://alkslayer.github.io/var-cvar-monte-carlo/

Portfel 100 000 PLN po równo w 5 spółkach z GPW (PKO BP, PZU, CD Projekt, KGHM, Orlen). Liczę VaR i CVaR metodą historyczną, parametryczną (normalną) i Monte Carlo (normalny, t-Student, bootstrap), backtestuję je w oknie kroczącym i porównuję z GARCH. Dane dzienne 2015-01-05 – 2026-10-02.

## Jak uruchomić (Mac)

Wymagane: R ≥ 4.3 i pakiety:

```r
install.packages(c("quantmod", "xts", "zoo", "ggplot2", "scales", "patchwork",
                   "MASS", "rugarch", "knitr", "rmarkdown"))
```

Całość (ok. 1–2 min, każdy skrypt w czystej sesji `Rscript --vanilla`, na końcu render raportu):

```bash
cd ~/Desktop/"Zadyma ekonometryczna/Projekt 3"
Rscript run_all.R
```

Raport otworzysz poleceniem `open output/raport.html` (jeden plik HTML z osadzonymi wykresami).

Pojedynczy etap: `Rscript 04_monte_carlo.R`. Skrypty trzeba uruchamiać po kolei, bo każdy czyta wyniki poprzednich z `data/processed/` i `output/`. Dane są pobierane z Yahoo tylko raz (`data/raw/yahoo_ceny.rds`). Odświeżenie danych: `REFRESH=1 Rscript 01_dane.R`. Uwaga: po odświeżeniu zmienią się wyniki, a `DATA_DO` w `R/config.R` trzeba ustawić ręcznie.

### Quarto

Quarto nie jest w `PATH`, ale jest dołączone do RStudio (`/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto`, wersja 1.x). `run_all.R` znajduje je samo. Gdyby go nie było, skrypt renderuje raport przez `rmarkdown` w uproszczonym wyglądzie. Żeby mieć `quarto` w terminalu:

```bash
brew install --cask quarto
```

Druga opcja to dopisać do `~/.zshrc`: `export PATH="/Applications/RStudio.app/Contents/Resources/app/quarto/bin:$PATH"`.

## Struktura

```
Projekt 3/
  R/config.R          parametry: tickery, wagi, daty, poziomy ufności, N symulacji, okno backtestu
  R/dane.R            pobieranie z Yahoo, kontrola jakości, czyszczenie, zwroty log
  R/ryzyko.R          VaR/CVaR z próby, normalny, t-Student; statystyki; kontrola spójności
  R/mc.R              symulacje MC (Cholesky, t wielowymiarowy, bootstrap), MLE stopni swobody
  R/backtest.R        prognozy kroczące, test Kupca, test Christoffersena
  R/wykresy.R         paleta i motyw ggplot2
  01_dane.R           pobranie, kontrola jakości, zwroty           -> data/, output/tables/qc_*
  02_eksploracja.R    rozkłady, QQ, kurtoza, ACF, korelacje         -> output/figures/02_*
  03_var_klasyczny.R  VaR/CVaR historyczny i normalny
  04_monte_carlo.R    MC: normalny, t, bootstrap; zbieżność; kontrola z wzorami
  05_porownanie.R     tabela zbiorcza, wrażliwość na alfa i horyzont
  06_backtest.R       backtest 500 dni, Kupiec, Christoffersen, marzec 2020, test look-ahead
  07_garch.R          GARCH(1,1) norm/t w backteście, warianty odporności, prognoza bieżąca
  run_all.R           uruchamia 01-07 i renderuje raport
  raport.qmd          raport (Quarto) -> output/raport.html
  docs/index.html     kopia raportu dla GitHub Pages
  data/raw/           surowe notowania z Yahoo (rds, z datą pobrania w atrybucie)
  data/processed/     zwroty log (rds + csv)
  output/figures/     wykresy PNG
  output/tables/      tabele CSV
```

## Podjęte decyzje (co i dlaczego)

| Decyzja | Wybór | Dlaczego |
|---|---|---|
| Aktywa | PKO.WA, PZU.WA, CDR.WA, KGH.WA, PKN.WA | Zweryfikowane w Yahoo: wszystkie mają dane od 2015-01-01, ten sam kalendarz i brak dziur poza wspólnymi (opis niżej). Plan B (ETF-y USD) nie był potrzebny. Odrzucone alternatywy: DNP.WA (dane od 2017), ALE.WA (od 2020), SPL.WA (Yahoo zwraca błąd). |
| Koniec próby | 2026-10-02 | Od 5 do 9.10.2026 Yahoo podaje płaskie ceny z zerowym wolumenem dla wszystkich spółek, a 12.10 to niepełna sesja z brakiem Close. |
| Ceny | Adjusted (dywidendy + splity) | Wymóg zadania. Bez korekty dzień odcięcia dywidendy wyglądałby jak strata. |
| Zwrot portfela | r_p = w'r (średnia ważona zwrotów log), wagi stałe 20% | Prosta agregacja w czasie (sumy) i spójność wszystkich metod. Różnica wobec dokładnego zwrotu portfela: VaR 99% 1d wyższy o 39 PLN (ok. 1%), czyli przybliżenie jest lekko konserwatywne. |
| Strata | L = V(1 − e^r), VaR i CVaR jako dodatnie kwoty w PLN | Dokładna strata kwotowa dla zwrotu log. |
| Średnia | Uwzględniona (średnia z próby) | Spójnie we wszystkich metodach. Przy 1 dniu ma znikomy wpływ, przy 10 dniach obniża VaR o ok. 0,6% wartości portfela. |
| Historyczny 10 dni | Zachodzące na siebie 10-dniowe okna | ~2900 obserwacji zamiast 293 rozłącznych. Dla porównania liczę też regułę √10 (`output/tables/03_var_klasyczny.csv`). |
| CVaR normalny | Wzór dokładny dla straty lognormalnej | Pozwala sprawdzić MC normalny co do złotówki. |
| MC t-Student | Wielowymiarowy t ze wspólnym ν, kowariancja = Σ z próby, ν z MLE (profil wiarygodności) | Wspólny mianownik daje zależność w ogonach (krach dotyka wszystkich naraz). ν = 5,97. Portfel ma wtedy dokładnie rozkład t, co daje analityczną kontrolę MC. |
| Bootstrap | iid po całych dniach | Zachowuje korelacje i empiryczne ogony. Nie zachowuje skupisk zmienności (to jest testowane w GARCH). |
| N symulacji | 100 000, `set.seed(20261013)` | Z wykresu zbieżności: względny błąd standardowy VaR 99% to 0,5–0,8%, CVaR 99% ok. 0,5–1,2% (zależnie od metody). |
| Backtest | VaR 1-dniowy 95% i 99%, okno 500 dni, prognoza na dzień t tylko z dni t−500…t−1 | Wymóg zadania. Look-ahead sprawdzany automatycznie (zepsucie przyszłych danych nie zmienia prognoz; okna GARCH kończą się przed dniem prognozy). |
| MC-t w backteście | ν reestymowane co 20 dni, 20 000 symulacji dziennie | Kompromis czasu i dokładności (backtest liczy się w ok. 7 s). |
| GARCH | sGARCH(1,1), stała w średniej, innowacje normalne i t, reestymacja co 20 dni w oknie 500 (`ugarchroll`) | Standardowy punkt odniesienia, ustalony przed obejrzeniem wyników. Estymacja na zwrotach ×100, bo na surowych optymalizator nie zbiegał się w 10 oknach. Warianty GJR i bez stałej pokazane jako test odporności. |
| Testy | Kupiec (POF), Christoffersen (niezależność i łączny) | Własna implementacja, sprawdzona z `rugarch::VaRTest` (te same statystyki LR). |
| Stres | Marzec 2020 (22 sesje) | Wymóg zadania. Dodatkowo podział przekroczeń na lata. |

## Dane: co znalazłem i jak to obsłużyłem

- **Zakres:** Yahoo zwraca 2995 wierszy na spółkę (2015-01-01 – 2026-10-12). Po czyszczeniu zostaje 2939 dni cen i 2938 zwrotów (2015-01-05 – 2026-10-02).
- **Zerowy wolumen:** 55 dni wspólnych dla wszystkich 5 spółek, z ceną przepisaną z poprzedniego dnia (sprawdzone w kodzie: zmiana ceny = 0). Do 2020 r. to głównie święta GPW, które Yahoo wstawiało jako sztuczne wiersze, plus pojedyncze dziury (np. 30.01.2018) i tydzień 5–9.10.2026. Usunąłem je. Gdyby zostały, zaniżyłyby zmienność sztucznymi zerami, a w pojedynczych dziurach zwrot z dwóch dni trafia do jednego dnia.
- **Braki:** 1 brak (12.10.2026, niepełna sesja), poza próbą.
- **Dni bez notowań:** najdłuższa przerwa to 6 dni kalendarzowych. Wszystkie przerwy powyżej 4 dni to Wielkanoc i Boże Narodzenie (`output/tables/qc_luki.csv`).
- **Splity:** Yahoo zgłasza tylko split PZU 1:10 (30.11.2015). Szereg Close jest już skorygowany (brak skoku).
- **Dywidendy:** 38 dni, w których zwrot z Adjusted różni się od zwrotu z Close o więcej niż 1 pp (PZU i Orlen co roku, PKO BP, KGHM). To korekty dywidendowe, poprawnie ujęte w Adjusted.
- **Obserwacje odstające (> 6 sd):** 9 dni, wszystkie to prawdziwe zdarzenia: 12.03.2020 (COVID, PKO/PZU/KGHM), 24–25.02.2022 (inwazja na Ukrainę, PKO −18%), 09.11.2020 (Orlen +13%), CD Projekt 09.06.2015 (+19%) i 29.01.2021 (−18%), PKO 22.08.2025 (−13%). Zostawiłem je, bo to właśnie ryzyko, które mierzymy.
- **Weryfikacja zewnętrzna:** próba porównania ze stooq.pl nie powiodła się (strona blokuje automatyczne pobieranie). Dane pochodzą wyłącznie z Yahoo.

## Najważniejsze wyniki (portfel 100 000 PLN)

| Metoda | VaR 99% 1d | CVaR 99% 1d | VaR 99% 10d |
|---|---|---|---|
| Historyczna | 3 652 | 5 049 | 12 077 |
| Normalna | 3 395 | 3 887 | 9 996 |
| MC normalny | 3 416 | 3 895 | 10 068 |
| MC t-Student | 3 741 | 4 797 | 10 211 |
| MC bootstrap | 3 611 | 5 069 | 10 496 |

Backtest VaR 99% (2438 dni, oczekiwane 24,4 przekroczenia): normalny 42 (Kupiec odrzuca), historyczny 32, MC-t 32 (Kupiec OK, Christoffersen odrzuca: przekroczenia w skupiskach), GARCH-t 37, GJR-GARCH-t bez stałej 32 (przechodzi oba testy). Marzec 2020, VaR 99%: 4–5 przekroczeń w modelach statycznych, 3 w GARCH (oczekiwane 0,22).

## Kontrole jakości wbudowane w kod

- `sprawdz_spojnosc()`: zatrzymuje skrypt, jeśli VaR99 < VaR95 albo CVaR < VaR (03, 04, 05, 07).
- MC normalny (10⁶) vs wzór parametryczny i MC-t vs dokładny rozkład t: błąd < 1% VaR i < 1,5% CVaR, inaczej `stop` (04).
- Look-ahead: po zepsuciu danych od dnia 1500 prognozy na wcześniejsze dni są identyczne (06), a każde okno estymacji GARCH kończy się przed dniem prognozy (07).
- Liczby przekroczeń i statystyki LR zgadzają się z `rugarch::VaRTest` (06, 07).
- Zgodność dat i realizacji backtestu GARCH z 06.

## Znane ograniczenia

- Metody historyczna, normalna i MC zakładają niezależne dni o stałym rozkładzie. Kwadraty zwrotów są silnie autokorelowane, więc VaR 10-dniowy z tych metod jest zaniżony wobec historycznych okien 10-dniowych.
- Korelacja jest stała, a w okresach wysokiej zmienności rośnie (0,30 → 0,41). Nie ma modelu DCC ani kopuł.
- Wspólne ν dla wszystkich spółek, choć CD Projekt ma grubszy ogon (ν ≈ 3,5).
- Okna 10-dniowe się nakładają, więc efektywna próba jest ok. 10 razy mniejsza.
- Backtest dotyczy tylko VaR 1-dniowego. Nie ma backtestu CVaR (np. test Acerbi-Szekely).
- GJR-GARCH (najlepszy wariant) wybrany po obejrzeniu wyników, więc to wskazówka, a nie dowód.
- Dane tylko z Yahoo, bez kosztów transakcyjnych, płynności i podatków.
- Wyniki zależą od dnia pobrania danych (Yahoo może retroaktywnie zmienić Adjusted).
