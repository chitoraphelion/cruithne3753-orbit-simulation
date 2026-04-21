# 3753 Cruithne — численное моделирование

Воспроизведение результатов Wiegert, Innanen & Mikkola (1998, AJ 115, 2604)
для околоземного астероида 3753 Cruithne: подковообразная орбита в
коротирующем кадре, эволюция большой полуоси, устойчивость по ансамблям
клонов. Интегратор — REBOUND/WHFast; начальные состояния — JPL
Horizons/DE441 (TDB, ICRF, барицентрическая).

## Структура

```
src/cruithne/     — пакет (simulation, clones, analysis, kepler, frames, …)
scripts/          — пост-обработка и графики
report/           — Typst-отчёт (main.typ, sample.bib)
data/             — sim_main.npz, clones/*.csv, horizons_cache/
plots/generated/  — все рисунки (PNG)
```

## Требования

- Python ≥ 3.11
- [uv](https://docs.astral.sh/uv/) — менеджер окружения (`brew install uv`)
- [Typst](https://typst.app/) ≥ 0.13 (`brew install typst`)

## Установка

```sh
uv sync
uv pip install -e .
```

Первый запуск скачает эфемериды из JPL Horizons (закэшируются в
`data/horizons_cache/`).

## Быстрый старт (smoke-прогон, ~минута)

```sh
make smoke
```

Короткий прогон 1900–2100 гг., все основные графики и валидация против
Horizons. Клоны на коротком окне ±2000 лет — для проверки пайплайна.

## Полный прогон

```sh
make full
```

Выполняет:

1. `cruithne-main` — интегрирование 1900–3900 гг. (Δ𝑡 = 0.5 сут).
2. `cruithne-clones` — 64 клона × {10⁻⁵, 10⁻⁶} × ±13 kyr.
3. Все графики в `plots/generated/`.

На MacBook Pro M1 ≈ 40–60 минут.

## Отчёт

```sh
make report
```

Компилирует `report/main.pdf`.

## Отдельные шаги

```sh
uv run cruithne-main --fast                    # короткая интеграция
uv run cruithne-clones --years 13000 --sample-years 0.5   # клоны
uv run python scripts/plot_main_figures.py     # основные рисунки
uv run python scripts/plot_fig3_horseshoe.py   # подкова
uv run python scripts/plot_validation.py       # Horizons vs REBOUND
uv run python scripts/plot_clones.py           # гистограммы клонов
```

## Лицензия

Для учебных/исследовательских целей. Исходные эфемериды принадлежат
JPL/NASA, см. условия Horizons.
