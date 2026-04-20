# 3753 Cruithne — численное моделирование

Воспроизведение результатов Wiegert, Innanen & Mikkola (1998, AJ 115, 2604)
для околоземного астероида 3753 Cruithne: подковообразная орбита в
коротирующем кадре, эволюция большой полуоси, устойчивость по ансамблям клонов.

## Требования

- Python ≥ 3.11
- [uv](https://docs.astral.sh/uv/) для управления окружением
- [Typst](https://typst.app/) для сборки отчёта (`brew install typst`)

## Установка

```sh
uv sync
```

## Запуск

Быстрый smoke-прогон (короткий интервал, для проверки кода):

```sh
uv run cruithne-main --fast
```

Полный прогон (1900–3900 гг., ~минуты):

```sh
uv run cruithne-main
```

Ансамбль клонов (±13 kyr, 64 клона × 2 ε):

```sh
uv run cruithne-clones
```

Все графики сохраняются в `plots/generated/`.

## Сборка отчёта

```sh
cd report && typst compile main.typ
```
