.PHONY: install smoke full figures clones report clean

install:
	uv sync
	uv pip install -e .

smoke:
	uv run cruithne-main --fast
	uv run cruithne-clones --fast
	$(MAKE) figures

full:
	uv run cruithne-main
	uv run cruithne-clones --years 13000 --sample-years 0.5
	$(MAKE) figures

figures:
	uv run python scripts/plot_main_figures.py
	uv run python scripts/plot_fig3_horseshoe.py
	uv run python scripts/plot_validation.py
	uv run python scripts/plot_clones.py

clones:
	uv run cruithne-clones --years 13000 --sample-years 0.5
	uv run python scripts/plot_clones.py

report:
	cd report && typst compile --root .. main.typ

clean:
	rm -rf plots/generated/*.png report/main.pdf
