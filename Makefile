FC      = gfortran
BUILD   = build
SRC     = src
FFLAGS  = -O3 -fopenmp -std=f2018 -Wall -Wextra -J$(BUILD)

MODULES = kinds_mod config_mod map_mod energy_mod init_mod evolve_mod io_mod
OBJS    = $(addprefix $(BUILD)/, $(addsuffix .o, $(MODULES) main))
EXE     = $(BUILD)/kicked_top

.PHONY: all run clean figures docs

all: $(EXE)

$(EXE): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

$(BUILD)/%.o: $(SRC)/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD):
	mkdir -p $@

# module dependencies
$(BUILD)/config_mod.o: $(BUILD)/kinds_mod.o
$(BUILD)/map_mod.o:    $(BUILD)/kinds_mod.o
$(BUILD)/energy_mod.o: $(BUILD)/kinds_mod.o
$(BUILD)/init_mod.o:   $(BUILD)/kinds_mod.o $(BUILD)/config_mod.o
$(BUILD)/evolve_mod.o: $(BUILD)/kinds_mod.o $(BUILD)/config_mod.o $(BUILD)/map_mod.o
$(BUILD)/io_mod.o:     $(BUILD)/kinds_mod.o $(BUILD)/config_mod.o $(BUILD)/energy_mod.o
$(BUILD)/main.o:       $(addprefix $(BUILD)/, $(addsuffix .o, $(MODULES)))

run: $(EXE)
	./$(EXE) config.nml

# lecture notes
DOCS        = docs
FIG_CONFIGS = $(wildcard $(DOCS)/figures/configs/*.nml)

figures: $(EXE)
	for f in $(FIG_CONFIGS); do ./$(EXE) $$f || exit 1; done
	python3 $(DOCS)/figures/plot_figures.py

docs:
	cd $(DOCS) && latexmk -pdf -interaction=nonstopmode -halt-on-error kicked_top_notes.tex

clean:
	rm -rf $(BUILD)
	cd $(DOCS) && latexmk -c kicked_top_notes.tex
