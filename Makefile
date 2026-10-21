FC      = gfortran
BUILD   = build
SRC     = src
FFLAGS  = -O3 -fopenmp -std=f2018 -Wall -Wextra -J$(BUILD)
LAPACK  = -llapack -lblas

COMMON    = kinds_mod sphere_mod io_utils_mod
CLASSICAL = config_mod map_mod energy_mod init_mod evolve_mod lyapunov_mod io_mod main
QUANTUM   = quantum_config_mod spin_mod floquet_mod husimi_mod spectrum_mod quantum_io_mod quantum_main

objs = $(addprefix $(BUILD)/, $(addsuffix .o, $(1)))
COMMON_OBJS    = $(call objs, $(COMMON))
CLASSICAL_OBJS = $(call objs, $(CLASSICAL))
QUANTUM_OBJS   = $(call objs, $(QUANTUM))

EXE  = $(BUILD)/kicked_top
QEXE = $(BUILD)/kicked_top_quantum

.PHONY: all classical quantum run run-quantum clean figures docs

all: classical quantum

classical: $(EXE)

quantum: $(QEXE)

$(EXE): $(COMMON_OBJS) $(CLASSICAL_OBJS)
	$(FC) $(FFLAGS) -o $@ $^

$(QEXE): $(COMMON_OBJS) $(QUANTUM_OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LAPACK)

$(BUILD)/%.o: $(SRC)/common/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/%.o: $(SRC)/classical/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/%.o: $(SRC)/quantum/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD):
	mkdir -p $@

# module dependencies
$(BUILD)/sphere_mod.o: $(BUILD)/kinds_mod.o
$(CLASSICAL_OBJS) $(QUANTUM_OBJS): $(COMMON_OBJS)

$(BUILD)/map_mod.o $(BUILD)/energy_mod.o: $(BUILD)/config_mod.o
$(BUILD)/init_mod.o:     $(BUILD)/config_mod.o
$(BUILD)/evolve_mod.o:   $(BUILD)/config_mod.o $(BUILD)/map_mod.o
$(BUILD)/lyapunov_mod.o: $(BUILD)/config_mod.o $(BUILD)/map_mod.o
$(BUILD)/io_mod.o:       $(BUILD)/config_mod.o $(BUILD)/energy_mod.o
$(BUILD)/main.o:         $(filter-out $(BUILD)/main.o, $(CLASSICAL_OBJS))

$(BUILD)/floquet_mod.o $(BUILD)/husimi_mod.o: $(BUILD)/spin_mod.o
$(BUILD)/quantum_io_mod.o: $(BUILD)/quantum_config_mod.o $(BUILD)/spectrum_mod.o
$(BUILD)/quantum_main.o:   $(filter-out $(BUILD)/quantum_main.o, $(QUANTUM_OBJS))

run: $(EXE)
	./$(EXE) config.nml

run-quantum: $(QEXE)
	./$(QEXE) config_quantum.nml

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
