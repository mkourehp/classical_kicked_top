FC      = gfortran
BUILD   = build
SRC     = src
FFLAGS  = -O3 -fopenmp -std=f2018 -Wall -Wextra -J$(BUILD)
LAPACK  = -llapack -lblas

COMMON    = kinds_mod sphere_mod io_utils_mod
CLASSICAL = config_mod map_mod energy_mod init_mod evolve_mod lyapunov_mod io_mod main
QUANTUM   = quantum_config_mod spin_mod floquet_mod husimi_mod spectrum_mod quantum_io_mod quantum_main
CHAIN     = chain_state_mod chain_floquet_mod chain_observables_mod
REGIMES   = regimes_config_mod regimes_io_mod regimes_main

objs = $(addprefix $(BUILD)/, $(addsuffix .o, $(1)))
COMMON_OBJS    = $(call objs, $(COMMON))
CLASSICAL_OBJS = $(call objs, $(CLASSICAL))
QUANTUM_OBJS   = $(call objs, $(QUANTUM))
CHAIN_OBJS     = $(call objs, $(CHAIN))
REGIMES_OBJS   = $(call objs, $(REGIMES))
# kicked_top_regimes reuses the classical modules and the spin-j Husimi code
REGIMES_DEPS   = $(call objs, config_mod map_mod energy_mod init_mod evolve_mod io_mod spin_mod husimi_mod)

EXE  = $(BUILD)/kicked_top
QEXE = $(BUILD)/kicked_top_quantum
REXE = $(BUILD)/kicked_top_regimes

.PHONY: all classical quantum regimes run run-quantum run-regimes clean figures docs

all: classical quantum regimes

classical: $(EXE)

quantum: $(QEXE)

regimes: $(REXE)

$(EXE): $(COMMON_OBJS) $(CLASSICAL_OBJS)
	$(FC) $(FFLAGS) -o $@ $^

$(QEXE): $(COMMON_OBJS) $(QUANTUM_OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LAPACK)

$(REXE): $(COMMON_OBJS) $(REGIMES_DEPS) $(CHAIN_OBJS) $(REGIMES_OBJS)
	$(FC) $(FFLAGS) -o $@ $^

$(BUILD)/%.o: $(SRC)/common/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/%.o: $(SRC)/classical/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/%.o: $(SRC)/quantum/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/%.o: $(SRC)/chain/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/%.o: $(SRC)/regimes/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD):
	mkdir -p $@

# module dependencies
$(BUILD)/sphere_mod.o: $(BUILD)/kinds_mod.o
$(CLASSICAL_OBJS) $(QUANTUM_OBJS) $(CHAIN_OBJS) $(REGIMES_OBJS): $(COMMON_OBJS)

$(BUILD)/map_mod.o $(BUILD)/energy_mod.o: $(BUILD)/config_mod.o
$(BUILD)/init_mod.o:     $(BUILD)/config_mod.o
$(BUILD)/evolve_mod.o:   $(BUILD)/config_mod.o $(BUILD)/map_mod.o
$(BUILD)/lyapunov_mod.o: $(BUILD)/config_mod.o $(BUILD)/map_mod.o
$(BUILD)/io_mod.o:       $(BUILD)/config_mod.o $(BUILD)/energy_mod.o
$(BUILD)/main.o:         $(filter-out $(BUILD)/main.o, $(CLASSICAL_OBJS))

$(BUILD)/floquet_mod.o $(BUILD)/husimi_mod.o: $(BUILD)/spin_mod.o
$(BUILD)/quantum_io_mod.o: $(BUILD)/quantum_config_mod.o $(BUILD)/spectrum_mod.o
$(BUILD)/quantum_main.o:   $(filter-out $(BUILD)/quantum_main.o, $(QUANTUM_OBJS))

$(BUILD)/regimes_io_mod.o: $(BUILD)/regimes_config_mod.o
$(BUILD)/regimes_main.o:   $(REGIMES_DEPS) $(CHAIN_OBJS) $(filter-out $(BUILD)/regimes_main.o, $(REGIMES_OBJS))

run: $(EXE)
	./$(EXE) config.nml

run-quantum: $(QEXE)
	./$(QEXE) config_quantum.nml

run-regimes: $(REXE)
	./$(REXE) config_regimes.nml

# lecture notes
DOCS        = docs
FIG_CONFIGS  = $(wildcard $(DOCS)/figures/configs/*.nml)
QFIG_CONFIGS = $(wildcard $(DOCS)/figures/configs_quantum/*.nml)
RFIG_CONFIGS = $(wildcard $(DOCS)/figures/configs_regimes/*.nml)

figures: $(EXE) $(QEXE) $(REXE)
	for f in $(FIG_CONFIGS); do ./$(EXE) $$f || exit 1; done
	for f in $(QFIG_CONFIGS); do ./$(QEXE) $$f || exit 1; done
	for f in $(RFIG_CONFIGS); do ./$(REXE) $$f || exit 1; done
	python3 $(DOCS)/figures/plot_figures.py
	python3 $(DOCS)/figures/plot_three_regimes.py

docs:
	cd $(DOCS) && latexmk -pdf -interaction=nonstopmode -halt-on-error kicked_top_notes.tex three_regimes.tex

clean:
	rm -rf $(BUILD)
	cd $(DOCS) && latexmk -c kicked_top_notes.tex three_regimes.tex
