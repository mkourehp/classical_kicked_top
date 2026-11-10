FC      = gfortran
BUILD   = build
SRC     = src
FFLAGS  = -O3 -fopenmp -std=f2018 -Wall -Wextra -J$(BUILD)
LAPACK  = -llapack -lblas

COMMON    = kinds_mod sphere_mod io_utils_mod
CLASSICAL = config_mod map_mod energy_mod init_mod evolve_mod lyapunov_mod io_mod main
QUANTUM   = quantum_config_mod spin_mod floquet_mod husimi_mod spectrum_mod quantum_io_mod quantum_main
CHAIN     = chain_state_mod chain_floquet_mod chain_observables_mod
COMBINED  = combined_config_mod combined_io_mod chain_run_mod combined_main
SPIN_CHAIN = spin_chain_main

objs = $(addprefix $(BUILD)/, $(addsuffix .o, $(1)))
COMMON_OBJS    = $(call objs, $(COMMON))
CLASSICAL_OBJS = $(call objs, $(CLASSICAL))
QUANTUM_OBJS   = $(call objs, $(QUANTUM))
CHAIN_OBJS     = $(call objs, $(CHAIN))
COMBINED_OBJS   = $(call objs, $(COMBINED))
SPIN_CHAIN_OBJS = $(call objs, $(SPIN_CHAIN))
# kicked_top_spin_chain shares everything of kicked_top_combined except its driver and the classical part
SPIN_CHAIN_DEPS = $(call objs, spin_mod husimi_mod combined_config_mod combined_io_mod chain_run_mod)
# kicked_top_combined reuses the classical modules and the spin-j Husimi code
COMBINED_DEPS   = $(call objs, config_mod map_mod energy_mod init_mod evolve_mod io_mod spin_mod husimi_mod)

EXE  = $(BUILD)/kicked_top
QEXE = $(BUILD)/kicked_top_quantum
CEXE = $(BUILD)/kicked_top_combined
SEXE = $(BUILD)/kicked_top_spin_chain

# make: builds all programs. make classical | quantum | spin_chain | combined: builds and runs one.
.PHONY: all classical quantum spin_chain combined clean doc

all: $(EXE) $(QEXE) $(SEXE) $(CEXE)

classical: $(EXE)
	./$(EXE) config_classical.nml

quantum: $(QEXE)
	./$(QEXE) config_quantum.nml

spin_chain: $(SEXE)
	./$(SEXE) config_spin_chain.nml

combined: $(CEXE)
	./$(CEXE) config_combined.nml

$(EXE): $(COMMON_OBJS) $(CLASSICAL_OBJS)
	$(FC) $(FFLAGS) -o $@ $^

$(QEXE): $(COMMON_OBJS) $(QUANTUM_OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LAPACK)

$(CEXE): $(COMMON_OBJS) $(COMBINED_DEPS) $(CHAIN_OBJS) $(COMBINED_OBJS)
	$(FC) $(FFLAGS) -o $@ $^

$(SEXE): $(COMMON_OBJS) $(SPIN_CHAIN_DEPS) $(CHAIN_OBJS) $(SPIN_CHAIN_OBJS)
	$(FC) $(FFLAGS) -o $@ $^

$(BUILD)/%.o: $(SRC)/common/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/%.o: $(SRC)/classical/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/%.o: $(SRC)/quantum/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/%.o: $(SRC)/chain/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/%.o: $(SRC)/combined/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/%.o: $(SRC)/spin_chain/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD):
	mkdir -p $@

# module dependencies
$(BUILD)/sphere_mod.o: $(BUILD)/kinds_mod.o
$(CLASSICAL_OBJS) $(QUANTUM_OBJS) $(CHAIN_OBJS) $(COMBINED_OBJS) $(SPIN_CHAIN_OBJS): $(COMMON_OBJS)

$(BUILD)/map_mod.o $(BUILD)/energy_mod.o: $(BUILD)/config_mod.o
$(BUILD)/init_mod.o:     $(BUILD)/config_mod.o
$(BUILD)/evolve_mod.o:   $(BUILD)/config_mod.o $(BUILD)/map_mod.o
$(BUILD)/lyapunov_mod.o: $(BUILD)/config_mod.o $(BUILD)/map_mod.o
$(BUILD)/io_mod.o:       $(BUILD)/config_mod.o $(BUILD)/energy_mod.o
$(BUILD)/main.o:         $(filter-out $(BUILD)/main.o, $(CLASSICAL_OBJS))

$(BUILD)/floquet_mod.o $(BUILD)/husimi_mod.o: $(BUILD)/spin_mod.o
$(BUILD)/quantum_io_mod.o: $(BUILD)/quantum_config_mod.o $(BUILD)/spectrum_mod.o
$(BUILD)/quantum_main.o:   $(filter-out $(BUILD)/quantum_main.o, $(QUANTUM_OBJS))

$(BUILD)/combined_io_mod.o: $(BUILD)/combined_config_mod.o
$(BUILD)/chain_run_mod.o:  $(BUILD)/husimi_mod.o $(BUILD)/combined_io_mod.o $(CHAIN_OBJS)
$(BUILD)/combined_main.o:   $(COMBINED_DEPS) $(CHAIN_OBJS) $(filter-out $(BUILD)/combined_main.o, $(COMBINED_OBJS))
$(BUILD)/spin_chain_main.o: $(SPIN_CHAIN_DEPS) $(CHAIN_OBJS)

# reports: make doc <dir> [compute], <dir> a directory in docs/.
# Without compute, compiles every .tex in docs/<dir>, assuming its figures exist.
# With compute, first runs every namelist in docs/<dir>/figures/configs{,_quantum,_spin_chain,_combined}/
# with kicked_top, kicked_top_quantum, kicked_top_spin_chain and kicked_top_combined,
# then every docs/<dir>/figures/plot_*.py.
DOCS = docs
ifeq (doc,$(firstword $(MAKECMDGOALS)))
  REPORT  := $(word 2,$(MAKECMDGOALS))
  COMPUTE := $(filter compute,$(wordlist 3,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS)))
  ifneq ($(REPORT),)
    $(eval $(REPORT) compute: ; @:)
  endif
endif
REPORT_DIR = $(DOCS)/$(REPORT)
REPORT_FIG = $(REPORT_DIR)/figures

doc: $(if $(COMPUTE),$(EXE) $(QEXE) $(SEXE) $(CEXE))
	@test -n "$(REPORT)" && test -d "$(REPORT_DIR)" || { echo "usage: make doc <dir> [compute], <dir> a directory in $(DOCS)/"; exit 1; }
ifneq ($(COMPUTE),)
	for f in $(wildcard $(REPORT_FIG)/configs/*.nml); do ./$(EXE) $$f || exit 1; done
	for f in $(wildcard $(REPORT_FIG)/configs_quantum/*.nml); do ./$(QEXE) $$f || exit 1; done
	for f in $(wildcard $(REPORT_FIG)/configs_spin_chain/*.nml); do ./$(SEXE) $$f || exit 1; done
	for f in $(wildcard $(REPORT_FIG)/configs_combined/*.nml); do ./$(CEXE) $$f || exit 1; done
	for f in $(wildcard $(REPORT_FIG)/plot_*.py); do python3 $$f || exit 1; done
endif
	cd $(REPORT_DIR) && latexmk -pdf -interaction=nonstopmode -halt-on-error $(notdir $(wildcard $(REPORT_DIR)/*.tex))

clean:
	rm -rf $(BUILD)
	for d in $(wildcard $(DOCS)/*/); do (cd $$d && latexmk -c) || exit 1; done
