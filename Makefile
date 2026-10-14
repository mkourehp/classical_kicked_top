FC      = gfortran
BUILD   = build
SRC     = src
FFLAGS  = -O3 -fopenmp -std=f2018 -Wall -Wextra -J$(BUILD)

MODULES = kinds_mod config_mod map_mod init_mod evolve_mod io_mod
OBJS    = $(addprefix $(BUILD)/, $(addsuffix .o, $(MODULES) main))
EXE     = $(BUILD)/kicked_top

.PHONY: all run clean

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
$(BUILD)/init_mod.o:   $(BUILD)/kinds_mod.o $(BUILD)/config_mod.o
$(BUILD)/evolve_mod.o: $(BUILD)/kinds_mod.o $(BUILD)/config_mod.o $(BUILD)/map_mod.o
$(BUILD)/io_mod.o:     $(BUILD)/kinds_mod.o $(BUILD)/config_mod.o
$(BUILD)/main.o:       $(addprefix $(BUILD)/, $(addsuffix .o, $(MODULES)))

run: $(EXE)
	./$(EXE) config.nml

clean:
	rm -rf $(BUILD)
