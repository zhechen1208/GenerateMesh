# GenerateMesh Makefile
#   make          — release: mesh + beta
#   make mesh     — release: mesh only
#   make beta     — release: beta_function_wing only
#   make debug    — debug build (checks + backtrace)
#   make clean    — remove build artifacts

FC        = mpif90
GFC       = gfortran
FFLAGS   ?= -O2 -ffree-line-length-none
FFLAGS_D  = -g -O0 -fcheck=all -fbacktrace -ffree-line-length-none
BUILD     = build
SRC       = src
INC       = params.inc

MOD_OBJS  = $(BUILD)/shared_date_module.o $(BUILD)/tip_square_o_module.o
OBJS      = $(MOD_OBJS) $(BUILD)/main.o
TARGET    = $(BUILD)/mesh
BETA      = $(BUILD)/beta_function_wing

# default: both
all: $(TARGET) $(BETA)

# mesh
$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

mesh: $(TARGET)

debug: FFLAGS = $(FFLAGS_D)
debug: $(TARGET)

# beta_function_wing (not MPI, always debug flags)
$(BETA): $(SRC)/beta_function_wing.f90 | $(BUILD)
	$(GFC) $(FFLAGS_D) $< -o $@

beta: $(BETA)

$(BUILD)/shared_date_module.o: $(SRC)/shared_date_module.f90 $(INC) | $(BUILD)
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -I. -c -o $@ $<

$(BUILD)/tip_square_o_module.o: $(SRC)/tip_square_o_module.f90 $(BUILD)/shared_date_module.o | $(BUILD)
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -I. -c -o $@ $<

$(BUILD)/main.o: $(SRC)/main.f90 $(BUILD)/shared_date_module.o $(BUILD)/tip_square_o_module.o | $(BUILD)
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -I. -c -o $@ $<

$(INC):
	@echo "params.inc not found."
	@echo "  cp params.inc.template params.inc   # then edit"
	@echo "  bash scripts/beta.sh                 # or auto-generate"
	@exit 1

$(BUILD):
	mkdir -p $(BUILD) output

clean:
	rm -rf $(BUILD)/*.o $(BUILD)/*.mod $(BUILD)/mesh $(BUILD)/beta_function_wing

distclean: clean
	rm -rf $(BUILD) params.inc

.PHONY: all mesh beta debug clean distclean
