PLATFORM := xilinx_vck5000_gen4x8_qdma_2_202220_1
#TARGET := x86sim
TARGET := hw
# ---- Configuration (override on command line) ----
KERNEL ?= wfa
ROWS ?= 8
COLS ?= 50
BUF_SIZE ?= 2048
SCORE_BUF ?= 64
VERSION ?= v3m
ITERS ?= 1

# Kernel selection → compiler define
ifeq ($(KERNEL),wfa)
  KERNEL_DEF := USE_WFA
  STACKSIZE := 8000
else ifeq ($(KERNEL),bitap)
  KERNEL_DEF := USE_BITAP
  STACKSIZE := 6000
else ifeq ($(KERNEL),striped_swg)
  KERNEL_DEF := USE_STRIPED_SWG
  STACKSIZE := 6000
else ifeq ($(KERNEL),interseq_swg)
  KERNEL_DEF := USE_INTERSEQ_SWG
  STACKSIZE := 6000
else ifeq ($(KERNEL),scalar_swg)
  KERNEL_DEF := USE_SCALAR_SWG
  STACKSIZE := 6000

else
  KERNEL_DEF := USE_DUMMY
  STACKSIZE := 4000
endif

XPFM = $(shell platforminfo -p $(PLATFORM) --json="file")
#BUILD_DIR = build.$(KERNEL).$(VERSION).$(ROWS)x$(COLS).$(TARGET)
BUILD_DIR = build.$(TARGET)
SRC_DIR = $(shell readlink -f src/)
DATA_DIR = $(shell readlink -f data/)
SW_STRIPED_SRC    = $(shell readlink -f ../../sw_striped/aie/src/)
BITAP_SRC         = $(shell readlink -f ../../bitap/aie/src/)
WFA_SRC           = $(shell readlink -f ../../wfa/aie/src/)
INTERSEQ_SWG_SRC  = $(shell readlink -f ../../swg_interseq/aie/src/)
SCALAR_SWG_SRC  = $(shell readlink -f ../../swg_scalar/aie/src/)
COMMON            = $(shell readlink -f ../../common/)
AIE_FLAGS += --platform=$(XPFM)
AIE_FLAGS += -include="$(SRC_DIR)"
AIE_FLAGS += -include="$(SW_STRIPED_SRC)"
AIE_FLAGS += -include="$(BITAP_SRC)"
AIE_FLAGS += -include="$(WFA_SRC)"
AIE_FLAGS += -include="$(INTERSEQ_SWG_SRC)"
AIE_FLAGS += -include="$(SCALAR_SWG_SRC)"
AIE_FLAGS += -include="$(COMMON)"


GRAPH_CPP := $(SRC_DIR)/array_graph_v3m.cpp
DEPS := $(SRC_DIR)/array_graph_v3m.cpp $(SRC_DIR)/array_graph_v3m.hpp
DEPS += $(SRC_DIR)/tile_bottom_v3.hh $(SRC_DIR)/tile_bottom_v3.cc
DEPS += $(SRC_DIR)/tile_upper_v3.hh $(SRC_DIR)/tile_upper_v3.cc


all: $(BUILD_DIR)/libadf.a

$(BUILD_DIR)/kernel_config.h:
	@mkdir -p $(BUILD_DIR)
	@echo "// Auto-generated — do not edit" > $@
	@echo "#define NUM_ROWS $(ROWS)" >> $@
	@echo "#define NUM_COLS $(COLS)" >> $@
	@echo "#define BUF_SIZE $(BUF_SIZE)" >> $@
	@echo "#define SCORE_BUF $(SCORE_BUF)" >> $@
	@echo "#define $(KERNEL_DEF)" >> $@
	@echo "#define NUM_ITERS $(ITERS)" >> $@
	@echo "Generated $@ (KERNEL=$(KERNEL) $(VERSION) $(ROWS)x$(COLS) ITERS=$(ITERS))"

$(BUILD_DIR)/libadf.a: $(DEPS) $(BUILD_DIR)/kernel_config.h
	@echo "=== Building: KERNEL=$(KERNEL) $(VERSION) $(ROWS)x$(COLS) ITERS=$(ITERS) TARGET=$(TARGET) ==="
	cd $(BUILD_DIR); \
	aiecompiler -v --target=$(TARGET) \
		--stacksize=$(STACKSIZE) \
		-include="$(XILINX_VITIS)/aietools/include" \
		-include="$(SRC_DIR)" \
		-include="$(DATA_DIR)" \
		-include="$(BUILD_DIR)" \
		$(AIE_FLAGS) \
		$(GRAPH_CPP) \
		-workdir=work 2>&1 | tee aiecompiler.log

clean:
	rm -rf $(BUILD_DIR)

sim:
	@if [ $(TARGET) = "x86sim" ]; then \
		cd $(BUILD_DIR); \
		x86simulator --pkg-dir=work --i=.. ; \
	fi
	@if [ $(TARGET) = "hw" ]; then\
		cd $(BUILD_DIR); \
		aiesimulator --pkg-dir=work --i=.. --profile ; \
	fi
.PHONY: all clean sim
