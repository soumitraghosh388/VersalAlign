PLATFORM := xilinx_vck5000_gen4x8_qdma_2_202220_1
TARGET := hw

XPFM = $(shell platforminfo -p $(PLATFORM) --json="file")
BUILD_DIR = build.$(TARGET)
WORK_DIR = work
SRC_DIR = $(shell readlink -f src/)
DATA_DIR = $(shell readlink -f data/)

DEPS := $(SRC_DIR)/graph.cpp
DEPS += $(SRC_DIR)/kernels.hpp
DEPS += $(SRC_DIR)/graph.hpp
DEPS += $(SRC_DIR)/swg_kernel.cc
COMMON_DIR = $(shell readlink -f ../../common/)
AIE_FLAGS += --platform=$(XPFM)
AIE_FLAGS += -include="$(COMMON_DIR)"

all: $(BUILD_DIR)/libadf.a

$(BUILD_DIR)/libadf.a: $(DEPS)
	@mkdir -p $(BUILD_DIR);
	cd $(BUILD_DIR); \
	aiecompiler -v --target=$(TARGET) \
		--stacksize=6000 \
		-include="$(XILINX_VITIS)/aietools/include" \
		-include="$(SRC_DIR)" \
		-include="$(DATA_DIR)" \
		$(AIE_FLAGS) \
		$(SRC_DIR)/graph.cpp \
		-workdir=$(WORK_DIR) 2>&1 | tee aiecompiler.log

clean:
	rm -rf $(BUILD_DIR)

sim:
	@if [ $(TARGET) = "x86sim" ]; then \
		cd $(BUILD_DIR); \
		x86simulator --pkg-dir=$(WORK_DIR) --i=.. ; \
	fi
	@if [ $(TARGET) = "hw" ]; then \
		cd $(BUILD_DIR); \
		aiesimulator --pkg-dir=$(WORK_DIR) --i=.. --profile ; \
	fi
