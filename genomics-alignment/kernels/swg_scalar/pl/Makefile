PLATFORM := xilinx_vck5000_gen4x8_qdma_2_202220_1
TARGET := hw
FREQ := 300

BUILD_DIR = build.$(TARGET)
VPP_FLAGS = -t $(TARGET) --platform $(PLATFORM)
VPP_FLAGS += --kernel_frequency $(FREQ)

all: $(BUILD_DIR)/pl_mm2s.xo $(BUILD_DIR)/pl_s2mm.xo

$(BUILD_DIR)/pl_mm2s.xo: pl_mm2s.cpp
	v++ $(VPP_FLAGS) -k $(subst $(BUILD_DIR)/,,$(subst .xo,,$@)) -c -g -o $@ $^

$(BUILD_DIR)/pl_s2mm.xo: pl_s2mm.cpp
	v++ $(VPP_FLAGS) -k $(subst $(BUILD_DIR)/,,$(subst .xo,,$@)) -c -g -o $@ $^

clean:
	rm -rf $(BUILD_DIR) *.log *.jou .Xil/
