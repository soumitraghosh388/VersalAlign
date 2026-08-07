FLAGS := -Wall -c -g -fmessage-length=0
FLAGS += -std=c++17 -Wno-unknown-pragmas -Wno-unused-label
FLAGS += -Wno-int-to-pointer-cast

INCLUDES += -I$(XILINX_VITIS)/aietools/include
INCLUDES += -I$(XILINX_VITIS)/include
INCLUDES += -I$(XILINX_XRT)/include

LIBS := -lxrt_core -lxrt_coreutil -lgcc -lc
LIBS += -lxilinxopencl -lpthread -lrt -ldl -lcrypt -lstdc++ -lOpenCL
LIBS += -Wno-unused-label -Wno-narrowing -DVERBOSE
LIBS += -L$(XILINX_XRT)/lib/

BUILD_DIR = build
EXECUTABLE = bitap.exe

all: host

.PHONY: host
host: $(BUILD_DIR) $(EXECUTABLE)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR);

$(EXECUTABLE): $(BUILD_DIR)/host.o
	g++ $^ $(LIBS) -std=c++17 -o $@
	@echo "COMPLETE: Host application $@ created."

$(BUILD_DIR)/host.o: host.cpp
	g++ $(FLAGS) $(INCLUDES) -o $@ $^

clean:
	rm -rf $(BUILD_DIR) *~ .run/ *.csv *.run_summary .Xil/ *.log *.jou $(EXECUTABLE)
