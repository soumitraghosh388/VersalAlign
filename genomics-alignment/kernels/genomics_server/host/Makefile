FLAGS := -Wall -c -g -fmessage-length=0 -std=c++17 -Wno-unknown-pragmas -Wno-unused-label -Wno-int-to-pointer-cast
INCLUDES += -I$(XILINX_VITIS)/aietools/include -I$(XILINX_VITIS)/include -I$(XILINX_XRT)/include
LIBS := -lxrt_core -lxrt_coreutil -lgcc -lc -lxilinxopencl -lpthread -lrt -ldl -lcrypt -lstdc++ -lOpenCL -Wno-unused-label -Wno-narrowing -DVERBOSE -L$(XILINX_XRT)/lib/
BUILD_DIR = build
EXECUTABLE = genomics_server

all: host
.PHONY: host
host: $(BUILD_DIR) $(EXECUTABLE)
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)
$(EXECUTABLE): $(BUILD_DIR)/host_genomics_server.o
	g++ $^ $(LIBS) -std=c++17 -o $@
	@echo "COMPLETE: $@"
$(BUILD_DIR)/host_genomics_server.o: host_genomics_server.cpp
	g++ $(FLAGS) $(INCLUDES) -o $@ $^
clean:
	rm -rf $(BUILD_DIR) $(EXECUTABLE)
