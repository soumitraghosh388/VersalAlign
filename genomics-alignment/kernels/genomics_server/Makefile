# genomics_server — top-level build driver (VCK5000, Vitis 2023.1)
# RESIDENT-server variant of genomics_dma: one launch serves nbatch batches
# via a DDR doorbell (see SERVER_NOTES.md). AIE graph + caps story identical.
#
#   make aie        # compile AIE graph (hw target) -> aie/build.hw/libadf.a
#   make pl         # synthesize the two server kernels -> pl/build.hw/*.xo
#   make link       # v++ --link  (libadf + 2 .xo + config) -> .xsa
#   make package    # v++ --package -> genomics_server.xclbin
#   make host       # build the XRT host
#   make all        # aie + pl + link + package + host
#   make run        # run on the board
#
# Cheap pre-hw validation:
#   make x86sim     # functional AIE sim (fast)
#   make aiesim     # cycle-approximate AIE sim (validates timing/stream sizes)

PLATFORM := xilinx_vck5000_gen4x8_qdma_2_202220_1
BUILD_DIR = build.hw
XCLBIN = $(BUILD_DIR)/genomics_server.xclbin
XSA    = $(BUILD_DIR)/genomics_server.xsa
LIBADF = aie/build.hw/libadf.a
XOS    = pl/build.hw/dist_server.xo pl/build.hw/collect_server.xo

.PHONY: all aie pl link package host run x86sim aiesim clean

all: aie pl link package host

aie:
	$(MAKE) -C aie TARGET=hw

pl:
	$(MAKE) -C pl TARGET=hw

x86sim:
	$(MAKE) -C aie TARGET=x86sim
	$(MAKE) -C aie sim TARGET=x86sim

aiesim:
	$(MAKE) -C aie TARGET=hw
	$(MAKE) -C aie sim TARGET=hw

link: $(XSA)
$(XSA): $(LIBADF) $(XOS) hw_link/config.cfg
	@mkdir -p $(BUILD_DIR)
	v++ -l -t hw --platform $(PLATFORM) \
		--config hw_link/config.cfg \
		$(LIBADF) $(XOS) \
		-o $(XSA)

package: $(XCLBIN)
$(XCLBIN): $(XSA) $(LIBADF)
	v++ -p -t hw --platform $(PLATFORM) \
		$(XSA) $(LIBADF) \
		-o $(XCLBIN)

host:
	$(MAKE) -C host

run:
	./host/genomics_server $(XCLBIN) $(NBATCH)

clean:
	$(MAKE) -C aie clean; $(MAKE) -C pl clean; $(MAKE) -C host clean
	rm -rf $(BUILD_DIR) *.log *.jou .Xil/
