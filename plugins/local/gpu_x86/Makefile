# Makefile for module gpu_arch
# Generated automatically

MODULE_DIR = $(QP_ROOT)/plugins/local/gpu_x86
MODULE_NAME = gpu_arch

# Dependencies: 
DEPS = 

SOURCES = \
	$(QP_ROOT)/plugins/local/gpu_x86/gpu.c

gpu_arch.bin: IRPF90_temp/Makefile
	@cd IRPF90_temp && $(MAKE)

IRPF90_temp/Makefile: $(SOURCES) $(DEPS)
	@echo 'Running IRPF90 for $(MODULE_NAME)'
	@$(IRPF90) $(IRPF90_FLAGS)

clean:
	@rm -rf IRPF90_temp *.bin

.PHONY: clean
