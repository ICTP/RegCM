# Global Makefile for RegCM v4

include Makefile.inc

HOST = @host@

printall::
	@echo ' this is Makefile for Regcm4 package: supported actions are: '
	@echo '  --- make terrain: to compile terrain' 
	@echo '  --- make icbc: to compile icbc codes'
	@echo '  ----make regcm:  to compile regcm code' 
	@echo '  ----make postnc:  to compile C++ NetCDF postprocessing code' 
	@echo '  ----make postproc:  to compile F90 postprocessing code' 
	@echo '  ----make clean:  to clean all ' 
	@echo '  ----make all: to compile all the codes and tools'

all: terrain icbc regcm postnc postproc

terrain:
	$(MAKE) -C $(REGCM_BASE_DIR)/PreProc/Terrain

icbc:
	$(MAKE) -C $(REGCM_BASE_DIR)/PreProc/ICBC

regcm:
	$(MAKE) -C $(REGCM_BASE_DIR)/PreProc/ICBC all TARGET=$(TARGET)

postnc:
	$(MAKE) -C $(REGCM_BASE_DIR)/PostNc

postproc:
	$(MAKE) -C $(REGCM_BASE_DIR)//PostProc all TARGET=$(TARGET)

clean:
	$(MAKE) -C $(REGCM_BASE_DIR)/PreProc/Terrain clean 
	$(MAKE) -C $(REGCM_BASE_DIR)/PreProc/ICBC clean 
	$(MAKE) -C $(REGCM_BASE_DIR)/Main clean 
	$(MAKE) -C $(REGCM_BASE_DIR)/PostNc clean 
	$(MAKE) -C $(REGCM_BASE_DIR)/PostProc clean 

DATAclean:
	cd $(REGCM_BASE_DIR)/PreProc/DATA && datalink.py -c 
