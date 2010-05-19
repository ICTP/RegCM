# Global Makefile for RegCM v4

include Makefile.inc

HOST = @host@


printall:
	@echo ' this is Makefile for Regcm4 package: supported actions are: '
	@echo '  --- make terrain: to compile terrain' 
	@echo '  --- make icbc: to compile icbc codes'
	@echo '  ----make regcm:  to compile regcm code' 
	@echo '  ----make clean:  to clean all ' 
	@echo '  ----make all: to compile all the codes and tools'

all: terrain icbc regcm

terrain:
	#Move to the Terrain directory and run make there
	cd $(REGCM_BASE_DIR)/PreProc/Terrain && $(MAKE) 


icbc:
	#Move to the ICBC directory and run make there
	cd $(REGCM_BASE_DIR)/PreProc/ICBC && $(MAKE) 

regcm:
	# Change to the Main directory and run make there
	cd $(REGCM_BASE_DIR)/Main && $(MAKE) all TARGET=$(TARGET)

postproc:
	cd $(REGCM_BASE_DIR)/PostProc && $(MAKE) all TARGET=$(TARGET)


clean:
	cd $(REGCM_BASE_DIR)/PreProc/Terrain && $(MAKE) clean 
	cd $(REGCM_BASE_DIR)/PreProc/ICBC && $(MAKE) clean
	cd $(REGCM_BASE_DIR)/Main && $(MAKE) clean
	cd $(REGCM_BASE_DIR)/PostProc && $(MAKE) clean

DATAclean:
	cd $(REGCM_BASE_DIR)/PreProc/DATA && datalink.py -c 
