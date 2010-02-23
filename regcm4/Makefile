# Global Makefile for RegCM v4 - credit to UCSC guys

include Makefile.inc

TARGET = AMD64
HOST = @host@
REGCM_NEST_DIR = 
NENSEMBLES = 
REMOVALSTR = --remove

ifeq ($(NENSEMBLES),)
  ENSFLAG=
else
  ENSFLAG="--ensembles=$(NENSEMBLES)"
endif

#all: terrain icbc regcm
printall:
	@echo ' this is Makefile for Regcm4 package: supported actions are: '
	@echo '  --- make data: to create DATA link '
	@echo '  --- make terrain: to compile terrain' 
	@echo '  --- make run_terrain: to execute terrain code' 
	@echo '  --- make icbc: to compile icbc codes'
	@echo '  ----make run_icbc: to run icbc codes'
	@echo '  ----make regcm:  to compile regcm code' 
	@echo '  ' 

data: 
	#link the Data directory to PreProc/DATA
	cd $(REGCM_BASE_DIR)/PreProc/DATA && ./datalinker.py -d  $(REGCM_DATA_DIR) 

terrain:

	#link domain.param to Terrain
	-ln -sf $(REGCM_BASE_DIR)/domain.param $(REGCM_BASE_DIR)/PreProc/Terrain/domain.param
	#Move to the Terrain directory and run make there
	cd $(REGCM_BASE_DIR)/PreProc/Terrain && $(MAKE) 

icbc:
	#Move to the ICBC directory and run make there
	cd $(REGCM_BASE_DIR)/PreProc/ICBC && $(MAKE) all TARGET=$(TARGET)

regcm:
	# Change to the Main directory and run make there
	cd $(REGCM_BASE_DIR)/Main && $(MAKE) all TARGET=$(TARGET)

postproc:
	cd $(REGCM_BASE_DIR)/PostProc && $(MAKE) all TARGET=$(TARGET)


$(ICBC_DEPS): terrain	

$(REGCM_DEPS): icbc	

clean:
	cd $(REGCM_BASE_DIR)/PreProc/Terrain && $(MAKE) clean 
	cd $(REGCM_BASE_DIR)/PreProc/ICBC && $(MAKE) clean
	cd $(REGCM_BASE_DIR)/Main && $(MAKE) clean
	cd $(REGCM_BASE_DIR)/PostProc && $(MAKE) clean

#dataclean:
#	cd $(REGCM_BASE_DIR)/PreProc/Terrain && $(MAKE) dataclean 
#	cd $(REGCM_BASE_DIR)/PreProc/ICBC && $(MAKE) dataclean
#	-rm -rf Input.*
