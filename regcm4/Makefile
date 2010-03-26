# Global Makefile for RegCM v4

include Makefile.inc

HOST = @host@


#all: terrain icbc regcm
printall:
	@echo ' this is Makefile for Regcm4 package: supported actions are: '
	@echo '  --- make data: to create DATA link '
	@echo '  --- make terrain: to compile terrain' 
	@echo '  --- make run_terrain: to execute terrain code' 
	@echo '  --- make icbc: to compile icbc codes'
	@echo '  ----make run_icbc: to run icbc codes'
	@echo '  ----make regcm:  to compile regcm code' 
	@echo '  ----make check : to check your configuration (run this if any of the steps above fails)' 

check: 
	# check compilation flags and how MACROs are defined in Makefile.inc 
	utils/check_make.py 


data: 
	#link the Data directory to PreProc/DATA
	cd $(REGCM_BASE_DIR)/PreProc/DATA && ./datalinker.py -d  $(REGCM_DATA_DIR) 

terrain:

	#Move to the Terrain directory and run make there
	cd $(REGCM_BASE_DIR)/PreProc/Terrain && $(MAKE) 

run_terrain:
	# execute the terrain command 
	cd $(REGCM_BASE_DIR)/PreProc/Terrain && ./terrain  

icbc:
	#Move to the ICBC directory and run make there
	cd $(REGCM_BASE_DIR)/PreProc/ICBC && $(MAKE) 

run_icbc:
	# Run icbc.x script 
	cd $(REGCM_BASE_DIR)/PreProc/ICBC && csh ./icbc_v4.x

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

DATAclean:
	cd $(REGCM_BASE_DIR)/PreProc/DATA && datalink.py -c 
