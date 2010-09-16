# Global Makefile for RegCM v4

include Makefile.inc

HOST = @host@

printall::
	@echo ' this is Makefile for Regcm4 package: supported actions are: '
	@echo '  --- make terrain: to compile terrain' 
	@echo '  --- make icbc: to compile icbc codes'
	@echo '  --- make clm2rcm: to compile clm24cm code'
	@echo '  --- make regcm:  to compile main regcm code' 
	@echo '  --- make postproc:  to compile F90 postprocessing code' 
	@echo '  --- make clean:  to clean all ' 
	@echo '  --- make all: to compile all the codes and tools'

all: terrain icbc clm2rcm regcm postproc
	@echo "##############################################################"
	@echo "YOU HAVE DONE IT! YOU HAVE COMPILED THE MODEL!"
	@echo "LET'S START PLAYING WITH THE BEAST....                        "
	@echo "##############################################################"

terrain:
	@$(MAKE) -C $(REGCM_BASE_DIR)/PreProc/Terrain
	@echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
	@echo "SUCCESSFULLY CREATED TERRAIN BINARY FILE"
	@echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"

icbc:
	@$(MAKE) -C $(REGCM_BASE_DIR)/PreProc/ICBC
	@echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
	@echo "SUCCESSFULLY CREATED ICBC BINARY FILES"
	@echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"

clm2rcm:
ifeq ($(USECLM),1)
	@$(MAKE) -C $(REGCM_BASE_DIR)/PreProc/CLM
	@echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
	@echo "SUCCESSFULLY CREATED clm2rcm BINARY FILES"
	@echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
else
	@echo "CLM option is not activated. Skipping creation of clm2rcm."
endif

regcm:
	@$(MAKE) -C $(REGCM_BASE_DIR)/Main all TARGET=$(TARGET)
	@echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
	@echo "SUCCESSFULLY CREATED MODEL BINARY FILE"
	@echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"

postproc:
	@$(MAKE) -C $(REGCM_BASE_DIR)//PostProc all TARGET=$(TARGET)
	@echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
	@echo "SUCCESSFULLY CREATED F90 POSTPROC BINARY FILES"
	@echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"

clean:
	@$(MAKE) -C $(REGCM_BASE_DIR)/PreProc/Terrain clean 
	@$(MAKE) -C $(REGCM_BASE_DIR)/PreProc/ICBC clean 
	@$(MAKE) -C $(REGCM_BASE_DIR)/Main clean 
	@$(MAKE) -C $(REGCM_BASE_DIR)/PostProc clean 
	@rm -f Bin/*
	@echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
	@echo "SUCCESSFULLY CLEANED ALL BINARY FILES"
	@echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"

DATAclean:
	cd $(REGCM_BASE_DIR)/PreProc/DATA && datalink.py -c 
	@echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
	@echo "SUCCESSFULLY CLEANED ALL DATA FILES"
	@echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
