#----------------------------------------------------------------------------------------------#
#  Copyright (C) 2012 Daniel Bolaños - www.bltek.com - Boulder Language Technologies           #
#                                                                                              #
#  www.bavieca.org is the website of the Bavieca Speech Recognition Toolkit                    #
#                                                                                              #
#  Licensed under the Apache License, Version 2.0 (the "License");                             #
#  you may not use this file except in compliance with the License.                            #
#  You may obtain a copy of the License at                                                     #
#                                                                                              #
#          http://www.apache.org/licenses/LICENSE-2.0                                          #
#                                                                                              #
#  Unless required by applicable law or agreed to in writing, software                         #
#  distributed under the License is distributed on an "AS IS" BASIS,                           #
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.                    #
#  See the License for the specific language governing permissions and                         #
#  limitations under the License.                                                              #
# ---------------------------------------------------------------------------------------------#
  
include ../Makefile.defines

INC = -I../common/algebra -I../common/alignment -I../common/audio -I../common/base -I../common/decoding \
      -I../common/dynamicdecoder -I../common/estimation -I../common/config -I../common/hmm -I../common/io \
      -I../common/other -I../common/param -I../common/text -I../common/sadmodule -I../common/transform -I../common/vtlestimator \
      -I../common/wfsabuilder -I../common/wfsadecoder ${INCS_DIR_CBLAS} ${INCS_DIR_LAPACK} 

SRC_DIR = .
OBJ_DIR = ../../obj/$(ARCH)-$(OS)/api
LIB_DIR = ../../lib/$(ARCH)-$(OS)
OBJFILES_BASE = $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%,$(patsubst %.cpp,%.o,$(wildcard $(SRC_DIR)/*.cpp)))
OBJFILES_BASE_JAVA = $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%,$(patsubst %.cxx,%.o,$(wildcard $(SRC_DIR)/*.cxx)))
LIBS = -L../../lib/$(ARCH)-$(OS)/ $(LIBS_DIR_CBLAS) $(LIBS_DIR_LAPACK)

all: createDirectories libbaviecaapi.a libbaviecaapi.so libbaviecaapijni.so
	@echo $(OBJFILES_BASE)

java-api: createDirectories libbaviecaapijni.so
	@echo $(OBJFILES_BASE)

clean: 
	rm -rf $(OBJ_DIR)
	rm -rf $(LIB_DIR)/libbaviecaapi.a
	rm -rf $(LIB_DIR)/libbaviecaapi.so
	rm -rf $(LIB_DIR)/libbaviecaapijni.so

createDirectories: 
	(mkdir -p $(LIB_DIR))
	(mkdir -p $(OBJ_DIR))

# ----------------------------------------------
# create the library
#-----------------------------------------------

libbaviecaapi.a: $(OBJFILES_BASE)
	$(AR) $@ $?
	(mv libbaviecaapi.a $(LIB_DIR))

libbaviecaapi.so: $(OBJFILES_BASE)
	$(XCC) $(LIBS) -shared -o libbaviecaapi.so $(OBJFILES_BASE) -lcommon_pic
	(mv libbaviecaapi.so $(LIB_DIR))

libbaviecaapijni.so: $(OBJFILES_BASE) $(OBJFILES_BASE_JAVA)
	$(XCC) $(LIBS) -shared -o libbaviecaapijni.so $(OBJFILES_BASE) $(OBJFILES_BASE_JAVA) -lcommon_pic
	(mv libbaviecaapijni.so $(LIB_DIR))

# ----------------------------------------------
# create the object files from the source files
# ----------------------------------------------

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(XCC) $(CPPFLAGS_SHARED) $(INC) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cxx
	$(XCC) $(CPPFLAGS_SHARED) $(INC) $(INCS_DIR_JNI) -fno-strict-aliasing -c $< -o $@





