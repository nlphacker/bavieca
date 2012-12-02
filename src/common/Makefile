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

INC = -I./algebra -I./alignment -I./audio -I./base -I./decoding -I./dynamicdecoder -I./estimation -I./config -I./hmm \
-I./io -I./other -I./param -I./text -I./transform -I./wfsadecoder ${INCS_DIR_CBLAS} ${INCS_DIR_LAPACK} 

SRC_DIR = ../common
OBJ_DIR = ../../obj/$(ARCH)-$(OS)/common
LIB_DIR = ../../lib/$(ARCH)-$(OS)
OBJFILES_BASE = $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%,$(patsubst %.cpp,%.o,$(wildcard $(SRC_DIR)/*/*.cpp)))

all: createDirectories libcommon.a

clean: 
	rm -rf $(OBJ_DIR)
	rm -rf $(LIB_DIR)/libcommon.a

createDirectories: 
	(mkdir -p $(LIB_DIR))
	(mkdir -p $(OBJ_DIR)/algebra)
	(mkdir -p $(OBJ_DIR)/alignment)
	(mkdir -p $(OBJ_DIR)/audio)
	(mkdir -p $(OBJ_DIR)/base)
	(mkdir -p $(OBJ_DIR)/config)
	(mkdir -p $(OBJ_DIR)/decoding)
	(mkdir -p $(OBJ_DIR)/dynamicdecoder)
	(mkdir -p $(OBJ_DIR)/estimation)
	(mkdir -p $(OBJ_DIR)/hmm)
	(mkdir -p $(OBJ_DIR)/io)
	(mkdir -p $(OBJ_DIR)/other)
	(mkdir -p $(OBJ_DIR)/param)
	(mkdir -p $(OBJ_DIR)/sadmodule)
	(mkdir -p $(OBJ_DIR)/text)
	(mkdir -p $(OBJ_DIR)/transform)
	(mkdir -p $(OBJ_DIR)/vtlestimator)
	(mkdir -p $(OBJ_DIR)/wfsabuilder)
	(mkdir -p $(OBJ_DIR)/wfsadecoder)

# ----------------------------------------------
# create the library
#-----------------------------------------------

libcommon.a: $(OBJFILES_BASE)
	$(AR) $@ $?
	(mv libcommon.a $(LIB_DIR))

# ----------------------------------------------
# create the object files from the source files
# ----------------------------------------------

$(OBJ_DIR)/algebra/%.o: $(SRC_DIR)/algebra/%.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/alignment/%.o: $(SRC_DIR)/alignment/%.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/audio/%.o: $(SRC_DIR)/audio/%.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/base/%.o: $(SRC_DIR)/base/%.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/config/%.o: $(SRC_DIR)/config/%.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/decoding/%.o: $(SRC_DIR)/decoding/%.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/dynamicdecoder/%.o: $(SRC_DIR)/dynamicdecoder/%.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/estimation/%.o: $(SRC_DIR)/estimation/%.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/hmm/%.o: $(SRC_DIR)/hmm/%.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/io/%.o: $(SRC_DIR)/io/%.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/other/%.o: $(SRC_DIR)/other/%.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/param/%.o: $(SRC_DIR)/param/%.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/sadmodule/%.o: $(SRC_DIR)/sadmodule/%.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/text/%.o: $(SRC_DIR)/text/%.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/transform/%.o: $(SRC_DIR)/transform/%.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/vtlestimator/%.o: $(SRC_DIR)/vtlestimator/%.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/wfsabuilder/%.o: $(SRC_DIR)/wfsabuilder/%.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/wfsadecoder/%.o: $(SRC_DIR)/wfsadecoder/%.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@









