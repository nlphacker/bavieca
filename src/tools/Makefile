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

OBJ_DIR = ../../obj/$(ARCH)-$(OS)/tools
LIB_DIR = ../../lib/$(ARCH)-$(OS)
BIN_DIR = ../../bin/$(ARCH)-$(OS)

LIBS = -L../../lib/$(ARCH)-$(OS)/ $(LIBS_DIR_CBLAS) $(LIBS_DIR_LAPACK)  

all: createDirectories aligner contextclustering dtaccumulator dtestimator dynamicdecoder fmllrestimator gmmeditor \
     hldaestimator hmminitializer hmmx latticeeditor ldaestimator mapestimator mlaccumulator mlestimator mllrestimator \
     param paramx regtree sadmodule vtlestimator wfsabuilder wfsadecoder

clean: 
	rm -rf $(OBJ_DIR)
	rm -rf $(BIN_DIR)
	rm -rf ../../obj
	rm -rf ../../lib
	rm -rf ../../bin

createDirectories: 
	(mkdir -p $(OBJ_DIR))
	(mkdir -p $(BIN_DIR))

# ----------------------------------------------
# create the tools
#-----------------------------------------------

aligner: $(OBJ_DIR)/mainAligner.o
	g++ $(CPPFLAGS) $(LIBS) -o $(BIN_DIR)/aligner $(OBJ_DIR)/mainAligner.o -lcommon ${LIB_LAPACK} ${LIB_CBLAS}  

contextclustering: $(OBJ_DIR)/mainContextClustering.o
	g++ $(CPPFLAGS) $(LIBS) -o $(BIN_DIR)/contextclustering $(OBJ_DIR)/mainContextClustering.o -lcommon ${LIB_LAPACK} ${LIB_CBLAS} 

dtaccumulator: $(OBJ_DIR)/mainDTAccumulator.o
	g++ $(CPPFLAGS) $(LIBS) -o $(BIN_DIR)/dtaccumulator $(OBJ_DIR)/mainDTAccumulator.o -lcommon ${LIB_LAPACK} ${LIB_CBLAS} 

dtestimator: $(OBJ_DIR)/mainDTEstimator.o
	g++ $(CPPFLAGS) $(LIBS) -o $(BIN_DIR)/dtestimator $(OBJ_DIR)/mainDTEstimator.o -lcommon ${LIB_LAPACK} ${LIB_CBLAS} 

dynamicdecoder: $(OBJ_DIR)/mainDynamicDecoder.o
	g++ $(CPPFLAGS) $(LIBS) -o $(BIN_DIR)/dynamicdecoder $(OBJ_DIR)/mainDynamicDecoder.o -lcommon ${LIB_LAPACK} ${LIB_CBLAS} 

fmllrestimator: $(OBJ_DIR)/mainfMLLREstimator.o
	g++ $(CPPFLAGS) $(LIBS) -o $(BIN_DIR)/fmllrestimator $(OBJ_DIR)/mainfMLLREstimator.o -lcommon ${LIB_LAPACK} ${LIB_CBLAS} 

gmmeditor: $(OBJ_DIR)/mainGMMEditor.o
	g++ $(CPPFLAGS) $(LIBS) -o $(BIN_DIR)/gmmeditor $(OBJ_DIR)/mainGMMEditor.o -lcommon ${LIB_LAPACK} ${LIB_CBLAS}

hldaestimator: $(OBJ_DIR)/mainHLDAEstimator.o
	g++ $(CPPFLAGS) $(LIBS) -o $(BIN_DIR)/hldaestimator $(OBJ_DIR)/mainHLDAEstimator.o -lcommon ${LIB_LAPACK} ${LIB_CBLAS} 

hmminitializer: $(OBJ_DIR)/mainHMMInitializer.o
	g++ $(CPPFLAGS) $(LIBS) -o $(BIN_DIR)/hmminitializer $(OBJ_DIR)/mainHMMInitializer.o -lcommon ${LIB_LAPACK} ${LIB_CBLAS}

hmmx: $(OBJ_DIR)/mainHMMX.o
	g++ $(CPPFLAGS) $(LIBS) -o $(BIN_DIR)/hmmx $(OBJ_DIR)/mainHMMX.o -lcommon ${LIB_LAPACK} ${LIB_CBLAS}

latticeeditor: $(OBJ_DIR)/mainLatticeEditor.o
	$(XCC) $(CPPFLAGS) $(LIBS) -o $(BIN_DIR)/latticeeditor $(OBJ_DIR)/mainLatticeEditor.o -lcommon ${LIB_LAPACK} ${LIB_CBLAS}

ldaestimator: $(OBJ_DIR)/mainLDAEstimator.o
	$(XCC) $(CPPFLAGS) $(LIBS) -o $(BIN_DIR)/ldaestimator $(OBJ_DIR)/mainLDAEstimator.o -lcommon ${LIB_LAPACK} ${LIB_CBLAS}

mapestimator: $(OBJ_DIR)/mainMAPEstimator.o
	$(XCC) $(CPPFLAGS) $(LIBS) -o $(BIN_DIR)/mapestimator $(OBJ_DIR)/mainMAPEstimator.o -lcommon ${LIB_LAPACK} ${LIB_CBLAS}

mlaccumulator: $(OBJ_DIR)/mainMLAccumulator.o
	$(XCC) $(CPPFLAGS) $(LIBS) -o $(BIN_DIR)/mlaccumulator $(OBJ_DIR)/mainMLAccumulator.o -lcommon ${LIB_LAPACK} ${LIB_CBLAS} 

mlestimator: $(OBJ_DIR)/mainMLEstimator.o
	$(XCC) $(CPPFLAGS) $(LIBS) -o $(BIN_DIR)/mlestimator $(OBJ_DIR)/mainMLEstimator.o -lcommon ${LIB_LAPACK} ${LIB_CBLAS}

mllrestimator: $(OBJ_DIR)/mainMLLREstimator.o
	$(XCC) $(CPPFLAGS) $(LIBS) -o $(BIN_DIR)/mllrestimator $(OBJ_DIR)/mainMLLREstimator.o -lcommon ${LIB_LAPACK} ${LIB_CBLAS} 

param: $(OBJ_DIR)/mainParam.o
	$(XCC) $(CPPFLAGS) $(LIBS) -o $(BIN_DIR)/param $(OBJ_DIR)/mainParam.o -lcommon ${LIB_LAPACK} ${LIB_CBLAS} 

paramx: $(OBJ_DIR)/mainParamX.o
	$(XCC) $(CPPFLAGS) $(LIBS) -o $(BIN_DIR)/paramx $(OBJ_DIR)/mainParamX.o -lcommon ${LIB_LAPACK} ${LIB_CBLAS} 

regtree: $(OBJ_DIR)/mainRegTree.o
	$(XCC) $(CPPFLAGS) $(LIBS) -o $(BIN_DIR)/regtree $(OBJ_DIR)/mainRegTree.o -lcommon ${LIB_LAPACK} ${LIB_CBLAS}   

sadmodule: $(OBJ_DIR)/mainSADModule.o
	$(XCC) $(CPPFLAGS) $(LIBS) -o $(BIN_DIR)/sadmodule $(OBJ_DIR)/mainSADModule.o -lcommon ${LIB_LAPACK} ${LIB_CBLAS}

vtlestimator: $(OBJ_DIR)/mainVTLEstimator.o
	$(XCC) $(CPPFLAGS) $(LIBS) -o $(BIN_DIR)/vtlestimator $(OBJ_DIR)/mainVTLEstimator.o -lcommon ${LIB_LAPACK} ${LIB_CBLAS} 

wfsabuilder: $(OBJ_DIR)/mainWFSABuilder.o
	$(XCC) $(CPPFLAGS) $(LIBS) -o $(BIN_DIR)/wfsabuilder $(OBJ_DIR)/mainWFSABuilder.o -lcommon ${LIB_LAPACK} ${LIB_CBLAS} 

wfsadecoder: $(OBJ_DIR)/mainWFSADecoder.o
	$(XCC) $(CPPFLAGS) $(LIBS) -o $(BIN_DIR)/wfsadecoder $(OBJ_DIR)/mainWFSADecoder.o -lcommon ${LIB_LAPACK} ${LIB_CBLAS} 


# ----------------------------------------------
# create the object files from the source files
# ----------------------------------------------

$(OBJ_DIR)/mainAligner.o: ./aligner/mainAligner.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/mainContextClustering.o: ./contextclustering/mainContextClustering.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/mainDTAccumulator.o: ./dtaccumulator/mainDTAccumulator.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/mainDTEstimator.o: ./dtestimator/mainDTEstimator.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/mainDynamicDecoder.o: ./dynamicdecoder/mainDynamicDecoder.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/mainfMLLREstimator.o: ./fmllrestimator/mainfMLLREstimator.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/mainGMMEditor.o: ./gmmeditor/mainGMMEditor.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/mainHLDAEstimator.o: ./hldaestimator/mainHLDAEstimator.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/mainHMMInitializer.o: ./hmminitializer/mainHMMInitializer.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/mainHMMX.o: ./hmmx/mainHMMX.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/mainLatticeEditor.o: ./latticeeditor/mainLatticeEditor.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/mainLDAEstimator.o: ./ldaestimator/mainLDAEstimator.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/mainMAPEstimator.o: ./mapestimator/mainMAPEstimator.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/mainMLAccumulator.o: ./mlaccumulator/mainMLAccumulator.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/mainMLEstimator.o: ./mlestimator/mainMLEstimator.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/mainMLLREstimator.o: ./mllrestimator/mainMLLREstimator.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/mainParam.o: ./param/mainParam.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/mainParamX.o: ./paramx/mainParamX.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/mainRegTree.o: ./regtree/mainRegTree.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/mainSADModule.o: ./sadmodule/mainSADModule.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/mainVTLEstimator.o: ./vtlestimator/mainVTLEstimator.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/mainWFSADecoder.o: ./wfsadecoder/mainWFSADecoder.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@

$(OBJ_DIR)/mainWFSABuilder.o: ./wfsabuilder/mainWFSABuilder.cpp
	$(XCC) $(CPPFLAGS) $(INC) -c $< -o $@










