#!/bin/bash

echo "building baviecaAPI for java"

# generate java files from the API header file
mkdir -p ./blt/bavieca
swig -c++ -java -package blt.bavieca -outdir ./blt/bavieca -I../src/api BaviecaAPI.i 

# move the wrapper file to the C++ api folder
cp BaviecaAPI_wrap.cxx ../src/api/BaviecaAPI_wrap.cxx

# compile java files
javac ./blt/bavieca/*.java

# create a jar file containing all the classes
jar -cf baviecaAPI.jar ./blt/bavieca/*.class

echo "baviecaAPI for java was successfully built!"
echo "note: library for native API needs to be rebuilt with the wrapper"

