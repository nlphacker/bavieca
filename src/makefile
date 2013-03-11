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

include Makefile.defines
  
all:
	(cd common; $(MAKE) static)
	(cd common; $(MAKE) pic)
	(cd api; $(MAKE) all)
	(cd tools; $(MAKE) all)

cmd_tools:
	(cd common; $(MAKE) static)
	(cd tools; $(MAKE) all)

java-api: 
	(cd common; $(MAKE) pic)
	(cd api; $(MAKE) java-api)

clean:
	(cd common; $(MAKE) clean)
	(cd api; $(MAKE) clean)
	(cd tools; $(MAKE) clean)
