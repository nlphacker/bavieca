/*---------------------------------------------------------------------------------------------*
 * Copyright (C) 2012 Daniel Bolaños - www.bltek.com - Boulder Language Technologies           *
 *                                                                                             *
 * www.bavieca.org is the website of the Bavieca Speech Recognition Toolkit                    *
 *                                                                                             *
 * Licensed under the Apache License, Version 2.0 (the "License");                             *
 * you may not use this file except in compliance with the License.                            *
 * You may obtain a copy of the License at                                                     *
 *                                                                                             *
 *         http://www.apache.org/licenses/LICENSE-2.0                                          *
 *                                                                                             *
 * Unless required by applicable law or agreed to in writing, software                         *
 * distributed under the License is distributed on an "AS IS" BASIS,                           *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.                    *
 * See the License for the specific language governing permissions and                         *
 * limitations under the License.                                                              *
 *---------------------------------------------------------------------------------------------*/
 
import java.lang.*;
import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ShortBuffer;
import java.nio.ByteOrder;
import blt.bavieca.*;

public class exampleAlign {
	
	static {
		try {
			System.loadLibrary("baviecaapijni");
		} catch (UnsatisfiedLinkError e) {
			System.err.println("Native code library failed to load.\n" + e);
			System.exit(1);
		}
	}

	public static void main(String argv[]) {
	
		// check parameters
		if (argv.length != 3) {
			System.out.println("usage: exampleAlign [configurationFile] [rawAudioFile] [txtTranscriptionFile]");
			return;
		}

		System.out.println("[java] - example begins --------------");
		blt.bavieca.BaviecaAPI baviecaAPI = new blt.bavieca.BaviecaAPI(argv[0]);

		// INIT_SAD, INIT_ALIGNER, INIT_DECODER, INIT_ADAPTATION
		if (baviecaAPI.initialize((short)BaviecaAPI_SWIGConstants.INIT_ALIGNER) == true) {
			System.out.println("[java] library initialized successfully");
		}		
		
		try {
		
			FileInputStream fis;		
		
			// open the raw audio file (16KHz 16bits)
			String strFileRaw = argv[1];		
			fis = new FileInputStream(strFileRaw);  	
			int iBytes = fis.available();
			byte [] bytes = new byte[iBytes];		

			// read bytes of audio
			if (fis.read(bytes) != iBytes) {
				System.out.println("Error: unable to read from file");
				return;
			}
			System.out.println("[java] audio read: " + ((float)bytes.length)/32000.0 + " seconds");
			
			// turn bytes into shorts (little endian)
			short [] sSamples = new short[iBytes/2]; 
			ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN).asShortBuffer().get(sSamples);	
			      		
			// extract features
			long lFeatures [] = new long[1];
	     	float [] fFeatures = baviecaAPI.extractFeatures(sSamples,(long)sSamples.length,lFeatures);
			System.out.println("[java] features extracted: " + lFeatures[0]);	     	
      	      
      	// forced alignment (Viterbi alignment)
			BufferedReader br = new BufferedReader(new FileReader(argv[2]));
			StringBuffer str = new StringBuffer();
			String line = br.readLine();
			while (line != null) {
				str.append(line);
				str.append("\n");
				line = br.readLine();
			}
			String strWords = str.toString();
     		AlignmentI alignment = baviecaAPI.align(fFeatures,lFeatures[0],strWords,false);
     		if (alignment != null) {
     			for(int i=0 ; i < alignment.size() ; ++i) {
     				WordAlignmentI wordAlignment = alignment.getWordAlignment(i);
     				System.out.format("[java] word: (%5d) %-12s (%5d)%n",wordAlignment.getFrameStart(),wordAlignment.getWord(),wordAlignment.getFrameEnd());
     			}
     		} else {
     			System.out.println("[java] alignment is null");     		
     		}
     		alignment.delete();
     		
  	  		// close the file
  			fis.close();

    	} catch (EOFException eof) {
			System.out.println("EOF reached"); 
		} catch (IOException ioe) {
			System.out.println("IO error: " + ioe);
		}	
		
		// uninitialization
		baviecaAPI.uninitialize();
		baviecaAPI.delete();
		System.out.println("[java] library uninitialized successfully");
		System.out.println("[java] - example ends ----------------");
	}
}

