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

public class exampleDecode {
	
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
		if (argv.length != 2) {
			System.out.println("usage: exampleDecode [configurationFile] [rawAudioFile]");
			return;
		}

		System.out.println("[java] - example begins --------------");
		blt.bavieca.BaviecaAPI baviecaAPI = new blt.bavieca.BaviecaAPI(argv[0]);

		// INIT_SAD, INIT_ALIGNER, INIT_DECODER, INIT_ADAPTATION
		if (baviecaAPI.initialize((short)BaviecaAPI_SWIGConstants.INIT_DECODER) == true) {
			System.out.println("[java] library initialized successfully");
		}
		
		String strFileRaw = argv[1];
		
		FileInputStream fis;
		
		try {
		
			// open the raw audio file (16KHz 16bits)
			fis = new FileInputStream(strFileRaw);  			

			// signal the beginning of an utterance         
			baviecaAPI.decBeginUtterance();

			// simulate live audio (chunks of audio are read from a file and processed one by one)
			int iBytesChunk = 3200;		// one 10th of a second
			byte [] bytesChunk = new byte[iBytesChunk];
			short [] sSamplesChunk = new short[iBytesChunk/2];

			int iBytesRead = 0;
			int iAudioProcessed = 0;
			do {			
				
				iBytesRead = fis.read(bytesChunk);
				if (iBytesRead == 0) {
					break;
				}
			
				// turn bytes into shorts (little endian) 
				ByteBuffer.wrap(bytesChunk).order(ByteOrder.LITTLE_ENDIAN).asShortBuffer().get(sSamplesChunk);	
			      		
				long lFeatures [] = new long[1];				float [] fFeatures = baviecaAPI.extractFeatures(sSamplesChunk,(long)iBytesChunk/2,lFeatures);
				iAudioProcessed += bytesChunk.length;
				System.out.println("[java] audio processed: " + ((float)iAudioProcessed)/32000.0 + " seconds");      		
				baviecaAPI.decProcess(fFeatures,lFeatures[0]);   	  		

			} while(iBytesRead == iBytesChunk);
			System.out.println("[java] audio processed: " + ((float)iAudioProcessed)/32000.0 + " seconds");

			// get the hypothesis
  			HypothesisI hypothesis = baviecaAPI.decGetHypothesis();
  			if (hypothesis != null) {
  				System.out.println("[java] hypothesis contains " + hypothesis.size() + " words");
    			for(int i=0 ; i<hypothesis.size() ; ++i) {
 	  				WordHypothesisI wordHypothesisI = hypothesis.getWordHypothesis(i);
  					System.out.format("[java] hyp: (%5d) %-12s (%5d)%n",wordHypothesisI.getFrameStart(),
 						wordHypothesisI.getWord(),wordHypothesisI.getFrameEnd());
  				}
    			hypothesis.delete();
 	  		} else {
  				System.out.println("[java] null hypothesis!");	
  			}

			// signal the end of utterance      				baviecaAPI.decEndUtterance();			
	     	
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

