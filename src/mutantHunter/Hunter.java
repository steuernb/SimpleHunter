package mutantHunter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;
import java.util.zip.GZIPInputStream;

import support.MPileupLine;


public class Hunter {

	HashMap<String, HashMap<Integer,Snp>> snps;
	HashMap<String, HashMap<String,Vector<int[]>>> deletions;   // HashMap<scaffold, HashMap<accession,Vector<positions>>>
	
	public Hunter(File wildtypePileup, HashMap<String, File > mutantPileups)throws IOException{
		this.snps= new HashMap<String, HashMap<Integer,Snp>>();
		this.deletions = new HashMap<String, HashMap<String, Vector<int[]>>>();
	}
	
	
	
	
	/**
	 * 
	 *   read simultanously a set of pileup files and record snps.
	 * 
	 *  Still to do: call deletions
	 *  
	 * @param pileupFiles
	 * @param fastaIndex
	 * @param minCoverage
	 * @param minAlleleFrequency
	 * @throws IOException
	 */
	public void readPileups(HashMap<String,File> pileupFiles, File fastaIndex, int minCoverageForSNP, double minAlleleFrequency, int coverageThresholdForDeletion)throws IOException{
		
		int counta = 0;
		int numSNPs = 0;
		
		BufferedReader inFAIDX = new BufferedReader(new FileReader(fastaIndex));

		HashMap<String, String> currentPositions = new HashMap<String, String>(); // store the current input lines of pileups
		HashSet<String> visitedScaffolds = new HashSet<String>(); //store, which scaffolds have been visited
		
		
		
		//open BufferedReaders for every pileup file.
		HashMap<String,BufferedReader> ins = new HashMap<String,BufferedReader> ();
		for(Iterator<String> iterator = pileupFiles.keySet().iterator(); iterator.hasNext();) {
			String name = iterator.next();
			File file = pileupFiles.get(name);
			BufferedReader in;
			FileInputStream fis = new FileInputStream(file);
			byte[] bytes = new byte[2];
			fis.read(bytes);
			int head = ((int) bytes[0] & 0xff) | ((bytes[1] << 8) & 0xff00);
			boolean gzip = GZIPInputStream.GZIP_MAGIC == head;
			fis.close();
			if(gzip){
				in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
			}else{
				in = new BufferedReader((new FileReader(file)));
			}
			
			ins.put(name,in); //add the buffered reader
			currentPositions.put(name, in.readLine()); //add the first line to currentpositions
			
		}
		
		
		
		
		for (String inputline = inFAIDX.readLine(); inputline != null; inputline = inFAIDX.readLine()) {
			String[] faidxSplit = inputline.split("\t");
			String scaffold = faidxSplit[0];
			int length = Integer.parseInt(faidxSplit[1]);
			visitedScaffolds.add(scaffold);
			
			
			HashMap<String, int[]> coverages = new HashMap<String, int[]>();
			for( Iterator<String> iterator = ins.keySet().iterator(); iterator.hasNext();) {
				String accession = iterator.next();
				coverages.put(accession, new int[length]);
			}
			
					
			for(int position = 1; position<= length; position++)	{
				System.out.println(scaffold + "\t" + position);
				counta++;
				if(counta%100000==0) {
					//System.out.print(".");
					if(counta%10000000==0) {
					//	System.out.println(" " + counta + " at " + scaffold +". Found " + numSNPs + " SNPs." );
					}
				}
				
				
				//check if one pileup has a snp
				
				HashMap<String, MPileupLine> currentMPileup = new HashMap<String, MPileupLine>();
				
				for(Iterator<String> iterator = currentPositions.keySet().iterator(); iterator.hasNext() ;) {
					String accession = iterator.next();
					
					MPileupLine line = new MPileupLine(currentPositions.get(accession));
					String myscaffold = line.getChromosome();
					
					
					
					if( !visitedScaffolds.contains(myscaffold)) { //wait
						continue;
					}
					int pos = line.getPosition();
					if( scaffold.equalsIgnoreCase(myscaffold) &&  pos > position) { //wait
						continue;
					}
					
					
					while( visitedScaffolds.contains(myscaffold) && pos < position) { //sync
						String currentInputLine =  ins.get(accession).readLine();
						currentPositions.put(accession, currentInputLine);
						line = new MPileupLine(currentInputLine);
						myscaffold = line.getChromosome();
						pos = line.getPosition();
					}
					
					if( myscaffold.equalsIgnoreCase(scaffold) && pos == position) { //here we are. In sync and ready to investigate.
						System.out.println("\t" + accession + " sync");
						currentMPileup.put(accession, line);
					}
					
				}	
				
				//check for SNP
				for(Iterator<String> iterator = currentMPileup.keySet().iterator(); iterator.hasNext();) {
					String accession = iterator.next();
					MPileupLine line = currentMPileup.get(accession);
					int cov = line.getCoverage();
					if( cov < minCoverageForSNP) {
						continue;
					}
					double freq = line.getLargestAlternativeAlleleFrequency();
					if( freq >= minAlleleFrequency) { //getSNP!
						Snp snp = new Snp(currentMPileup, scaffold, position);
						
						if(!this.snps.containsKey(scaffold)) {
							this.snps.put(scaffold, new HashMap<Integer, Snp>());
						}
						this.snps.get(scaffold).put(position, snp);
						System.out.println("found SNP. " + scaffold + ": " + position + ", " + currentMPileup.size() + " accessions\t" + snp.getAccessions());
						break;
					}
				}
				
				
			}	
		}		
				
					

		inFAIDX.close();
		
		
		for(Iterator<String> iterator = ins.keySet().iterator(); iterator.hasNext();) {
			String accession = iterator.next();
			ins.get(accession).close();
		}
		
		
	}
	
	/**
	 * 
	 * a helper method to extract deletions in form of intervals from an int[] that holds coverage.
	 * 
	 * @param coverage
	 * @param coverageThreshold
	 * @return
	 */
	private Vector<int[]> callDeletions(int[] coverage, int coverageThreshold){
		Vector<int[]> deletions = new Vector<int[]>();
		
		int start = -1;
		for( int i = 0; i< coverage.length; i++) {
			if (coverage[i] < coverageThreshold) {
				if( start ==-1) {
					start = i;
				}
			}else {
				if( start >-1) {
					int end = i-1;
					if(end-start>0) {
						int[] deletion = {start+1, end+1};
						deletions.add(deletion);
						start = -1;
					}
				}
			}
			
			
		}
		
		return deletions;
	}
	
	
	
	/**
	 * 
	 * write a TSV table with every SNP. 
	 * 
	 * First column: Scaffold
	 * Second column: position
	 * third column: reference allele
	 * fourth column: alternative allele
	 * 
	 * every additional collumn is accession,alternative_allele_frequency_coverage
	 * 
	 * 
	 * @param outputFile
	 * @throws IOException
	 */
	public void writeSNPs(File outputFile)throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		
		for(Iterator<String> iterator = this.snps.keySet().iterator(); iterator.hasNext();) {
			String scaffold = iterator.next();
			for(Iterator<Integer> iterator2 = this.snps.get(scaffold).keySet().iterator(); iterator2.hasNext();) {
				int pos = iterator2.next();
				Snp snp = this.snps.get(scaffold).get(pos);
				out.write(scaffold + "\t" + pos + "\t" + snp.toString());
				out.newLine();
			}
		}
		
		out.close();
	}
	
	
	/**
	 * 
	 * load snps from a TSV table
	 * 
	 * First column: Scaffold
	 * Second column: position
	 * third column: reference allele
	 * fourth column: alternative allele
	 * 
	 * every additional collumn is accession,alternative_allele_frequency_coverage
	 
	 * @param inputFile
	 * @throws IOException
	 */
	public void loadSNPs(File inputFile)throws IOException{
		BufferedReader in = new BufferedReader(new FileReader(inputFile));

		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String[] split = inputline.split("\t");
			Snp snp = new Snp(split);
			
			int pos = Integer.parseInt(split[1]);
			if(!this.snps.containsKey(split[0])) {
				this.snps.put(split[0], new HashMap<Integer, Snp> ());
			}
			this.snps.get(split[0]).put(pos, snp);
		}

		in.close();
		
	}
	
	
	/**
	 * 
	 * write deletions to a TSV file. One row per scaffold and accession. First column is scaffold, second column is accession. Every following column is start and end of a deletion.
	 * 
	 * @param outputFile
	 * @throws IOException
	 */
	public void writeDeletions(File outputFile)throws IOException{
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		
		for(Iterator<String> iterator = this.deletions.keySet().iterator(); iterator.hasNext();) {
			String scaffold = iterator.next();
			for(Iterator<String> iterator2 = this.deletions.get(scaffold).keySet().iterator(); iterator2.hasNext();) {
				String accession = iterator2.next();
				out.write(scaffold + "\t" + accession);
				for(Iterator<int[]> iterator3 = this.deletions.get(scaffold).get(accession).iterator(); iterator3.hasNext();) {
					int[] deletion = iterator3.next();
					out.write("\t" + deletion[0] + "," + deletion[1]);
				}
				out.newLine();
			}
		}
		out.close();
	}
	
	
	/**
	 * 
	 * load deletions from TSV file. One row per scaffold and accession. First column is scaffold, second column is accession. Every following column is start and end of a deletion.
	 * 
	 * @param inputFile
	 * @throws IOException
	 */
	public void loadDeletions(File inputFile)throws IOException{
		BufferedReader in = new BufferedReader(new FileReader(inputFile));

		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String[] split = inputline.split("\t");
			if(!this.deletions.containsKey(split[0])) {
				this.deletions.put(split[0],new HashMap<String, Vector<int[]>>());
			}
			Vector<int[]> v = new Vector<int[]>();
			for( int i = 2; i< split.length; i++) {
				String[] split2 = split[i].split(",");
				int[] deletion = {Integer.parseInt(split2[0]), Integer.parseInt(split2[1])};
				v.add(deletion);
			}
			this.deletions.get(split[0]).put(split[1], v);
		}

		in.close();
		
	}
}
