package simpleHunter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Vector;

import support.commandLineInterface.CLI;

public class Hunter {

	Hashtable<String, Hashtable<String,Hashtable<Integer, SNP>>> mainTable;
	
	
	
	public Hunter(){
		mainTable = new Hashtable<String, Hashtable<String,Hashtable<Integer, SNP>>>();
	}
	
	


	
	
	
	public void readMutantFiles(Vector<File> snpFiles)throws IOException{
		
		for(Enumeration<File> enum_inputFiles = snpFiles.elements(); enum_inputFiles.hasMoreElements();){
			File file = enum_inputFiles.nextElement();
			
			String fileName = file.getName();
			BufferedReader in = new BufferedReader(new FileReader(file));
			String currentContig = "";
			Hashtable<Integer, SNP> h = new Hashtable<Integer, SNP>();
			for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
				if(inputline.startsWith("#")){
					continue;
				}
				String[] split = inputline.split("\t"); //contig	position	referenceAllele	alternativeAllele	refAlleleFreq largestAltAlleleFreq	Coverage
				
				if(!currentContig.equalsIgnoreCase(split[0])){
					if(!this.mainTable.containsKey(currentContig)){
						this.mainTable.put(currentContig, new Hashtable<String,Hashtable<Integer, SNP>>());
					}
					
					this.mainTable.get(currentContig).put(fileName, h);
					h = new Hashtable<Integer, SNP>();
					currentContig = split[0];
				}
				
				if( 	(split[2].equalsIgnoreCase("A") || split[2].equalsIgnoreCase("T") || split[2].equalsIgnoreCase("G") || split[2].equalsIgnoreCase("C") ) &&
						(split[3].equalsIgnoreCase("A") || split[3].equalsIgnoreCase("T") || split[3].equalsIgnoreCase("G") || split[3].equalsIgnoreCase("C") )){
					
					SNP snp = new SNP(fileName, currentContig, split[2], split[3], Integer.parseInt(split[1]), Integer.parseInt(split[6]), Double.parseDouble(split[4]), Double.parseDouble(split[5]));
					h.put(Integer.parseInt(split[1]), snp);
				}	
			}

			in.close();
		}
		
	}
	
	
	public String hunt(int windowSize, int minMutations, int minEMSMutations, int maxMutantsSharePosition, double maxRefAlleleFrequency, HashSet<String> wtFileNames,HashMap<String,Vector<int[]>> exons){
		
		String returnString = "";
		
		for(Enumeration<String> enumContig = mainTable.keys(); enumContig.hasMoreElements();){
			String contig = enumContig.nextElement();
		
			
			
			//mark unusable mutants
			HashSet<Integer> remove = this.filterForMultipleMutantLinesAtPosition(contig, maxMutantsSharePosition);
			
			
			
			
			
			
			//create a data structure that can quickly be searched with sliding window.
			Hashtable<String, Vector<Integer>> positions = new Hashtable<String,Vector<Integer>>();
			int maxPosition = 0; //need that for later
			for(Enumeration<String> enumFileNames = mainTable.get(contig).keys(); enumFileNames.hasMoreElements(); ){
				String fileName = enumFileNames.nextElement();
				if( wtFileNames.contains(fileName)){
					continue;
				}
				
				Vector<Integer> v = new Vector<Integer>();
				for(Enumeration<Integer> enumPositions = mainTable.get(contig).get(fileName).keys(); enumPositions.hasMoreElements();){
					Integer position = enumPositions.nextElement();
					if( !remove.contains(position)  &&  mainTable.get(contig).get(fileName).get(position).getReferenceAlleleFrequency() <=maxRefAlleleFrequency ){
						v.add(position);
						if(position>maxPosition){
							maxPosition = position;
						}
					}
					
				}
				Collections.sort(v);
				if(v.size()>0){	
					positions.put(fileName, v);
				}
			}
			
			
			
			if(positions.size() < minMutations){
				continue;
			}
			
			
		
			
			
			
			
			Hashtable<String, Integer> indices = new Hashtable<String, Integer>();
			for(Enumeration<String> enumFileNames = positions.keys(); enumFileNames.hasMoreElements();){
				String fileName = enumFileNames.nextElement();
				indices.put(fileName, 0);
			}
			
			while(true){
				int minPosition = maxPosition+1; // only to initialize;
				String fileNameWithLowestPosition = "";
				for(Enumeration<String> enumFileNames = positions.keys(); enumFileNames.hasMoreElements();){
					String fileName = enumFileNames.nextElement();
					if(positions.get(fileName).size() > indices.get(fileName)){
						int position = positions.get(fileName).get(indices.get(fileName));
						if(position < minPosition){
							minPosition = position;
							fileNameWithLowestPosition = fileName;
						}
					}
				}
				
				if(fileNameWithLowestPosition.equalsIgnoreCase("")){ //stop condition
					break;
				}
				
				int numSNPs = 0; 
				int numEMSSNPs = 0;
				
				String report = contig + "\t" + contig + ":"+ Math.max(0,minPosition-1000) +"-" + (minPosition+ windowSize + 1000) ;
				Vector<String> myFileNames = new Vector<String>();
				for(Enumeration<String> enumFileNames = positions.keys(); enumFileNames.hasMoreElements();){
					myFileNames.add(enumFileNames.nextElement());
				}
				Collections.sort(myFileNames);
				for(Enumeration<String> enumFileNames = myFileNames.elements(); enumFileNames.hasMoreElements();){
					String fileName = enumFileNames.nextElement();
					
					int index = indices.get(fileName).intValue();
					if(positions.get(fileName).size()>index){
						int position = positions.get(fileName).get(indices.get(fileName));
						if( position >= minPosition && position <= minPosition + windowSize){
							numSNPs++;
							
							String coding = "";
							if(this.isInExon(exons, contig, position)) {
								coding = "-coding-";
							}
							
							SNP snp = mainTable.get(contig).get(fileName).get(position);
							report = report + "\t" + fileName +"(" +  position+")"+coding+":"+snp.getReferenceAllele()+"->"+snp.getLargestAlternativeAllele();
							
							if(snp.isEMSMutation()){
								numEMSSNPs++;
							}
						}
					}
					
				}
				
				
				if(numSNPs >= minMutations && numEMSSNPs >= minEMSMutations ){
					returnString = returnString + "\n" + report;
					
				}
				
				
				int index = indices.get(fileNameWithLowestPosition).intValue();
				index ++;
				indices.put(fileNameWithLowestPosition, index);
				
			}
			
			
			
			
			
			
			
		}
		return returnString;
	}
	
	private boolean isInExon(HashMap<String,Vector<int[]>> exons, String chr, int pos) {
	
		for(Iterator<int[]> iterator= exons.get(chr).iterator(); iterator.hasNext();) {
			int[] a = iterator.next();
			if(a[0] <=pos && a[1] >=pos) {
				return true;
			}
			
		}
		
		return false;
		
	}
	
	
	private HashSet<Integer> filterForMultipleMutantLinesAtPosition(String contig, int maxMutantsSharePosition){
		//mark unusable mutants
		Hashtable<Integer, Integer> wasThere = new Hashtable<Integer, Integer>();
		HashSet<Integer> remove = new HashSet<Integer>();
		for(Enumeration<String> enumFileNames = mainTable.get(contig).keys(); enumFileNames.hasMoreElements();){
			String fileName = enumFileNames.nextElement();
			for(Enumeration<Integer> enumPositions = mainTable.get(contig).get(fileName).keys();enumPositions.hasMoreElements();){
				Integer position = enumPositions.nextElement();
				int num = 0;
				if(wasThere.containsKey(position)){
					num = wasThere.get(position).intValue();
				}
				num++;
				wasThere.put(position, num);
			}
		}
		
		for(Enumeration<Integer> enumPositions = wasThere.keys(); enumPositions.hasMoreElements();){
			Integer position = enumPositions.nextElement();
			if( wasThere.get(position).intValue() > maxMutantsSharePosition){
				remove.add(position);
			}else{
				
			}
		}
		return remove;
	}
	
	
	
	
	public HashMap<String,Vector<int[]>> loadExons(File inputGFF)throws IOException{
		
		HashMap<String,Vector<int[]>> h = new HashMap<String,Vector<int[]>>();
		
		
		BufferedReader in = new BufferedReader(new FileReader(inputGFF));

		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {

			if(inputline.startsWith("#")) {
				continue;
			}
			String[] split = inputline.split("\t");
			if(split[2].equalsIgnoreCase("exon")) { 
				
				String chr  = split[0];
				int[] pos = {Integer.parseInt(split[3]) , Integer.parseInt(split[4])};
				if(!h.containsKey(chr)) {
					h.put(chr, new Vector<int[]>());
				}
				h.get(chr).add(pos);
			}
			
		}

		in.close();
		
		return h;
		
	}
	
	
	
	
	public static void main(String[] args) {
		
		
		
		
		
		String help = "java -jar SimpleHunter.jar -w wt.snps.txt -i m1.snps.txt [mn.snps.txt]+ -o outputFile -l window_length -m minNumberOfMutants -a maxRefAlleleFrequency -n maxNumberMutantsSharingPosition -g exons.gff";
		
				
		try {
			CLI cli = new CLI();
			cli.parseOptions(args);
			
			if(!cli.hasOption("i")) {
				System.out.println("Please provide mutant files");
				System.out.println("\n"+help);
				
			}else if(!cli.hasOption("w")) {
				System.out.println("Please provide wt file");
				System.out.println("\n"+help);
			}else {
				
				File wtFile = new File(cli.getArg("w"));
				Vector<String> mutantFileNames = cli.getArgs("i");
				
				Vector<File> inputFiles = new Vector<File>();
				inputFiles.add(wtFile);
				for(Iterator<String> iterator = mutantFileNames.iterator(); iterator.hasNext();) {
					inputFiles.add(new File(iterator.next()));
				}
		
				
				HashSet<String> wtFileNames = new HashSet<String>();
				wtFileNames.add(new File(cli.getArg("w")).getName());
				
				int windowSize = 5000;
				int minMutations = 2;
				int minEMSMutations = 0;
				double maxRefAlleleFrequency = 0.7;
				int maxMutantsSharePosition = 1;
				
				if( cli.hasOption("l")) {
					windowSize = Integer.parseInt(cli.getArg("l"));
				}
				if( cli.hasOption("m")) {
					minMutations = Integer.parseInt(cli.getArg("m"));
				}
				
				if( cli.hasOption("a")) {
					maxRefAlleleFrequency = Double.parseDouble(cli.getArg("a"));
				}
				if( cli.hasOption("n")) {
					maxMutantsSharePosition = Integer.parseInt(cli.getArg("n"));
				}
				
				
				
				Hunter hunter = new Hunter();	
				
				HashMap<String,Vector<int[]>> exons = new HashMap<String,Vector<int[]>>();
				
				if(cli.hasOption("g")) {
					exons = hunter.loadExons(new File(cli.getArg("g")));
				}
			
			
				
				hunter.readMutantFiles(inputFiles);
	
				String s = hunter.hunt(windowSize, minMutations, minEMSMutations, maxMutantsSharePosition, maxRefAlleleFrequency, wtFileNames, exons);	
				
				
				if( cli.hasOption("o")) {
					BufferedWriter out = new BufferedWriter(new FileWriter(cli.getArg("o")));
					out.write(s);
					out.close();
				}else {
					System.out.println(s);
				}
				
				
			}
			
			
			
			
			
			
			
				
				
				
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println(help);
		}		
		
		
		
		
		
		
		
		
	}
	
	
	
	
}
