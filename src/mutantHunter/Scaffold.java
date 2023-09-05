package mutantHunter;

import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;
/**
 * @deprecated
 * @author steuernb
 *
 */
public class Scaffold {

	String id;
	
	HashMap<String, HashMap<Integer, Snp> >  snps;   // <accession, <position<SNP>>
	
	HashMap<String, Vector<int[]>> deletions;
	
	static final String wildtype = "wildtype";
	
	public Scaffold(String id) {
		
		this.id = id;
		
		snps = new HashMap<String, HashMap<Integer, Snp> >();
		deletions = new HashMap<String, Vector<int[]>>();
		
	}
	
	
	public void addVariation(String line, int position, Snp var) {
		
		if( !snps.containsKey(line)) {
			snps.put(line, new HashMap<Integer, Snp>());
		}
		snps.get(line).put(position, var);
	}
	
	public void addDeletion(String line, int[] deletion) {
		
		if( !deletions.containsKey(line)) {
			deletions.put(line, new Vector<int[]>());
		}
		deletions.get(line).addElement(deletion);
		
	}
	
	/**
	 * 
	 * remove every position that is shared by more than the allowed number of mutations that may share a SNP.
	 * 
	 * Note: the wildtype is not counted into the allowed number of mutants here. But positions that are removed are also removed from wildtype/
	 * 
	 * @param maxNumberOfMutantsSharingASnp
	 * 			The maximum number of mutants that may share a position having a SNP.
	 */
	public void filterSNPsForNumberOfMutants(int maxNumberOfMutantsSharingASnp) {
		HashSet<Integer> positions = new HashSet<Integer>();
		for(Iterator<String> iterator = this.snps.keySet().iterator(); iterator.hasNext();) {
			String accession = iterator.next();
			if(!accession.equalsIgnoreCase(wildtype)) {
				for( Iterator<Integer> iterator2 = snps.get(accession).keySet().iterator(); iterator2.hasNext();) {
					int pos = iterator2.next();
					positions.add(pos);
				}
				
			}
			
			
		}
		
		for(Iterator<Integer> iterator = positions.iterator(); iterator.hasNext();){
			int pos = iterator.next();
			int num = 0;
			for( Iterator<String> iterator2 = snps.keySet().iterator(); iterator2.hasNext();) {
				String accession = iterator2.next();
				if(!accession.equalsIgnoreCase(wildtype)) {
					if(snps.get(accession).containsKey(pos) ) {
						num++;
					}
				}
			}
			
			if( num > maxNumberOfMutantsSharingASnp ) {
				for( Iterator<String> iterator2 = snps.keySet().iterator(); iterator2.hasNext();) {
					String accession = iterator2.next();
					snps.get(accession).remove(pos);
				}
			}
			
			
		}
	
	}
	
	
	
	
	public void filterSNPsForWildtype() {
		for(Iterator<Integer> iterator = this.snps.get(wildtype).keySet().iterator(); iterator.hasNext();) {
			int pos = iterator.next();
			for(Iterator<String> iterator2 = this.snps.keySet().iterator(); iterator2.hasNext();) {
				String accession = iterator2.next();
				if( !accession.equalsIgnoreCase(wildtype)) {
					this.snps.get(accession).remove(pos);
				}
				
			}
		}
	}
	
	/*
	public void getCandidates(int intervalLength, int minMutants ) {
		
		Vector<String[]> positions = new Vector<String[]>();
		
		for(Iterator<String> iterator1 = snps.keySet().iterator(); iterator1.hasNext();) {
			String line = iterator1.next();	
			
			if( line.equalsIgnoreCase(wildtype)) {
				continue;
			}
			
			for( Iterator<Integer> iterator2 = snps.get(line).keySet().iterator(); iterator2.hasNext();) {
				int position = iterator2.next();
				Snp snp = snps.get(line).get(position);
				String[] a = { position + "", line , snp.getReferenceAllele() + "," + snp.getAlternativeAllele() + "," + snp.getCoverage()   };
				positions.add(a);
			}
			
			
		}
		
		for( Iterator<String> iterator1 = deletions.keySet().iterator(); iterator1.hasNext();) {
			String line = iterator1.next();
			if( line.equalsIgnoreCase(wildtype)) {
				continue;
			}
			for(Iterator<int[]> iterator2 = deletions.get(line).iterator(); iterator2.hasNext();) {
				int[] deletion_interval = iterator2.next();
				String[] a = { deletion_interval[0]+"", line, "deletion", deletion_interval[1] +"" };
				positions.add(a);
				
				for( int i = deletion_interval[0]; i<= deletion_interval[1]; i++ ) {
					
				}
			}
		}
		
		Collections.sort(positions, new Comparator<String[]>() {
			public int compare(String[] s1, String[] s2) {
				int i1 = Integer.parseInt( s1[0] );
				int i2 = Integer.parseInt( s2[0] );
				
				if( i1 < i2) {return -1;}
				if( i2 < i1) {return 1;}
				return 0;
			}
		});
		
		Vector<String[]> myInterval =  new Vector<String[]>();
		myInterval.add(positions.firstElement());
		int startPosition = Integer.parseInt(myInterval.firstElement()[0]); 
		int[] currentInterval = {startPosition, startPosition + intervalLength};
		
		for( Iterator<String[]> iterator = positions.iterator(); iterator.hasNext();) {
			
			String[] a = iterator.next();
			int pos = Integer.parseInt(a[0]);
			if( pos <= currentInterval[1]) {
				
				myInterval.add(a);
				
			}else {
				
				//evaluate interval
				String s = "Interval " + this.id + ":" + currentInterval[0] + "-" + currentInterval[1];
				HashSet<String> mutants = new HashSet<String>();//figure out how many mutants are mutated in this interval
				for(Iterator<String[]> iteratorInterval = myInterval.iterator(); iteratorInterval.hasNext(); ) {
					String[] variation = iteratorInterval.next();
					mutants.add(variation[1]);
					if( variation[2].equalsIgnoreCase("deletion")) {
						s = s + "\n\t" + variation[1] + "\tdeletion:" + variation[0] + "-" + variation[3];  
					}else {
						s = s + "\n\t" + variation[1] + "\tSNP:" + variation[0] + ";" + variation[2];
					}
				}
				
				if( mutants.size() >= minMutants) {//enough mutants mutated in this interval?
					
					System.out.println(s);
				}
				
				//remove elements downstream of pos - intervalLength
				int newIntervalStart = pos;
				
				for( int i = myInterval.size()-1; i>=0;i--) {
					String[] variation = myInterval.get(i);
					int variationPosition = 0;
					if( variation[2].equalsIgnoreCase("deletion")) {
						variationPosition = Integer.parseInt(variation[3]);
					}else {
						variationPosition = Integer.parseInt(variation[0]);
					}
					if(variationPosition < pos - intervalLength) {
						myInterval.remove(i);
					}else {
						if( newIntervalStart > variationPosition) {
							newIntervalStart = variationPosition;
						}
					}
				}
				currentInterval[0] = newIntervalStart;
				currentInterval[1] = newIntervalStart + intervalLength;
				
			}
			
		}
		
		
		
	}
	
	
	*/
	
	
}
