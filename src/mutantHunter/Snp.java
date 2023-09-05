package mutantHunter;

import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

import javax.sound.sampled.Line;

import support.MPileupLine;


/**
 * 
 * @author steuernb
 * @deprecated
 *
 */
public class Snp {

	char referenceAllele;
	char alternativeAllele;
	
	public final char[] alleleIndex = {'A', 'T','G','C'};
	
	HashMap<String, int[]> alleles;
	
	

	public Snp(HashMap<String, MPileupLine> h, String scaffold, int position) {
		
		this.alleles = new HashMap<String, int[]>();
		
		
		double largestAlternativeAlleleFrequency = 0;
		char largestAlteranativeAllele = ' ';
		
		for(Iterator<String> iterator = h.keySet().iterator(); iterator.hasNext();) {
			String accession = iterator.next();
			MPileupLine line = h.get(accession);
			
			if( !line.getChromosome().equalsIgnoreCase(scaffold) || line.getPosition() != position) {
				System.out.println("SNP " + scaffold + " " + position + ":out of sync!");
			}
			this.referenceAllele = line.getRefBase();
			int[] alleleCounts = new int[4];
			int coverage = 0;
			for(Iterator<String> iterator2 = line.getAlleles().keySet().iterator(); iterator2.hasNext();) {
				String allele = iterator2.next();
				int count = line.getAlleles().get(allele);
				if( allele.equalsIgnoreCase("A")) {alleleCounts[0] = count;coverage = coverage + count;}
				if( allele.equalsIgnoreCase("T")) {alleleCounts[1] = count;coverage = coverage + count;}
				if( allele.equalsIgnoreCase("G")) {alleleCounts[2] = count;coverage = coverage + count;}
				if( allele.equalsIgnoreCase("C")) {alleleCounts[3] = count;coverage = coverage + count;}
				
			}
			
			for( int i = 0; i< alleleCounts.length; i++) {
				double freq = (double)alleleCounts[i]/(double)coverage;
				if( freq >  largestAlternativeAlleleFrequency && alleleIndex[i] != this.referenceAllele) {
					largestAlternativeAlleleFrequency = freq;
					largestAlteranativeAllele = alleleIndex[i];
					
				}
				
			}
			
			this.alleles.put(accession, alleleCounts);
			
		}
		
		this.alternativeAllele = largestAlteranativeAllele;
	}
	
	public Snp(String[] split) {
		this.alleles = new HashMap<String, int[]>();
		
		this.setReferenceAllele(split[2].toCharArray()[0]);
		this.setAlternativeAllele(split[3].toCharArray()[0]);
		
		for( int i = 4; i< split.length; i++) {
			String[] alleleString = split[i].split(",");
			
			int[] allele = {Integer.parseInt(alleleString[1]),Integer.parseInt(alleleString[2]),Integer.parseInt(alleleString[3]),Integer.parseInt(alleleString[4]) };
			this.alleles.put(alleleString[0],allele);
			
		}
	}
	
	public String toString() {
		String s  = referenceAllele + "\t" + alternativeAllele;
		
		for(Iterator<String> iterator = this.alleles.keySet().iterator(); iterator.hasNext();) {
			String accession = iterator.next();
			int[] allele = this.alleles.get(accession);
			s = s + "\t" + accession+","+allele[0] +","+allele[1] +","+allele[2] +","+allele[3];
		}
		
		return s;
	}
	
	public char getReferenceAllele() {
		return referenceAllele;
	}

	

	public void setReferenceAllele(char referenceAllele) {
		this.referenceAllele = referenceAllele;
	}


	public char getAlternativeAllele() {
		return alternativeAllele;
	}


	public void setAlternativeAllele(char alternativeAllele) {
		this.alternativeAllele = alternativeAllele;
	}


	
	public double getAlleleFrequency(String accession) {
		int coverage = 0;
		int alleleCount = 0;
		int[] allele = this.alleles.get(accession);
		for( int i = 0; i< allele.length; i++) {
			coverage = coverage + allele[i];
			if( alleleIndex[i] == this.alternativeAllele) {
				alleleCount = allele[i];
			}
		}
		return (double) alleleCount /  (double) coverage;
	}
	
	
	
	
	public int getCoverage(String accession) {
		int[] allele = this.alleles.get(accession);
		return allele[0] + allele[1] + allele[2] + allele[3];
	}
	
	
	
	public boolean isTransition() {
		
		if( referenceAllele == 'G' && alternativeAllele == 'A') {
			return true;
		}
		if( referenceAllele == 'C' && alternativeAllele == 'T') {
			return true;
		}
		return false;
		
	}
	
	public String getAccessions() {
		String s = "";
		Vector<String> v= new Vector<String>();
		v.addAll(alleles.keySet());
		Collections.sort(v);
		for(Iterator<String> iterator = v.iterator(); iterator.hasNext();) {
			s = s + "," + iterator.next();
		}
		return s.substring(1);
	}
	
}
