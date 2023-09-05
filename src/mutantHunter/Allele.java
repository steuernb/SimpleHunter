package mutantHunter;

public class Allele implements Comparable<Allele>{

	int id;
	int variation_id;
	String nucleotideSequence;
	boolean isRef;
	
	public Allele(int id, int variation_id, String nucleotideSequence, boolean isRef) {
		super();
		this.id = id;
		this.variation_id = variation_id;
		this.nucleotideSequence = nucleotideSequence;
		this.isRef = isRef;
	}
	

	public Allele(int id, int variation_id, String nucleotideSequence) {
		super();
		this.id = id;
		this.variation_id = variation_id;
		this.nucleotideSequence = nucleotideSequence;
		this.isRef = false;
	}


	public boolean isRef() {
		return isRef;
	}


	public void setRef(boolean isRef) {
		this.isRef = isRef;
	}


	public int getId() {
		return id;
	}


	public int getVariation_id() {
		return variation_id;
	}


	public String getNucleotideSequence() {
		return this.nucleotideSequence;
	}
	
	
	public int compareTo(Allele allele) {
		return this.nucleotideSequence.compareTo(allele.getNucleotideSequence());
	}
	
	public boolean isSameAllele(Allele allele){
		
		if( this.nucleotideSequence.equalsIgnoreCase(allele.getNucleotideSequence())) {
			return true;
		}
		
		return false;
	}
	
	
	
	
	/**
	 * 
	 * Check if an allele is the same as this one. Heterozygosity is generously regarded (see below).
	 * 
	 * The current encoding of a "heterozygous allele" is for example HET(A/T)
	 * If a het allele is compared to any allele with this method any match will return true: 
	 * A to HET(A/T) is true
	 * HET(A/C) to HET(A/T) is true
	 * C to HET(A/T) is false
	 * 
	 * 
	 * @param allele
	 * 			The allele this one is compared to.
	 * @return
	 */
	public boolean isSameAlleleIgnoreHet(Allele allele) {
		
		String[] myalleles1 = {this.nucleotideSequence};
		String[] myalleles2 = {allele.getNucleotideSequence()};
		
		if( this.nucleotideSequence.startsWith("HET")) {
			 
			myalleles1 = this.nucleotideSequence.substring(4, this.nucleotideSequence.length()-1).split("/");
		
		}
		
		
		if( allele.nucleotideSequence.startsWith("HET")) {
		
			myalleles2 = allele.getNucleotideSequence().substring(4, allele.getNucleotideSequence().length()-1).split("/");
		
		}
		
		for( int i = 0; i< myalleles1.length; i++) {
			String a1 = myalleles1[i];
			for( int j = 0; j< myalleles2.length; j++) {
				String a2 = myalleles2[j];
				if(a1.equalsIgnoreCase(a2)) {
					return true;
				}
			}
		}
		
		
		
		return false;
	}
	
}
