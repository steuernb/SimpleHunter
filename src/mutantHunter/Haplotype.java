package mutantHunter;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

/**
 * 
 * A haplotype is a specific set of alleles for a given set of variations
 * 
 * @author steuernb
 *
 */
public class Haplotype implements Comparable<Haplotype> {

	int id;
	
	HashMap<Integer, Integer> alleles; //<variation_id, allele_id>
	
	Vector<Integer> accessions;
	
	
	public Haplotype(int id) {
		this.id = id;
		this.alleles = new HashMap<Integer, Integer>();
		this.accessions = new Vector<Integer>();
	}
	
	
	
	/**
	 * 
	 * this returns 0 only if both haplotypes contain the same variations and the allele for each variation is the same.
	 * Otherwise,
	 * 			1. The haplotype with more variations comes first
	 * 			2. This haplotype comes first if the other does not have one of the variations this one has (yep, that's random)
	 * 			3. The haplotype with more accessions
	 * 			4. The other haplotype (in case of same number of accessions)
	 * 
	 * 
	 */
	public int compareTo(Haplotype haplotype) {
		
		if( this.alleles.size() > haplotype.getAlleles().size()) {
			//System.out.println("number of alleles does not match");
			return -1;
		}
		if( this.alleles.size() < haplotype.getAlleles().size()) {
			//System.out.println("number of alleles does not match");
			return 1;
		}
		//String s1 = "";
		//String s2 = "";
		
		for( Iterator<Integer> iterator = this.alleles.keySet().iterator(); iterator.hasNext();) {
			
			int variation_id = iterator.next();
			//s1 = s1 + ":" + this.alleles.get(variation_id);
			//s2 = s2 + ":" + haplotype.getAlleles().get(variation_id);
			
			if( haplotype.getAlleles().containsKey(variation_id)) {
				if( this.alleles.get(variation_id).intValue() !=  haplotype.getAlleles().get(variation_id).intValue()) {
					//System.out.println("allele id does not match: " + this.alleles.get(variation_id) +" " + haplotype.getAlleles().get(variation_id));
					if(this.getAccessions().size() > haplotype.getAccessions().size()) {
						return -1;
					}else {
						return 1;
					}
					
				}
				
			}else {
				//System.out.println(" variation id is missing in other haplotype");
				return -1;
			}
			
			
			
		}
		/*
		System.out.println(s1.substring(1));
		System.out.println(s2.substring(1));
		System.out.println("equal");
		System.out.println();
		*/
		
		return 0;
	}


	public HashMap<Integer, Integer> getAlleles() {
		return alleles;
	}


	public int getAllele_id(int variation_id) {
		return this.alleles.get(variation_id);
	}
	
	
	
	/**
	 * 
	 * Set the allele for a specific variation.
	 * 
	 * @param variation_id
	 * 			The id for the variation (A Haplotype object only holds the IDs of variations. The actual Variation objects are held by an HaplotypeCaller object).
	 * @param allele_id
	 * 			The id for the allele (A Haplotype object only holds the IDs of alleles. The actual Allele objects are held by a Variation object).
	 */
	public void setAllele_id(int variation_id, int allele_id) {
		this.alleles.put(variation_id, allele_id);
	}
	
	
	public Vector<Integer> getAccessions() {
		return accessions;
	}
	
	public void addAccession(int accession_id){
		this.accessions.add(accession_id);
	}
	
	public void addAccesions(Vector<Integer> accession_ids) {
		this.accessions.addAll(accession_ids);
	}
	
	/**
	 * 
	 * Add a variation with a specific allele to this haplotype
	 * 
	 * @param allele
	 */
	public void addAllele(Allele allele) {
		if( 	this.alleles.containsKey(allele.getVariation_id()) &&
				this.alleles.get(allele.getVariation_id()) != null &&	
				this.alleles.get(allele.getVariation_id()).intValue() != allele.getId()) {
			
			System.out.println("DEBUG: This haplotype already has an allele for variation " + allele.getVariation_id());
		
		}
		this.alleles.put(allele.getVariation_id(), allele.getId());
	}
	
	
	/**
	 * Add a variation with a specific allele to this haplotype
	 * 
	 * @param variation_id
	 * @param allele_id
	 */
	public void addAllele(int variation_id, int allele_id) {
		if( 	this.alleles.containsKey(variation_id) &&
				this.alleles.get(variation_id) != null &&	
				this.alleles.get(variation_id).intValue() != allele_id) {
			
			System.out.println("DEBUG: This haplotype already has an allele for variation " + variation_id);
		
		}
		
		
		this.alleles.put(variation_id, allele_id);
	}
	
	
	public int getID() {
		return this.id;
	}
	
	
	
	/**
	 * 
	 * Get the distance between two haplotypes. Either this is a variation missing in one of them or it is a different allele for a variation.
	 * 
	 * 
	 * @param haplotype
	 * @return
	 */
	public int getDistanceTo(Haplotype haplotype) {
		int distance = 0;
		
		for( Iterator<Integer> iterator = this.alleles.keySet().iterator(); iterator.hasNext();) {
			int variation_id = iterator.next();
			if(!haplotype.getAlleles().containsKey(variation_id) || this.alleles.get(variation_id).intValue() != haplotype.getAlleles().get(variation_id).intValue() ) {
				distance ++;
			}
			
		}
		
		for(Iterator<Integer> iterator = haplotype.getAlleles().keySet().iterator(); iterator.hasNext();) {
			int variation_id = iterator.next();
			if( !this.alleles.containsKey(variation_id)) {
				distance ++;
			}
		}
		
		return distance;
	}
	
	
	public int getDistanceToAllowHet(Haplotype haplotype, HashMap<Integer, Variation> variations) {
		int distance=0;
		
		for( Iterator<Integer> iterator = this.alleles.keySet().iterator(); iterator.hasNext();) {
			int variation_id = iterator.next();
			
			if(!haplotype.getAlleles().containsKey(variation_id) ) {
				
				distance ++;
				
			}else {
				Allele allele1 = variations.get(variation_id).getAllele(this.getAllele_id(variation_id));
				Allele allele2 = variations.get(variation_id).getAllele(haplotype.getAllele_id(variation_id));
				
				if(!allele1.isSameAlleleIgnoreHet(allele2)) {
					distance ++;
				}
			}
		}
		
		
		for(Iterator<Integer> iterator = haplotype.getAlleles().keySet().iterator(); iterator.hasNext();) {
			int variation_id = iterator.next();
			if( !this.alleles.containsKey(variation_id)) {
				distance ++;
			}
		}
		
		return distance;
	}
	
	
	
	
	
	
	/**
	 * 
	 * Check if any of the alleles of this haplotype is not resolved, i.e. is of the format HET(A/T)
	 * 
	 * 
	 * 
	 * @param variations
	 * @return
	 */
	public boolean hasHetAllele(HashMap<Integer, Variation> variations) {
		
		for(Iterator<Integer> iterator = this.alleles.keySet().iterator(); iterator.hasNext();) {
			int variation_id = iterator.next();
			int allele_id = this.alleles.get(variation_id);
			if( variations.get(variation_id).getAllele(allele_id).getNucleotideSequence().startsWith("HET")) {
				return true;
			}
		}
		
		
		return false;
	}
	
	
	
}
