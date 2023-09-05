package mutantHunter;

import java.util.HashMap;
import java.util.Iterator;

import support.vcf.VCF_Entry;

public class Variation {

	int id;
	String chromosome;
	int position;
	HashMap<Integer, Allele> alleles;
	
	
	
	
	
	public Variation(int id, String chromosome, int position) {
		super();
		this.id = id;
		this.chromosome = chromosome;
		this.position = position;
		this.alleles = new HashMap<Integer, Allele>();
	}	
	
	public Variation(int id, VCF_Entry vcf_entry, ID_Generator allele_ids) {
		super();
		this.id = id;
		this.chromosome = vcf_entry.getChromosome();
		this.position = vcf_entry.getPosition();
		this.alleles = new HashMap<Integer, Allele>();
		
		Allele refAllele = new Allele(allele_ids.next(),id, vcf_entry.getRef(), true);
		this.addAllele(refAllele);
		String[] split = vcf_entry.getAlt().split(",");
		for( int i = 0; i< split.length; i++) {
			this.addAllele(new Allele(allele_ids.next(), id  ,split[i]));
		}
		
		
		
	}	
	
	public void addAllele(Allele allele) {
		this.alleles.put(allele.getId(), allele);
	}
	
	
	/**
	 * 
	 * Get the id of a specific allele.
	 * 
	 * @param nucleotideSequence
	 * @return
	 */
	public int getAlleleID(String nucleotideSequence) {
		for(Iterator<Integer> iterator = this.alleles.keySet().iterator(); iterator.hasNext();) {
			int allele_id = iterator.next();
			Allele allele = this.alleles.get(allele_id);
			
			if( allele.getNucleotideSequence().equalsIgnoreCase(nucleotideSequence)) {
				return allele_id;
			}
		}
		
		
		
		return -1;
	}
	
	
	public int getID() {
		return this.id;
	}

	

	public String getChromosome() {
		return chromosome;
	}

	public int getPosition() {
		return position;
	}
	
	public Allele getAllele(int allele_id) {
		return this.alleles.get(allele_id);
	}
	
	
	
}
