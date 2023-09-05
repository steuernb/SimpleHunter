package mutantHunter;

public class Accession {

	int id;
	
	String sample_name;  //this should correspond to the column header in the vcf file. So that is the readgroup id as well.
	String description;
	
	int vcfColumnNumber;
	
	
	public Accession(int id, String sample_id) {
		this.id = id;
		this.sample_name = sample_id;
	}
	
	public int getID() {
		return this.id;
	}
	
	public String getSample_id() {
		return sample_name;
	}
	
	
	public String getDescription() {
		return this.description;
	}
	
	public void setDescription(String description) {
		this.description = description;
	}
	
	public void setVcfColumnNumber(int vcfColumnNumber) {
		this.vcfColumnNumber = vcfColumnNumber;
	}
	public int getVcfColumnNumber() {
		return this.vcfColumnNumber;
	}
	
	
}
