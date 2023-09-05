package simpleHunter;

public class SNP {

	
	String mutantLine;
	String contig;
	String referenceAllele;
	String largestAlternativeAllele;
	int position;
	int coverage;
	
	double referenceAlleleFrequency;
	double largestAlternativeAlleleFrequency;
	
	
	public SNP(String mutantLine, String contig, String referenceAllele,
			String largestAlternativeAllele, int position, int coverage,
			double referenceAlleleFrequency,
			double largestAlternativeAlleleFrequency) {
		super();
		this.mutantLine = mutantLine;
		this.contig = contig;
		this.referenceAllele = referenceAllele;
		this.largestAlternativeAllele = largestAlternativeAllele;
		this.position = position;
		this.coverage = coverage;
		this.referenceAlleleFrequency = referenceAlleleFrequency;
		this.largestAlternativeAlleleFrequency = largestAlternativeAlleleFrequency;
	}


	public String getMutantLine() {
		return mutantLine;
	}


	public String getContig() {
		return contig;
	}


	public String getReferenceAllele() {
		return referenceAllele;
	}


	public String getLargestAlternativeAllele() {
		return largestAlternativeAllele;
	}


	public int getPosition() {
		return position;
	}


	public int getCoverage() {
		return coverage;
	}


	public double getReferenceAlleleFrequency() {
		return referenceAlleleFrequency;
	}


	public double getLargestAlternativeAlleleFrequency() {
		return largestAlternativeAlleleFrequency;
	}
	
	
	
	public boolean isEMSMutation(){
		if (this.getReferenceAllele().equalsIgnoreCase("G") && this.getLargestAlternativeAllele().equalsIgnoreCase("A")){
			return true;
		}
		if (this.getReferenceAllele().equalsIgnoreCase("C") && this.getLargestAlternativeAllele().equalsIgnoreCase("T")){
			return true;
		}
		return false;
	}
	
}
