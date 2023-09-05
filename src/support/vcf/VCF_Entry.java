package support.vcf;

import java.util.HashMap;
import java.util.Iterator;

public class VCF_Entry {

	
	//CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	
	String chrom;
	int pos;
	String id;
	String ref;
	String alt;
	double qual;
	String filter;
	
	String infoString;
	String formatString;
	
	HashMap<String, String> info;
	HashMap<String, Integer> format;
	
	
	HashMap<String, HashMap<Integer, String>> samples; //the format fields for particular samples. <field<column_number,value>>
	
	
	public VCF_Entry(String inputLine) {
		
		String[] splitInputString = inputLine.split("\t");
		this.chrom = splitInputString[0];
		this.pos = Integer.parseInt(splitInputString[1]);
		this.id = splitInputString[2];
		this.ref = splitInputString[3];
		this.alt = splitInputString[4];
		this.qual = Double.parseDouble(splitInputString[5]);
		this.filter = splitInputString[6];
		
		this.infoString = splitInputString[7];
		this.formatString = splitInputString[8];
		
		
		this.initializeSamples(splitInputString);
		
		
		
	}
	
	private void initializeSamples(String[] splitInputString) {
		//initialize the samples HashMap
				this.samples= new HashMap<String, HashMap<Integer, String>>();
				//fill it with the format fields
				for(Iterator<String> iterator = this.getFormat().keySet().iterator(); iterator.hasNext();) {
					String formatField = iterator.next();
					this.samples.put( formatField, new HashMap<Integer, String>() );
				}
				
				//add the format values for each sample
				for( int i = 9;i< splitInputString.length; i++) {
					if(splitInputString[i].equalsIgnoreCase(".")) {
						continue;
					}else {
						String[] sample = splitInputString[i].split(":");
						for(Iterator<String> iterator = this.getFormat().keySet().iterator(); iterator.hasNext();) {
							String formatField = iterator.next();
							this.samples.get( formatField ).put(i, sample[this.getFormat().get(formatField)]);
						}
					}
					
					
				}
	}
	
	public HashMap<String, Integer> getFormat(){
		if(this.format == null) {
			String[] split = formatString.split(":");
			this.format= new HashMap<String, Integer>();
			for( int i = 0; i< split.length; i++) {
				this.format.put(split[i], i);
			}
		}
		return this.format;
	}
	
	
	
	
	public String getFormatField(int sampleColumnNumber, String formatHeader) {
		return this.samples.get(formatHeader).get(sampleColumnNumber);
	}
	
	
	/**
	 * 
	 * 	Get the key/value entries of the info column of this vcf entry.
	 *  To safe compute time, this HashMap will only be created upon first call. 
	 * 
	 * @return
	 * 		Return the full HashMap with all info entries. The IDs are the keys and the entries are the values. Description of the fields are maintained in the class VCF_Reader
	 * 
	 */
	public HashMap<String, String> getInfo(){
		if( this.info == null) {
			this.info = new HashMap<String, String>();
			String[] split= infoString.split(";");
			for(int i = 0; i< split.length; i++) {
				String[] split2 = split[i].split("=");
				String key = split2[0];
				String value = split2[1];
				this.info.put(key, value);
			}
			
			
			
		}
		return this.info;
	}
	
	/**
	 * 
	 * Get the value of a specific field from the info column. If the key does not exist this returns null.
	 * 
	 * @param key
	 * 		The ID of the field that is requested.
	 * @return
	 * 		Value of the field from the info column specified with ID
	 */
	public String getInfo(String id) {
		if( this.getInfo().containsKey(id)) {
			return this.info.get(id);
		}else {
			return null;
		}
		
		
	}
	
	
	
	
	
	public String getChromosome() {
		return this.chrom;
	}
	
	public int getPosition() {
		return this.pos;
	}

	public String getChrom() {
		return chrom;
	}

	public int getPos() {
		return pos;
	}

	public String getId() {
		return id;
	}

	public String getRef() {
		return ref;
	}

	public String getAlt() {
		return alt;
	}

	public double getQual() {
		return qual;
	}

	public String getFilter() {
		return filter;
	}

	public String getInfoString() {
		return infoString;
	}

	public HashMap<String, HashMap<Integer, String>> getSamples() {
		return samples;
	}

	public int getNumberOfAltAlleles() {
		return this.alt.split(",").length;
	}
	
	
	public boolean isHomozygous(int sampleColumnNumber) {
		String gt = this.samples.get("GT").get(sampleColumnNumber);
		if(gt == null) {
			throw new NullPointerException("no entry in column " + sampleColumnNumber);
		}
		if( gt.substring(0, 1).equalsIgnoreCase(gt.substring(2)) ) {
			return true;
		}else {
			return false;
		}
		
	}
	
	
	
	/**
	 * 
	 * get the allele for a specific column/sample
	 * 
	 * 
	 * 
	 * @param sampleColumnNumber
	 * @return
	 */
	public String getAllele(int sampleColumnNumber) {
		String gt = this.samples.get("GT").get(sampleColumnNumber);
		if(gt == null) {
			return "(NA)";
		}
		
		HashMap<Integer, String> alleles = this.getAlleles();
		boolean isHomozygous = this.isHomozygous(sampleColumnNumber);
		
		if( ! isHomozygous) {
			return "HET(" + alleles.get(Integer.parseInt(gt.substring(0, 1))) + "/" +  alleles.get(Integer.parseInt(gt.substring(2) ) ) +")"  ;
		}else {
			return alleles.get(Integer.parseInt(gt.substring(2)));
		}
			
			
		
		
		
		
	}
	
	
	public HashMap<Integer, String> getAlleles(){
		HashMap<Integer, String> alleles = new HashMap<Integer, String>();
		alleles.put(0, this.getRef());
		
		String[] split = this.getAlt().split(",");
		for( int i = 0; i< split.length; i++) {
			alleles.put( (i+1), split[i]);
		}
		return alleles;
	}
	
	
	
	
	
	
}
