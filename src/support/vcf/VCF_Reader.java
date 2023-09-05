package support.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;

public class VCF_Reader {

	
	BufferedReader in; //The object reading the vcf
	String lastLine;
	
	HashMap<String, String> infoFields;
	HashMap<String, String> formatFields;
	
	HashMap<Integer, String> samples; //the keys are the column numbers, the values the sample ids
	
	
	
	public VCF_Reader(File inputVCFFile)throws IOException{
		this.initialize(inputVCFFile);
	}
	
	
	private void initialize(File inputFile)throws IOException{
		
		
		
		//check if the fastq is gzipped
		FileInputStream fis = new FileInputStream(inputFile);
		byte[] bytes = new byte[2];
		fis.read(bytes);
		int head = ((int) bytes[0] & 0xff) | ((bytes[1] << 8) & 0xff00);
		boolean gzip = GZIPInputStream.GZIP_MAGIC == head;
		fis.close();

		
		if(gzip){
			this.in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputFile))));
		}else{
			this.in = new BufferedReader((new FileReader(inputFile)));
		}
		
		this.infoFields = new HashMap<String, String>();
		this.formatFields = new HashMap<String, String>();
		this.samples = new HashMap<Integer,String>();
		
		
		this.lastLine = in.readLine();
		
		while(this.lastLine.startsWith("#")) {
			
			
			
			// record the info fields
			if( this.lastLine.startsWith("##INFO")) {
				//##INFO=<ID=LEN,Number=A,Type=Integer,Description="allele length">
				String id = this.lastLine.split("ID=")[1].split(",")[0];
				String description = this.lastLine.split("Description=")[1].split(",")[0];
				infoFields.put(id, description);
			}
			
			
			// record the format fields
			if( this.lastLine.startsWith("##FORMAT")) {
				String id = this.lastLine.split("ID=")[1].split(",")[0];
				String description = this.lastLine.split("Description=")[1].split(",")[0];
				infoFields.put(id, description);
			}
			
			
			//read the sample names and record the columns they are in.
			if( this.lastLine.startsWith("#CHROM")) {
				String[] split = this.lastLine.split("\t");
				for( int i = 9; i< split.length; i++ ) {
					this.samples.put(i, split[i]);
				}
			}
			
			
			this.lastLine = in.readLine();
		}
		
		
		
		
	}
	
	
	
	
	
	/**
	 * 
	 * Get the header of a sample column; The counting starts at 0, so the first sample column is 9.
	 * 
	 * @param columnNumber
	 * @return
	 */
	public String getSample(int columnNumber) {
		return this.samples.get(columnNumber);
	}
	
	
	/**
	 * 
	 * 
	 * 
	 * @return
	 * @throws IOException
	 */
	public VCF_Entry next()throws IOException {
		
		if(this.lastLine == null) {
			this.in.close(); //automatically close the reading stream once every line is read.
			return null;
		}
		
		VCF_Entry entry = new VCF_Entry(this.lastLine);
		this.lastLine = this.in.readLine();
		return entry;
	}
	
	
	/**
	 * 
	 * "manually" close the reading stream. Note that if method next() returns null, the stream is closed automatically.
	 * 
	 * @throws IOException
	 */
	
	public void close() throws IOException{
		this.in.close();
	}
	public HashMap< Integer,String> getSamples(){
		return this.samples;
	}
	
	
	
	
}
