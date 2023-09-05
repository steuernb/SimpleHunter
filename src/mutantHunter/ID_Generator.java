package mutantHunter;

public class ID_Generator {

	int last_id;
	
	
	public ID_Generator() {
		last_id = 0;
	}
	
	public int next() {
		last_id++;
		return last_id;
	}
	
}
