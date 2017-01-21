import java.util.*;
import java.io.*;

/**
 *  This class loads a Bayes net from a data file.  A description of
 *  the topology of the Bayes net is stored in public fields, while
 *  the conditional probability tables are accessed using a public
 *  method.  Internally, variables are identified with the integers 0,
 *  1, 2, ..., up to (but not including) the number of variables.  The
 *  values that each variable can take are also given as integers in a
 *  similar manner.  */
public class BayesNet {

    /** The number of variables. */
    public int numVariables = 0;

    /** The names of each of the variables.  <tt>varName[i]</tt> is the
     *   name of variable <tt>i</tt>.
     */
    public String varName[];

    /** The number of values that each variable can take.
     *  <tt>numValues[i]</tt> is the number of values that variable
     *  <tt>i</tt> can take.
     */
    public int numValues[];

    /** The names of the values of each variable.
     *  <tt>valueNames[i][j]</tt> is the name of variable <tt>i</tt>'s
     *  <tt>j</tt>-th value.
     */
    public String valueNames[][];

    /** The number of parents of each variable.
     *  <tt>numParents[i]</tt> is the number of parents of variable
     *  <tt>i</tt>.
     */
    public int numParents[];

    /** The parents of each of the variables.  <tt>parents[i][]</tt>
     *  is a list of the parents of variable <tt>i</tt>.
     */
    public int parents[][];


    /** This method can be used to access the conditional probability
     *  table of each of the variables.  The method returns the
     *  conditional probability of a specified variable <tt>x</tt>
     *  taking a particular value, given the settings of variable
     *  <tt>x</tt>'s parents.  The first argument to this method is
     *  simply the variable <tt>x</tt> whose (conditional) probability
     *  is to be looked up.  The second argument <tt>vals[]</tt> is an
     *  integer array of length equal to the total number of
     *  variables, where <tt>vals[i]</tt> is a provided value for
     *  variable <tt>i</tt>.  Values must be provided in this array
     *  for <tt>x</tt> and for all of <tt>x</tt>'s parents; all other
     *  values in this array will be ignored.  The value returned by
     *  the method is the conditional probability of variable
     *  <tt>x</tt> taking the value <tt>vals[x]</tt>, given that all
     *  of <tt>x</tt>'s parents are set to the values specified in the
     *  <tt>vals[]</tt> array.
     */
    public double getCondProb(int x, int vals[]) {
	int idx = 0;
	for (int i = 0; i < numParents[x]; i++) {
	    int p = parents[x][i];
	    if (vals[p] < 0 || vals[p] >= numValues[p])
		throw new
		    IndexOutOfBoundsException("value of variable " + p +
					      " out of bounds");
	    idx = (idx * numValues[p]) + vals[p];
	}
	return cpt[x][idx][vals[x]];
    }


    /* Minimum value of any probability.  All values in the datafile
     * less than this value are replaced by it.
     */
    private final double MIN_PROB_VALUE = 1e-50;

    /**
     *  A constructor for constructing a Bayes net by reading in a
     *  description from the given data file.
     */
    public BayesNet(String filename)
	throws FileNotFoundException, IOException {

	in = new LineNumberReader(new FileReader(filename));
	this.filename = filename;

	building_cpt = false;
	VarRecord vr = null;

	String line;
	HashMap<String,VarRecord> var_table = new HashMap<String,VarRecord>();

	while((line = in.readLine()) != null) {
	    line = line.trim( );
	    if (line.equals(""))
		continue;

	    String words[] = line.split("\\s+");

	    if (words[0].equals("variable")) {
		if (building_cpt) finish_cpt(vr);
		if (words.length == 1)
		    throw_exception("expected variable name");
		String new_var = words[1];
		if (words.length == 2)
		    throw_exception("variable \"" + new_var +
				    "\" must have at least one value");
		vr = new VarRecord();
		vr.name = new_var;
		vr.num_values = words.length - 2;
		vr.values = new String[vr.num_values];
		for(int i = 0; i < vr.num_values; i++)
		    vr.values[i] = words[i+2];
		vr.id = numVariables++;
		if (var_table.containsKey(new_var))
		    throw_exception("variable \"" + new_var +
				    "\" appears second time");
		var_table.put(new_var, vr);
	    } else if (words[0].equals("probability")) {
		if (building_cpt) finish_cpt(vr);
		String var = null;
		if (words.length == 1
		    || !var_table.containsKey(var = words[1]))
		    throw_exception("expected variable name");
		vr = var_table.get(var);
		if (vr.num_parents >= 0)
		    throw_exception("trying to build probability table for "+
				    "variable \"" + var + "\" a second time");
		if (words.length == 2) {
		    vr.num_parents = 0;
		    vr.parents = new VarRecord[0];
		} else {
		    if (words.length == 3
			|| !words[2].equals("|"))
			throw_exception("expected '|' following variable name");
		    vr.num_parents = words.length - 3;
		    vr.parents = new VarRecord[vr.num_parents];
		    for (int i = 0; i < vr.num_parents; i++) {
			String w = words[i+3];
			if (!var_table.containsKey(w))
			    throw_exception("variable \"" + w +
					    "\" not yet declared");
			vr.parents[i] = var_table.get(w);
		    }
		}
		int num_settings = 1;
		for (int i = 0; i < vr.num_parents; i++)
		    num_settings *= vr.parents[i].num_values;
		vr.cpt = new double[num_settings][];
		building_cpt = true;
	    } else if (building_cpt) {
		if (words.length != vr.num_parents + vr.num_values)
		    throw_exception("badly formatted data");
		int idx = 0;
		int wc = 0;
		for (int i = 0; i < vr.num_parents; i++) {
		    String w = words[wc++];
		    int par_nvals = vr.parents[i].num_values;
		    idx *= par_nvals;
		    int j;
		    for (j = 0; j < par_nvals; j++) {
			if (vr.parents[i].values[j].equals(w))
			    break;
		    }
		    if (j >= par_nvals)
			throw_exception("\"" + w + "\" is not a declared " +
					"value of \"" + vr.parents[i].name
					+ "\"");
		    idx += j;
		}
		vr.cpt[idx] = new double[vr.num_values];
		double sum = 0.0;
		double d = 0.0;
		for (int i = 0; i < vr.num_values; i++) {
		    try {
			d = Double.parseDouble(words[wc++]);
		    } catch (NumberFormatException e) {
			throw new NumberFormatException(
				    "expected double in column " + wc +
				    " at line " + in.getLineNumber() +
				    " in file " + filename);
		    }
		    if (d < 0.0)
			throw_exception("probabilities cannot be negative");
		    if (d < MIN_PROB_VALUE)
			d = MIN_PROB_VALUE;
		    sum += vr.cpt[idx][i] = d;
		}
		if (sum == 0.0)
		    throw_exception("probabilities must sum to positive number");
		for (int i = 0; i < vr.num_values; i++)
		    vr.cpt[idx][i] /= sum;
	    } else
		throw_exception("badly formatted data");
	}
	if (building_cpt) finish_cpt(vr);

	varName = new String[numVariables];
	numValues = new int[numVariables];
	valueNames = new String[numVariables][];
	numParents = new int[numVariables];
	parents = new int[numVariables][];
	cpt = new double[numVariables][][];

	Iterator it = var_table.keySet().iterator();
	while (it.hasNext()) {
	    vr = var_table.get(it.next());
	    int i = vr.id;
	    if (vr.num_parents < 0)
		throw new RuntimeException(
                       "variable \"" + vr.name + "\" was never " +
		       "assigned probabilities");
	    varName[i] = vr.name;
	    numValues[i] = vr.num_values;
	    valueNames[i] = vr.values;
	    numParents[i] = vr.num_parents;
	    parents[i] = new int[vr.num_parents];
	    for(int j = 0; j < vr.num_parents; j++)
		parents[i][j] = vr.parents[j].id;
	    cpt[i] = vr.cpt;
	}

	in.close();

    }

    // private stuff

    private double cpt[][][];

    private boolean building_cpt = false;

    private void finish_cpt(VarRecord vr) {
	for (int i = 0; i < vr.cpt.length; i++)
	    if (vr.cpt[i] == null)
		throw new RuntimeException(
				"probability table for variable \"" +
				vr.name + "\" not completely specified " +
				"in file " + filename);
	building_cpt = false;
    }

    private void throw_exception(String err) {
	err = err + " at line " + in.getLineNumber() + " in file " + filename;
	throw new RuntimeException(err);
    }

    private class VarRecord {
	private String name;
	private int num_values;
	private String[] values;
	private int id;
	private int num_parents = -1;
	private VarRecord parents[];
	private double cpt[][];
    }

    private String filename;
    private LineNumberReader in;

}
