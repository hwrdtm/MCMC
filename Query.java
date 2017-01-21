import java.io.*;

/**
 *  This class can be used for reading in queries either from files or
 *  from standard input.  The query itself is stored in public fields.
 */
public class Query {

    /** The variable whose value is being queried. */
    public int queryVar;

    /** An array specifying the settings of the evidence variables.
     *  More specifically, <tt>evidence[i]</tt> gives the value that
     *  variable <tt>i</tt> is to be set to if it is one of the
     *  evidence variables for this query; otherwise, if variable
     *  <tt>i</tt> is not an evidence variable, then
     *  <tt>evidence[i]</tt> is set to -1.
     */
    public int evidence[] = null;

    /** Use this constructor to read a query from standard input. */
    public Query(BayesNet bn)
	throws IOException {
	this(bn, new InputStreamReader(System.in), "<stdin>");
    }

    /** Use this constructor to read a query from a data file. */
    public Query(BayesNet bn, String filename)
	throws FileNotFoundException, IOException {
	this(bn, new FileReader(filename), filename);
    }

    // a private constructor
    private Query(BayesNet bn, Reader in_rdr, String filename)
	throws IOException {
	in = new LineNumberReader(in_rdr);
	this.filename = filename;
	String line;

	while((line = in.readLine()) != null) {
	    line = line.trim( );
	    if (line.equals(""))
		continue;

	    String words[] = line.split("\\s+");

	    String w = words[0];

	    int v = 0;
	    while (v < bn.numVariables && !w.equals(bn.varName[v]))
		v++;

	    if (v >= bn.numVariables)
		throw_exception("expected variable");

	    if (evidence == null) {
		if (words.length != 1)
		    throw_exception("badly formatted line");
		queryVar = v;
		evidence = new int[bn.numVariables];
		for (int i = 0; i < bn.numVariables; i++)
		    evidence[i] = -1;
	    } else {
		if (words.length != 2)
		    throw_exception("badly formatted line");

		w = words[1];
		int i = 0;
		while (i < bn.numValues[v] && !w.equals(bn.valueNames[v][i]))
		    i++;
		if (i >= bn.numValues[v])
		    throw_exception("unrecognized variable value");
		evidence[v] = i;
	    }
	}

	in.close();
    }

    // private stuff
    private String filename;
    private LineNumberReader in;

    private void throw_exception(String err) {
	err = err + " at line " + in.getLineNumber() + " in file " + filename;
	throw new RuntimeException(err);
    }

}
