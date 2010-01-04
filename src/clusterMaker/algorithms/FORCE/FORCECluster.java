package clusterMaker.algorithms.FORCE;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import java.util.ArrayList;
import java.util.Date;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import javax.swing.JOptionPane;
import javax.swing.JPanel;

import cytoscape.Cytoscape;
import cytoscape.CyEdge;
import cytoscape.CyNode;
import cytoscape.data.CyAttributes;
import cytoscape.groups.CyGroup;
import cytoscape.groups.CyGroupManager;
import cytoscape.layout.Tunable;
import cytoscape.logger.CyLogger;
import cytoscape.task.TaskMonitor;

import clusterMaker.algorithms.ClusterAlgorithm;
import clusterMaker.algorithms.AbstractClusterAlgorithm;
import clusterMaker.ui.ClusterViz;

import de.layclust.layout.forcend.FORCEnDLayoutConfig;
import de.layclust.taskmanaging.TaskConfig;


public class FORCECluster extends AbstractClusterAlgorithm implements ActionListener {
	private CyLogger logger = null;
	private	TaskMonitor monitor = null;

	private CyAttributes cyNodeAttributes;
	private CyAttributes cyEdgeAttributes;

	private List<CyNode> nodes;
	private List<CyEdge> edges;
	
	private Hashtable<String, Double> similaritiesForGivenEdges = new Hashtable<String, Double>(); // key: s#t
	private Hashtable<String, Double> normalizedSimilaritiesForGivenEdges = new Hashtable<String, Double>(); // key: s#t
	
	// Tunable values
	private double threshold;
	private double maxSim = Double.NEGATIVE_INFINITY;
	private double minSim = Double.POSITIVE_INFINITY;
	private boolean isDistanceFunction = false;
	private boolean mergeSimilar = false;
	private double normalizedThreshold = 0;
	private double mergeThreshold = 0;
	private String edgeAttribute = null;
	private boolean evolutionaryTraining = false;
	private int dimension = 0;
	private int generations = 0;
	private double minThreshold = 0;
	private double maxThreshold = 0;
	private double stepSize = 0;

	private String groupAttributeCC = "__FORCEccGroups";
	private String groupAttribute = "__FORCEGroups";

	private String attributeNameCluster = "FORCECluster";
	private String attributeNameCC = "FORCE_connected_component";


	/**
 	 * Main constructor -- calls constructor of the superclass and initializes
 	 * all of our tunables.
 	 */
	public FORCECluster() {
		super();
		logger = CyLogger.getLogger(FORCECluster.class);
		initializeProperties();
	}

	public String getShortName() {return "force";};
	public String getName() {return "FORCE cluster";};

	public JPanel getSettingsPanel() {
		return clusterProperties.getTunablePanel();
	}

	/**
 	 * At this point, we don't have a specific visualizer, although
 	 * there might be some value at some juncture to use a heatmap to
 	 * view the resulting matrix.
 	 */
	public ClusterViz getVisualizer() {
		return null;
	}

	public void updateSettings() {
		updateSettings(false);
	}

	/**
 	 * Update all of our tunables
 	 */
	public void updateSettings(boolean force) {
		clusterProperties.updateValues();
		super.updateSettings(force);

		Tunable t = clusterProperties.get("edgeAttribute");
		if ((t != null) && (t.valueChanged() || force))
			edgeAttribute = (String) t.getValue();

		t = clusterProperties.get("weightType");
		if ((t != null) && (t.valueChanged() || force)) {
			int val = ((Integer) t.getValue()).intValue();
			if (val == 0) 
				this.isDistanceFunction = false;
			else
				this.isDistanceFunction = true;
		}

		t = clusterProperties.get("threshold");
		if ((t != null) && (t.valueChanged() || force))
			threshold = ((Double) t.getValue()).doubleValue();

		t = clusterProperties.get("dimension");
		if ((t != null) && (t.valueChanged() || force))
			dimension = ((Integer) t.getValue()).intValue();

		t = clusterProperties.get("iterations");
		if ((t != null) && (t.valueChanged() || force))
			FORCEnDLayoutConfig.iterations = ((Integer) t.getValue()).intValue();

		t = clusterProperties.get("attractionFactor");
		if ((t != null) && (t.valueChanged() || force))
			FORCEnDLayoutConfig.attractionFactor = ((Double) t.getValue()).doubleValue();

		t = clusterProperties.get("repulsionFactor");
		if ((t != null) && (t.valueChanged() || force))
			FORCEnDLayoutConfig.repulsionFactor = ((Double) t.getValue()).doubleValue();

		t = clusterProperties.get("mergeSimilar");
		if ((t != null) && (t.valueChanged() || force))
			mergeSimilar = ((Boolean) t.getValue()).booleanValue();

		t = clusterProperties.get("mergeThreshold");
		if ((t != null) && (t.valueChanged() || force))
			mergeThreshold = ((Integer) t.getValue()).intValue();

		t = clusterProperties.get("evolutionaryTraining");
		if ((t != null) && (t.valueChanged() || force))
			evolutionaryTraining = ((Boolean) t.getValue()).booleanValue();

		t = clusterProperties.get("generations");
		if ((t != null) && (t.valueChanged() || force))
			generations = ((Integer) t.getValue()).intValue();

		t = clusterProperties.get("componentAttribute");
		if ((t != null) && (t.valueChanged() || force))
			attributeNameCC = (String)t.getValue();

		t = clusterProperties.get("clusterAttribute");
		if ((t != null) && (t.valueChanged() || force))
			attributeNameCluster = (String)t.getValue();

		t = clusterProperties.get("minThreshold");
		if ((t != null) && (t.valueChanged() || force))
			minThreshold = ((Double) t.getValue()).doubleValue();

		t = clusterProperties.get("maxThreshold");
		if ((t != null) && (t.valueChanged() || force))
			maxThreshold = ((Double) t.getValue()).doubleValue();

		t = clusterProperties.get("stepsize");
		if ((t != null) && (t.valueChanged() || force))
			stepSize = ((Double) t.getValue()).doubleValue();

	}

	/**
 	 * Perform the actual clustering.  For FORCE, there are really
 	 * two steps:
 	 * 	1) Assign all of the connected components
 	 * 	2) Do the FORCE clustering.
 	 *
 	 * There is also an optional approach called evolutionary parameter
 	 * tuning, which takes a lot longer and is probably less relevant for
 	 * the Cytoscape integration.
 	 *
 	 * @param monitor the TaskMonitor to use
 	 */
	public void doCluster(TaskMonitor monitor) {
		this.monitor = monitor;

		updateSettings();

		cyNodeAttributes = Cytoscape.getNodeAttributes();
		cyEdgeAttributes = Cytoscape.getEdgeAttributes();

		cyNodeAttributes.deleteAttribute(this.attributeNameCluster);

		nodes = Cytoscape.getCurrentNetwork().nodesList();
		edges = Cytoscape.getCurrentNetwork().edgesList();

		try {
			// Calculate all of the connected components
			runCC();

			// Cluster
			runFORCE();
		} catch (IOException e) {
			StackTraceElement stack[] = e.getStackTrace();
			logger.error("IOException: "+e.getMessage()+" at: ");
			for (int i = 0; i < stack.length; i++) {
				logger.error("      "+stack[i].toString());
			}
		}
	}

	/**
 	 * This is called if the user requests evolutionary parameter
 	 * tuning.
 	 */
	public void actionPerformed(ActionEvent e) {
	}

	/**
 	 * Initialize all of our tunables
 	 */
	protected void initializeProperties() {
		super.initializeProperties();
		// Attributes group
		clusterProperties.add(new Tunable("attributesGroup", 
		                                  "Attributes",
		                                  Tunable.GROUP, new Integer(3)));
		{
			// 	Edge weight attribute (Edge attribute combo-box)
    	clusterProperties.add(new Tunable("edgeAttribute",
	                                      "Edge weight attribute:",
 	                                     Tunable.EDGEATTRIBUTE, "",
 	                                     (Object)null, (Object)null, 
			                                  Tunable.NUMERICATTRIBUTE));
	
			// 	Edge weights correspond to: Similarity|Distance combo box
			String[] weightTypes = { "Similarity", "Distance" };
 	   clusterProperties.add(new Tunable("weightType",
 	                                     "Edge weights correspond to:",
 	                                     Tunable.LIST, new Integer(0),
 	                                     (Object)weightTypes, (Object)null, 0));
			// 	Threshold integer
 	   clusterProperties.add(new Tunable("threshold",
 	                                     "Threshold:",
 	                                     Tunable.DOUBLE, new Double(10)));
		}

		// Advanced attributes group
		clusterProperties.add(new Tunable("advancedAttributesGroup", 
		                                  "Advanced Attributes",
		                                  Tunable.GROUP, new Integer(4),
		                                  new Boolean(true), null, Tunable.COLLAPSABLE));
		{
			// 	Layouter options group
			clusterProperties.add(new Tunable("layouterOptionsGroup", 
			                                  "Layouter Options",
			                                  Tunable.GROUP, new Integer(4)));
			{
				// 		Dimension (slider 2-6)
		    clusterProperties.add(new Tunable("dimension",
   		 	                                  "Dimension:",
     			                                 Tunable.INTEGER, new Integer(3),
		   			                               new Integer(2), new Integer(6), Tunable.USESLIDER));
		
				// 		Iterations (slider 50-250)
 		   	clusterProperties.add(new Tunable("iterations",
   		 	                                  "Iterations:",
   		 	                                  Tunable.INTEGER, new Integer(FORCEnDLayoutConfig.iterations),
		 			                                 new Integer(50), new Integer(250), Tunable.USESLIDER));
				// 		Attraction factor (double)
   		 	clusterProperties.add(new Tunable("attractionFactor",
     			                                 "Attraction Factor:",
       			                               Tunable.DOUBLE, new Double(FORCEnDLayoutConfig.attractionFactor)));
				// 		Repulsion factor (double)
   		 	clusterProperties.add(new Tunable("repulsionFactor",
     			                                 "Repulsion Factor:",
       			                               Tunable.DOUBLE, new Double(FORCEnDLayoutConfig.repulsionFactor)));
			}
		
			//	Merge nodes
			clusterProperties.add(new Tunable("mergeNodesGroup", 
			                                  "Merge nodes",
			                                  Tunable.GROUP, new Integer(2)));
			{
				//		Merge very similar nodes to one (boolean)
				clusterProperties.add(new Tunable("mergeSimilar", 
				                                  "Merge very similar nodes to one?",
		 		                                  Tunable.BOOLEAN, new Boolean(false)));
				//		Threshold (integer)
				clusterProperties.add(new Tunable("mergeThreshold", 
				                                  "Threshold:",
		 		                                  Tunable.INTEGER, new Integer(100)));
			}

			//	Parameter training
			clusterProperties.add(new Tunable("parameterGroup", 
			                                  "Parameter training",
			                                  Tunable.GROUP, new Integer(2)));
			{
				//		Evolutionary parameter training? (boolean)
				clusterProperties.add(new Tunable("evolutionaryTraining", 
				                                  "Evolutionary parmeter training?",
		 		                                  Tunable.BOOLEAN, new Boolean(false)));
				//		Generations (1-10)
				clusterProperties.add(new Tunable("generations", 
				                                  "Generations:",
		 		                                  Tunable.INTEGER, new Integer(3),
				                                  new Integer(1), new Integer(10), Tunable.USESLIDER));
			}

			//	Cytoscape
			clusterProperties.add(new Tunable("cytoscapeGroup", 
			                                  "Cytoscape",
			                                  Tunable.GROUP, new Integer(2)));
			{
				//		Node attribute connected component (string)
				clusterProperties.add(new Tunable("componentAttribute", 
			 		                                "Node attribute connected component:",
			 		                                Tunable.STRING, attributeNameCC));
				//		Node attribute FORCE cluster (string)
				clusterProperties.add(new Tunable("clusterAttribute", 
			 		                                "Node attribute FORCE cluster:",
			 		                                Tunable.STRING, attributeNameCluster));
			}
		}
		// Determine optimal threshold
		clusterProperties.add(new Tunable("optimalGroup", 
		                                  "Determine optimal threshold",
		                                  Tunable.GROUP, new Integer(5)));
		{
			//	 Gold standard attribute (node attribute)
    	clusterProperties.add(new Tunable("goldStandardAttributes",
	                                      "Gold standard attribute:",
 	                                     Tunable.NODEATTRIBUTE, "",
 	                                     (Object)null, (Object)null, 
			                                  Tunable.NUMERICATTRIBUTE));
			//	 Minimal threshold (integer)
			clusterProperties.add(new Tunable("minThreshold", 
			                                  "Minimal Threshold:",
	 		                                  Tunable.DOUBLE, new Double(10)));
			//	 Maximal threshold (integer)
			clusterProperties.add(new Tunable("maxThreshold", 
			                                  "Maximal Threshold:",
	 		                                  Tunable.DOUBLE, new Double(100)));
			// 	 Stepsize (float)
			clusterProperties.add(new Tunable("stepsize", 
			                                  "Stepsize:",
	 		                                  Tunable.DOUBLE, new Double(0.5)));
			// 	 Run Comparison (button)
			clusterProperties.add(new Tunable("runComparison", 
			                                  "",
	 		                                  Tunable.BUTTON, "Run Comparison", 
			                                  this, null, 0));
		}
	}

	private boolean runCC() {

		cyNodeAttributes = Cytoscape.getNodeAttributes();
		cyEdgeAttributes = Cytoscape.getEdgeAttributes();
		
		int step = 0;
		int max = this.edges.size()*2 + this.nodes.size();
		monitor.setStatus("Assigning connected components...");
		monitor.setPercentCompleted(0);
		
		if (edgeAttribute == null) {
			JOptionPane.showMessageDialog(Cytoscape.getDesktop(), 
				"ERROR - You must select an edge attribute! ");
			return false;
		}
		byte type = cyEdgeAttributes.getType(edgeAttribute);
		
		for (CyEdge e: edges) {
			
			String sourceID = e.getSource().getIdentifier();
			String targetID = e.getTarget().getIdentifier();
			
			// THIS IS A STRANGE BUG IN CYTOSCAPE; it always give u 2 additional nodes :-(
			if (!sourceID.equalsIgnoreCase("Source") && 
			    !sourceID.equalsIgnoreCase("Target") && 
			    !targetID.equalsIgnoreCase("Source") && 
			    !targetID.equalsIgnoreCase("Target")) {
				double sim = 0;
				if (type == CyAttributes.TYPE_INTEGER) {
					try {
						sim = (double) cyEdgeAttributes.getIntegerAttribute(e.getIdentifier(), edgeAttribute);
					} catch(NullPointerException e1) {
						continue;
					}
				} else {
					try {
						sim = cyEdgeAttributes.getDoubleAttribute(e.getIdentifier(), edgeAttribute);
					} catch(NullPointerException e1) {
						continue;
					}
				}
				
				similaritiesForGivenEdges.put((sourceID + "#" + targetID), sim);
				
				if (sim < this.minSim) {
					this.minSim = sim;
				}
				if (sim > this.maxSim) {
					this.maxSim = sim;
				}
			}
			
			step++;
			if ((step%10) == 0) monitor.setPercentCompleted((step/max)*100);
		}
		
		if (!this.isDistanceFunction) { // is sim function
			this.normalizedThreshold = this.getNormalizedValue(this.minSim, this.maxSim, this.threshold);
		} else {
			this.normalizedThreshold = 100-this.getNormalizedValue(this.minSim, this.maxSim, this.threshold);
		}
		
		for (String key: similaritiesForGivenEdges.keySet()) {

			double sd = similaritiesForGivenEdges.get(key);
			double s;
			if (!this.isDistanceFunction) { // is sim function
				s = this.getNormalizedValue(this.minSim, this.maxSim, sd);
			} else {
				s = 100-this.getNormalizedValue(this.minSim, this.maxSim, sd);
			}
			
			normalizedSimilaritiesForGivenEdges.put(key, s);
			
			step++;
			if ((step%10) == 0) monitor.setPercentCompleted((step/max)*100);
		}
				
		Hashtable<String, Boolean> already = new Hashtable<String, Boolean>();
		
		int clusterNr = 0;
		
		String networkID = Cytoscape.getCurrentNetwork().getIdentifier();
		CyAttributes netAttributes = Cytoscape.getNetworkAttributes();

		if (netAttributes.hasAttribute(networkID, groupAttributeCC)) {
			List<String> groupList = (List<String>)netAttributes.getListAttribute(networkID, groupAttributeCC);
			for (String groupName: groupList) {
				CyGroup group = CyGroupManager.findGroup(groupName);
				if (group != null)
					CyGroupManager.removeGroup(group);
			}
		}
		
		List<String>groupList = new ArrayList<String>();
		for (CyNode node:nodes) {
			String id = node.getIdentifier();
			if (!already.containsKey(id)) {
				Vector<CyNode> nodesInThisCluster = new Vector<CyNode>();
				assignNodeToCluster(already, node, id, clusterNr, nodesInThisCluster);
				clusterNr++;
				String groupName = this.attributeNameCC +"_" + clusterNr;
				// groupList.add(groupName);
				
				CyGroup newgroup = CyGroupManager.createGroup(groupName, nodesInThisCluster, null);
				if (newgroup != null) {
					CyGroupManager.setGroupViewer(newgroup, "metaNode", Cytoscape.getCurrentNetworkView(), true);
					groupList.add(groupName);
				}
				
			}
			
			step++;
			if ((step%10) == 0) monitor.setPercentCompleted((step/max)*100);
		}
		
		netAttributes.setListAttribute(networkID, groupAttributeCC, groupList);
		return true;
	}

	/**
	 * Execute the FORCE algorithm on the current dataset.  Note that this does
	 * assume that the connected components have already been calculated.
	 */
	private void runFORCE() throws IOException {
		// Have we already calculated the connected components?
		// No, calculate them
		
		Vector<Vector<CyNode>> connectedComponents = new Vector<Vector<CyNode>>();
		
		for (CyNode n: nodes) {
			int ccNr = cyNodeAttributes.getIntegerAttribute(n.getIdentifier(), attributeNameCC);
			try {
				connectedComponents.get(ccNr).add(n);
			} catch (Exception e) {
				Vector<CyNode> v = new Vector<CyNode>();
				v.add(n);
				connectedComponents.add(ccNr, v);
			}
		}
		
		monitor.setPercentCompleted(0);
		monitor.setStatus("FORCE");
		
		Date date = new Date(System.currentTimeMillis());
		String dateTimeStr = date.toString().replaceAll(" ", "_");
		dateTimeStr = dateTimeStr.replaceAll(":", "-");
		
		File cmSubTempDir = File.createTempFile("cm_"+dateTimeStr, null);
		// This creates the temp file, but we want a directory
		cmSubTempDir.delete(); // Delete the tmp file

		// Now make a directory
		boolean suc = cmSubTempDir.mkdir();
		if (!suc) {
			throw new IOException("Can't write to temp directory: "+cmSubTempDir.toString());
		}
		
		for (int i = 0; i < connectedComponents.size(); i++) {
			monitor.setPercentCompleted((i/connectedComponents.size()*2)*100);
			monitor.setStatus("Writing temp file nr. " + i + " of " + connectedComponents.size());

			
			Vector<CyNode> cc = connectedComponents.get(i);
			
			if(mergeSimilar){
				writeCCtoTempDirWithMergedNodes(cmSubTempDir, cc, i);
			}else{
				writeCCtoTempDir(cmSubTempDir, cc, i);
			}
			
			
		}
		
		monitor.setStatus("Running FORCE clustering (might take a while)...");
		
		File resultsFileName = new File(cmSubTempDir, "_results.txt");

		if(evolutionaryTraining) {
			String[] args = {"-i",cmSubTempDir.toString(),
			                 "-o",resultsFileName.toString(),
			                 "-cf", "FALSE",
			                 "-ld", "" + dimension,
			                 "-fi", "" + FORCEnDLayoutConfig.iterations,
			                 "-lp", TaskConfig.parameterTrainingClass,
			                 "-lpn", ""+generations,
			                 "-fa", ""+FORCEnDLayoutConfig.attractionFactor,
			                 "-fr", ""+FORCEnDLayoutConfig.repulsionFactor};
			FORCEnD_ACC.main(args);
			TaskConfig.doLayoutParameterTraining = false;
		} else {
			String[] args = {"-i", cmSubTempDir.toString(),
			                 "-o", resultsFileName.toString(),
			                 "-cf", "FALSE",
			                 "-ld", "" + dimension,
			                 "-fi", "" + FORCEnDLayoutConfig.iterations,
			                 "-fa", "" + FORCEnDLayoutConfig.attractionFactor,
			                 "-fr","" + FORCEnDLayoutConfig.repulsionFactor};
			FORCEnD_ACC.main(args);
		}
		
		
		readFORCEresults(resultsFileName);
		
		deleteDirectory(cmSubTempDir);
		resultsFileName.delete();
		
	}


	/**
 	 * As part of the FORCE algorithm, we write out information into a temporary file.
 	 * This routine writes that file with the merged nodes.
 	 *
 	 * @param cmTempDir a File representing the temporary directory
 	 * @param cc the Vector of CyNodes
 	 * @param ccNr the component number of this file
 	 */
	private void writeCCtoTempDirWithMergedNodes(File cmTempDir, 
	                                             Vector<CyNode> cc, 
	                                             int ccNr) throws IOException {
		
		double normalizedUpperBound = 0;
		
		if (!isDistanceFunction) { // is sim function
			normalizedUpperBound = getNormalizedValue(minSim, maxSim, mergeThreshold);
		} else {
			normalizedUpperBound = 100-getNormalizedValue(minSim, maxSim, mergeThreshold);
		}
		
		// find connected components with upperBound to merge them to one node
		Vector<Vector<CyNode>> upperBoundMergedNodes = new Vector<Vector<CyNode>>();
		
		Hashtable<String, Boolean> already = new Hashtable<String, Boolean>();
		
		for (int j = 0; j < cc.size(); j++) {
			
			CyNode n = cc.get(j);
			
			if(!already.containsKey(n.getIdentifier())){
			
				Vector<CyNode> v = new Vector<CyNode>();
				
				findMergeNodes(n,cc,v,already,normalizedUpperBound);
				upperBoundMergedNodes.add(v);
				
			}
		}

		File tmpFile = new File(cmTempDir, "costmatrix_nr_"+ccNr+"_size+"+cc.size()+".rcm");
		BufferedWriter bw = new BufferedWriter(new FileWriter(tmpFile));
		
		bw.write("0");
		bw.newLine();
		bw.write("" + upperBoundMergedNodes.size());
		bw.newLine();
		for (int j = 0; j < upperBoundMergedNodes.size(); j++) {
			Vector<CyNode> v = upperBoundMergedNodes.get(j);
			for (int k = 0; k < v.size()-1; k++) {
				bw.write(v.get(k).getIdentifier() + "\t");
			}
			bw.write(v.get(v.size()-1).getIdentifier());
			bw.newLine();
		}
		for (int j = 0; j < upperBoundMergedNodes.size(); j++) {
			Vector<CyNode> v = upperBoundMergedNodes.get(j);
			for (int k = j+1; k < upperBoundMergedNodes.size(); k++) {			
				Vector<CyNode> v2 = upperBoundMergedNodes.get(k);
				double sim = calculateSimilarityForMergeNodes(v,v2);
				bw.write(Double.toString(sim));
				if(k<upperBoundMergedNodes.size()-1){
					bw.write("\t");
				}else if(j<upperBoundMergedNodes.size()-1){
					bw.write("\n");
				}
			}	
		}
		bw.close();
	}
		
	/**
 	 * As part of the FORCE algorithm, we write out information into a temporary file.
 	 * This routine writes that file without mering the nodes.
 	 *
 	 * @param cmTempDir a File representing the temporary directory
 	 * @param cc the Vector of CyNodes
 	 * @param ccNr the component number of this file
 	 */
	private void writeCCtoTempDir(File cmTempDir, 
	                              Vector<CyNode> cc, int ccNr) throws IOException {

		File tmpFile = new File(cmTempDir, "costmatrix_nr_"+ccNr+"_size+"+cc.size()+".cm");
		BufferedWriter bw = new BufferedWriter(new FileWriter(tmpFile));
		
		bw.write("" + cc.size());
		bw.newLine();
		
		for (int i = 0; i < cc.size(); i++) {
			CyNode n = cc.get(i);
			bw.write(n.getIdentifier());
			bw.newLine();
		}
		
		for (int i = 0; i < cc.size(); i++) {
			String s = cc.get(i).getIdentifier();
			for (int j = i+1; j < cc.size(); j++) {
				String t = cc.get(j).getIdentifier();
				
				double cost = -this.normalizedThreshold;
				if (this.normalizedSimilaritiesForGivenEdges.containsKey(s + "#" + t)) {
					cost = this.normalizedSimilaritiesForGivenEdges.get(s + "#" + t) - this.normalizedThreshold;
				} else if (this.normalizedSimilaritiesForGivenEdges.containsKey(t + "#" + s)) {
					cost = this.normalizedSimilaritiesForGivenEdges.get(t + "#" + s) - this.normalizedThreshold;
				}
				
				if (j != cc.size()-1) {
					bw.write(cost + "\t");
				} else {
					bw.write("" + cost);
				}
			}
			
			if (i != cc.size()-1) {
				bw.newLine();
			}
		}
		bw.close();
	}

	private double getNormalizedValue(double min, double max, double value) {
		double span = max-min;
		return ((value-min)/span)*100;
	}

	private boolean deleteDirectory(File path) {
		if (path.exists()) {
			File[] files = path.listFiles();
			for (int i = 0; i < files.length; i++) {
				if (files[i].isDirectory()) {
					deleteDirectory(files[i]);
				} else {
					files[i].delete();
				}
			}
		}
		return (path.delete());
	}


	private double calculateSimilarityForMergeNodes(Vector<CyNode> v, Vector<CyNode> v2) {
		
		double value = 0;
		
		for (int i = 0; i < v.size(); i++) {
			CyNode n = v.get(i);
			for (int j = 0; j < v2.size(); j++) {
				CyNode n2 = v2.get(j);
				String key = n.getIdentifier() + "#" + n2.getIdentifier();
				String keyI = n2.getIdentifier() + "#" + n.getIdentifier();
				if(this.normalizedSimilaritiesForGivenEdges.containsKey(key)){
					value+=(this.normalizedSimilaritiesForGivenEdges.get(key)-this.threshold);
				}else if(this.normalizedSimilaritiesForGivenEdges.containsKey(keyI)){
					value+=(this.normalizedSimilaritiesForGivenEdges.get(keyI)-this.threshold);
				}
			}
		}
		return value;
	}

	private void findMergeNodes(CyNode n, Vector<CyNode> cc, 
	                            Vector<CyNode> v, Hashtable<String, Boolean> already, 
	                            double normalizedUpperBound) {
		
		v.add(n);
		already.put(n.getIdentifier(), true);
		
		for (int i = 0; i < cc.size(); i++) {
			
			CyNode n2 = cc.get(i);
			
			String key = n.getIdentifier() + "#" + n2.getIdentifier();
			String keyI = n2.getIdentifier() + "#" + n.getIdentifier();
			
			if(!already.containsKey(n2.getIdentifier())){
				if(this.normalizedSimilaritiesForGivenEdges.containsKey(key)){
					if(this.normalizedSimilaritiesForGivenEdges.get(key)<normalizedUpperBound){
						findMergeNodes(n2, cc, v, already, normalizedUpperBound);
					}
				} else if(this.normalizedSimilaritiesForGivenEdges.containsKey(keyI)){
					if(this.normalizedSimilaritiesForGivenEdges.get(keyI)<normalizedUpperBound){
						findMergeNodes(n2, cc, v, already, normalizedUpperBound);
					}
				}
			}
		}
	}

	
	private void readFORCEresults(File resultsFileName) throws IOException {
		
		Hashtable<String, Integer> clusterForGivenNode = new Hashtable<String, Integer>();
		
		BufferedReader br = new BufferedReader(new FileReader(resultsFileName));
		
		String line;
		while ((line=br.readLine()) != null) {
			
			String[] d = line.split("\t");
			
			clusterForGivenNode.put(d[0].trim(), Integer.parseInt(d[1].trim()));
		}
		
		br.close();
		
		Hashtable<Integer, Vector<CyNode>> nodeListForGivenClusterNumber = new Hashtable<Integer, Vector<CyNode>>();
		
		for (CyNode n: nodes) {
			
			int clusterNr =  clusterForGivenNode.get(n.getIdentifier());
			cyNodeAttributes.setAttribute(n.getIdentifier(), this.attributeNameCluster, clusterNr);
			
			Vector<CyNode> v = new Vector<CyNode>();
			if (nodeListForGivenClusterNumber.containsKey(clusterNr)) {
				v = nodeListForGivenClusterNumber.get(clusterNr);
			}
			v.add(n);
			nodeListForGivenClusterNumber.put(clusterNr, v);
		}
		
		// See if we already have groups defined (from a previous run?)
		CyAttributes netAttributes = Cytoscape.getNetworkAttributes();
		String networkID = Cytoscape.getCurrentNetwork().getIdentifier();
		if (netAttributes.hasAttribute(networkID, this.groupAttribute)) {
			List<String> groupList = (List<String>)netAttributes.getListAttribute(networkID, this.groupAttribute);
			for (String groupName: groupList) {
				CyGroup group = CyGroupManager.findGroup(groupName);
				if (group != null)
					CyGroupManager.removeGroup(group);
			}
		}
		
		// Now, create the groups
		List<String>groupList = new ArrayList<String>();
		for (Integer clusterNumber: nodeListForGivenClusterNumber.keySet()) {
			String groupName = this.attributeNameCluster+"_"+clusterNumber;
			List<CyNode>nodeList = nodeListForGivenClusterNumber.get(clusterNumber);
			// Create the group
			CyGroup newgroup = CyGroupManager.createGroup(groupName, nodeList, null);
			if (newgroup != null) {
				// Now tell the metanode viewer about it
				CyGroupManager.setGroupViewer(newgroup, "metaNode", Cytoscape.getCurrentNetworkView(), true);
				groupList.add(groupName);
			}
		}

		// Save the network attribute so we remember which groups are ours
		netAttributes.setListAttribute(networkID, this.groupAttribute, groupList);
	}

	private void assignNodeToCluster(Hashtable<String, Boolean> already, 
	                                 CyNode node, String id, int clusterNr, 
	                                 Vector<CyNode> nodesInThisCluster) {
		
		already.put(id, true);
		
		cyNodeAttributes.setAttribute(id, this.attributeNameCC, clusterNr);
		nodesInThisCluster.add(node);
		
		for (CyNode node2: nodes) {
			String id2 = node2.getIdentifier();
			
			if (!id.equalsIgnoreCase("Source") && 
			    !id.equalsIgnoreCase("Target") && 
			    !id2.equalsIgnoreCase("Source") && 
			    !id2.equalsIgnoreCase("Target")) {
				
				if (!already.containsKey(id2)) {
					String key = id + "#" + id2;
					String keyI = id2 + "#" + id;
					
					double sim = this.normalizedThreshold -1;
					if (this.normalizedSimilaritiesForGivenEdges.containsKey(key)) {
						sim = this.normalizedSimilaritiesForGivenEdges.get(key);
					} else if (this.normalizedSimilaritiesForGivenEdges.containsKey(keyI)) {
						sim = this.normalizedSimilaritiesForGivenEdges.get(keyI);
					}
					// System.out.println("sim: " + sim);
					if (sim >= this.normalizedThreshold) {
						assignNodeToCluster(already, node2, id2, clusterNr, nodesInThisCluster);
					}
				}
			}
		}
	}
}
