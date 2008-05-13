/* vim: set ts=2: */
/**
 * Copyright (c) 2008 The Regents of the University of California.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions, and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above
 *      copyright notice, this list of conditions, and the following
 *      disclaimer in the documentation and/or other materials provided
 *      with the distribution.
 *   3. Redistributions must acknowledge that this software was
 *      originally developed by the UCSF Computer Graphics Laboratory
 *      under support by the NIH National Center for Research Resources,
 *      grant P41-RR01081.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
 * OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
package clusterMaker.algorithms.hierarchical;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.awt.GridLayout;
import javax.swing.JPanel;

// Cytoscape imports
import cytoscape.Cytoscape;
import cytoscape.data.CyAttributes;
import cytoscape.layout.Tunable;
import cytoscape.task.TaskMonitor;
import clusterMaker.algorithms.AbstractClusterAlgorithm;

// clusterMaker imports

public class HierarchicalCluster extends AbstractClusterAlgorithm {
	/**
	 * Linkage types
	 */
	ClusterMethod[] linkageTypes = { ClusterMethod.AVERAGE_LINKAGE,
	                                 ClusterMethod.SINGLE_LINKAGE,
	                                 ClusterMethod.MAXIMUM_LINKAGE,
	                                 ClusterMethod.CENTROID_LINKAGE };

	DistanceMetric[] distanceTypes = { DistanceMetric.EUCLIDEAN,
	                                   DistanceMetric.CITYBLOCK,
	                                   DistanceMetric.CORRELATION,
	                                   DistanceMetric.ABS_CORRELATION,
	                                   DistanceMetric.UNCENTERED_CORRELATION,
	                                   DistanceMetric.ABS_UNCENTERED_CORRELATION,
	                                   DistanceMetric.SPEARMANS_RANK,
	                                   DistanceMetric.KENDALLS_TAU };
	String[] attributeArray = new String[1];

	ClusterMethod clusterMethod =  ClusterMethod.AVERAGE_LINKAGE;
	DistanceMetric distanceMetric = DistanceMetric.EUCLIDEAN;
	boolean transposeMatrix = false;
	String dataAttributes = null;
	TaskMonitor monitor = null;

	public HierarchicalCluster() {
		super();
		initializeProperties();
	}

	public String getShortName() {return "hierarchical";};
	public String getName() {return "Hierarchical cluster";};

	public JPanel getSettingsPanel() {
		// Everytime we ask for the panel, we want to update our attributes
		Tunable attributeTunable = clusterProperties.get("attributeList");
		attributeArray = getAllAttributes();
		attributeTunable.setLowerBound((Object)attributeArray);

		return clusterProperties.getTunablePanel();
	}

	protected void initializeProperties() {
		super.initializeProperties();

		/**
		 * Tuning values
		 */

		// The linkage to use
		clusterProperties.add(new Tunable("linkage",
		                                  "Linkage",
		                                  Tunable.LIST, new Integer(0),
		                                  (Object)linkageTypes, (Object)null, 0));

		// The distance metric to use
		clusterProperties.add(new Tunable("dMetric",
		                                  "Distance Metric",
		                                  Tunable.LIST, new Integer(0),
		                                  (Object)distanceTypes, (Object)null, 0));

		clusterProperties.add(new Tunable("attributeListGroup",
		                                  "Source for array data",
		                                  Tunable.GROUP, new Integer(1)));

		// The attribute to use to get the weights
		attributeArray = getAllAttributes();
		clusterProperties.add(new Tunable("attributeList",
		                                  "Array sources",
		                                  Tunable.LIST, "",
		                                  (Object)attributeArray, (Object)null, Tunable.MULTISELECT));

		// Whether or not to transpose the matrix
		clusterProperties.add(new Tunable("transposeMatrix",
		                                  "Transpose matrix (if using node attributes)",
		                                  Tunable.BOOLEAN, new Boolean(false)));

		clusterProperties.initializeProperties();
		updateSettings(true);
	}

	public void updateSettings() {
		updateSettings(false);
	}

	public void updateSettings(boolean force) {
		clusterProperties.updateValues();
		super.updateSettings(force);

		Tunable t = clusterProperties.get("linkage");
		if ((t != null) && (t.valueChanged() || force))
			clusterMethod = linkageTypes[((Integer) t.getValue()).intValue()];

		t = clusterProperties.get("dMetric");
		if ((t != null) && (t.valueChanged() || force))
			distanceMetric = distanceTypes[((Integer) t.getValue()).intValue()];

		t = clusterProperties.get("transposeMatrix");
		if ((t != null) && (t.valueChanged() || force))
			transposeMatrix = ((Boolean) t.getValue()).booleanValue();

		t = clusterProperties.get("attributeList");
		if ((t != null) && (t.valueChanged() || force)) {
			dataAttributes = (String) t.getValue();
		}
	}

	public void doCluster(TaskMonitor monitor) {
		this.monitor = monitor;
		// Sanity check all of our settings
		System.out.println("Performing hierarchical cluster with method: "+clusterMethod+" using "+distanceMetric+" and attributes: "+dataAttributes);
		// Get our attributes we're going to use for the cluster
		String attributeArray[] = getAttributeArray(dataAttributes);
		// To make debugging easier, sort the attribute array
		Arrays.sort(attributeArray);

		// Cluster!
		EisenCluster.cluster(attributeArray, distanceMetric, clusterMethod, transposeMatrix);
	}

	private void getAttributesList(List<String>attributeList, CyAttributes attributes, 
	                              String prefix) {
		String[] names = attributes.getAttributeNames();
		for (int i = 0; i < names.length; i++) {
			if (attributes.getType(names[i]) == CyAttributes.TYPE_FLOATING ||
			    attributes.getType(names[i]) == CyAttributes.TYPE_INTEGER) {
				attributeList.add(prefix+names[i]);
			}
		}
	}

	private String[] getAllAttributes() {
		// Create the list by combining node and edge attributes into a single list
		List<String> attributeList = new ArrayList<String>();
		getAttributesList(attributeList, Cytoscape.getNodeAttributes(),"node.");
		getAttributesList(attributeList, Cytoscape.getEdgeAttributes(),"edge.");
		return attributeList.toArray(attributeArray);
	}

	private String[] getAttributeArray(String dataAttributes) {
		String indices[] = dataAttributes.split(",");
		String selectedAttributes[] = new String[indices.length];
		for (int i = 0; i < indices.length; i++) {
			selectedAttributes[i] = attributeArray[Integer.parseInt(indices[i])];
		}
		return selectedAttributes;
	}
}