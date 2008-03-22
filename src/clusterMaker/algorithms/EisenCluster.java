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
package clusterMaker.algorithms;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.awt.GridLayout;
import javax.swing.JPanel;

// Cytoscape imports
import cytoscape.Cytoscape;
import cytoscape.data.CyAttributes;
import cytoscape.layout.Tunable;
import cytoscape.task.TaskMonitor;

// clusterMaker imports

public class EisenCluster {
	final static int IS = 0;
	final static int JS = 1;

	public static String cluster(String weightAttributes[], DistanceMetric metric, 
	                      ClusterMethod clusterMethod, boolean transpose) {

		String keyword = "GENE";
		if (transpose) keyword = "ARRY";

		// Create the matrix
		Matrix matrix = new Matrix(weightAttributes, transpose);

		// Create a weight vector of all ones (we don't use individual weighting, yet)
		matrix.setUniformWeights();

		// Cluster
		TreeNode[] nodeList = treeCluster(matrix, metric, clusterMethod);
		if (nodeList == null || nodeList.length == 0) System.out.println("treeCluster returned empty tree!");

		if (metric == DistanceMetric.EUCLIDEAN || metric == DistanceMetric.CITYBLOCK) {
			// Normalize distances to between 0 and 1
			double scale = 0.0;
			for (int node = 0; node < nodeList.length; node++) {
				if (nodeList[node].getDistance() > scale) scale = nodeList[node].getDistance();
			}
			if (scale != 0.0) {
				for (int node = 0; node < nodeList.length; node++) {
					double dist = nodeList[node].getDistance();
					nodeList[node].setDistance(dist/scale);
				}
			}
		}

		// Join the nodes
		double[] nodeOrder = new double[nodeList.length];
		int[] nodeCounts = new int[nodeList.length];
		String[] nodeID = new String[nodeList.length];

		for (int node = 0; node < nodeList.length; node++) {
			int min1 = nodeList[node].getLeft();
			int min2 = nodeList[node].getRight();

			double order1;
			double order2;
			int counts1;
			int counts2;
			String ID1;
			String ID2;
			nodeID[node] = "GROUP"+(node+1)+"X";
			if (min1 < 0) {
				int index1 = -min1-1;
				order1 = nodeOrder[index1];
				counts1 = nodeCounts[index1];
				ID1 = nodeID[index1];
				nodeList[node].setDistance(Math.max(nodeList[node].getDistance(), nodeList[index1].getDistance()));
			} else {
				order1 = min1;
				counts1 = 1;
				// ID1 = keyword+min1+"X"; // Shouldn't this be the name of the gene/condition?
				ID1 = matrix.getRowLabel(min1);
			}

			if (min2 < 0) {
				int index2 = -min2-1;
				order2 = nodeOrder[index2];
				counts2 = nodeCounts[index2];
				ID2 = nodeID[index2];
				nodeList[node].setDistance(Math.max(nodeList[node].getDistance(), nodeList[index2].getDistance()));
			} else {
				order2 = min2;
				counts2 = 1;
				// ID2 = keyword+min2+"X"; // Shouldn't this be the name of the gene/condition?
				ID2 = matrix.getRowLabel(min2);
			}

			System.out.println(nodeID[node]+"\t"+ID1+"\t"+ID2+"\t"+(1.0-nodeList[node].getDistance()));

			nodeCounts[node] = counts1 + counts2;
			nodeOrder[node] = (counts1*order1 + counts2*order2) / (counts1 + counts2);
		}

		// Now sort based on tree structure
		TreeSort(nodeList.length, nodeOrder, nodeCounts, nodeList);

		// Update the network attribute "HierarchicalCluster" and make it hidden
		return "Complete";
	}

	private static TreeNode[] treeCluster(Matrix matrix, DistanceMetric metric, ClusterMethod clusterMethod) { 

		matrix.printMatrix();
		double[][] distanceMatrix = matrix.getDistanceMatrix(metric);
		TreeNode[] result = null;

		switch (clusterMethod) {
			case SINGLE_LINKAGE:
				System.out.println("Calculating single linkage hierarchical cluster");
				result = pslCluster(matrix, distanceMatrix, metric);
				break;

			case MAXIMUM_LINKAGE:
				System.out.println("Calculating maximum linkage hierarchical cluster");
				result = pmlcluster(matrix.nRows(), distanceMatrix);
				break;

			case AVERAGE_LINKAGE:
				System.out.println("Calculating average linkage hierarchical cluster");
				result = palcluster(matrix.nRows(), distanceMatrix);
				break;

			case CENTROID_LINKAGE:
				System.out.println("Calculating centroid linkage hierarchical cluster");
				result = pclcluster(matrix, distanceMatrix, metric);
				break;
		}
		return result;
	}

	/**
 	 * The pslcluster routine performs single-linkage hierarchical clustering, using
 	 * either the distance matrix directly, if available, or by calculating the
 	 * distances from the data array. This implementation is based on the SLINK
 	 * algorithm, described in:
 	 * Sibson, R. (1973). SLINK: An optimally efficient algorithm for the single-link
 	 * cluster method. The Computer Journal, 16(1): 30-34.
 	 * The output of this algorithm is identical to conventional single-linkage
 	 * hierarchical clustering, but is much more memory-efficient and faster. Hence,
 	 * it can be applied to large data sets, for which the conventional single-
 	 * linkage algorithm fails due to lack of memory.
 	 *
 	 * @param matrix the data matrix containing the data and labels
 	 * @param distanceMatrix the distances that will be used to actually do the clustering.
 	 * @param metric the distance metric to be used.
 	 * @return the array of TreeNode's that describe the hierarchical clustering solution, or null if
 	 * it it files for some reason.
 	 **/

	private static TreeNode[] pslCluster(Matrix matrix, double[][] distanceMatrix, DistanceMetric metric) {
		int nRows = matrix.nRows();
		int nNodes = nRows-1;

		int[] vector = new int[nNodes];
		TreeNode[] nodeList = new TreeNode[nNodes]; 
		// Initialize
		for (int i = 0; i < nNodes; i++) {
			vector[i] = i;
			nodeList[i] = new TreeNode(Double.MAX_VALUE);
		}

		int k = 0;
		double[] temp = new double[nNodes];

		for (int row = 0; row < nRows; row++) {
			if (distanceMatrix != null) {
				for (int j = 0; j < row; j++) temp[j] = distanceMatrix[row][j];
			} else {
				for (int j = 0; j < row; j++)
					temp[j] = metric.getMetric(matrix, matrix, matrix.getWeights(), row, j);
			}
			for (int j = 0; j < row; j++) {
				k = vector[j];
				if (nodeList[j].getDistance() >= temp[j]) {
					if (nodeList[j].getDistance() < temp[k]) temp[k] = nodeList[j].getDistance();
					nodeList[j].setDistance(temp[j]);
					vector[j] = row;
				} else if (temp[j] < temp[k]) temp[k] = temp[j];
			}
			for (int j = 0; j < row; j++) {
				if (vector[j] == nNodes || nodeList[j].getDistance() >= nodeList[vector[j]].getDistance()) vector[j] = row;
			}
		}


		for (int row = 0; row < nNodes; row++)
			nodeList[row].setLeft(row);

		Arrays.sort(nodeList, new NodeComparator());

		int[] index = new int[nRows];
		for (int i = 0; i < nRows; i++) index[i] = i;
		for (int i = 0; i < nNodes; i++) {
			int j = nodeList[i].getLeft();
			k = vector[j];
			nodeList[i].setLeft(index[j]);
			nodeList[i].setRight(index[k]);
			index[k] = -i-1;
		}

		return nodeList;
	}

	/**
 	 * The pclcluster routine performs clustering, using pairwise centroid-linking
 	 * on a given set of gene expression data, using the distrance metric given by metric.
 	 *
 	 * @param matrix the data matrix containing the data and labels
 	 * @param distanceMatrix the distances that will be used to actually do the clustering.
 	 * @param metric the distance metric to be used.
 	 * @return the array of TreeNode's that describe the hierarchical clustering solution, or null if
 	 * it it files for some reason.
 	 **/
	private static TreeNode[] pclcluster(Matrix matrix, double[][] distanceMatrix, DistanceMetric metric) {
		int nRows = matrix.nRows();
		int nColumns = matrix.nColumns();
		int nNodes = nRows-1;
		double mask[][] = new double[matrix.nRows()][matrix.nColumns()];

		TreeNode[] nodeList = new TreeNode[nNodes]; 

		// Initialize
		Matrix newData = new Matrix(matrix);
		// System.out.println("New matrix: ");
		// newData.printMatrix();

		int distID[] = new int[nRows];
		for (int row = 0; row < nRows; row++) {
			distID[row] = row;
			for (int col = 0; col < nColumns; col++) {
				if (newData.hasValue(row, col))
					mask[row][col] = 1.0;
				else
					mask[row][col] = 0.0;
			}
			if (row < nNodes)
				nodeList[row] = new TreeNode(Double.MAX_VALUE);
		}

		int pair[] = new int[2];

		for (int inode = 0; inode < nNodes; inode++) {
			// find the pair with the shortest distance
			pair[IS] = 1; pair[JS] = 0;
			double distance = findClosestPair(nRows-inode, distanceMatrix, pair);
			nodeList[inode].setDistance(distance);

			int is = pair[IS];
			int js = pair[JS];
			nodeList[inode].setLeft(distID[js]);
			nodeList[inode].setRight(distID[is]);
	
			// make node js the new node
			for (int col = 0; col < nColumns; col++) {
				double jsValue = newData.doubleValue(js, col);
				double isValue = newData.doubleValue(is, col);
				double newValue = 0.0;
				if (newData.hasValue(js,col)) newValue = jsValue * mask[js][col];
				if (newData.hasValue(is,col)) newValue += isValue * mask[is][col];

				if (newData.hasValue(js,col) || newData.hasValue(is,col)) {
					newData.setValue(js, col, newValue);
				}
				mask[js][col] += mask[is][col];
				if (mask[js][col] != 0) {
					newData.setValue(js, col, newValue / mask[js][col]);
				}
			}

			for (int col = 0; col < nColumns; col++) {
				mask[is][col] = mask[nNodes-inode][col];
				newData.setValue(is, col, newData.getValue(nNodes-inode, col));
			}

			// Fix the distances
			distID[is] = distID[nNodes-inode];
			for (int i = 0; i < is; i++) {
				distanceMatrix[is][i] = distanceMatrix[nNodes-inode][i];
			}

			for (int i = is+1; i < nNodes-inode; i++) {
				distanceMatrix[i][is] = distanceMatrix[nNodes-inode][i];
			}

			distID[js] = -inode-1;
			for (int i = 0; i < js; i++) {
				distanceMatrix[js][i] = metric.getMetric(newData, newData, newData.getWeights(), js, i);
			}
			for (int i = js+1; i < nNodes-inode; i++) {
				distanceMatrix[i][js] = metric.getMetric(newData, newData, newData.getWeights(), js, i);
			}
		}

		return nodeList;
	}

	/**
	 * The pmlcluster routine performs clustering using pairwise maximum- (complete-)
	 * linking on the given distance matrix.
	 * 
	 * @param nRows The number of rows to be clustered
	 * @param distanceMatrix The distance matrix, with rows rows, each row being filled up to the
	 * diagonal. The elements on the diagonal are not used, as they are assumed to be
	 * zero. The distance matrix will be modified by this routine.
	 * @return the array of TreeNode's that describe the hierarchical clustering solution, or null if
	 * it fails for some reason.
	 */
	private static TreeNode[] pmlcluster(int nRows, double[][] distanceMatrix) {
		int[] clusterID = new int[nRows];
		TreeNode[] nodeList = new TreeNode[nRows-1]; 

		for (int j = 0; j < nRows; j++) {
			clusterID[j] = j;
		}

		int pair[] = new int[2];
		for (int n = nRows; n > 1; n--) {
			pair[0] = 1; pair[1] = 2;
			if (nodeList[nRows-n] == null)
				nodeList[nRows-n] = new TreeNode(Double.MAX_VALUE);
			nodeList[nRows-n].setDistance(findClosestPair(n, distanceMatrix, pair));
			int is = pair[0];
			int js = pair[1];

			// Fix the distances
			for (int j = 0; j < js; j++)
				distanceMatrix[js][j] = Math.max(distanceMatrix[is][j],distanceMatrix[js][j]);
			for (int j = js+1; j < is; j++)
				distanceMatrix[j][js] = Math.max(distanceMatrix[is][j],distanceMatrix[j][js]);
			for (int j = is+1; j < n; j++)
				distanceMatrix[j][js] = Math.max(distanceMatrix[j][is],distanceMatrix[j][js]);
			for (int j = 0; j < is; j++)
				distanceMatrix[is][j] = distanceMatrix[n-1][j];
			for (int j = is+1; j < n-1; j++)
				distanceMatrix[j][is] = distanceMatrix[n-1][j];

			// Update cluster IDs
			nodeList[nRows-n].setLeft(clusterID[is]);
			nodeList[nRows-n].setRight(clusterID[js]);
			clusterID[js] = n-nRows-1;
			clusterID[is] = clusterID[n-1];
		}
		return nodeList;
	}

	/**
	 * The pmlcluster routine performs clustering using pairwise average
	 * linking on the given distance matrix.
	 * 
	 * @param nRows The number of rows to be clustered
	 * @param distanceMatrix The distance matrix, with rows rows, each row being filled up to the
	 * diagonal. The elements on the diagonal are not used, as they are assumed to be
	 * zero. The distance matrix will be modified by this routine.
	 * @return the array of TreeNode's that describe the hierarchical clustering solution, or null if
	 * it fails for some reason.
	 */
	private static TreeNode[] palcluster(int nRows, double[][] distanceMatrix) {
		int[] clusterID = new int[nRows];
		int[] number = new int[nRows];
		TreeNode[] nodeList = new TreeNode[nRows-1]; 

		// Setup a list specifying to which cluster a gene belongs, and keep track
		// of the number of elements in each cluster (needed to calculate the
		// average).
		for (int j = 0; j < nRows; j++) {
			number[j] = 1;
			clusterID[j] = j;
		}

		int pair[] = new int[2];
		for (int n = nRows; n > 1; n--) {
			int sum = 0;
			pair[IS] = 1; pair[JS] = 0;
			if (nodeList[nRows-n] == null)
				nodeList[nRows-n] = new TreeNode(Double.MAX_VALUE);
			double distance = findClosestPair(n, distanceMatrix, pair);
			nodeList[nRows-n].setDistance(distance);

			// Save result
			int is = pair[IS];
			int js = pair[JS];
			nodeList[nRows-n].setLeft(clusterID[is]);
			nodeList[nRows-n].setRight(clusterID[js]);

			// Fix the distances
			sum = number[is] + number[js];
			for (int j = 0; j < js; j++) {
				distanceMatrix[js][j] = (distanceMatrix[is][j]*(double)number[is] + distanceMatrix[js][j]*(double)number[js])/(double)sum;
			}

			for (int j = js+1; j < is; j++) {
				distanceMatrix[j][js] = (distanceMatrix[is][j]*(double)number[is] + distanceMatrix[j][js]*(double)number[js])/(double)sum;
			}

			for (int j = is+1; j < n; j++) {
				distanceMatrix[j][js] = (distanceMatrix[j][is]*(double)number[is] + distanceMatrix[j][js]*(double)number[js])/(double)sum;
			}

			for (int j = 0; j < is; j++) {
				distanceMatrix[is][j] = distanceMatrix[n-1][j];
			}
			for (int j = is+1; j < n-1; j++) {
				distanceMatrix[j][is] = distanceMatrix[n-1][j];
			}

			// Update number of elements in the clusters
			number[js] = sum;
			number[is] = number[n-1];

			// Update cluster IDs
			clusterID[js] = n-nRows-1;
			clusterID[is] = clusterID[n-1];
		}
		return nodeList;
	}

	private static void getAttributesList(List<String>attributeList, CyAttributes attributes, 
	                              String prefix) {
		String[] names = attributes.getAttributeNames();
		for (int i = 0; i < names.length; i++) {
			if (attributes.getType(names[i]) == CyAttributes.TYPE_FLOATING ||
			    attributes.getType(names[i]) == CyAttributes.TYPE_INTEGER) {
				attributeList.add(prefix+names[i]);
			}
		}
	}

	private static String[] getAllAttributes() {
		// Create the list by combining node and edge attributes into a single list
		List<String> attributeList = new ArrayList<String>();
		attributeList.add("-- select attribute --");
		getAttributesList(attributeList, Cytoscape.getNodeAttributes(),"node.");
		getAttributesList(attributeList, Cytoscape.getEdgeAttributes(),"edge.");
		return (String[])attributeList.toArray();
	}
		
	/**
 	 * This function searches the distance matrix to find the pair with the shortest
 	 * distance between them. The indices of the pair are returned in ip and jp; the
 	 * distance itself is returned by the function.
 	 *
 	 * n          (input) int
 	 * The number of elements in the distance matrix.
 	 *
 	 * distanceMatrix (input) double[][]
 	 * A ragged array containing the distance matrix. The number of columns in each
 	 * row is one less than the row index.
 	 *
 	 * pair         (output) int[2]
 	 * An array with two values representing the first and second indices of the pair
 	 * with the shortest distance.
 	 */
	private static double findClosestPair(int n, double[][] distanceMatrix, int[] pair) {
		int ip = 1;
		int jp = 0;
		double temp;
		double distance = distanceMatrix[1][0];
		for (int i = 1; i < n; i++) {
			for (int j = 0; j < i; j++) {
				temp = distanceMatrix[i][j];
				if (temp < distance) {
					distance = temp;
					ip = i;
					jp = j;
				}
			}
		}
		pair[IS] = ip;
		pair[JS] = jp;
		return distance;
	}

	private static void TreeSort(int nNodes, double nodeOrder[], int nodeCounts[], TreeNode nodeList[]) {
		int nElements = nNodes+1;
		double newOrder[] = new double[nElements];
		int clusterID[] = new int[nElements];
		
		for (int i = 0; i < nElements; i++) clusterID[i] = i;

		for (int i = 0; i < nNodes; i++) {
			int i1 = nodeList[i].getLeft();
			int i2 = nodeList[i].getRight();
		}
	}
}
