package clusterMaker.algorithms.attributeClusterers.hopach.types;

public enum SplitCost {
	//MEAN_SPLIT_SILHOUETTE("Mean split silhouette"),
	//MEDIAN_SPLIT_SILHOUETTE("Median split silhouette"),
	AVERAGE_SPLIT_SILHOUETTE("Average split silhouette"),
	//MEAN_SILHOUETTE("Mean silhouette"),
	//MEDIAN_SILHOUETTE("Median silhouette"),
	AVERAGE_SILHOUETTE("Average silhouette");

	private String name;

	SplitCost(String name) {
		this.name = name;
	}

	public String toString() {
		return this.name;
	}

}
