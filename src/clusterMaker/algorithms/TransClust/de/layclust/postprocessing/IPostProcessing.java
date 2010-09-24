/*
 *  Created on 3. December 2007
 */
package clusterMaker.algorithms.TransClust.de.layclust.postprocessing;

import clusterMaker.algorithms.TransClust.de.layclust.datastructure.ConnectedComponent;


/**
 * @author sita
 *
 */
public interface IPostProcessing {
	
	public void initPostProcessing(ConnectedComponent cc);

	public void run();
	
}
