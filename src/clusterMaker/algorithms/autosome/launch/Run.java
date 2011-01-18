/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package clusterMaker.algorithms.autosome.launch;

import clusterMaker.algorithms.autosome.cluststruct.*;
import java.io.*;
import java.util.*;
import clusterMaker.algorithms.autosome.mapping.som.*;
import clusterMaker.algorithms.autosome.mapping.cartogram.*;
import clusterMaker.algorithms.autosome.mapping.sammonmapping.*;
import clusterMaker.algorithms.autosome.clustering.*;
import clusterMaker.algorithms.autosome.clustering.agglomerative.*;
import clusterMaker.algorithms.autosome.clustering.kmeans.*;
import clusterMaker.algorithms.autosome.clustering.mst.*;
//import view.view3d.*;
import java.text.*;
import clusterMaker.algorithms.autosome.view.view2d.viewer2D;
import javax.swing.*;
import cytoscape.task.TaskMonitor;
//import org.jvnet.substance.*;

/**
 *  launch AutoSOME clustering suite
 * @author Aaron
 */
public class Run {

    private static boolean substance = true;
    private boolean fillMissing = false;
    public boolean openingFile=false;
    private static javax.swing.JFrame splash;
    private boolean negativeLog=false;
    private launchSplash ls;
    private TaskMonitor monitor;


    public Run(boolean splash){
        if(splash){
            ls = new launchSplash();
            Thread t = new Thread(ls);
            t.start();
            try{
                t.join();
            }catch(Exception err){};
        }
    }

    public Run(){};
    
    public static void main(String[] args){
         //s.printDefaultSettings();
            //System.exit(0);
            try {
            // Set System L&F
                UIManager.setLookAndFeel(
                    UIManager.getSystemLookAndFeelClassName());
            }
            catch (UnsupportedLookAndFeelException e) {
            // handle exception
            }
            catch (ClassNotFoundException e) {
            // handle exception
            }
            catch (InstantiationException e) {
            // handle exception
            }
            catch (IllegalAccessException e) {
            // handle exception
            }
javax.swing.JFrame.setDefaultLookAndFeelDecorated(true);

Run r = new Run(true);

        printBanner();    
        Settings s = new Settings();
 
        if(args.length == 0) {
           java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                try{
                runGUI = true;
                viewer2D v2d = new viewer2D();
                javax.swing.ImageIcon icon = new javax.swing.ImageIcon(getClass().getResource("/imgs/icon.jpg"));
                v2d.setIconImage(icon.getImage());
                v2d.setVisible(true);
                }catch(Exception err){};
            }
           });

        }else{
            s.setParams(args); //set user parameters in 's'

                if(s.invokeViewerOnly) {
                    viewer2D v2d = new viewer2D();
                    v2d.setVisible(true);
                    v2d.invokeClusterViewer_From_trIDFile(new File(s.inputFile),s);
                }

                new Run().run(s, true, new viewer2D());

        
        }
    }




    public clusterRun run(Settings s, boolean readInput, viewer2D jpb){
        
        long t = System.currentTimeMillis();
        //this.jpb = jpb; //set progress bar

        if(readInput){
            jpb.jLabel7.setText("Preprocessing");
            jpb.jProgressBar1.setIndeterminate(true);
            s.input = getInput(s); //read input file and load into dataItem struct
            if(s.input==null) return null;
        }else{
            for(int i = 0; i < s.input.length; i++){
                s.input[i].setValue(i, s.input[i].getOriginalValues()[i]);
            }
        }
            //for(int k = 0; k < s.input.length; k++) System.out.println(s.input[k]);
             if(s.distMatrix) {

                if(runGUI) {
                    jpb.jLabel7.setText("Distance Matrix");
                    jpb.jProgressBar1.setIndeterminate(true);
                }

                if((s.unitVar || s.scale>0 || s.logNorm)) {
                    if(s.unitVar && s.logNorm){
                        s.unitVar=false;
                        s.input = norm(s); //normalize input
                        s.unitVar=true; s.logNorm=false;
                        s.input = norm(s);
                        s.logNorm=true;
                    }else   s.input = norm(s);
                }
                if(s.medCenter||s.medCenterCol) s.input = medCenter(s);
                if(s.sumSqrRows||s.sumSqrCol) s.input = getSumSqr(s);
                if(negativeLog && s.logNorm) return null;
                //if(s.medCenter) s.input = medCenter(s);
                s.input = new makeDistMatrix().getDistMatrix(s);
                if(s.benchmark) s.convertDMLabels();

                if(runGUI) {
                    jpb.jLabel7.setText("Clustering");
                    jpb.jProgressBar1.setIndeterminate(false);
                }
             }
       
            s.inputSize = s.input.length;
       // }
        if(!s.distMatrix){
            jpb.jLabel7.setText("Preprocessing");
            jpb.jProgressBar1.setIndeterminate(true);
            if((s.unitVar || s.scale>0 || s.logNorm)) {
                if(s.unitVar && s.logNorm){
                    s.unitVar=false;
                    s.input = norm(s); //normalize input
                    s.unitVar=true; s.logNorm=false;
                    s.input = norm(s);
                    s.logNorm=true;
            }else   s.input = norm(s);
        }
                if(s.medCenter||s.medCenterCol) s.input = medCenter(s);
                if(s.sumSqrRows||s.sumSqrCol) s.input = getSumSqr(s);
        
        }else if(s.unitVarAfterDM) s.input = norm(s);

        

        if(runGUI) {
             jpb.jProgressBar1.setIndeterminate(false);
             jpb.jLabel7.setText("Clustering");
             jpb.jProgressBar1.setIndeterminate(false);
        }
        
        if(negativeLog && s.logNorm) return null;
       
        s.setInputMinMax(); //store minimum and maximum values in input

        if(s.Pearson) s.setCenter(); //compute average of each attribute
        if(s.writeTemp) {
            new File(s.outputDirectory+s.getFolderDivider()+s.getName()+"_temp").mkdir();
        }

        
        invokeAutoSOME(s);    
        s.runTime = getrunTime(t);
       // System.out.println(">Running time: "+s.runTime);
        try{
             return doOutput(s, clusterRuns);
        }catch(Exception err){};

        return new clusterRun();
        
    }


    public clusterRun runAutoSOMEBasic(Settings s, TaskMonitor monitor){
            long t = System.currentTimeMillis();
            this.monitor=monitor;

            if(s.fillMissing){
                s.input=addMissing(s.input, s);
            }

             if(s.distMatrix) {


                if((s.unitVar || s.scale>0 || s.logNorm)) {
                    if(s.unitVar && s.logNorm){
                        s.unitVar=false;
                        s.input = norm(s); //normalize input
                        s.unitVar=true; s.logNorm=false;
                        s.input = norm(s);
                        s.logNorm=true;
                    }else   s.input = norm(s);
                }
                if(s.medCenter||s.medCenterCol) s.input = medCenter(s);
                if(s.sumSqrRows||s.sumSqrCol) s.input = getSumSqr(s);
                if(negativeLog && s.logNorm) return null;
                //if(s.medCenter) s.input = medCenter(s);
                s.input = new makeDistMatrix().getDistMatrix(s);
                if(s.benchmark) s.convertDMLabels();


             }

            s.inputSize = s.input.length;
       // }
        if(!s.distMatrix){

            if((s.unitVar || s.scale>0 || s.logNorm)) {
                if(s.unitVar && s.logNorm){
                    s.unitVar=false;
                    s.input = norm(s); //normalize input
                    s.unitVar=true; s.logNorm=false;
                    s.input = norm(s);
                    s.logNorm=true;
                }else   s.input = norm(s);
            }
            if(s.medCenter||s.medCenterCol) s.input = medCenter(s);
            if(s.sumSqrRows||s.sumSqrCol) s.input = getSumSqr(s);

        }else if(s.unitVarAfterDM) s.input = norm(s);


        if(negativeLog && s.logNorm) return null;

        s.setInputMinMax(); //store minimum and maximum values in input

        if(s.Pearson) s.setCenter(); //compute average of each attribute

        invokeAutoSOME(s);
        s.runTime = getrunTime(t);
      //  System.out.println(">Running time: "+s.runTime);
        try{
             return doOutput(s, clusterRuns);
        }catch(Exception err){};

        return new clusterRun();
    }
    
    
    
    public dataItem[] getInput(Settings s){
        
        dataItem[] data = new dataItem[1]; //store spreadsheet data in dataItem
        int maxNum = 0; // maximum cluster number counter
          int rowLength = 0; //count number of columns for each row to validate consistent input

        if(s.outputDirectory.equals("C:")) s.outputDirectory=new File(s.inputFile).getParent();
        //long t = System.currentTimeMillis();
        try{
            BufferedReader bf = new BufferedReader(new FileReader(s.inputFile));
            
            String line = new String();
            ArrayList a = new ArrayList();
            int lineCount = 0;
            int rowCount = 0;

            ///////additional input format////////
            String GEOsamples = new String();
            if(s.GEOformat){
                 while((line = bf.readLine())!= null){
                     if(line.length()==0) continue;
                     line = line.replaceAll("\"", "");
                     String[] tokens = line.split("\t");
                     if(tokens.length==0) continue;
                     if(tokens[0].equals("!Sample_title")) GEOsamples=line;
                     if(tokens[0].equals("ID_REF")){
                         if(GEOsamples.length()==0) GEOsamples=line;
                         break;
                     }
                 }
            }
            
      
            //////////////////////////////////////////////////


            while((line = bf.readLine())!= null){
                if(line.length()!=0){
                    lineCount++;
                    if(lineCount%1000==0) {
                        clusterMaker.algorithms.autosome.fileio.File_Open.jLabel15.setText(String.valueOf(lineCount));
                    }
                    if(s.GEOformat) {
                        if(lineCount==1) line=GEOsamples;
                        //fix this
                       // line = line.replaceAll("\"", "");
                        if(line.contains("\"")){
                        int pos = 0;
                        StringBuilder sb = new StringBuilder(line);
                        while((pos = sb.indexOf("\"", pos)) != -1) {
                            sb.insert(pos + 1, "");
                            pos += 2;
                        }
                        line = sb.toString();
                        }

                    }
                    //fix this
                   // line = line.replaceAll("\\t\\t", "\t?\t");
                    if(line.contains("\\t\\t")){
                        int pos = 0;
                        StringBuilder sb = new StringBuilder(line);
                        while((pos = sb.indexOf("\\t\\t", pos)) != -1) {
                            sb.insert(pos + 1, "\t?\t");
                            pos += 2;
                        }
                    line = sb.toString();
                    }              
                    String[] tokens = line.split("\t"); //try parsing by tab
                    

                    if(tokens.length == 1) {
                        line = line.replaceAll(",,", ",?,");
                        if(line.contains(",,")){
                        int pos = 0;
                        StringBuilder sb = new StringBuilder(line);
                        while((pos = sb.indexOf(",,", pos)) != -1) {
                            sb.insert(pos + 1, ",?,");
                            pos += 2;
                        }
                        line = sb.toString();
                        }
                        tokens = line.split(",");
                    } //else try by comma
                    if(tokens.length == 1) {
                        line = line.replaceAll("  ", " ? ");
                        if(line.contains("  ")){
                        int pos = 0;
                        StringBuilder sb = new StringBuilder(line);
                        while((pos = sb.indexOf("  ", pos)) != -1) {
                            sb.insert(pos + 1, " ? ");
                            pos += 2;
                        }
                        line = sb.toString();
                        }
                        tokens = line.split(" ");
                    } //else try by space


                    if(s.GEOformat) if(tokens[0].equals("!series_matrix_table_end")) continue;
                    if(s.PCLformat && tokens.length>0) if(tokens[0].equals("EWEIGHT")) {
                        s.EWEIGHT = new double[tokens.length-s.startData];
                        for(int i = s.startData; i < tokens.length; i++) {
                            s.EWEIGHT[i-s.startData] = Double.valueOf(tokens[i]);
                        }
                        continue;
                    }
                    if(tokens.length == 1 && !s.GEOformat) {
                     //   System.err.println("Error: Input file cannot be parsed correctly. Only tab-,comma-,space-delimited files can be read.");
                        if(!runGUI) System.exit(1);
                        else{
                            clusterMaker.algorithms.autosome.fileio.File_Open.ofile.kill();
                           
                            
                            return null;
                        }

                        }
                    //autocheck for column header
                    if(lineCount==1){
                        try{
                            Double.parseDouble(tokens[s.startData]);
                        }catch(NumberFormatException e){
                            s.readColumns=true;
                        }
                    }
                    ///////////////////////////
                    if(s.readColumns && lineCount == 1){
                       
                        if(s.PCLformat){
                            if(tokens[2].equals("GWEIGHT")) s.startData=3;
                            else s.startData=2;
                        }
                        s.columnHeaders = new String[tokens.length];
                        for(int i = 0; i < tokens.length; i++) {
                            tokens[i] = tokens[i].replace(",",";");
                            s.columnHeaders[i] = tokens[i];
                        }
                        continue;
                    }
                    
                    String identifier = tokens[s.PCLformat ? 0 : (s.startData-1)];
                    //if identifier is number, find maximum number (might be number of clusters in dataset)
                    try{
                        int num = Integer.parseInt(identifier);
                        if(maxNum < num)
                            maxNum = num;
                    
                    }catch(NumberFormatException err){
                        
                    };
                    
                    //load dataItem struct with input file information   
                    float[] values = new float[tokens.length-(s.startData )];
                    if(s.readColumns){
                        if(values.length<s.columnHeaders.length-s.startData){
                            float[] temp = new float[s.columnHeaders.length-s.startData];
                            for(int i = 0; i < temp.length; i++){
                                if(i<values.length) temp[i] = values[i];
                                else temp[i] = -99999999;
                            }
                            values=temp;
                        }
                    }
                    if(rowLength==0) rowLength=values.length;
                    if(s.GEOformat && values.length==0) continue;
                    else if (values.length!= rowLength && ((s.GEOformat&&values.length!=1) || !s.GEOformat)){
                    
                        if(runGUI){
                            clusterMaker.algorithms.autosome.fileio.File_Open.ofile.kill();
                                       
                                        }
                            

                       //   System.out.println(">Improper format: One or more rows of different lengths detected.");
                          return null;
                    }
                    for(int i = s.startData; i < tokens.length; i++) {
                            if(tokens[i].equals("?") || tokens[i].equals("NA") || tokens[i].equals("") || tokens[i].equals("null")){
                                {
                                  if(!fillMissing){
                                    if(runGUI && openingFile){
                                       clusterMaker.algorithms.autosome.fileio.File_Open.ofile.kill();
                                        
                                        }
                                        
                                        fillMissing=true;
                                       // if(openingFile) System.out.println(">One or more missing values detected. Will replace with "+((s.mvMedian) ? "median" : "mean")+" value of each "+((s.mvCol) ? "column" : "row")+".");

                                  }
                                }
                            }
                        try{
                            values[i-s.startData] = Float.valueOf(((tokens[i].equals("?")) || tokens[i].equals("")) ? "-99999999" : tokens[i]);
                        }catch(NumberFormatException err){
    
                        if(!fillMissing){
                            if(runGUI && openingFile){
                                clusterMaker.algorithms.autosome.fileio.File_Open.ofile.kill();
                                
                            }// if(openingFile)System.out.println("Improper format: Non-numerical data detected: \""+tokens[i]+"\". Will treat non-numbers as missing values.");
                             
                            fillMissing=true;

                        }
                             values[i-s.startData] = Float.valueOf("-99999999");
                        };
                    }
               
                    identifier = identifier.replace(" ","_");
                    identifier = identifier.replace(",",";");
                    dataItem d = new dataItem(values, identifier);
  
                    if(s.startData > 1){
                        StringBuilder sb = new StringBuilder();
                        for(int i = 0; i < s.startData; i++)
                            sb.append(tokens[i]+"\t");
                        d.setDesc(sb.toString());
                    }
                    if(s.kept.size()>0){
                       // System.out.println(a.size()+" "+s.kept.size()+" "+s.kept.containsKey(a.size()));
                        if(!s.kept.containsKey(rowCount++)) continue;
                    }
                    //a.add(getBytes(d));
                   // Object[] o = new Object[]{values, identifier.getBytes()};
                    //System.out.println(identifier);
                    a.add(d);
                }
            }
            
            data = new dataItem[a.size()];
            for(int i = 0; i < data.length; i++) data[i] = (dataItem)a.get(i);//toObject((byte[])a.remove(0));

            if(fillMissing){
                data=addMissing(data, s);
            }

        }catch(IOException err){//System.err.println("Error Reading File: "+s.inputFile);
            if(!runGUI) System.exit(1);
            else{
                           if(clusterMaker.algorithms.autosome.fileio.File_Open.ofile!=null) clusterMaker.algorithms.autosome.fileio.File_Open.ofile.kill();
                           
                          
                           return null;

            }
        }
        
        if(maxNum > 0 && (s.batch || s.known_clusters == 0)) s.known_clusters = maxNum;
        
        return data;
    }

     


    //allocate threads to ensemble process and start AutoSOME
    private void invokeAutoSOME(Settings s){


        if(!s.batch) {
           // System.out.println(">Input Size: "+s.inputSize+"\n>Attributes: "+s.input[0].getValues().length+"\n");
        }
        
        if(!s.doSM && s.verbose) {
            // if(s.som_gridSize == 0) System.out.println(">SOM grid size: "+((int)Math.min(s.som_maxGrid, Math.max(s.som_minGrid, Math.sqrt(s.input.length*2)))));
          //   else System.out.println(">SOM grid size: "+s.som_gridSize);
        }
                
        if(s.ensemble_runs>1) {
            clusterRuns = new ArrayList(); //store output of each cluster run
            // if(!s.batch) System.out.println(">Running Ensemble Clustering\n\n...computing clusters\n\n          |100%");
             monitor.setStatus("Clustering "+s.input.length+" rows by "+s.input[0].getValues().length+" columns");
        }

        long maxMemory = Runtime.getRuntime().maxMemory()/(1024*1024);

        if(!s.writeTemp && ((s.input.length>=100000) || ((s.input.length >= 50000 || (s.input.length * s.input[0].getValues().length >= 1200000)) && s.ensemble_runs>500 && maxMemory<=1600))){

            if(runGUI){
                                        
             }

          //  System.out.println(">You may need to write intermediate data to disk to avoid running out of memory (AutoSOME will hang if this happens).\nTo do this, go to Advanced Fields, then select the 'Write Ensemble Runs to Disk' checkbox.");

        }

        



        threads = new Thread[s.threads];
        
        for(int q = 0; q < threads.length; q++){
              threads[q] = new Thread(new runAutoSOME(s, q*s.ensemble_runs/(threads.length), (q+1)*s.ensemble_runs/(threads.length)));
              threads[q].start();
        }
        try{
             for(int q = 0; q < threads.length; q++){
                   threads[q].join();
             }
             
             if(progressCount < 10 && !s.batch && s.ensemble_runs>1){
                 for(; progressCount < 10; progressCount++) {
                     //System.out.print("*");
                     int currVal = 0;//jpb.jProgressBar1.getValue();
                    
                 }
             }
             //if(!s.batch) System.out.println("\n");
            
        }catch(Exception err){};
         if(s.ensemble_runs>1 && !s.batch) ensemble = (new Ensemble(clusterRuns, true, 0, (!s.doKmeans && !s.doHierarchical) ? true : false, s, monitor)).run();

    }
    
    
     //run an instance of AutoSOME in single thread
     public class runAutoSOME implements Runnable{
        Settings s;
        int start;
        int end;
        
        public runAutoSOME(Settings s, int start, int end) {
            this.s = s;
            this.start = start;
            this.end = end;
        }
        public void run(){
            
            for(;start<end;start++){
                try{
                    
                    //mapping to lower dimensional space
                    Object[] results = doMapping(s);

                    storeMapping.add(results);

                    float[][] coors = (float[][]) results[0];
                    ArrayList dataLabels = (ArrayList) results[1];
                    //clustering
                    clusterRun cr = doClustering(s, coors, dataLabels);
                    if(start>=10) cr.DEC = null;
                    if(!s.writeTemp) clusterRuns.add(cr);
                    else {
                        if(clusterRuns.size()==0) clusterRuns.add(cr);
                        else clusterRuns.add("1");
                        writeTemp(cr, clusterRuns.size(), s);
                    }

                    totalProgress++;
                   // System.out.println(((int)(100*(double)start/s.ensemble_runs))+" "+jpb.jProgressBar1.getValue());
                    //jpb.jProgressBar1.setValue((int)(100*((double)totalProgress/s.ensemble_runs)));
                    monitor.setPercentCompleted((int)(100*((double)totalProgress/s.ensemble_runs)));
                    if(!s.batch) if(s.ensemble_runs>1 && Math.floor((runCount++) %((double)(s.ensemble_runs)/10)) == 0) {
                    //    System.out.print("*");
                        
                        progressCount++;
                    }
                }catch(Exception err){};
            }
        }
    }
     
     
     
    private Object[] doMapping(Settings s){
        
                ArrayList dataLabels = new ArrayList();
                float[][] coors = new float[1][1];
                
                if(!s.doSM) {
                        //invoke Self-Organizing Map
                        SOM som  = new SOM(s);
                        som.run();
                        s.som_gridSize = som.getGridSize();
                        Object[] info = som.getDEInfo(3);
                        dataLabels = (ArrayList)info[2];
                        //invoke Density-equalizing Cartogram
                        if(s.doCart){
                            coors = (new DEC((ArrayList[])info[0], (float[])info[1], (float[])info[3], (ArrayList[])info[4],
                            true, false,
                            s.de_resolution, s.de_resolution)).makeCartogram();
                        }else{
                            ArrayList[] somCoors = (ArrayList[])info[0];
                            coors = new float[somCoors.length*4][3];
                            for(int j = 0, itor = 0; j < somCoors.length; j++){
                                    float[][] f = (float[][]) somCoors[j].get(0);
                                    for(int p = 0; p < f.length; p++){
                                        coors[itor][0] = f[p][0];
                                        coors[itor][1] = f[p][1];
                                        coors[itor++][2] = f[p][2];
                                        //System.out.println(f[p][0]+" "+f[p][1]);
                                    }
                            }
                        }
                    }  
                    else {
                        //invoke Sammon Mapping
                        sammonMapping sm = new sammonMapping(new javax.swing.JProgressBar(), s.sm_iters, 2);
                        coors = sm.run(s.input);
                        ArrayList ids = new ArrayList();
   
                        for(int j = 0; j < s.input.length; j++) {
                            ids.add(new int[]{j,j});
                            //System.out.println(s.input[j].getIdentity()+"\t"+coors[j][0]+"\t"+coors[j][1]);
                        }
                        dataLabels = ids;
                    }
                
        Object[] results = {coors,dataLabels};
        
        return results;
    }
    
    
    public clusterRun doClustering(Settings s, float[][] coors, ArrayList dataLabels){
        
        
                    MSTCluster mst = new MSTCluster((s.ensemble_runs>1) ? true : false);
                    mst.run(new javax.swing.JProgressBar(), coors, dataLabels, s.mst_pval,s.mst_MC, true, s);
                    
                    clusterRun cr = mst.getClusterRun();
                    cr.setInputFile(new String(s.inputFile));
                    
                    //System.out.println(cr.nodes.length);
                    if(s.doKmeans || s.doHierarchical){
                        
                        getClusters gc = new getClusters(cr,s);
                        gc.findClusters(true);
                        cr.c = gc.getClust();

                        if(s.doKmeans){
                            cr.c = new runKMeans().kdesom(cr, s.known_clusters, true, s);
                        }else{
                            if(s.hierarchical_choice == 1){
                                 cr.c = new agglomerative().runAgg_Mapping(cr, s.known_clusters, 1, true, s);
                            }
                            if(s.hierarchical_choice == 2){
                                 cr.c = new agglomerative().runAgg_Mapping(cr, s.known_clusters, 2, true, s);
                            }
                            if(s.hierarchical_choice == 3){
                                 cr.c = new agglomerative().runAgg_Mapping(cr, s.known_clusters, 3, true, s);
                            }
                            if(s.hierarchical_choice == 4){
                                 cr.c = new agglomerative().runAgg_Mapping(cr, s.known_clusters, 4, true, s);
                            }

                        }
                        
          
                    }
            return cr;
    }
     
     
    //output of clustering
    public clusterRun doOutput(Settings s, ArrayList clusterRuns){
        
        if(clusterRuns.isEmpty()) return null;
        clusterRun cr = (clusterRun) clusterRuns.get(0);
        if(clusterRuns.size() > 1) cr.c=ensemble.c;
        getClusters clust = new getClusters(cr,s);
        clusterMaker.algorithms.autosome.fileio.printClusters pc = new clusterMaker.algorithms.autosome.fileio.printClusters();
        
        if(clusterRuns.size() == 1 && !s.doKmeans && !s.doHierarchical){       
            clust.findClusters(true);
            cr.c = clust.getClust();
        }
      
       
        
        if(s.confidence && s.ensemble_runs > 1 && clusterRuns.size() > 1){
            sortCluster sc = new sortCluster();
            for(int i = 0; i < cr.c.length; i++) {
                cr.c[i] = sc.sortConf(cr.c[i],s);
                for(int j = 0; j < cr.c[i].ids.size(); j++){
                    int id = Integer.valueOf(cr.c[i].ids.get(j).toString());
                    s.input[id].setConf(Integer.valueOf(cr.c[i].confidence.get(j).toString()));
                }
            }
        }if(s.ensemble_runs == 1 || clusterRuns.size() == 1) s.confidence = false;
        

        Arrays.sort(cr.c);
 
        
        if(s.display2D){
           // System.out.println("...Launching Cluster Viewer");
            new clusterMaker.algorithms.autosome.view.view2d.viewer2D(cr.c,s).setVisible(true);
        }

        
        if(s.htmlOut || s.textOut) {
            s.add=s.getName()+"_E"+s.ensemble_runs+"_Pval"+s.mst_pval+((s.printRowsCols==1) ? "_rows" : ((s.printRowsCols==2) ? "_columns" : ""));
            pc.printClusters(s.add, cr.c, s);
        }
        
        if(s.benchmark || s.trBM) {
       
            clust.getClusterValidity(false);
        }
       // if(s.benchmark && !s.batch) System.out.println("\n>Benchmarking Results:\nF: "+cr.Fmeasure+"\tP: "+cr.Precision+"\tR: "+cr.Recall+"\tNMI: "+cr.NMI+"\tcR: "+cr.adjRand+"\n");
        
        
        return cr;
    }
    
     
    
     
    private static void printBanner(){
     /*  System.out.println("===============================================\n" +
                           "         AutoSOME Clustering Suite 2.0\n\n" +
                           "   Automatic clustering of density-equalized\n" +
                           "         Self-Organizing Map Ensembles\n\n" +
                           "        BMC Bioinformatics 2010, 11:117\n\n"+
                           "    University of California, Santa Barbara\n" +
                           "===============================================\n");
          */
    }
    
    
     //format running time and return time string
     public String getrunTime(long t){
        DecimalFormat Format = new DecimalFormat();
        Format.setMinimumFractionDigits(0);
        DecimalFormat Format2 = new DecimalFormat("#");
        Format2.setMinimumFractionDigits(2);
        
        double runT = System.currentTimeMillis()-t;
        String rt = (runT < 1000) ? String.valueOf(Format.format(runT)).concat(" msecs")
                  : (runT < 60000)? String.valueOf(Format2.format(runT/1000)).concat(" secs") :
                                    String.valueOf(Format2.format(runT/60000)).concat(" mins");

        return rt;
    }
     
    
   
    
    
    //median center normalization
    public dataItem[] medCenter(Settings s){
      if(s.medCenter){
        for(int i = 0; i < s.input.length; i++){
            float[] sort = new float[s.input[i].getValues().length];
            for(int j = 0; j < s.input[i].getValues().length; j++){
                sort[j] = s.input[i].getValues()[j];
            }
            Arrays.sort(sort); 
            float median = (sort.length%2!=0) ? sort[(int)Math.floor(sort.length/2)]
                                            : (sort[sort.length/2]
                                              +sort[(int)Math.floor(sort.length/2)])/2;
            for(int j = 0; j < s.input[i].getValues().length; j++){
                s.input[i].getValues()[j] -= median;                
            }
        }
      }
      if(s.medCenterCol){
        for(int i = 0; i < s.input[0].getValues().length; i++){
            float[] sort = new float[s.input.length];
            for(int j = 0; j < s.input.length; j++){
                sort[j] = s.input[j].getValues()[i];
            }
            Arrays.sort(sort);
            float median = (sort.length%2!=0) ? sort[(int)Math.floor(sort.length/2)]
                                            : (sort[sort.length/2]
                                              +sort[(int)Math.floor(sort.length/2)])/2;
            for(int j = 0; j < s.input.length; j++){
                s.input[j].getValues()[i] -= median;
            }
        }
      }
        return s.input;
    }

    //normalize rows/columns to 1
    public dataItem[] getSumSqr(Settings s){
    
        if(s.sumSqrRows){
          for(int i = 0; i < s.input.length; i++){
              float sumSqr=0;
              for(int j = 0; j < s.input[i].getValues().length; j++) sumSqr += (s.input[i].getValues()[j]*s.input[i].getValues()[j]);
              double sqrRoot = Math.sqrt(1/(double)sumSqr);
              if(Double.isNaN(sqrRoot)) sqrRoot=0;

              for(int j = 0; j < s.input[i].getValues().length; j++) {
                  s.input[i].getValues()[j] *= sqrRoot;
              }
              
          }
        }

        if(s.sumSqrCol){
            for(int i = 0; i < s.input[0].getValues().length; i++){
                float sumSqr=0;
                for(int j = 0; j < s.input.length; j++) sumSqr += (s.input[j].getValues()[i]*s.input[j].getValues()[i]);
                double sqrRoot = Math.sqrt(1/(double)sumSqr);
                if(Double.isNaN(sqrRoot)) sqrRoot=0;
                for(int j = 0; j < s.input.length; j++) s.input[j].getValues()[i] *= sqrRoot;
            }
        }

        return s.input;
    }
    
    
    //normalize input
    public dataItem[] norm(Settings s){
        //normalize input to unit variance, logarithm base 2, or range [0,s.scale]
        dataItem[] d = s.input;
        if(s.unitVar || s.scale > 0 || s.logNorm){
            for(int i = 0; i < d[0].getValues().length; i++){
                float sum = 0;    
                float max = Float.MIN_VALUE;
                float min = Float.MAX_VALUE;
                
                for(int j = 0; j < d.length; j++){
                    sum += d[j].getValues()[i];
                    if(d[j].getValues()[i] > max) max = d[j].getValues()[i];
                    if(d[j].getValues()[i] < min) min = d[j].getValues()[i];
                }
                
                float ave = sum / d.length;
                float x = 0;
                if(s.unitVar){
                    for(int j = 0; j < d.length; j++){
                        x += Math.pow(d[j].getValues()[i]-ave,2);
                    }
                }
                float stdev = (float)Math.sqrt(x/((float)d.length-1));
                if(stdev == 0) stdev = 1;
                
                for(int j = 0; j < d.length; j++){
                    if(s.unitVar) d[j].getValues()[i] = (d[j].getValues()[i]-ave)/stdev; //if unit var
                    else if(s.logNorm) {
                        if(d[j].getValues()[i]<=0){
                            /*if(runGUI){
                                        
                                        }

                            monitor.setStatus(">Improper number(s) detected. Cannot take logarithm of: "+d[j].getValues()[i]);
                            negativeLog=true;
                            return d;*/
                            d[j].getValues()[i]=0.00000001f;
                        }
                        d[j].getValues()[i] = (float)(Math.log10(d[j].getValues()[i]) / Math.log10(2)); //if log normalization
                        //if( d[j].getValues()[i] < 0)  d[j].getValues()[i]  = 0;
                    }
                    else d[j].getValues()[i] = s.scale*(d[j].getValues()[i] - min)/(max-min); //if scale
                }                
            }
        }
        
        return d;
    }
     




     private dataItem[] addMissing(dataItem[] d, Settings s){

         if(s.mvCol){
         for(int i = 0; i < d[0].getValues().length; i++){
            ArrayList notMissing = new ArrayList();
            for(int j = 0; j < d.length; j++){
                if(d[j].getValues()[i]!=-99999999) notMissing.add(d[j].getValues()[i]);
            }
            float[] sort = new float[notMissing.size()];
            float sum = 0;
            float median = 0;
            if(sort.length>0){
            for(int l = 0; l < sort.length; l++) {
                sort[l] = Float.valueOf(notMissing.get(l).toString());
                sum+=sort[l];
            }
            Arrays.sort(sort);
            median = (sort.length%2!=0) ? sort[(int)Math.floor(sort.length/2)]
                                            : (sort[sort.length/2]
                                              +sort[(int)Math.floor(sort.length/2)])/2;
            }
            for(int j = 0; j < d.length; j++){
                if(d[j].getValues()[i]==-99999999) {
                    d[j].getValues()[i] = (s.mvMedian) ? median : (sort.length>0) ? (sum/sort.length) : 0;
                   // if(d[j].getValues()[i]<-100) System.out.println(median+" "+sum+" "+sort.length);
                    d[j].setOriginalValue(i,d[j].getValues()[i]);
                }
            }
        }
         }else{
            for(int i = 0; i < d.length; i++){
            ArrayList notMissing = new ArrayList();
            for(int j = 0; j < d[i].getValues().length; j++){
                if(d[i].getValues()[j]!=-99999999) notMissing.add(d[i].getValues()[j]);
            }
            float[] sort = new float[notMissing.size()];
            float sum = 0;
            float median=0;
            if(sort.length>0){
            for(int l = 0; l < sort.length; l++) {
                sort[l] = Float.valueOf(notMissing.get(l).toString());
                sum+=sort[l];
            }
            Arrays.sort(sort);
            median = (sort.length%2!=0) ? sort[(int)Math.floor(sort.length/2)]
                                            : (sort[sort.length/2]
                                              +sort[(int)Math.floor(sort.length/2)])/2;
            }
            for(int j = 0; j < d[i].getValues().length; j++){
                if(d[i].getValues()[j]==-99999999) {
                    d[i].getValues()[j] = (s.mvMedian) ? median : (sort.length>0) ? (sum/sort.length) : 0;
                    // if(d[j].getValues()[i]<-100) System.out.println(median+" "+sum+" "+sort.length);
                    d[i].setOriginalValue(j,d[i].getValues()[j]);
                }
            }
           
        }
         }

      
         return d;
     }
     
     //write clusterRuns to disk to save memory
     private void writeTemp(clusterRun cr, int start, Settings s){
         try{
             FileOutputStream fos = new FileOutputStream(s.outputDirectory+s.getFolderDivider()+s.getName()+"_temp"+s.getFolderDivider()+start);
             ObjectOutputStream oos = new ObjectOutputStream(fos);
             oos.writeObject(cr);
             oos.close();
         }catch(IOException err){};
     }

     public void kill(){
         if(threads==null) return;
         for(int i = 0; i < threads.length; i++){
             threads[i].stop();
         }
     }

     public clusterRun runNewPValue(ArrayList clusterRuns, Settings s, viewer2D v2d){
         ensemble = (new Ensemble(clusterRuns, true, 0, (!s.doKmeans && !s.doHierarchical) ? true : false, s, monitor)).run();
         return doOutput(s, clusterRuns);
     }

     public ArrayList getMappingArrayList() {return storeMapping;}


     public class launchSplash implements Runnable{
        public void run(){
            splash = new javax.swing.JFrame();
            splash.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
            splash.setUndecorated(true);
            JLabel jLabel1 = new javax.swing.JLabel();
            javax.swing.ImageIcon icon = new javax.swing.ImageIcon(getClass().getResource("/imgs/autosomebanner.jpg"));
            jLabel1.setIcon(icon);
            splash.add(jLabel1);
            splash.pack();
            splash.setLocationRelativeTo(null);
            splash.setAlwaysOnTop(true);
            splash.setVisible(true);
            appInit();
        }
        public void appInit(){
              try
            {
                Thread.sleep(2200);
            }
            catch (InterruptedException ex)
            {
                // ignore it
         //   }
        }
        splash.setVisible(false);
        }
    }

     
    private int runCount = 0; //how many clusters runs so far?
    private int progressCount = 0;
    private int totalProgress = 0;
    private ArrayList clusterRuns = new ArrayList(); //store all cluster runs
    private clusterRun ensemble; //combined clustering from ensemble
    private Thread[] threads; //all AutoSOME threads
    private ArrayList storeMapping = new ArrayList(); //store mapping results before MST clustering
    private static boolean runGUI = false;
}

