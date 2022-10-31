package bsm;


import javax.swing.*;

import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.*;
import java.text.*;
import java.io.*;

/**
 * This class implements the core methods for learning the BSM models
 */

public class LearningCase {

	JTextArea casetext;

	static final String SZDELIM = "|;,";
	int knn;
	double similarity;
	int numrows;
	int nmingenes;
	
	/**
	 * Object containing all of the Regulator-Gene binding data.
	 */
	List<String> genelist;
	String[] reglist;
	List<RegulatorBindingData> pathwayData;
	PPIInteractionData ppiData;
	String[][] mintftable;
	String[][] minpathwaytable;
	String[][] minmirnatable;
	HashMap<Integer, String[][]> tfenrichment;
	HashMap<Integer, String[][]> pathwayenrichment;
	HashMap<Integer, String[][]> mirnaenrichment;
	
	Double[][] knnweight;
	Integer[][] knnindex;
	double[][] Pmatrix;
	List<Integer> CSOs;
	List<Integer> Outliers;
	List<Integer> Others;
	HashMap<Integer, List<Integer>> resultlist;
	int ModuleNum;
	double NDmean;
	double NDsigma;
	NormalDistribution Normalcorr;
	double Sumcorr;   
	double Sumcorrsq;
	NumberFormat nf4;
	double dlasterror;
	double currenerror;
	List<BSM_DataSet> theDataSet; 

	public LearningCase(List<String> genelist, int[][] measuredtable, List<BSM_DataSet> theDataSet, 
			List<RegulatorBindingData> pathwayData, PPIInteractionData ppiData, String[] totalregnames,
			int[][] totalreg2geneindex, int[][] totalgene2regindex,
			boolean busego, String szknn, int similarity, String mingenes,
			JTextArea casetext, JButton endSearchButton) throws Exception {
	
		nf4 = NumberFormat.getInstance(Locale.ENGLISH); //nf3 transform double to 4 decimal places
		nf4.setMinimumFractionDigits(6);
		nf4.setMaximumFractionDigits(6);
		
		this.genelist = genelist;
		this.reglist = totalregnames;
		this.theDataSet = theDataSet;
		this.knn = Integer.parseInt(szknn);
		this.nmingenes = Integer.parseInt(mingenes);
		this.pathwayData = pathwayData;
		this.ppiData = ppiData;
		this.casetext = casetext;
		this.numrows = genelist.size();
		this.similarity = similarity;
		
	    casetext.append("Conbined object number: "+numrows+"\n");
	    casetext.paintImmediately(casetext.getBounds());
	    
		int[][] state = new int[numrows][theDataSet.size()];
		Integer[][] totalknnindex = new Integer[numrows][knn*theDataSet.size()];
		double[][] totalknnweight = new double[numrows][knn*theDataSet.size()];
	    
	    if(theDataSet != null && theDataSet.size()>0) {
	    	for(int index=0;index<theDataSet.size();index++) {
	    		BSM_DataSet data = theDataSet.get(index);
	    		 GetCorrMatrix gcm = new GetCorrMatrix();
	    		    gcm.multithread(data.controlnorm, data.numrows, similarity, 200);
	    			double[][] corrmatrix = gcm.sum;
	    			for(int i=0;i<data.numrows;i++){
	    				for(int j=0;j<i;j++){
	    					corrmatrix[i][j] = corrmatrix[j][i];
	    				}
	    			}
	    			/*
	    			GetpdfMatrix gpf = new GetpdfMatrix();
	    			gpf.multithread(corrmatrix, Normalcorr, trainPearson, data.numrows, 100);
	    			double[][] pdfmatrix = gpf.pdftable;
	    			for(int i=0;i<data.numrows;i++){
	    				for(int j=0;j<i;j++){
	    					pdfmatrix[i][j] = pdfmatrix[j][i];
	    				}
	    			}
	    			*/
	    			GetKNN2 getknn = new GetKNN2();
	    			getknn.multithread(corrmatrix, data.numrows, 200, knn);
	    			double[][] knnvalue = getknn.knnvalue;
	    			int[][] knnindex = getknn.knnindex;
	    			double[][] knnweight = getknn.knnweight;
	    			double[] density = getknn.density;
	    			
	    			//成为CSO有两个条件：1）在其他数据中既不是空也不是outlier
	    			//成为outlier条件：1）在鉴定到的所有数据中是outlier

	    			for(int i=0;i<data.numrows;i++){
	    				for(int j=0;j<knn;j++){
	    					totalknnindex[data.this2all[i]][index*knn+j] = knnindex[i][j];
	    					totalknnweight[data.this2all[i]][index*knn+j] = knnweight[i][j];
	    				}
	    				if(density[i] > 0){
	    					boolean  iscso = true;
	    					//boolean isoutlier = true;
	    					for(int j=0;j<knn;j++){
	    						if(density[i] < density[knnindex[i][j]]) {
	    							iscso = false;
	    							break;
	    						}
	    					}
	    					if(iscso) state[data.this2all[i]][index] = 1;
	    					//if(isoutlier) state[data.this2all[i]][index] = -1;
	    				} else {
	    					state[data.this2all[i]][index] = -1;
	    				}
	    			}
	    	}
	    } else {
	    	throw new IllegalArgumentException("No expression data file found.");
	    }
	    
	    knnweight = new Double[numrows][];
		knnindex = new Integer[numrows][];
		for(int i=0;i<numrows;i++) {
			List<Integer> r1 = new ArrayList<Integer>();
			List<Double> r2 = new ArrayList<Double>();
			for(int j=0;j<totalknnindex[i].length;j++) {
				if(totalknnindex[i][j] != null) {
					r1.add(totalknnindex[i][j]);
					r2.add(totalknnweight[i][j]);
				}
			}
			knnindex[i] = r1.toArray(new Integer[r1.size()]);
			knnweight[i] = r2.toArray(new Double[r2.size()]);
		}
		
	    int[] allstate = new int[numrows];
	    for(int i=0;i<state.length;i++) {
	    	boolean isout1 = false;
	    	boolean isout2 = true;
	    	for(int j=0;j<state[i].length;j++) {
	    		BSM_DataSet data = theDataSet.get(j);
	    		if(state[i][j]==-1) {
	    			isout1 = true;
	    		} else if(data.all2this[i] != null){
	    			isout2 = false;
	    		}
	    	}
	    	if(isout1 == true && isout2 == true) {
	    		allstate[i] = -1;
	    	}
	    }
	    for(int i=0;i<state.length;i++) {
	    	boolean iscso1 = false;
	    	boolean iscso2 = true;
	    	for(int j=0;j<state[i].length;j++) {
	    		BSM_DataSet data = theDataSet.get(j);
	    		if(state[i][j]==1) {
	    			iscso1 = true;
	    		} else if(data.all2this[i] == null || state[i][j]==-1){
	    			iscso2 = false;
	    		}
	    	}
	    	if(iscso1 == true && iscso2 == true) {
	    		allstate[i] = 1;
	    	}
	    }
	    CSOs = new ArrayList<Integer>();
		Outliers = new ArrayList<Integer>();
		Others = new ArrayList<Integer>();
	    for(int i=0;i<numrows;i++) {
	    	if(allstate[i]==1) {
	    		CSOs.add(i);
	    	} else if(allstate[i]==-1) {
	    		Outliers.add(i);
	    	} else {
	    		Others.add(i);
	    	}
	    }
	    casetext.append("Center:"+CSOs.size()+" Outlier:"+Outliers.size()+ " Others:"+Others.size()+"\n");
	    casetext.paintImmediately(casetext.getBounds());
	    
	    for(int i=0;i<numrows;i++) {
	    	if(busego && totalreg2geneindex != null && totalgene2regindex != null) {
	    		boolean constanta = false;
	    		double[] dist = GettfTable(i, knnindex[i], totalreg2geneindex, totalgene2regindex, constanta);
	    		double[] prior;
	    		if(constanta) {
	    			prior = dist;
	    		} else {
	    			prior = new double[dist.length];
	    			double sum = 0;
		    		for(int j=0;j<dist.length;j++) {
		    			prior[j] = Math.exp(dist[j]);
		    			sum += prior[j];
		    		}
		    		for(int j=0;j<prior.length;j++) {
		    			prior[j] /= sum;
		    		}
	    		}
	    		for(int j=0;j<knnweight[i].length;j++) {
	    			knnweight[i][j] = prior[j]*knnweight[i][j];
	    		}
	    	}
	    	if(busego && ppiData != null) {
	    		boolean constanta = false;
	    		double[] dist = GetppiTable(i, knnindex[i], ppiData.ppiindex, ppiData.ppivalue, constanta);
	    		double[] prior;
	    		if(constanta) {
	    			prior = dist;
	    		} else {
	    			prior = new double[dist.length];
	    			double sum = 0;
		    		for(int j=0;j<dist.length;j++) {
		    			prior[j] = Math.exp(dist[j]);
		    			sum += prior[j];
		    		}
		    		for(int j=0;j<prior.length;j++) {
		    			prior[j] /= sum;
		    		}
	    		}
	    		for(int j=0;j<knnweight[i].length;j++) {
	    			knnweight[i][j] = prior[j]*knnweight[i][j];
	    		}
	    	}
	    	
	    	 double sum = 0;
	    	 for(int j=0;j<knnweight[i].length;j++) sum += knnweight[i][j];
	    	 for(int j=0;j<knnweight[i].length;j++) knnweight[i][j] /= sum;
	    }
	    
	    
		ModuleNum = CSOs.size()+1;
		InitePmatrix();
		searchstage1();
		
		/*
	    try{
			BufferedWriter outXml2 = new BufferedWriter(new FileWriter("D:/Pmatrix.txt"));
	    	for(int i=0;i<Pmatrix.length;i++){
	    		for(int j=0;j<Pmatrix[i].length-1;j++){
	    			outXml2.write(Pmatrix[i][j] + "\t");
	    		}
	    		outXml2.write(Pmatrix[i][Pmatrix[i].length-1] + "\n");
	    	}
			outXml2.flush(); 
			outXml2.close();
			System.out.println("DONE");
			
	    }catch (Exception e){
			System.out.println("FALSE"); 
		e.printStackTrace(); 
		}
	    */
		
		assigngene(Pmatrix);
		String[][] resulttable = new String[numrows][2];
		int count = 0;
		Iterator<Map.Entry<Integer, List<Integer>>> entries = resultlist.entrySet().iterator();
		while(entries.hasNext()){
			Map.Entry<Integer, List<Integer>> entry = entries.next();
			int key = entry.getKey();
			List<Integer> set = entry.getValue();
			if(key == ModuleNum-1){
				for(int i=0;i<set.size();i++){
					resulttable[count][0] = genelist.get(set.get(i));
					resulttable[count][1] = "Outlier";
					count++;
				}
			} else {
				for(int i=0;i<set.size();i++){
					resulttable[count][0] = genelist.get(set.get(i));
					resulttable[count][1] = key+"";
					count++;
				}
			}
		}
		
		try{
			File filePath = new File("Results");
			if (!filePath.exists()) {
				filePath.mkdir();
			}
	    	BufferedWriter outXml1 = new BufferedWriter(new FileWriter("Results/Final modules.txt"));
	    	outXml1.write("Name" + "\t" + "Module" + "\n");
	    	for(int i=0;i<resulttable.length;i++){
	    		for(int j=0;j<resulttable[i].length-1;j++){
	    			outXml1.write(resulttable[i][j] + "\t");
	    		}
	    		outXml1.write(resulttable[i][resulttable[i].length-1] + "\n");
	    	}
			outXml1.flush(); 
			outXml1.close();
			System.out.println("DONE");
			
	    } catch (Exception e) {
			System.out.println("FALSE"); 
		e.printStackTrace(); 
		}
		
		if(pathwayData != null && pathwayData.size() > 0) {
			File filePath = new File("Results/Enrichment results");
			if (!filePath.exists()) {
				filePath.mkdir();
			}
			for(int reg=0;reg<pathwayData.size();reg++) {
				RegulatorBindingData tfData = pathwayData.get(reg);
				HashMap<Integer, String[][]> tfenrichment = new HashMap<Integer, String[][]>();
				String[][] mintftable = new String[tfData.regNames.length+1][2];
				enrichmentanalysis(tfData, tfenrichment, mintftable);
				try{
					File filePath2 = new File(filePath+"/Enrichment "+reg);
					if (!filePath2.exists()) {
						filePath2.mkdir();
					}
					Iterator<Map.Entry<Integer, String[][]>> entries2 = tfenrichment.entrySet().iterator();
					while(entries2.hasNext()){
						Map.Entry<Integer, String[][]> entry = entries2.next();
						int key = entry.getKey();
						String[][] set = entry.getValue();
						BufferedWriter outXml1 = new BufferedWriter(new FileWriter(filePath2+"/Cluster"+key+".txt"));
				    	for(int i=0;i<set.length;i++){
				    		for(int j=0;j<set[i].length-1;j++){
				    			outXml1.write(set[i][j] + "\t");
				    		}
				    		outXml1.write(set[i][set[i].length-1] + "\n");
				    	}
						outXml1.flush(); 
						outXml1.close();
						System.out.println("DONE");
					}
					
			    	BufferedWriter outXml1 = new BufferedWriter(new FileWriter(filePath2+"/Min pvalues.txt"));
			    	for(int i=0;i<mintftable.length;i++){
			    		for(int j=0;j<mintftable[i].length-1;j++){
			    			outXml1.write(mintftable[i][j] + "\t");
			    		}
			    		outXml1.write(mintftable[i][mintftable[i].length-1] + "\n");
			    	}
					outXml1.flush(); 
					outXml1.close();
					System.out.println("DONE");
					
			    }catch (Exception e) {
					System.out.println("FALSE"); 
				e.printStackTrace(); 
				}
			}
			
		}
		
	}
	
	private double[] GettfTable(int i, Integer[] knnindex, int[][] totalreg2geneindex, 
			int[][] totalgene2regindex, boolean constanta) {
		double[] knnvalue = new double[knnindex.length];
			if(totalgene2regindex[i] != null && totalgene2regindex[i].length>0) {
				for(int j=0;j<totalgene2regindex[i].length;j++){
					int nReg = totalgene2regindex[i][j];
					int[] list1 = totalreg2geneindex[nReg];
					if(list1.length > 1) {
						for(int m=0;m<list1.length;m++) {
							if(list1[m] != i) {
								for(int n=0;n<knnindex.length;n++) {
									if(list1[m] == knnindex[n]) {
										knnvalue[n] += 1/(list1.length-1);
										break;
									}
								}
							}
						}
					}
				}
				double sum=0;
				for(int j=0;j<knnvalue.length;j++) sum += knnvalue[j];
				if(sum > 0) {
					for(int j=0;j<knnvalue.length;j++) knnvalue[j] /= sum;
				} else {
					constanta = true;
					knnvalue = CONSTANTA(knnindex.length);
				}
			} else {
				constanta = true;
				knnvalue = CONSTANTA(knnindex.length);
			}
			return knnvalue;
	}
	
	private double[] GetppiTable(int i, Integer[] knnindex, HashMap<Integer, List<Integer>> ppiindex,
			HashMap<Integer, List<Double>> ppivalue, boolean constanta) {
		double[] knnvalue = new double[knnindex.length];
		if(ppiindex.get(i) != null && ppivalue.get(i) != null) {
			List<Integer> protein = ppiindex.get(i);
			List<Double> value = ppivalue.get(i);
				for(int j=0;j<protein.size();j++){
					for(int n=0;n<knnindex.length;n++) {
						if(protein.get(j) == knnindex[n]) {
							knnvalue[n] += Math.abs(value.get(j));
							break;
						}
					}
				}
				double sum=0;
				for(int j=0;j<knnvalue.length;j++) sum += knnvalue[j];
				if(sum > 0) {
					for(int j=0;j<knnvalue.length;j++) knnvalue[j] /= sum;
				} else {
					constanta = true;
					knnvalue = CONSTANTA(knnindex.length);
				}
		} else {
			constanta = true;
			knnvalue = CONSTANTA(knnindex.length);
		}
		return knnvalue;
	}
	
	
	public double[] CONSTANTA(int nclass) {
		double[] ptran = new double[nclass];
		for(int i=0;i<ptran.length;i++){
			ptran[i] = (double)1/nclass;
		}
		return ptran;
	}
	
	private void InitePmatrix(){
		Pmatrix = new double[numrows][ModuleNum];
		for(int i=0;i<CSOs.size();i++){
			Pmatrix[CSOs.get(i)][i] = 1;
		}
		for(int i=0;i<Outliers.size();i++){
			Pmatrix[Outliers.get(i)][ModuleNum-1] = 1;
		}
		for(int i=0;i<Others.size();i++){
			for(int j=0;j<ModuleNum;j++){
				Pmatrix[Others.get(i)][j] = (double) 1/ModuleNum;
			}
		}
	}
	private double GetError(double[][] Pmatrix, double[][] updatedPmatrix){
		double error = 0;
		for(int gene:Others){
			for(int i=0;i<ModuleNum;i++){
				error += Math.pow((Pmatrix[gene][i] - updatedPmatrix[gene][i]), 2);
			}
		}
		return error;
	}
	private void assigngene(double[][] Pmatrix){
		HashMap<Integer, List<Integer>> genelist = new HashMap<Integer, List<Integer>>();
		for(int i=0;i<numrows;i++){
			double maxp = 0;
			int maxcso = -1;
			for(int j=0;j<ModuleNum;j++){
				if(maxp < Pmatrix[i][j]){
					maxp = Pmatrix[i][j];
					maxcso = j;
				}
			}
			if(genelist.get(maxcso) != null){
				List<Integer> set = genelist.get(maxcso);
				set.add(i);
				genelist.put(maxcso, set);
			}else{
				List<Integer> set = new ArrayList<Integer>();
				set.add(i);
				genelist.put(maxcso, set);
			}
		}
		
		resultlist = new HashMap<Integer, List<Integer>>();
		Iterator<Map.Entry<Integer, List<Integer>>> entries = genelist.entrySet().iterator();
		while(entries.hasNext()){
			Map.Entry<Integer, List<Integer>> entry = entries.next();
			int key = entry.getKey();
			List<Integer> set = entry.getValue();
			if(set.size() < nmingenes){
				if(resultlist.get(ModuleNum-1) != null){
					List<Integer> outlier = resultlist.get(ModuleNum-1);
					outlier.addAll(set);
					resultlist.put(ModuleNum-1, outlier);
				} else {
					resultlist.put(ModuleNum-1, set);
				}
			} else {
				resultlist.put(key, set);
			}
		}
	}
	
	private void assigngene(double[][] Pmatrix, double threshold){
		HashMap<Integer, List<Integer>> genelist = new HashMap<Integer, List<Integer>>();
		for(int i=0;i<numrows;i++){
			for(int j=0;j<ModuleNum;j++){
				if(Pmatrix[i][j] > threshold){
					if(genelist.get(j) != null){
						List<Integer> set = genelist.get(j);
						set.add(i);
						genelist.put(j, set);
					}else{
						List<Integer> set = new ArrayList<Integer>();
						set.add(i);
						genelist.put(j, set);
					}
				}
			}
		}
		resultlist = new HashMap<Integer, List<Integer>>();
		if(genelist.get(ModuleNum-1) != null){
			resultlist.put(ModuleNum-1, genelist.get(ModuleNum-1));		
		}
		Iterator<Map.Entry<Integer, List<Integer>>> entries = genelist.entrySet().iterator();
		while(entries.hasNext()){
			Map.Entry<Integer, List<Integer>> entry = entries.next();
			int key = entry.getKey();
			List<Integer> set = entry.getValue();
			if(key != ModuleNum-1){
				if(set.size() < nmingenes){
					if(resultlist.get(ModuleNum-1) != null){
						List<Integer> outlier = resultlist.get(ModuleNum-1);
						outlier.addAll(set);
						resultlist.put(ModuleNum-1, outlier);		
					}else{
						resultlist.put(ModuleNum-1, set);
					}
				}else{
					resultlist.put(key, set);
				}
			}
		}
	}
	/**
	 * First step of tree searching��add and delete path
	 */
	public void searchstage1() throws Exception {
		UpdatePmatrix up = new UpdatePmatrix();
		up.multithread(Pmatrix, Others, knnweight, knnindex, numrows, ModuleNum, 100);
		double[][] updatedPmatrix = up.newpmatrix;
		for(int cso:CSOs) updatedPmatrix[cso] = Pmatrix[cso];
		for(int out:Outliers) updatedPmatrix[out] = Pmatrix[out];
		dlasterror = GetError(Pmatrix, updatedPmatrix);
		Pmatrix = updatedPmatrix;
		
		boolean bendsearchlocal;
		boolean keepsearching;
		do {
			keepsearching = false;
			up.multithread(Pmatrix, Others, knnweight, knnindex, numrows, ModuleNum, 100);
			double[][] currentmatrix = up.newpmatrix;
			for(int cso:CSOs) currentmatrix[cso] = Pmatrix[cso];
			for(int out:Outliers) currentmatrix[out] = Pmatrix[out];
			double currenterror = GetError(Pmatrix, currentmatrix);	
			
			bendsearchlocal = Main_interface2.bendsearch; //bend searching
			casetext.append(" Previous error: "+dlasterror+"; Current error: "+currenterror+"\n");
			casetext.paintImmediately(casetext.getBounds());
			
			if(currenterror < dlasterror){
				keepsearching = true;
				Pmatrix = currentmatrix;
				dlasterror = currenterror;
			}
		}while(!bendsearchlocal && keepsearching);
	}
	
	public void enrichmentanalysis(RegulatorBindingData bindingdata,
			HashMap<Integer, String[][]> enrichment, String[][] mintable){
		EnrichmentAnalysis EA = new EnrichmentAnalysis();
	    Iterator<Map.Entry<Integer, List<Integer>>> entries = resultlist.entrySet().iterator();
		while(entries.hasNext()){
			Map.Entry<Integer, List<Integer>> entry = entries.next();
			int key = entry.getKey();
			if(key != ModuleNum-1){
				List<Integer> set = entry.getValue();
				String[][] table = EA.enrichment(set, bindingdata.regNames,
						bindingdata.reg2GeneDedupliIndex, numrows);
				enrichment.put(key, table);
			}	
		}
		mintable[0][0] = "Name";
		mintable[0][1] = "Pvalue";
	    for(int i=0;i<bindingdata.regNames.length;i++){
	    	double minpvalue = 1;
	    	for(String[][] table:enrichment.values()){
		    	double pval = Double.parseDouble(table[i+1][3]);
		    	if(minpvalue > pval){
		    		minpvalue = pval;
		    	}
		    }
	    	mintable[i+1][0] = bindingdata.regNames[i];
	    	mintable[i+1][1] = minpvalue+"";
	    }
	}
}
