package bsm;

import org.apache.commons.math3.distribution.NormalDistribution;

public class GetpdfMatrix extends Thread {
	private int startNum;
	double[][] corrmatrix;
	int numrows;
	int segment;
	NormalDistribution normal;
	double traincorr; 
	public static double[][] pdftable;
	public GetpdfMatrix(double[][] corrmatrix, NormalDistribution normal, 
			double traincorr, int numrows, int segment, int startNum) {
		this.startNum = startNum;
		this.corrmatrix = corrmatrix;
		this.segment = segment;
		this.numrows = numrows;
		this.normal = normal;
		this.traincorr = traincorr;
	}
	public GetpdfMatrix() {
		// TODO Auto-generated constructor stub
	}
	public synchronized void add(double[][] sub, int startnum) {//被设置为线程安全的！该方法不能被多个线程同时执行！
		//System.out.println(sub.length);
		for(int i=startnum;i<Math.min(startnum+segment, numrows);i++){
			pdftable[i] = sub[i-startnum];
		}
	}
	public void run() {
		double[][] sub = new double[Math.min(segment, numrows-startNum)][numrows];
		for (int i = startNum; i < Math.min(startNum+segment, numrows); i++){
			for(int j=i+1;j<numrows;j++){
				sub[i-startNum][j] = getProportion(normal, corrmatrix[i][j], traincorr);
			}
		}
		add(sub, startNum);
	}
	public void multithread(double[][] corrmatrix, NormalDistribution normal, 
			double traincorr, int numrows, int segment) {
		int divide = numrows / segment;
		pdftable = new double[numrows][numrows];
		
		Thread[] threadList = new Thread[divide+1];
		for (int i = 0; i < divide+1; i++) {
			threadList[i] = new GetpdfMatrix(corrmatrix, normal, traincorr, numrows, segment, segment * i);
			/* 这里将分别传入每组的初始数据 */
			threadList[i].start();
		}
		for (int i = 0; i < divide+1; i++) {
			try {
				threadList[i].join();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	private double getProportion(NormalDistribution normal, double corr, double traincorr){
		double pdf = 0;
		if(Math.abs(corr) > 1){
			throw new IllegalArgumentException("The correlation coefficient should be between -1 and 1.");
		}else if(Math.abs(corr) > traincorr){
			double pdfposi = normal.cumulativeProbability(Math.abs(corr));
			double pdfnega = normal.cumulativeProbability(-Math.abs(corr));
			double pdfdelta = normal.cumulativeProbability(traincorr)
					- normal.cumulativeProbability(-traincorr);
			double P1 = normal.cumulativeProbability(1);
			double P2 = normal.cumulativeProbability(-1);
			pdf = (pdfposi - pdfnega - pdfdelta)/(P1-P2-pdfdelta);
			//result = 1 / (1 + Math.exp(-2*Math.log(fdr)));
		}
		return pdf;
	}
}


