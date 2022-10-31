package bsm;

import java.util.List;

public class UpdatePmatrix extends Thread {
	private int startNum;
	private double[][] pmatrix;
	private Double[][] knnweight;
	private Integer[][] knnindex;
	List<Integer> genes;
	private int nummodule;
	private int segment;
	public static double[][] newpmatrix;
	public UpdatePmatrix(double[][] pmatrix, List<Integer> genes, Double[][] knnweight, 
			Integer[][] knnindex, int nummodule, int segment, int startNum) {
		this.startNum = startNum;
		this.pmatrix = pmatrix;
		this.knnweight = knnweight;
		this.knnindex = knnindex;
		this.genes = genes;
		this.segment = segment;
		this.nummodule = nummodule;
	}

	public UpdatePmatrix() {
		// TODO Auto-generated constructor stub
	}

	public void run() {
		for(int i=startNum;i<Math.min(startNum+segment, genes.size());i++){
			int totalindex = genes.get(i);
			Integer[] neibours = knnindex[totalindex];
			Double[] weights = knnweight[totalindex];
			for(int j=0;j<neibours.length;j++){
				if(weights[j] > 0){
					for(int m=0;m<nummodule;m++) newpmatrix[totalindex][m] += weights[j] * pmatrix[neibours[j]][m];
				}
			}
		}
	}
	
	public void multithread(double[][] pmatrix,  List<Integer> genes, 
			Double[][] knnweight, Integer[][] knnindex, int numrows, 
			int nummodule, int segment) {
		
		int divide = genes.size()/segment;
		newpmatrix = new double[numrows][nummodule];
		
		Thread[] threadList = new Thread[divide+1];
		for (int i = 0; i < divide+1; i++) {
			threadList[i] = new UpdatePmatrix(pmatrix, genes, knnweight, 
					knnindex, nummodule, segment, segment * i);
			/* ���ｫ�ֱ���ÿ��ĳ�ʼ���� */
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
	/*
	public static void main(String[] args) {
		double[][] a = null;
		multithread(a);
	}
	*/
}


