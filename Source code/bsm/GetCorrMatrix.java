package bsm;

import bsm.core.Similarity;

public class GetCorrMatrix extends Thread {
	private int startNum;
	double[][] data;
	int numrows;
	int similarity;
	int segment;
	public static double[][] sum;
	public GetCorrMatrix(double[][] data, int numrows, int similarity, int segment, int startNum) {
		this.startNum = startNum;
		this.data = data;
		this.segment = segment;
		this.numrows = numrows;
		this.similarity = similarity;
	}
	public GetCorrMatrix() {
		// TODO Auto-generated constructor stub
	}
	public synchronized void add(double[][] sub, int startnum) {//������Ϊ�̰߳�ȫ�ģ��÷������ܱ�����߳�ͬʱִ�У�
		//System.out.println(sub.length);
		for(int i=startnum;i<Math.min(startnum+segment, numrows);i++){
			sum[i] = sub[i-startnum];
		}
	}
	public void run() {
		Similarity pr = new Similarity();
		double[][] sub = new double[Math.min(segment, numrows-startNum)][numrows];
		double[] exp1, exp2;
		for (int i = startNum; i < Math.min(startNum+segment, numrows); i++){
			exp1 = data[i];
			for(int j=i+1;j<numrows;j++){
				exp2 = data[j];
				sub[i-startNum][j] = pr.measure(exp1, exp2, similarity);
			}
		}
		add(sub, startNum);
	}
	public void multithread(double[][] data, int numrows, int similarity, int segment) {
		int divide = numrows/segment;
		sum = new double[numrows][numrows];
		
		Thread[] threadList = new Thread[divide+1];
		for (int i = 0; i < divide+1; i++) {
			threadList[i] = new GetCorrMatrix(data, numrows, similarity, segment, segment * i);
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


