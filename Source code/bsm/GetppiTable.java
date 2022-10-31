package bsm;

import java.util.HashMap;
import java.util.List;

public class GetppiTable extends Thread {
	private int startNum;
	PPIInteractionData ppiData;
	int numrows;
	int segment;
	public static double[][] coregulator;
	public GetppiTable(PPIInteractionData ppiData, int numrows, int segment, int startNum) {
		this.startNum = startNum;
		this.ppiData = ppiData;
		this.segment = segment;
		this.numrows = numrows;
	}
	public GetppiTable() {
		// TODO Auto-generated constructor stub
	}
	public synchronized void add(double[][] sub, int startnum) {//������Ϊ�̰߳�ȫ�ģ��÷������ܱ�����߳�ͬʱִ�У�
		for(int i=startnum;i<Math.min(startnum+segment, numrows);i++){
			coregulator[i] = sub[i-startnum];
		}
	}
	public void run() {
		double[][] sub = new double[Math.min(segment, numrows-startNum)][numrows];
		makeCoregulatorIndex(ppiData.ppiindex, ppiData.ppivalue, 
				sub, startNum, Math.min(startNum+segment, numrows));
		
		for(int i=0;i<sub.length;i++){
			double sum = 0;
			for(int j=0;j<numrows;j++){
				sum += sub[i][j];
			}
			if(sum > 0){
				double sum2 = 0;
				for(int j=0;j<numrows;j++){
					double v = Math.exp(sub[i][j]/sum);
					sub[i][j] = v;
					sum2 += v;
				}
				sum2--;
				for(int j=0;j<numrows;j++){
					sub[i][j] /= sum2;
				}
			}else{
				for(int j=0;j<numrows;j++){
					sub[i][j] = (double) 1/(numrows-1);
				}
			}
		}
		add(sub, startNum);
	}
	
	public void multithread(PPIInteractionData ppiData, int numrows, int segment) {
		int divide = numrows/segment;
		coregulator = new double[numrows][numrows];
		
		Thread[] threadList = new Thread[divide+1];
		for (int i = 0; i < divide+1; i++) {
			threadList[i] = new GetppiTable(ppiData, numrows, segment, segment * i);
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
	private void makeCoregulatorIndex(HashMap<Integer, List<Integer>> valIndex, 
			HashMap<Integer, List<Double>> val, double[][] coregulation, int start, int end) {	
		for(int i=start;i<end;i++){
			List<Integer> tglist = valIndex.get(i);
			List<Double> tgvalue = val.get(i);
			if(tglist != null && tglist.size() > 0){
				for(int p=0;p<tglist.size();p++){
					coregulation[i-start][tglist.get(p)] += Math.abs(tgvalue.get(p)) / tglist.size();
				}
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


