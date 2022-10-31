package bsm;


public class GetKNN2 extends Thread {
	public static double[][] knnvalue;
	public static double[][] knnweight;
	public static int[][] knnindex;
	public static double[] density;
	private int startNum;
	double[][] pdftable;
	int numrows;
	int segment;
	int k;
	public GetKNN2(double[][] pdftable, int numrows, int segment, int k, int startNum) {
		this.startNum = startNum;
		this.pdftable = pdftable;
		this.segment = segment;
		this.numrows = numrows;
		this.k = k;
	}

	public GetKNN2() {
		// TODO Auto-generated constructor stub
	}

	public void run() {	
		for (int index = startNum; index < Math.min(startNum+segment, numrows); index++){
			double[] corrvalue  = new double[k];
			int[] corrindex  = new int[k];
			if(k > index){
				for(int i=0;i<index;i++){
					corrindex[i] = i;
					corrvalue[i] = pdftable[index][i];
				}
				for(int i=index;i<k;i++){
					corrindex[i] = i+1;
					corrvalue[i] = pdftable[index][i+1];
				}
				BubbleSort(corrvalue, corrindex);
				/******************************************/
				for (int i = k+1; i < numrows; i++){
					double pdf = pdftable[index][i];
					if(pdf > corrvalue[k-1]){
						for(int j=0;j<k;j++){
							if(pdf > corrvalue[j]){
								for(int p=k-1;p>j;p--){
									corrvalue[p] = corrvalue[p-1];
									corrindex[p] = corrindex[p-1];
								}
								corrvalue[j] = pdf;
								corrindex[j] = i;
								break;
							}
						}
					}
				}
			}else{
				for(int i=0;i<k;i++){
					corrindex[i] = i;
					corrvalue[i] = pdftable[index][i];
				}
				BubbleSort(corrvalue, corrindex);
				/******************************************/
				for (int i = k; i < index; i++){
					double pdf = pdftable[index][i];
					if(pdf > corrvalue[k-1]){
						for(int j=0;j<k;j++){
							if(pdf > corrvalue[j]){
								for(int p=k-1;p>j;p--){
									corrvalue[p] = corrvalue[p-1];
									corrindex[p] = corrindex[p-1];
								}
								corrvalue[j] = pdf;
								corrindex[j] = i;
								break;
							}
						}
					}
				}
				for (int i = index+1; i < numrows; i++){
					double pdf = pdftable[index][i];
					if(pdf > corrvalue[k-1]){
						for(int j=0;j<k;j++){
							if(pdf > corrvalue[j]){
								for(int p=k-1;p>j;p--){
									corrvalue[p] = corrvalue[p-1];
									corrindex[p] = corrindex[p-1];
								}
								corrvalue[j] = pdf;
								corrindex[j] = i;
								break;
							}
						}
					}
				}
			}
			double sum = 0;
			for(int i=0;i<k;i++) sum += corrvalue[i];
			if(sum > 0){
				density[index] = sum;
				for(int i=0;i<k;i++) knnweight[index][i] = corrvalue[i]/sum;
			}
			knnvalue[index] = corrvalue;
			knnindex[index] = corrindex;
		}
	}
	/*
	public void run() {	
		for (int index = startNum; index < Math.min(startNum+segment, numrows); index++){
			double[] corrvalue  = new double[k];
			int[] corrindex  = new int[k];
			
			if(k > index){
				for(int i=0;i<index;i++){
					corrindex[i] = i;
					corrvalue[i] = pdftable[i][index];
				}
				for(int i=index;i<k;i++){
					corrindex[i] = i+1;
					corrvalue[i] = pdftable[index][i+1];
				}
				BubbleSort(corrvalue, corrindex);
				for (int i = k+1; i < numrows; i++){
					double corr = pdftable[index][i];
					if(Math.abs(corr) > Math.abs(corrvalue[k-1])){
						for(int j=0;j<k;j++){
							if(Math.abs(corr) > Math.abs(corrvalue[j])){
								for(int p=k-1;p>j;p--){
									corrvalue[p] = corrvalue[p-1];
									corrindex[p] = corrindex[p-1];
								}
								corrvalue[j] = corr;
								corrindex[j] = i;
								break;
							}
						}
					}
				}
			}else{
				for(int i=0;i<k;i++){
					corrindex[i] = i;
					corrvalue[i] = pdftable[i][index];
				}
				BubbleSort(corrvalue, corrindex);
				for (int i = k; i < index; i++){
					double corr = pdftable[i][index];
					if(Math.abs(corr) > Math.abs(corrvalue[k-1])){
						for(int j=0;j<k;j++){
							if(Math.abs(corr) > Math.abs(corrvalue[j])){
								for(int p=k-1;p>j;p--){
									corrvalue[p] = corrvalue[p-1];
									corrindex[p] = corrindex[p-1];
								}
								corrvalue[j] = corr;
								corrindex[j] = i;
								break;
							}
						}
					}
				}
				for (int i = index+1; i < numrows; i++){
					double corr = pdftable[index][i];
					if(Math.abs(corr) > Math.abs(corrvalue[k-1])){
						for(int j=0;j<k;j++){
							if(Math.abs(corr) > Math.abs(corrvalue[j])){
								for(int p=k-1;p>j;p--){
									corrvalue[p] = corrvalue[p-1];
									corrindex[p] = corrindex[p-1];
								}
								corrvalue[j] = corr;
								corrindex[j] = i;
								break;
							}
						}
					}
				}
			}
			knnvalue[index] = corrvalue;
			knnindex[index] = corrindex;
		}
	}
	*/
	public void multithread(double[][] pdftable, int numrows, int segment, int k) {
		int divide = numrows/segment;
		knnvalue = new double[numrows][k];
		knnindex = new int[numrows][k];
		knnweight = new double[numrows][k];
		density = new double[numrows];
		
		Thread[] threadList = new Thread[divide+1];
		for (int i = 0; i < divide+1; i++) {
			threadList[i] = new GetKNN2(pdftable, numrows, segment, k, segment * i);
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
	
	static void BubbleSort(double[] value, int[] index) //反向冒泡
 	{
 		 int low = 0;   
 		 int high= value.length - 1; //设置变量的初始值  
 		 double tmp;
 		 int temp2;
 		 int j;  
 		 while (low < high) {
 		     for (j= low; j< high; ++j) //正向冒泡,找到最大者  
 		         if (Math.abs(value[j]) < Math.abs(value[j+1])) {  
 		             tmp = value[j]; value[j]=value[j+1];value[j+1]=tmp;
 		            temp2 = index[j]; index[j]=index[j+1];index[j+1]=temp2;
 		         }  
 		     --high;                 //修改high值, 前移一位  
 		     for ( j=high; j>low; --j) //反向冒泡,找到最小者  
 		         if (Math.abs(value[j]) > Math.abs(value[j-1])) {  
 		             tmp = value[j]; value[j]=value[j-1];value[j-1]=tmp;
 		            temp2 = index[j]; index[j]=index[j-1];index[j-1]=temp2;
 		         }
 		     ++low; //修改low值,后移一位  
 		 }
 	}
	
	/*
	public static void main(String[] args) {
		double[][] a = null;
		multithread(a);
	}
	*/
}


