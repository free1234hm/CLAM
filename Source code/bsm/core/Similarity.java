package bsm.core;

import java.util.Map.Entry;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

public class Similarity {
	
	public double measure(double[] vectorA, double[] vectorB, int sm) {
		double score = 0;
		if(sm == 0) {
			//System.out.println("Euclidean distance");
			double dis = 0;
			for (int i=0;i<vectorA.length;i++) {
				dis += (vectorA[i]-vectorB[i])*(vectorA[i]-vectorB[i]);
			}
			score = Math.sqrt(dis);
		} else if(sm == 1) {
			//System.out.println("Pearson correlation");
			int m = vectorA.length;
		    double[][] aa = new double[m][2];
		    for (int i = 0; i < m; i++) {
		           aa[i][0] = vectorA[i];
		           aa[i][1] = vectorB[i];
		    }
		    RealMatrix matrix = new Array2DRowRealMatrix(aa, false);
		    PearsonsCorrelation pc = new PearsonsCorrelation(matrix);
		    score = pc.getCorrelationMatrix().getEntry(0, 1);
		} else if(sm == 2) {
			//System.out.println("|Pearson correlation|");
			int m = vectorA.length;
		    double[][] aa = new double[m][2];
		    for (int i = 0; i < m; i++) {
		           aa[i][0] = vectorA[i];
		           aa[i][1] = vectorB[i];
		    }
		    RealMatrix matrix = new Array2DRowRealMatrix(aa, false);
		    PearsonsCorrelation pc = new PearsonsCorrelation(matrix);
		    score = Math.abs(pc.getCorrelationMatrix().getEntry(0, 1));
		} else if(sm == 3) {
			//System.out.println("Cosine similarity");
			double dotProduct = 0.0;
		    double normA = 0.0;
		    double normB = 0.0;
		    for (int i = 0; i < vectorA.length; i++) {
		        dotProduct += vectorA[i] * vectorB[i];
		        normA += Math.pow(vectorA[i], 2);
		        normB += Math.pow(vectorB[i], 2);
		    }   
		    score = dotProduct / (Math.sqrt(normA) * Math.sqrt(normB));
		} else if(sm == 4) {
			//System.out.println("Mutual information");
			JointProbabilityState state = new JointProbabilityState(vectorA,vectorB);
		    double jointValue, firstValue, secondValue;
		    double mutualInformation = 0.0;
		    for (Entry<Pair<Integer,Integer>,Double> e : state.jointProbMap.entrySet())
		    {
		      jointValue = e.getValue();
		      firstValue = state.firstProbMap.get(e.getKey().a);
		      secondValue = state.secondProbMap.get(e.getKey().b);

		      if ((jointValue > 0) && (firstValue > 0) && (secondValue > 0))
		      {
		        mutualInformation += jointValue * Math.log((jointValue / firstValue) / secondValue);
		      }
		    }
		    mutualInformation /= Math.log(Entropy.LOG_BASE);
		    score = mutualInformation;
		}
		return score;
	}
	
	public double euclidean(double[] vectorA, double[] vectorB) {
		double dis = 0;
		for (int i=0;i<vectorA.length;i++) {
			dis += (vectorA[i]-vectorB[i])*(vectorA[i]-vectorB[i]);
		}
		return Math.sqrt(dis);
	}
	
	public double cosineSimilarity(double[] vectorA, double[] vectorB) {
	    double dotProduct = 0.0;
	    double normA = 0.0;
	    double normB = 0.0;
	    for (int i = 0; i < vectorA.length; i++) {
	        dotProduct += vectorA[i] * vectorB[i];
	        normA += Math.pow(vectorA[i], 2);
	        normB += Math.pow(vectorB[i], 2);
	    }   
	    return dotProduct / (Math.sqrt(normA) * Math.sqrt(normB));
	}
	
	public double pearsonr(double[] x, double[] mean) {
		int m = x.length;
	    double[][] aa = new double[m][2];
	    for (int i = 0; i < m; i++) {
	           aa[i][0] = x[i];
	           aa[i][1] = mean[i];
	    }
	    RealMatrix matrix = new Array2DRowRealMatrix(aa, false);
	    PearsonsCorrelation pc = new PearsonsCorrelation(matrix);
	    double r = pc.getCorrelationMatrix().getEntry(0, 1);
		return r;
	}
	
	public double MutualInformation(double[] firstVector, double[] secondVector){
	    JointProbabilityState state = new JointProbabilityState(firstVector,secondVector);
	    double jointValue, firstValue, secondValue;
	    double mutualInformation = 0.0;
	    for (Entry<Pair<Integer,Integer>,Double> e : state.jointProbMap.entrySet())
	    {
	      jointValue = e.getValue();
	      firstValue = state.firstProbMap.get(e.getKey().a);
	      secondValue = state.secondProbMap.get(e.getKey().b);

	      if ((jointValue > 0) && (firstValue > 0) && (secondValue > 0))
	      {
	        mutualInformation += jointValue * Math.log((jointValue / firstValue) / secondValue);
	      }
	    }
	    mutualInformation /= Math.log(Entropy.LOG_BASE);
	    return mutualInformation;
	  }
	
}
