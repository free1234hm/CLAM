package bsm;

import java.util.List;
import bsm.core.StatUtil;

public class EnrichmentAnalysis {

	public String[][] enrichment(List<Integer> genelist, String[] regNames,
			int[][] bindingdata, int numgenes){
		String[][] table = new String[regNames.length+1][4];
		table[0][0] = "Name";
		table[0][1] = "Gene Number";
		table[0][2] = "Overlaps";
		table[0][3] = "Pvalue";
		for(int i=0;i<regNames.length;i++){
			int[] tfgene = bindingdata[i];
			if(tfgene != null && tfgene.length > 0){
				int all_overlap = 0;
				for(int tg:tfgene){
					if(genelist.contains(tg)){
						all_overlap++;
					}
				}
				double pvalue = StatUtil.hypergeometrictail(all_overlap-1, tfgene.length, 
						numgenes-tfgene.length, genelist.size());
				table[i+1][0] = regNames[i];
				table[i+1][1] = tfgene.length+"";
				table[i+1][2] = all_overlap+"";
				table[i+1][3] = pvalue+"";
			}else{
				table[i+1][0] = regNames[i];
				table[i+1][1] = 0+"";
				table[i+1][2] = 0+"";
				table[i+1][3] = 1+"";
			}
		}
		return table;
	}
	
}
