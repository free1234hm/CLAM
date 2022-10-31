package bsm;


import bsm.core.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.util.*;
import java.util.List;

/**
 * Class implementing the main input interface
 */
public class Main_interface2 extends JFrame implements ActionListener {
	
	static final String SZDELIM = "|;,";
	static final boolean BDEBUG = false;
	
	
	List<BSM_DataSet> theDataSet;
	// Main
	static String szoptionalFile = " Optional";
	static String szDataFileDEF = "";
	static boolean bspotcheckDEF = false;
	static String szGeneAnnotationFileDEF = "";
	static String szCrossRefFileDEF = "";
	static int ndbDEF = 0;
	static int nstaticsourceDEF = 0;
	static int numchildDEF = 100;
	static int nummissing = 0;

	// Repeat
	static Vector<String> vRepeatFilesDEF0 = new Vector<String>();
	static Vector<String> vRepeatFilesDEF1 = new Vector<String>();
	static Vector<String> vRepeatFilesDEF2 = new Vector<String>();
	ListDialog theRepeatList0, theRepeatList2;
	static Color lightBlue = new Color(190, 255, 190);
	Vector<String> knownmodulefiles, expressionfiles;
	static String szinteractionval;

	// Search Options
	static double dCONVERGENCEDEF = 0.01;
	static double dMinScoreDEF = 0.0;
	static double dDELAYPATHDEF = .15;
	static double dDMERGEPATHDEF = .15;
	static double dPRUNEPATHDEF = .15;
	static double dTHRESHOLD = 0.000001;
	static double dPvalueshold = 0.01;

	// Filtering
	JRadioButton log1, log2, go1, go2;
	static int nnormalizeDEF = 1;
	static int ngoDEF = 0;
	static int ntfDEF = 1;
	static int nmirnaDEF = 1;
	ButtonGroup normGroup = new ButtonGroup();
	ButtonGroup goGroup = new ButtonGroup();
	ButtonGroup tfGroup = new ButtonGroup();
	ButtonGroup mirnaGroup = new ButtonGroup();
	
	static boolean busego;
	
	String labels1[] = {" Min value", " Mean value", " Zero"};
	JComboBox comboBox1 = new JComboBox(labels1);
	String labels3[] = {" Max value", " Mean value"};
	JComboBox comboBox3 = new JComboBox(labels3);
	String labels4[] = {" Euclidean distance"," Pearson correlation"," |Pearson correlation|", 
			" Cosine similarity", " Mutual information"};
	JComboBox comboBox4 = new JComboBox(labels4);
	JLabel filterlable;
	static int missinter = 0; // Set the missing value as
	static int filterduplicates = 0; // Set the duplicate genes as
	static int similarity = 1;
	static String szPrefilteredDEF = "";
	static int nMaxMissingDEF = 0;
	static double dMinExpressionDEF = 1;
	static double dMinCorrelationRepeatsDEF = 0;

	// GUI
	static long s1;
	static boolean bendsearch = false;
	static String szsavedval;
	static String missing;
	static boolean btakelog;
	static String szknn;
	static String mingenes;

	int npathwaysourcecb;
	int ntfsourcecb;
	int nppisourcecb;
	String szuserFileField1;
	String szuserFileField2;
	String szuserFileField3;
	
	static JFileChooser theChooser = new JFileChooser();

	// Regulator Scoring GUI
	JButton regScoreHButton = new JButton(Util.createImageIcon("Help16.gif"));
	JTextField regScoreField;
	static int NUMCOLS = 42;
	// Strings for the labels
	static Color gray = new Color(235,235,235);
	static Color defaultColor;
	static String[] staticsourceArray1 = { "User provided" };
	static String[] staticsourceArray2 = { "User provided" };
	static String[] staticsourceArray3 = { "User provided" };

	JTextField goField;
	JTextField extraField;
	JTextField xrefField;
	JTextField categoryIDField;
	JTextField taxonField;
	JTextField evidenceField;
	JButton infoButton = new JButton(Util.createImageIcon("About16.gif"));
	static JFileChooser fc0 = new JFileChooser(new File("Known modules"));
	static JFileChooser fc1 = new JFileChooser(new File("Interaction files"));
	static JFileChooser fc2 = new JFileChooser(new File("Expression data (sample)"));
	

	static Container contentPane;
	JPanel textpanel;
	JPanel panel4;

	JButton transcriptomebutton= new JButton();
	JButton interactionbutton= new JButton();
	JButton knownmodulebutton= new JButton();
	JTextArea runningtext;
	JButton Run = new JButton();
	JButton currentButton = new JButton();
	JButton endSearchButton = new JButton();
	JTextField knownmoduleField, interactionField, trancriptomeField;
	final JSpinner j18, j20, j23;

	
	public Main_interface2() throws FileNotFoundException, IOException {
		super("CLAM v.1.0");

		contentPane = getContentPane();
		BoxLayout layout = new BoxLayout(contentPane, BoxLayout.Y_AXIS);
		contentPane.setLayout(layout);

		JPanel panel1=new JPanel();
		panel1.setBorder(new TitledBorder(null,"Load Data",TitledBorder.LEFT,TitledBorder.TOP));
		panel1.setLayout(new BorderLayout());
		
        /*****************************************************************/
		
		JLabel j5=new JLabel("Load known modules:");
		j5.setBounds(15,20,500,35);
		knownmoduleField=new JTextField(szoptionalFile, JLabel.TRAILING);
		knownmoduleField.setBounds(150,23,400,25);
		knownmoduleField.setBorder(BorderFactory.createLineBorder(Color.black));
		knownmoduleField.setOpaque(true);
		knownmoduleField.setBackground(Color.white);
		knownmodulebutton.setText("Load File");
		knownmodulebutton.setHideActionText(true);
		knownmodulebutton.addActionListener(this);
		knownmodulebutton.setBounds(555,23,90,25);
		
		JLabel j6=new JLabel("Load interaction file:");
		j6.setBounds(15,50,500,35);
		interactionField=new JTextField(szoptionalFile, JLabel.TRAILING);
		interactionField.setBounds(150,53,400,25);
		interactionField.setBorder(BorderFactory.createLineBorder(Color.black));
		interactionField.setOpaque(true);
		interactionField.setBackground(Color.white);
		interactionbutton.setText("Load File");
		interactionbutton.setHideActionText(true);
		interactionbutton.addActionListener(this);
		interactionbutton.setBounds(555,53,90,25);
		
		
		JLabel j7=new JLabel("Load expression files:");
		j7.setBounds(15,80,300,35);
		trancriptomeField=new JTextField(szDataFileDEF, JLabel.TRAILING);
		trancriptomeField.setBounds(150,83,400,25);
		trancriptomeField.setBorder(BorderFactory.createLineBorder(Color.black));
		trancriptomeField.setOpaque(true);
		trancriptomeField.setBackground(Color.white);
		transcriptomebutton.setText("Load File");
		transcriptomebutton.setHideActionText(true);
		transcriptomebutton.addActionListener(this);
		transcriptomebutton.setBounds(555,83,90,25);
		
		/*****************************************************************/

		panel1.setLayout(null); 
		panel1.add(j5);
		panel1.add(j6);
		panel1.add(j7);
		panel1.add(knownmoduleField);
		panel1.add(knownmodulebutton);
		panel1.add(interactionField);
		panel1.add(trancriptomeField);
		panel1.add(interactionbutton);
		panel1.add(transcriptomebutton);
		
        JScrollPane   scrollpanel1   =   new   JScrollPane(panel1, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollpanel1.setBounds(100, 100, 745, 65);
		panel1.setPreferredSize(new Dimension(scrollpanel1.getWidth() - 50, scrollpanel1.getHeight()*2));
		
		/*******************************************************/
		
		JPanel panel3=new JPanel();
		JScrollPane   scrollpanel2   =   new   JScrollPane(panel3, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollpanel2.setBounds(100, 100, 745, 100);
        panel3.setPreferredSize(new Dimension(scrollpanel2.getWidth() - 50, scrollpanel2.getHeight()*2));
		panel3.setBorder(new TitledBorder(null,"Set Parameters",TitledBorder.LEFT,TitledBorder.TOP));
		panel3.setLayout(new BorderLayout());
		
		JPanel panel31=new JPanel();
		panel31.setBorder(new TitledBorder(null,"Data Preprocessing",TitledBorder.LEFT,TitledBorder.TOP));
		panel31.setBounds(5,20,365,175);
		panel31.setLayout(null);
		
		JLabel log=new JLabel("Log normalize data:");
		log.setBounds(15,20,300,35);
	    log1 = new JRadioButton("Yes");
		log2 = new JRadioButton("No");
		log1.setBounds(140, 25, 50, 22);
		log2.setBounds(190, 25, 50, 22);
		if (nnormalizeDEF == 0) {
			log1.setSelected(true);
		} else if (nnormalizeDEF == 1) {
			log2.setSelected(true);
		}
		
		JLabel mm=new JLabel("Max number of missing values:");
		mm.setBounds(15,55,300,35);
		SpinnerNumberModel misscontrol = new SpinnerNumberModel(new Integer(nummissing), new Integer(0), null, new Integer(1));
		j18 = new JSpinner(misscontrol);
		j18.setBounds(210,60,45,22);
		
		JLabel inter=new JLabel("Set missing values as:");
		inter.setBounds(15,95,300,35);
		comboBox1.setBounds(155,100,100,21);
	    if (missinter == 0) {
	    	comboBox1.setSelectedIndex(0);
		} else if(missinter == 1){
			comboBox1.setSelectedIndex(1);
		} else{
			comboBox1.setSelectedIndex(2);
		}
	    
	    JLabel duplicates=new JLabel("Set duplicate genes as:");
		duplicates.setBounds(15,135,300,35);
	    comboBox3.setBounds(155,140,100,21);
	    if (filterduplicates == 0) {
	    	comboBox3.setSelectedIndex(0);
		} else if(filterduplicates == 1){
			comboBox3.setSelectedIndex(1);
		}
		
		panel31.add(log);
		panel31.add(log1);
		panel31.add(log2);
		normGroup.add(log1);
		normGroup.add(log2);
		panel31.add(mm);
		panel31.add(j18);
		panel31.add(inter);
		panel31.add(comboBox1);
		panel31.add(duplicates);
		panel31.add(comboBox3);
		
		/*******************************************************/
	    JPanel panel32=new JPanel();
		panel32.setBorder(new TitledBorder(null,"Module Detecting",TitledBorder.LEFT,TitledBorder.TOP));
		panel32.setBounds(370,20,365,175);
		panel32.setLayout(null);
		
		JLabel minir=new JLabel("Similarity measure:");
		minir.setBounds(15,20,300,35);
		comboBox4.setBounds(145,25,150,21);
	    if (similarity == 0) {
	    	comboBox4.setSelectedIndex(0);
		} else if(similarity == 1){
			comboBox4.setSelectedIndex(1);
		}
		
		JLabel knn=new JLabel("K-nearest neighbors:");
		knn.setBounds(15,55,200,35);
		SpinnerNumberModel con = new SpinnerNumberModel(10, 2, null, 1);
		j20 = new JSpinner(con);
		j20.setBounds(145,60,50,22);
		
		JLabel ming=new JLabel("Min number of module genes:");
		ming.setBounds(15,90,300,35);
		SpinnerNumberModel minigene = new SpinnerNumberModel(new Integer(5), new Integer(1), null, new Integer(1));
		j23 = new JSpinner(minigene);
		j23.setBounds(195,95,50,22);
		
		JLabel j3=new JLabel("Use prior knowledge to assist");
		JLabel j4=new JLabel("in module detection:");
		j3.setBounds(15,123,500,35);
		j4.setBounds(15,138,500,35);
		go1 = new JRadioButton("Yes");
		go2 = new JRadioButton("No");
		go1.setBounds(195, 140, 50, 22);
		go2.setBounds(245, 140, 50, 22);
		if (ngoDEF == 0) {
			go1.setSelected(true);
		} else if (ngoDEF == 1) {
			go2.setSelected(true);
		}
		

		panel32.add(knn);
		panel32.add(minir);
		panel32.add(ming);
		panel32.add(j20);
		panel32.add(comboBox4);
		panel32.add(j23);
		panel32.add(go1);
		panel32.add(go2);
		goGroup.add(go1);
		goGroup.add(go2);	
		panel32.add(j3);
		panel32.add(j4);

		panel3.setLayout(null);
		panel3.add(panel31); 
		panel3.add(panel32);
		/********************************************************/
		textpanel=new JPanel();
		textpanel.setPreferredSize(new Dimension(745, 150));
		textpanel.setBorder(new TitledBorder(null,"Search Gene Modules",TitledBorder.LEFT,TitledBorder.TOP));
		BoxLayout layout2 = new BoxLayout(textpanel, BoxLayout.X_AXIS);
		textpanel.setLayout(layout2);
		
		panel4=new JPanel();
		panel4.setLayout(new BorderLayout());
		runningtext=new JTextArea();
		runningtext.setLineWrap(true);
		JScrollPane sp1=new JScrollPane(runningtext);
		panel4.add(sp1);
		textpanel.add(panel4);

		Run.setText("Run");
		Run.setHideActionText(true);
		Run.addActionListener(this);
		endSearchButton.setText("Stop Searching");
		endSearchButton.setEnabled(false);
		endSearchButton.addActionListener(this);
		infoButton.addActionListener(this);
		infoButton.setBackground(gray);
		
		JPanel panel6=new JPanel();
		panel6.setPreferredSize(new Dimension(745, 50));
		panel6.add(Run);
		panel6.add(endSearchButton);
		panel6.add(infoButton);
		contentPane.add(scrollpanel1);
		contentPane.add(scrollpanel2);
		contentPane.add(textpanel);
		contentPane.add(panel6);
		
		theRepeatList0 = new ListDialog(this, Main_interface2.vRepeatFilesDEF0,
				knownmodulebutton, Main_interface2.fc0);
		theRepeatList2 = new ListDialog(this, Main_interface2.vRepeatFilesDEF2,
				transcriptomebutton, Main_interface2.fc2);
	}

	/**
	 * Checks if the two data sets have the same number of rows, time points,
	 * and the gene name matches.
	 */
	public static void errorcheck(BSM_DataSet theDataSet1,
			BSM_DataSet theOtherSet) {

		if (theDataSet1.numcols != theOtherSet.numcols) {
			throw new IllegalArgumentException(
					"Repeat data set must have same "
							+ "number of columns as original, expecting "
							+ theDataSet1.numcols + " found "
							+ theOtherSet.numcols + " in the repeat");
		} else if (theDataSet1.numrows != theOtherSet.numrows) {
			throw new IllegalArgumentException(
					"Repeat data set must have same "
							+ "number of spots as the original, expecting "
							+ theDataSet1.numrows + " found "
							+ theOtherSet.numrows + " in the repeat");
		} else {
			for (int nrow = 0; nrow < theDataSet1.numrows; nrow++) {
				if (!theDataSet1.genenames[nrow]
						.equals(theOtherSet.genenames[nrow])) {
					throw new IllegalArgumentException("In row " + nrow
							+ " of the repeat set " + "expecting gene symbol "
							+ theDataSet1.genenames[nrow] + " found "
							+ theOtherSet.genenames[nrow]);
				} else if (!theDataSet1.probenames[nrow]
						.equals(theOtherSet.probenames[nrow])) {
					throw new IllegalArgumentException("In row " + nrow
							+ " of the repeat set " + "expecting gene symbol "
							+ theDataSet1.probenames[nrow] + " found "
							+ theOtherSet.probenames[nrow]);
				}
			}
		}
	}

	/**
	 * Checks if origcols and nrepeat cols are the same value, the length of
	 * origgenes and repeatgenes is the same, and the gene names are the same
	 */
	public static void errorcheck(String[] origgenes, String[] repeatgenes,
			int norigcols, int nrepeatcols) {
		if (norigcols != nrepeatcols) {
			throw new IllegalArgumentException(
					"Repeat data set must have same "
							+ "number of columns as original, expecting "
							+ norigcols + " found " + nrepeatcols
							+ " in the repeat");
		} else if (origgenes.length != repeatgenes.length) {
			throw new IllegalArgumentException(
					"Repeat data set must have same "
							+ "number of spots as the original, expecting "
							+ origgenes.length + " found " + repeatgenes.length
							+ " in the repeat");
		} else {
			for (int nrow = 0; nrow < origgenes.length; nrow++) {
				if (!origgenes[nrow].equals(repeatgenes[nrow])) {
					throw new IllegalArgumentException("In row " + nrow
							+ " of the repeat set " + "expecting gene symbol "
							+ origgenes[nrow] + " found " + repeatgenes[nrow]);
				} else if (!origgenes[nrow].equals(repeatgenes[nrow])) {
					throw new IllegalArgumentException("In row " + nrow
							+ " of the repeat set " + "expecting gene symbol "
							+ origgenes[nrow] + " found " + repeatgenes[nrow]);
				}
			}
		}
	}

	private BSM_DataSet buildset(String szexp1val, 
			Integer maxcase, Integer missinter, boolean btakelog) 
			throws Exception {
		BSM_DataSet theDataSet1 = new BSM_DataSet(szexp1val, maxcase, missinter, btakelog);
		theDataSet1 = new BSM_DataSet(theDataSet1.filterMissing());	
		if(filterduplicates==0){
			theDataSet1 = new BSM_DataSet(theDataSet1.maxAndFilterDuplicates());
		}else if(filterduplicates==1){
			theDataSet1 = new BSM_DataSet(theDataSet1.averageAndFilterDuplicates());
		}
		return theDataSet1;	
	}
	
	/**
	 * A control method that handles the response for when the execute button on
	 * the interface is pressed including building the data set, running the
	 * TSMiner modeling procedure, and displaying the results
	 */
	public void clusterscript(Vector<String> knownmodules, String ppifile, Vector<String> theData, 
			String missing, int missinter, boolean btakelog) throws Exception {
		/*
		String aa1 = "D:/Personal issues/CLAM/CLAM tool/CLAM v3.0/Expression data (sample)/Sample data 1.txt";
		String aa2 = "D:/Personal issues/CLAM/CLAM tool/CLAM v3.0/Expression data (sample)/Sample data 2.txt";
		String aa3 = "D:/Personal issues/CLAM/CLAM tool/CLAM v3.0/Expression data (sample)/Sample data 3.txt";
		theData.add(aa1);
		theData.add(aa2);
		theData.add(aa3);
		*/
		/********************************************Expression Data*********************************************/
		List<String> finallist = new ArrayList<String>();
				if (theData != null && theData.size()>0) {
					this.theDataSet = new ArrayList<BSM_DataSet>();
					for(String file:theData) {
						if(new File(file).exists()) {
							BSM_DataSet dataSet = buildset(file, Integer.parseInt(missing), missinter, btakelog);
							dataSet.normalization(finallist);
							theDataSet.add(dataSet);
						} else {
							throw new IllegalArgumentException("The expression data file '" + file+ "' cannot be found.");
						}
					}
				}
				
				int[][] measuredtable = new int[finallist.size()][theData.size()];
				if(theDataSet != null && theDataSet.size()>0) {
					for(int i=0;i<theDataSet.size();i++) {
						BSM_DataSet data = theDataSet.get(i);
						data.filtergenes(finallist, measuredtable, i);
					}
				}
				
				if(theDataSet == null || theDataSet.size()==0) {
					throw new IllegalArgumentException("No expression data file uploaded.");
				}
				System.out.println(finallist.size());
				
				/********************************************Interaction Data********************************************/
				/*
				String save1 = "D:/Personal issues/CLAM/CLAM tool/CLAM v3.0/Known modules/KEGG pathways (hsa).txt";
				String save2 = "D:/Personal issues/CLAM/CLAM tool/CLAM v3.0/Known modules/miRNA-targets (hsa).txt";
				String save3 = "D:/Personal issues/CLAM/CLAM tool/CLAM v3.0/Known modules/TF-targets (hsa).txt";
				knownmodules.add(save1);
				knownmodules.add(save2);
				knownmodules.add(save3);
				ppifile = "D:/Personal issues/CLAM/CLAM tool/CLAM v3.0/Interaction files/Protein interactions (hsa).txt";
				*/
				
				List<RegulatorBindingData> KnownModules = new ArrayList<RegulatorBindingData>();
				int totalregs = 0;
				if(knownmodules.size()>0) {
					try {
						for(String file:knownmodules) {
							RegulatorBindingData pathwayData = new RegulatorBindingData(file,
									finallist, SZDELIM, null, null);
							pathwayData.deDuplication();
							KnownModules.add(pathwayData);
							totalregs += pathwayData.numberRegs;
						}
					} catch (Exception e) {
						throw new IllegalArgumentException("Read known module file error, please check the knwon module files.");
				    }
				}
				
				PPIInteractionData ppiData = null;
				String[] totalregnames = null;
				int[][] reg2geneindex = null;
				int[][] gene2regindex = null;
				
				if(busego) {
					totalregnames = new String[totalregs];
					reg2geneindex = new int[totalregs][];
					gene2regindex = new int[finallist.size()][];
					HashMap<String, Integer> reg2index = new HashMap<String, Integer>();
					int count=0;
					for(RegulatorBindingData pathwayData:KnownModules) {
						for(int i=0;i<pathwayData.numberRegs;i++) {
							totalregnames[count] = pathwayData.regNames[i];
							reg2index.put(pathwayData.regNames[i], count);
							reg2geneindex[count] = pathwayData.reg2GeneDedupliIndex[i];
							count++;
						}
					}
					for(RegulatorBindingData pathwayData:KnownModules) {
						for(int i=0;i<pathwayData.numberGenes;i++) {
							if(pathwayData.gene2RegDedupliIndex[i] != null && pathwayData.gene2RegDedupliIndex[i].length>0) {
								if(gene2regindex[i] == null) {
									gene2regindex[i] = new int[pathwayData.gene2RegDedupliIndex[i].length];
									for(int j=0;j<pathwayData.gene2RegDedupliIndex[i].length;j++) {
										gene2regindex[i][j] = reg2index.get(pathwayData.regNames[pathwayData.gene2RegDedupliIndex[i][j]]);
									}
								} else {
									int[] vv = new int[gene2regindex[i].length+pathwayData.gene2RegDedupliIndex[i].length];
									for(int j=0;j<gene2regindex[i].length;j++) {
										vv[j] = gene2regindex[i][j];
									}
									for(int j=0;j<pathwayData.gene2RegDedupliIndex[i].length;j++) {
										vv[gene2regindex[i].length+j] = reg2index.get(pathwayData.regNames
												[pathwayData.gene2RegDedupliIndex[i][j]]);
									}
									gene2regindex[i] = vv;
								}
							}
						}
					}
					if(ppifile != null && ppifile.length()>0 && !ppifile.equalsIgnoreCase(" Optional")) {
						try {
							ppiData = new PPIInteractionData(ppifile, finallist);
						}catch (Exception e) {
			            	throw new IllegalArgumentException("Read file "+ ppifile +" error.");
			            }
					}
				}
				
				/*****************************************************************************************************/
				
				bendsearch = false;
				try {
					LearningCase caseModules = new LearningCase(finallist, measuredtable, theDataSet, KnownModules, 
							ppiData, totalregnames, reg2geneindex, gene2regindex, busego, szknn, similarity, 
							mingenes, runningtext, endSearchButton);
		            } catch (Exception e) {
						e.printStackTrace();
				}
		
		System.exit(0);
		long e1 = System.currentTimeMillis();
		System.out.println("Time: " + (e1 - s1) + "ms");
	}
	
	
	public List<List<String>> readFile(String filename){
        List<List<String>> value = new ArrayList<List<String>>();
		try {
			BufferedReader bufferedReader1 = new BufferedReader(new FileReader(filename));   
			String lineTxt = null;
            while((lineTxt = bufferedReader1.readLine()) != null){
                value.add(Arrays.asList(lineTxt.split("\t")));
            }
            bufferedReader1.close();  
                
    } catch (Exception e) {
    	throw new IllegalArgumentException("Read saved model file error.");
    }
		return value;
	}

	/**
	 * define the button methods
	 */
	public void actionPerformed(ActionEvent e) {
		Object esource = e.getSource();

		if (esource == endSearchButton) {
			bendsearch = true;
			runningtext.append("End Search Requested. Search Will End Soon..."+"\n");
			runningtext.paintImmediately(runningtext.getBounds());
			endSearchButton.setEnabled(false);
		} else if (esource == transcriptomebutton) {
			theRepeatList2.setLocation(this.getX() + 75, this.getY() + 100);
			theRepeatList2.setVisible(true);
		}  else if (esource == knownmodulebutton) {
			theRepeatList0.setLocation(this.getX() + 75, this.getY() + 100);
			theRepeatList0.setVisible(true);
		}  else if (esource == interactionbutton) {
			int returnVal = fc1.showOpenDialog(this);
			if (returnVal == JFileChooser.APPROVE_OPTION) {
				File file = fc1.getSelectedFile();
				interactionField.setText(file.getAbsolutePath());
			}
		} else if (esource == Run) {
			s1 = System.currentTimeMillis();
			busego = go1.isSelected();
			expressionfiles = theRepeatList2.data;
			knownmodulefiles = theRepeatList0.data;
			szinteractionval = interactionField.getText();
			missing = j18.getValue().toString();
			szknn = j20.getValue().toString();
			mingenes = j23.getValue().toString();
			
			int temp1=comboBox4.getSelectedIndex();
			if(temp1==0){
				similarity=0;
			}else if(temp1==1){
				similarity=1;
			}else if(temp1==2){
				similarity=2;
			}else if(temp1==3){
				similarity=3;
			}else if(temp1==4){
				similarity=4;
			}
			
			int temp2=comboBox1.getSelectedIndex();
			if(temp2==0){
				missinter=0;
			}else if(temp2==1){
				missinter=1;
			}else{
				missinter=2;
			}
			int temp3=comboBox3.getSelectedIndex();
			if(temp3==0){
				filterduplicates=0;
			}else if(temp3==1){
				filterduplicates=1;
			}
			
			btakelog = log1.isSelected();

			this.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));

			final JFrame fframe = this;

			Runnable clusterrun = new Runnable() {
				public void run() {
					Run.setEnabled(false);
					try {
						clusterscript(knownmodulefiles, szinteractionval, expressionfiles, missing, missinter, btakelog);
					} catch (IllegalArgumentException iex) {
						final IllegalArgumentException fiex = iex;
						iex.printStackTrace(System.out);

						javax.swing.SwingUtilities.invokeLater(new Runnable() {
							public void run() {
								JOptionPane.showMessageDialog(fframe, fiex.getMessage(), "Error",JOptionPane.ERROR_MESSAGE);
							}
						});
					} catch (Exception ex) {
						final Exception fex = ex;

						javax.swing.SwingUtilities.invokeLater(new Runnable() {
							public void run() {
								JOptionPane.showMessageDialog(fframe, fex.toString(), "Exception thrown",JOptionPane.ERROR_MESSAGE);
								fex.printStackTrace(System.out);
							}
						});
					}
					Run.setEnabled(true); 
					fframe.setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
				}
			};
			(new Thread(clusterrun)).start();
			
		} else if (esource == infoButton) {//The Information and Help button
			String szMessage = "This is version 1.0.0 of CLAM.\n\n"
					+ "The CREAM is available under a GPL v3.0 license.\n"
					+ "Any questions or bugs found should "
					+ "be emailed to free1234hm@163.com.";

			Util.renderDialog(this, szMessage, 50, 100, "Information");
		}
	}
	
	/**
	 * Create the GUI and show it. For thread safety, this method should be
	 * invoked from the event-dispatching thread.
	 */
	private static void createAndShowGUI() throws FileNotFoundException,IOException {
		// Make sure we have nice window decorations.
		// JFrame.setDefaultLookAndFeelDecorated(true);
		// Create and set up the window.
		JFrame frame = new Main_interface2();
		frame.setLocation(10, 25);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		// Display the window.
		frame.pack();
		frame.setVisible(true);
	}

	/**
	 * The main method which when executed will have the input interface created
	 */
	public static void main(String[] args) throws Exception {
			javax.swing.SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					try {
						createAndShowGUI();
					} catch (FileNotFoundException ex) {
						ex.printStackTrace(System.out);
					} catch (IOException ex) {
						ex.printStackTrace(System.out);
					}
				}
			});
	}
}
