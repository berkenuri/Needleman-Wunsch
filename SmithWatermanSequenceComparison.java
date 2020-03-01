import java.io.*;
import java.lang.*;
import java.util.*;

import Jama.*;


public class SmithWatermanSequenceComparison{

	static double[] gap_penalty = null;
	static double[][] weight_matrix = null;
	static Matrix score = null;

	public static void main(String[] args) throws Exception {

				// no config file specified
				if (args.length == 0){
            System.out.println("Usage: Smith-WatermanSequenceComparison <configFile>!");
            System.exit(-1);
        }

        // Load config file
        Properties properties = new Properties();

        InputStream instream = null;

        try{
            instream = new FileInputStream(args[0]);
        }
        catch(Exception e) {
            System.out.println("Invalid configuration file");
            System.exit(-1);
        }

				// load in config file
				properties.load(instream);

        String sequenceInputFile = properties.getProperty("sequenceInputFile");

        if (sequenceInputFile.equals(null)){
            System.out.println("Invalid sequence input file");
            System.exit(-1);
        }

        String sequenceA_number = properties.getProperty("sequenceA");
        String sequenceB_number = properties.getProperty("sequenceB");

        String[] sequences = parseSequenceFile(sequenceInputFile, sequenceA_number, sequenceB_number);
        String sequenceA = sequences[0];
        String sequenceB = sequences[1];

        String weightMatrixFile = properties.getProperty("weightMatrixFile");
        String gapPenaltyFile = properties.getProperty("gapPenaltyFile");

        int maxColumns = 0;

        try{
            maxColumns = Integer.parseInt(properties.getProperty("maxColumns"));
        }
        catch(Exception e){
            System.out.println("Invalid max columns value");
            System.exit(-1);
        }

        boolean writeToFile = false;

				try{
        		writeToFile= Boolean.parseBoolean(properties.getProperty("writeToFile"));
        }
        catch(Exception e){
            System.out.println("Improper boolean value");
            System.exit(-1);
        }
        String outputFile = properties.getProperty("outputFile");

				// calculate matrices and gap penalty
				gap_penalty = parseGapPenalty(gapPenaltyFile);
        weight_matrix = parseWeightMatrix(weightMatrixFile);
        score = calculateScoreMatrix(sequenceA, sequenceB);

        String[] alignments = findAlignment(score, sequenceA, sequenceB);

				if (writeToFile) {
					writeFile(outputFile, alignments, maxColumns);
				}
				else {

        	System.out.println();
        	System.out.println();

        	if (alignments[0].length() <= maxColumns){
            System.out.println(alignments[0]);
            System.out.println(alignments[1]);
        }
        else {
            int i = 0;
            int divisions = alignments[0].length()/maxColumns;

            while (i < divisions ){
                System.out.println(alignments[0].substring(i*maxColumns, (i+1)*maxColumns));
                System.out.println(alignments[1].substring(i*maxColumns, (i+1)*maxColumns));
                System.out.println();

                i++;
            }

            if(i*maxColumns < alignments[0].length()){
                System.out.println(alignments[0].substring(i * maxColumns));
                System.out.println(alignments[1].substring(i * maxColumns));
            }
        }

				System.out.println();
			}
}

	public static void print2D(double[][] arr)
    {
        // Loop through all rows
        for (double[] row : arr)

            // converting each row as String
            // and then printing in a separate line
            System.out.println(Arrays.toString(row));
    }

	// get the index number depending on the nucleotide type
	public static int getVal(char element) {

		switch (element)
		{
				case 'A': return 0;
				case 'C': return 1;
				case 'G': return 2;
				case 'T': return 3;
		}

		return -1;
	}

	// write the alignments on the specified file
	public static void writeFile(String file_name, String[] alignments, int maxColumns) {

		File file = new File(file_name);

		try {
      FileWriter fr = new FileWriter(file);
      fr.write("\n");
      fr.write("\n");

			// write all of the sequence if the sequence is shorter than the max column number
			if (alignments[0].length() <= maxColumns){
					fr.write(alignments[0] + "\n");
					fr.write(alignments[1] + "\n");
			}
			// divide the sequence into multiple lines if its longer than the max column number
			else{
					int i = 0;
					int divisions = alignments[0].length()/maxColumns;

					while (i < divisions ) {
							fr.write(alignments[0].substring(i*maxColumns, (i+1)*maxColumns) + "\n");
							fr.write(alignments[1].substring(i*maxColumns, (i+1)*maxColumns) + "\n");
							fr.write("\n");

							i++;
					}

					if(i*maxColumns < alignments[0].length()) {
							fr.write(alignments[0].substring(i * maxColumns) + "\n");
							fr.write(alignments[1].substring(i * maxColumns) + "\n");
					}
			}

			fr.close();

    } catch (Exception e) {
    	System.out.println("Improper write file");
        System.exit(-1);
    }

	}

	// parses the specific sequences in a file depending on their sequence number
	public static String[] parseSequenceFile(String file_name, String seqNumber_A, String seqNumber_B) {

		String sequenceA = "";
		String sequenceB = "";
		String[] sequences = new String[2];

		Scanner sc;
		try {

			sc = new Scanner(new File(file_name));
			String line  = "" ;

			// start reading the file
			while(sc.hasNext()) {

				line = sc.nextLine();

				if (line.contains("sequence")) {

					String[] temp = line.trim().split(":");

					// checks the sequence number is equal to the parameter sequence number
					if (seqNumber_A.equals(temp[1])) {
						sequenceA = temp[2]; //adds the sequence part after its number

						String seqA_line = "";

						while(sc.hasNext()) {
							seqA_line = sc.nextLine();

							// if the sequence is still going on continues adding the next lines
							if (!seqA_line.isEmpty() )
								sequenceA += seqA_line;
							else
								break; //stops adding to sequence if the line is empty

						}
					}

					// checks the sequence number is equal to the second parameter sequence number
					if (seqNumber_B.equals(temp[1])) {
						sequenceB = temp[2]; //adds the sequence part after its number

						String seqB_line = "";

						while(sc.hasNext()) {
							seqB_line = sc.nextLine();

							// if the sequence is still going on continues adding the next lines
							if (!seqB_line.isEmpty())
								sequenceB += seqB_line;
							else
								break; //stops adding to sequence if the line is empty

						}
					}

				}

			}

			// get rid of the X at the end of the sequence
      if (sequenceA.charAt(sequenceA.length()-1) == 'X')
        	sequenceA = sequenceA.substring(0, sequenceA.length() -1);

      for (int i = 0; i < sequenceA.length(); i++) {
          if (getVal(sequenceA.charAt(i)) == -1){
              System.out.println("Invalid characters in sequence A");
              System.exit(-1);
          }
      }


			sequences[0] = sequenceA;

			// get rid of the X at the end of the sequence
      if (sequenceB.charAt(sequenceB.length()-1) == 'X')
          sequenceB = sequenceB.substring(0, sequenceB.length() -1);

      for (int i = 0; i < sequenceB.length(); i++){
          if (getVal(sequenceB.charAt(i)) == -1){
              System.out.println("Invalid characters in sequence B");
              System.exit(-1);
          }
      }

      sequences[1] = sequenceB;

			sc.close();

			return sequences;

		} catch (Exception e) {
						System.out.println("Improper entry in sequence input file");
						System.exit(-1);
		}

		return null;
	}

	// parses the file and returns an array of gap penalties
	public static double[] parseGapPenalty(String file_name) {

		Scanner sc;

		try {

			sc = new Scanner(new File(file_name));

			int size = 0;

			String input = "" ;
			String line = "" ;

			// start reading the file
			while(sc.hasNext()) {

				line = sc.nextLine();

				// ignore the empty lines and comments
				if (line.length() == 0 || line.contains("#"))
					continue;

				input += line; // add the other lines into a string
	    }

			// split the string and store the elements in an array
			String[] temp = input.trim().split(" ");
			size = temp.length;

			double[] penalty_array = new double [size];

			// store the elements in the array into a double array of gap penalties
			for (int i = 0; i < size; i++) {
				penalty_array[i] = Double.parseDouble(temp[i]);
			}

	    sc.close();

	    return penalty_array;

		} catch (Exception e) {
            System.out.println("Improper penalty file");
            System.exit(-1);
		}

		return null;
	}

	//parses the file and returns a 2d array of the weight matrix
	public static double[][] parseWeightMatrix(String file_name) {

		Scanner sc;

		try {

			sc = new Scanner(new File(file_name));

			int row = 0;
			int col = 0;
			String input = "" ;
			String line = "" ;

			// start reading the file
			while(sc.hasNext()) {

				line = sc.nextLine();

				// ignore the empty lines and comments
				if (line.length() == 0 || line.contains("#"))
					continue;

        String[] temp = line.trim().split(" ");

				input += " "+ line; // store the elements in a string
				row++; //count the rows
	    }


			String[] temp = line.trim().split(" ");
			col = temp.length; // count the columns

			// if the matrix is asymmetrix, exit
      if (row != col || row != 4){
          System.out.println("Invalid matrix input file");
          System.exit(-1);
      }

			double[][] w_matrix = new double [row][col];

			// seperate the elements and store them in an array
			String[] temp2 = input.trim().split(" ");
			int count = 0;

			for (int i = 0; i < row; i++) {
				for (int j = 0; j < col; j++) {

					// assign the elements in the array into the 2d weight matrix
					w_matrix[i][j] = Double.parseDouble(temp2[count]);
					count++;
	      }
	    }

	    sc.close();

	    return w_matrix;

		} catch (Exception e) {
			System.out.println("Improper Weight Matrix file");
            System.exit(-1);
		}

		return null;

	}

	// Get maximum value of the choices
	public static double getMax( double[] choices )
	{
		// start with the first value
		double maximum = choices[0];

		for (int i=1; i < choices.length; i++) {
	        if (choices[i] > maximum) {
	            maximum = choices[i];   // new maximum
	        }
	    }

		return maximum;
	}

	// calculates the score matrix of the parameter sequences
	public static Matrix calculateScoreMatrix ( String seq1, String seq2  ) {

		double choices[] = new double[3];

		//matriix dimensions are 1 more then the sequence lengths
		score = new Matrix(seq1.length() + 1, seq2.length() + 1);

    score.set(0, 0, 0.0);

		// first column and row are summation of gap penalties
    for(int x = 1; x < seq1.length() + 1; x++)
			score.set(x, 0, score.get(x-1, 0) + gap_penalty[getVal(seq1.charAt(x-1))]);

		for(int y = 1; y < seq2.length() + 1; y++)
			score.set(0, y, score.get(0, y-1) + gap_penalty[getVal(seq2.charAt(y-1))]);

		// fill the rest of the matrix
		for(int x = 1; x <= seq1.length(); x++) {
			for(int y = 1; y <= seq2.length(); y++) {

				// if the last sequence elements match
        choices[0] = score.get(x - 1, y - 1) +
						weight_matrix[ getVal(seq1.charAt(x-1)) ][ getVal(seq2.charAt(y-1)) ];
				// if the last sequence elements don't match and there is a gap at sequence 1
				choices[1] = score.get(x - 1, y) + gap_penalty[ getVal(seq1.charAt(x-1)) ];
				// if the last sequence elements don't match and there is a gap at sequence 2
				choices[2] = score.get(x, y - 1) + gap_penalty[ getVal(seq2.charAt(y-1)) ];

				// get the max of the choices
				score.set(x, y, getMax(choices));
			}
		}

		return score;
	}

    public static String[] findAlignment(Matrix fMatrix, String A, String B){
        String alignmentA = "";
        String alignmentB = "";
        int i = A.length(); // start the index at the end
        int j = B.length(); // start the index at the end


        while (i > 0 || j > 0){
            char A_i;
            char B_j;

            if (i == 0)
                A_i = '-'; // reached the beginnıng of A
            else
                A_i = A.charAt(i-1);

            if (j == 0)
                B_j = '-'; // reached the beginnıng of B
            else
                B_j = B.charAt(j-1);

						// Go up first
						if ((i > 0) &&
						((fMatrix.get(i, j) == (fMatrix.get(i-1, j) + gap_penalty[getVal(A_i)])))){
                alignmentA = A_i + alignmentA;
                alignmentB = "-" + alignmentB;

                i--;
            }
						// Diagonal second (match)
            else if ((i> 0) && (j > 0) &&
						((fMatrix.get(i, j) == (fMatrix.get(i-1, j-1) + weight_matrix[getVal(A_i)][getVal(B_j)])) )){
                alignmentA = A_i + alignmentA;
                alignmentB = B_j + alignmentB;

                i--;
                j--;
            }
						// back third
            if (j > 0 && (fMatrix.get(i, j) == (fMatrix.get(i, j-1) + gap_penalty[getVal(B_j)]))){
                alignmentA = "-" + alignmentA;
                alignmentB = B_j + alignmentB;

                j--;
            }

        }

        String[] alignments = {alignmentA, alignmentB};

        return alignments;
    }

}
