#!/usr/bin/env python3
# =====================================================================
# Title             : Christian_Tran_R11641653_assignment_4.py
# Description       : This program reads command line arguments for input FASTA files, aligns multiple sequences from score matrix and outputs a FASTA style alignment using Feng Doolittle
# Author            : Christian Tran (R11641653)
# Date              : 10/19/2023
# Usage             : python3 Christian_Tran_R11641653_assignment_4.py -i <input file path for .fna> -o <output file path for .fna> -s <path for score matrix>
# Notes             : This program requires compatibility for Python 3.11.3 and has to accept command line arguments -i and -o respectively for HPCC grading script
# Python Version    : Python 3.11
# =====================================================================

import argparse #import argparse for command line arguments 
from operator import itemgetter #import the itemgetter to specify which element in the list of lists to base the sorting off of
import math #used for mathematical operations (like getting natural log in this case)
import numpy as np #used for matrix manipulation

# Function Name : readFastA
# Parameters    : 
#                   inputFile (File type) - Accepts a input file (in this case .fna) to read and sort
# Returns       : list_of_lists[] - mainly used to write the sorted sequences back to another .fna file
# Purpose       : Reads the given inputFile (from the -i (inputpath) command), sorts the sequences based on the length from longest to shortest, then returns the sorted list to main
def readFastA(inputFile):
    list_of_lists = [] #Combined list of all the elements together
    title_sequence = [] #list for the titles of each sequence starting with >
    full_sequence_list = [] #list for the full combined sequence
    sequence = "" #temp string for reading full sequences
    j = 0 #Counter for the sequence element

    #for loop to read the inputFile given and appends them to their respective list depending on if its a sequence title, or part of the actual sequence
    for line in inputFile: #For each line in the file given
        if line[0] == ">": #if the first character in the line starts with a > then:
            title_sequence.append(line.strip()) #Add the title of the sequence to the title_sequence list
            if sequence: #if the sequence HAS something, then add it to the full_sequence_list before getting to the next > symbol
                full_sequence_list.append(sequence)
            sequence = "" #empty the list whenever we reach a new sequence that starts with the character >
        else: #if the first character in a line doesn't start with >, then we can assume that its part of a sequence
            sequence = sequence + line.strip() #until we get to the next sequence title, add lines to make the full sequence
    if sequence: #at the end of the file, if there still is something inside this sequence variable empty it into the full_sequence_list
        full_sequence_list.append(sequence)
        sequence = "" #empty the list
    
    #for loop to create the list of lists to compare and sort the elements
    for element in title_sequence: 
        list_of_lists.append([element, full_sequence_list[j], len(full_sequence_list[j])])#gets the length of the sequence
        j = j + 1 #adds 1 to the counter to iterate through the full_sequence_list
    
    #Sorts the list of lists based off of the length of the sequence string value in list index 2
    list_of_lists.sort(key = itemgetter(2), reverse = True) 

    #Removes the length value from the list of lists
    for h in list_of_lists: 
        h.pop(2) #pops element 2 out of list

    return list_of_lists

#   Function Name: readScoreMatrix
#   Parameters: 
#               ScoreMatrixFile - file to read score matrix from
#   Returns: A set of tuples for keys (letter pairs) assigned to their corresponding score
#   Purpose: Reads the score matrix file using nested for loops and yield calls to generate a set of tuples to be turned into a dictionary
def readScoreMatrix(ScoreMatrixFile):
  #note: Use .rstrip() (gets rid of all trailing characters) and .split() (makes a list based on the separated elements)
  #List comprehension: A way to apply the given condition to all of the elements in a list. In this case, we want all of the lines in the file
  #                    to have no trailing characters, but also in it's own list, so a list for each row with each score as elements of the list.
  inputFile_lines = (line.rstrip().split() for line in ScoreMatrixFile) #Creates a list of the rows/lines in the file with list comprehension
  
  Score_header = next(inputFile_lines) #iterate the rows once to get the next row, which in this case is the header letters (A, B, C, etc.)

  #for each row in the list of lines from the file that are .rstripped and .split
  for row in inputFile_lines:
    row_letter = row[0] #get the index of the letter
    #Set zip(Score_header, to current whole row) 
    #The zip takes the individual score header (ex: A) and "zips" it in a tuple pairs, so in format (A, 1)
    #Goes through whole list line/row by row and creates tuples of tuples based on the row_letter and col_letter, and goes through all possibilities and assigns them at the current score position
    for col_letter, Score_num in zip(Score_header, row[1:]):
      yield ((row_letter, col_letter), Score_num) #yield a set of tuples as the two letter key values and their score to be used in dict() to turn into dictionary


#   Function Name: needlemanWunsch
#   Parameters: 
#               seq1 - The first sequence to align
#               seq2 - The second sequence to align with the first sequence
#               score_dict - The score dictionary from the score matrix file; used to look up matching score values
#   Returns: aligned sequences and the optimal alignment score (score to the bottom right corner) as a tuple
#   Purpose: Implements the needleman-wunsch algorithm to calculate the optimal score and alignment for 2 given sequences
def needlemanWunsch(seq1, seq2, score_dict):
    
    #initialize empty matrix
    Initial_Matrix = [] 

    #Initialize Traceback matrix
    Traceback_Matrix = []

    #Alignment score variable
    optimal_alignment_score = 0

    #finished alignment variable strings
    aligned_seq1 = ""
    aligned_seq2 = ""
    
    #Get gap penalty from score_dict
    #Note: originally I thought that the gap score came from element key tuple ("-", "-"), but that seems wrong so I am testing it with the letter tuple for the new gap value which might be correct (letter, "-") or ("-", letter)
    #note: getting a dict element from the key is a string, make sure to convert it to an integer for calculation purposes
    #gap_pen = -int(score_dict.get(("-","-"))) #Remember to turn the gap score negative per the algorithm requirment
    #gap_pen = -2 #Test gap score
    #gap_pen = int(score_dict.get()) #Unfinished test gap score
    #Note: changed gap score to (letter and "-") or vice versa, and commented out the ("-", "-") gap score to see if this works.

    #print("CURRENT GAP PENALTY: "+str(gap_pen)) #Test to display current gap penalty

    #Fill empty matrix with zeroes to start off, also fill the traceback matrix with lists and empty strings
    for row in range(len(seq1) + 1): #each row made from the amount of len in seq 1
        Initial_Matrix.append([]) #add empty list inside list to make a list of list (matrix)
        Traceback_Matrix.append([]) #add empty lists to traceback matrix also
        # print(Initial_Matrix)  #used to test and see where the added lists were going
        for col in range(len(seq2) + 1): #each column made from the amount of len in seq 2
          #  print(str(Initial_Matrix[-1])) #used to test and track where the zeroes were going
            Initial_Matrix[-1].append(0) #add zero to each column in each row
            Traceback_Matrix[-1].append("?") #add a character random into the traceback
            #note: we append to [-1] because we want to add the zeroes to the list we JUST added in the row for loop

    #Fill out the first row and first column with default needleman_wunsch values by taking multiples of the gap score
    #From lecture 05 needleman-wunsch algorithm slides:
    for row_iterator in range(len(seq1) + 1): #Fill out score matrix rows with multiples of gap penalty
        Initial_Matrix[row_iterator][0] = row_iterator * (int(score_dict.get((seq1[row_iterator-1], "-"))))

    for row_iterator in range(len(seq1) + 1): #Fill out traceback matrix rows with "up"
        Traceback_Matrix[row_iterator][0] = "up"

    for col_iterator in range(len(seq2) + 1): #Fill out score matrix columns with multiples of gap penalty
        Initial_Matrix[0][col_iterator] = col_iterator * (int(score_dict.get(("-", seq2[col_iterator-1]))))    

    for col_iterator in range(len(seq2) + 1): #Fill out traceback matrix rows with "left"
        Traceback_Matrix[0][col_iterator] = "left"

    #Note: fill out the Traceback[0,0] with DONE to keep accuracy and to mark where the traceback ends
    Traceback_Matrix[0][0] = "DONE"

    #use the needleman wunsch rules to find the greatest value of the score matrix and fill out all cells
    for row_iterator in range(1, len(seq1) + 1): #Start at range 1 to the end of seq1 length because the first row and column do NOT need to be recalculated due to having default values
        for col_iterator in range(1, len(seq2) + 1): #Start at range 1 because at position (0,0) there is a 0, so theres a small shift in the table by 1 for the default values of the starting score matrix
            #Algorithm pseudocode is on Slides 05 - pairwise alignment slide 26 for Needleman-wunsch rules

            diag_score = Initial_Matrix[row_iterator - 1][col_iterator-1] + int(score_dict.get((seq1[row_iterator-1],seq2[col_iterator-1]))) #add diagonal score from score matrix
            left_score = Initial_Matrix[row_iterator][col_iterator-1] - abs(int(score_dict.get(("-", seq2[col_iterator-1])))) #its initially negative, but to follow the subtract rule, we absolute value it back to positive
            up_score = Initial_Matrix[row_iterator - 1][col_iterator] - abs(int(score_dict.get((seq1[row_iterator-1], "-")))) #same thing with up calculation
            #note: I mixed up the left and up scores because I wasn't paying attention. This should be fixed after I changed it (10/10/23)

            max_value = max(diag_score, up_score, left_score) #finds the max/largest value out of the 3 scores, compares scores for max one 

            Initial_Matrix[row_iterator][col_iterator] = max_value #adds the max/largest value out of the 3 and stores it as the main score in this current position

            #Input Traceback directions into traceback matrix
            if(diag_score == max_value):
                Traceback_Matrix[row_iterator][col_iterator] = "diag"
            elif(up_score == max_value):
                Traceback_Matrix[row_iterator][col_iterator] = "up"
            elif(left_score == max_value):
                Traceback_Matrix[row_iterator][col_iterator] = "left"

            #Test display for the score tracking
            #print("DIAG SCORE: "+str(diag_score))
            #print("LEFT SCORE: "+str(left_score))
            #print("UP SCORE:   "+str(up_score))
            #print("MAX SCORE FOUND: "+str(max_value))
            #print("ROW: "+str(row_iterator))
            #print("COL: "+str(col_iterator))
            #print("GAP SCORE UP: "+str(score_dict.get(seq1[row_iterator-1], "-")))
            #print("GAP SCORE LEFT: "+str(score_dict.get("-", seq2[col_iterator-1])))
            #print("\n\n\n")

            

    #Test display for both score and traceback matrix:
    #displayMatrix(Initial_Matrix)
    #print("\n\n")
    #displayMatrix(Traceback_Matrix)

    #Assign the optimal alignment score to the lower right most cell score
    optimal_alignment_score = Initial_Matrix[len(seq1)][len(seq2)]


    #Iterate through the traceback matrix and combine sequence depending on traceback direction
    #Note: Changelog (10/12/23): 
    #I originally used a nested for loop to go through the traceback matrix, which did not work.
    #I went through and saw that I got the score and the traceback matrix right, its just that
    #I wasn't traversing through the traceback matrix properly. Whenever there was a diag, up, or left
    #operation, the for loop would decrement (since I started in the bottom right) along with the 
    #movement operations, adding an extra -1 to the row and column iterator which made the traversal
    #incorrect. I changed it to a while loop to cause the loop to not decrement the row and column iterator
    #like a for loop, and to ONLY CHANGE POSITION WHEN IT CHECKS SPECIFIC DIRECTION. No other thing can
    #cause the row and column iterators to change value except those.

    #Another note: I do NOT have to redefine/reassign col_iterator and row_iterator because they both
    #have a starting value of len(row_iterator) and len(col_iterator) (they both start at the end, bottom right) 
    #due to the previous nested for loop from the score matrix. I am simply just reusing those length values 
    #and going backwards from bottom right to top left. I verified this method works on (10/12/23)
    while col_iterator > -1 and row_iterator > -1:
        #print("INNER LOOP: Curr row: "+str(row_iterator))
        #print("INNER LOOP: Curr col: "+str(col_iterator))
        #Start at bottom right and go to DONE 
        if(Traceback_Matrix[row_iterator][col_iterator] == "DONE"): #Check if reached DONE, then exit nested for loops
           # print("DONE DETECTED")
            break #Exit inner for loop
                
        elif(Traceback_Matrix[row_iterator][col_iterator] == "diag"): #if the current traceback is diagonal, then align both sequence letters together and move diagonal
            aligned_seq1 += seq1[row_iterator-1] #align seq1 letter with sequence 2
            aligned_seq2 += seq2[col_iterator-1]
            row_iterator -= 1 #move diagonally backwards since we started at the bottom right and moving to the upper left
            col_iterator -= 1
          #  print("DIAG CHANGE: Curr row: "+str(row_iterator))
          #  print("DIAG CHANGE: Curr col: "+str(col_iterator))
          #  print("diag DETECTED")
            continue

        elif(Traceback_Matrix[row_iterator][col_iterator] == "up"): #if the current traceback is up, then align just seq1 and gap the up letter from seq2      
            aligned_seq1 += seq1[row_iterator-1] #align just seq1 current letter
            aligned_seq2 += "-" #gap the up from seq2
            row_iterator -= 1 #move up on the traceback matrix
            #col_iterator += 1 #stay on the same column
          #  print("UP CHANGE: Curr row: "+str(row_iterator))
          #  print("UP CHANGE: Curr col: "+str(col_iterator))
          #  print("up DETECTED")
            continue
               
        elif(Traceback_Matrix[row_iterator][col_iterator] =="left"): #if the current traceback is left, then align just seq2 and gap the left letter from seq1    
            aligned_seq1 += "-" #gap the left from seq1
            aligned_seq2 += seq2[col_iterator-1] #align just seq2 current letter
            col_iterator -= 1 #move left on the traceback matrix
            #row_iterator += 1 #stay on same row
           # print("LEFT CHANGE: Curr row: "+str(row_iterator))
           # print("LEFT CHANGE: Curr col: "+str(col_iterator))
           # print("left DETECTED")
            continue
              
    #reverse alignment strings to get the correct alignment since it aligns starting at the bottom right like we do when we're drawing it out, from right to left
    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]

    return [aligned_seq1, aligned_seq2, optimal_alignment_score] #return finished aligned sequences as a LIST because I cannot change values in a tuple



#   Function Name: displayMatrix
#   Parameteres: 
#               matrix - The input matrix to be displayed
#   Returns: nothing
#   Purpose: Used mainly for test purposes, this operation displays a given list of lists in matrix style
def displayMatrix(matrix):
    for i in matrix:
        print(str(i))


#   Function Name: convertAlignToDist
#   Parameters: 
#               AlignmentScore - Given an alignment score, convert this score into a distance score and return the value as a distance value
#               AlignmentScoreMatrix - Given a diagonal matrix of alignment scores between all pairs of sequences, use this to get the max/min alignment score for calculations
#      
#   Returns: distanceScore- return the distance score of the converted Alignment Score
#   Purpose: Used to convert a given alignment score to a distance score, then returns the distance score value
def convertAlignToDist(AlignmentScore, subMatrix, AlignmentLen):
    #Variable to represent the distance score
    distance_score = 0

    #   1. Normalize the score first
    #           S(n) = (S(a) - S(min)) / (S(max) - S(min))

    lowest_score_sub = subMatrix[min(subMatrix, key = subMatrix.get)] #Finds the LOWEST score possible in the subMatrix (nucleotide.mtx, BLOSUM50.mtx, etc.)

    highest_score_sub = subMatrix[max(subMatrix, key = subMatrix.get)] #Finds the HIGHEST score possible in the subMatrix (nucleotide.mtx, BLOSUM50.mtx, etc.)
    
    s_min = lowest_score_sub * AlignmentLen #S(min) = lowest possible score in sub matrix and multiply it by the length of the current alignment

    s_max = highest_score_sub * AlignmentLen #S(max) = highest possible score in sub matrix and multiply it by the length of the current alignment

    s_normalized = (AlignmentScore - s_min) / (s_max - s_min) #complete normalization equation as shown above
    
    #   2. After normalizing, convert to distance using natural log and scaling factor (D = -k * ln(S(n)))
    k = 10  #scaling factor k in the equation

    distance_score = -(k) * math.log(s_normalized) #Use the convert to distance score equation from Slide deck 06 slide 21

    #After these 2 steps, the given alignment scores should now be some decimal distance score

    return distance_score #After calculating, return distance score








#   Function Name: buildDistanceMatrix
#   Parameters: 
#               sequences - given a list of sequences, use that to generate the distance matrix by first getting a diagonal matrix of all the alignment scores to each other, then convert those scores to distance scores
#               subMatrix - given substitution matrix (in this case, our subMatrix is a dictionary) use this to compute the needleman-wunsch initial alignment scores to each other
#   Returns: distance_score_matrix - A matrix that has all of the alignment scores converted into distance scores with the normalizing and computing equations (slide deck 06 slide 20 and 21)
#   Purpose: Used to generate a distance matrix based off of all the combinations of alignments, and returns the distance matrix
def buildDistanceMatrix(sequences, subMatrix):
    alignment_score_matrix = [] #List representing the matrix from the alignment scores

    distance_score_matrix = [] #List representing the converted alignment score matrix to distance scores based on the method given in the lecture 06 slides - 20 and 21??

    distance_score_dict = {} #A dictionary associating whatever sequence combination is associated with a distance score


    #Establish an empty matrix of 0 alignment scores and distance scores, where the length of the matrix is the total number of sequences from sequences list/variable
    for i in range(len(sequences)): #Create a matrix with size N of number of sequences for row and column to represent the matricies
        alignment_score_matrix.append([]) #add empty list inside list to make a list of list (matrix) for alignment score matrix
        distance_score_matrix.append([]) #add empty lists to distance score matrix also
        
        for j in range(len(sequences)): 
            alignment_score_matrix[-1].append(0) #We use position [-1] because on every iteration through the first for loop, we are adding a new list to the end of the list of lists, and we want to specify to fill that list with zeroes
            distance_score_matrix[-1].append(0) #Same thing for the distance matrix
            
    #Test display to see the initial matricies
    #displayMatrix(alignment_score_matrix)  
    #print("\n\n")
    #displayMatrix(distance_score_matrix)  

    #Convert all key values into integer values so that we can find the min and max score value during normalization
    for k in subMatrix:
        subMatrix[k] = int(subMatrix[k])

    #Test print to see if the conversion worked
    #print(str(subMatrix))

    #Use a for loop to call needleman wunsch for every possible combination of sequences and perform
    #pairwise alignments on all of the combination of sequences and assignment to the alignment_score_matrix accordingly  
    for i in range(len(sequences)): #Iterate to the total length of the number of sequences from the sequences list 
        #Note: all sequences in the sequences list of lists are in position 1, because position 0 is the header for each sequence
        sequence_1  = sequences[i][1] #get sequence 1 to compare to every other sequence

        #Test display for the current sequence. Un-comment if you want to see the process of comparing each sequence to each other
        #print("\n\n------CURRENT SEQUENCE BEING COMPARED TO OTHER SEQUENCES "+str(sequences[i]))
        #print("current num of i == "+str(i))

        for j in range (len(sequences)): #Iterate to the total length of the number of sequences to compare each combination of sequences together
            #note: we iterate by i + 1(current sequence 1 + 1) because since its a diagonal matrix, we do not care about the other scores on the bottom half since they're the same score mirrored over the middle.
            sequence_2 = sequences[j][1] #get sequence 2 to compare to sequence 1, where sequence 1 has to be compared to ALL sequences in order to iterate to the next sequence combination
            
            #Run the needleman-wunsch algorithm and assign the returned alignment score tuple to needleman_list
            needleman_list = needlemanWunsch(sequence_1, sequence_2, subMatrix)

            #note: I used to have an if statement here that compared if the header of the sequences were the same (if the sequence was aligning to itself, then make it 0), but it was un-needed since I changed the range to iterate by i+1
            #Assign alignment score for this current alignment with sequence[i] and sequence[j] to the alignment score matrix at this current position
            alignment_score_matrix[i][j] = needleman_list[2] #Set the current position of the alignment score matrix to the needleman-wunsch score for the current sequence alignment 
            
            #Use the alignment_score_matrix and convert all those relevant scores to distance scores
            #Take the current position alignment score from the alignment_score_matrix and convert it to a distance score, then set the current position on the distance matrix to whatever that new distance score is
            distance_score_matrix[i][j] = convertAlignToDist(alignment_score_matrix[i][j], subMatrix, len(needleman_list[1])) #call convertAlignToDist function
            #note: sending the length of the current combination of alignment to convert the current cell into a distance score.
            #      needleman_list in format = (aligned seq 1, aligned seq 2, alignment_score)
            #      So pass the "length" of either aligned seq 1 or aligned seq 2. In this case, we are pasing however long aligned seq 2 is 
            #      also pasing the subMatrix (ex: nucleotide.mtx, BLOSUM50.mtx, etc.) because thats where we get the lowest and max possible score for normalizing 

            #This creates a dictionary in format:
            #   (seq 1, seq 2) = distance score
            #This is to make sure that we can look up the associations to keep track of what distance score is attached to what combination of sequences
            distance_score_dict[(sequences[i][0],sequences[j][0])] = distance_score_matrix[i][j]

            #Test display for the current sequence. Un-comment if you want to see the process of comparing each sequence to each other
            #print("\n\nSequence being compared to the current sequence: "+str(sequences[i]) +" TO "+ str(sequences[j]))
            #print(str(needleman_list))
            #print("current num of J == "+str(j))

            #print("THIS IS A TEST OF THE needleman_list: ")
            #print(str(needleman_list[1]))
            #print(str(len(needleman_list[1])))
            
    #Test display to see diagonal matrix after scoring.
    #print("\n\n\n")
    #displayMatrix(alignment_score_matrix) #display diagonal alignment score matrix before converting values to distances
    #print("\n\n")
    #displayMatrix(distance_score_matrix) #display diagonal distance score matrix after conversion
    #print("\n\n\n")
    #print(distance_score_dict)

    #Diagonal matrix should just have the top right corner scores, everything else INCLUDING THE MIDDLE SCORES are disregarded as 0.
    #use test display to check if the alignment score matrix was generated correctly.

    #After generating the alignment_score_matrix, we need to convert all these scores into distance scores by:
    #   1. Normalizing the score
    #   2. Converting normalized alignment score to a distance score
    #   3. place the distance score back into the new distance score matrix
    #   4. repeat for all relevant values (i+1)

    #There used to be another nested for loop, but I could just run the convertAlignToDist in the for loop above. (I changed it 10/21/23)
    #       distance_score_matrix[i][j] = convertAlignToDist(alignment_score_matrix[i][j], alignment_score_matrix, len()) #call convertAlignToDist function
            #THIS IS WHERE I STOPPED AT (10/19/23) ASK TOMORROW ABOUT HOW TO MEET WITH PROFESSOR DURING OFFICE HOURS AND FIND OUT HOW TO NORMALIZE THE SCORES

            #Note: I have met with Eric Rees in class and asked him a bunch of questions on (10/20/23).
            #Calculating distance:
            #   Normalizing Score:
            #       S(n) = (S(a)- S(min)) / (S(max) - S(min))
            #           Where: 
            #               S(min) = lowest possible score in the substitution matrix (from score matrix or whatever matrix was "given" to us (ex: nucleotide or BLOSUM)) *  length
            #                   length = the total length of the alignment sequence (after running NW) (ex: algined seq B = ---ACT length = 6)
            #                       note: the length from both sequences after alignment should be equal always (Ex: if I use NW on seq A and seq B, the resulting aligned sequences A and B should be the same length)
            #               S(max) = same calculation as S(min), but using the highest score in the substitution matrix from the "given" matrix
            #   Converting Normalized Score to distance:
            #       We just use natural log on the normalized score multiplied by -k where k is a scaling factor of our choosing.
            #       k's default value is 1 but we do this to try to get rid of as any decimal places as possible, so we can choose the scaling factor 
            #       to be anything, doesn't really matter.
            #       (k can be -10, -100, -1000, etc.) (Eric Rees recommended -10, but he said it doesn't matter)
            #       D = -k * ln(S(n))
            #           where D is the distance score
            #           where k is the scaling factor
            #   Then after all that, we should have a distance score.

    #Make sure that distance matrix has zeroes going through the middle
    distance_score_matrix = np.array(distance_score_matrix)

    np.fill_diagonal(distance_score_matrix, 0)
    
    return distance_score_matrix #After generating the distance matrix, return it to be used in constructing a guide tree with fitch-margoliash clustering algorithm
   
    #note: with distance_score_dict we can use this to lookup what are the distance scores for a given key value pair (where the key value pair is the desired combination of sequences you want to look up)



#   Function Name: constructGuideTree
#   Parameters: 
#               sequences - given a list of sequences 
#               subMatrix - given substitution matrix (in this case, our subMatrix is a dictionary)
#               D - The distance matrix from the buildDistanceMatrix function represented as a dictionary to look up values from based off of a given key
#   Returns: guide_tree - returns a guide tree based on the distance matrix
#   Purpose: Used to generate a guide tree with the fitch-margoliash algorithm from the given distance matrix
def constructGuideTree(sequences, subMatrix, distanceMatrix):
    #List of lists to represent the guide tree? idk what to use to represent the guide tree, maybe a stack **
    #STOPPED HERE (10/22/23) confused on the scope the algorithm since the matrix is changing size/dictionary is changing size
    guide_tree = []  #Represents a the guide_tree in a stack???
    
    temp_seqs = list(sequences) #This is used to iterate through the matrix and keep track of position values

    changed_num_distance_matrix = np.array(distanceMatrix)#converts the distance matrix into a numpy matrix to improve matrix manipulation
    
    original_reference_matrix = np.array(distanceMatrix) #Uses this untouched matrix to lookup values for averaging purposes

    #Met with Eric Rees today (10/23/23) and he said that he didn't use a dictionary, he just used a list of lists
    flag = True #flag is true by default for while loop

    loop_counter = 1 #This just keeps a number count of how many iterations it took to complete the clustering fitch_margoliash

    while flag:
        #********
        # How this works as far as I know:
        # 1. Find the smallest/shortest distance in the distance matrix.
        # 2. After finding the shortest distance, assign those two sequences found at those positions to be merged (positions of sequences is determined by temp_seqs)
        # 3. Delete the rows and columns that involve the two closest sequences (seq_1 and seq_2), this reduces the size of the overall matrix because we are merging two together
        # 4. Remove the closest sequences (seq_1 and seq_2) from the temp_seqs, and append them to the end of the temp_seqs list, changing the overall position of the
        #    sequences by 1. 
        #    Ex:
        #       before:
        #       length of temp_seqs = 5
        #       remember, the length of temp_seqs determines the size of the matrix. This is a starting matrix of 5 by 5 
        #       temp_seqs = [>A, >B, >C, >D, >E]  -> I identified >C and >D to be merged together since they are the closest distance
        #
        #       after:
        #       length of temp_seqs = 4
        #       remember, the length of temp_seqs determins the size of the next matrix. Since we merged them, the new matrix is going to be 4 by 4
        #       note: the new position of >E is 2 instead of 4 from the previous list above
        #       temp_seqs = [>A, >B, >E, (>C, >D)] -> Representing 
        #       
        # 5. add the two closest sequences being merged into the guide_tree
        # 6. Generate a new matrix to be used on the next loop, with the new matrix being the length of the temp_seqs after merging, so in the example above
        #    the original matrix (5 by 5 matrix) turns into a 4 by 4 matrix.
        # 7. Copy the values from the old matrix (after the rows and columns have been deleted from step 3) to the new matrix. There will be an extra new row and 
        #    new column that is empty at the end of the matrix to represent the newly merged values. The merged values need to be calculated with the fitch cluster method
        #    in the next step.
        # 8. To accurately fill the merged column values, I took the two sequences to be merged and compared their values to the rest of the temp_seqs by iterating through a for loop.
        #    Ex from above:
        #    temp_seqs = [>A, >B, >E, (>C, >D)] -> >C and >D are the closest, so average their values to get the merged column value (>C>D):
        #                   note: in this case/iteration, the number of sequences added is 2.
        #                   (>C to >A + >D to >A) / number of sequences added = (>C>D) to >A
        #                   (>C to >B + >D to >B) / number of sequences added = (>C>D) to >B
        #                   (>C to >E + >D to >E) / number of sequences added = (>C>D) to >E
        #   note: I also included multiple if statements to calculate if its merging a MSA to MSA, sequence to MSA, MSA to sequence, or both just sequences.
        #   Ex again:
        #   temp_seqs = [>A, >B, >E, (>C, >D)] -> in this case, lets say that the closest distance value is >E to (>C, >D):
        #   note: this means we have to merge sequence >E to (>C, >D)
        #         same idea, merge >E and change temp_seqs to change the length of the new matrix
        #
        #   the new temp_seqs becomes:
        #   temp_seqs = [>A, >B, (>E, >C, >D)] -> again, notice how the temp_seqs changes length, and indexes of the letters.
        #   note: then after merging, we just calculate the values AGAIN but change the number of sequences to be divided because we are averaging  
        #             note: in this case, we are averaging again but since we are merging 3 sequences together this time, we divide by 3
        #             (>E to >A + >C to >A + >D to A)  / number of sequences added = (>E>C>D) to >A    
        #             (>E to >B + >C to >B + >D to B)  / number of sequences added = (>E>C>D) to >B    
        #and so on, and so on. The algorithm keeps count of how many sequences we are adding together so that we can divide properly.
        # 9. (How do we get individual values if the new matrix keeps changing?) Thats what the purpose of these two variables are for:
        #           original_reference_matrix - This is an untouched version of the original distance matrix.
        #           sequences - This is an untouched version of the temp_seqs list. We use this in combination of the original_reference_matrix to find values.
        #           so, how it works:
        #           sequences = [>A, >B, >C, >D, >E] <- Not changing this list, using it for indexing into the original_reference_matrix
        #           positions = [0,  1,   2,  3,  4] <- not an actual variable in this program, just used to show where all the letter's indexes are
        #                                         A  B  C  D  E
        #           original_reference_matrix = A[0  1  2  3  4 ]      <- not changing this variable because we use this original matrix to look up specific                   
        #                                       B[1  0  7  8  9 ]         combination values.             
        #                                       C[2  7  0  13 14]         Ex: if we wanted to find one of the combinations above:
        #                                       D[3  8  13 0  19]                >E to >A -> we use positions: (4 to 0) = (4,0)  
        #                                       E[4  9  14 19 0 ]  or we can do: >A to >E -> we use positions: (0 to 4) = (0,4)
        #
        #
        # based on the temp_seqs list, since the index of the letters change, we can still find the values we are looking for just like how we used the original
        # matrix in the hand-written examples to find values and average them out when we were merging sequences together:
        # Ex: merge sequences -> temp_seqs = [>A, >B, >C, (>D, >E)]
        #       so we need to calculate:
        #                       (>D to >A + >E to >A) / 2 sequences added together = (>D>E) to >A
        #                       (>D to >B + >E to >B) / 2 sequences added together = (>D>E) to >B
        #                       (>D to >C + >E to >C) / 2 sequences added together = (>D>E) to >C
        #                       so we use sequences list to look up values in the original_reference_matrix:
        #                       find:
        #                        index     =   0   1   2     3
        #                        temp_seqs = [>A, >B, >C, (>D, >E)]
        #                       lets say, the next merge is >C to (>D to >E) where seq_1 = >C and seq_2 = (>D, >E)
        #                       
        #                       this is the new temp_seqs:
        #                         index     =  0    1       2
        #                         temp_seqs = [>A, >B, (>C, >D, >E)]
        #                       
        #                         index     =  0    1   2   3   4
        #                         sequences = [>A, >B, >C, >D, >E]

        #                           lets say, row is equal to 1.
        #
        #                        sequences.index(temp_seqs[row])
        #                         sequences.index(temp_seqs[1])
        #                           sequences.index(>B) = 1
        #                       
        #                        sequences.index(seq_1)
        #                           sequences.index(>C) = 2
        #
        #                           
        #                         original_reference_matrix[sequences.index(seq_1)][sequences.index(temp_seqs[row])]
        #                           original_reference_matrix[2][1]   = >C to >B value from original matrix 
        #     
        #                   and so on, and so on until we are done.                      
        #
        #                                                       
        #      for row in range(len(temp_seqs)) # 4 length of temp_seqs                
        #                                                                             >C    to        >A and then >B                                                  (>D, >E)     to      >A and then >B
        #            new_matrix[row][-1] = (original_reference_matrix[sequences.index(seq_1)][sequences.index(temp_seqs[row])] + original_reference_matrix[sequences.index(seq_2)][sequences.index(temp_seqs[row])] / number of total sequences added
        #
        #
        # 10. Set the new_matrix to the changed_num_matrix so that on the next iteration of the while loop, it will use the updated matrix
        #               
        # 11. repeat previous steps to complete the guide tree.
        # 12. After completing the guide tree, run MSA
        #********    
        #Test display for each iteration of the while loop
        print("\n\nTHIS IS ITERATION: "+str(loop_counter) +"---------------------------------------------------------")
        print("Starting List:")
        print(changed_num_distance_matrix)
        print("\n\n")
        print("Original Sequence list:")
        displayMatrix(sequences)
        print("\nTemp Sequences:")
        displayMatrix(temp_seqs)
        

        temp_list = [] #A temporary list to store the combined values after averaging, this will determine the size of the new matrix.


        default_min_num = changed_num_distance_matrix[0][1] #default min value to compare
    
        #finds the smallest number in the distance matrix
        for row in range(len(temp_seqs)):
            for col in range(row + 1, len(temp_seqs)):
                default_min_num = min(default_min_num, changed_num_distance_matrix[row][col]) 
        
        #Test display
        #print("SMALLEST NUM FOUND IN CURRENT MATRIX:")
        #print("number: "+str(default_min_num))

        #Finds the positon of the smallest number in the distance matrix and appends the sequences to the guide tree
        for row in range(len(temp_seqs)): 
            for col in range(row + 1, len(temp_seqs)): 
                if (changed_num_distance_matrix[row][col] == default_min_num):
                        seq_1 = temp_seqs[row] #Store the closest sequence names in seq_1 and seq_2
                        seq_2 = temp_seqs[col] 
                       # seq_1_index = row #To keep track of where the closest sequences are in the matrix
                       # seq_2_index = col #To keep track of where the closest sequences are in the matrix
                        break
            else:
                continue #If the inner loop doesn't break, keep going through the rows of the 2d matrix
            break #If the inner loop DOES break, then break out of the nested for loops
        
      
        
        #remove the rows and columns at the position of the lowest distance score to merge, merge their values.
        #Note: This deletion block gets rid of the rows and columns that we are merging, and will place the merged value at the end of the matrix.
        changed_num_distance_matrix = np.delete(changed_num_distance_matrix, row, 0) #Note: the third value, is 0 = row and 1 = column axis
        changed_num_distance_matrix = np.delete(changed_num_distance_matrix, row, 1)
        changed_num_distance_matrix = np.delete(changed_num_distance_matrix, col-1, 0) #-1 to make up for the change in matrix size from previous delete
        changed_num_distance_matrix = np.delete(changed_num_distance_matrix, col-1, 1)

        #Remove the closest found sequences from temp_seqs to be appended at the end as a new, 1 combined sequence

        temp_seqs.remove(seq_1) #remove the 2 sequences row and column from the sequence list because we are going to combine them, this will make all of the 
        temp_seqs.remove(seq_2) #other sequences the new indexes for the new matrix
        
        #Combine and merge the sequences. Just add on the string to an existing sequence/MSA
        temp_seqs.append(seq_1+seq_2) #Add combined sequence to the end of the temp_sequence list for indexing 

        #ATTEMPTING TO CONSTRUCT GUIDE TREE BASED ON WHETHER THE SEQUENCES BEING MERGED ARE MSA OR SEQ OR BOTH MSA OR BOTH SEQUENCES
        if (len(seq_1) <= 2 and len(seq_2) > 2): #If seq_2 is an MSA and seq_1 is just one sequence (sequence vs MSA)
                #print("ACTIVATED IF WHERE SEQ 2 IS AN MSA")
                guide_tree.append([sequences.index(seq_1), []])
                #print("SEQ 1 (non-MSA) = "+ str(seq_1))
                temp_list = []
                temp_index = []
                count = 0

                #print("SEQ 2 (MSA) = ")
                for i in seq_2: #Add ALL of the sequences in the MSA in seq_2
                    
                    temp_list.append(i) #Add first two elements
                    #print(temp_list)
                    count = count + 1
                    if count == 2:
                        #print("IF TEMP LIST ACTIVATED COUNT == 2: ")
                        guide_tree[-1][1].append(sequences.index(temp_list)) #append every sequence's index located in the MSA to the end of the empty list in the tree
                        print("THIS IS THE CURRENT GUIDE TREE")
                        print(guide_tree[-1])
                        print(temp_list)
                        count = 0
                        temp_list.clear()

                print("THIS IS THE GUIDE TREE PRINT ----1 " + str(guide_tree[-1]))
                print("SEQ_1 = "+str(seq_1))
                
                #Check to see if combined sequence has already been added to the list
                #if not, add the MSA to the list to index properly
                if guide_tree[-1][1] in sequences:
                    print("") #print empty string if this happens for flag
                    
                    #This is a very special case, where if the sequences being appended into the sequences list
                    #already exists in the sequences, and IF the first element is NOT found in the guide_tree[-1][1] list, then 
                    #create a new temp list, store the guide_tree[-1][1] elements into the list, then append that to the sequences list
                    check_if_list = False

                    for i in guide_tree[-1][1]:
                        if i == guide_tree[-1][0]:
                            check_if_list = True #element found

                    #print("MSA NOT FOUND IN SEQUENCES")
                    #if len(guide_tree[-1][1]) >= 2 and (guide_tree[-1]):
                      #  guide_tree[-1][1].insert(0, guide_tree[-1][0])
                      #  sequences.append(guide_tree[-1][1] )
                    if check_if_list == False:
                        temp_seq_index_list = []
                        temp_seq_index_list = list(guide_tree[-1][1])
                        temp_seq_index_list.insert(0, guide_tree[-1][0])
                        sequences.append(temp_seq_index_list)

                else:
                    #print("ELSE ACTIVATED")

                    
    
                    sequences.append(guide_tree[-1][1]) #Add the finished entry into the guide tree for indexing
                    #temp_index = list(guide_tree[-1][1]) #Store the combined list to be removed after finding index

                
            
                guide_tree[-1][1] = int(sequences.index(guide_tree[-1][1])) #Search for an existing combination
                


        elif(len(seq_2) <= 2 and len(seq_1) > 2): #if seq_1 is an MSA and seq_2 is just one sequence (MSA vs Sequence)
                #print("ACTIVATED IF WHERE SEQ 1 IS MSA")
                guide_tree.append([sequences.index(seq_2), []])
                #print("SEQ 2 (non-MSA) = "+ str(seq_2))
                temp_list = []
                count = 0

                #print("SEQ 1 (MSA) = ")
                for i in seq_1: #Add ALL of the sequences in the MSA in seq_1
                    
                    temp_list.append(i) #Add first two elements
                   #print(temp_list)
                    count = count + 1
                    if count == 2:
                        #print("IF TEMP LIST ACTIVATED COUNT == 2: ")
                        guide_tree[-1][1].append(sequences.index(temp_list)) #append every sequence's index located in the MSA to the end of the empty list in the tree
                        #print(temp_list)
                        count = 0
                        temp_list.clear()

                #Check to see if combined sequence has already been added to the list
                #if not, add the MSA to the list to index properly
               #Check to see if combined sequence has already been added to the list
                #if not, add the MSA to the list to index properly
                if guide_tree[-1][1] in sequences:
                    print("")
                    #print("MSA NOT FOUND IN SEQUENCES")
                else:
                    #print("ELSE ACTIVATED")
                    #sequences.append(list(guide_tree[-1][1])) #Add the finished entry into the guide tree for indexing
                    temp_index = list(guide_tree[-1][1]) #Store the combined list to be removed after finding index
                #After building the MSA with indexes, find the MSA in the sequences list

               
                guide_tree[-1][1] = int(sequences.index(guide_tree[-1][1])) #Search for an existing combination
          
            
        elif(len(seq_2) > 2 and len(seq_1) > 2): #MSA vs MSA
                #print("ACTIVATED IF SEQ 1 and SEQ 2 ARE BOTH MSA")
                guide_tree.append([[], []]) #This time, create 2 lists to represent both MSA's
                temp_list = []
                count = 0

                #print("SEQ 1 (MSA) = ")
                for i in seq_1: #Add ALL of the sequences in the MSA in seq_1
                    #print("SEQ 1 ADDING MSA")
                    temp_list.append(i) #Add first two elements
                    #print(temp_list)
                    count = count + 1
                    if count == 2:
                        #print("IF TEMP LIST ACTIVATED COUNT == 2: ")
                        guide_tree[-1][0].append(sequences.index(temp_list)) #Append every sequence to the first MSA empty list
                        #print(temp_list)
                        count = 0
                        temp_list.clear()

                #Check to see if combined sequence has already been added to the list
                #if not, add the MSA to the list to index properly
                if guide_tree[-1][0] in sequences:
                    print("")
                    #print("MSA NOT FOUND IN SEQUENCES")
                else:
                    #print("ELSE ACTIVATED")
                    #sequences.append(list(guide_tree[-1][1])) #Add the finished entry into the guide tree for indexing
                    temp_index = list(guide_tree[-1][0]) #Store the combined list to be removed after finding index
                #After building the MSA with indexes, find the MSA in the sequences list

                
                guide_tree[-1][0] = int(sequences.index(guide_tree[-1][0])) #Search for an existing combination
                #Then we repeat the process above and do it for the 2nd MSA

                #print("SEQ 2 (MSA) = ")
                for i in seq_2: #Add ALL of the sequences in the MSA in seq_2
                    #print("SEQ 2 ADDING MSA")
                    temp_list.append(i) #Add first two elements
                    #print(temp_list)
                    count = count + 1
                    if count == 2:
                        #print("IF TEMP LIST ACTIVATED COUNT == 2: ")
                        guide_tree[-1][1].append(sequences.index(temp_list)) #Append every sequence to the second MSA empty list
                        #print(temp_list)
                        count = 0
                        temp_list.clear()
                
                #Check to see if combined sequence has already been added to the list
                #if not, add the MSA to the list to index properly
                if guide_tree[-1][1] in sequences:
                    print("")
                    #print("MSA NOT FOUND IN SEQUENCES")
                else:
                    #print("ELSE ACTIVATED")
                    #sequences.append(list(guide_tree[-1][1])) #Add the finished entry into the guide tree for indexing
                    temp_index = list(guide_tree[-1][1]) #Store the combined list to be removed after finding index
                #After building the MSA with indexes, find the MSA in the sequences list
                
                
                guide_tree[-1][1] = int(sequences.index(guide_tree[-1][1])) #Search for an existing combination

                #After finding the correct guide tree, replace the MSA with the previous guide_tree iteration

                #print("GUIDE TREE MSA TO MSA TEST PRINT")
                #print(guide_tree)
                #displayMatrix(sequences)


        else: #Sequence Vs Sequence
                #print("ACTIVATED ELSE WHERE NONE OF SEQ 1 OR SEQ 2 ARE MSAs")
                #print("SEQ 1 = "+str(seq_1))
                #print("SEQ 2 = "+str(seq_2))
                guide_tree.append([sequences.index(seq_1) , sequences.index(seq_2)]) #Append the two sequences to the tree as one unit. 
                sequences.append([sequences.index(seq_1) , sequences.index(seq_2)]) #Append the two sequences to the sequences list to keep track of index
                
                #print("GUIDE TREE SO FAR: ")
                #print(guide_tree)
                #print("SEQUENCES LIST SO FAR:")
                #print(sequences)

        #Because of the initial merging, the code above appends 2 of the first merging sequence into the sequences list, but thats okay since the tree program still finds
        #the correct indexes for all of the values (10/26/23)
        
        

        #With this new distance matrix, generate a new matrix while calculating the values out for combined sequences:
        new_matrix = []#represents a new matrix from the previous matrix, plus the newly combined sequences
        #Generate new matrix with the dimensions of the new temp_sequences length
        #Note: THIS FOR LOOP DETERMINES THE SIZE OF THE NEXT MATRIX BY INITIALIZING A NEW MATRIX BASED ON THE NEW LENGTH OF temp_seqs
        for row in range(len(temp_seqs)):
            new_matrix.append([])
            for col in range(len(temp_seqs)):
                new_matrix[-1].append(0)

        #copy previous values from changed_num_distance_matrix to new generated matrix
        for row in range(len(changed_num_distance_matrix)):
            for col in range(len(changed_num_distance_matrix)):
                new_matrix[row][col] = changed_num_distance_matrix[row][col]

        #print("\n\n")
        #print("New temp sequences")
        #print(temp_seqs)
        #print("\n\n")

        #This for loop uses the clustering. For each iteration, it fills in the new merged values of the matrix based on what sequences are being merged.
        #THIS IS THE FITCH-MARGOLIASH FOR LOOP CLUSTERING 
        for row in range(len(temp_seqs)-1):
            #print("\n\nSTART OF NEW Sequence COMBO CALC: COMBINE SEQUENCES: "+str(seq_1)+" AND "+str(seq_2))
            #print("THIS IS VALUE 1 (Seq 1): ")
            #print("Value 1 = "+str(original_reference_matrix[sequences.index(seq_1)][sequences.index(temp_seqs[row])]))
            #print("\n\n\n")
            #print(seq_1)
            #print(len(seq_1))
            #print(seq_2)
            #print(len(seq_2))

            #print("THIS is VALUE 2 (Seq 2): ")
            #print("Value 2 = "+str(original_reference_matrix[sequences.index(seq_2)][sequences.index(temp_seqs[row])]))
            if (len(seq_1) <= 2 and len(seq_2) > 2): #If seq_2 is an MSA and seq_1 is just one sequence (sequence vs MSA)
                print("ACTIVATED IF WHERE SEQ 2 IS AN MSA")

                if len(temp_seqs[row]) > 2: #if we are adding up merged row:
                    #print("ADDING MERGED ROW ACTIVATED")
                    #print("SEQ 1: "+str(seq_1))
                    #print("SEQ 2: "+str(seq_2))
                    temp_list = [] #to temporarily store the read sequences from the MSA
                    temp_list_seq_2 = [] #to temporarily store the read sequences from MSA seq_2
                    count2 = 0 # To keep track of the sequences added to the temp_list_seq_2
                    count = 0 #To keep track of the sequences added to the temp_list
                    operation_counter = 0 #to divide by how many operations we had to do 

                    for i in temp_seqs[row]: #Add ALL of the sequences from seq_1 (sequence) to MSA (temp_seqs[row])
                    
                        temp_list.append(i) #Add first two elements
                        count = count + 1
                        if count == 2:
                            #print("IF TEMP LIST ACTIVATED COUNT == 2 FOR SEQ_1: ")
                            #print(temp_list)
                            #print("compare to Seq_1")
                            #print(seq_1)
                            operation_counter = operation_counter + 1 #Add 1 to the operation counter
                            new_matrix[row][-1] = new_matrix[row][-1] + original_reference_matrix[sequences.index(seq_1)][sequences.index(temp_list)]
                            count = 0 #reset the count to store new operation
                            temp_list.clear() #clear the temp_list to store another sequence from the MSA
                    

                    #Do it again, but this time, seq_2 is an MSA to an MSA (temp_seqs[row])
                    for i in temp_seqs[row]: #Add ALL of the sequences from seq_2 (MSA) to MSA (temp_seqs[row])
                        temp_list.append(i) #Add first two elements
                        count = count + 1
                        if count == 2:
                            #print("IF TEMP LIST ACTIVATED COUNT == 2: ")
                            #print(temp_list)
                            #print("CURRENT ROW:")
                            #print(temp_seqs[row])
                            
                            for j in seq_2:  #Since seq_2 is also an MSA that we have to add, make sure to grab from seq_2 as well
                                temp_list_seq_2.append(j)
                                count2 = count2 + 1
                                if count2 == 2:
                                    #print("IF TEMP LIST ACTIVATED COUNT INNER**** == 2: ")
                                    #print(temp_list_seq_2)
                                    #print("CURRENT seq_2:")
                                    #print(seq_2)
                                    operation_counter = operation_counter + 1 #Add 1 to the operation counter
                                    new_matrix[row][-1] = new_matrix[row][-1] + original_reference_matrix[sequences.index(temp_list_seq_2)][sequences.index(temp_list)]
                                    count2 = 0 #clear seq_2 count2
                                    temp_list_seq_2.clear() #clear seq_2 temp_list_seq_2
                                    
                            count = 0 #clear temp_seqs[row] count
                            temp_list.clear() #clear temp_seqs[row] temp_list

                    #print("TOTAL OPERATIONS = " +str(operation_counter))

                    #After calculating what seq_1 (a sequence) to MSA (temp_seqs[row]) and seq_2 (an MSA) to MSA (temp_seqs[row]), divide the cell by the total number of operations
                    new_matrix[row][-1] = new_matrix[row][-1] / operation_counter #divide the current cell with the total number of operations

                    new_matrix[-1][row] = new_matrix[row][-1] #Mirror distance values
                

                else:
                    temp_list = []
                    count = 0 #To keep track of the sequences added to the temp_list
                    total_seqs = 0 #Total sequences counter
                    operation_counter = 0 #to divide by how many operations we had to do 
                    for i in seq_2: #Add ALL of the sequences in the MSA in seq_2
                    
                        temp_list.append(i) #Add first two elements
                        #print(temp_list)
                        count = count + 1
                        if count == 2:
                            total_seqs = total_seqs + 1 #add 1 to the total_seqs since 1 sequence has length 2
                            #print("IF TEMP LIST ACTIVATED COUNT == 2: ")
                            #print(temp_list)
                            #print("CURRENT ROW:")
                            #print(temp_seqs[row])
                            #print("VALUE FOUND " + str(original_reference_matrix[sequences.index(temp_list)][sequences.index(temp_seqs[row])]))
                            new_matrix[row][-1] = new_matrix[row][-1] + original_reference_matrix[sequences.index(temp_list)][sequences.index(temp_seqs[row])]
                            count = 0
                            temp_list.clear()
                    
                    new_matrix[row][-1] = new_matrix[row][-1] + original_reference_matrix[sequences.index(seq_1)][sequences.index(temp_seqs[row])] #Add the last sequence 1 value to combine them
                
                    #print("VALUE FOUND "+ str(original_reference_matrix[sequences.index(seq_1)][sequences.index(temp_seqs[row])]))

                    total_seqs = total_seqs + 1 #because we just added the other sequence to the MSA

                    #print("TOTAL_SEQS "+str(total_seqs))
                

                    new_matrix[row][-1] = new_matrix[row][-1] / total_seqs

                    new_matrix[-1][row] = new_matrix[row][-1] #Mirror distance values

                    #print(new_matrix[row][-1])

            elif(len(seq_2) <= 2 and len(seq_1) > 2): #if seq_1 is an MSA and seq_2 is just one sequence (MSA vs Sequence)
                #print("ACTIVATED IF WHERE SEQ 1 IS MSA")
                temp_list = []
                count = 0
                total_seqs = 0
                for i in seq_1: #Add ALL of the sequences in the MSA in seq_1
                    
                    temp_list.append(i) #Add first two elements
                    #print(temp_list)
                    count = count + 1
                    if count == 2:
                        total_seqs = total_seqs + 1 #add 1 to the total_seqs since 1 sequence has length 2
                        #print("IF TEMP LIST ACTIVATED COUNT == 2: ")
                        #print(temp_list)
                        #print("CURRENT ROW:")
                        #print(temp_seqs[row])
                        #print("VALUE FOUND " + str(original_reference_matrix[sequences.index(temp_list)][sequences.index(temp_seqs[row])]))
                        new_matrix[row][-1] = new_matrix[row][-1] + original_reference_matrix[sequences.index(temp_list)][sequences.index(temp_seqs[row])]
                        count = 0
                        temp_list.clear()
                    
                new_matrix[row][-1] = new_matrix[row][-1] + original_reference_matrix[sequences.index(seq_2)][sequences.index(temp_seqs[row])] #Add the last sequence 1 value to combine them
                
                #print("VALUE FOUND "+ str(original_reference_matrix[sequences.index(seq_2)][sequences.index(temp_seqs[row])]))

                total_seqs = total_seqs + 1 #because we just added the other sequence to the MSA

                #print("TOTAL_SEQS "+str(total_seqs))

                new_matrix[row][-1] = new_matrix[row][-1] / total_seqs

                new_matrix[-1][row] = new_matrix[row][-1] #Mirror distance values

                #print(new_matrix[row][-1])
            
            elif(len(seq_2) > 2 and len(seq_1) > 2): #MSA vs MSA
                print("ACTIVATED IF SEQ 1 and SEQ 2 ARE BOTH MSA")
                temp_list = []
                count = 0
                total_seqs = 0
                for i in seq_1: #Add ALL of the sequences in the MSA in seq_1
                    #print("SEQ 1 ADDING MSA")
                    temp_list.append(i) #Add first two elements
                    #print(temp_list)
                    count = count + 1
                    if count == 2:
                        total_seqs = total_seqs + 1 #add 1 to the total_seqs since 1 sequence has length 2
                        #print("IF TEMP LIST ACTIVATED COUNT == 2: ")
                        #print(temp_list)
                        #print("CURRENT ROW:")
                        #print(temp_seqs[row])
                        #print("VALUE FOUND " + str(original_reference_matrix[sequences.index(temp_list)][sequences.index(temp_seqs[row])]))
                        new_matrix[row][-1] = new_matrix[row][-1] + original_reference_matrix[sequences.index(temp_list)][sequences.index(temp_seqs[row])]
                        count = 0
                        temp_list.clear()
                
                for i in seq_2: #Add ALL of the sequences in the MSA in seq_2
                    #print("SEQ 2 ADDING MSA")
                    temp_list.append(i) #Add first two elements
                    #print(temp_list)
                    count = count + 1
                    if count == 2:
                        total_seqs = total_seqs + 1 #add 1 to the total_seqs since 1 sequence has length 2
                        #print("IF TEMP LIST ACTIVATED COUNT == 2: ")
                        #print(temp_list)
                        #print("CURRENT ROW:")
                        #print(temp_seqs[row])
                        #print("VALUE FOUND " + str(original_reference_matrix[sequences.index(temp_list)][sequences.index(temp_seqs[row])]))
                        new_matrix[row][-1] = new_matrix[row][-1] + original_reference_matrix[sequences.index(temp_list)][sequences.index(temp_seqs[row])]
                        count = 0
                        temp_list.clear()
                
                #Test for total sequences
                #print("TOTAL SEQUENCES: "+str(total_seqs))
                
                new_matrix[row][-1] = new_matrix[row][-1] / total_seqs

                new_matrix[-1][row] = new_matrix[row][-1] #Mirror distance values


                


            else: #Sequence Vs Sequence
                print("SEQUENCE VS SEQUENCE ACTIVATED")
                #print(str(seq_1) + " TO " +str(temp_seqs[row]))
                #print("PLUS")
                #print(str(seq_2) + " TO " +str(temp_seqs[row]))
                #print("ACTIVATED ELSE WHERE NONE OF SEQ 1 OR SEQ 2 ARE MSAs")
                #print("\n\nPART OF THE TEMP_SEQS[ROW]")
                #print(temp_seqs[row])

                if len(temp_seqs[row]) > 2: #If the row we are currently at a merged row
                    #print("----------ACTIVATED IF STATEMENT---------")
                    #print("TEMP SEQ LIST: ")
                    #displayMatrix(temp_seqs)
                    #print("TEMP SEQ[ROW]: ")
                    #print(temp_seqs[row])

                    counter = 0 #to keep track of how many sequences there are
                    temp_seq = [] #to keep track of what sequence we are finding the index for 
                    operation_counter = 0 #to divide by how many operations we had to do 

                    #Add all combinations in seq_1 to everything in current msa
                    for i in temp_seqs[row]:
                        counter = counter + 1
                        temp_seq.append(i) #Add the element into the temp_seq list
                        if counter == 2:
                            operation_counter = operation_counter + 1
                            new_matrix[row][-1] = new_matrix[row][-1] + original_reference_matrix[sequences.index(seq_1)][sequences.index(temp_seq)]
                            counter = 0
                            temp_seq.clear()

                    #Do this operation again, and add every combination of seq_2 with the current msa
                    for i in temp_seqs[row]:
                        counter = counter + 1
                        temp_seq.append(i) #Add the element into the temp_seq list
                        if counter == 2:
                            operation_counter = operation_counter + 1
                            new_matrix[row][-1] = new_matrix[row][-1] + original_reference_matrix[sequences.index(seq_2)][sequences.index(temp_seq)]
                            counter = 0
                            temp_seq.clear()
                    
                    #print("SEQUENCES TO BE MERGED")
                    #print(seq_1)
                    #print(seq_2)
                    #print(temp_seqs[row])
                    #print(operation_counter)

                    #Get rid of this line below
                    new_matrix[row][-1] = new_matrix[row][-1] / operation_counter #Calculates the last columns of the diagonal matrix (The combined sequence column)
                    new_matrix[-1][row] = new_matrix[row][-1] #Mirror values
                else:
                    new_matrix[row][-1] = (original_reference_matrix[sequences.index(seq_1)][sequences.index(temp_seqs[row])] + original_reference_matrix[sequences.index(seq_2)][sequences.index(temp_seqs[row])]) / 2 #Calculates the last columns of the diagonal matrix (The combined sequence column)
                    new_matrix[-1][row] = new_matrix[row][-1] #Mirror values


        new_matrix = np.array(new_matrix)# turn new_matrix into a np.array
        
        changed_num_distance_matrix = new_matrix # use the new_matrix as the next changed_num_distance
       
        print("\n\nNEW MATRIX")
        print(new_matrix)
        print("\n\nTEMP SEQS")
        displayMatrix(temp_seqs)
        print("\n\nSEQUENCES LIST")
        displayMatrix(sequences)

        #print("Seq 1 closest")
        #print(seq_1)
        #print("Seq 2 closest")
        #print(seq_2)
        #print("\n\n")
        print("Guide tree:")
        print(guide_tree)

        #Combine the two closest sequences and add them to the end of the row and col, and the temp seqs for indexing
        #Add another row and column for the new combined entry to the changed_num_distance_matrix
        

        loop_counter = loop_counter + 1 #Add one to the loop counter to keep track of the iterations for test display
        
    
        #If the length of the sequences is less than 2 sequences, stop the while loop.
        if len(temp_seqs) < 2:
            flag = False

    
    
    #Test to see if the index of sequences is correct
    #print("\n\nFINAL SEQUENCES:")
    #print(sequences)

    #The purpose of this next while loop is to replace all the merged sequences that we used to create guide tree with the correct tree combination for indexing through the
    #sequences list properly when creating the MSA

    #guide_tree_index = 0 #Start of the guide tree

    #starting_index = sequences.index(guide_tree[0]) #Start of the guide tree off of sequences list

    #flag2 = True #2nd while loop flag

    #loop_counter2 = 1 #2nd while loop iteration counter for testing purposes
    
    #while flag2: #REMOVE THIS WHILE LOOP AFTER
        
        #print("while loop activated: "+ str(loop_counter2))
    
    #    sequences[starting_index] = guide_tree[guide_tree_index]

    #    guide_tree_index = guide_tree_index + 1

    #    starting_index = starting_index + 1

    #    loop_counter2 = loop_counter2 + 1

        #print(sequences)

    #    if guide_tree_index >= len(guide_tree) - 1: #Check if the index has reached the guide tree length - 1 limit
    #        flag2 = False

    #After generating the guide tree, return it as T= in the main function
    return guide_tree



#   Function Name: getMSA
#   Parameters: 
#               sequences - given a list of sequences 
#               subMatrix - given substitution matrix (in this case, our subMatrix is a dictionary)
#               T - The generated guide tree to show us the order in which to perform pairwise alignment
#   Returns: MSA - returns MSA after performing it from the guide tree
#   Purpose: Used to perform pairwise alignments based on the order of the given tree and return it as a MSA
def getMSA(indexSequencesWithT, referenceSequences, subMatrix, T):
    # I am going to try to read through the guide tree in order by using the indexes from the tree and the indexes from the sequences list to find the
    # correct order to perform needleman wunsch.
    MSA = [] 

    #Test display before converting all the sequences found in the indexSequencesWithT list into index values
    #print("SEQUENCES LIST WITH INDEXES FOR EVERYTHING")
    #print(indexSequencesWithT)
    #print("SEQUENCES LIST CONTAINING ORIGINAL SEQUENCES")
    #print(referenceSequences)
    

    #First, convert the sequences in indexSequencesWithT to just their indexes from referenceSequences
    end_index = indexSequencesWithT.index(T[0]) #Defines where the end of the iteration will be, which is the start of wherever the guide tree starts in this list

    for i in range(end_index): #From 0 (start) to end index (wherever the guide tree starts in the indexSequencesWithT list)
        indexSequencesWithT[i] = referenceSequences.index(indexSequencesWithT[i]) #Replace the actual sequence in the indexSequencesWithT list with their original index

    #Test display to check if the sequences were converted into index
    #print("\n\n")
    #print(indexSequencesWithT)
    #print(referenceSequences)
    #print("\n\nGUIDE TREE TO ITERATE THROUGH: ")
    #print(T)

    #Now we have two lists:
    # one list that contains the actual sequences: referenceSequences
    # one list that contains the indexes for the entire guide tree to read through: indexSequencesWithT

    #I guess the idea is:
    #I want to go through the entire guide tree by looking for it's indexes for everything that I need, and perform NW on them with the MSA vs MSA or MSA vs Seq rules etc.
    #Iterate through the guide tree and perform NW based on index values from indexSequencesWithT
    #Note: all guide tree entries will be in this format: [seq1 or MSA1 , seq2 or MSA2] 
    #                                                          i[0]     ,     i[1]
    for i in T:

        #Find out what operation we are doing based off of the guide tree's index value pairs
        #Test display to show what current pair operation we are on:
        print("\nCURRENT OPERATION")
        print(i)
        print("\n")
        


        #IF the 1st number in the guide tree is a sequence (an int) and the 2nd number is a sequence (an int): sequence vs sequence
        if isinstance(indexSequencesWithT[i[0]], int) and isinstance(indexSequencesWithT[i[1]], int):
            print("\nguide tree says: ITS A SEQUENCE VS A SEQUENCE") 
            #If it is a sequence vs a sequence, just align the first sequence with the second sequence via needleman-wunsch (from slides)
            #Run standard needleman wunsch if its a sequence with a sequence
            #Test display to check if the two sequences were pulled correctly 
            #print("\nSeq from 1st num")
            #print(referenceSequences[indexSequencesWithT[i[0]]][1]) #Pull sequence 1 from referenceSequences
            #print("\nSeq from 2nd num")
            #print(referenceSequences[indexSequencesWithT[i[1]]][1]) #Pull sequence 2 from referenceSequences

            #Run needleman-wunsch on the two sequences found from the indexes from the guide tree and assign tuple to needleman_list
            needleman_list = needlemanWunsch(referenceSequences[indexSequencesWithT[i[0]]][1], referenceSequences[indexSequencesWithT[i[1]]][1], subMatrix)

            #Test display for needleman_list
            #print(needleman_list)

            #note: the tuple is returned as: (seq_1 aligned sequence, seq_2 aligned sequence, optimal alignment score)

            #REMEMBER TO REPLACE ALL GAPS (-) WITH X'S. ALWAYS REPLACE GAPS WITH X AFTER EACH ALIGNMENT

            #go through both sequences and replace each gap with X. Because I stored the output of needleman wunsch in a tuple, we need to assign 
            #the newly edited sequences to another variable, then replace the existing sequences with that new variable.
            needleman_list[0] = needleman_list[0].replace("-", "X") #Replace gaps with X for aligned sequence 1

            needleman_list[1] = needleman_list[1].replace("-", "X") #Replace gaps with X for aligned sequence 2
            
            #print(needleman_list)
            #Test to see if the gaps were replaced with X's successfully
            #print(new_aligned_seq1)
            #print(new_aligned_seq2)

            #Reassign the newly aligned sequences from needleman tuple back into the referenceSequences
            referenceSequences[indexSequencesWithT[i[0]]][1] = needleman_list[0] #reassign newly aligned seq_1 back to reference sequence list to be used with other guide tree operations
            referenceSequences[indexSequencesWithT[i[1]]][1] = needleman_list[1] #same thing to aligned seq_2, reassign back to reference sequence list

            #Test display to check if the newly aligned sequences replaced the old sequences in the reference sequences list
            displayMatrix(referenceSequences)

            #And thats all you do for this case in the slides (no need to insert X's or anything because we are not copying gaps in sequence vs sequence)

        #IF the 1st number in the guide tree is a sequence (an int) and the 2nd number is an MSA (a list): sequence vs. MSA
        elif isinstance(indexSequencesWithT[i[0]], int) and len(indexSequencesWithT[i[1]]) >= 2:
            print("\nguide tree says: ITS A SEQUENCE VS AN MSA")
            #BEFORE WE EVEN START TO ALIGN, WE HAVE TO TURN ALL OF THE GAPS IN EACH SEQUENCE INVOLVED INTO AN X
            #If it is a sequence vs an MSA, we have to align whatever the 1st number's sequence is to EVERY sequence found in the 2nd guide tree number list using needleman-wunsch.
            #From all of the alignment scores, pick the BEST alignment score out of all of the combinations and copy all the gaps from 
            #So, because in this case, the 1st number in our guide tree represents a single sequence against a MSA (a list), we want to find all the possible aligments combinations
            #       from the 1st number sequence to every sequence in the 2nd number guide tree MSA, and pick the alignment from 1st number sequence to whatever MSA had the highest
            #       alignment score, and copy all the gaps from that sequence over to all of the other sequences in the MSA
            
            #Test display to check if the two sequences were pulled correctly 
            #print("\nSeq from 1st num")
            #print(referenceSequences[indexSequencesWithT[i[0]]][1]) #Pull sequence 1 from referenceSequences
            seq_1 = referenceSequences[indexSequencesWithT[i[0]]][1] #Sequence 1 to go into the needleman-wunsch
            #print("\nMSA from 2nd num")
            #print(indexSequencesWithT[i[1]]) #pull the MSA (2nd number in the guide tree has index for MSA in indexSequencesWithT) 
            seq_2_MSA = indexSequencesWithT[i[1]] #Assign seq_2 MSA for organization
            #print(seq_2_MSA) #Test print
            #Temp list of lists to hold all of the sequence alignments
            temp_aligned_list = [] # variable to hold all the combined alginment tuples
            
            #Align X with EACH Y (so every sequence in the MSA) and find the best score.
            #Store all the alignments into temp_aligned_list
            for j in seq_2_MSA:
                needleman_list = needlemanWunsch(seq_1, referenceSequences[j][1], subMatrix) #Run needleman-wunsch over EVERY combination of X with all Y's 
                temp_aligned_list.append(needleman_list)

            #Test display to see if the alignments were output correctly
            print("TEMP ALIGNED LIST")
            displayMatrix(temp_aligned_list)
            print("Seq_2_MSA: ")
            print(seq_2_MSA)
            print("\n\n")
            
            #Find the alignment that has the best optimal alignment score:
            #Assign an initial current max number to compare to:
            curr_max = temp_aligned_list[0][2] #Setting it to the first alignment's optimal score as a initial max value to compare
            #compare the curr_max of the temp_aligned list (lists of alignments with all possible Y's to X) and find the highest alignment score
            for b in temp_aligned_list:
                if b[2] > curr_max:
                    curr_max = b[2]

            #Test display to see if the curr_max was obtained
            #print(curr_max)

            
            best_score_pos = 0 #A variable that holds the index of the best score for the temp_aligned list

            #Get the best score's index from temp_aligned_list
            for c in range(len(temp_aligned_list)):
                if temp_aligned_list[c][2] == curr_max:
                    best_score_pos = c
                    break #once it finds the index of the best score, stop the loop
                    
            #Once we obtained the best score out of all the alignments, copy all of the gaps from the highest score from temp_aligned_list[highest score position][2nd sequence]
            #Since we want to inject all the gaps within the MSA. Do not mess with the first seq_1 because its just a sequence, not an MSA like seq_2_MSA
            
            #Test display
            #print(best_score_pos)
            

            #Replace the value of the seq_1 that had the best_score_pos
            referenceSequences[indexSequencesWithT[i[0]]][1] = temp_aligned_list[best_score_pos][0]

            referenceSequences[indexSequencesWithT[i[0]]][1] = referenceSequences[indexSequencesWithT[i[0]]][1].replace("-", "X")


            #For every other sequence that ISNT the best_score_pos, insert gaps from best_score_pos
            ## FIX THIS PART, WE NEED TO ADD THE GAP INSERTION FUNCTIONALITY FOR THIS TO WORK ## (10/29/23)
            for g in range(len(seq_2_MSA)):
                
                #print("CURRENT POSITION G -------------------")
                #print(g)
                #print(seq_2_MSA[g])
                #print("Current sequence that I am copying best sequence gaps over to: ")
                #print(referenceSequences[seq_2_MSA[g]])
                #print("Sequence to copy gaps: ")
                #print(temp_aligned_list[best_score_pos][1])


                #for each other sequence, copy all the gaps from highest scoring sequence to their other sequences
                letter_iterator = -1
                temp_build_new_seq = ""
                referenceSequences_length = len(referenceSequences[seq_2_MSA[g]][1]) #length of the current copy we are doing 
                
                
                #note: temp_aligned_list[best_score_pos][1] = best score sequence 1 to copy to everyone else
                #note: referenceSequences[seq_2_MSA[g]][1] = the current sequence to copy gaps to
                #THIS INNER FOR LOOP RE-BUILDS THE COPYING SEQUENCE BACKWARDS
                for letter_pos in range(len(temp_aligned_list[best_score_pos][1]) - 1, 0 - 1, -1): #Added -1 to add an offset to iterate through the letters backwards
                    if temp_aligned_list[best_score_pos][1][letter_pos] == "-":
                        #print("ADD A GAP")
                        temp_build_new_seq = temp_build_new_seq + "-"
                        continue

                    temp_build_new_seq = temp_build_new_seq + referenceSequences[seq_2_MSA[g]][1][letter_iterator] #Build the other sequence, but add a gap whenever it comes to it

                    letter_iterator = letter_iterator - 1

                    #print(str(letter_pos) +" --- " + str(temp_aligned_list[best_score_pos][1][letter_pos] + " --- " + str(temp_build_new_seq)))

                    if letter_iterator >= referenceSequences_length:
                        break
                
                #print("FINISHED SEQUENCE BUILT: ")
                #print(temp_build_new_seq[::-1])

                #Assign the new copied gap sequence into the main referenceSequences list
                referenceSequences[seq_2_MSA[g]][1] = temp_build_new_seq[::-1]

                #optional, ###I do not know if I need to make the copied gaps into X's --------------------------------------
                #Confirmed that this is what you do, you're supposed to replace the new gapped sequences with X's (10/31/23)
                referenceSequences[seq_2_MSA[g]][1] = referenceSequences[seq_2_MSA[g]][1].replace("-", "X")
                #The real question is, after I copy the gaps, do I make them into X's or not? idk

            #After copying over gaps to each other sequence in the MSA, replace gaps with X in all alignments for MSA
            for e in temp_aligned_list:
                for f in range(len(e) - 1):
                    e[f] = e[f].replace("-", "X") 

            #print(temp_aligned_list)

            #replace the value of the seq_2_MSA that had the best_score_pos
            #print(best_score_pos)
            #print(indexSequencesWithT[i[1]][best_score_pos])
            #print(referenceSequences[indexSequencesWithT[i[1]][best_score_pos]])
            #print(referenceSequences[indexSequencesWithT[i[1]][best_score_pos]][1])
            referenceSequences[indexSequencesWithT[i[1]][best_score_pos]][1] = temp_aligned_list[best_score_pos][1]

            #Test display to check if the newly aligned sequences replaced the old sequences in the reference sequences list
            displayMatrix(referenceSequences)

            #Read and store all the sequences from the MSA into the temp_seq_list.
            #This is complicated. Will stop here for the time being (10/27/23) may ask him how to track our way throughout this MSA to get all of the sequences


        #IF the 1st number in the guide tree is an MSA (a list) and the 2nd number is a sequence (an int): MSA vs sequence
        #elif len(indexSequencesWithT[i[0]]) >= 2 and isinstance(indexSequencesWithT[i[1]], int):
        #    print("guide tree says: ITS AN MSA VS A SEQUENCE")



        #IF the 1st number in the guide tree is an MSA (a list) and the 2nd number is an MSA (a list): MSA vs MSA
        elif len(indexSequencesWithT[i[0]]) >= 2 and len(indexSequencesWithT[i[1]]) >= 2: 
            print("\nguide tree says: ITS AN MSA VS AN MSA")

            #Get the MSA's from the guide tree
            seq_1_MSA = indexSequencesWithT[i[0]] #Assign seq_1 MSA for organization

            seq_2_MSA = indexSequencesWithT[i[1]] #Assign seq_2 MSA for organization

            #In MSA to MSA, every sequence in seq_1_MSA is aligned with every sequence in sequence_2_MSA
            #Ex:      MSA_1 = (AB)    MSA_2 = (CDE)
            #        NW ->         A to C   some_score...
            #        NW ->         A to D   some_score...
            #        NW ->         A to E   some_score...
            #        NW ->         B to C   some_score...
            #        NW ->         B to D   some_score...
            #        NW ->         B to E   some_score...
            #Same idea, as the seq vs MSA, but when we choose the highest scoring optimal alignment, we have to now copy gaps to both MSA.
            #   so lets say in the above example, A to E was the highest/best. Then in MSA_1, we would have to copy gaps from A to B,
            #   and in MSA_2 you would have to copy gaps from E to C and D.
            
            #Temp list of lists to hold all of the sequences alignments
            temp_aligned_list = []

            #perform needleman wunsch for all combinations from seq_1_MSA to seq_2_MSA
            for i in seq_1_MSA: #Pull sequence index from seq_1_MSA
                for j in seq_2_MSA: #Pull sequence index from seq_2_MSA
                        needleman_list = needlemanWunsch(referenceSequences[i][1], referenceSequences[j][1], subMatrix)
                        temp_aligned_list.append(needleman_list)

            #Test display
            #print("\n\nTEMP ALIGNED LIST: ")
            #displayMatrix(temp_aligned_list)

            #Find the alignment that has the best optimal alignment score:
            #Assign an initial current max number to compare to:
            curr_max = temp_aligned_list[0][2] #Setting it to the first alignment's optimal score as a initial max value to compare
            #compare the curr_max of the temp_aligned list (lists of alignments with all possible Y's to X) and find the highest alignment score
            for b in temp_aligned_list:
                if b[2] > curr_max:
                    curr_max = b[2]

            
            best_score_pos = 0 #A variable that holds the index of the best score for the temp_aligned list
            #Get the best score's index from temp_aligned_list
            for c in range(len(temp_aligned_list)):
                if temp_aligned_list[c][2] == curr_max:
                    best_score_pos = c
                    break #once it finds the index of the best score, stop the loop


            #Once we obtained the best score out of all the alignments, copy all of the gaps from the highest score from temp_aligned_list[highest score position][2nd sequence]
            #Since we want to inject all the gaps within the MSA. Do not mess with the first seq_1 because its just a sequence, not an MSA like seq_2_MSA

            #For every other sequence that ISNT the best_score_pos, insert gaps from best_score_pos
            ## FIX THIS PART, WE NEED TO ADD THE GAP INSERTION FUNCTIONALITY FOR THIS TO WORK ## (10/29/23)
            ## I THINK I FIXED THE GAP COPYING PART, CHECK AND VERIFY WITH GRADER SCRIPT ## (10/30/23)
            #Copying gaps from MSA 1 best scoring to other original sequences in MSA 1
            for g in range(len(seq_1_MSA)):
                
                #print("CURRENT POSITION G -------------------")
                #print(g)
                #print(seq_1_MSA[g])
                #print("Current sequence that I am copying best sequence gaps over to: ")
                #print(referenceSequences[seq_1_MSA[g]])
                #print("Sequence to copy gaps: ")
                #print(temp_aligned_list[best_score_pos][0])


                #for each other sequence, copy all the gaps from highest scoring sequence to their other sequences
                letter_iterator = -1
                temp_build_new_seq = ""
                referenceSequences_length = len(referenceSequences[seq_1_MSA[g]][1]) #length of the current copy we are doing 
                
                
                #note: temp_aligned_list[best_score_pos][0] = best score sequence 1 to copy to everyone else
                #note: referenceSequences[seq_1_MSA[g]][1] = the current sequence to copy gaps to
                #THIS INNER FOR LOOP RE-BUILDS THE COPYING SEQUENCE BACKWARDS
                for letter_pos in range(len(temp_aligned_list[best_score_pos][0]) - 1, 0 - 1, -1): #Added -1 to add an offset to iterate through the letters backwards
                    if temp_aligned_list[best_score_pos][0][letter_pos] == "-":
                        #print("ADD A GAP")
                        temp_build_new_seq = temp_build_new_seq + "-"
                        continue

                    temp_build_new_seq = temp_build_new_seq + referenceSequences[seq_1_MSA[g]][1][letter_iterator] #Build the other sequence, but add a gap whenever it comes to it

                    letter_iterator = letter_iterator - 1

                    #print(str(letter_pos) +" --- " + str(temp_aligned_list[best_score_pos][0][letter_pos] + " --- " + str(temp_build_new_seq)))

                    if letter_iterator >= referenceSequences_length:
                        break
                
                #print("FINISHED SEQUENCE BUILT: ")
                #print(temp_build_new_seq[::-1])

                #Assign the new copied gap sequence into the main referenceSequences list
                referenceSequences[seq_1_MSA[g]][1] = temp_build_new_seq[::-1]

                #optional, ###I do not know if I need to make the copied gaps into X's --------------------------------------
                #Confirmed that this is what you do, you're supposed to replace the new gapped sequences with X's (10/31/23)
                referenceSequences[seq_1_MSA[g]][1] = referenceSequences[seq_1_MSA[g]][1].replace("-", "X")                                                               
                #The real question is, after I copy the gaps, do I make them into X's or not? idk

            #Copying gaps from MSA 2 best scoring to other original sequences in MSA 2
            for g in range(len(seq_2_MSA)):
                
                #print("CURRENT POSITION G -------------------")
                #print(g)
                #print(seq_2_MSA[g])
                #print("Current sequence that I am copying best sequence gaps over to: ")
                #print(referenceSequences[seq_2_MSA[g]])
                #print("Sequence to copy gaps: ")
                #print(temp_aligned_list[best_score_pos][1])


                #for each other sequence, copy all the gaps from highest scoring sequence to their other sequences
                letter_iterator = -1
                temp_build_new_seq = ""
                referenceSequences_length = len(referenceSequences[seq_2_MSA[g]][1]) #length of the current copy we are doing 
                
                
                #note: temp_aligned_list[best_score_pos][1] = best score sequence 1 to copy to everyone else
                #note: referenceSequences[seq_2_MSA[g]][1] = the current sequence to copy gaps to
                #THIS INNER FOR LOOP RE-BUILDS THE COPYING SEQUENCE BACKWARDS
                for letter_pos in range(len(temp_aligned_list[best_score_pos][1]) - 1, 0 - 1, -1): #Added -1 to add an offset to iterate through the letters backwards
                    if temp_aligned_list[best_score_pos][1][letter_pos] == "-":
                        #print("ADD A GAP")
                        temp_build_new_seq = temp_build_new_seq + "-"
                        continue

                    temp_build_new_seq = temp_build_new_seq + referenceSequences[seq_2_MSA[g]][1][letter_iterator] #Build the other sequence, but add a gap whenever it comes to it

                    letter_iterator = letter_iterator - 1

                    #print(str(letter_pos) +" --- " + str(temp_aligned_list[best_score_pos][1][letter_pos] + " --- " + str(temp_build_new_seq)))

                    if letter_iterator >= referenceSequences_length:
                        break
                
                #print("FINISHED SEQUENCE BUILT: ")
                #print(temp_build_new_seq[::-1])

                #Assign the new copied gap sequence into the main referenceSequences list
                referenceSequences[seq_2_MSA[g]][1] = temp_build_new_seq[::-1]

                #optional, ###I do not know if I need to make the copied gaps into X's --------------------------------------
                #Confirmed that this is what you do, you're supposed to replace the new gapped sequences with X's (10/31/23)
                referenceSequences[seq_2_MSA[g]][1] = referenceSequences[seq_2_MSA[g]][1].replace("-", "X")                                                               
                #The real question is, after I copy the gaps, do I make them into X's or not? idk
                        

            #After copying over gaps to each other sequence in the MSA, replace gaps with X in all alignments for MSA
            for e in temp_aligned_list:
                for f in range(len(e) - 1):
                    e[f] = e[f].replace("-", "X") 


            #Test displays
            #print("\n\nTHIS IS THE TEMP ALIGNED LIST")
            #print(temp_aligned_list[best_score_pos])
            #print("\n\nIndexSequencesWithT")
            #print(indexSequencesWithT)
            #print(indexSequencesWithT[8][1])
            #displayMatrix(referenceSequences)
            #print(seq_1_MSA)
            #print(seq_2_MSA)

            #Get the best_score_sequences and replace their original sequences
            temp_seqCombo_index = [] #List to keep track of what sequence got compared to what
            for f in range(len(seq_1_MSA)):
                for g in range(len(seq_2_MSA)):
                    temp_seqCombo_index.append((f, g)) #make an index combination list
                    #print(str(f) + " ---- " + str(g)) 

            #Test display for combination list
            #print(temp_seqCombo_index)
            #print(temp_seqCombo_index[best_score_pos][0])
            print(seq_1_MSA[temp_seqCombo_index[best_score_pos][0]])
            print(seq_2_MSA[temp_seqCombo_index[best_score_pos][1]])
            
            #Replace the value of the seq_1_MSA that had the best_score_pos
            referenceSequences[seq_1_MSA[temp_seqCombo_index[best_score_pos][0]]][1] = temp_aligned_list[best_score_pos][0]
            
            #replace the value of the seq_2_MSA that had the best_score_pos
            referenceSequences[seq_2_MSA[temp_seqCombo_index[best_score_pos][1]]][1] = temp_aligned_list[best_score_pos][1]

            #Test display to check if the newly aligned sequences replaced the old sequences in the reference sequences list
            displayMatrix(referenceSequences)
            
        
       
    #After doing the guide tree, turn all the X back into gaps
    for n in referenceSequences:
        n[1] = n[1].replace("X", "-")

    #print(referenceSequences)


    #Then assign MSA to referenceSequences, and return the MSA
    MSA = list(referenceSequences)

    return MSA #return the finished MSA alignment sequences list????/



#   Function Name: insDash
#   Parameters: 
#               givenStr - The given string I want to insert a gap into (in this case, the sequence)
#               insertIndex - The position from the best scoring alignment out of the MSA to insert gap into (Copying gaps over to the other sequences)
#   Returns: sequence with gap inserted at position
#   Purpose: Used as a utility function to insert a gap
def insDash(givenStr, insertIndex):
    return givenStr[:insertIndex] + "-" +givenStr[insertIndex:]
    
    


#   Function Name: sumOfPairs
#   Parameters: 
#               MSA - list of lists containing the multiple sequence alignment
#               subMatrix - used for scoring based on .mtx file
#   Returns: final_score - The sum of all the columns in a aligned MSA sequence scores with all possible combinations added together to get the total MSA score
#   Purpose: Used to compute the final MSA score
def sumOfPairs(MSA, subMatrix):
    final_score = 0 #Variable to hold the final score

    #Test display to see subMatrix
    #print(subMatrix)

    #For every combination of all sequences, sum up all the columns based on their score from subMatrix (.mtx file) and add them all together to get the final score
    for i in range(len(MSA)): #Iterate to the total length of the number of sequences from the sequences list 
        #Note: all sequences in the sequences list of lists are in position 1, because position 0 is the header for each sequence
        sequence_1  = MSA[i][1] #get sequence 1 to compare to every other sequence

        #Test display for the current sequence. Un-comment if you want to see the process of comparing each sequence to each other
        #print("\n\n------CURRENT SEQUENCE BEING COMPARED TO OTHER SEQUENCES "+str(MSA[i]))

        for j in range (i + 1,len(MSA)): #Iterate to the total length of the number of sequences to compare each combination of sequences together. i + 1 offset because we do not need to compare C to A if we already compared A to C
            sequence_2 = MSA[j][1] #get sequence 2 to compare to sequence 1, where sequence 1 has to be compared to ALL sequences in order to iterate to the next sequence combination
            
            #Test display to see if we are getting sequence 2 
            #print(MSA[j]) 

            #compare columns of current sequence 1 to columns of current sequence 2 and add them to the final score by summing their values up based on subMatrix
            for letter_column in range(len(sequence_1)):
                final_score = final_score + (int(subMatrix.get((sequence_1[letter_column], sequence_2[letter_column]))))
                #Test to see what letters (columns) are being scored together
                #print(str(sequence_1[letter_column]) + "    VS      "+str(sequence_2[letter_column]) + "    SCORE FOR THIS COMBINATION = "+ str(int(subMatrix.get((sequence_1[letter_column], sequence_2[letter_column])))))

    #Test print to see the final_score
    #print(final_score)

    return final_score



def main():
    parser = argparse.ArgumentParser() #Establish the parser
    sorted_sequences = [] #List to represent the finished, sorted list_of_lists from the readFastA function
    score_dict = {} #List to represent the scores from each letter in the dictionary from the score matrix BLOSUM50 or BLOSUM62
    distanceMatrix = [] #List to represent the distance matrix from running needleman wunsch for all possible combinations of sequences
    # Adding an argument/option to the command line
    # -i is used to specify the file path to the input FASTA file
    # called in this format:
    # python3 <python file name> -i <input file path>
    parser.add_argument("-i", "--inputpath", help = "This option retrieves the file path to the starting cellular matrix", type = argparse.FileType('r'))
    
    # -o is used to specify the file path to the output file
    # I think that he wants us to have the capability to write out the output to an output file
    # called in similar format to our -i:
    # python3 <python file name> -o <output file path> 
    parser.add_argument("-o", "--outputpath", help = "This option retrieves the file path for the final output file", type = argparse.FileType('w'))
    #Note: I verified that option -o specifies the and writes to a new file with the output, but you can specify the filepath to where the output goes. Requires this code to be ran in the HPCC grading script (verified on 9/18/2023)
    

    # -s is used to specify the file path to the score matrix
    # calledin format:
    # python3 <python file name> -s <file path to the score matrix>
    parser.add_argument("-s", "--scorepath", help = "This operation retrieves the file path for the score matrix", type = argparse.FileType('r'))

    args = parser.parse_args() #runs parser and parses the arguments from the command line input

    #Calls the readFastA function to read all sequences, sort them, then return the sorted list of lists
    sorted_sequences = readFastA(args.inputpath.readlines()) 

    sorted_sequences2 = list(sorted_sequences)
    #Reminder:
    #The first line of the output printed to the screen during execution MUST be in the format:
    # "Assignmetn 1 :: R#"
    #Where R# is the TTU R#. 
    #This is REQUIRED for the grading script
    #Displays the required 1st line printed output fromt he program output specifications:
    print("Assignment 1 :: R11641653")

    #note: BLOSUM and the nucleotide.mtx files are also known as substitution matricies, so I am renaming the variables for clarity
    SUBSTITUTION_MATRIX_FILE = open(args.scorepath.name) #Opens the Score Matrix file to be used in readScoreMatrix


    sub_dict = dict(readScoreMatrix(SUBSTITUTION_MATRIX_FILE)) #Calls the readScoreMatrix and generates a dictionary based on the yield tuples of letter pairs WITH their given matching score
    

    #print(sorted_sequences)

    #build the distance diagonal matrix first by passing a list of all the sequences from the .fna file, and the substitution matrix (BLOSUM, or nucleotide.mtx or etc.)
    #Note: From my last assignment 2, I made the substitution matrix a dictionary so that I can find whatever match for each nucleotide combination I needed (ex: what score is (A with C) = 12 or something defined in the given .mtx file)
    D = buildDistanceMatrix(sorted_sequences, sub_dict) #From the assignment_4 hints, where D is the distance matrix from the function
    #Returns a tuple of the distance matrix in D[0] and a dictionary for the distance matrix in D[1]

    #Test to see if the distance matrix (D) was generated properly
    #print("THIS IS THE D (DISTANCE MATRIX AFTER BUILD DISTANCE MATRIX FUNCTION)")
    #print(D)

    
 

    #After getting the distance matrix, construct the guide tree of of that distance matrix
    T = constructGuideTree(sorted_sequences, sub_dict, D)

    #(10/24/23) I think I finished the constructGuideTree. I didn't use a dictionary at all, just raw values from the matrix. 
    #           I also put a lot of explaination on how it works for myself

    #(10/27/23) I fixed the guide tree output and now it works based off of index values from the sorted_sequences (sequences) list. 
    #           this will mainly be used to index through the guide tree and do the Needleman-wunsch in the correct order. 
    #           Note: I did not fix the sequences just yet. The guide tree is correct but the sequences inserted MSA are wrong. I need to fix

    #(10/27/23) Okay, I think the construct guide tree works now, now with the correct guide tree construction (I think...)
    #           Work on building the MSA to check order and the alignment sum of pairs score later.
    #           The sorted_sequences list should be modified to have the guide tree in there in order to index through it

    #Test display to also show the added sequence indexes into the original list to index through the guide tree properly 
    print("FINAL SEQUENCES LIST USED FOR INDEXING THROUGH THE GUIDE TREE: ")
    displayMatrix(sorted_sequences)

    #DISPLAY TEST FOR GUIDE TREE (T)
    print("\nFINAL GUIDE TREE GENERATED:")
    for i in T: 
       print(i)

    
    
    #After constructing the guide tree, do the MSA in the order that the tree tells us to
    MSA = getMSA(sorted_sequences, sorted_sequences2, sub_dict, T)

    #(10/29/23) I think I got the MSA to work at least for the sequence vs MSA part. I still need to build MSA to MSA part
    #Work on sum of pairs for now

    #The output MSA looks incorrect, so many gaps. Ask Mr. Rees about it when I get the chance. For now, I am stopping (10/24/23)

    #DISPLAY TEST FOR MSA 
    print("\nFINAL MSA GENERATED: ")
    for i in MSA:
        print(i)


    #After generating the MSA, get the final score with sumOfPairs:
    final_score = sumOfPairs(MSA, sub_dict)


    #After getting the final score of the MSA, add the final score to all the sequence headers in MSA
    #Writes out the output to the specified output location and file name from -o argument from command line, adding a new line for each sequence
    for elements in MSA:
        for h in elements:
            if h[0] == ">": #If h is a sequence header add the final_score to the header
                args.outputpath.write(h+"; score="+str(final_score)+ "\n")

            else: #If h is not a sequence header and is a sequence, then write out like normal
                args.outputpath.write(h + "\n")

    #**THIS WAS FOR ASSIGNMENT 2
    #sequence_1 = sorted_sequences[0][1] #Gets the second sequence from FASTA 

    #sequence_2 = sorted_sequences[1][1] #Gets the first sequence from FASTA


    #Test display to check for sequences
    #print("INPUT SEQUENCES: \n")
    #print("INPUT 1: "+sequence_1+"\n")
    #print("INPUT 2: "+sequence_2+"\n")
    #print(sorted_sequences)
    #print(sub_dict)
    
   
    #run needleman-wunsch algorithm and align the two given sequences with the given score matrix, then return a tuple of the aligned matricies
    #aligned_sequences = needlemanWunsch(sequence_1, sequence_2, sub_dict)

    #note: sequence_2 goes first because the "readFastA" sorts the 

    #Test display for tuple output check
    #print(aligned_sequences)
    #print(aligned_sequences[0])
    #print(aligned_sequences[1])

    #Take the output tuple of "aligned_sequences" and reassign them back into the list with appended score values and new aligned sequences
    #sorted_sequences[0][1] = aligned_sequences[0] #Reassign the sequence 1 array with the aligned version of it back to position 0 since the algorithm always puts sequence/longest sequence at position 0
    #sorted_sequences[1][1] = aligned_sequences[1] #Reassign the sequence 2 array with the aligned version of it back to position 1 since the algorithm has the shortest sequence at position 1
    #sorted_sequences[1][0] = sorted_sequences[1][0]+"; score="+str(aligned_sequences[2]) #append the optimal alignment score to the end of the their sequence headers
    #sorted_sequences[0][0] = sorted_sequences[0][0]+"; score="+str(aligned_sequences[2]) #append the optimal alignment score to the end of the their sequence headers
    
    #Test display to see sorted_sequences
    #print(sorted_sequences)
    #**THIS WAS FOR ASSIGNMENT 2

    


#Runs the Main() function
if __name__ == '__main__': 
    main()