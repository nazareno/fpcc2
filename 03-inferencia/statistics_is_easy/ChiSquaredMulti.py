#!/usr/bin/python

######################################
# Multi-Variable Chi-Squared Significance Test
# From: Statistics is Easy! By Dennis Shasha and Manda Wilson
# 
# Assuming that health does not influence wealth, tests to see the probability
# of getting a chi-squared by change alone greater than or equal to the one observed.
# Uses shuffling to get a distribution to compare 
# to the observed chi-squared.
# 
# Author: Manda Wilson
#
# Example of FASTA formatted input file: 
# >sick
# 20 18 8
# >healthy
# 24 24 16
#
# Pseudocode:
#
# 1. Calculate chi-squared for the observed values (in our example it is 0.97).  
#    a. Calculate the expected counts for each category (it is row total / total * column total)
#       Example: for the sick/poor category we expected (46 / 110) * 44 = 18.48 
#    b. For each category (in this example sick/poor, sick/middle, sick/rich, healthy/poor, 
#       healthy/middle, and healthy/rich are the categories):
#       i. Subtract the expected count from the observed count
#       ii. Square the result of step (1bi)
#       iii. Divide the result of step (1bii) by the expected count
#    c. Sum the results of step (1b), this is r, our observed chi-squared value
#
# 2. Set a counter to 0, this will count the number of times we get a chi-squared
#    greater than or equal to 0.97.  
#
# 3. Do the following 10,000 times:
#    a. Create a new matrix of counts, preserving the marginals from our original matrix
#       i.  While there are more rows:
#	       I. While there are more columns:
#		      * If we ARE NOT at the last element in the row 
#             * If we ARE NOT at the last element in the column
#                - Pick a random number between 0
#                   and min(available_in_this_row, 
#                      available_in_this_column) (inclusive)
#                - Store this value in our new maxtrix in the 
#                     current row, current column position
#             * If we ARE at the last element in the column
#                - Store whatever is available_in_this_column 
#                   in the current row, current column position
#             * If we ARE at the last element in the row
#                - Store whatever is available_in_this_row 
#                   in the current row, current column position
#             * Subtract whatever is stored in the current row, 
#                current column position from the available_in_this_column
#                as well as from available_in_this_row
#
#    b. Calculate chi-squared on the results from step (3a), just as we did in step (1).
#    c. If the result from step (3b) is greater than or equal to our observed chi-squared (0.97), 
#       increment our counter from step (2). 
#
# 4. counter / 10,000 equals the probability of getting a chi-squared greater tha or equal to
#    0.97, if wealth has no influence on health
#
######################################

import random

######################################
#
# Adjustable variables
#
######################################

input_file = "chisquaredmulti.vals"

######################################
#
# Subroutines
#
######################################

# takes a list of values, plus the number of rows and columns
# For an m x n matrix, values must be ordered in a 1 dimensional array:
# row 1, columns 1 - n, row 2, columns 1 - n, ... row m, columns 1 - n
# For example, for health (sick/healthy) and wealth (poor, middle, rich),
# (2 variables with 2 and 3 categories respectively)
# If input is:
#	Poor  	Middle  	Rich
# Sick 	20 	18 	8
# Healthy 	24 	24 	16
# We expect the following 1 dimensional array
# [20, 18, 8, 24, 24, 16]
# shuffles values, preserving row and column totals
def shuffle(observed, num_rows, num_cols):

	# this array will store values for the first variable
	# in our example, it will store health values
	# 0 = sick, 1 = healthy
	row_vals = []
	
	# this array will store values for the second variable
	# in our example, it will store wealth values
	# 0 = poor, 1 = middle, 2 = rich
	col_vals = []

	for r in range(num_rows):
        	for c in range(num_cols):
			count_in_cell = int(observed[(r * num_cols) + c])
			# there must be a better way...
			for i in range(count_in_cell):
				row_vals.append(r)
				col_vals.append(c)

	# now shuffle one variable (breaking any association between the two) 
	# in our example we are shuffling wealth lables
	random.shuffle(col_vals)

	# reassemble counts like original one-dimensional array
	new_counts = []

	for r in range(num_rows):
                for c in range(num_cols):
			new_counts.append(0)

	for i in range(len(row_vals)):
		# use values in row_vals and col_vals to insert into this array
		new_counts[(row_vals[i] * num_cols) + col_vals[i]] += 1
	return new_counts
				
def chisquared(expected, observed):
	count = len(expected)
	total = 0
	for i in range (count):
		total += ((observed[i] - expected[i])**2) / float(expected[i])
	return total

######################################
#
# Computations
#
######################################

observed = []
column_totals = []
row_totals = []

# file must be in FASTA format
infile=open(input_file)
for line in infile:
	if not line.isspace() and not line.startswith('>'):
		# this is one row
		row = map(float,line.split())
		row_totals.append(sum(row))
		if len(column_totals) == 0:
			column_totals = row
		else: 
			for c in range(len(column_totals)):
				column_totals[c] += row[c]
		observed += row
infile.close()

total = sum(column_totals)
expected = []

# calculate expected values based on row and column totals
for r in range(len(row_totals)):
	for c in range(len(column_totals)):
		expected.append((row_totals[r] / float(total) * column_totals[c]))

observed_chi_squared = chisquared(expected, observed)

count = 0
num_runs = 10000

for i in range(num_runs):
	shuffled_observed = shuffle(observed, len(row_totals), len(column_totals))
	chi_squared = chisquared(expected, shuffled_observed)
	if (chi_squared >= observed_chi_squared):
		count = count + 1

######################################
#
# Output
#
######################################

print "Observed chi-squared: %.2f" % observed_chi_squared
print count, "out of", num_runs, "experiments had a chi-squared greater than or equal to %.2f" % observed_chi_squared
print "Probability that chance alone gave us a chi-squared greater than or equal to", 
print "%.2f" % observed_chi_squared, "is", (count / float(num_runs))

	
