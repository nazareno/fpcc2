#!/usr/bin/python

######################################
# Fisher's Exact Test - Significance Test
# From: Statistics is Easy! By Dennis Shasha and Manda Wilson
# 
# Assuming that our tea taster can not identify which cups of tea were milk
# first and which were not, tests to see the probability of getting an outcome 
# as extreme or more extreme than the observed one (i.e. more correct than in the 
# observed case) by chance alone.  
#
# This is a one-tailed test.
#
# Author: Manda Wilson
#
# Example of FASTA formatted input file: 
# >tea taster claimed milk first
# 3 1
# >tea taster claimed tea first
# 2 4
#
# The above generalizes to:
#
# >variable 1
# a b
# >variable 2
# c d
#
# Pseudocode:
#
# 1. Calculate Fisher's Exact Test for the observed values (in our example it is 0.2619).
#    a. Calculate a + b, c + d, a + c, b + d, and a + b + c + d  
#    b. Set observed_prob_tail = 0
#    c. Compute all possible outcomes (matricies) 
#       i. For a' = 0 to a + b + c + d
#          A. b' = (a + b) - a' (if this is impossible, skip to next a')
#          B. c' = (a + c) - a' (if this is impossible, skip to next a')
#          C. d' = (c + d) - c' (if this is impossible, skip to next a')
#          D. If this matrix is "as extreme", or "more extreme" than the observed compute the
#             probability of getting this matrix by chance and add this probability
#             to observed_prob_tail.  In this case, the matrix is "more extreme" if 
#             a' + d' > a + b, i.e., the tea taster was correct more often than observed.
#             If the test is a two-tailed test, "more extreme" matricies include
#             those where a' + d' > a + b as well as those where b' + c' > a + b
#             (i.e. the matrix is more unbalanced than the observed).
#             - Compute the probability of getting the matrix
#               a' b'
#               c' d'
#               prob = ((a' + b')!(c' + d')!(a' + c')!(b' + d')! / (a'!b'!c'!d'!n'!)
#               Where n = a' + b' + c' + d'
#               Example: ((3 + 1)!(2 + 4)!(3 + 2)!(1 + 4)!) / (3!1!2!4!10!) = 10/42
#             - observed_prob_tail += probability from step (1ciD)
#
# 2. Set a counter to 0, this will count the number of times we get a fisher's exact test 
#    less than or equal to 0.2619 (less than because a smaller value is more unlikely).
#
# 3. Do the following 10,000 times:
#    a. Create a new matrix of counts, preserving the marginals from our original matrix
#       i.  While there are more rows:
#              I. While there are more columns:
#                     * If we ARE NOT at the last element in the row
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
#    b. Calculate Fisher's Exact Test on the results from step (3a), just as we did in step (1).
#    c. If the result from step (3b) is less than or equal to our observed Fisher's Exact Test (0.2619),
#       increment our counter from step (2).
#
# 4. counter / 10,000 equals the probability of getting a Fisher's Exact Test less than
#    or equal to our observed probability (0.2619)
#
######################################

import random

######################################
#
# Adjustable variables
#
######################################

input_file = "fishersexact.vals"

######################################
#
# Subroutines
#
######################################

# takes a list of values, plus row and column totals
# For an m x n matrix, values must be ordered:
# row 1, columns 1 - n, row 2, columns 1 - n, ... row m, columns 1 - n
# shuffles values, preserving row and column totals
def shuffle(orig_row_totals, orig_column_totals):
        # start with an empty m x n matrix
        # fill it with randomly generated counts
        # that preserve the row and column totals

        current_row = 0
        # keeps track of the available values for this row
        available_row_vals = orig_row_totals[:]
        # keeps track of the available values for this column
        available_column_vals = orig_column_totals[:]

        new_counts = []

        while current_row < len(orig_row_totals):
                current_column = 0
                while current_column < len(orig_column_totals):
                        # if we are not at the last element in either the row or the column
                        if current_row < len(orig_row_totals) - 1:
                                if current_column < len(orig_column_totals) - 1:
                                        # get a random number between 0 and
                                        # min(available_row_vals[current_row],
                                        #               available_column_vals[current_column])
                                        max_val = min(available_row_vals[current_row], available_column_vals[current_column])
                                        new_val = random.randint(0, max_val)
                                        # put this value in the new matrix
                                        new_counts.append(new_val)
                                else: # we are at the last column, this value must be whatever is available for this row
                                        new_counts.append(available_row_vals[current_row])
                        else: # we are at the last row, this value must be whatever is available for this column
                                new_counts.append(available_column_vals[current_column])
                        # remove this amount from both the available row and column values
                        available_row_vals[current_row] -= new_counts[-1]
                        available_column_vals[current_column] -= new_counts[-1]
                        current_column += 1
                current_row += 1
        return new_counts

def factorial(n):
	prod = 1
	for i in range(n):
		prod = prod * (i + 1)
	return prod

def prob_of_matrix(a, b, c, d):
	return (factorial(a + b) * factorial(c + d) * factorial(a + c) * factorial(b + d)) / float(factorial(a) * factorial(b) * factorial(c) * factorial(d) * factorial(a + b + c + d))

def fishers_exact_test(a, b, c, d):
	# now we have to figure out possible outcomes
	# that are more extreme than ours
	# and sum the probability of each
	# this is the part of the code that should be tailored
	# for your definition of "more extreme"
	# here we are doing a one-tailed test
	# where "more extreme" means, more correct answers than
	# what was observed 
	# this translates to any matrix where a + d is larger than ours

	prob_tail = 0.0 

	for a_prime in range(a + b + c + d + 1):
        	b_prime = a + b - a_prime
        	if b_prime >= 0: 
                	c_prime = a + c - a_prime
                	if c_prime >= 0:
                        	d_prime = c + d - c_prime
                        	if d_prime >= 0:
                                	# this matrix is valid
                                	# now check if it is "more extreme" than the observed
                                	if a_prime + d_prime >= a + d:
                                        	prob_tail = prob_tail + prob_of_matrix(a_prime, b_prime, c_prime, d_prime)
	return prob_tail

######################################
#
# Computations
#
######################################

# set to invalid counts
a = -1
b = -1
c = -1
d = -1

# file must be in FASTA format
infile=open(input_file)
for line in infile:
	if not line.isspace() and not line.startswith('>'):
		# this is one row
		if a == -1 and b == -1:
			(a, b) = map(int,line.split())
		else:
			(c, d) = map(int,line.split())
infile.close()

observed_prob_tail = fishers_exact_test(a, b, c, d)

count = 0
num_runs = 10000

for i in range(num_runs):
        [a_prime, b_prime, c_prime, d_prime] = shuffle([a + b, c + d], [a + c, b + d])
        prob_tail = fishers_exact_test(a_prime, b_prime, c_prime, d_prime)
        if (prob_tail <= observed_prob_tail):
                count = count + 1

######################################
#
# Output
#
######################################

print "Observed Fisher's Exact Test: %.4f" % observed_prob_tail
print count, "out of 10000 experiments had a Fisher's Exact Test less than or equal to %.4f" % observed_prob_tail
print "Probability that chance alone gave us a Fisher's Exact Test",
print "of %.4f" % observed_prob_tail, "or less is", (count / float(num_runs))
