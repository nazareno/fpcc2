#!/usr/bin/python

######################################
# Difference between Two Means Significance Test
# From: Statistics is Easy! By Dennis Shasha and Manda Wilson
# 
# Assuming that there is no significant difference in the means of 
# the two samples, tests to see the probability of getting a difference 
# greater than or equal to the observed difference in the means by chance 
# alone.  Uses shuffling & bootstrapping to get a distribution to compare 
# to the observed statistic.
# 
# Author: Manda Wilson
#
# Example of FASTA formatted input file (this is only for 2 groups): 
# >placebo_vals
# 54 51 58 44 55 52 42 47 58 46
# >drug_vals
# 54 73 53 70 73 68 52 65 65
#
# Pseudocode:
#
# 1. Measure the difference between the two group means.  The difference in means is measured
#    by (sum(grpA) / len(grpA)) - (sum(grpB) / len(grpB)).  In this example the difference between
#    the two group means is 12.97.
#
# 2. Set a counter to 0, this will count the number of times we get a difference
#    between the means greater than or equal to 12.97.  
#
# 3. Do the following 10,000 times:
#    a. Shuffle the original measurements.  To do this:
#       i. put the values from all the groups into one array but remembering the start and end indexes of each group
#       ii. shuffle the values in the array, effectively reassigning the values to different groups
#    b. Measure the difference between the two group means, just as we did in step (1).
#    c. If the difference from step (3b) is greater than or equal to 12.97, increment our counter 
#       from step (2). Note: if our original difference between the means were a negative value 
#       we would check for values less than or equal to that value.
#
# 4. counter / 10,000 equals the probability of getting our observed difference of two means greater than
#    or equal to 12.97, if there is in fact no significant difference.
#
######################################

import random

######################################
#
# Adjustable variables
#
######################################

input_file = "Diff2Mean.vals"

######################################
#
# Subroutines
#
######################################

# takes a list of groups (two or more)
# pools all values, shuffles them, and makes new groups 
# of same size as original groups
# returns these new groups
# example of shuffle with more than two groups: OneWayAnovaSig.py
def shuffle(grps):
	num_grps = len(grps)
	pool = []
	
	# pool all values
	for i in range(num_grps):
		pool.extend(grps[i])
	# mix them up
	random.shuffle(pool)
	# reassign to groups of same size as original groups
	new_grps = []
	start_index = 0
	end_index = 0
	for i in range(num_grps):
		end_index = start_index + len(grps[i])
		new_grps.append(pool[start_index:end_index])
		start_index = end_index
	return new_grps

# subtracts group a mean from group b mean and returns result
def meandiff(grpA, grpB):
	return sum(grpB) / float(len(grpB)) - sum(grpA) / float(len(grpA))

######################################
#
# Computations
#
######################################

# list of lists
samples = [] 
a = 0
b = 1

# file must be in FASTA format
infile=open(input_file)
for line in infile:
	if line.startswith('>'):
		# start of new sample
		samples.append([])
	elif not line.isspace():
		# line must contain values for previous sample
		samples[len(samples) - 1] += map(float,line.split())
infile.close()

observed_mean_diff = meandiff(samples[a], samples[b])

count = 0
num_shuffles = 10000

for i in range(num_shuffles):
	new_samples = shuffle(samples)
	mean_diff = meandiff(new_samples[a], new_samples[b])
	# if the observed difference is negative, look for differences that are smaller
	# if the observed difference is positive, look for differences that are greater
	if observed_mean_diff < 0 and mean_diff <= observed_mean_diff:
		count = count + 1
	elif observed_mean_diff >= 0 and mean_diff >= observed_mean_diff:
		count = count + 1

######################################
#
# Output
#
######################################

print "Observed difference of two means: %.2f" % observed_mean_diff 
print count, "out of", num_shuffles, "experiments had a difference of two means",
if observed_mean_diff < 0:
	print "less than or equal to",
else:
	print "greater than or equal to",
print "%.2f" % observed_mean_diff, "."
print "The chance of getting a difference of two means",
if observed_mean_diff < 0:
	print "less than or equal to",
else:
	print "greater than or equal to",
print "%.2f" % observed_mean_diff, "is", (count / float(num_shuffles)), "."
