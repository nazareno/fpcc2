#!/usr/bin/python

######################################
# Difference between Two Means Confidence Interval
# From: Statistics is Easy! By Dennis Shasha and Manda Wilson
# 
# Uses shuffling & bootstrapping to get a 90% confidence interval for the difference between two means.
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
#    by sum(grpA) / len(grpA) - sum(grpB) / len(grpB).  In this example the difference between
#    the two group means is 12.97.
#
# 2. Do the following 10,000 times:
#    a. For each sample we have get a bootstrap sample:
#       i. Create a new array of the same size as the original sample
#       ii. Fill the array with randomly picked values from the original sample (randomly picked with replacement)
#    b. Measure the difference between the two bootstrap group means, just as we did in step (1)
#       with the original samples.
#
# 3. Sort the differences computed in step (2).
# 
# 4. Compute the size of each interval tail.  If we want a 90% confidence interval, then 1 - 0.9 yields the
#    portion of the interval in the tails.  We divide this by 2 to get the size of each tail, in this case 0.05.
#
# 5. Compute the upper and lower bounds.  To get the lower bound we multiply the tail size by the number of 
#    bootstraps we ran and round this value up to the nearest integer (to make sure this value maps
#    to an actual boostrap we ran).  In this case we multiple 0.05 by 10,000 and get 500, which rounds up to 500.
#    
#    To compute the upper bound we subtract the tail size from 1 and then multiply that by the number of bootstraps
#    we ran.  Round the result of the last step down to the nearest integer.  Again, this is to ensure this value
#    maps to an actual bootstrap we ran.  We round the lower bound up and the upper bound down to reduce the
#    confidence interval size, so that we can still say we have as much confidence in this result.
#
# 6. The bootstrap values at the lower bound and upper bound give us our confidence interval.
#
###################################### 

import random
import math
import sys

######################################
#
# Adjustable variables
#
######################################

input_file = "Diff2Mean.vals"
conf_interval = 0.9

######################################
#
# Subroutines
#
######################################

# x is an array of sample values
# returns a new array with randomly picked 
# (with replacement) values from x
def bootstrap(x):
	samp_x = []
	for i in range(len(x)):
		samp_x.append(random.choice(x))
	return samp_x

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

num_resamples = 10000   # number of times we will resample from our original samples
out = []                # will store results of each time we resample

for i in range(num_resamples):
	# get bootstrap samples for each of our groups
	# then compute our statistic of interest
	# append statistic to out
	bootstrap_samples = []  # list of lists
	for sample in samples:
		bootstrap_samples.append(bootstrap(sample))
	# now we have a list of bootstrap samples, run meandiff
	out.append(meandiff(bootstrap_samples[a], bootstrap_samples[b]))

out.sort()

tails = (1 - conf_interval) / 2

# in case our lower and upper bounds are not integers,
# we decrease the range (the values we include in our interval),
# so that we can keep the same level of confidence
lower_bound = int(math.ceil(num_resamples * tails))
upper_bound = int(math.floor(num_resamples * (1 - tails)))

######################################
#
# Output
#
######################################

# print observed value and then confidence interval
print "Observed difference between the means: %.2f" % observed_mean_diff
print "We have", conf_interval * 100, "% confidence that the true difference between the means",
print "is between: %.2f" % out[lower_bound], "and %.2f" % out[upper_bound]
