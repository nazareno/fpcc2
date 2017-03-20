#!/usr/bin/python

######################################
# One-Way ANOVA Significance Test 
# From: Statistics is Easy! By Dennis Shasha and Manda Wilson
# 
# Assuming that there is no difference in the three drugs used, tests to see the probability
# of getting a f-statistic by chance alone greater than or equal to the one observed.
# Uses shuffling to get a distribution to compare 
# to the observed f-statistic.
# 
# Author: Manda Wilson
#
# Example of FASTA formatted input file: 
# >grp_a
# 45 44 34 33 45 46 34
# >grp_b
# 34 34 50 49 48 39 45
# >grp_c
# 24 34 23 25 36 28 33 29
#
# Pseudocode:
#
# 1. Calculate f-statistic for the observed values (in our example it is 11.27).  
#    a. Initialize within sum of squares (wss), total sum (ts),
#       total count (tc), and between sum of squares (bss) to 0
#    b. For each group:
#	i. Within group mean (wgm) = sum of group values / number of values in group
#	ii. Store the number of values in group
#	iii. Sum the following: for each value in the group
#            I. Subtract the value from wgm
#            II. Square the result of step (1biiiI)
#            III.  Add the result of step (1biiiII) to wss
#    c. Total mean (tm) = ts / tc
#    d. For each group: 
#       i. Subtract the wgm from tm
#       ii. Square the result of step (1di)
#       iii. Multiply the result of step (1dii) by 
#            the number of values in that group (stored in step (1bii))
#       iv. Add the result of step (1diii) to bss
#    e. Between degrees of freedom (bdf) = number of groups - 1
#    f. Within degrees of freedome (wdf) = tc - number of groups
#    g. Within group variance (wgv) = wss / wdf
#    h. Between group variance (bgv) = bss / bdf
#    i. f-statistic = wgv / bgv
#
# 2. Set a counter to 0, this will count the number of times we get a f-statistic 
#    greater than or equal to 11.27.  
#
# 3. Do the following 10,000 times:
#    a. Shuffle the observed values. To do this:
#       i. Put the values from all the groups into one array
#       ii. Shuffle the pooled values
#       iii. Reassign the pooled values to groups of the same size as the original groups
#    b. Calculate f-statistic on the results from step (3a), just as we did in step (1)
#    c. If the result from step (3b) is greater than or equal to our observed f-statistic (11.27), 
#       increment our counter from step (2). 
#
# 3. counter / 10,000 equals the probability of getting a f-statistic greater than or equal to
#    11.27, assuming there is no difference between the groups
#
######################################

import random

######################################
#
# Adjustable variables
#
######################################

input_file = 'OneWayAnova.vals'

######################################
#
# Subroutines
#
######################################

# a list of lists
def shuffle(grps):
	num_grps = len(grps)
	pool = []
	
	# throw all values together
	for i in range(num_grps):
		pool.extend(grps[i])
	# mix them up
	random.shuffle(pool)
	# reassign to groups
	new_grps = []
	start_index = 0
	end_index = 0
	for i in range(num_grps):
		end_index = start_index + len(grps[i])
		new_grps.append(pool[start_index:end_index])
		start_index = end_index
	return new_grps


def sumofsq(vals, mean):
	# the sum of squares for a group is calculated 
	# by getting the sum of the squared difference 
	# of each value and the mean of the group that value belongs to
	count = len(vals)
	total = 0
	for i in range (count):
		diff_sq = (vals[i] - mean)**2
		total += diff_sq
	return total

def weightedsumofsq(vals, weights, mean):
	count = len(vals)
	total = 0
	for i in range (count):
		diff_sq = (vals[i] - mean)**2
		total += weights[i] * diff_sq
	return total

# expects list of lists
def onewayanova(grps):
	num_grps = len(grps)
	within_group_means = []
	grp_counts = []
	within_ss = 0
	total_sum = 0
	total_count = 0
	for i in range (num_grps):
		grp = grps[i]
		grp_count = len(grp)
		grp_sum = sum(grp)
		within_group_mean = grp_sum / float(grp_count)

		grp_counts.append(grp_count)

		total_count += grp_count	# add to total number of vals
		total_sum += grp_sum

		within_group_means.append(within_group_mean) 

		# to get the within group sum of squares:
		# sum the following: for every element in the overall group
		# subtract that element's value from that element's group mean
		# square the difference
		# get within group sum of squares
		# this is calculated by summing each group's sum of squares
		within_ss += sumofsq(grp, within_group_mean) 
		
	total_mean = total_sum / total_count 
	
	# to get the between group sum of squares:
	# sum the following: for every element in the overall group
	# subtract that element's group mean from the overall group mean
	# square the difference
	# grp_counts are used as weights
	between_ss = weightedsumofsq(within_group_means, grp_counts, total_mean)

	# now we want to find out how different the groups are between each other
	# compared to how much the values vary within the groups
	# if all groups vary a lot within themselves, 
	# and there is no significant difference
	# between the groups, then we expect the differences between 
	# the groups to vary by about the same amount
	# so lets get the ratio of the between group variance and 
	# the within group variance
	# if the ratio is 1, then there is no difference between the groups
	# if it is significantly larger than one, 
	# then there is a significant difference between the groups
	# remember: even if the groups are significantly different, 
	# we still won't know which groups are different

	# the between group degrees of freedom
	# is equal to the number of groups - 1
	# this is because once we know the number of groups - 1, 
	# we know the last group
	between_df = len(grp_counts) - 1

	# the within group degrees of freedom
	# is equal to the total number of values minus the number of groups
	# this is because for each group, once we know the count - 1 values, 
	# we know the last value for that group
	# so we lose the number of groups * 1 degrees of freedom
	within_df = total_count - num_grps

	within_var = within_ss / within_df
	between_var = between_ss / between_df

	f_stat = between_var / within_var

	return f_stat

######################################
#
# Computations
#
######################################

# list of lists
samples = [] 

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

observed_f_statistic = onewayanova(samples)

count = 0
num_shuffles = 10000

for i in range(num_shuffles):
	new_samples = shuffle(samples)
	f_statistic = onewayanova(new_samples)
	if (f_statistic >= observed_f_statistic):
		count = count + 1

######################################
#
# Output
#
######################################

print "Observed F-statistic: %.2f" % observed_f_statistic 
print count, "out of 10000 experiments had a F-statistic greater than or equal to %.2f" % observed_f_statistic
print "Probability that chance alone gave us a F-statistic", 
print "of %.2f or more" % observed_f_statistic, "is", (count / float(num_shuffles))
