#!/usr/bin/python

######################################
# Two-Way ANOVA Significance Test 
# From: Statistics is Easy! By Dennis Shasha and Manda Wilson
# 
# Assuming that there is no interaction between the three drugs used and sex, tests to see the probability
# of getting a f-statistic by chance alone greater than or equal to the one observed.
# Uses shuffling to get a distribution to compare 
# to the observed f-statistic.
# 
# Author: Manda Wilson
#
# Example of FASTA formatted input file: 
# >grp 1 1
# 15 12 13 16
# >grp 2 1
# 19 17 16 15
# >grp 3 1
# 14 13 12 17
# >grp 1 2
# 13 13 12 11
# >grp 2 2
# 13 11 11 17
# >grp 3 2
# 11 12 10
#
# Note that the input file is a little different from our other input files.  
# On each descriptive line (beginning with '>') there is the group name (which is ignored)
# followed by a whitespace character, followed by a number, then another whitespace character,
# then another number.  These numbers represent indexes in a two-dimensional matrix.  So ">grp 1 1"
# represents the element (a list) at row 1, column 1, ">grp 2 1" represents the element (a list) 
# at row 2, column 1, etc.  This is how we seperate the groups by two factors.  The rows are 
# different categories within one factor, the columns different categories within another factor.
# 
# Pseudocode:
#
# 1. Calculate f-statistic for the observed values (in our example it is 0.93).  
#    a. Initialize within sum of squares (wss), total sum (ts), and
#       total count (tc) to 0
#    b. For each row: (here we loop through the categories of factor a)
#       For each column: (here we loop through the categories of factor b)
#       i. Within group mean (wgm) = sum of group values / number of values in group
#	   ii. Store the number of values in group
#	   iii. Sum the following: for each value in the group
#            I. Subtract the value from wgm
#            II. Square the result of step (1biiiI)
#            III.  Add the result of step (1biiiII) to wss
#    c. Total mean (tm) = ts / tc
#    d. Mean group means (mgm) = sum(wgm) / num_groups
#    e. Total sum of squares (tss) = 
#	   i. Sum the following: for each value (include all values)
#          I. Subtract the value from tm
#          II. Square the result of step (1eiI)
#		 III.  Add the result of step (1eiII) to tss
#    d. Factor a sum of squares (fass) = 
#       i. Sum the following: for each category in factor a
#          I. Subtract the sum of values in this category from tm
#          II. Square the result of step (1diI)
#          III. Multiply the result of step (1diII) by 
#               the number of values in that category
#    e. Factor b sum of squares (fbss) = 
#       i. Sum the following: for each category in factor b
#          I. Subtract the sum of values in this category from tm
#          II. Square the result of step (1eiI)
#          III. Multiply the result of step (1eiII) by 
#               the number of values in that category
#    f. Between group sum of squares (bss) = tss - wss
#    g. Factor sum of squares (fss) = bss - (fass + fbss)
#    h. Between degrees of freedom (bdf) = number of groups - 1
#    i. Within degrees of freedom (wdf) = tc - number of groups
#    j. Factor a degrees of freedom (fadf) = number of categories in factor_a - 1
#    k. Factor b degrees of freedom (fbdf) = number of categories in factor_b - 1
#    l. Interaction degrees of freedom (idf) = bdf - fadf - fbdf
#    m. Within group variance (wgv) = wss / wdf
#    n. Between group variance (bgv) = bss / bdf
#    o. Interaction variance (iv) = fss / idf
#    i. f-statistic = iv / wgv
#
# 2. Set a counter to 0, this will count the number of times we get a f-statistic 
#    greater than or equal to our observed f-statistic (in our example 0.93).  
#
# 3. Do the following 10,000 times:
#    a. Shuffle the observed values. To do this:
#       i. Put the values from all the groups into one array
#       ii. Shuffle the pooled values
#       iii. Reassign the pooled values to groups of the same size as the original groups
#    b. Calculate f-statistic on the results from step (3a), just as we did in step (1)
#    c. If the result from step (3b) is greater than or equal to our observed f-statistic (0.93), 
#       increment our counter from step (2). 
#
# 3. counter / 10,000 equals the probability of getting a f-statistic greater than or equal to
#    our observed f-stat (0.93), assuming there is no difference between the groups
#
######################################

import random

######################################
#
# Adjustable variables
#
######################################

input_file = 'TwoWayAnova.vals'

######################################
#
# Subroutines
#
######################################

# a list of lists
def shuffle(grps):
	num_rows = len(grps)
	num_cols = 0
	pool = []
	
	# throw all values together
	for r in range(num_rows):
		num_cols = len(grps[r])
		for c in range(num_cols):
			pool.extend(grps[r][c])
	# mix them up
	random.shuffle(pool)
	# reassign to groups
	new_grps = []
	start_index = 0
	end_index = 0
	for r in range(num_rows):
		new_grps.append([]);
		for c in range(num_cols):
			new_grps[r].append([])	
			end_index = start_index + len(grps[r][c])
			new_grps[r][c] = pool[start_index:end_index]
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
def twowayanova(grps):
	num_rows = len(grps)	# equals the number of categories in factor a
	num_cols = len(grps[0])  # equals the number of categories in factor b
	num_grps = num_rows * num_cols
	within_group_means = []
	grp_counts = []
	within_ss = 0
	total_sum = 0
	total_count = 0
	factor_a = []  # will have num_rows categories, each cateory has num_cols lists added to it
	factor_b = []  # will have num_cols categories, each cateory has num_rows lists added to it
	all_vals = []
	
	for r in range (num_rows):		# looping through factor a categories
		if len(factor_a) <= r:		# if this category was not initialized yet, initialize it
			factor_a.append([])
		for c in range (num_cols):    # looping through factor b categories
			if len(factor_b) <= c:   # if this category was not initialized yet, initialize it
				factor_b.append([])
			grp = grps[r][c]
			all_vals += grp
			factor_a[r] += grp
			factor_b[c] += grp
			grp_count = float(len(grp))
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

	total_mean = total_sum / float(total_count) 
	mean_grp_means = sum(within_group_means) / float(num_grps)
	
	# total sum of squares 
	total_ss = sumofsq(all_vals, total_mean)

	factor_a_ss = weightedsumofsq(map(lambda x: sum(x) / float(len(x)), factor_a), map(lambda x: float(len(x)), factor_a), total_mean)
	factor_b_ss = weightedsumofsq(map(lambda x: sum(x) / float(len(x)), factor_b), map(lambda x: float(len(x)), factor_b), total_mean)
	
	# get between group sum of squares
	# NOTE: to compute the between group sum of squares:
	# sum the following: for every element in the overall group
	# subtract that element's group mean from the overall group mean
	# square the difference
	between_ss = total_ss - within_ss
	# NOTE: we could have done the following
	# between_ss = weightedsumofsq(within_group_means, grp_counts, total_mean)
	
	factor_ss = between_ss - (factor_a_ss + factor_b_ss)

	# now we want to find out how different the groups are between each other
	# compared to how much the values vary within the groups
	# if all groups vary a lot within themselves, and there is no significant difference
	# between the groups, then we expect the differences between the groups to vary by about the same amount
	# so lets get the ratio of the between group variance and the within group variance
	# if the ratio is 1, then there is no difference between the groups
	# if it is significantly larger than one, then there is a significant difference between the groups
	# remember: even if the groups are significantly different, we still won't know which groups are different

	# the between group degrees of freedom
	# is equal to the number of groups - 1
	# this is because once we know the number of groups - 1, we know the last group
	between_df = num_grps - 1

	# the within group degress of freedom
	# is equal to the total number of values minus the number of groups
	# this is because for each group, once we know the count - 1 values, we know the last value for that group
	# so we lose the number of groups * 1 degrees of freedom
	within_df = total_count - num_grps

	# the interaction degrees of freedom
	# is equal to between group degrees of freedom
	# minus factor a degrees of freedom
	# minus factor b degrees of freedom
	factor_a_df = num_rows - 1    # number of categories in factor_a - 1
	factor_b_df = num_cols - 1    # number of categories in factor_a - 1
	interaction_df = between_df - factor_a_df - factor_b_df

	within_var = within_ss / within_df
	between_var = between_ss / between_df
	interaction_var = factor_ss / interaction_df

	# return our f-statistic
	return interaction_var / within_var

######################################
#
# Computations
#
######################################

# list of lists
samples = [] 

# file must be in FASTA format
infile=open(input_file)
r = 0
c = 0
for line in infile:
	if line.startswith('>'):
		# read in index
		tokens = line.split()    # split on whitespace
		c = int(tokens[-1]) - 1  # column index is last token on line, we subtract one b/c python indexes start at 0
		r = int(tokens[-2]) - 1  # row index is second to last token on line, we subtract one b/c python indexes start at 0
		# make sure we have enough rows in our samples matrix
		while r >= len(samples):
			samples.append([])
		# now make sure that row has enough columns
		while c >= len(samples[r]):
			samples[r].append([])
	elif not line.isspace():
		# line must contain values for previous sample
		samples[r][c] += map(float,line.split())
infile.close()

observed_f_statistic = twowayanova(samples)

count = 0
num_shuffles = 10000

for i in range(num_shuffles):
	new_samples = shuffle(samples)
	f_statistic = twowayanova(new_samples)
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
