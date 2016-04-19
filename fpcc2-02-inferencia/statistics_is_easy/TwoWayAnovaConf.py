#!/usr/bin/python

######################################
# Two-Way Confidence Interval
# From: Statistics is Easy! By Dennis Shasha and Manda Wilson
# 
# Uses shuffling & bootstrapping to get a 90% confidence interval for the f-statistic.
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
# Included in the code, but NOT in the pseudocode,
# is the bias-corrected confidence interval.  
# See the chapter titled "Bias Corrected Confidence Intervals" in Statistics is Easy!
# for more information on bias-corrected cofidence intervals.
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
# 2. Do the following 10,000 times:
#    a. For each sample we have get a bootstrap sample:
#       i. Create a new array of the same size as the original sample
#       ii. Fill the array with randomly picked values from the original sample (randomly picked with replacement)
#    b. Calculate the f-statistic on the bootstrap samples, just as we did in step (1)
#       with the original samples.
#
# 3. Sort the f-statistics computed in step (2).
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

input_file = 'TwoWayAnova.vals'

######################################
#
# Subroutines
#
######################################

# maps proportion of values above statistic
# to number of standard deviations above statistic
# keys will be index / 100 \:[0-9]\.[0-9][0-9]\,
area_to_sd_map = [0.0000, 0.0040, 0.0080, 0.0120, 0.0160, 0.0199, 0.0239, 0.0279, 0.0319, 0.0359, 0.0398, 0.0438, 0.0478, 0.0517, 0.0557, 0.0596, 0.0636, 0.0675, 0.0714, 0.0753, 0.0793, 0.0832, 0.0871, 0.0910, 0.0948, 0.0987, 0.1026, 0.1064, 0.1103, 0.1141, 0.1179, 0.1217, 0.1255, 0.1293, 0.1331, 0.1368, 0.1406, 0.1443, 0.1480, 0.1517, 0.1554, 0.1591, 0.1628, 0.1664, 0.1700, 0.1736, 0.1772, 0.1808, 0.1844, 0.1879, 0.1915, 0.1950, 0.1985, 0.2019, 0.2054, 0.2088, 0.2123, 0.2157, 0.2190, 0.2224, 0.2257, 0.2291, 0.2324, 0.2357, 0.2389, 0.2422, 0.2454, 0.2486, 0.2517, 0.2549, 0.2580, 0.2611, 0.2642, 0.2673, 0.2704, 0.2734, 0.2764, 0.2794, 0.2823, 0.2852, 0.2881, 0.2910, 0.2939, 0.2967, 0.2995, 0.3023, 0.3051, 0.3078, 0.3106, 0.3133, 0.3159, 0.3186, 0.3212, 0.3238, 0.3264, 0.3289, 0.3315, 0.3340, 0.3365, 0.3389, 0.3413, 0.3438, 0.3461, 0.3485, 0.3508, 0.3531, 0.3554, 0.3577, 0.3599, 0.3621, 0.3643, 0.3665, 0.3686, 0.3708, 0.3729, 0.3749, 0.3770, 0.3790, 0.3810, 0.3830, 0.3849, 0.3869, 0.3888, 0.3907, 0.3925, 0.3944, 0.3962, 0.3980, 0.3997, 0.4015, 0.4032, 0.4049, 0.4066, 0.4082, 0.4099, 0.4115, 0.4131, 0.4147, 0.4162, 0.4177, 0.4192, 0.4207, 0.4222, 0.4236, 0.4251, 0.4265, 0.4279, 0.4292, 0.4306, 0.4319, 0.4332, 0.4345, 0.4357, 0.4370, 0.4382, 0.4394, 0.4406, 0.4418, 0.4429, 0.4441, 0.4452, 0.4463, 0.4474, 0.4484, 0.4495, 0.4505, 0.4515, 0.4525, 0.4535, 0.4545, 0.4554, 0.4564, 0.4573, 0.4582, 0.4591, 0.4599, 0.4608, 0.4616, 0.4625, 0.4633, 0.4641, 0.4649, 0.4656, 0.4664, 0.4671, 0.4678, 0.4686, 0.4693, 0.4699, 0.4706, 0.4713, 0.4719, 0.4726, 0.4732, 0.4738, 0.4744, 0.4750, 0.4756, 0.4761, 0.4767, 0.4772, 0.4778, 0.4783, 0.4788, 0.4793, 0.4798, 0.4803, 0.4808, 0.4812, 0.4817, 0.4821, 0.4826, 0.4830, 0.4834, 0.4838, 0.4842, 0.4846, 0.4850, 0.4854, 0.4857, 0.4861, 0.4864, 0.4868, 0.4871, 0.4875, 0.4878, 0.4881, 0.4884, 0.4887, 0.4890, 0.4893, 0.4896, 0.4898, 0.4901, 0.4904, 0.4906, 0.4909, 0.4911, 0.4913, 0.4916, 0.4918, 0.4920, 0.4922, 0.4925, 0.4927, 0.4929, 0.4931, 0.4932, 0.4934, 0.4936, 0.4938, 0.4940, 0.4941, 0.4943, 0.4945, 0.4946, 0.4948, 0.4949, 0.4951, 0.4952, 0.4953, 0.4955, 0.4956, 0.4957, 0.4959, 0.4960, 0.4961, 0.4962, 0.4963, 0.4964, 0.4965, 0.4966, 0.4967, 0.4968, 0.4969, 0.4970, 0.4971, 0.4972, 0.4973, 0.4974, 0.4974, 0.4975, 0.4976, 0.4977, 0.4977, 0.4978, 0.4979, 0.4979, 0.4980, 0.4981, 0.4981, 0.4982, 0.4982, 0.4983, 0.4984, 0.4984, 0.4985, 0.4985, 0.4986, 0.4986, 0.4987, 0.4987, 0.4987, 0.4988, 0.4988, 0.4989, 0.4989, 0.4989, 0.4990, 0.4990]

def sd_to_area(sd):
	sign = 1
	if sd < 0:
		sign = -1
	sd = math.fabs(sd)  # get the absolute value of sd
	index = int(sd * 100)
	if len(area_to_sd_map) <= index:
		return sign * area_to_sd_map[-1] # return last element in array
	if index == (sd * 100):
		return sign * area_to_sd_map[index]
	return sign * (area_to_sd_map[index] + area_to_sd_map[index + 1]) / 2

def area_to_sd(area):
	sign = 1
	if area < 0:
		sign = -1
	area = math.fabs(area)
	for a in range(len(area_to_sd_map)):
		if area == area_to_sd_map[a]:
			return sign * a / 100
		if 0 < a and area_to_sd_map[a - 1] < area and area < area_to_sd_map[a]:
			# our area is between this value and the previous
			# for simplicity, we will just take the sd half way between a - 1 and a
			return sign * (a - .5) / 100
	return sign * (len(area_to_sd_map) - 1) / 100

def bootstrap(x):
        samp_x = []
        for i in range(len(x)):
                samp_x.append(random.choice(x))
        return samp_x

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

num_resamples = 10000    # number of times we will resample from our original samples
num_below_observed = 0   # count the number of bootstrap values below the observed sample statistic
out = []				# will store results of each time we resample

for i in range(num_resamples):
	# get bootstrap samples for each of our groups
	# then compute our statistic of interest
	# append statistic to out
	bootstrap_samples = []  # list of lists
	r = 0
	c = 0
	for r in range(len(samples)):
		bootstrap_samples.append([])
		for c in range(len(samples[r])):
			bootstrap_samples[r].append(bootstrap(samples[r][c]))
	# now we have a list of new samples, run onewayanova
	boot_f_statistic = twowayanova(bootstrap_samples)
	if boot_f_statistic < observed_f_statistic:
		num_below_observed += 1
	out.append(boot_f_statistic)

out.sort()

# standard confidence interval computations
conf_interval = 0.9
tails = (1 - conf_interval) / 2

# in case our lower and upper bounds are not integers,
# we decrease the range (the values we include in our interval),
# so that we can keep the same level of confidence
lower_bound = int(math.ceil(num_resamples * tails))
upper_bound = int(math.floor(num_resamples * (1 - tails)))

# bias-corrected confidence interval computations
p = num_below_observed / float(num_resamples)	# proportion of bootstrap values below the observed value

dist_from_center = p - .5	# if this is negative, the original is below the center, if positive, it is above
z_0 = area_to_sd(dist_from_center)

# now we want to find the proportion that should be between the mean and one of the tails
tail_sds = area_to_sd(conf_interval / 2) 
z_alpha_over_2 = 0 - tail_sds
z_1_minus_alpha_over_2 = tail_sds

# in case our lower and upper bounds are not integers,
# we decrease the range (the values we include in our interval),
# so that we can keep the same level of confidence
bias_corr_lower_bound = int(math.ceil(num_resamples * (0.5 + sd_to_area(z_alpha_over_2 + (2 * z_0)))))
bias_corr_upper_bound =  int(math.floor(num_resamples * (0.5 + sd_to_area(z_1_minus_alpha_over_2 + (2 * z_0)))))

######################################
#
# Output
#
######################################

# print observed value and then confidence interval
print "Observed F-statistic: %.2f" % observed_f_statistic
print "We have", conf_interval * 100, "% confidence that the true F-statistic",
print "is between: %.2f" % out[int(lower_bound)], "and %.2f" % out[int(upper_bound)]

print "***** Bias Corrected Confidence Interval *****"
print "We have", conf_interval * 100, "% confidence that the true F-statistic",
print "is between: %.2f" % out[bias_corr_lower_bound], "and %.2f" % out[bias_corr_upper_bound]
