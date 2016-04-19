#!/usr/bin/python

######################################
# Regression Confidence Interval Example 
# From: Statistics is Easy! By Dennis Shasha and Manda Wilson
# 
# Uses shuffling & bootstrapping to get a 90% confidence interval for the slope of the regression line.
# 
# Author: Manda Wilson
#
# Example of FASTA formatted input file: 
# >grp_x
# 1350 1510 1420 1210 1250 1300 1580 1310 1290 1320 1490 1200 1360
# >grp_y
# 3.6 3.8 3.7 3.3 3.9 3.4 3.8 3.7 3.5 3.4 3.8 3.0 3.1
# 
# Included in the code, but NOT in the pseudocode,
# is the bias-corrected confidence interval.  
# See the chapter titled "Bias Corrected Confidence Intervals" in Statistics is Easy!
# for more information on bias-corrected cofidence intervals.
#
# Pseudocode:
#
# 1. Calculate the regression line for the observed values (in our example it is y' = 0.0014 x +  1.6333).  
#    a. mean_x = sum of x values / number of x values
#    b. mean_y = sum of y values / number of y values
#    c. Calculate the sum of products:
#	i. Initialize total to 0
#	ii. For each pair of x, y values:
#           I. product = (x - mean_x) * (y - mean_y) 
#           II. total += product
#    d. Calculate the sum of squares for the x values: 
#       i. Initialize total to 0 
#       ii. For each x value:
#           I. diff_sq = (x - mean_x)^2
#           II. total += diff_sq
#    e. b = sum of products / sum of squares for the x values 
#    f. a = mean_y - (b * mean_x)
#    g. Line of best fit is: y' = bx + a
#
# 2. Do the following 10,000 times:
#    a. 
#       i. Create two new arrays of the same size as the original sample arrays (one for x and one for y)
#       ii. Fill the new arrays with randomly picked pairs from the original arrays (randomly picked with replacement)
#    b. Compute the regression line on the bootstrap samples, just as we did in step (1)
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

######################################
#
# Adjustable variables
#
######################################

input_file = 'Correlation.vals'
conf_interval = 0.9

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

# x, y are arrays of sample values
# returns two new arrays with randomly picked
# (with replacement) pairs from x and y
def bootstrap(x, y):
        samp_x = []
	samp_y = []
	num_samples = len(x)
        for i in range(num_samples):
		index = random.randint(0, num_samples - 1)
	        samp_x.append(x[index])
        	samp_y.append(y[index])
	return (samp_x, samp_y)

def sumofsq(vals, mean):
	count = len(vals)
	total = 0
	for i in range (count):
		diff_sq = (vals[i] - mean)**2
		total += diff_sq
	return total

def sumofproducts(x_vals, y_vals,  mean_x, mean_y):
	count = len(x_vals)
	total = 0
	for i in range (count):
		product = (x_vals[i] - mean_x) * (y_vals[i] - mean_y)
		total += product
	return total

def regressionline(grp_x, grp_y):
	sum_x = sum(grp_x)
	sum_y = sum(grp_y)

	count_x = float(len(grp_x))
	count_y = float(len(grp_y))

	mean_x = sum_x / count_x
	mean_y = sum_y / count_y

	# get the sum of products
	sum_of_prod = sumofproducts(grp_x, grp_y, mean_x, mean_y)

	# get the sum of squares for x
	sum_of_sq_x = sumofsq(grp_x, mean_x)

	b = sum_of_prod / sum_of_sq_x
	a = mean_y - (b * mean_x)
	return (a, b)

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

grp_x = samples[0]
grp_y = samples[1]

(observed_a, observed_b) = regressionline(grp_x, grp_y)

num_resamples = 10000   # number of times we will resample from our original samples
num_below_observed = 0   # count the number of bootstrap values below the observed sample statistic
out = []                # will store results of each time we resample

for i in range(num_resamples):
	# bootstrap - then compute our statistic of interest
	# keep pairs together
	# append statistic to out
	(boot_x, boot_y) = bootstrap(grp_x, grp_y)
	# now we have a list of bootstrap samples, run regressionline
	(boot_a, boot_b) = regressionline(boot_x, boot_y)
	if boot_b < observed_b:
		num_below_observed += 1
	out.append(boot_b)

out.sort()

# standard confidence interval computations
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

sign = "+"

if (observed_a < 0):
        sign = "-"

observed_a = abs(observed_a)

print "Line of best fit for observed data: ",
print "y' = %.4f" % observed_b, "x", sign, " %.4f" % observed_a
print "We have", conf_interval * 100, "% confidence that the true slope",
print "is between: %.4f" % out[lower_bound], "and %.4f" % out[upper_bound]

print "***** Bias Corrected Confidence Interval *****"
print "We have", conf_interval * 100, "% confidence that the true slope",
print "is between: %.4f" % out[bias_corr_lower_bound], "and %.4f" % out[bias_corr_upper_bound]
