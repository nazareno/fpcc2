#!/usr/bin/python

######################################
# Regression Significance Test 
# From: Statistics is Easy! By Dennis Shasha and Manda Wilson
# 
# Assuming that x is not a good predictor of y, tests to see the probability
# of getting a slope by chance alone greater than or equal to the one observed (less
# than or equal to the one observed if the observed slope is negative).
# Uses shuffling to get a distribution to compare 
# to the observed f-statistic.
# 
# Author: Manda Wilson
#
# Example of FASTA formatted input file: 
# >grp_x
# 1350 1510 1420 1210 1250 1300 1580 1310 1290 1320 1490 1200 1360
# >grp_y
# 3.6 3.8 3.7 3.3 3.9 3.4 3.8 3.7 3.5 3.4 3.8 3.0 3.1
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
# 2. Set a counter to 0, this will count the number of times we get a slope (b) 
#    greater than or equal to 0.0014.  Note: if you have a negative slope, count
#    the number of times you get a slope less than or equal to the original negative slope.  
#
# 3. Do the following 10,000 times:
#    a. Shuffle the y values. 
#    b. Calculate regression line on the results from step (3a), just as we did in step (1)
#    c. If the slope (b)  from step (3b) is greater than or equal to our observed slope (0.0014), 
#       increment our counter from step (2). 
#
# 3. counter / 10,000 equals the probability of getting a slope greater than
#    or equal to 0.0014, assuming x does not predict y
#
######################################

import random

######################################
#
# Adjustable variables
#
######################################

input_file = 'Correlation.vals'

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

count = 0
num_shuffles = 10000

for i in range(num_shuffles):
        new_y_values = shuffle([grp_y])[0]
        (a, b) = regressionline(grp_x, new_y_values)
        if ((observed_b >= 0 and b > observed_b) or (observed_b < 0 and b < observed_b)):
                count = count + 1

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
print count, "out of", num_shuffles, "experiments had a slope",
if observed_b < 0:
        print "less than or equal to",
else:
        print "greater than or equal to",
print "%.4f" % observed_b, "."
print "The chance of getting a slope",
if observed_b < 0:
        print "less than or equal to",
else:
        print "greater than or equal to",
print "%.4f" % observed_b, "is", (count / float(num_shuffles)), "."

