import math

# Frequencies:
p_c1_loss  = 0.03  # category 1
p_c1_gain  = 0.02  # category 2
p_c4_loss  = 0.06  # category 3
p_c4_gain  = 0.04  # category 4
aneuploid_fraction = (p_c1_loss + p_c1_gain + p_c4_loss + p_c4_gain)
p_normal = 1.0 - aneuploid_fraction
p = [p_normal, p_c1_loss, p_c1_gain, p_c4_loss, p_c4_gain]

# For weighting detection among all aneuploids:
# => among aneuploids, 3/15 = 0.2 are chr1-loss, 2/15 = 0.1333 are chr1-gain, 6/15=0.4 are chr4-loss, 4/15=0.2667 are chr4-gain
frac_c1_loss = p_c1_loss/aneuploid_fraction
frac_c1_gain = p_c1_gain/aneuploid_fraction
frac_c4_loss = p_c4_loss/aneuploid_fraction
frac_c4_gain = p_c4_gain/aneuploid_fraction

def factorial(n):
    return math.factorial(n)

def multinomial_coefficient(counts):
    """
    counts = (x0, x1, x2, x3, x4)
    returns n! / (x0! * x1! * x2! * x3! * x4!)
    """
    n = sum(counts)
    num = factorial(n)
    for c in counts:
        num //= factorial(c)
    return num

def pool_negative_prob(n):
    """
    Probability a pool of size n is 'negative':
      => #chr1-loss == #chr1-gain  AND  #chr4-loss == #chr4-gain
    """
    neg_sum = 0.0
    # x0= #normal, x1= #c1Loss, x2= #c1Gain, x3= #c4Loss, x4= #c4Gain
    for x1 in range(n+1):
        for x2 in range(n+1 - x1):
            if x2 != x1:
                continue
            for x3 in range(n+1 - x1 - x2):
                for x4 in range(n+1 - x1 - x2 - x3):
                    if x4 != x3:
                        continue
                    x0 = n - (x1 + x2 + x3 + x4)
                    # Probability of exactly (x0,x1,x2,x3,x4):
                    coeff = multinomial_coefficient([x0, x1, x2, x3, x4])
                    prob = (coeff *
                            (p[0]**x0) *
                            (p[1]**x1) *
                            (p[2]**x2) *
                            (p[3]**x3) *
                            (p[4]**x4))
                    neg_sum += prob
    return neg_sum

def detection_prob_for_aneuploid(n, which_aneu):
    """
    Probability that an individual with 'which_aneu' (1..4) is detected
    in a 2-stage scheme with pool size n (given the occupant is of that type).
    We'll enumerate all ways the other n-1 can be distributed,
    and check whether the pool ends up flagged 'positive' or 'negative'.
    """
    missed_sum = 0.0  # sum of probabilities for "pool is negative" with occupant
    # We'll do a nested loop for the other (n-1) individuals:
    for x0_ in range(n):
        for x1_ in range(n - x0_):
            for x2_ in range(n - x0_ - x1_):
                for x3_ in range(n - x0_ - x1_ - x2_):
                    x4_ = (n - 1) - (x0_ + x1_ + x2_ + x3_)
                    if x4_ < 0:
                        continue
                    
                    # Probability of exactly x0_,x1_,x2_,x3_,x4_ among the OTHER n-1
                    # (each category has prob p[k])
                    # => Multinomial over n-1
                    cnts_other = (x0_, x1_, x2_, x3_, x4_)
                    sum_cnts = sum(cnts_other)
                    if sum_cnts != (n - 1):
                        continue
                    # multinomial probability for these n-1 draws:
                    coeff = multinomial_coefficient(cnts_other)
                    prob_others = coeff
                    for k, cval in enumerate(cnts_other):
                        prob_others *= (p[k]**cval)
                    
                    if prob_others == 0:
                        continue
                    
                    # Now add the occupant's category
                    x0 = x0_
                    x1 = x1_
                    x2 = x2_
                    x3 = x3_
                    x4 = x4_
                    
                    if which_aneu == 1: 
                        x1 += 1
                    elif which_aneu == 2:
                        x2 += 1
                    elif which_aneu == 3:
                        x3 += 1
                    elif which_aneu == 4:
                        x4 += 1
                    
                    # Check if final pool is negative:
                    # negative => #c1Loss == #c1Gain AND #c4Loss == #c4Gain
                    if (x1 == x2) and (x3 == x4):
                        missed_sum += prob_others
    # Probability occupant is missed:
    p_missed = missed_sum
    # => detection probability = 1 - p_missed
    return 1.0 - p_missed

# We only need detection rates, so let's compute and print for n=2,3,5,10.
pool_sizes = [2, 3, 5, 10]
print("POOLING CALCULATIONS FOR ANEUPLOID DETECTION")
print("---------------------------------------------")
print(f"{'Pool Size':>8} | {'Prob(negative)':>13} | {'d(c1Loss)':>9} | {'d(c1Gain)':>9} | {'d(c4Loss)':>9} | {'d(c4Gain)':>9} | {'WeightedDet':>11}")

for n in pool_sizes:
    pNeg = pool_negative_prob(n)
    pPos = 1 - pNeg
    
    d1L = detection_prob_for_aneuploid(n, 1)  # chr1 loss
    d1G = detection_prob_for_aneuploid(n, 2)  # chr1 gain
    d4L = detection_prob_for_aneuploid(n, 3)  # chr4 loss
    d4G = detection_prob_for_aneuploid(n, 4)  # chr4 gain
    
    # Weighted average detection among aneuploids:
    #   among aneuploids: 3/15 = c1Loss, 2/15 = c1Gain, 6/15 = c4Loss, 4/15 = c4Gain
    frac_c1L = 3/15
    frac_c1G = 2/15
    frac_c4L = 6/15
    frac_c4G = 4/15
    
    overall_det = (frac_c1L*d1L
                   + frac_c1G*d1G
                   + frac_c4L*d4L
                   + frac_c4G*d4G)
    
    print(f"{n:8d} | {pNeg:13.4f} | {d1L:9.3f} | {d1G:9.3f} | {d4L:9.3f} | {d4G:9.3f} | {overall_det:11.3f}")