import math

# === Adjustable Variables ===
num_targets = 5                             # Number of targets in multiplex
detection_efficiency_percent = 50        # e.g., 50 means only half the copies are effectively detected
haploid_genome_fg = 150.0 / 977.0 * 1e3     # A. thaliana haploid genome in fg
num_droplets = 15_000                       # Typical droplet count

# === 1) Compute the Optimal λ per Target ===
# That maximizes single-positive fraction
optimal_lambda = -math.log((num_targets - 1) / num_targets)

# === 2) Theoretical Input for 100% Detection ===
# Haploid genome copies needed = optimal_lambda * num_droplets
genomes_needed = optimal_lambda * num_droplets
theoretical_input_ng = (genomes_needed * haploid_genome_fg) / 1e6

# === 3) Adjust for Real-World Detection Efficiency
detection_efficiency = detection_efficiency_percent / 100.0
adjusted_input_ng = theoretical_input_ng / detection_efficiency

# === 4) Estimate Droplet Categories ===
# Using Poisson logic for multiplex
lambda_per_target = optimal_lambda

# Probability of a droplet being negative for all targets
p_negative = math.exp(-num_targets * lambda_per_target)
negative_droplets = num_droplets * p_negative

# Probability of being positive for exactly one target
p_single_positive = num_targets * math.exp(-(num_targets - 1) * lambda_per_target) * (1 - math.exp(-lambda_per_target))
single_positive_droplets = num_droplets * p_single_positive
single_positive_per_target = single_positive_droplets / num_targets


# Everything else is multi-positive
multiple_positive_droplets = num_droplets - single_positive_droplets - negative_droplets

# === Output ===
print("\n=== ddPCR Multiplex Calculations ===")
print("\nOptimal Conditions:")
print(f" - Number of targets:                   {num_targets}")
print(f" - Optimal λ per target:                {optimal_lambda:.2f}")
print(f" - Genome copies needed:                {genomes_needed:.0f}")
print("\nInput adjustment:")
print(f" - Theoretical input:                   {theoretical_input_ng:.2f} ng")
print(f" - Detection efficiency:                {detection_efficiency_percent:.1f}%")
print(f" - Adjusted input:                      {adjusted_input_ng:.2f} ng")
print("\nDroplet distribution:")
print(f" - Singe-positive droplets (per target):{single_positive_per_target:.0f}")
print(f" - Target-negative droplets:            {negative_droplets:.0f}")
print(f" - Multi-positive droplets:             {multiple_positive_droplets:.0f}")
print ("")