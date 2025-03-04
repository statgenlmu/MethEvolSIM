import argparse
import random
import csv
import os
import shutil


parser = argparse.ArgumentParser(description='Sample recombination rates for simulations')
parser.add_argument('--nsim', type=int, required=True, help='Number of simulations')
parser.add_argument('--seed', type=int, default=123456, help='Random seed')
parser.add_argument('--output', type=str, default='recombination_rates.csv', help='Output CSV file')
parser.add_argument('--lowerbound', type=float, required=True, help='Lower bound of uniform distribution')
parser.add_argument('--upperbound', type=float, required=True, help='Upper bound of uniform distribution')
parser.add_argument('--outdir', type=str, required=True, help='Output directory')

args = parser.parse_args()

def recreate_directory(directory_path):
    if os.path.exists(directory_path):
        shutil.rmtree(directory_path)
    os.makedirs(directory_path)


# Set the random seed for reproducibility
#random.seed(args.seed)

# Sample recombination rates
sim_ids = [str(i + 1).zfill(len(str(args.nsim))) for i in range(args.nsim)]
recombination_rates = []

recreate_directory(args.outdir)

for sim_id in sim_ids:
    recrate = random.uniform(args.lowerbound, args.upperbound)
    recombination_rates.append({'sim_id': sim_id, 'recrate': recrate})

    # Save to CSV
with open(args.output, 'w', newline='') as csvfile:
    fieldnames = ['sim_id', 'recrate']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(recombination_rates)












