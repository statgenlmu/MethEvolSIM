import msprime
import argparse
import os
import pandas as pd
import shutil

parser = argparse.ArgumentParser(description='Simulate ARGs with msprime')


parser.add_argument('--samples', type=int, help='number of samples')
parser.add_argument('--ploidy', type=int, help='value of ploidy')
parser.add_argument('--populationsize', type=int, help='size of the population')
parser.add_argument('--sequencelength', type=int, help='length of the sequence')
parser.add_argument('--output', type=str, help='output location')
parser.add_argument('--nsim', type=int, help='number of simulations')
parser.add_argument('--recrates', type=str, help='Path to recombination rates CSV')

def recreate_directory(directory_path):
    if os.path.exists(directory_path):
        shutil.rmtree(directory_path)
    os.makedirs(directory_path)



args = parser.parse_args()

samples = args.samples
ploidy = args.ploidy
population_size = args.populationsize
sequence_length = args.sequencelength
output_dir = args.output
nsim = args.nsim
recrates = args.recrates

recreate_directory(os.path.join(output_dir, "nodes"))
recreate_directory(os.path.join(output_dir, "edges"))
nodes_dir = os.path.join(output_dir, "nodes")
edges_dir = os.path.join(output_dir, "edges")
os.makedirs(nodes_dir, exist_ok=True)
os.makedirs(edges_dir, exist_ok=True)

rates_df = pd.read_csv(recrates)

num_digits = len(str(args.nsim))+1
for sim_id, rec_rate in zip(rates_df["sim_id"], rates_df["recrate"]):
    nodes_file_name = f"nodes_{str(sim_id).zfill(num_digits)}.csv"
    edges_file_name = f"edges_{str(sim_id).zfill(num_digits)}.csv"
    


    ts = msprime.sim_ancestry(
        samples=args.samples,
        ploidy=args.ploidy,
        population_size=args.populationsize,
        recombination_rate=rec_rate,
        sequence_length=args.sequencelength,
        #random_seed=sim_id,
        additional_nodes=msprime.NodeType.RECOMBINANT,
        coalescing_segments_only=False
        #end_time=160 #this always simulates two generations
    )
    tables = ts.tables

    nodes_df = pd.DataFrame({
        "id": range(tables.nodes.num_rows),
        "time": tables.nodes.time,
        "flags": tables.nodes.flags,
        "population": tables.nodes.population,
        "individual": tables.nodes.individual,
    })
    edges_df = pd.DataFrame({
        "id": range(tables.edges.num_rows),
        "left": tables.edges.left,
        "right": tables.edges.right,
        "parent": tables.edges.parent,
        "child": tables.edges.child,
    })

    nodes_df.to_csv(os.path.join(nodes_dir, nodes_file_name), index=False)
    edges_df.to_csv(os.path.join(edges_dir, edges_file_name), index=False)


