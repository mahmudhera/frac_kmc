import argparse
import subprocess
import time

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filename", help="name of file")
    args = parser.parse_args()
    filename = args.filename

    # run sourmash and frackmc for different kmer sizes and scaled values
    records_file = open("records", "w")
    records_file.write("kmer_size,scaled,sourmash_time,frackmc_time\n")

    scaled_values = [1, 10, 100, 1000, 10000, 100000]
    kmer_sizes = [21, 31, 41, 51, 61, 71]
    seed = 42
    for scaled in scaled_values:
        for kmer_size in kmer_sizes:
            print(f"Running k={kmer_size}, scaled={scaled}")
            sourmash_sketch_name = f"{filename}_sketch_sm_k_{kmer_size}_s_{scaled}.sig"
            frackmc_sketch_name =  f"{filename}_sketch_fk_k_{kmer_size}_s_{scaled}.sig"
            
            # run sourmash
            start_time = time.time()
            cmd = f"sourmash sketch dna -p k={kmer_size},scaled={scaled} -o {sourmash_sketch_name} {filename}"
            subprocess.run(cmd.split(' '))
            end_time = time.time()
            sourmash_time = end_time - start_time
            
            # run frackmc
            start_time = time.time()
            cmd = f"fracKmcSketch {filename} {frackmc_sketch_name} --ksize {kmer_size} --scaled {scaled} --seed {seed}"
            subprocess.run(cmd.split(' '))
            end_time = time.time()
            frackmc_time = end_time - start_time
            
            # record results
            records_file.write(f"{kmer_size},{scaled},{sourmash_time},{frackmc_time}\n")
            
    records_file.close()

if __name__ == "__main__":
    main()
