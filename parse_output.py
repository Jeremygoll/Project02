# a python function that opens a folder of files and puts them into a csv

import os
import csv
import sys

def main():
    # get the folder name
    folder = sys.argv[1]
    # get the output file
    output = sys.argv[2]
    # get the list of files
    files = os.listdir(folder)
    # open the output file
    with open(output, 'w') as csvfile:
        # create a csv writer
        writer = csv.writer(csvfile)
        # write the header
        writer.writerow(['Collected By', 'Hardware', 'shared_or_distributed', 'World Size', 'iterations', 'threads', 'tasks', 'time'])
        # loop over the files
        for file in files:
            # open the file
            with open(folder + '/' + file, 'r') as f:
                # read the content
                collected_by = "Seth and Jeremy"
                hardware = "mucluster"

                filename_split = file.split('-')
                filename_split[-1] = filename_split[-1][:-4]
                is_distributed = len(filename_split) == 4

                shared_or_distributed = "distributed" if is_distributed else "shared"
                world_size = filename_split[1]
                iterations = filename_split[3] if is_distributed else filename_split[2]
                threads = filename_split[0]
                tasks = filename_split[2] if is_distributed else 1
                time = f.readline()[5:-5].strip()

            # write the content to the csv
            writer.writerow([collected_by, hardware, shared_or_distributed, world_size, iterations, threads, tasks, time])


if __name__ == '__main__':
    main()