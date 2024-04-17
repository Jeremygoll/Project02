# things to find crossover in
# number of threads
# number of nodes
# world size

# TODO: make a python script that puts the outputs into a csv file

# shared
for threads in $(seq 1 1 8); do
    for world_size in 8 16 32 64 128 256; do
        for iter in 1 2 3; do
            echo $threads threads, $world_size world_size, $iter iteration
            sbatch -o /home/colemans/project02/output/$threads-$world_size-$iter.out \
                -e /home/colemans/project02/errors/$threads-$world_size-$iter.err \
                --export=threads=$threads,filename=big_$world_size.txt \
                --cpus-per-task=$threads ./shared_sbatch.sh
        done
    done
done

# distributed
for threads in $(seq 1 1 8); do
    for world_size in 8 16 32 64 128 256; do
        for nodes in 1 2 3 4; do
            for iter in 1 2 3; do
                echo $threads threads, $world_size world_size, $nodes nodes, $iter iteration
                sbatch -o /home/colemans/project02/output/$threads-$world_size-$nodes-$iter.out \
                    -e /home/colemans/project02/errors/$threads-$world_size-$nodes-$iter.err \
                    --export=threads=$threads,filename=big_$world_size.txt,nodes=$nodes \
                    --nodes=$nodes --cpus-per-task=$threads ./distributed_sbatch.sh
            done
        done
    done
done