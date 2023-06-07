for algo in {0..3}; do
    for n_bodies in {20,100,500,1000}; do
        for n_threads in {1,2,4,6}; do
            ./main $algo $n_bodies $n_threads
        done
    done
done