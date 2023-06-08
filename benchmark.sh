make main
for n_threas in 1 2 4 8; do
    for n_bodies in 5 50 500 5000; do
        for algo in 0 1 2 3 4; do
            ./main $algo $n_threas $n_bodies
        done
    done
done