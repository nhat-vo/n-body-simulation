make benchmark
for n_threads in 2 4 8; do
	for n_bodies in 5 50 500 5000; do
		for algo in 1 2 4; do
			./main $algo $n_threads $n_bodies
		done
		echo ""
	done
	echo ""
done
