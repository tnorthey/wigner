#!/bin/bash

# usage ./workflow 1000 1 0; run generate_pool with N=1000
# usage ./workflow 1000 0 1; run fit_pool with N=1000
# usage ./workflow 1000 1 1; run both with N=1000

./clear_reset.sh

N=$1
if [ $2 = 1 ]; then
  echo "generate_pool.py ..."
  python3 generate_pool.py "initconds.xyz" $N
  cp pool.npz "results/pool_$N.npz"
fi

if [ $3 = 1 ]; then
  echo "fit_pool.py ..."
  python3 fit_pool.py "results/pool_$N.npz"
fi

