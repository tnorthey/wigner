#!/bin/bash

#n=100000
./clear_reset.sh
#echo "converting $n samples to xyz..."
#./initconds2xyz.sh $n
#echo "creating xyz_traj"
#./create_xyz_traj.sh
echo "calculating distances, and x-ray"
./run_get_distances_xray.sh

