#$ -S /bin/bash
##cd /home/zsu/NASAIRAD/opt_stratafiedsampling
cp e0m0/e0m0.e* /home/zsu/NASAIRAD/opt_stratafiedsampling/microns100_results/opt/level2stratafiedsampling2
cp e0m0/e0m0.o* /home/zsu/NASAIRAD/opt_stratafiedsampling/microns100_results/opt/level2stratafiedsampling2

for x in 1 2 4 8
  do
    for y in {0..4..1}
      do
        cp e${x}m${y}/e${x}m${y}.e* /home/zsu/NASAIRAD/opt_stratafiedsampling/microns100_results/opt/level2stratafiedsampling2
        cp e${x}m${y}/e${x}m${y}.o* /home/zsu/NASAIRAD/opt_stratafiedsampling/microns100_results/opt/level2stratafiedsampling2

      done
  done
