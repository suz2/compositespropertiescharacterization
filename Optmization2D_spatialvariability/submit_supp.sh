#$ -S /bin/bash
cd e0_cpu1
##qsub -l h=node009 initial.sh
for x in 0 1 2 4 8
      do
	for z in 1
	  do
	    cd ../e${x}_cpu${z}
	    qsub -l h=node009 initial.sh        
	  done

        for z in 2 
          do
            cd ../e${x}_cpu${z}
            qsub -l h=node010 initial.sh
          done

        for z in 3
          do
            cd ../e${x}_cpu${z}
            qsub -l h=node011 initial.sh
          done

        for z in 4
          do
            cd ../e${x}_cpu${z}
            qsub -l h=node012 initial.sh
          done

  done
##cd ..
