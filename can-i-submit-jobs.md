# Can I submit jobs onto the cluster?

## Connect to the cluster

You need a cluster account for doing this. Please request one via Helpdesk - IT <ithelpdesk@cruk.cam.ac.uk>.

Log in onto the cluster head node
```
ssh my_username@clust1-headnode.cri.camres.org
```

Go to your scratch space
```
cd /scratcha/xxlab/my_username
```

## Create a very simple job in a shell script

We are going to submit a very simple job to the cluster using the `echo` command and a shell script `job.sh`.

Create a `job.sh` file in your scratch space using `nano` and containing these lines below. Do replace `/scratcha/xxlab/my_username` by your scratch space
```
#!/bin/sh
#SBATCH --partition general
#SBATCH --mem 512
#SBATCH --job-name hello_world
#SBATCH --output /scratcha/xxlab/my_username/hello_world.%j.out

echo Hello World from the cluster!
```

## Submit your job

Submit your job to the cluster using the command `sbatch` followed by the name of the script
```
sbatch job.sh
```
