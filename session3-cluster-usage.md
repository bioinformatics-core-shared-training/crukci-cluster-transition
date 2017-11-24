# Session 3: Usage of cluster


## Learning Objectives

- Submit jobs to the cluster using SLURM
- Get cluster file systems mounted on your local macOS
- Run quality control on the cluster


## Submit a job

We did submit a very simple job to the cluster using the `echo` command and a shell script `job.sh` containing specific Slurm instructions. Please check the instructions on [Can I submit jobs onto the cluster?](can-i-submit-jobs.md) and run it again to make sure you can submit jobs to the cluster and understand how to do it.

In the `job.sh` script, we have specific `SBATCH` instructions:
- `--partition`: select a specific queue to submit the job
- `--mem`: specify the real memory required e.g. `2GB` is `2048`
- `--job-name`: specify a name for the job allocation
- `--output`: connect the batch script's standard output directly to the file name specified. By default both standard output and standard error are directed to the same file. The filename is `hello_world.%j.out` where the `%j` is replaced by the job ID.

See [sbatch](http://slurm.schedmd.com/sbatch.html) man page for all the options and explanation on submitting a batch script to Slurm.

> :computer: **EXERCISE** Go to your Terminal window, or open a new one and go to `session1-data/nelle/`.
>
> - Copy `molecules/` onto your scratch space on the cluster using `scp -r`
> - Submit `middle.sh` script to Slurm to extract lines 20-23 of `octane.pdb`
> - Write a loop to submit jobs for all PDB files by modifying `job.sh` to take one command line argument which will be given at each step of the loop
>
> :tada: Congratulations! :thumbsup: You did it! :wink:


## Get cluster file systems mounted on your macOS

[FUSE for macOS](https://osxfuse.github.io/) allows you to extend macOS's native file handling capabilities via third-party file systems. Combining with [SSHFS](https://github.com/osxfuse/osxfuse/wiki/SSHFS), it gets file system accessible via `ssh` mounted directly on your macOS computer.

Install:
- Latest release of [FUSE for macOS](https://github.com/osxfuse/osxfuse/releases) by downloading the `.dmg` file
- Latest release of [SSHFS](https://github.com/osxfuse/sshfs/releases) by downloading the `.pkg` file

Open a Terminal window, and create a mount point in your home directory `mnt/scratchb` for example:
```
cd ~                    # go to home directory
mkdir -p mnt/scratchb   # create intermediate directories as required with -p option
cd mnt/scratchb         # go to this directory
pwd                     # return current working directory name
```

Mount your cluster `/mnt/scratchb/my_group/my_username` onto your local machine using:
```
sshfs my_username@clust1-headnode.cri.camres.org:/mnt/scratchb/my_group/my_username /Users/my_username/mnt/scratchb
```
To list all mounted file system, use `mount`:
```
mount
```
To unmount `scratchb`, use:
```
umount /Users/my_username/mnt/scratchb
```

Aliases can be created in your `~/.profile` to save you time and to avoid typing complex commands every time you need them. You could enter these lines at the beginning of your `~/.profile` file using your preferred editor [atom](https://atom.io/) which could be launch from the command line by typing `atom`:
```
atom ~/.profile
```
Add enter these lines into the file:
```
### Aliases
alias mntclustsb='sshfs pajon01@clust1-headnode.cri.camres.org:/mnt/scratchb/my_group/my_username /Users/my_username/mnt/scratchb'
alias umntclustsb='umount /Users/my_username/mnt/scratchb'
```
Save and open a new Terminal window, your new aliases are now available as 'new' commands! :thumbsup:


## FastQ quality control

We are now going back to our own data and check read quality using FastQC. Here we will use the command line version of FastQC to check what kind of options this tool has. This command will display all the parameters you can use when running `FastQC`.
```
/home/bioinformatics/software/fastqc/fastqc-v0.11.5/fastqc --help
```
Now to run FastQC on your downloaded sequencing data in `SLX-ID` folder, you will need to run this command but we are not running this command directly, we will submit this job to the cluster:
```
/home/bioinformatics/software/fastqc/fastqc-v0.11.5/fastqc -o /scratcha/xxlab/my_username/SLX-ID/ --noextract -f fastq my_file.fastq.gz
```
If you need to get sequencing data again, check the steps from [Session 1: Shell - Getting sequencing: Using CRUKCI infrastructure data](session1-shell.md#using-crukci-infrastructure)

The options are:
- `-o FOLDER_NAME`: you can define this way in which folder you want FastQC to output its results
- `–noextrac`t: tells FastQC not to uncompress the output file after creating it
- `-f fastq`: defines that the input file is a fastq file (valid formats are bam, sam, bam_mapped, sam_mapped and fastq)

> :computer: **EXERCISE** Go to your Terminal window, or open a new one and log in onto the cluster head node.
>
> - Navigate to your project data.
> - Create a `job.sh` to run FastQC
> - Send job to the cluster and wait for the results
> - Update your `README.txt` file with what you've done
> - View the html report in a web browser, you may have to copy back this file on your own computer to be able to view it using the `scp` command or mount your scratch space using `sshfs`
>
> :tada: Congratulations! :thumbsup: You did it! :wink:

When running jobs on the cluster, you may wish to keep an eye on the output by using `tail -f` which appends data as the file grows:
```
tail -f my_output_file_name.out
```
Check your job is still running and in the queue using [squeue](http://slurm.schedmd.com/squeue.html), you can combine it with `grep my_username` to only extract information about your jobs:
```
squeue | grep my_username
```
You can display information of all your submitted jobs using [sacct](http://slurm.schedmd.com/sacct.html) but most importantly you need to check that your job has completed using:
```
sacct -j JobID
```
If you wish to kill your job before it completes, run [scancel](http://slurm.schedmd.com/scancel.html) using:
```
scancel JobID
```


## Take home message: everyday cluster commands

```
sbatch job.sh
sacct -j JobID
tail -f my_output_file_name.JobID.out
```

## Reference materials

- [CRUK Bioinformatics Autumn School 2017: Functional Genomics](https://bioinformatics-core-shared-training.github.io/cruk-autumn-school-2017/)
  - [Quality control and artefact removal](https://bioinformatics-core-shared-training.github.io/cruk-autumn-school-2017/Introduction/SS_DB/Materials/Lectures/Lecture2_qualityControl_artefactRemoval_DB.pdf)
