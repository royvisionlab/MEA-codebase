# The spike sorting pipeline


## Using SBATCH jobs on Hyak
```console
[user@klone-login01 utilities]$ cd /gscratch/retina/GitRepos/manookin-lab/MEA/src/utilities
```

```console
[user@klone-login01 utilities]$ sbatch --export=EXP='20221228C',CHUNK='chunk1',FILES='data028 data029 data030 data031 data032' prepare_data.slurm
```

### Running Kilosort
```console
[user@klone-login01 utilities]$ sbatch --export=EXP='20221228C',CHUNK='chunk1' sort_kilosort.slurm
```

### Run YASS
To run YASS, you will need to create a subdirectory with the name of the chunk and move the .bin, config.yaml, and geom.txt files to that directory.
```console
[user@klone-login01 utilities]$ mkdir /gscratch/scrubbed/retina/data/sorted/20221228C/chunk1
[user@klone-login01 utilities]$ mv /gscratch/scrubbed/retina/data/sorted/20221228C/chunk1.bin /gscratch/scrubbed/retina/data/sorted/20221228C/chunk1
[user@klone-login01 utilities]$ cp /gscratch/scrubbed/retina/data/sorted/20221228C/config.yaml /gscratch/scrubbed/retina/data/sorted/20221228C/chunk1
[user@klone-login01 utilities]$ cp /gscratch/scrubbed/retina/data/sorted/20221228C/geom.txt /gscratch/scrubbed/retina/data/sorted/20221228C/chunk1
```

Check that the config file has the proper .bin name (e.g., chunk1.bin) and edit if necessary. 
```console
[user@klone-login01 utilities]$ nano /gscratch/scrubbed/retina/data/sorted/20221228C/chunk1/config.yaml
```
Exit the editor using Ctrl-X and select 'y' to save any changes made.

Now, you're ready to spin up a YASS batch process.

```console
[user@klone-login01 utilities]$ sbatch --export=EXP='20221228C',CHUNK='chunk1' sort_yass.slurm
```
