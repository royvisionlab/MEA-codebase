
## Kilosort pipeline using BASH shell

### Step 1: Prepare the data for spike sorting

```console
bash prepare_data.sh 20240130C chunk4 data015 data016 data017 data018 data019 data020 data021 data022 data023 data024 data027
```

### Step 2: Run spike sorting in Kilosort2

```console
matlab -nodisplay
run_kilosort('20240130C','chunk4')
```

### Step 3: Process the spikes from Kilosort
```console
bash process_spikes.sh 20240130C chunk4 kilosort2 16
```

### Step 4: Run deduplication
```console
bash deduplication.sh 20240130C chunk4 kilosort2 12 data019
```

### Step 4: Run the STA analysis
```console
cd ../analysis/protocol/typing/
bash analyze_chunk.sh 20240130C chunk4 kilosort2 data019
```

## Kilosort pipeline using SBATCH/SLURM

### Step 1: Prepare the data for spike sorting.

```console
[user@klone-login01 analysis]$ sbatch --export=EXP='20221006C',CHUNK='chunk1' sort_kilosort.slurm
```

### Step 2: Run spike sorting in Kilosort2

```console
[user@klone-login01 analysis]$ sbatch --export=EXP='20221006C',CHUNK='chunk1',ALG='kilosort2' process_spikes.slurm
```

### Step 3: Run the STA analysis on the FastNoise data files.

```console
[user@klone-login01 analysis]$ sbatch --export=EXPERIMENT_DATE='20221006C',CHUNK_NAME='chunk1',SORT_ALGORITHM='kilosort2',DATA_FILES='data001 data009' analyze_chunk.slurm
```

