# MEA: Analysis tools for data collected on UC Santa Cruz multi-electrode arrays.

### [Install Kilosort on your server](KilosortConfig.md)

## Set up paths for Matlab-based versions of Kilosort (version <= 3)
If you intend to run Matlab-based versions of kilosort (version<= 3), you will need to copy the the kilosort_paths.m file to your user directory and modify that file with the correct paths to the directories indicated in the file. You can find the path to your user directory by typing ***userpath*** in the Matlab console window. I recommend that you add path to the required NPY-Matlab and pipeline_utilities directories in your startup.m file located in the same (user root) directory. You can edit the startup.m file by typing the following into your Matlab console:

```console
edit startup
```

Then add the following lines to the file with the correct paths set and save the file.

```console
addpath('/home/mike/Documents/git_repos/npy-matlab/npy-matlab/');
addpath('/home/mike/Documents/git_repos/manookin-lab/MEA/src/pipeline_utilities/');
```

The NPY-Matlab repository is available on GitHub. Change to the directory where you want to store the repository and type the following line.

```console
git clone https://github.com/kwikteam/npy-matlab.git
```


## Set the proper paths to your data and to the Vision.jar file in the .bashrc file.
```console
sudo nano ~/.bashrc

# Paths for analysis pipeline.
export TEMPORARY_SORT_PATH='/data/data/'
export LITKE_PATH='/data/data/'
export SORTED_SPIKE_PATH='/data/data/sorted'
export KILOSORT_TTL_PATH='/data/data/sorted'
export RAW_DATA_PATH='/data/data/raw'
export VISIONPATH='/home/mike/Documents/git_repos/manookin-lab/MEA/src/Vision7_for_2015DAQ/Vision.jar'
export LAB_NAME='RiekeManookin' # Your lab name goes here.
```

## For Hyak (Manookin/Rieke labs)
Each user who wants to run the pipeline code must edit their .bashrc file as follows.
```console
nano ~/.bashrc

# Paths for analysis pipeline.
export LITKE_PATH='/gscratch/scrubbed/retina/data/'
export SORTED_SPIKE_PATH='/gscratch/retina/data/sorted'
export KILOSORT_TTL_PATH='/gscratch/retina/data/sorted'
export RAW_DATA_PATH='/gscratch/scrubbed/retina/data/raw'
export VISIONPATH="/gscratch/retina/GitRepos/manookin-lab/MEA/src/Vision7_for_2015DAQ/"
export LAB_NAME='RiekeManookin' # Your lab name goes here.
```


Ctrl-X to save the file and 'Y'

Test it out
```console
source ~/.bashrc
echo $TEMPORARY_SORT_PATH
```


I highly recommend using conda environments for setting up these analysis tools. This makes it possible to delete/upgrade/reinstall environments if needed. Conda can be installed on your computer using the following links: [Anaconda](https://docs.anaconda.com/free/anaconda/install/mac-os/) or [Miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/macos.html), which is a lightweight version of Anaconda.

You will need a valid gcc compiler to build/run the code in this repository. After installing conda, you need to create a conda environment for your analysis pipeline (NEVER use the base environment for this purpose as you would need to completely uninstall and reinstall conda if something goes awry!). In the steps below, we will create an environment named "kilosort" for our analysis. You can name the environment whatever you like. You can take one of the following steps to install the environment.

## Step 1a: Create an environment directly using the .yml file.
Navigate to the ./MEA/hyak/containers/ folder. This will vary depending on your current path in the terminal.
```console
cd ./MEA/hyak/containers/
```

Install the environment using the provided .yml file.
```console
conda env create -f kilosort.yml
```

Skip to [Step 2: Install the pipeline utilities](#step-2-install-the-pipeline-utilities).

## Step 1b: Install the environment and components manually. (Do not run this step if you used Step 1a to create the environment.)
```console
conda create --name kilosort python=3.9
```

Activate the conda environment.
```console
conda activate kilosort
```

Install the prerequisite libraries. 
Note: Scipy version 1.13 has a serious bug. You need to downgrade to 1.12.

```console
conda install conda-forge::hdf5storage
```
```console
pip install progressbar2
```
```console
conda install conda-forge::cython
```
```console
pip install opencv-python
```
```console
pip install -U matplotlib
```
```console
pip install -U scikit-learn
```
```console
pip install pandas
```
```console
pip install pybind11
```
```console
conda install pytorch pytorch-cuda=11.8 -c pytorch -c nvidia
```
```console
conda install conda-forge::tqdm
```
```console
pip install -Iv scipy==1.12.0
```
```console
pip install -Iv numpy==1.25.0
```

## Step 2: Install the pipeline utilities.
If you haven't already done so, activate your conda environment.

```console
conda activate kilosort
```

```console
 cd ../../../chichilnisky-lab/artificial-retina-software-pipeline/utilities/
```
 ```console
 pip install .
```

```console
 cd ../../vision-convert
```
 ```console
 pip install .
```

## Step 3: Install the analysis tools from the MEA repository.
```console
cd ../../../Manookin-Lab/MEA/src/analysis
```
```console
python setup.py build_ext -i
```
```console
pip install .
```

# Kilosort setup.

## Kilosort v4 (Python-based)

```console
pip install -e ".[gui]"
```

## Kilosort <= v3 (Matlab-based)

