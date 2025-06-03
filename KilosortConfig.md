# Kilosort Setup

## Install the Kilosort versions that you plan to use.

## Kilosort v4 (Python-based)

### Option 1: Install the development version.
```console
pip install -e ".[gui]"
```

### Option 2: Install the development version.

```console
pip install -e ".[gui]"
```

## Kilosort <= v3 (Matlab-based)


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

