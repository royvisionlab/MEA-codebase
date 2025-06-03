# Set up the tools on your Mac.

I highly recommend using conda environments for setting up these analysis tools. This makes it possible to delete/upgrade/reinstall environments if needed. Conda can be installed on your computer using the following links: [Anaconda](https://docs.anaconda.com/free/anaconda/install/mac-os/) or [Miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/macos.html), which is a lightweight version of Anaconda.

You will need a valid gcc compiler to build/run the code in this repository. On MacOS, you can install the necessary tools by opening a terminal session and typing the following command:
```console
xcode-select --install
```

After installing conda, you need to create a conda environment for your analysis pipeline (NEVER use the base environment for this purpose as you would need to completely uninstall and reinstall conda if something goes awry!). In the steps below, we will create an environment named "mea" for our analysis. You can name the environment whatever you like. You can take one of the following steps to install the environment.

## Step 1a: Create an environment directly using the .yml file.
Navigate to the ./MEA/environments folder. This will vary depending on your current path in the terminal.
```console
cd ./MEA/environments
```

Install the environment using the provided .yml file.
```console
conda env create -f mac_mea.yml
```

Skip to [Step 2: Install the pipeline utilities](#step-2-install-the-pipeline-utilities).

## Step 1b: Install the environment and components manually. (Do not run this step if you used Step 1a to create the environment.)
```console
conda create --name mea python=3.9
```

Activate the conda environment.
```console
conda activate mea
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
conda install pytorch::pytorch torchvision torchaudio -c pytorch
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
conda activate mea
```

```console
 cd ../../../chichilnisky-lab/artificial-retina-software-pipeline/utilities/
```
 ```console
 pip install .
```

## Step 3: Install the analysis tools from the MEA repository.
```console
cd ../../../Manookin-Lab/MEA/src/analysis
```
```console
pip install -e .
```
```console
pip install .
```
