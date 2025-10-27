# KakuEtAl2025
Data analysis code for [**Kaku, Liu et al. A brainstem map of orofacial rhythms**](https://www.biorxiv.org/content/10.1101/2025.01.27.635041v1)

## Installation
#### NWB to MAT converter (MATLAB >=2024b required)
- Install nwb2mat using [setup.exe](https://github.com/NuoLiLabBCM/KakuEtAl2025/tree/main/nwb2mat). The default installation directory is `C:\Program Files\nwb2mat`.
- Add `your\MATLAB\Path\bin\win64` to PATH in environment variables. ([Help](https://www.java.com/en/download/help/path.html))


#### Environment
The environment can be installed using requirements.yml. 

`conda env create -f requiements.yml`




## Usage
### Step 1
Download the NWB dataset from [DANDI archive](https://dandiarchive.org/dandiset/001619).

### Step 2
Convert NWB files to MAT files for processing using the nwb2mat from command line. In the command prompt use `"C:\Program Files\nwb2mat\nwb2mat.exe" -d path\to\folder\containing\NWB`

### Step 3
Run process_mat.py (Change `input_dir` to the folder containing MAT files generated in step 2). This will create and save MAT files with the suffix `'_processed'` containing key analysis results from the [publication](https://www.biorxiv.org/content/10.1101/2025.01.27.635041v1).

### Step 4
Run generate_figures.py and generate_activity_maps.m to plot key figures from the [publication](https://www.biorxiv.org/content/10.1101/2025.01.27.635041v1).

Download TIFF required for plotting activity maps [here](https://drive.google.com/file/d/18zMYyCuro2lHrILHIM2xHU-94VT9LCMF/view?usp=drive_link).