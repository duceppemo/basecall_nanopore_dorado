# Basecall Nanopore Dorado

## Description
This is a wrapper script to streamline nanopore basecalling using `dorado`. It will also perform read QC (pycoQC) and filtering (Filtlong).

Under the hood, it's using `dorado_basecall_server` + `ont_basecall_client` rather than the standalone `dorado` (https://github.com/nanoporetech/dorado). The rationale behind using the basecall server is that can take both fast5 and pod5 as input, performs sample demultuplexing and outputs fastq with "complete" headers (needed for pycoQC).

## Important notes
- This script will only work in "Linux".
- You must have conda and dorado basecall server pre-installed.
- You need to have dorado install location added to your `.bashrc`.
- The other dependencies will be automatically installed via conda during runtime, the first time.
- PycoQC integration is not implemented yet.
- Only minimal read filtering is done (remove bottom 5%).

## Installation
- Miniconda installation (say "yes" to when asked automatically load conda on startup):
```bash
# https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html
# https://docs.anaconda.com/miniconda/
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh

# Close and open a new terminal after installation.
# You should see "(base)" at the beginning of your terminal line. 
```
- Dorado Basecall Server installation:
```bash
# Create folder to hold program
[ -d $HOME/prog ] || mkdir -p $HOME/prog
cd $HOME/prog

# Get and extract dorado basecall server
wget https://cdn.oxfordnanoportal.com/software/analysis/ont-dorado-server_7.3.9_linux64.tar.gz  # Get dorado basecall server
tar zxvf ont-dorado-server_7.3.9_linux64.tar.gz  # Extract archive
rm ont-dorado-server_7.3.9_linux64.tar.gz  # Remove archive

# Make sure we can call dorado from any location from terminal
echo "export PATH=\$PATH:\$HOME/prog/ont-dorado-server/bin" | tee -a $HOME/.bashrc
source $HOME/.bashrc  # Apply changes
```
* Pipeline installation:
```bash
# Create folder to hold program
[ -d $HOME/prog ] || mkdir -p $HOME/prog
cd $HOME/prog

# Clone repo
git clone https://github.com/duceppemo/basecall_nanopore_dorado

# Test pipeline
cd basecall_nanopore_dorado
python basecall_nanopore_dorado.py -h
```
## Usage
```commandline
usage: python basecall_nanopore_dorado.py [-h] -i /path/to/input_folder/ -o /path/to/output_folder/ [-s {minion,promethion}] [-c dna_r9.4.1_450bps_sup.cfg] [-f FLO-MIN106] [-l SQK-LSK109] [-b EXP-NBD104 [EXP-NBD104 ...]] [-d /path/to/barcode_description.tsv] [--min-qscore MIN_QSCORE]
                                          [--port PORT] [-r] [-t 24] [-g "cuda:0"] [-p 2] [-m 114] [-v]

Basecall Nanopore raw data to fastq using Dorado Basecall Server.

options:
  -h, --help            show this help message and exit
  -i /path/to/input_folder/, --input /path/to/input_folder/
                        Folder that contains the fast5 files. Mandatory.
  -o /path/to/output_folder/, --output /path/to/output_folder/
                        Folder to hold the result files. Mandatory.
  -s {minion,promethion}, --sequencer {minion,promethion}
                        Sequencer used. "minion" includes all sequencers except "promethion". Optional. Default is "minion".
  -c dna_r9.4.1_450bps_sup.cfg, --config dna_r9.4.1_450bps_sup.cfg
                        Dorado config file. Typically found in dorado installation folder and ending with ".cfg". This is the preferred methods over choosing a "library kit/flowcell" combination. Both methods are incompatible. Optional.
  -f FLO-MIN106, --flowcell FLO-MIN106
                        Flowcell type used for sequencing. Optional.
  -l SQK-LSK109, --library-kit SQK-LSK109
                        Library kit used. Optional.
  -b EXP-NBD104 [EXP-NBD104 ...], --barcode-kit EXP-NBD104 [EXP-NBD104 ...]
                        Barcoding kit(s). Not using this option will not perform barcode splitting. For multiple barcoding kits, use double quotes and space like this: "EXP-NBD104 EXP-NBA114". Optional
  -d /path/to/barcode_description.tsv, --description /path/to/barcode_description.tsv
                        Tab-separated file with two columns with barcode assignments. First column contains barcode names [barcode01, barcode02, etc.]. Second column contains sample name. Avoid using special characters. Sample file in data folder. Optional.
  --min-qscore MIN_QSCORE
                        Minimum acceptable qscore for a read to be filtered into the PASS folder. Accepted values: [0 .. 30]. Default 10. Optional.
  --port PORT           Port for basecalling service. Default 5555. Optional.
  -r, --recursive       Look for pod5 or fast5 recursively. Optional
  -t 24, --threads 24   Number of threads. Default is maximum available(24). Optional.
  -g "cuda:0", --gpu "cuda:0"
                        GPU device to use. Typically use "cuda:0". Default is "auto". Optional.
  -p 2, --parallel 2    Number of samples to process in parallel for trimming and filtering. Default is 2. Optional.
  -m 114, --memory 114  Memory in GB. Default is 85% of total memory (114). Optional.
  -v, --version         show program's version number and exit
```
## Sample description file
- When the `--description` argument is being used, demultiplexed fastq will be renamed. It's important to know that any barcode detected that are not listed in your `sample description` file will be deleted. The `unclassified` folder is being kept, but not used at the filtering step.
- Look at  the file `sample_description.tsv` file in the `/data` folder for a template for the `--description` argument.
- File format is two tab-separated colums. First column must be the barcode numbers from `barcode01` up to `barcode96`. Second colum is your prefered sample name. Don't use spaces or any special characters (`_` and `-` are OK).
```text
barcode01	my_sample1
barcode02	my_sample2
barcode03	my_sample3
barcode04	my_sample4
```

## Examples
- If you're having trouble getting a working value for the config file name, barcoding kit, library kit or flowcell, you can have a peak at the `kits.py` file content, which contains all the valid names. This file needs to be updated as kits come and go. Please file an issue is you specific kit is not listed and should be. Note that Dorado is not compatible with some older kits.

### Fast5 generated from a multiplexed run with a R9.4.1 flowcell:
```bash
python basecalle_nanopore_dorado.py \
    --input "/data/run1/fast5" \
    --output "$HOME/analyses/run1_dorado" \
    --config "dna_r9.4.1_450bps_sup.cfg" \
    --description "/data/run1/bc.txt" \
    --barcode-kit "EXP-NBD104" \
    --recursive
```
### Pod5 generated from a multiplexed run with a R10.4.1 flowcell and v14 chemistry:
```bash
python basecalle_nanopore_dorado.py \
    --input "/data/run2/pod5" \
    --output "$HOME/analyses/run2_dorado" \
    --config "dna_r10.4.1_e8.2_400bps_5khz_sup.cfg" \
    --description "/data/run2/bc.txt" \
    --barcode-kit "SQK-NBD114-24" \
    --recursive
```
