import os
from argparse import ArgumentParser
from multiprocessing import cpu_count
from psutil import virtual_memory
from basecall_nanopore_dorado_methods import Methods
from importlib.resources import files
from kits import Kits
from conda_methods import CondaMethods

# pandas=2.2.2 psutil=5.9.8


__author__ = 'duceppemo'
__version__ = '0.1'


# TODO: create log of versions and time stamps for QA

class Basecaller(object):
    def __init__(self, args):
        # I/O
        self.input = os.path.abspath(args.input)
        self.output_folder = os.path.abspath(args.output)

        # Performance
        self.cpu = args.threads
        self.parallel = args.parallel
        self.mem = args.memory

        # Dorado related
        self.gpu = args.gpu
        self.description = args.description
        self.barcode_kit = args.barcode_kit
        if self.barcode_kit:
            self.barcode_kit = self.barcode_kit[0].split()
        self.sequencer = args.sequencer
        self.config = args.config
        self.flowcell = args.flowcell
        self.library_kit = args.library_kit
        self.min_qscore = args.min_qscore
        self.port = args.port
        self.recursive = args.recursive
        self.workflows = str(files('data').joinpath('workflows.tsv'))

        # Data
        self.sample_dict = dict()

        # Conda
        self.pycoQC_env_path = ''
        self.nbc_env_path = ''

        # Run
        self.run()

    def run(self):

        ##################
        #
        # Checks
        #
        ##################

        print('Checking a few things...')

        # Check if number of CPU and memory requested are valid
        self.cpu, self.parallel = Methods.check_requested_cpus(self.cpu, self.parallel)
        self.mem = Methods.check_requested_mem(self.mem)

        # Check I/O
        Methods.check_input_folder(self.input)
        Methods.check_raw_exist(self.input)
        Methods.check_config(self.config, self.flowcell, self.sequencer, self.library_kit)
        if self.barcode_kit:
            Methods.check_barcode(self.barcode_kit, self.description)

        # Check software
        Methods.check_dorado_installed()
        Methods.check_conda_installed()

        # Check environments
        if not CondaMethods.is_conda_env_installed('pycoQC'):
            CondaMethods.install_pycoQC_env()
        self.pycoQC_env_path = CondaMethods.get_conda_env_path('pycoQC')

        if not CondaMethods.is_conda_env_installed('nbc'):
            CondaMethods.install_nbc_env()
        self.pycoQC_env_path = CondaMethods.get_conda_env_path('nbc')

        print('\tAll checks passed')

        ##################
        #
        # Preparing outputs
        #
        ##################

        # Step completion report files
        done_basecalling = self.output_folder + '/done_basecalling'
        done_qc = self.output_folder + '/done_QC'
        done_filtering = self.output_folder + '/done_filtering'

        # Output folders to create
        basecalled_folder = self.output_folder + '/1_basecalled/'
        qc_folder = self.output_folder + '/2_qc/'
        filtered_folder = self.output_folder + '/3_filtered/'

        # Create output folder
        Methods.make_folder(self.output_folder)

        ##################
        #
        # 1- Basecalling
        #
        ##################

        if not os.path.exists(done_basecalling):
            # Retrieve proper configuration file
            if not self.config:
                dorado_conf = Methods.get_dorado_config(self.flowcell, self.library_kit, self.sequencer, self.workflows)
            else:
                dorado_conf = self.config

            # Basecall
            # Methods.is_basecall_server_running()
            print('Starting basecalling server...')
            p = Methods.start_dorado_basecall_server(basecalled_folder, self.config, self.port, self.gpu)

            print('Basecalling with Dorado')
            Methods.run_dorado(self.input, basecalled_folder, dorado_conf, self.recursive,
                               self.gpu, self.barcode_kit, self.min_qscore, self.port)
            # Terminate the basecalling server
            p.terminate()

            # Merge all fastq per barcode, if more than one file present
            Methods.merge_rename_fastq(basecalled_folder, self.barcode_kit)

            if self.description:
                sample_dict = Methods.parse_samples(self.description)
                Methods.rename_barcode(sample_dict, basecalled_folder)  # Also remove extra barcode folders

            # Create "done" file for resuming purposes
            Methods.flag_done(done_basecalling)
        else:
            print('Skipping basecalling. Already done.')

        # Update sample_dict after extracting, only keep "pass" files
        self.sample_dict['basecalled'] = Methods.get_files(basecalled_folder, 'pass.fastq.gz')

        # Remove "unclassified" for next step if barcodes used
        if self.barcode_kit:
            self.sample_dict['basecalled'].pop('unclassified')

        ##################
        #
        # 2- QC
        #
        ##################

        if not os.path.exists(done_qc):
            print('Performing read QC with PycoQC...')
            Methods.fastq_to_seq_summary(basecalled_folder, qc_folder, self.cpu, 'pycoQC')
            Methods.run_pycoQC(basecalled_folder, qc_folder)
            Methods.flag_done(done_qc)
        else:
            print('Skipping QC. Already done.')

        ##################
        #
        # 3- Filter reads
        #
        ##################

        # Get reference size
        if not os.path.exists(done_filtering):
            print('Filtering lower quality reads with Filtlong...')
            Methods.run_filtlong_parallel(self.sample_dict['basecalled'], filtered_folder, self.parallel, 'nbc')
            Methods.flag_done(done_filtering)
        else:
            print('Skipping filtering. Already done.')

        # Update sample_dict after trimming
        self.sample_dict['filtered'] = Methods.get_files(filtered_folder, '.fastq.gz')

        ##################
        #
        # Done
        #
        ##################

        print('DONE!')


if __name__ == "__main__":
    max_cpu = cpu_count()
    max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB

    parser = ArgumentParser(description='Basecall Nanopore raw data to fastq using Dorado Basecall Server.')
    parser.add_argument('-i', '--input', metavar='/path/to/input_folder/',
                        required=True, type=str,
                        help='Folder that contains the fast5 files. Mandatory.')
    parser.add_argument('-o', '--output', metavar='/path/to/output_folder/',
                        required=True, type=str,
                        help='Folder to hold the result files. Mandatory.')
    parser.add_argument('-s', '--sequencer',
                        required=False, type=str,
                        choices=['minion', 'promethion'],
                        help='Sequencer used. "minion" includes all sequencers except "promethion". '
                             'Optional. Default is "minion".')
    parser.add_argument('-c', '--config', metavar='dna_r9.4.1_450bps_sup.cfg',
                        required=False, type=str,
                        help='Guppy config file. Typically found in "/opt/ont/guppy/data" ending with ".cfg". '
                             'This is the prefered methods over choosing a "library kit/flowcell" combination. '
                             'Both methods are uncompatible. Optional.')
    parser.add_argument('-f', '--flowcell', metavar='FLO-MIN106',
                        required=False, type=str,
                        help='Flowcell type used for sequencing. Optional.')
    parser.add_argument('-l', '--library-kit', metavar='SQK-LSK109',
                        required=False, type=str,
                        help='Library kit used. Optional.')
    parser.add_argument('-b', '--barcode-kit', metavar='EXP-NBD104',
                        required=False, type=str, nargs='+',
                        help='Barcoding kit(s) used. Use "unknown" if you know barcodes were used, but do not know '
                             'which kit. Not using this option will not perform barcode splitting. For multiple '
                             'barcoding kits, use double quotes and space like this: "EXP-NBD104 EXP-NBA114". Optional')
    parser.add_argument('-d', '--description', metavar='/path/to/barcode_description.tsv',
                        required=False, type=str,
                        help='Tab-separated file with two columns with barcode assignments. '
                             'First column contains barcode names [barcode01, barcode02, etc.]. '
                             'Second column contains sample name. Avoid using special character. '
                             'Sample file in data folder. Optional.')
    parser.add_argument('--min-qscore', type=int, default=10, required=False,
                        help='Minimum acceptable qscore for a read to be filtered into the PASS folder.	'
                             'Accepted values: [0 .. 30]. Default 10. Optional.')
    parser.add_argument('--port', type=int, default=5555, required=False,
                        help='Port for basecalling service. Default 5555. Optional.')
    parser.add_argument('-r', '--recursive',
                        action='store_true',
                        help='Look for fast5 recursively. Useful if fast5 are in multiple sub-folders. Optional')
    parser.add_argument('-t', '--threads', metavar=str(max_cpu),
                        required=False, type=int, default=max_cpu,
                        help='Number of threads. Default is maximum available({}). Optional.'.format(max_cpu))
    parser.add_argument('-g', '--gpu', metavar='"cuda:0"',
                        required=False, type=str, default='auto',
                        help='GPU device to use. Typically use "cuda:0". Default is "auto". Optional.')
    parser.add_argument('-p', '--parallel', metavar='2',
                        required=False, type=int, default=2,
                        help='Number of samples to process in parallel for trimming and filtering. '
                             'Default is 2. Optional.')
    parser.add_argument('-m', '--memory', metavar=str(max_mem),
                        required=False, type=int, default=max_mem,
                        help='Memory in GB. Default is 85%% of total memory ({}). Optional.'.format(max_mem))
    parser.add_argument('-v', '--version', action='version',
                        version=f'{os.path.basename(__file__)}: version {__version__}')

    # Get the arguments into an object
    arguments = parser.parse_args()

    Basecaller(arguments)
