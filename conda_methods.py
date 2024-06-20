import os
import sys
import subprocess


class CondaMethods(object):
    @staticmethod
    def is_conda_installed():
        return os.path.exists(os.path.join(sys.prefix, 'conda-meta'))

    @staticmethod
    def is_conda_env_installed(env):
        cmd = ['conda', 'info', '--envs']
        return env in subprocess.check_output(cmd).decode()

    @staticmethod
    def is_env_activated(env):
        return env in os.environ['CONDA_DEFAULT_ENV']

    @staticmethod
    def get_conda_env_path(env):
        cmd = ['conda', 'info', '--envs']
        env_dict = dict()
        env_out = subprocess.check_output(cmd).decode().split('\n')
        for line in env_out:
            if line.startswith('#') or line == '':
                continue
            e, p = line.split(' ', maxsplit=1)
            e = e.strip()
            p = p.strip()
            env_dict[e] = p

        return env_dict[env]

    @staticmethod
    def install_pycoQC_env():
        cmd = ['conda', 'create', '-y', '-n', 'pycoQC', 'pycoQC=3.0.0']
        subprocess.run(cmd)

    @staticmethod
    def install_nbc_env():
        cmd = ['conda', 'create', '-y', '-n', 'nbc', 'filtlong', 'pigz']
        subprocess.run(cmd)

