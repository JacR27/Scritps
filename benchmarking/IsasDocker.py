#!/usr/bin/env python

import argparse
import os.path
import subprocess
import tempfile
import shutil
import string

__author__ = 'pcherng'

def GetOptions():
    """Parses options from command-line arguments

    :rtype : object
    :return: Application options
    """
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-r', required=True, action='store', dest='run_folder',
                            help='input run folder path')
    arg_parser.add_argument('-a', required=False, action='store', dest='analysis_folder',
                            help='output folder path')
    arg_parser.add_argument('-g', required=False, action='store', dest='genome_folder',
                            default='/illumina/development/Isis/Genomes', help='genome folder path')
    docker_image_group = arg_parser.add_mutually_exclusive_group(required=True)
    docker_image_group.add_argument('-i', action='store', dest='isis_docker_name',
                                    help='isis docker image name')
    docker_image_group.add_argument('-t', action='store', dest='isis_docker_tar',
                                    help='isis docker image tar path')
    arg_parser.add_argument('--extra-paths', nargs='*', required=False, action='store', dest='extra_paths',
                            help='extra docker run command mount paths')
    arg_parser.add_argument('--extra-docker-args', nargs='*', required=False, action='store', dest='extra_docker_run_args',
                            help='extra docker run command args')
    arg_parser.add_argument('--extra-isis-args', required=False, action='store', dest='extra_isis_args',
                            help='extra docker run command args')
    arg_parser.add_argument('--use-scratch', required=False, action='store_true', dest='use_scratch',
                            help='toggles use of scratch directory')
    arg_parser.add_argument('--scratch-dir', required=False, action='store', dest='scratch_folder',
                            default='/scratch', help='scratch folder path')
    arg_parser_options = arg_parser.parse_args()
    if arg_parser_options.use_scratch and not arg_parser_options.analysis_folder:
        arg_parser.error("-a option must be used with --use-scratch option")
    if arg_parser_options.genome_folder == None:
        arg_parser_options.genome_folder = "/illumina/development/Isis/Genomes"
    return arg_parser_options


def RunBashCommand(command):
    """Runs a bash command

    :rtype : None
    :param command: List of command arguments to run with bash
    """
    bash_command = ['bash', '-c', " ".join(command)]
    print 'Running command: ' + ' '.join(bash_command)
    subprocess.check_call(bash_command)


def LoadDockerImage(isis_docker_tar):
    """Loads a Docker image tar file

    :rtype : None
    :param isis_docker_tar: Path to the Isis Docker tar file
    """
    docker_import_command = ['sudo', 'docker', 'load', '--input="%s"' % isis_docker_tar]
    RunBashCommand(docker_import_command)


def ConvertIsisTarFilenameToDockerImageName(isis_docker_tar_path):
    """Converts an Isis Docker image tar filename to a Docker image name

    Only works with Isis tar files that follow the correct naming convention:
     Isis Docker tar filename - <repository_name>-<version_number>.tar
     Isis Docker image name - illumina/<repository_name>:<version_number>

    :rtype : string
    :param isis_docker_tar_path: Path to the Isis Docker tar file
    :return: Name of Isis Docker image
    """
    tar_file_name = os.path.basename(isis_docker_tar_path)
    tar_file_name_without_extension = os.path.splitext(tar_file_name)[0]
    dashPos = tar_file_name_without_extension.rfind("-")
    repo_name = tar_file_name_without_extension[:dashPos]
    tag = tar_file_name_without_extension[dashPos + 1:]
    image_name = 'illumina/%s:%s' % (repo_name, tag)
    return image_name


def RunDockerImage(isis_docker_name, run_folder, analysis_folder, genomes_folder, extra_paths, extra_isis_args,
                     extra_docker_run_args):
    """Starts a Docker container and runs Isis

    :rtype : None
    :param isis_docker_name: Name of the Isis Docker image
    :param run_folder: Path to the input Illumina run folder
    :param analysis_folder: Path to the output analysis folder (optional)
    :param genomes_folder: Path to the input top-level genomes folder (optional)
    :param extra_paths: List of extra paths to be read-only mounted to the Docker container
    :param extra_args: Extra arguments to pass to the "docker run" command
    """
    docker_run_command = ['sudo', 'docker', 'run', '-t'] 
    docker_run_command.append("-v /etc/localtime:/etc/localtime:ro") # to use local time zone for logging
    if analysis_folder:
        docker_run_command.append('-v "%s":"%s":ro' % (run_folder, run_folder))
        docker_run_command.append('-v "%s":"%s":rw' % (analysis_folder, analysis_folder))
    else:
        docker_run_command.append('-v "%s":"%s":rw' % (run_folder, run_folder))

    if genomes_folder:
        docker_run_command.append('-v "%s":/illumina/development/Isis/Genomes:ro' % genomes_folder)

    if extra_paths:
        for extra_path in extra_paths:
            docker_run_command.append('-v "%s":"%s":ro' % (extra_path, extra_path))

    if extra_docker_run_args:
        for extra_docker_run_arg in extra_docker_run_args:
            docker_run_command.append(extra_docker_run_arg)

    docker_run_command.extend([isis_docker_name, 'mono', '/opt/illumina/Isis/Isis.exe', '-r', '"%s"' % run_folder])

    if analysis_folder:
        docker_run_command.extend(['-a', '"%s"' % analysis_folder])

    if extra_isis_args:
        docker_run_command.append(extra_isis_args)

    RunBashCommand(docker_run_command)


def CopyFolder(source_path, destination_path):
    """Recursively copies a folder

    :rtype : None
    :param source_path: Path to the source folder to be copied
    :param destination_path: Path to the destination folder where the source folder be copied to
    """
    RunBashCommand(["rsync", "-avuh", os.path.join(source_path, ''), destination_path])


def GetTempFolder(parent_folder):
    """Creates a temp folder

    :rtype : string
    :param parent_folder: Path to the parent folder where the temp folder will be created
    :return: Path to the newly created temp folder
    """
    temp_folder_path = tempfile.mkdtemp(dir=parent_folder)
    print 'Created temp folder: %s' % temp_folder_path
    return temp_folder_path


def DeleteFolder(folder_path):
    """Recursively deletes a folder

    :rtype : None
    :param folder_path: Path to the folder to be deleted
    """
    print 'Deleting folder: %s' % folder_path
    shutil.rmtree(folder_path)


if __name__ == '__main__':
    options = GetOptions()

    if options.isis_docker_tar:
        LoadDockerImage(os.path.abspath(options.isis_docker_tar))
        docker_image_name = ConvertIsisTarFilenameToDockerImageName(options.isis_docker_tar)
    else:
        docker_image_name = options.isis_docker_name

    # When scratch mode is enabled:
    # -Create a temp folder in the scratch folder
    # -Run Isis with -a set to the temp folder
    # -Copy the temp folder to the desired analysis folder
    # -Delete the temp folder
    if options.use_scratch:
        scratch_analysis_folder = GetTempFolder(os.path.abspath(options.scratch_folder))
        RunDockerImage(docker_image_name,
                         os.path.abspath(options.run_folder),
                         scratch_analysis_folder,
                         os.path.abspath(options.genome_folder) if options.genome_folder else None,
                         [os.path.abspath(path) for path in options.extra_paths] if options.extra_paths else None,
                         options.extra_isis_args,
                         options.extra_docker_run_args)
        CopyFolder(scratch_analysis_folder, os.path.abspath(options.analysis_folder))
        DeleteFolder(scratch_analysis_folder)
    # When scratch mode is disabled, just run Isis with outputs set to desired analysis folder or default location
    else:
        RunDockerImage(docker_image_name,
                         os.path.abspath(options.run_folder),
                         os.path.abspath(options.analysis_folder) if options.analysis_folder else None,
                         os.path.abspath(options.genome_folder) if options.genome_folder else None,
                         [os.path.abspath(path) for path in options.extra_paths] if options.extra_paths else None,
                         options.extra_isis_args,
                         options.extra_docker_run_args)
