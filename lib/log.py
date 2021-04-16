"""Common logging functions."""

import os
import signal
import subprocess
import sys
import tempfile
from datetime import datetime

import psutil

from . import db, util

DEBUG = 10
INFO = 20
ERROR = 30
FATAL = 40

LEVEL = {
    'debug': DEBUG,
    'info': INFO,
    'error': ERROR,
    'fatal': FATAL,
}


class Logger:
    """Common logging functions."""

    def __init__(self, log_file, log_level):
        self.log_file = log_file
        self.log_level = LEVEL[log_level]
        self.file_handle = open(log_file, 'a') if log_file else None

    def __del__(self):
        if self.file_handle:
            self.file_handle.close()

    def header(self):
        """Log header information."""
        self.info('#' * 80)
        self.info('aTRAM version: {}'.format(db.ATRAM_VERSION))
        self.info('Python version: {}'.format(' '.join(sys.version.split())))
        self.info(' '.join(sys.argv[:]))

    def subcommand(self, cmd, temp_dir, timeout=None):
        """
        Call a subprocess and log the output.

        Note: stdout=PIPE is blocking and large logs cause a hang.
        So we don't use it.
        """
        self.debug(cmd)

        error = None

        with tempfile.NamedTemporaryFile(mode='w', dir=temp_dir) as log_output:
            with subprocess.Popen(
                    cmd, shell=True, stdout=log_output, stderr=log_output) as proc:
                try:
                    proc.communicate(timeout=timeout)

                # Catch any error and kill all child processes
                except Exception as err:  # pylint: disable=broad-except
                    pid = psutil.Process(proc.pid)

                    error = err
                    self.error('Exception: {}'.format(err))

                    wait = 5

                    killed, alive = util.kill_proc_tree(pid, timeout=wait)
                    self.error('Killing processes with SIGTERM.')
                    self.error('Processes still alive: {}'.format(len(alive)))
                    self.error('Processes killed: {}'.format(len(killed)))

                    if alive:
                        killed, alive = util.kill_proc_tree(
                            pid, timeout=wait, sig=signal.SIGKILL)
                        self.error('Killing processes with SIGKILL.')
                        self.error('Processes still alive: {}'.format(len(alive)))
                        self.error('Processes killed: {}'.format(len(killed)))

                # On success or failure log what we can
                finally:
                    with open(log_output.name) as log_input:
                        for line in log_input:
                            line = line.strip()
                            if line:
                                self.debug(line)
                    if error:
                        raise error

    def _output(self, msg, level):
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')
        msg = '({:8d}) {} {:<5}: {}'.format(os.getpid(), timestamp, level, msg)
        print(msg)
        if self.file_handle:
            self.file_handle.write(msg)
            self.file_handle.write('\n')
            self.file_handle.flush()

    def debug(self, msg):
        """Log an info message."""
        if self.log_level <= DEBUG:
            self._output(msg, 'DEBUG')

    def info(self, msg):
        """Log an info message."""
        if self.log_level <= INFO:
            self._output(msg, 'INFO')

    def error(self, msg):
        """Log an error message."""
        if self.log_level <= ERROR:
            self._output(msg, 'ERROR')

    def fatal(self, msg):
        """Log an error message and exit."""
        if self.log_level <= FATAL:
            self._output(msg, 'FATAL')
        sys.exit(1)
