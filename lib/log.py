"""Common logging functions."""

from datetime import datetime
import subprocess
import sys
import tempfile

from . import db


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
            try:
                subprocess.check_call(
                    cmd,
                    shell=True,
                    timeout=timeout,
                    stdout=log_output,
                    stderr=log_output)
            except Exception as err:  # pylint: disable=broad-except
                self.error('Exception: {}'.format(err))
                error = err
            finally:
                with open(log_output.name) as log_input:
                    for line in log_input:
                        line = line.strip()
                        if line:
                            self.debug(line)
                if error:
                    raise error

    def _output(self, msg, level):
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        msg = '{} {:<5}: {}'.format(timestamp, level, msg)
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
