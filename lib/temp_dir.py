"""Handle creation and deletion of temporary files."""

from os import rmdir
from tempfile import mkdtemp


class TempDir:
    """Handle creation and deletion of temporary files."""

    def __init__(self, dir_name, prefix, keep=False):
        """Setup the kind of temporary directory from the arguments."""
        self.dir_name = dir_name
        self.prefix = prefix
        self.keep = keep
        self.temp_dir = None

    def __enter__(self):
        """Create the temporary directory."""
        self.temp_dir = mkdtemp(prefix=self.prefix, dir=self.dir_name)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Destroy the temporary directory."""
        if not self.keep:
            rmdir(self.temp_dir)
