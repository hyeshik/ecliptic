__all__ = ['TemporaryDirectory', 'DeleteOnError']

import tempfile
import shutil
import os

class TemporaryDirectory(object):
    def __init__(self, dir='.'):
        self.dir = dir
        self.path = None

    def __enter__(self):
        self.path = tempfile.mkdtemp(dir=self.dir)
        return self.path

    def __exit__(self, type, value, traceback):
        if self.path is not None:
            shutil.rmtree(self.path)

class DeleteOnError(object):
    def __init__(self, path, opener=None):
        self.path = path
        self.opener = opener

    def __enter__(self):
        if self.opener is None:
            return open(self.path, 'wb')
        else:
            return self.opener(self.path, 'wb')

    def __exit__(self, type, value, traceback):
        if type is not None and os.path.exists(self.path):
            os.unlink(self.path)

