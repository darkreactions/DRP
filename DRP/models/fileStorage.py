"""Implements a special class of FileSystemStorage which keeps the file name the same."""

from django.core.files.storage import FileSystemStorage
import os

class OverwriteStorage(FileSystemStorage):
    """Overwrites a file if it already exists, rather than trying to uniqueise the new filename."""

    def get_available_name(self, name, max_length=None): 
        if max_length and len(name) > max_length:
            raise ValueError('File name too long for secure storage.')
        if self.exists(name):
            os.remove(os.path.join(self.location, name))
        return name
