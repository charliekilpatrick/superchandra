"""
File I/O helpers for the superchandra pipeline.

Provides decompression and similar utilities (e.g. gunzip in place).
"""
import os


def gunzip(filepath):
    """
    Decompress a gzipped file in place and return the path without .gz.

    Runs the system ``gunzip -f`` command; the original .gz file is removed.

    Parameters
    ----------
    filepath : str
        Path to a .gz file.

    Returns
    -------
    str
        Path to the decompressed file (same basename without .gz extension).

    Notes
    -----
    Requires the ``gunzip`` command to be available on the system.
    """
    newfile = filepath.replace('.gz', '')
    os.system('gunzip -f ' + filepath)
    return newfile
