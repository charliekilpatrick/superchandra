"""
Centralized logging for the superchandra package.

Provides get_logger() and configure_logging() so all pipeline output uses
standard logging (info, warning, error) instead of print. The package root
logger name is ``superchandra``; child loggers are ``superchandra.module``.
"""
import logging
import sys

# Root logger name for the package; child loggers use superchandra.module.
PACKAGE_LOG_NAME = "superchandra"


def get_logger(name=None):
    """
    Return a logger for the superchandra package.

    Parameters
    ----------
    name : str, optional
        If None, returns the package root logger ("superchandra").
        If provided, returns a child logger "superchandra.name" (e.g. "superchandra.superchandra").

    Returns
    -------
    logging.Logger
        Logger instance.
    """
    if name is None:
        return logging.getLogger(PACKAGE_LOG_NAME)
    if name.startswith(PACKAGE_LOG_NAME + "."):
        return logging.getLogger(name)
    return logging.getLogger(PACKAGE_LOG_NAME + "." + name)


def configure_logging(level=logging.INFO, stream=None):
    """
    Configure the package root logger with a stream handler and level.

    Idempotent: if the root package logger already has a handler, this does
    not add another. Clear handlers on the logger first if you need to replace.

    Parameters
    ----------
    level : int, optional
        Logging level (e.g. logging.INFO, logging.DEBUG). Default is INFO.
    stream : file-like, optional
        Stream for the handler (e.g. sys.stdout, sys.stderr).
        Default is sys.stderr.

    Returns
    -------
    logging.Logger
        The package root logger (superchandra).
    """
    if stream is None:
        stream = sys.stderr
    log = logging.getLogger(PACKAGE_LOG_NAME)
    log.setLevel(level)
    # Avoid duplicate handlers when used as a library.
    if not log.handlers:
        handler = logging.StreamHandler(stream)
        handler.setLevel(level)
        handler.setFormatter(logging.Formatter("%(message)s"))
        log.addHandler(handler)
    return log
