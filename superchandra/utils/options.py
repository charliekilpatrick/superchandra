"""
Command-line option definitions for the superchandra CLI.

Defines positional ra/dec and optional flags (--model, --download, --before,
--after, --extinction, --redshift, --radius, etc.) used by the pipeline
entry point.
"""
import argparse


def get_parser(parser=None, usage=None):
    """
    Return an ArgumentParser with all superchandra CLI options.

    If an existing parser is passed, options are added to it (e.g. for
    embedding in another CLI). Otherwise a new parser is created.

    Parameters
    ----------
    parser : argparse.ArgumentParser, optional
        Parser to extend; if None, a new one is created.
    usage : str, optional
        Usage string for the new parser (e.g. "superchandra ra dec").

    Returns
    -------
    argparse.ArgumentParser
        Parser with positional ra, dec and all superchandra options added.
    """
    if parser is None:
        parser = argparse.ArgumentParser(usage=usage, conflict_handler="resolve")
    parser.add_argument("ra", default=None, type=str, help="RA")
    parser.add_argument("dec", default=None, type=str, help="DEC")
    parser.add_argument(
        "--makeclean",
        default=False,
        action="store_true",
        help="Clean up all output files from previous runs then exit.",
    )
    parser.add_argument(
        "--download",
        default=False,
        action="store_true",
        help="Download the raw data files given input ra and dec.",
    )
    parser.add_argument(
        "--before",
        default=None,
        type=str,
        metavar="YYYY-MM-DD",
        help="Date after which we should reject all Chandra/ACIS observations for reduction.",
    )
    parser.add_argument(
        "--after",
        default=None,
        type=str,
        metavar="YYYY-MM-DD",
        help="Date before which we should reject all Chandra/ACIS observations for reduction.",
    )
    parser.add_argument(
        "--clobber",
        default=False,
        action="store_true",
        help="Overwrite files when using download mode.",
    )
    parser.add_argument(
        "--extinction",
        default=0.0,
        type=float,
        help="Host extinction (Av) for the observed SN. MW-like dust is assumed (Güver & Özel 2009).",
    )
    parser.add_argument(
        "--redshift",
        "-z",
        default=0.0,
        type=float,
        help="Host redshift for evaluating extinction.",
    )
    parser.add_argument(
        "--quick",
        "-q",
        default=False,
        action="store_true",
        help="Only calculate the first model in a given set.",
    )
    parser.add_argument(
        "--model",
        "-m",
        default="sss",
        type=str,
        help="Model spectrum to use (options are sss|bh|frb).",
    )
    parser.add_argument(
        "--search_radius",
        "-s",
        default=0.1667,
        type=float,
        help="Override search radius with input value (units are degrees).",
    )
    parser.add_argument(
        "--radius",
        "-r",
        default=None,
        type=float,
        help="Radius (in arcsec) for performing aperture photometry.",
    )
    return parser
