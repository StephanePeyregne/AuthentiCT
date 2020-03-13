__all__ = ['CIGAR_REGEX']

import re

CIGAR_REGEX = re.compile(r"(\d+)([MIDNSHPX=])", flags=re.ASCII)
