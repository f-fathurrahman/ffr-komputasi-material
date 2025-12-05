# -*- coding: utf-8 -*
#    my_boltztrap2, a program for interpolating band structures and calculating
#    semi-classical transport coefficients.
#    Copyright (C) 2017-2025 Georg K. H. Madsen <georg.madsen@tuwien.ac.at>
#    Copyright (C) 2017-2025 Jesús Carrete <jesus.carrete.montana@tuwien.ac.at>
#    Copyright (C) 2017-2025 Matthieu J. Verstraete <matthieu.verstraete@ulg.ac.be>
#    Copyright (C) 2018-2019 Genadi Naydenov <gan503@york.ac.uk>
#    Copyright (C) 2020 Gavin Woolman <gwoolma2@staffmail.ed.ac.uk>
#    Copyright (C) 2020 Roman Kempt <roman.kempt@tu-dresden.de>
#    Copyright (C) 2022 Robert Stanton <stantor@clarkson.edu>
#    Copyright (C) 2024 Haoyu (Daniel) Yang <yanghaoyu97@outlook.com>
#
#    This file is part of my_boltztrap2.
#
#    my_boltztrap2 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    my_boltztrap2 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with my_boltztrap2. If not, see <http://www.gnu.org/licenses/>.

import contextlib
import io
import logging
import os
import os.path
import sys
import time

# Miscellaneous code used in several parts of the package.

logger = logging.getLogger(__name__)


def ffloat(arg):
    """Try to parse a string from Fortran code as a floating-point number."""
    return float(arg.lower().replace("d", "e"))


def info(*args):
    """Print something to a ioStringIO and log its contents as info."""
    with io.StringIO("") as message:
        print(*args, file=message)
        for l in message.getvalue().splitlines():
            logger.info(l)


def warning(*args):
    """Print something to a ioStringIO and log its contents as a warning."""
    with io.StringIO("") as message:
        print(*args, file=message)
        for l in message.getvalue().splitlines():
            logger.warning(l)


def lexit(message, code=1):
    """Log an error message and exit."""
    logger.error(message)
    sys.exit(code)


@contextlib.contextmanager
def dir_context(dirname):
    """Create a context to run code in a different directory."""
    cwd = os.getcwd()
    try:
        os.chdir(dirname)
        yield
    finally:
        os.chdir(cwd)


class TimerContext:
    """Simple stopwatch-line context manager class."""

    def __enter__(self):
        """Start counting the time."""
        self.t0 = time.time()
        return self

    def __exit__(self, *args):
        """Do nothing."""
        pass

    def get_deltat(self):
        """Return the time elapsed since the object was created."""
        return time.time() - self.t0
