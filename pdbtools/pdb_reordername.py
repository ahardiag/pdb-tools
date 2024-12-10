#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2018 João Pedro Rodrigues
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Reorder atoms based on the first occurrence of residue atom names in the PDB file.

Residues names are matched *without* taking into consideration spaces.

Combine with pdb_reatom and pdb_sort to rearrange the records and renumber the atom after lines have been exchanged.

Usage:
    python pdb_reordername.py -<option> <pdb file>

Example:
    python pdb_reordername.py -ALA 1CTF.pdb

This program is part of the `pdb-tools` suite of utilities and should not be
distributed isolatedly. The `pdb-tools` were created to quickly manipulate PDB
files using the terminal, and can be used sequentially, with one tool streaming
data to another. They are based on old FORTRAN77 code that was taking too much
effort to maintain and compile. RIP.
"""

import os
import sys

__author__ = ["Joao Rodrigues", "Joao M.C. Teixeira"]
__email__ = ["j.p.g.l.m.rodrigues@gmail.com", "joaomcteixeira@gmail.com"]


def check_input(args):
    """Checks whether to read from stdin/file and validates user input/options.
    """

    # Defaults
    option = ''
    fh = sys.stdin  # file handle

    if not len(args):
        # Reading from pipe with default option
        if sys.stdin.isatty():
            sys.stderr.write(__doc__)
            sys.exit(1)

    elif len(args) == 1:
        # One of two options: option & Pipe OR file & default option
        if args[0].startswith('-'):
            option = args[0][1:]
            if sys.stdin.isatty():  # ensure the PDB data is streamed in
                emsg = 'ERROR!! No data to process!\n'
                sys.stderr.write(emsg)
                sys.stderr.write(__doc__)
                sys.exit(1)

        else:
            if not os.path.isfile(args[0]):
                emsg = 'ERROR!! File not found or not readable: \'{}\'\n'
                sys.stderr.write(emsg.format(args[0]))
                sys.stderr.write(__doc__)
                sys.exit(1)

            fh = open(args[0], 'r')

    elif len(args) == 2:
        # Two options: option & File
        if not args[0].startswith('-'):
            emsg = 'ERROR! First argument is not an option: \'{}\'\n'
            sys.stderr.write(emsg.format(args[0]))
            sys.stderr.write(__doc__)
            sys.exit(1)

        if not os.path.isfile(args[1]):
            emsg = 'ERROR!! File not found or not readable: \'{}\'\n'
            sys.stderr.write(emsg.format(args[1]))
            sys.stderr.write(__doc__)
            sys.exit(1)

        option = args[0][1:]
        fh = open(args[1], 'r')

    else:  # Whatever ...
        sys.stderr.write(__doc__)
        sys.exit(1)

    # Validate option
    option_set = set([o.upper().strip() for o in option.split(',') if o.strip()])
    if not option_set:
        sys.stderr.write('ERROR!! Residue name set cannot be empty\n')
        sys.stderr.write(__doc__)
        sys.exit(1)
    else:
        for elem in option_set:
            if len(elem) > 3:
                emsg = 'ERROR!! Residue name is invalid: \'{}\'\n'
                sys.stderr.write(emsg.format(elem))
                sys.stderr.write(__doc__)
                sys.exit(1)

    return (fh, option_set)


def run(fhandle, resname_set):
    """
    Reorder lines for residues with non-unique atom name order.
    Maintains all non-coordinate lines (e.g., CRYST1, REMARK, etc.).
    Brute-force algorithm.

    Parameters
    ----------
    fhandle : iterable
        A line-by-line iterator of the original PDB file.
    resname_set : set, list, or tuple
        The name of the residues with non-unique atom name order.

    Yields
    ------
    str
        The reordered PDB lines.
    """
    records = ('ATOM', 'HETATM', 'ANISOU', 'TER')

    ref_atom_order = {}
    residue_lines = {}
    not_residues_lines = []
    buffer_line = None  # To hold the line when the atom names repeat

    # Parse the reference atom order and yield non-record lines
    for line in fhandle:
        if line.startswith(records):
            resname = line[17:20].strip()
            if resname in resname_set:
                atom_name = line[12:16].strip()
                if resname not in ref_atom_order:
                    ref_atom_order[resname] = []
                if atom_name not in ref_atom_order[resname]:
                    ref_atom_order[resname].append(atom_name)
                    yield line
                else:
                    # Stop after parsing the first residue and move to reordering
                    buffer_line = line
                    break
        else:
            # Yield non-coordinate lines (e.g., CRYST1, REMARK)
            yield line

    # Continue processing from the line after the break
    if buffer_line:
        # Process the line where the break occurred
        if buffer_line.startswith(records):
            resname = buffer_line[17:20].strip()
            resid = buffer_line[22:26].strip()
            if resname in resname_set:
                if resid not in residue_lines:
                    residue_lines[resid] = []
                residue_lines[resid].append(buffer_line)

    # Continue reading from fhandle without converting to a list
    for line in fhandle:
        if line.startswith(records):
            resname = line[17:20].strip()
            resid = line[22:26].strip()
            if resname in resname_set:
                if resid not in residue_lines:
                    residue_lines[resid] = []
                residue_lines[resid].append(line)
            else:
                # Collect lines not matching resname_set
                not_residues_lines.append(line)
        else:
            # Collect remaining non-coordinate lines
            not_residues_lines.append(line)

    # Reorder and yield the grouped residue lines
    for resid, lines in residue_lines.items():
        resname = lines[0][17:20].strip()
        if resname in ref_atom_order:
            atom_order = ref_atom_order[resname]
            sorted_lines = []
            for atom_name in atom_order:
                sorted_lines.extend(
                    [line for line in lines if line[12:16].strip() == atom_name]
                )
            for line in sorted_lines:
                yield line
        else:
            for line in lines:
                yield line

    # Yield any collected non-residues lines
    for line in not_residues_lines:
        yield line


filter_residue_by_name = run


def main():
    # Check Input
    pdbfh, resname_set = check_input(sys.argv[1:])

    # Do the job
    new_pdb = run(pdbfh, resname_set)

    try:
        _buffer = []
        _buffer_size = 5000  # write N lines at a time
        for lineno, line in enumerate(new_pdb):
            if not (lineno % _buffer_size):
                sys.stdout.write(''.join(_buffer))
                _buffer = []
            _buffer.append(line)

        sys.stdout.write(''.join(_buffer))
        sys.stdout.flush()
    except IOError:
        # This is here to catch Broken Pipes
        # for example to use 'head' or 'tail' without
        # the error message showing up
        pass

    # last line of the script
    # We can close it even if it is sys.stdin
    pdbfh.close()
    sys.exit(0)


if __name__ == '__main__':
    main()
