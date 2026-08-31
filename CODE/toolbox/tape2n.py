"""Reader for RCN's TAPE2N, the unformatted file RCN2 takes as input.

`cmp_num.py` tokenizes text.  TAPE2N is not text: it is Fortran
sequential unformatted, one record per configuration, and each record
mixes integers, `character*6`/`character*8` strings and doubles in a
single `write` list.  Nothing generic can split that, so until this
module existed `cmp_num.py` could only answer for TAPE2N with an md5,
which says "different" and stops -- exactly the answer that stops being
useful once `PLAN_fundconst.md` Stage B deliberately moves numbers.

This module mirrors the one statement that writes the file,
`FOR/RCN36K.F` in `subroutine outpt`:

    write (2) nnn,l,wwnl,ncspvs,vprm,conf,nconf,iz,z,ion,irel,hxid,
   1  mesh,c,idb,exf,corrf,kut,npar,nlbcdu,nnlz,nkkk,ee,ee8,norb,
   2  r,ru,((pnl(i,m),i=1,kmsh),m=ns,ncspvs),iw6

It is the only `write` to unit 2 in the program, so this is the whole
format.  `FIELDS` below is that list, in that order, and is the thing to
edit if the statement ever changes -- keep the two in step, because a
reader that has silently drifted from the writer is worse than no reader
at all.

Types come from the routine's declarations, not from guesswork.  RCN is
`implicit real*8 (a-h,o-z)` with the default `implicit integer (i-n)`,
so `nnn`, `l`, `ncspvs`, `nconf`, `iz`, `ion`, `irel`, `mesh`, `idb`,
`kut`, `npar`, `nnlz`, `nkkk`, `norb` and `iw6` are `integer*4` and
everything else in the list is `real*8`, apart from the two character
arrays: `conf` is `character*6` dimensioned 3, and `nlbcdu` is
`character*8` dimensioned `ko`.  `kmsh=1801` and `ko=20` are the
`parameter` values in that routine.

The record is self-checking in a useful way: `norb` is written inside
it, before `r`, and fixes the length of the `pnl` block that follows.
The total then has to come out exactly equal to the record-length marker
gfortran wrote at the front.  `read_records` asserts that, so a layout
that has drifted from the Fortran fails loudly on the first record
rather than returning plausible nonsense.

Byte order and marker width are the native little-endian 4-byte form
that gfortran writes on every platform this package is built for.  A
record longer than 2 GB would be split into subrecords with negative
markers; no TAPE2N comes close, and `read_records` says so rather than
mis-reading one if it ever happens.
"""

import struct

#: The `parameter` values in `subroutine outpt`.
KMSH = 1801
KO = 20

#: The write list, in order.  Each entry is (name, kind, count), where
#: kind is 'i' for integer*4, 'd' for real*8, or an int for a character
#: string of that many bytes.  `count` is either an integer or the name
#: of an earlier field whose value gives the count -- `pnl` is the only
#: one of those, and it is `kmsh*norb`.
FIELDS = (
    ('nnn',    'i', KO),
    ('l',      'i', KO),
    ('wwnl',   'd', KO),
    ('ncspvs', 'i', 1),
    ('vprm',   'd', 50),
    ('conf',    6,  3),
    ('nconf',  'i', 1),
    ('iz',     'i', 1),
    ('z',      'd', 1),
    ('ion',    'i', 1),
    ('irel',   'i', 1),
    ('hxid',   'd', 1),
    ('mesh',   'i', 1),
    ('c',      'd', 1),
    ('idb',    'i', 1),
    ('exf',    'd', 1),
    ('corrf',  'd', 1),
    ('kut',    'i', 1),
    ('npar',   'i', 1),
    ('nlbcdu',  8,  KO),
    ('nnlz',   'i', KO),
    ('nkkk',   'i', KO),
    ('ee',     'd', KO),
    ('ee8',    'd', 1),
    ('norb',   'i', 1),
    ('r',      'd', KMSH),
    ('ru',     'd', KMSH),
    ('pnl',    'd', 'kmsh*norb'),
    ('iw6',    'i', 1),
)

#: The fields that carry numbers a tolerance comparison applies to.
#: Integers and strings are structure: they must match exactly, and a
#: difference in one of them is not a rounding shift.
REAL_FIELDS = tuple(n for n, k, _ in FIELDS if k == 'd')
INT_FIELDS = tuple(n for n, k, _ in FIELDS if k == 'i')
CHAR_FIELDS = tuple(n for n, k, _ in FIELDS if isinstance(k, int))


class Tape2nError(Exception):
    """The file is not a TAPE2N this reader understands."""


def _unpack(kind, count, blob, pos):
    width = 4 if kind == 'i' else 8 if kind == 'd' else kind
    size = width * count
    if len(blob) - pos < size:
        raise Tape2nError('record ends inside a field: %d bytes wanted at '
                          'offset %d, %d left' % (size, pos, len(blob) - pos))
    if kind == 'i':
        values = list(struct.unpack('<%di' % count, blob[pos:pos + size]))
    elif kind == 'd':
        values = list(struct.unpack('<%dd' % count, blob[pos:pos + size]))
    else:                                   # character*kind, `count` of them
        chunk = blob[pos:pos + size]
        values = [chunk[i * kind:(i + 1) * kind].decode('latin-1')
                  for i in range(count)]
    return values, pos + size


def read_records(path):
    """Read every record of `path` and return them as a list of dicts.

    A scalar field comes back as a scalar, an array field as a list.
    `pnl` comes back flat, `kmsh*norb` long, in the order the implied-do
    writes it: all of `i` for one `m`, then the next `m`.
    """
    records = []
    with open(path, 'rb') as handle:
        blob = handle.read()

    pos = 0
    total = len(blob)
    while pos < total:
        if total - pos < 8:
            raise Tape2nError('trailing bytes at %d: not a whole record' % pos)
        (length,) = struct.unpack('<i', blob[pos:pos + 4])
        if length < 0:
            raise Tape2nError('record at %d is a subrecord (length marker '
                              '%d); this reader does not join subrecords, '
                              'and no TAPE2N should be large enough to need '
                              'them' % (pos, length))
        body = blob[pos + 4:pos + 4 + length]
        if len(body) != length:
            raise Tape2nError('record at %d claims %d bytes and the file has '
                              '%d' % (pos, length, len(body)))
        (tail,) = struct.unpack('<i', blob[pos + 4 + length:pos + 8 + length])
        if tail != length:
            raise Tape2nError('record at %d: head marker %d, tail marker %d'
                              % (pos, length, tail))

        record = {}
        at = 0
        for name, kind, count in FIELDS:
            if count == 'kmsh*norb':
                count = KMSH * record['norb']
                if count < 0:
                    raise Tape2nError('record at %d has norb=%d'
                                      % (pos, record['norb']))
            values, at = _unpack(kind, count, body, at)
            record[name] = values[0] if count == 1 else values
        if at != length:
            raise Tape2nError('record at %d: layout consumed %d bytes of %d. '
                              'FIELDS no longer matches the write statement '
                              'in RCN36K.F.' % (pos, at, length))
        records.append(record)
        pos += 8 + length

    if not records:
        raise Tape2nError('%s holds no records' % path)
    return records


def flatten(record):
    """Yield (label, value) for every real*8 number in one record.

    The label is what a report quotes, so it carries the index: `r[17]`,
    `pnl[3421]`.  Scalars carry no index.
    """
    for name in REAL_FIELDS:
        value = record[name]
        if isinstance(value, list):
            for index, item in enumerate(value):
                yield '%s[%d]' % (name, index), item
        else:
            yield name, value


def structure(record):
    """The part of a record that must match exactly: integers and text.

    Returned as a tuple so that two records can be compared with `==`.
    `norb` is in here as well as being the length of `pnl`, which is the
    point: if it moved, the numeric comparison of `pnl` is meaningless
    and the structural check has to fail first.
    """
    return tuple((name, tuple(record[name])
                  if isinstance(record[name], list) else record[name])
                 for name in INT_FIELDS + CHAR_FIELDS)
