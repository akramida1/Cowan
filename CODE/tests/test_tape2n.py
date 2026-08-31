"""toolbox/tape2n.py: the record layout of RCN's unformatted TAPE2N.

The reader mirrors one `write` statement in `FOR/RCN36K.F`, and the way
it can fail is by drifting from it -- a field added, reordered, or given
a different type -- while still returning numbers that look plausible.
So the tests here are mostly about the self-check that makes drift loud:
the record carries its own length marker, and the layout has to consume
exactly that many bytes.

They run against the Sn VIII TAPE2N in `WORK/`, which ships as inputs
only, so they skip rather than fail when the case has not been run.
"""

import struct

import pytest

from conftest import WORK
from toolbox import cmpnum, tape2n


def tape():
    path = WORK / 'TAPE2N'
    if not path.is_file():
        pytest.skip('WORK/TAPE2N is not present: run the Sn VIII case first')
    return path


# ------------------------------------------------------------- the layout

def test_the_layout_consumes_exactly_the_record_gfortran_wrote():
    """The one check that catches drift from the write statement.

    gfortran brackets each record with its length; `read_records` raises
    if the fields in `FIELDS` do not add up to it.  Reading the real file
    without raising *is* the assertion -- but assert the count too, so
    that a reader which somehow returned nothing cannot pass.
    """
    records = tape2n.read_records(tape())
    assert records


def test_a_record_describes_the_atom_the_deck_asked_for():
    """Sn VIII, 4d7.  If a field had shifted by so much as one integer,
    these would come back as noise rather than as the deck's own values."""
    first = tape2n.read_records(tape())[0]
    assert first['iz'] == 50
    assert first['ion'] == 7
    assert ''.join(first['conf']).split() == ['Sn7+', '4d7']
    assert 1 <= first['ncspvs'] <= tape2n.KO
    assert 1 <= first['norb'] <= first['ncspvs']
    assert first['nlbcdu'][0].strip() == '1S'


def test_the_orbital_labels_agree_with_the_quantum_numbers():
    """`nnlz` is `100*n+l` and `nlbcdu` is the same orbital spelled out.

    Two independent fields at opposite ends of the record: they can only
    agree if everything between them is being read at the right offset.
    """
    first = tape2n.read_records(tape())[0]
    letters = 'SPDFGHIK'
    for m in range(first['ncspvs']):
        n, l = divmod(first['nnlz'][m], 100)
        assert first['nnn'][m] == n
        assert first['l'][m] == l
        assert first['nlbcdu'][m].strip() == '%d%s' % (n, letters[l])


def test_the_radial_mesh_is_monotonic_and_positive():
    """`r` is the radial mesh.  A field read at the wrong offset would
    put integers or text into it, and it would not come back ordered.

    Only the first `mesh` points of the `kmsh`-long array are the mesh;
    the rest are zero, and that is the array's own statement of where it
    ends.
    """
    first = tape2n.read_records(tape())[0]
    assert len(first['r']) == tape2n.KMSH
    used = first['r'][:first['mesh']]
    assert used[0] > 0.0
    assert all(b > a for a, b in zip(used, used[1:]))
    assert all(value == 0.0 for value in first['r'][first['mesh']:])


def test_pnl_is_kmsh_long_for_every_orbital_written():
    first = tape2n.read_records(tape())[0]
    assert len(first['pnl']) == tape2n.KMSH * first['norb']


def test_flatten_covers_every_real_and_nothing_else():
    first = tape2n.read_records(tape())[0]
    labels = [label for label, _ in tape2n.flatten(first)]
    assert len(labels) == len(set(labels))
    expected = sum(len(first[name]) if isinstance(first[name], list) else 1
                   for name in tape2n.REAL_FIELDS)
    assert len(labels) == expected
    assert all(isinstance(value, float) for _, value in tape2n.flatten(first))


def test_structure_holds_the_integers_and_the_text_and_no_reals():
    first = tape2n.read_records(tape())[0]
    names = [name for name, _ in tape2n.structure(first)]
    assert set(names) == set(tape2n.INT_FIELDS) | set(tape2n.CHAR_FIELDS)
    assert not set(names) & set(tape2n.REAL_FIELDS)


# ------------------------------------------------------- refusing bad files

def test_a_truncated_record_is_refused_rather_than_guessed_at(tmp_path):
    blob = tape().read_bytes()
    broken = tmp_path / 'TAPE2N'
    broken.write_bytes(blob[:len(blob) // 2])
    with pytest.raises(tape2n.Tape2nError):
        tape2n.read_records(broken)


def test_a_record_whose_markers_disagree_is_refused(tmp_path):
    """The head and tail markers must match.  They are what makes the
    length claim trustworthy in the first place."""
    blob = bytearray(tape().read_bytes())
    (length,) = struct.unpack('<i', bytes(blob[:4]))
    blob[4 + length:8 + length] = struct.pack('<i', length - 8)
    broken = tmp_path / 'TAPE2N'
    broken.write_bytes(bytes(blob))
    with pytest.raises(tape2n.Tape2nError):
        tape2n.read_records(broken)


def test_a_layout_that_no_longer_matches_the_writer_is_refused(monkeypatch):
    """Drop a field from FIELDS and the byte count stops adding up.

    This is the failure the whole design is built around: a reader that
    has drifted from `RCN36K.F` must fail loudly rather than return
    numbers read at the wrong offset.
    """
    monkeypatch.setattr(tape2n, 'FIELDS', tape2n.FIELDS[:-1])
    with pytest.raises(tape2n.Tape2nError) as raised:
        tape2n.read_records(tape())
    assert 'FIELDS' in str(raised.value)


def test_a_file_that_is_not_a_tape2n_is_refused(tmp_path):
    junk = tmp_path / 'TAPE2N'
    junk.write_bytes(b'not a fortran record at all')
    with pytest.raises(tape2n.Tape2nError):
        tape2n.read_records(junk)


# --------------------------------------------------- through the comparator

def test_cmp_num_reads_tape2n_by_layout_instead_of_by_md5():
    """A TAPE2N pair now gets the same treatment as a text file: a count
    of numbers compared and a deviation, not just "different"."""
    path = tape()
    report = cmpnum.compare_files(path, path)
    assert report.binary and report.layout == 'TAPE2N'
    assert report.records and report.compared
    assert not report.structural and not report.offenders and report.ok


def test_a_changed_real_in_tape2n_is_numeric_not_structural(tmp_path):
    """Move one `ee` by a part in 1e6 -- what a constants change does --
    and it must be reported as a deviation, with the record intact."""
    blob = bytearray(tape().read_bytes())
    offset = 4 + 4 * tape2n.KO * 2          # past nnn and l, into wwnl
    (value,) = struct.unpack('<d', bytes(blob[offset:offset + 8]))
    blob[offset:offset + 8] = struct.pack('<d', value * (1.0 + 1e-6))
    moved = tmp_path / 'TAPE2N'
    moved.write_bytes(bytes(blob))
    report = cmpnum.compare_files(tape(), moved)
    assert not report.structural
    assert len(report.offenders) == 1
    assert report.offenders[0][6].startswith('wwnl[')
    assert 5e-7 < report.largest < 2e-6


def test_a_changed_integer_in_tape2n_is_structural_not_numeric(tmp_path):
    """`iz` is the atomic number.  No tolerance applies to it, and a
    file whose integers moved is not the same calculation."""
    blob = bytearray(tape().read_bytes())
    records = tape2n.read_records(tape())
    where = blob.find(struct.pack('<i', records[0]['iz']),
                      4 + 4 * tape2n.KO * 2 + 8 * tape2n.KO)
    assert where > 0
    blob[where:where + 4] = struct.pack('<i', records[0]['iz'] + 1)
    moved = tmp_path / 'TAPE2N'
    moved.write_bytes(bytes(blob))
    report = cmpnum.compare_files(tape(), moved)
    assert report.structural and not report.ok


def test_a_tape2n_that_will_not_parse_is_structural_not_a_silent_md5(tmp_path):
    """Falling back to md5 here would report "not quantifiable" for what
    is really a broken reader.  It has to say which."""
    good = tmp_path / 'a'
    good.mkdir()
    bad = tmp_path / 'b'
    bad.mkdir()
    (good / 'TAPE2N').write_bytes(tape().read_bytes())
    (bad / 'TAPE2N').write_bytes(b'\0' * 64)
    report = cmpnum.compare_files(good / 'TAPE2N', bad / 'TAPE2N')
    assert report.layout == 'TAPE2N'
    assert report.structural and not report.ok
    assert 'RCN36K.F' in report.structural[0][1]
