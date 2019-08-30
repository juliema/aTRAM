"""Testing functions in lib/blast."""

import lib.blast as blast


# for (ends, clamp) in [('mixed_ends', ''), ('end_1', '1'),
#                       ('end_2', '2'), ('single_ends', '')]:


def test_parse_fasta_title_01():
    """It handles empty strings."""
    actual_seq_name, actual_seq_end = blast.parse_fasta_title('', '', '')
    assert actual_seq_name == ''
    assert actual_seq_end == ''


def test_parse_fasta_title_02():
    """It handles a 1 or 2 at the end of the title."""
    seq_name, seq_end = blast.parse_fasta_title(
        'title1 1 after', 'end_1', '1')
    assert seq_name == 'title1'
    assert seq_end == '1'


def test_parse_fasta_title_03():
    """It handles a 1 or 2 at the end of the title."""
    seq_name, seq_end = blast.parse_fasta_title(
        'title1/2 after', 'end_2', '2')
    assert seq_name == 'title1'
    assert seq_end == '2'


def test_parse_fasta_title_04():
    """It handles a slash delimited end."""
    seq_name, seq_end = blast.parse_fasta_title(
        'title/2 after', 'end_2', '2')
    assert seq_name == 'title'
    assert seq_end == '2'


def test_parse_fasta_title_05():
    """It handles an underscore delimited end."""
    seq_name, seq_end = blast.parse_fasta_title(
        'title_1', 'end_1', '1')
    assert seq_name == 'title'
    assert seq_end == '1'


def test_parse_fasta_title_06():
    """It handles a dot delimited end."""
    seq_name, seq_end = blast.parse_fasta_title(
        'title.1 after', 'end_1', '1')
    assert seq_name == 'title'
    assert seq_end == '1'


def test_parse_fasta_title_07():
    """It handles mixed ends with no sequence end."""
    seq_name, seq_end = blast.parse_fasta_title(
        'title', 'mixed_ends', '')
    assert seq_name == 'title'
    assert seq_end == ''


def test_parse_fasta_title_08():
    """It handles mixed ends with a delimited sequence end."""
    seq_name, seq_end = blast.parse_fasta_title(
        'title_1', 'mixed_ends', '')
    assert seq_name == 'title'
    assert seq_end == '1'


def test_parse_fasta_title_09():
    """It handles mixed ends with a space delimited sequence end."""
    seq_name, seq_end = blast.parse_fasta_title(
        'title 2 after', 'mixed_ends', '')
    assert seq_name == 'title'
    assert seq_end == '2'


def test_parse_fasta_title_10():
    """It handles single ends with no sequence end."""
    seq_name, seq_end = blast.parse_fasta_title(
        'title after', 'single_ends', '')
    assert seq_name == 'title'
    assert seq_end == ''


def test_parse_fasta_title_11():
    """It handles single ends with a delimited sequence end."""
    seq_name, seq_end = blast.parse_fasta_title(
        'title_1', 'single_ends', '')
    assert seq_name == 'title_1'
    assert seq_end == ''


def test_parse_fasta_title_12():
    """It handles single ends with a space delimited sequence end."""
    seq_name, seq_end = blast.parse_fasta_title(
        'title 2 words', 'single_ends', '')
    assert seq_name == 'title 2'
    assert seq_end == ''
