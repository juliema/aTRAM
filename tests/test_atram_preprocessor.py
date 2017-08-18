"""Testing functions in atram_preprocessor."""

# pylama: ignore=D103

from os.path import join
import atram_preprocessor
import lib.db as db
import lib.blast as blast
import tests.mock as mock


def test_load_seq(monkeypatch):
    mock.history = []
    db.BATCH_SIZE = 5

    mock.mock(monkeypatch, db, 'insert_sequences_batch')

    file_1 = join('tests', 'data', 'load_seq1.txt')
    file_2 = join('tests', 'data', 'load_seq2.txt')
    atram_preprocessor.load_seqs('connection', [file_1, file_2])

    assert mock.history == [
        {'module': 'lib.db', 'func': 'insert_sequences_batch',
         'db_conn': 'connection',
         'batch': [
             ('seq1', '1', 'AAAAAAAAAA'),
             ('seq2', '1', 'AAAAAAAAAAGGGGGGGGGG'),
             ('seq3', '1', 'AAAAAAAAAA'),
             ('seq4', '1', 'AAAAAAAAAA'),
             ('seq5/3', '', 'AAAAAAAAAAGGGGGGGGGG')]},
        {'module': 'lib.db', 'func': 'insert_sequences_batch',
         'db_conn': 'connection',
         'batch': [
            ('seq1', '2', 'AAAAAAAAAA'),
            ('seq2', '2', 'AAAAAAAAAAGGGGGGGGGG'),
            ('seq3', '2', 'AAAAAAAAAA'),
            ('seq4', '2', 'AAAAAAAAAAGGGGGGGGGG')]},
        {'module': 'lib.db', 'func': 'insert_sequences_batch',
         'db_conn': 'connection',
         'batch': [
            ('seq6', '1', 'TTTTTTTTTT'),
            ('seq7', '1', 'TTTTTTTTTTCCCCCCCCCC'),
            ('seq8/a.1 suffix', '', 'TTTTTTTTTT'),
            ('seq8', '2', 'TTTTTTTTTTCCCCCCCCCC')]}]


def test_assign_seqs_to_shards(monkeypatch):
    mock.mock(monkeypatch, db, 'get_two_sequences', returns=[('seq1', 'seq2'),
                                                             ('seq3', 'seq3')])
    mock.mock(monkeypatch, db, 'get_sequence_count', returns=100)

    shard_list = atram_preprocessor.assign_seqs_to_shards(True, 3)

    # A list of pairs of LIMIT and OFFSET pairs for queries that build shards
    # The values will depend on what is returned from get_two_sequences
    assert shard_list == [(34, 0), (32, 34), (34, 66)]


def test_create_one_blast_shard(monkeypatch):
    mock.history = []
    mock.mock(monkeypatch, blast, 'create_db')
    mock.mock(monkeypatch, blast, 'shard_path', returns='shard/path')
    mock.mock(monkeypatch, atram_preprocessor, 'fill_blast_fasta')

    args = {'blast_db': 'my_blast_db', 'temp_dir': 'my_temp_dir'}
    shard_params = ['limit', 'offset']

    atram_preprocessor.create_one_blast_shard(args, shard_params, 11)

    history = mock.filter('lib.blast', 'shard_path')
    assert history == [{'blast_db': 'my_blast_db',
                        'shard_index': 11}]

    history = mock.filter('atram_preprocessor', 'fill_blast_fasta')
    assert history == [{'blast_db': 'my_blast_db',
                        'fasta_path': 'my_temp_dir/pyt_011.fasta',
                        'shard_params': ['limit', 'offset']}]

    history = mock.filter('lib.blast', 'create_db')
    assert history == [{'fasta_file': 'my_temp_dir/pyt_011.fasta',
                        'shard_path': 'shard/path',
                        'temp_dir': 'my_temp_dir'}]
