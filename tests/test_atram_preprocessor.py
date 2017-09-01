"""Testing functions in atram_preprocessor."""

# pylama: ignore=D103

from os.path import join
import tempfile
import atram_preprocessor
import lib.db as db
import lib.log as log
import lib.blast as blast
import tests.mock as mock


def test_preprocess():
    args = {
        'log_file': 'log_file_arg',
        'blast_db': 'blast_db_arg',
        'sra_files': 'sra_files_arg',
        'shard_count': 'shard_count_arg'
    }

    mock.it(log, 'setup')
    mock.context(db, 'connect', 'my_connection')
    mock.it(db, 'create_metadata_table')
    mock.it(db, 'create_sequences_table')
    mock.it(atram_preprocessor, 'load_seqs')
    mock.it(log, 'info')
    mock.it(db, 'create_sequences_index')
    mock.it(atram_preprocessor, 'assign_seqs_to_shards', 'shard_list')
    mock.it(atram_preprocessor, 'create_all_blast_shards')

    atram_preprocessor.preprocess(args)

    expect = [
        {
            'module': 'lib.log',
            'func': 'setup',
            'blast_db': 'blast_db_arg',
            'log_file': 'log_file_arg',
        }, {
            'module': 'lib.db',
            'blast_db': 'blast_db_arg',
            'func': 'connect',
        }, {
            'module': 'lib.db',
            'func': 'create_metadata_table',
            'db_conn': 'my_connection',
        }, {
            'module': 'lib.db',
            'func': 'create_sequences_table',
            'db_conn': 'my_connection',
        }, {
            'module': 'atram_preprocessor',
            'func': 'load_seqs',
            'db_conn': 'my_connection',
            'sra_files': 'sra_files_arg',
        }, {
            'module': 'lib.log',
            'func': 'info',
            'msg': 'Creating an index for the sequence table',
        }, {
            'db_conn': 'my_connection',
            'func': 'create_sequences_index',
            'module': 'lib.db'
        }, {
            'module': 'atram_preprocessor',
            'db_conn': 'my_connection',
            'func': 'assign_seqs_to_shards',
            'shard_count': 'shard_count_arg',
        }, {
            'module': 'atram_preprocessor',
            'func': 'create_all_blast_shards',
            'shard_list': 'shard_list',
            'args': {
                'blast_db': 'blast_db_arg',
                'log_file': 'log_file_arg',
                'shard_count': 'shard_count_arg',
                'sra_files': 'sra_files_arg'},
        }]
    assert expect == mock.history


def test_load_seqs():
    db.BATCH_SIZE = 5

    mock.it(log, 'info')
    mock.it(db, 'insert_sequences_batch')

    file_1 = join('tests', 'data', 'load_seq1.txt')
    file_2 = join('tests', 'data', 'load_seq2.txt')
    atram_preprocessor.load_seqs('connection', [file_1, file_2])

    expect = [
        {'module': 'lib.log',
         'func': 'info',
         'msg': 'Loading "tests/data/load_seq1.txt" into sqlite database'},
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
        {'module': 'lib.log',
         'func': 'info',
         'msg': 'Loading "tests/data/load_seq2.txt" into sqlite database'},
        {'module': 'lib.db', 'func': 'insert_sequences_batch',
         'db_conn': 'connection',
         'batch': [
            ('seq6', '1', 'TTTTTTTTTT'),
            ('seq7', '1', 'TTTTTTTTTTCCCCCCCCCC'),
            ('seq8/a.1 suffix', '', 'TTTTTTTTTT'),
            ('seq8', '2', 'TTTTTTTTTTCCCCCCCCCC')]}]
    assert expect == mock.history


def test_assign_seqs_to_shards():
    mock.it(log, 'info')
    mock.it(db, 'get_sequence_count', returns=100)
    mock.it(db, 'get_shard_cut_pair', returns=[('seq1', 'seq2'),
                                               ('seq3', 'seq3')])

    # We are breaking the database into three shards. The endpoints and lengths
    # of the shards will be adjusted based upon where pairs of sequences call.
    shard_list = atram_preprocessor.assign_seqs_to_shards(True, 3)

    # A list of pairs of (LIMIT, OFFSET) for queries that build shards
    # The LIMIT/lengths are calculated from current & previous offsets
    # The OFFSET will depend on what is returned from get_shard_cut_pair
    #   1) The first pair always has offset 0
    #   2) Because the sequences for the first pair are different the offset
    #      will be pushed forward one from 33 to 34
    #   3) Because the sequences are the same the offset will stay at 66.

    assert [(34, 0), (32, 34), (34, 66)] == shard_list

    expect = [{'msg': 'Assigning sequences to shards'}]
    assert expect == mock.filter('lib.log', 'info')


def test_create_one_blast_shard():
    mock.it(blast, 'create_db')
    mock.it(blast, 'shard_path', returns='shard/path')
    mock.it(atram_preprocessor, 'fill_blast_fasta')

    args = {'blast_db': 'my_blast_db', 'temp_dir': 'my_temp_dir'}
    shard_params = ['limit', 'offset']

    atram_preprocessor.create_one_blast_shard(args, shard_params, 11)

    expect = [{'blast_db': 'my_blast_db', 'shard_index': 11}]
    assert expect == mock.filter('lib.blast', 'shard_path')

    expect = [{'blast_db': 'my_blast_db',
               'fasta_path': 'my_temp_dir/pyt_011.fasta',
               'shard_params': ['limit', 'offset']}]
    assert expect == mock.filter('atram_preprocessor', 'fill_blast_fasta')

    expect = [{'fasta_file': 'my_temp_dir/pyt_011.fasta',
               'shard_path': 'shard/path',
               'temp_dir': 'my_temp_dir'}]
    assert expect == mock.filter('lib.blast', 'create_db')


def test_fill_blast_fasta():
    mock.context(db, 'connect', 'my_connection')

    # We want to return an array of tuples as one item
    mock.it(db, 'get_sequences_in_shard', [[
        ('seq1', '1', 'AAAAAAAAAA'),
        ('seq1', '2', 'CCCCCCCCCC'),
        ('seq2', '', 'GGGGGGGGGG'),
        ('seq3', '1', 'TTTTTTTTTT')]])

    with tempfile.TemporaryDirectory(prefix='test_') as temp_dir:
        blast_db = 'test_blast_db'
        fasta_path = join(temp_dir, 'test_output.fasta')
        shard_params = (100, 200)  # limit and offset

        atram_preprocessor.fill_blast_fasta(blast_db, fasta_path, shard_params)

        with open(fasta_path) as test_file:
            expect = (
                '>seq1/1\n'
                'AAAAAAAAAA\n'
                '>seq1/2\n'
                'CCCCCCCCCC\n'
                '>seq2\n'
                'GGGGGGGGGG\n'
                '>seq3/1\n'
                'TTTTTTTTTT\n')
            assert expect == test_file.read()

    expect = [
        {'blast_db': 'test_blast_db',
         'func': 'connect',
         'module': 'lib.db'},
        {'db_conn': 'my_connection',
         'func': 'get_sequences_in_shard',
         'limit': 100,
         'module': 'lib.db',
         'offset': 200}]
    assert expect == mock.history
