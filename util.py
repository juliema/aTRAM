import datetime


LOG_FILE = None
COMPLEMENT = str.maketrans('ACGTUWSMKRYBDHVNXacgtuwsmkrybdhvnx',
                           'TGCAAWSKMYRVHDBNXtgcaawskmyrvhdbnx')


def open_log_file(file_prefix):
    global LOG_FILE
    log_file = '{}log.txt'.format(file_prefix)
    LOG_FILE = open(log_file, 'w')


def log(s):
    global LOG_FILE
    LOG_FILE.write('{}: {}\n'.format(datetime.datetime.now().isoformat(' '), s))


def close_log_file():
    log('Done')
    LOG_FILE.close()


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]
