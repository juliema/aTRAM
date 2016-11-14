import datetime


LOG_FILE = None


def open_log_file(args):
    global LOG_FILE
    log_file = '{}log.txt'.format(args.out)
    LOG_FILE = open(log_file, 'w')


def log(s):
    global LOG_FILE
    LOG_FILE.write('{}: {}\n'.format(datetime.datetime.now().isoformat(' '), s))


def close_log_file():
    log('Done')
    LOG_FILE.close()
