## Preprocessor arguments:

- `-h, --help`

  - List arguments

- `--version`

  - Give version number

- `-b DB, --blast-db DB, --output DB, --db DB`

  - This is the prefix of all of the blast database files so you can identify different blast database sets and so they can be stored together without resorting to subdirectories. You may include a directory as part of the prefix. The default is `atram_<today's date>`.

- `--cpus CPUS, --processes CPUS, --max-processes CPUS`

  - We default to a number that will not use up all of your cores. You may set this number to use more preocessors, but be aware that aTRAM uses a lot of temporary disk space (usually in /tmp), so you should balance the increased parallelization with the increase in temporary disk space. We should also note that you can use the `--temp-dir` option (below) to use a higher capacity disk.

- `-t DIR, --temp-dir DIR`

  - Place temporary files in this directory. All files will be deleted after aTRAM completes. The directory must exist.

- `-l LOG_FILE, --log-file LOG_FILE`

  - Log file (full path). The default is to use the DB and program name to come up with a name like so: `<DB>_atram_preprocessor.log`

- `-s SHARDS, --shards SHARDS, --number SHARDS`

  - Number of blast DB shards to create. The default is to have each shard contain roughly 250MB of sequence data. Since this parameter affects parallelization and I/O load, it may be worth experimenting with.

