# Troubleshooting

- [Checking the system configuration](#Checking-the-system-configuration)
- [Log file](#Log-file)
- [Saving temporary data](#Saving-temporary-data)
- [Debugging assembler issues](#Debugging-assembler-issues)

## Checking the system configuration

aTRAM requires some external programs and a few set of pip installed Python
modules. We have provided a utility to check that you have the required
programs installed.

```bash
./util_check_requirements.py
```

## Log file

aTRAM log files contain more information than the screen output. This should be
your first stop for debugging aTRAM issues. The log file concatenates so you
need to scroll to the end of it to see the latest error message.

If you don't specify a location for the log file via the
`--log-file=/path/to/log_file.log` then the default location will be in the
same directory as your database files given by the `--blast-db` option.

## Saving temporary data

It is possible to save the temporary data from an aTRAM run. To do this you
need to add two arguments to atram.py or to atram_preprocessor.py:

- `--temp-dir=/path/to/existing/directory`
This will create subdirectories underneath the given directory.

- `--keep-temp-dir`
This instructs aTRAM to not delete the data after it is done. You can then
examine the directory contents or use it for
[further analysis](#Debugging-assembler-issues).

## Debugging assembler issues

*Advanced debugging technique.*

The log file will almost always have a meaningful error message but in rare
cases it can be useful to see how the assembler itself behaves. If aTRAM
crashes while running an assembler we can debug the actual assembler error if
we kept the [temporary data](#Saving-temporary-data). What I do is take the
assembler command that is typically displayed within single quotes `'` in the
aTRAM output and run it directly from the command line.

**TODO Change to use tutorial arguments**

For instance if I ran atram with the following arguments.
```bash
./atram.py \
  --query=query/Phum.PHUM003340-PA.pep.fasta \
  --blast-db=db/ptgor \
  --output-prefix=output/ptgor \
  --assembler=velvet
```

I add the temporary directory arguments and rerun atram.py:
```bash
./atram.py \
  --query=query/Phum.PHUM003340-PA.pep.fasta \
  --blast-db=db/ptgor \
  --output-prefix=output/ptgor \
  --assembler=velvet \
  --temp-dir=temp \
  --keep-temp-dir
```

And I see the following error message.
```
2019-09-13 15:06:17 ERROR: The assembler failed with error: Command
'velveth temp/atram_dt0jmqy0/ptgor_Phum.PHUM003340-PA.pep.fasta_01_7xy5cnle 31
-fasta
-shortPaired '/home/user/work/aTRAM/temp/atram_dt0jmqy0/ptgor_Phum.PHUM003340-PA.pep.fasta_01_7xy5cnle/paired_1.fasta'
'/home/user/work/aTRAM/temp/atram_dt0jmqy0/ptgor_Phum.PHUM003340-PA.pep.fasta_01_7xy5cnle/paired_2.fasta'
-short '/home/user/work/aTRAM/temp/atram_dt0jmqy0/ptgor_Phum.PHUM003340-PA.pep.fasta_01_7xy5cnle/single_1.fasta'' returned non-zero exit status 1.
```

I pick up everything between the outermost single quotes
(The assembler command itself):
```
velveth temp/atram_dt0jmqy0/ptgor_Phum.PHUM003340-PA.pep.fasta_01_7xy5cnle 31
-fasta
-shortPaired '/home/user/work/aTRAM/temp/atram_dt0jmqy0/ptgor_Phum.PHUM003340-PA.pep.fasta_01_7xy5cnle/paired_1.fasta'
'/home/user/work/aTRAM/temp/atram_dt0jmqy0/ptgor_Phum.PHUM003340-PA.pep.fasta_01_7xy5cnle/paired_2.fasta'
-short '/home/user/work/aTRAM/temp/atram_dt0jmqy0/ptgor_Phum.PHUM003340-PA.pep.fasta_01_7xy5cnle/single_1.fasta'
```
And then paste it into the command line to see what happens to the assembler.
