# file_def.py : A python tool to generate files definition in pipeline configuration file

I hate to write the files configuration hash table by hand, especially when there are a lot of samples. So I wrote this tools to help me generate the files hash table.

```
usage: file_def.py [-h] -i [INPUT] -f [FILE_PATTERN] [-n [NAME_PATTERN]] [-a] [-d USE_DIR_NAME] [-r] [-p] [-x] [--meta_file [META_FILE]]
                   [--meta_file_index META_FILE_INDEX] [--meta_name_index META_NAME_INDEX] [--meta_file_delimiter [META_FILE_DELIMITER]]

Get file definition

options:
  -h, --help            show this help message and exit
  -i [INPUT], --input [INPUT]
                        Input source folder (default: None)
  -f [FILE_PATTERN], --file_pattern [FILE_PATTERN]
                        Input file pattern (default: None)
  -n [NAME_PATTERN], --name_pattern [NAME_PATTERN]
                        Input name pattern (default: None)
  -a, --add_zero        Filling digits with zero (default: False)
  -d USE_DIR_NAME, --use_dir_name USE_DIR_NAME
                        Use X level dir name as sample name, 0 means use file name (default: 0)
  -r, --recursive_dir   Find file in recursive dir (default: False)
  -p, --add_prefix_P    Add P as prefix of sample name (default: False)
  -x, --add_prefix_X    Add X as prefix of sample name (default: False)
  --meta_file [META_FILE]
                        Input meta file (default: None)
  --meta_file_index META_FILE_INDEX
                        Input index corresponding to real file name (default: None)
  --meta_name_index META_NAME_INDEX
                        Input index corresponding to expected file name (default: None)
  --meta_file_delimiter [META_FILE_DELIMITER]
                        Input meta file delimiter (default: None)

```

## Application 1: get sample name from file name

```{shell}
file_def.py -i /data/cqs/example_data/exomeseq -f "*.fastq.gz" -n "(S.)"
```

It will generate such hash table. You can copy and paste it into configuration file.

```
  files => {
    'S1' => [ '/data/cqs/example_data/exomeseq/S1-AGTGTTGC-ATGTAACG_S329_R1_001.fastq.gz', '/data/cqs/example_data/exomeseq/S1-AGTGTTGC-ATGTAACG_S329_R2_001.fastq.gz' ],
    'S2' => [ '/data/cqs/example_data/exomeseq/S2-TTACCTGG-GATTCTGA_S330_R1_001.fastq.gz', '/data/cqs/example_data/exomeseq/S2-TTACCTGG-GATTCTGA_S330_R2_001.fastq.gz' ],
    'S3' => [ '/data/cqs/example_data/exomeseq/S3-TCTATCCT-GAGAGGTT_S331_R1_001.fastq.gz', '/data/cqs/example_data/exomeseq/S3-TCTATCCT-GAGAGGTT_S331_R2_001.fastq.gz' ],
    'S4' => [ '/data/cqs/example_data/exomeseq/S4-TTCTACAT-TTGTATCA_S332_R1_001.fastq.gz', '/data/cqs/example_data/exomeseq/S4-TTCTACAT-TTGTATCA_S332_R2_001.fastq.gz' ],
  },
```

By change the name pattern, you can get another hash table

```{shell}
file_def.py -i /data/cqs/example_data/exomeseq -f "*.fastq.gz" -n "(S...)_R"
```

The result:

```
  files => {
    'S329' => [ '/data/cqs/example_data/exomeseq/S1-AGTGTTGC-ATGTAACG_S329_R1_001.fastq.gz', '/data/cqs/example_data/exomeseq/S1-AGTGTTGC-ATGTAACG_S329_R2_001.fastq.gz' ],
    'S330' => [ '/data/cqs/example_data/exomeseq/S2-TTACCTGG-GATTCTGA_S330_R1_001.fastq.gz', '/data/cqs/example_data/exomeseq/S2-TTACCTGG-GATTCTGA_S330_R2_001.fastq.gz' ],
    'S331' => [ '/data/cqs/example_data/exomeseq/S3-TCTATCCT-GAGAGGTT_S331_R1_001.fastq.gz', '/data/cqs/example_data/exomeseq/S3-TCTATCCT-GAGAGGTT_S331_R2_001.fastq.gz' ],
    'S332' => [ '/data/cqs/example_data/exomeseq/S4-TTCTACAT-TTGTATCA_S332_R1_001.fastq.gz', '/data/cqs/example_data/exomeseq/S4-TTCTACAT-TTGTATCA_S332_R2_001.fastq.gz' ],
  },
```

If you want only the first read, change the file_pattern.

```{shell}
file_def.py -i /data/cqs/example_data/exomeseq -f "*R1*.fastq.gz" -n "(S...)_R"
```

The result:

```
  files => {
    'S329' => [ '/data/cqs/example_data/exomeseq/S1-AGTGTTGC-ATGTAACG_S329_R1_001.fastq.gz' ],
    'S330' => [ '/data/cqs/example_data/exomeseq/S2-TTACCTGG-GATTCTGA_S330_R1_001.fastq.gz' ],
    'S331' => [ '/data/cqs/example_data/exomeseq/S3-TCTATCCT-GAGAGGTT_S331_R1_001.fastq.gz' ],
    'S332' => [ '/data/cqs/example_data/exomeseq/S4-TTCTACAT-TTGTATCA_S332_R1_001.fastq.gz' ],
  },

```

Let's try another folder.

```{shell}
file_def.py -i /data/cqs/example_data/chipseq -f "*.gz" -n "(.+)_R1"
```

The result:

```
  files => {
    'batch01_S1' => [ '/data/cqs/example_data/chipseq/batch01_S1_R1_001.fastq.gz' ],
    'batch01_S10' => [ '/data/cqs/example_data/chipseq/batch01_S10_R1_001.fastq.gz' ],
    'batch01_S2' => [ '/data/cqs/example_data/chipseq/batch01_S2_R1_001.fastq.gz' ],
    'batch01_S5' => [ '/data/cqs/example_data/chipseq/batch01_S5_R1_001.fastq.gz' ],
    'batch01_S6' => [ '/data/cqs/example_data/chipseq/batch01_S6_R1_001.fastq.gz' ],
    'batch01_S7' => [ '/data/cqs/example_data/chipseq/batch01_S7_R1_001.fastq.gz' ],
    'batch02_S1' => [ '/data/cqs/example_data/chipseq/batch02_S1_R1_001.fastq.gz' ],
    'batch02_S3' => [ '/data/cqs/example_data/chipseq/batch02_S3_R1_001.fastq.gz' ],
    'batch02_S4' => [ '/data/cqs/example_data/chipseq/batch02_S4_R1_001.fastq.gz' ],
    'batch02_S5' => [ '/data/cqs/example_data/chipseq/batch02_S5_R1_001.fastq.gz' ],
    'batch02_S7' => [ '/data/cqs/example_data/chipseq/batch02_S7_R1_001.fastq.gz' ],
    'batch02_S8' => [ '/data/cqs/example_data/chipseq/batch02_S8_R1_001.fastq.gz' ],
  },
```

You may notice that th batch01_S10 was in the second line. Personally, I don't like it. It will also occur in such way in the table and figure when we perform analysis. So I will like to add additional '0' in other sample names to make them in order.

```{shell}
file_def.py -i /data/cqs/example_data/chipseq -f "*.gz" -n "(.+)_R1" -a
```

The result:

```
  files => {
    'batch01_S01' => [ '/data/cqs/example_data/chipseq/batch01_S1_R1_001.fastq.gz' ],
    'batch01_S02' => [ '/data/cqs/example_data/chipseq/batch01_S2_R1_001.fastq.gz' ],
    'batch01_S05' => [ '/data/cqs/example_data/chipseq/batch01_S5_R1_001.fastq.gz' ],
    'batch01_S06' => [ '/data/cqs/example_data/chipseq/batch01_S6_R1_001.fastq.gz' ],
    'batch01_S07' => [ '/data/cqs/example_data/chipseq/batch01_S7_R1_001.fastq.gz' ],
    'batch01_S10' => [ '/data/cqs/example_data/chipseq/batch01_S10_R1_001.fastq.gz' ],
    'batch02_S01' => [ '/data/cqs/example_data/chipseq/batch02_S1_R1_001.fastq.gz' ],
    'batch02_S03' => [ '/data/cqs/example_data/chipseq/batch02_S3_R1_001.fastq.gz' ],
    'batch02_S04' => [ '/data/cqs/example_data/chipseq/batch02_S4_R1_001.fastq.gz' ],
    'batch02_S05' => [ '/data/cqs/example_data/chipseq/batch02_S5_R1_001.fastq.gz' ],
    'batch02_S07' => [ '/data/cqs/example_data/chipseq/batch02_S7_R1_001.fastq.gz' ],
    'batch02_S08' => [ '/data/cqs/example_data/chipseq/batch02_S8_R1_001.fastq.gz' ],
  },
```

## Application 2: get sample name from folder name

For single cell data, usually the file name is fixed, such as 'filtered_feature_bc_matrix.h5'. We need to use folder name as sample name. Here '-r' means to look at sub folders and '-d 1' means using first parent folder name as sample name.

```{shell}
file_def.py -i /data/cqs/example_data/singlecellseq -f "filtered_feature_bc_matrix.h5" -r -d 1
```

The result:

```
  files => {
    '8822_PL_1' => [ '/data/cqs/example_data/singlecellseq/8822-PL-1/filtered_feature_bc_matrix.h5' ],
    '8822_PL_2' => [ '/data/cqs/example_data/singlecellseq/8822-PL-2/filtered_feature_bc_matrix.h5' ],
    '8822_PL_3' => [ '/data/cqs/example_data/singlecellseq/8822-PL-3/filtered_feature_bc_matrix.h5' ],
  },
```

## Application 3: replace sample name from meta file

Sometimes, PI will provide us the preferred sample names in the meta file. We can link that meta file with real data.

Content of sample.csv.

```text
ID,sample
8822-PL-1,patient1
8822-PL-2,patient2
8822-PL-3,patient3
```

Command:

```shell
file_def.py -i /data/cqs/example_data/singlecellseq -f "filtered_feature_bc_matrix.h5" -r -d 1 --meta_file /data/cqs/example_data/singlecellseq/sample.csv --meta_file_index 0 --meta_name_index 1 --meta_file_delimiter ','
```

The result:

```text
  files => {
    'patient1' => [ '/data/cqs/example_data/singlecellseq/8822-PL-1/filtered_feature_bc_matrix.h5' ],
    'patient2' => [ '/data/cqs/example_data/singlecellseq/8822-PL-2/filtered_feature_bc_matrix.h5' ],
    'patient3' => [ '/data/cqs/example_data/singlecellseq/8822-PL-3/filtered_feature_bc_matrix.h5' ],
  },
```
