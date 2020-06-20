cd /scratch/jbrown_lab/shengq2/projects/crc-test

java -Dconfig.file=/home/zhaos/source/perl_cqs/test/cromwell/cromwell.examples.local.conf \
    -jar /scratch/cqs/zhaos/test/cromwell/cromwell-47.jar \
    run /scratch/cqs_share/softwares/ngsperl/lib/WDL/CRC/crc.wdl \
    --inputs /scratch/cqs_share/softwares/ngsperl/lib/WDL/CRC/crc.inputs.json \
    --options /home/zhaos/source/perl_cqs/workflow/cromwell.options.json


#caper run /home/shengq2/program/projects/jonathan_brown/20200228_norvatis/crc.wdl -i /home/shengq2/program/projects/jonathan_brown/20200228_norvatis/crc.inputs.json --singularity /scratch/cqs_share/softwares/singularity/novartis.sif

