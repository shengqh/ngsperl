#cd /scratch/jbrown_lab/shengq2/projects/crc-test

# java -Dconfig.file=/home/zhaos/source/perl_cqs/test/cromwell/cromwell.examples.local.conf \
#     -jar /scratch/cqs/zhaos/test/cromwell/cromwell-47.jar \
#     run /scratch/cqs_share/softwares/ngsperl/lib/WDL/CRC/crc.wdl \
#     --inputs /scratch/cqs_share/softwares/ngsperl/lib/WDL/CRC/crc.inputs.json \
#     --options /home/zhaos/source/perl_cqs/workflow/cromwell.options.json


cd /scratch/jbrown_lab/shengq2/projects/20200409_chipseq_4615_human_encode_afterTrimming/wdl_active_gene

if [[ ! -s shengqh.bioinfo.novartis.simg ]]; then
  ln -s /scratch/cqs_share/softwares/singularity/novartis.simg shengqh.bioinfo.novartis.simg 
fi

java -Dconfig.file=/home/zhaos/source/perl_cqs/test/cromwell/cromwell.examples.local.conf \
    -jar /scratch/cqs/zhaos/test/cromwell/cromwell-47.jar \
    run /scratch/cqs_share/softwares/ngsperl/lib/WDL/CRC/crc_pipeline.wdl \
    --inputs H3K27ac_25mg.json \
    --options /home/zhaos/source/perl_cqs/workflow/cromwell.options.json


#caper run /home/shengq2/program/projects/jonathan_brown/20200228_norvatis/crc.wdl -i /home/shengq2/program/projects/jonathan_brown/20200228_norvatis/crc.inputs.json --singularity /scratch/cqs_share/softwares/singularity/novartis.sif

