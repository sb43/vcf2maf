bsub -oo logs/pindel_proct%I.log  -q normal -J 'UNM[1-83]' -P analysis-cgp -n 1 -R'select[mem>=500] span[hosts=1] rusage[mem=500]' -M 500  'perl farm_idx_exec.pl pindel_protected.cmd  $LSB_JOBINDEX'
bsub -oo logs/caveman_proct%I.log  -q normal -J 'UNM[1-83]' -P analysis-cgp -n 1 -R'select[mem>=500] span[hosts=1] rusage[mem=500]' -M 500  'perl farm_idx_exec.pl caveman_protected.cmd  $LSB_JOBINDEX'
bsub -oo logs/pindel_som%I.log  -q normal -J 'UNM[1-83]' -P analysis-cgp -n 1 -R'select[mem>=500] span[hosts=1] rusage[mem=500]' -M 500  'perl farm_idx_exec.pl pindel_somatic.cmd  $LSB_JOBINDEX'
bsub -oo logs/caveman_som%I.log  -q normal -J 'UNM[1-83]' -P analysis-cgp -n 1 -R'select[mem>=500] span[hosts=1] rusage[mem=500]' -M 500  'perl farm_idx_exec.pl caveman_somatic.cmd  $LSB_JOBINDEX'

# testing
#bsub -oo logs/test%I.log  -q normal -J 'UNM[1-5]' -P analysis-cgp -n 1 -R'select[mem>=500] span[hosts=1] rusage[mem=500]' -M 500  'perl farm_idx_exec.pl test.cmd  $LSB_JOBINDEX'
