The pipeline used to design MIP for somatic mutation

SM_MIP.py --help
usage: SM_MIP.py [-h] --input_file INPUT_FILE --output_file OUTPUT_FILE
                 [--tm_min TM_MIN] [--tm_max TM_MAX]
                 [--min_hom_length MIN_HOM_LENGTH]
                 [--max_hom_length MAX_HOM_LENGTH]
                 [--gc_threshold_min GC_THRESHOLD_MIN]
                 [--gc_threshold_max GC_THRESHOLD_MAX]

Design MIPs for somatic mutations (SM)

optional arguments:
  -h, --help            show this help message and exit
  --input_file INPUT_FILE
                        SM input file
  --output_file OUTPUT_FILE
                        Output file of designed MIPs for SM
  --tm_min TM_MIN       MIPs with Tms lower than this are ignored (default 55)
  --tm_max TM_MAX       MIPs with Tms higher than this are ignored (default
                        62)
  --min_hom_length MIN_HOM_LENGTH
                        Min hom length (default 17)
  --max_hom_length MAX_HOM_LENGTH
                        Max hom length (default 28)
  --gc_threshold_min GC_THRESHOLD_MIN
                        Low GC content threshold (default .15)
  --gc_threshold_max GC_THRESHOLD_MAX
