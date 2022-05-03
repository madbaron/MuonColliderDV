import os
import logging
from optparse import OptionParser

logging.basicConfig(
	format='%(asctime)s %(levelname)s %(message)s', datefmt='%H:%M:%S')
log = logging.getLogger("myLogger")
log.setLevel(logging.INFO)
log.info('Loading options ...')

parser = OptionParser()
parser.add_option("--type", help="type of job", default="sim")
parser.add_option("--min", help="min event", default="0")
parser.add_option("--max", help="max events", default="1000")
parser.add_option("--step", help="which run", default="1")
parser.add_option("--outFolder", help="expected output folder",
                  default="/nfs/dust/atlas/user/fmeloni/MuonCollider/Data/sbottom")
(options, args) = parser.parse_args()

log.info('Resubmitting ...')

for ievt in range(int(options.min), int(options.max), int(options.step)):

    do_resub = False

    if options.type == "sim":
        path = "sim/sbottom_sim_"
        job_path = "scripts_sim/sim"+str(ievt).zfill(6)+".submit"
        size_threshold = 1000000
    elif options.type == "overlay":
        path = "overlay/bb_BIB_"
        job_path = "scripts_overlay/overlay_"+str(ievt).zfill(6)+".submit"
        size_threshold = 100000000
    elif options.type == "digi":
        path = "digi/bb_digi_"
        job_path = "scripts_digi/digi_"+str(ievt).zfill(6)+".submit"
        size_threshold = 100000000
    elif options.type == "reco":
        path = "reco/bb_reco_"
        job_path = "scripts_reco/reco_"+str(ievt).zfill(6)+".submit"
        size_threshold = 200000
    elif options.type == "plots":
        path = "histos/bb_histos_"
        job_path = "scripts_plots/histos_"+str(ievt).zfill(6)+".submit"
        size_threshold = 8000

    expected_output = options.outFolder + "/" + path + str(ievt) + ".slcio"
    if options.type == "plots":
        expected_output = options.outFolder + "/" + path + str(ievt) + ".root"
    log.debug("Checking " + expected_output)

    if not os.path.isfile(expected_output):
        do_resub = True
        log.info("Couldn't find " + expected_output)
    else:
        if os.path.getsize(expected_output) < size_threshold:
            do_resub = True
            log.info("File is too small " + expected_output)
            os.remove(expected_output)

    #if do_resub:
    #    log.info("Re-submitting")
    #    command = "condor_submit " + job_path
    #    os.system(command)

log.info('Resubmission done')
