from optparse import OptionParser
import logging

logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', datefmt='%H:%M:%S')
log = logging.getLogger("mySubmissionLogger")
log.setLevel(logging.INFO)
log.info('Loading options ...')

parser = OptionParser()
parser.add_option("--outDir", help="dir to store the output", default="/nfs/dust/atlas/user/fmeloni/MuonCollider/MuonColliderDV/Batch/scripts_reco/")
parser.add_option("--process", help="process", default="sbottom")
parser.add_option("--step", help="step", default="1")
parser.add_option("--min", help="min event", default="0")
parser.add_option("--max", help="max event", default="10")
parser.add_option("--nologs", help="logs to /dev/null", action='store_true', default=False)
(options, args) = parser.parse_args()

log.info('Preparing ...')

process_list = ["sbottom"]
if options.process not in process_list:
	log.error('Unsupported process!')
	exit()

logDir = "/nfs/dust/atlas/user/fmeloni/MuonCollider/Data/logs/"

for ievt in range(int(options.min),int(options.max),int(options.step)):

	startevent = ievt
	
	filename = "reco_"+str(options.process)+"_"+str(startevent).zfill(6)+".submit"

	fsub = open(options.outDir+filename,"w+")
	
	fsub.write("executable = /nfs/dust/atlas/user/fmeloni/MuonCollider/MuonColliderDV/Batch/RECO_Job_EVENT.sh \n")
	fsub.write("arguments = "+str(startevent)+" "+options.process+" \n")
	fsub.write("universe = vanilla"+" \n")
	if options.nologs:
		fsub.write("output = /dev/null \n")
		fsub.write("error = /dev/null \n")
		fsub.write("log = /dev/null \n")
	else:
		fsub.write("output = "+logDir+"reco_"+str(options.process)+"_"+str(startevent).zfill(6)+"_$(ClusterId).$(ProcId).out"+" \n")
		fsub.write("error = "+logDir+"reco_"+str(options.process)+"_"+str(startevent).zfill(6)+"_$(ClusterId).$(ProcId).err"+" \n")
		fsub.write("log = "+logDir+"reco_"+str(options.process)+"_"+str(startevent).zfill(6)+"_$(ClusterId).$(ProcId).log"+" \n")
	fsub.write("requirements = OpSysAndVer == \"CentOS7\" \n")
	fsub.write("RequestMemory = 2000 \n")
	fsub.write("on_exit_hold = (ExitBySignal == True) || (ExitStatus != 0) \n") 
	fsub.write("+MySingularityImage = \"/cvmfs/unpacked.cern.ch/registry.hub.docker.com/infnpd/mucoll-ilc-framework:1.6-centos8\" \n") 
	fsub.write("+MySingularityArgs = \"--no-home -B /nfs/dust/atlas/user/fmeloni/MuonCollider/MuonColliderDV:/code -B /nfs/dust/atlas/user/fmeloni/MuonCollider/Data:/data\" \n") 
	fsub.write("transfer_executable = False"+" \n")
	fsub.write("should_transfer_files = False"+" \n")
	fsub.write("periodic_release = ((JobStatus == 5) && (time() - EnteredCurrentStatus) > 240)"+" \n")
	fsub.write("queue"+" \n\n")
	fsub.close()

log.info('Scripts prepared')






