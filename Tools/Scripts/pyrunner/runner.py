#!/usr/bin/env python3

import sys
import os
import argparse
import string
from datetime import timedelta, datetime

outnwf = 0
static_done = False
dynpft = False
extension = ""
only_model = False
only_preproc = False
no_preproc = False
only_postproc = False
no_postproc = False
keep_files = False

parser = argparse.ArgumentParser(description = "Run a RegCM simulation",
        epilog = "Built to use SLURM.", allow_abbrev = True)
parser.add_argument("namelist", help = "Namelist file")
parser.add_argument("-c","--config", metavar = "regcm_config.yaml",
        default = "regcm_config.yaml",help = "Configuration file")
parser.add_argument("-i","--init", metavar = "YYYYMMDDHH",
        help = "Initial date in YYYYMMDDHH format")
parser.add_argument("-s","--start", metavar = "YYYYMMDDHH",
        help = "Start date in YYYYMMDDHH format")
parser.add_argument("-e","--end", metavar = "YYYYMMDDHH",
        help = "End date in YYYYMMDDHH format")
parser.add_argument("-d","--directory", metavar = "DIRECTORY",
        help = "Input/output directory", default = os.getenv("REGCM_IODIR"))
parser.add_argument("-m", "--makedir",
        help = "Create input/output directory (if not exist)",
        action = "store_true")
parser.add_argument("-opre","--only-preproc",
        help = "Run only the pre-processing (no model or postproc)",
        action = "store_true")
parser.add_argument("-npre","--no-preproc",
        help = "Don't run the pre-processing",
        action = "store_true")
parser.add_argument("-omdl","--only-model",
        help = "Run only the model stage (no preproc or postproc)",
        action = "store_true")
parser.add_argument("-opst","--only-postproc",
        help = "Run only the post-processing (no preproc or model)",
        action = "store_true")
parser.add_argument("-npst","--no-postproc",
        help = "Don't run the post-processing",
        action = "store_true")
parser.add_argument("-keep","--keep-files",
        help = "Do not remove original output file (enable for nesting into)",
        action = "store_true")
args = parser.parse_args()

def parse_step(t1,step):
    from dateutil.relativedelta import relativedelta
    global outnwf
    numbers_list = [char for char in step if char.isdigit()]
    alphabets_list = [char for char in step if char.isalpha()]
    units = ''.join(alphabets_list)
    values = ''.join(numbers_list)
    if "y" in units:
        outnwf = 0
        return t1 + relativedelta(years = int(values))
    if "m" in units:
        outnwf = 0
        return t1 + relativedelta(months = int(values))
    if "d" in units:
        outnwf = 1
        return t1 + relativedelta(days = int(values))
    return None

def runner(args):
    "Run the regcm model using a namelist from start_date to end_date"
    import f90nml
    import yaml
    global extension , keep_files
    global only_model, only_preproc, only_postproc
    global no_preproc, no_postproc

    only_preproc = args.only_preproc
    no_preproc = args.no_preproc
    no_postproc = args.no_postproc
    only_model = args.only_model
    only_postproc = args.only_postproc
    keep_files = args.keep_files

    if not args.directory:
        print("An output directory must be defined.")
        print("Either define it as argument or as REGCM_IODIR environment.")
        sys.exit(1)
    iodir = os.path.abspath(args.directory)

    namelist_file = args.namelist
    try:
        namelist = f90nml.read(namelist_file)
    except:
        print("Error. First argument not a fortran namelist.")
        sys.exit(1)

    domain = namelist["terrainparam"]["domname"]
    gcmdata = namelist["globdatparam"]["dattyp"].strip( )

    try:
        mdate0 = args.init
        try:
            date0 = datetime.strptime(mdate0,"%Y%m%d%H")
        except:
            print("Init not conforming: requested format is YYYYMMDDHH")
            sys.exit(1)
    except:
        mdate0 = str(namelist["restartparam"]["mdate0"])
        date0 = datetime.strptime(mdate0,"%Y%m%d%H")
    try:
        mdate1 = args.start
        try:
            date1 = datetime.strptime(mdate1,"%Y%m%d%H")
        except:
            print("Start not conforming: requested format is YYYYMMDDHH")
            sys.exit(1)
    except:
        mdate1 = str(namelist["restartparam"]["mdate1"])
        date1 = datetime.strptime(mdate1,"%Y%m%d%H")
    try:
        mdate2 = args.end
        try:
            date2 = datetime.strptime(mdate2,"%Y%m%d%H")
        except:
            print("End not conforming: requested format is YYYYMMDDHH")
            sys.exit(1)
    except:
        mdate2 = str(namelist["restartparam"]["mdate2"])
        date2 = datetime.strptime(mdate2,"%Y%m%d%H")

    oneday = timedelta(days = 1).total_seconds()
    period = timedelta(days = (date2-date1).days).total_seconds( )
    if period < oneday:
        print("Final date is less than one day after initial date")
        sys.exit(-1)
    date2 = date1 + timedelta(days = period//oneday)

    try:
        with open(args.config,"r") as f:
            config = yaml.safe_load(f)
    except:
        print("Cannot open configuration file " + args.config)
        sys.exit(1)

    print("Init at  : " + str(date0))
    print("Start at : " + str(date1))
    print("Stop at  : " + str(date2))

    if not os.path.isdir(iodir):
        if args.makedir:
            os.mkdir(iodir)
        else:
            print("Asked not to create the output directory.")
            print("Bailing out.")
            sys.exit(1)
    if not os.path.isdir(iodir):
        print("IODIR error : " + iodir + " is not an existing directory.")
        sys.exit(1)

    icbcpath = os.path.join(iodir,gcmdata,domain,"input")
    outpath = os.path.join(iodir,gcmdata,domain,"output")

    extension = "_".join(config["RegCM_Conf"].translate(
                str.maketrans("", "", string.whitespace)).split(","))
    dynpft = "DYNPFT" in extension

    namelist["terrainparam"]["dirter"] = icbcpath
    namelist["globdatparam"]["dirglob"] = icbcpath
    namelist["outparam"]["dirout"] = outpath
    namelist["terrainparam"]["inpter"] = config["RegCM_Data"]
    namelist["globdatparam"]["inpglob"] = config["RegCM_Data"]
    namelist["physicsparam"]["radclimpath"] = os.path.join(
           config["RegCM_Data"],"MERRA2","OPPMONTH")

    os.makedirs(icbcpath, exist_ok = True)
    os.makedirs(outpath, exist_ok = True)
    os.makedirs(os.path.join(iodir,"namelists"), exist_ok = True)
    os.makedirs(os.path.join(iodir,"logs"), exist_ok = True)
    os.makedirs(os.path.join(iodir,"jobs"), exist_ok = True)
    os.makedirs(os.path.join(iodir,gcmdata,domain,"postproc"), exist_ok = True)

    if date0 != date1:
        static_done = True

    mdate0 = datetime.strftime(date0,"%Y%m%d%H")
    temp1 = date1
    temp2 = parse_step(date1, config["SLURM"]["Step"])
    namelist["outparam"]["outnwf"] = outnwf
    icbc_jobid = 0
    model_jobid = 0
    post_jobid = 0
    with open(os.path.join(iodir,
           namelist["terrainparam"]["domname"]+"_namelist.yaml"),"w") as f:
        yaml.dump(namelist.todict( ),f)

    namelist["restartparam"]["mdate0"] = int(mdate0)

    while temp2 < date2:
        mdate1 = datetime.strftime(temp1,"%Y%m%d%H")
        mdate2 = datetime.strftime(temp2,"%Y%m%d%H")
        if mdate0 == mdate1:
            namelist["restartparam"]["ifrest"] = False
        else:
            namelist["restartparam"]["ifrest"] = True
        namelist["restartparam"]["mdate1"] = int(mdate1)
        namelist["restartparam"]["mdate2"] = int(mdate2)
        namelist["globdatparam"]["gdate1"] = int(mdate1)
        namelist["globdatparam"]["gdate2"] = int(mdate2)
        icbc_jobid = icbcrun(config,namelist,mdate1,mdate2,iodir,icbc_jobid)
        model_jobid = regcmrun(config,namelist,mdate1,mdate2,iodir,
                icbc_jobid,model_jobid)
        post_jobid = postrun(config,namelist,mdate1,mdate2,iodir,
                model_jobid,post_jobid)
        temp1 = temp2
        temp2 = parse_step(temp1, config["SLURM"]["Step"])

    mdate1 = datetime.strftime(temp1,"%Y%m%d%H")
    mdate2 = datetime.strftime(date2,"%Y%m%d%H")
    if mdate0 == mdate1:
        namelist["restartparam"]["ifrest"] = False
    else:
        namelist["restartparam"]["ifrest"] = True
    namelist["restartparam"]["mdate1"] = int(mdate1)
    namelist["restartparam"]["mdate2"] = int(mdate2)
    namelist["globdatparam"]["gdate1"] = int(mdate1)
    namelist["globdatparam"]["gdate2"] = int(mdate2)
    icbc_jobid = icbcrun(config,namelist,mdate1,mdate2,iodir,icbc_jobid)
    model_jobid = regcmrun(config,namelist,mdate1,mdate2,iodir,
            icbc_jobid,model_jobid)
    post_jobid = postrun(config,namelist,mdate1,mdate2,iodir,
            model_jobid,post_jobid)

def icbcrun(c,n,d1,d2,io,jid):
    import string
    from simple_slurm import Slurm
    global static_done, dynpft, extension, only_model , no_preproc
    if no_preproc:
        return 0
    if only_model or only_postproc:
        print("No preprocessing requested. Assuming ICBC ready.")
        return 0
    tempfile = os.path.join(io,"namelists",n["terrainparam"]["domname"]+
           "_namelist." + d1 + "." + d2 + ".in")
    n.write(tempfile,force = True)
    mail_type = "NONE"
    if c["SLURM"]["Mail"]:
        mail_type = "ALL"
    partition = c["SLURM"]["ICBCPartition"]
    job_name = n["terrainparam"]["domname"] + "_ICBC"
    account = c["SLURM"]["Account"]
    ntasks_per_node = c["SLURM"]["ICBCNtasks"]
    time = timedelta(hours = c["SLURM"]["Time"])
    mail_user = c["SLURM"]["Email"]
    nodes = c["SLURM"]["ICBCNodes"]
    output = os.path.join(io,"logs",f"{Slurm.JOB_NAME}_{Slurm.JOB_ID}.out")
    error = os.path.join(io,"logs",f"{Slurm.JOB_NAME}_{Slurm.JOB_ID}.err")
    if jid > 0:
        slurm = Slurm(
                partition = partition, job_name = job_name,
                account = account, ntasks_per_node = ntasks_per_node,
                time = time, mail_type = mail_type, mail_user = mail_user,
                output = output, error = error,
                nodes = nodes, dependency = dict(afterok = jid))
    else:
        slurm = Slurm(
                partition = partition, job_name = job_name,
                account = account, ntasks_per_node = ntasks_per_node,
                time = time, mail_type = mail_type, mail_user = mail_user,
                output = output, error = error, nodes = nodes)
    if "Extra" in c["SLURM"]:
        for item in c["SLURM"]["Extra"].split(','):
            slurm.add_cmd("#SBATCH "+item)
    slurm.add_cmd("module purge")
    slurm.add_cmd("source " + c["RegCM_Env"])
    slurm.add_cmd("export OMP_NUM_THREADS=" + str(c["SLURM"]["ICBCOMP"]))
    if not static_done:
        thebin = "terrain" + extension
        slurm.add_cmd(os.path.join(c["RegCM_Path"], thebin) + " " + tempfile)
    if "CLM45" in extension:
        if (not static_done) or dynpft:
            thebin = "mksurfdata" + extension
            slurm.add_cmd(os.path.join(c["RegCM_Path"],thebin) + " " + tempfile)
    thebin = "sst" + extension
    slurm.add_cmd(os.path.join(c["RegCM_Path"], thebin) + " " + tempfile)
    thebin = "icbc" + extension
    job_id = slurm.sbatch(os.path.join(c["RegCM_Path"], thebin) +
                " " + tempfile, verbose = True, shell = "/bin/bash")
    static_done = True
    with open(os.path.join(io,"jobs","ICBC." + str(job_id) + ".job"),"w") as f:
        f.write(str(slurm))
    return job_id

def regcmrun(c,n,d1,d2,io,jid,mjid):
    import string
    from simple_slurm import Slurm
    global extension, only_preproc, only_postproc
    if only_preproc or only_postproc:
        print("No RegCM model run, assuming model output ready.")
        return 0
    tempfile = os.path.join(io,"namelists",n["terrainparam"]["domname"]+
           "_namelist." + d1 + "." + d2 + ".in")
    n.write(tempfile,force = True)
    mail_type = "NONE"
    if c["SLURM"]["Mail"]:
        mail_type = "ALL"
    partition = c["SLURM"]["ModelPartition"]
    job_name = n["terrainparam"]["domname"] + "_MODEL"
    account = c["SLURM"]["Account"]
    ntasks_per_node = c["SLURM"]["ModelNtasks"]
    time = timedelta(hours = c["SLURM"]["Time"])
    mail_user = c["SLURM"]["Email"]
    nodes = c["SLURM"]["ModelNodes"]
    output = os.path.join(io,"logs",f"{Slurm.JOB_NAME}_{Slurm.JOB_ID}.out")
    error = os.path.join(io,"logs",f"{Slurm.JOB_NAME}_{Slurm.JOB_ID}.err")
    if mjid > 0:
        if jid > 0:
            slurm = Slurm(
                partition = partition, job_name = job_name,
                account = account, ntasks_per_node = ntasks_per_node,
                time = time, mail_type = mail_type, mail_user = mail_user,
                output = output, error = error,
                nodes = nodes, dependency = dict(afterok = mjid, after = jid))
        else:
            slurm = Slurm(
                partition = partition, job_name = job_name,
                account = account, ntasks_per_node = ntasks_per_node,
                time = time, mail_type = mail_type, mail_user = mail_user,
                output = output, error = error,
                nodes = nodes, dependency = dict(afterok = mjid))
    else:
        if jid > 0:
            slurm = Slurm(
                partition = partition, job_name = job_name,
                account = account, ntasks_per_node = ntasks_per_node,
                time = time, mail_type = mail_type, mail_user = mail_user,
                output = output, error = error,
                nodes = nodes, dependency = dict(afterok = jid))
        else:
            slurm = Slurm(
                partition = partition, job_name = job_name,
                account = account, ntasks_per_node = ntasks_per_node,
                time = time, mail_type = mail_type, mail_user = mail_user,
                output = output, error = error, nodes = nodes)
    if "Extra" in c["SLURM"]:
        for item in c["SLURM"]["Extra"].split(','):
            slurm.add_cmd("#SBATCH "+item)
    slurm.add_cmd("module purge")
    slurm.add_cmd("source " + c["RegCM_Env"])
    slurm.add_cmd("export OMP_NUM_THREADS=" + str(c["SLURM"]["ModelOMP"]))
    modelbin = ("regcm" + "MPI" + extension)
    job_id = slurm.sbatch("mpirun " + os.path.join(c["RegCM_Path"], modelbin) +
                  " " + tempfile, verbose = True, shell = "/bin/bash")
    with open(os.path.join(io,"jobs","MODEL." + str(job_id) + ".job"),"w") as f:
        f.write(str(slurm))
    return job_id

def postrun(c,n,d1,d2,io,mjid,jid):
    import string
    import glob
    from simple_slurm import Slurm
    from dateutil.rrule import rrule, MONTHLY, DAILY
    global extension, only_preproc, only_model , keep_files
    global no_postrpc
    if no_postproc:
        return 0
    if only_preproc or only_model:
        print("No postprocessing requested. Assuming products ready.")
        return 0
    mail_type = "NONE"
    if c["SLURM"]["Mail"]:
        mail_type = "ALL"
    partition = c["SLURM"]["PostPartition"]
    job_name = n["terrainparam"]["domname"] + "_POST"
    account = c["SLURM"]["Account"]
    ntasks_per_node = c["SLURM"]["PostNtasks"]
    time = timedelta(hours = c["SLURM"]["Time"])
    mail_user = c["SLURM"]["Email"]
    nodes = c["SLURM"]["PostNodes"]
    output = os.path.join(io,"logs",f"{Slurm.JOB_NAME}_{Slurm.JOB_ID}.out")
    error = os.path.join(io,"logs",f"{Slurm.JOB_NAME}_{Slurm.JOB_ID}.err")
    if mjid > 0:
        if jid > 0:
            slurm = Slurm(
                partition = partition, job_name = job_name,
                account = account, ntasks_per_node = ntasks_per_node,
                time = time, mail_type = mail_type, mail_user = mail_user,
                output = output, error = error,
                nodes = nodes, dependency = dict(afterok = mjid, after = jid))
        else:
            slurm = Slurm(
                partition = partition, job_name = job_name,
                account = account, ntasks_per_node = ntasks_per_node,
                time = time, mail_type = mail_type, mail_user = mail_user,
                output = output, error = error,
                nodes = nodes, dependency = dict(afterok = mjid))
    else:
        if jid > 0:
            slurm = Slurm(
                partition = partition, job_name = job_name,
                account = account, ntasks_per_node = ntasks_per_node,
                time = time, mail_type = mail_type, mail_user = mail_user,
                output = output, error = error,
                nodes = nodes, dependency = dict(afterok = jid))
        else:
            slurm = Slurm(
                partition = partition, job_name = job_name,
                account = account, ntasks_per_node = ntasks_per_node,
                time = time, mail_type = mail_type, mail_user = mail_user,
                output = output, error = error, nodes = nodes)
    if "Extra" in c["SLURM"]:
        for item in c["SLURM"]["Extra"].split(','):
            slurm.add_cmd("#SBATCH "+item)
    slurm.add_cmd("module purge")
    slurm.add_cmd("source " + c["Pycordex_Env"])
    slurm.add_cmd("export OMP_NUM_THREADS=" + str(c["SLURM"]["PostOMP"]))
    opath = n["outparam"]["dirout"]
    domain = n["terrainparam"]["domname"]
    gcmdata = n["globdatparam"]["dattyp"].strip( )
    postpath = os.path.join(io,gcmdata,domain,"postproc")
    command = (os.path.join(c["Pycordex_Path"],"pycordexer.py") + 
           " -m " + c["CORDEX"]["email"] + 
           " -d " + c["CORDEX"]["domain"] +
           " -g " + c["CORDEX"]["global"] + 
           " -e " + c["CORDEX"]["experiment"] +
           " -b " + c["CORDEX"]["ensemble"] + 
           " -n " + c["CORDEX"]["notes"] +
           " -o " + postpath +
           " -p " + repr(ntasks_per_node) +
           " -r " + repr(c["CORDEX"]["regcm_release"]) +
           " --regcm-version-id " + repr(c["CORDEX"]["regcm_version_id"]))
    fbase = os.path.join(opath,n["terrainparam"]["domname"])
    dd1 = datetime.strptime(d1,"%Y%m%d%H")
    dd2 = datetime.strptime(d2,"%Y%m%d%H")
    if "outnwf" in n["outparam"]:
        dayfrq = int(n["outparam"]["outnwf"])
        if dayfrq > 0:
            dates = [dt for dt in rrule(DAILY, dtstart=dd1, until=dd2)]
            if dayfrq > 1:
                dates = dates[::dayfrq]
        else:
            dates = [dt for dt in rrule(MONTHLY, dtstart=dd1, until=dd2)]
    else:
        dates = [dt for dt in rrule(MONTHLY, dtstart=dd1, until=dd2)]
    dd1 = list(x.strftime("%Y%m%d%H") for x in dates[:-1])
    atmfiles = list(fbase+"_ATM."+x+'.nc' for x in dd1)
    radfiles = list(fbase+"_RAD."+x+'.nc' for x in dd1)
    srffiles = list(fbase+"_SRF."+x+'.nc' for x in dd1)
    stsfiles = list(fbase+"_STS."+x+'.nc' for x in dd1)
    slurm.add_cmd("set -e")
    slurm.add_cmd("{")
    # Produce cordex variables
    for ff in srffiles:
        slurm.add_cmd(command + " " + ff + " " + c["CORDEX"]["core_srfvars"])
    for ff in stsfiles:
        slurm.add_cmd(command + " " + ff + " " + c["CORDEX"]["core_stsvars"])
    for ff in srffiles:
        slurm.add_cmd(command + " " + ff + " " + c["CORDEX"]["tier1_srfvars"])
    for ff in radfiles:
        slurm.add_cmd(command + " " + ff + " " + c["CORDEX"]["tier1_radvars"])
    for ff in stsfiles:
        slurm.add_cmd(command + " " + ff + " " + c["CORDEX"]["tier1_stsvars"])
    for ff in atmfiles:
        slurm.add_cmd(command + " " + ff + " " + c["CORDEX"]["tier1_atmvars"])
    for ff in srffiles:
        slurm.add_cmd(command + " " + ff + " " + c["CORDEX"]["tier2_srfvars"])
    for ff in radfiles:
        slurm.add_cmd(command + " " + ff + " " + c["CORDEX"]["tier2_radvars"])
    for ff in stsfiles:
        slurm.add_cmd(command + " " + ff + " " + c["CORDEX"]["tier2_stsvars"])
    slurm.add_cmd("}")
    if not keep_files:
        slurm.add_cmd("rm -f "+(" ".join(str(x) for x in srffiles)))
        slurm.add_cmd("rm -f "+(" ".join(str(x) for x in atmfiles)))
        slurm.add_cmd("rm -f "+(" ".join(str(x) for x in radfiles)))
        slurm.add_cmd("rm -f "+(" ".join(str(x) for x in stsfiles)))
    job_id = slurm.sbatch("echo Done.", verbose = True, shell = "/bin/bash")
    with open(os.path.join(io,"jobs","POST." + str(job_id) + ".job"),"w") as f:
        f.write(str(slurm))
    return job_id

if __name__ == "__main__":
    runner(args)
