from ..utils import util, program_caller, files, parallel_task_manager
from ..gen_tree_inference import orthologues


def CanRunBLAST():
    if parallel_task_manager.CanRunCommand("makeblastdb -help") and parallel_task_manager.CanRunCommand("blastp -help"):
        return True
    else:
        print("ERROR: Cannot run BLAST+")
        program_caller.ProgramCaller.PrintDependencyCheckFailure("makeblastdb -help\nblastp -help")
        print("Please check BLAST+ is installed and that the executables are in the system path\n")
        return False

def CanRunMCL():
    command = "mcl -h"
    if parallel_task_manager.CanRunCommand(command):
        return True
    else:
        print("ERROR: Cannot run MCL with the command \"%s\"" % command)
        program_caller.ProgramCaller.PrintDependencyCheckFailure(command)
        print("Please check MCL is installed and in the system path. See information above.\n")
        return False

def CanRunASTRAL():
    d_deps_test = files.FileHandler.GetDependenciesCheckDir()
    fn_test = d_deps_test + "astral.input.nwk"
    with open(fn_test, 'w') as outfile:
        for _ in range(20):
            outfile.write("(((A:1,B:1):2,(C:1,D:1):2,E:3):4);i\n")
    fn_output = d_deps_test + "astral.output.nwk"
    cmd = " ".join(["astral-pro", "-i", fn_test, "-o", fn_output, "-t", "2"])
    if parallel_task_manager.CanRunCommand(cmd, qAllowStderr=True, qRequireStdout=False, qCheckReturnCode=True):
        return True
    else:
        print("ERROR: Cannot run astral-pro")
        program_caller.ProgramCaller.PrintDependencyCheckFailure(cmd)
        print("Please check astral-pro is installed and that the executables are in the system path\n")
        return False

def CheckDependencies(options, user_specified_m, prog_caller, dirForTempFiles):
    util.PrintUnderline("Checking required programs are installed")
    if not user_specified_m:
        print("Running with the recommended MSA tree inference by default. To revert to legacy method use '-M dendroblast'.\n")
    if options.qStartFromFasta or options.qFastAdd:
        if options.search_program == "blast":
            if not CanRunBLAST(): util.Fail()
        else:
            d_deps_check = files.FileHandler.GetDependenciesCheckDir()
            success, stdout, stderr, cmd = prog_caller.TestSearchMethod(options.search_program, d_deps_check)
            if not success:
                print("\nERROR: Cannot run %s" % options.search_program)
                prog_caller.PrintDependencyCheckFailure(cmd)
                print("Please check %s is installed and that the executables are in the system path\n" % options.search_program)
                util.Fail()
    if (options.qStartFromFasta or options.qStartFromBlast) and not CanRunMCL():
        util.Fail()
    if not (options.qStopAfterPrepare or options.qStopAfterSeqs or options.qStopAfterGroups or options.qStartFromTrees):
        if not orthologues.CanRunOrthologueDependencies(dirForTempFiles, 
                                                            options.qMSATrees, 
                                                            options.qPhyldog, 
                                                            options.qStopAfterTrees, 
                                                            options.msa_program, 
                                                            options.tree_program, 
                                                            options.recon_method,
                                                            prog_caller, 
                                                            options.qStopAfterAlignments):
            print("Dependencies have been met for inference of orthogroups but not for the subsequent orthologue inference.")
            print("Either install the required dependencies or use the option '-og' to stop the analysis after the inference of orthogroups.\n")
            util.Fail()
    if options.qFastAdd:
        if not CanRunASTRAL():
            util.Fail()