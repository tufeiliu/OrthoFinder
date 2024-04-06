#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2014 David Emms
#
# This program (OrthoFinder) is distributed under the terms of the GNU General Public License v3
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#  
#  When publishing work that uses OrthoFinder please cite:
#      Emms, D.M. and Kelly, S. (2015) OrthoFinder: solving fundamental biases in whole genome comparisons dramatically 
#      improves orthogroup inference accuracy, Genome Biology 16:157
#
# For any enquiries send an email to David Emms
# david_emms@hotmail.com

from __future__ import absolute_import

# first import parallel task manager to minimise RAM overhead for small processes
import multiprocessing as mp                    # optional  (problems on OpenBSD)
import platform                                 # Y

if __name__ == "__main__":
    if platform.system() == "Darwin":
        # https://github.com/davidemms/OrthoFinder/issues/570
        # https://github.com/davidemms/OrthoFinder/issues/663
        mp.set_start_method('fork')
    # else:
    #     # Should be more RAM efficient than fork and the time penalty
    #     # should be very small as we never try to create many processes
    #     mp.set_start_method('spawn')


from ..utils import parallel_task_manager, blast_file_processor, \
    files, util, matrices, program_caller, split_ortholog_files, fasta_processor

ptm_initialised = parallel_task_manager.ParallelTaskManager_singleton()

import os                                       # Y
os.environ["OPENBLAS_NUM_THREADS"] = "1"    # fix issue with numpy/openblas. Will mean that single threaded options aren't automatically parallelised 

import sys                                      # Y
import copy                                     # Y
import subprocess                               # Y
import glob                                     # Y
import shutil                                   # Y
import time                                     # Y
import itertools                                # Y
import datetime                                 # Y
from collections import Counter                 # Y
from scipy.optimize import curve_fit            # install
import numpy as np                              # install
import csv                                      # Y
import scipy.sparse as sparse                   # install
import os.path                                  # Y
import numpy.core.numeric as numeric            # install
from collections import defaultdict             # Y
import xml.etree.ElementTree as ET              # Y
from xml.etree.ElementTree import SubElement    # Y
from xml.dom import minidom                     # Y
try: 
    import queue
except ImportError:
    import Queue as queue                       # Y
import warnings                                 # Y

from ..orthogroups import gathering, orthogroups_set
from ..orthogroups import accelerate as acc

from ..tools import astral, tree, mcl, trees_msa
from ..gen_tree_inference import trees2ologs_of, orthologues
from . import process_args, check_dependencies, run_commands, species_info
from orthofinder import orphan_genes_version, __version__


# Get directory containing script/bundle
if getattr(sys, 'frozen', False):
    __location__ = os.path.split(sys.executable)[0]
else:
    __location__ = os.path.realpath(os.path.join(os.getcwd(), 
                                    os.path.dirname(os.path.dirname(__file__))))

configfile_location = os.path.join(__location__, "run")  
max_int = sys.maxsize
ok = False
while not ok:
    try:
        csv.field_size_limit(max_int)
        ok = True
    except OverflowError:
        max_int = int(max_int/10)
sys.setrecursionlimit(10**6)

# uncomment to get round problem with python multiprocessing library that can set all cpu affinities to a single cpu
# This can cause use of only a limited number of cpus in other cases so it has been commented out
# if sys.platform.startswith("linux"):
#     with open(os.devnull, "w") as f:
#         subprocess.call("taskset -p 0xffffffffffff %d" % os.getpid(), shell=True, stdout=f) 

my_env = os.environ.copy()
# use orthofinder supplied executables by preference
my_env['PATH'] = os.path.join(__location__, 'bin:') + my_env['PATH']
# Fix LD_LIBRARY_PATH when using pyinstaller 
if getattr(sys, 'frozen', False):
    if 'LD_LIBRARY_PATH_ORIG' in my_env:
        my_env['LD_LIBRARY_PATH'] = my_env['LD_LIBRARY_PATH_ORIG']  
    else:
        my_env['LD_LIBRARY_PATH'] = ''  
    if 'DYLD_LIBRARY_PATH_ORIG' in my_env:
        my_env['DYLD_LIBRARY_PATH'] = my_env['DYLD_LIBRARY_PATH_ORIG']  
    else:
        my_env['DYLD_LIBRARY_PATH'] = ''



"""
OrthoFinder
-------------------------------------------------------------------------------
"""   
    
def GetProgramCaller():
    config_file = os.path.join(configfile_location, 'config.json') 
    pc = program_caller.ProgramCaller(config_file if os.path.exists(config_file) else None)
    config_file_user = os.path.expanduser("~/config_orthofinder_user.json")
    if os.path.exists(config_file_user):
        pc_user = program_caller.ProgramCaller(config_file_user)
        pc.Add(pc_user)
    return pc
    
"""
Main
-------------------------------------------------------------------------------
"""   

# 9
def GetOrthologues(speciesInfoObj, options, prog_caller, i_og_restart=0):
    util.PrintUnderline("Analysing Orthogroups", True)
    orthologues.OrthologuesWorkflow(speciesInfoObj.speciesToUse,
                                    speciesInfoObj.nSpAll,
                                    prog_caller,
                                    options.msa_program,
                                    options.tree_program,
                                    options.recon_method,
                                    options.nBlast,
                                    options.nProcessAlg,
                                    options.qDoubleBlast,
                                    options.qAddSpeciesToIDs,
                                    options.qTrim,
                                    options.fewer_open_files,
                                    options.speciesTreeFN,
                                    options.qStopAfterSeqs,
                                    options.qStopAfterAlignments,
                                    options.qStopAfterTrees,
                                    options.qMSATrees,
                                    options.qPhyldog,
                                    options.name,
                                    options.qSplitParaClades,
                                    options.save_space,
                                    root_from_previous = False,
                                    i_og_restart=i_og_restart)
    util.PrintTime("Done orthologues")


def BetweenCoreOrthogroupsWorkflow(continuationDir, speciesInfoObj, seqsInfo, options, prog_caller, speciesNamesDict, results_files, q_hogs):
    """
    Infer clade-specific orthogroups for the new species clades
    n_unassigned: List[int] - number of unassigned genes per species
    """
    # Get current orthogroups - original orthogroups plus genes assigned to them
    if q_hogs:
        ogs = acc.read_hogs(continuationDir, "N0")
    else:
        ogs = acc.get_original_orthogroups()
    i_og_restart = 0
    ogs_new_species, _ = acc.assign_genes(results_files)
    clustersFilename_pairs = acc.write_all_orthogroups(ogs, ogs_new_species, [])  # this updates ogs

    ogSet = orthogroups_set.OrthoGroupsSet(files.FileHandler.GetWorkingDirectory1_Read(), speciesInfoObj.speciesToUse,
                                       speciesInfoObj.nSpAll, options.qAddSpeciesToIDs, idExtractor=util.FirstWordExtractor)

    if options.qStopAfterGroups and options.speciesTreeFN is None:
        # Can't infer clade-specific groups
        print("\nSpecies tree required for clade-speicfic orthogroups - skipping")
        return clustersFilename_pairs, i_og_restart

    n_unassigned = acc.write_unassigned_fasta(ogs, None, speciesInfoObj)

    # Get/Infer species tree
    if options.speciesTreeFN is None:
        # Infer gene trees
        # We write orthogroup & stats results files in the following code, which we should avoid & only do once all OGs are done.
        gathering.post_clustering_orthogroups(clustersFilename_pairs, speciesInfoObj, seqsInfo, speciesNamesDict, options,
                                              speciesXML=None, q_incremental=True)

        orthologues.InferGeneAndSpeciesTrees(ogSet,
                                              prog_caller, options.msa_program, options.tree_program,
                                              options.nBlast, options.nProcessAlg, options.qDoubleBlast,
                                              options.qAddSpeciesToIDs,
                                              options.qTrim, userSpeciesTree=None, qStopAfterSeqs=False,
                                              qStopAfterAlign=False, qMSA=options.qMSATrees,
                                              qPhyldog=False, results_name=options.name,
                                              root_from_previous=True)

        # Infer species tree
        astral_fn = files.FileHandler.GetAstralFilename()
        astral.create_input_file(files.FileHandler.GetOGsTreeDir(), astral_fn)
        species_tree_unrooted_fn = files.FileHandler.GetSpeciesTreeUnrootedFN()
        parallel_task_manager.RunCommand(astral.get_astral_command(astral_fn, species_tree_unrooted_fn, options.nBlast))

        # Root it
        core_rooted_species_tree = tree.Tree(files.FileHandler.GetCoreSpeciesTreeIDsRootedFN(), format=1)
        species_to_speices_map = lambda x: x
        rooted_species_tree_ids, qHaveSupport = trees2ologs_of.CheckAndRootTree(species_tree_unrooted_fn, core_rooted_species_tree, species_to_speices_map)
        rooted_species_tree_fn = files.FileHandler.GetSpeciesTreeIDsRootedFN()
        rooted_species_tree_ids.write(outfile=rooted_species_tree_fn)

        spTreeUnrootedFN = files.FileHandler.GetSpeciesTreeResultsFN(None, True)
        util.RenameTreeTaxa(rooted_species_tree_ids, spTreeUnrootedFN, ogSet.SpeciesDict(), qSupport=qHaveSupport, qFixNegatives=True)

        labeled_tree_fn = files.FileHandler.GetSpeciesTreeResultsNodeLabelsFN()
        util.RenameTreeTaxa(rooted_species_tree_ids, labeled_tree_fn, ogSet.SpeciesDict(), qSupport=False, qFixNegatives=True, label='N')
        i_og_restart = len(ogs)  # Need to process the clade-specific orthogroups only
    else:
        util.PrintUnderline("Using user-supplied species tree")
        spTreeFN_ids = files.FileHandler.GetSpeciesTreeUnrootedFN()
        orthologues.ConvertUserSpeciesTree(options.speciesTreeFN, ogSet.SpeciesDict(), spTreeFN_ids)
        rooted_species_tree_fn = spTreeFN_ids

    # Identify clades for clade-specific orthogroup inference
    iSpeciesCore = set(speciesInfoObj.get_original_species())
    species_clades = acc.get_new_species_clades(rooted_species_tree_fn, iSpeciesCore)
    util.PrintUnderline("Identifying clade-specific orthogroups for the following clades:", qHeavy=True)
    species_dict = ogSet.SpeciesDict()
    for i, clade in enumerate(species_clades):
        print(str(i) + ": " + ", ".join([species_dict[str(isp)] for isp in clade]))
    print("")

    # Clade-specific orthogroup inference
    run_commands.CreateSearchDatabases(speciesInfoObj, options, prog_caller, q_unassigned_genes=True)
    # provide list of clades, only run these searches (and only if the fasta files are non-empty)
    run_commands.RunSearch(options, speciesInfoObj, seqsInfo, prog_caller, n_genes_per_species=n_unassigned, species_clades=species_clades)
    # process the results files - only if they are present and non-empty
    options.v2_scores = True
    n_clades = len(species_clades)
    clustersFilename_pairs_unassigned_all = []
    for i_clade, clade in enumerate(species_clades):
        util.PrintUnderline("OrthoFinder clutering on new species clade %d of %d" % (i_clade+1, n_clades))
        print(str(i) + ": " + ", ".join([species_dict[str(isp)] for isp in clade]))
        speciesInfo_clade = copy.deepcopy(speciesInfoObj)
        speciesInfo_clade.speciesToUse = clade
        seqsInfo_clade = util.SeqsInfoRecompute(seqsInfo, clade)
        clustersFilename_pairs_unassigned = gathering.DoOrthogroups(options, speciesInfo_clade, seqsInfo_clade,
                                                         speciesNamesDict, speciesXML=None, i_unassigned=i_clade)
        clustersFilename_pairs_unassigned_all.append(clustersFilename_pairs_unassigned)
    ogs_clade_specific_list = [mcl.GetPredictedOGs(filename) for filename in clustersFilename_pairs_unassigned_all]

    # OGs have had assigned genes added to them already
    clustersFilename_pairs = acc.write_all_orthogroups(ogs, {}, ogs_clade_specific_list)
    return clustersFilename_pairs, i_og_restart

def GetOrthologues_FromTrees(options):
    orthologues.OrthologuesFromTrees(options.recon_method, options.nBlast, options.nProcessAlg, options.speciesTreeFN,
                                     options.qAddSpeciesToIDs, options.qSplitParaClades, options.fewer_open_files)

def main(args=None): 
    try:
        if args is None:
            args = sys.argv[1:]

        if args[0] == "-v" or args[0] == "--version":
            print(f"orthofinder:v{__version__}")
            sys.exit(0)

        else: 

            print("")
            print(("OrthoFinder version %s Copyright (C) 2014 David Emms\n" % __version__))

        # Create PTM right at start
            ptm_initialised = parallel_task_manager.ParallelTaskManager_singleton()
            prog_caller = GetProgramCaller()
            
            options, fastaDir, continuationDir, resultsDir_nonDefault, pickleDir_nonDefault, user_specified_M = process_args.ProcessArgs(prog_caller, args)
            
            files.InitialiseFileHandler(options, fastaDir, continuationDir, resultsDir_nonDefault, pickleDir_nonDefault)     
            print("Results directory: %s" % files.FileHandler.GetResultsDirectory1())

            check_dependencies.CheckDependencies(options, user_specified_M, prog_caller, files.FileHandler.GetWorkingDirectory1_Read()[0])
                
            # if using previous Trees etc., check these are all present - Job for orthologues
            if options.qStartFromBlast and options.qStartFromFasta:
                # 0. Check Files
                speciesInfoObj, speciesToUse_names = species_info.ProcessPreviousFiles(files.FileHandler.GetWorkingDirectory1_Read(), options.qDoubleBlast)
                print("\nAdding new species in %s to existing analysis in %s" % (fastaDir, continuationDir))
                # 3. 
                speciesInfoObj = fasta_processor.ProcessesNewFasta(fastaDir, options.dna, speciesInfoObj, speciesToUse_names)
                files.FileHandler.LogSpecies()
                options = process_args.CheckOptions(options, speciesInfoObj.speciesToUse)
                # 4.
                seqsInfo = util.GetSeqsInfo(files.FileHandler.GetWorkingDirectory1_Read(), speciesInfoObj.speciesToUse, speciesInfoObj.nSpAll)
                # 5.
                speciesXML = species_info.GetXMLSpeciesInfo(speciesInfoObj, options) if options.speciesXMLInfoFN else None
                # 6.    
                util.PrintUnderline("Dividing up work for BLAST for parallel processing")
                run_commands.CreateSearchDatabases(speciesInfoObj, options, prog_caller)
                # 7.  
                run_commands.RunSearch(options, speciesInfoObj, seqsInfo, prog_caller)
                # 8.
                speciesNamesDict = species_info.SpeciesNameDict(files.FileHandler.GetSpeciesIDsFN())
                gathering.DoOrthogroups(options, speciesInfoObj, seqsInfo, speciesNamesDict, speciesXML)
                # 9.
                if not options.qStopAfterGroups:
                    GetOrthologues(speciesInfoObj, options, prog_caller)  

            elif options.qStartFromFasta:
                # 3. 
                speciesInfoObj = None
                speciesInfoObj = fasta_processor.ProcessesNewFasta(fastaDir, options.dna)
                files.FileHandler.LogSpecies()
                options = process_args.CheckOptions(options, speciesInfoObj.speciesToUse)
                # 4
                seqsInfo = util.GetSeqsInfo(files.FileHandler.GetWorkingDirectory1_Read(), speciesInfoObj.speciesToUse, speciesInfoObj.nSpAll)
                # 5.
                speciesXML = species_info.GetXMLSpeciesInfo(speciesInfoObj, options) if options.speciesXMLInfoFN else None
                # 6.    
                util.PrintUnderline("Dividing up work for BLAST for parallel processing")
                run_commands.CreateSearchDatabases(speciesInfoObj, options, prog_caller)
                # 7. 
                run_commands.RunSearch(options, speciesInfoObj, seqsInfo, prog_caller)
                # 8.
                speciesNamesDict = species_info.SpeciesNameDict(files.FileHandler.GetSpeciesIDsFN())
                gathering.DoOrthogroups(options, speciesInfoObj, seqsInfo, speciesNamesDict, speciesXML)
                # 9. 
                if not options.qStopAfterGroups:
                    GetOrthologues(speciesInfoObj, options, prog_caller)

            elif options.qStartFromBlast:
                # 0.
                speciesInfoObj, _ = species_info.ProcessPreviousFiles(files.FileHandler.GetWorkingDirectory1_Read(), options.qDoubleBlast)
                files.FileHandler.LogSpecies()
                print("Using previously calculated BLAST results in %s" % (files.FileHandler.GetWorkingDirectory1_Read()[0]))
                options = process_args.CheckOptions(options, speciesInfoObj.speciesToUse)
                # 4.
                seqsInfo = util.GetSeqsInfo(files.FileHandler.GetWorkingDirectory1_Read(), speciesInfoObj.speciesToUse, speciesInfoObj.nSpAll)
                # 5.
                speciesXML = species_info.GetXMLSpeciesInfo(speciesInfoObj, options) if options.speciesXMLInfoFN else None
                # 8
                speciesNamesDict = species_info.SpeciesNameDict(files.FileHandler.GetSpeciesIDsFN())
                gathering.DoOrthogroups(options, speciesInfoObj, seqsInfo, speciesNamesDict, speciesXML)
                # 9
                if not options.qStopAfterGroups:
                    GetOrthologues(speciesInfoObj, options, prog_caller)

            elif options.qStartFromGroups:
                # 0.
                check_blast = not options.qMSATrees
                speciesInfoObj, _ = species_info.ProcessPreviousFiles(continuationDir, options.qDoubleBlast, check_blast=check_blast)
                files.FileHandler.LogSpecies()
                options = process_args.CheckOptions(options, speciesInfoObj.speciesToUse)
                # 9
                GetOrthologues(speciesInfoObj, options, prog_caller)

            elif options.qStartFromTrees:
                speciesInfoObj, _ = species_info.ProcessPreviousFiles(files.FileHandler.GetWorkingDirectory1_Read(), options.qDoubleBlast, check_blast=False)
                files.FileHandler.LogSpecies()
                options = process_args.CheckOptions(options, speciesInfoObj.speciesToUse)
                GetOrthologues_FromTrees(options)

            elif options.qFastAdd:
                # Prepare previous directory as database
                speciesInfoObj, speciesToUse_names = species_info.ProcessPreviousFiles(files.FileHandler.GetWorkingDirectory1_Read(), options.qDoubleBlast, check_blast=False)
                # Check previous directory has been done with MSA trees
                if not acc.check_for_orthoxcelerate(continuationDir, speciesInfoObj):
                    util.Fail()
                util.PrintUnderline("Creating orthogroup profiles")
                wd_list = files.FileHandler.GetWorkingDirectory1_Read()
                fn_diamond_db, q_hogs = acc.prepare_accelerate_database(continuationDir, wd_list, speciesInfoObj.nSpAll)
                print("\nAdding new species in %s to existing analysis in %s" % (fastaDir, continuationDir))
                speciesInfoObj = fasta_processor.ProcessesNewFasta(fastaDir, options.dna, speciesInfoObj, speciesToUse_names)

                options = process_args.CheckOptions(options, speciesInfoObj.speciesToUse)
                seqsInfo = util.GetSeqsInfo(files.FileHandler.GetWorkingDirectory1_Read(), speciesInfoObj.speciesToUse, speciesInfoObj.nSpAll)
                # Add genes to orthogroups
                results_files = acc.RunSearch(options, speciesInfoObj, fn_diamond_db, prog_caller)

                # Clade-specific genes
                speciesNamesDict = species_info.SpeciesNameDict(files.FileHandler.GetSpeciesIDsFN())
                # if orphan_genes_version == 1:
                #     # v1 - This is unsuitable, it does an all-v-all search of all unassigned genes. Although these should have
                #     # been depleted of all genes that are not clade-specific, the resulting search still takes too long.
                #     # clade_specific_orthogroups_v1 function is now inside the src/orthofinder/legacy/utils/clade_specific_orthogroups.py
                #     clustersFilename_pairs, i_og_restart = clade_specific_orthogroups_v1(speciesInfoObj, seqsInfo, options, prog_caller, speciesNamesDict, results_files, q_hogs)
                #     raise Exception("If q_hjogs then should be reading the N0.tsv file, not the original clusters")
                #     gathering.post_clustering_orthogroups(clustersFilename_pairs, speciesInfoObj, seqsInfo, speciesNamesDict, options, speciesXML=None)
                if orphan_genes_version == 2:
                    # v2 - Infer rooted species tree from new rooted gene trees, identify new species-clades & search within these
                    clustersFilename_pairs, i_og_restart = BetweenCoreOrthogroupsWorkflow(continuationDir, speciesInfoObj, seqsInfo, options, prog_caller, speciesNamesDict, results_files, q_hogs)

                    # Infer clade-specific orthogroup gene trees
                    gathering.post_clustering_orthogroups(clustersFilename_pairs, speciesInfoObj, seqsInfo,
                                                        speciesNamesDict, options, speciesXML=None)
                    if options.speciesTreeFN is None:
                        # No user species tree, use the one we've just inferred
                        options.speciesTreeFN = files.FileHandler.GetSpeciesTreeResultsFN(None, True)
                if not options.qStopAfterGroups:
                    GetOrthologues(speciesInfoObj, options, prog_caller, i_og_restart)
            else:
                raise NotImplementedError
                ptm = parallel_task_manager.ParallelTaskManager_singleton()
                ptm.Stop()
            if not options.save_space and not options.qFastAdd:
                # split up the orthologs into one file per species-pair
                split_ortholog_files.split_ortholog_files(files.FileHandler.GetOrthologuesDirectory())
            d_results = os.path.normpath(files.FileHandler.GetResultsDirectory1()) + os.path.sep
            print("\nResults:\n    %s" % d_results)
            util.PrintCitation(d_results)
            files.FileHandler.WriteToLog("OrthoFinder run completed\n", True)

    except Exception as e:
        print(str(e))
        parallel_task_manager.print_traceback(e)
        ptm = parallel_task_manager.ParallelTaskManager_singleton()
        ptm.Stop()
        raise

    except KeyboardInterrupt:
        print("\nProgram terminated by user.")
        sys.exit(1)
    
    finally:
        ptm = parallel_task_manager.ParallelTaskManager_singleton()
        ptm.Stop()
