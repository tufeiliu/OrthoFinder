from ..utils import util, files
import itertools

def GetOrderedSearchCommands(seqsInfo, speciesInfoObj, 
                             options, prog_caller, 
                             n_genes_per_species=None, 
                             q_new_species_unassigned_genes=False):

    """ Using the nSeq1 x nSeq2 as a rough estimate of the amount of work required for a given species-pair, returns the commands 
    ordered so that the commands predicted to take the longest come first. This allows the load to be balanced better when processing 
    the BLAST commands.
    n_genes_per_species: List[int] - number of genes per species
    """
    iSpeciesPrevious = list(range(speciesInfoObj.iFirstNewSpecies))
    iSpeciesNew = list(range(speciesInfoObj.iFirstNewSpecies, speciesInfoObj.nSpAll))
    if q_new_species_unassigned_genes:
        if n_genes_per_species is not None:
            exclude = {isp for isp, n_genes in enumerate(n_genes_per_species) if n_genes == 0}
            iSpeciesNew = list(set(iSpeciesNew).difference(exclude))
        speciesPairs = [(i, j) for i, j in itertools.product(iSpeciesNew, iSpeciesNew)]
        taskSizes = [n_genes_per_species[i]*n_genes_per_species[j] for i,j in speciesPairs]
    else:
        speciesPairs = [(i, j) for i, j in itertools.product(iSpeciesNew, iSpeciesNew) if (options.qDoubleBlast or i <=j)] + \
                       [(i, j) for i, j in itertools.product(iSpeciesNew, iSpeciesPrevious) if (options.qDoubleBlast or i <=j)] + \
                       [(i, j) for i, j in itertools.product(iSpeciesPrevious, iSpeciesNew) if (options.qDoubleBlast or i <=j)]
        taskSizes = [seqsInfo.nSeqsPerSpecies[i]*seqsInfo.nSeqsPerSpecies[j] for i,j in speciesPairs]
    taskSizes, speciesPairs = util.SortArrayPairByFirst(taskSizes, speciesPairs, True)
    if options.search_program == "blast":
        commands = [" ".join(["blastp", "-outfmt", "6", "-evalue", "0.001",
                              "-query", files.FileHandler.GetSpeciesUnassignedFastaFN(iFasta) if q_new_species_unassigned_genes else files.FileHandler.GetSpeciesFastaFN(iFasta),
                              "-db", files.FileHandler.GetSpeciesDatabaseN(iDB),
                              "-out", files.FileHandler.GetBlastResultsFN(iFasta, iDB, qForCreation=True)]) for iFasta, iDB in speciesPairs]
    else:
        commands = [prog_caller.GetSearchMethodCommand_Search(options.search_program,
                        files.FileHandler.GetSpeciesUnassignedFastaFN(iFasta) if q_new_species_unassigned_genes else files.FileHandler.GetSpeciesFastaFN(iFasta),
                        files.FileHandler.GetSpeciesDatabaseN(iDB, options.search_program),
                        files.FileHandler.GetBlastResultsFN(iFasta, iDB, qForCreation=True),
                        scorematrix=options.score_matrix,
                        gapopen=options.gapopen,
                        gapextend=options.gapextend) 
                        for iFasta, iDB in speciesPairs]
    return commands


def GetOrderedSearchCommands_clades(seqsInfo, speciesInfoObj, 
                                    options, prog_caller,
                                    n_genes_per_species, 
                                    species_clades):
    """
    Search all species
    """
    exclude = {isp for isp, n_genes in enumerate(n_genes_per_species) if n_genes == 0}
    speciesPairs = []
    for clade in species_clades:
        clade = list(set(clade).difference(exclude))
        speciesPairs.extend([(i, j) for i, j in itertools.product(clade, clade)])
    if options.search_program == "blast":
        commands = [" ".join(["blastp", "-outfmt", "6", "-evalue", "0.001",
                              "-query", files.FileHandler.GetSpeciesUnassignedFastaFN(iFasta),
                              "-db", files.FileHandler.GetSpeciesDatabaseN(iDB),
                              "-out", files.FileHandler.GetBlastResultsFN(iFasta, iDB, qForCreation=True)]) for iFasta, iDB in speciesPairs]
    else:
        commands = [prog_caller.GetSearchMethodCommand_Search(options.search_program,
                        files.FileHandler.GetSpeciesUnassignedFastaFN(iFasta),
                        files.FileHandler.GetSpeciesDatabaseN(iDB, options.search_program),
                        files.FileHandler.GetBlastResultsFN(iFasta, iDB, qForCreation=True),
                        scorematrix=options.score_matrix,
                        gapopen=options.gapopen,
                        gapextend=options.gapextend) 
                        for iFasta, iDB in speciesPairs]
    return commands