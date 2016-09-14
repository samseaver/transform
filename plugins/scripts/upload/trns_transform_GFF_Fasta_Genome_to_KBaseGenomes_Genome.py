#!/usr/bin/env python
from __future__ import print_function

#GFF3 format
#http://www.sequenceontology.org/gff3.shtml
#http://gmod.org/wiki/GFF3

# Standard imports
import sys,os,time,datetime,re
import itertools,hashlib,logging
import gzip,shutil

# Defining pythonic stderr
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# 3rd party imports
import simplejson
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable

codon_table = CodonTable.ambiguous_generic_by_name["Standard"]

# KBase imports
import biokbase.workspace.client 
import biokbase.Transform.script_utils as script_utils
import trns_transform_FASTA_DNA_Assembly_to_KBaseGenomeAnnotations_Assembly as assembly

def remove_last(substr, str):
    length = len(substr)
    index = str.rfind(substr)
    if(index >= 0):
        str = str[0:index] + str[index+length:]
    return str

def remove_version(old_id,version):
    new_id = old_id
    new_id = remove_last(version,new_id)

    #strip period or underscore
    if(new_id.startswith(".") or new_id.startswith("_") or new_id.startswith("-")):
        new_id = new_id[1:]
    if(new_id.endswith(".") or new_id.endswith("_") or new_id.endswith("-")):
        new_id = new_id[:-1]        
    new_id = new_id.replace("..",".")
    new_id = new_id.replace("__","_")
    new_id = new_id.replace("-.",".")

    return new_id

def sanitize_identifiers(feature,text):
    for cleanup in text.split(","):
        feature["ID"]=remove_version(feature["ID"],cleanup)
        if("Parent" in feature):
            feature["Parent"]=remove_version(feature["Parent"],cleanup)
    return feature

def convert_ftr_object(old_ftr,contig):
    new_ftr = dict()
    new_ftr["id"] = old_ftr["ID"]

    #GFF use 1-based integers
    substr_start = old_ftr["start"]-1 
    substr_end = old_ftr["end"]
    if(old_ftr["strand"] == "-"):
        substr_start = old_ftr["end"]-1
        substr_end = old_ftr["start"]

    dna_sequence = Seq(contig[old_ftr["start"]-1:old_ftr["end"]], IUPAC.ambiguous_dna)

    #reverse complement
    if(old_ftr["strand"] == "-"):
        dna_sequence = dna_sequence.reverse_complement()

    new_ftr["dna_sequence"]=str(dna_sequence).upper()
    new_ftr["dna_sequence_length"]=len(dna_sequence)
    new_ftr["md5"]=hashlib.md5(str(dna_sequence)).hexdigest()
    new_ftr["location"] = [[old_ftr["contig"],old_ftr["start"],old_ftr["strand"],len(dna_sequence)]]
    return new_ftr
    
def upload_genome(input_gff_file=None, input_fasta_file=None, workspace_name=None,
                  shock_service_url=None, handle_service_url=None, workspace_service_url=None,
                  taxon_reference = None, source=None, release=None, core_genome_name=None, genome_type=None,scientific_name=None,
                  level=logging.INFO, logger=None, id_cleanup=None, annotations=None, proteins=None):

    ws_client = biokbase.workspace.client.Workspace(workspace_service_url)

    #Retrieve current object, and use references from it for uploaded objects, to avoid redundancy
    original_genome = ws_client.get_objects2( { 'ignoreErrors':1, 'objects': [ {'workspace':workspace_name,'name':core_genome_name} ] } )['data'][0]['data']
    assembly_ref = None
    gff_handle_ref = None
    if(original_genome != None):
        if('assembly_ref' in original_genome):
            assembly_ref=original_genome['assembly_ref']
        if('gff_handle_ref' in original_genome):
            gff_handle_ref=original_genome['gff_handle_ref']

    #Get GO OntologyDictionary
    go_ontology = ws_client.get_object( {'workspace':'KBaseOntology', 'id':'gene_ontology'} )['data']['term_hash']

    #No time string stored in GFF
    #Fasta file headers have time strings
    time_string = str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y_%m_%d_%H_%M_%S'))

    print(assembly_ref,gff_handle_ref)

    if(assembly_ref == None):
        logger.info("Uploading Assembly")
        assembly_name = "%s_assembly" % (core_genome_name)
        assembly_ref = "%s/%s" % (workspace_name,assembly_name)
        input_directory = "/".join(input_fasta_file.split("/")[0:-1])
    
        #Files should be gzipped
        with gzip.open(input_fasta_file, 'rb') as f_in:
            with open(input_fasta_file[0:-3], 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        try:
            assembly.upload_assembly(shock_service_url = shock_service_url,
                                     handle_service_url = handle_service_url,
                                     workspace_service_url = workspace_service_url,
                                     input_directory = input_directory,
                                     workspace_name = workspace_name,
                                     assembly_name = assembly_name,
                                     source = source,
                                     date_string = time_string,
                                     taxon_reference = taxon_reference,
                                     logger = logger)
        except Exception, e: 
            logger.exception(e) 
            sys.exit(1) 

        logger.info("Assembly Uploaded as "+assembly_ref)
        os.remove(input_fasta_file[0:-3])

    ##########################################
    #Reading in Fasta file, Code taken from https://www.biostars.org/p/710/
    ##########################################
    logger.info("Reading FASTA file.") 

    contigs_sequences = dict()
    dna_size=0
    gc_count=0

    input_file_handle = gzip.open(input_fasta_file,'rb')
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in itertools.groupby(input_file_handle, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())

        try:
            fasta_header,fasta_description = header.split(' ',1)
        except:
            fasta_header = header
            fasta_description = None

        seq=seq.upper()

        #count gc content
        gc_count+=seq.count("G")
        gc_count+=seq.count("C")

        #count dna size
        dna_size+=len(seq)

        contigs_sequences[fasta_header]= {'description':fasta_description,'sequence':seq}

    gc_percent = float(gc_count) / float(dna_size)

    logger.info("Reading GFF file.") 

    header = list()
    feature_list = dict()

    gff_file_handle = gzip.open(input_gff_file, 'rb')
    current_line = gff_file_handle.readline()
    while ( current_line != '' ):
        current_line=current_line.strip()
        
        if(current_line.startswith("##")):
            header.append(current_line)
        else:
            contig_id, source, feature_type, start, end, score, strand, phase, attributes = current_line.split('\t')
            if(contig_id not in contigs_sequences):
                logger.warn("Missing contig: "+contig_id)

            if(contig_id not in feature_list):
                feature_list[contig_id]=list()

            feature = {'type':feature_type,'start':int(start),'end':int(end),'score':score,'strand':strand,'phase':phase}
            for attribute in attributes.split(";"):
                key, value = attribute.split("=")
                feature[key]=value

            #Append contig identifier
            feature["contig"]=contig_id
            feature_list[contig_id].append(feature)

        current_line = gff_file_handle.readline()

    #New code inserted to better handle feature identifiers
    #Start by extracting and group them first
    original_features_identifiers_dict = dict()
    original_features_identifiers_list = list()
    original_features_identifiers_count = dict()
    original_CDS_count=dict()
    original_features_parents_dict = dict()
    for contig in sorted(feature_list):
        for feature in feature_list[contig]:
            #We're only considering gene, mRNA, and CDS for brevity's sake
            if(feature["type"] not in ("gene","mRNA","CDS")):
                continue

            if("ID" not in feature):
                if("Name" in feature):
                    feature["ID"]=feature["Name"]
                else:
                    feature["ID"]=feature["Parent"]+"."+feature["type"]

                    #if CDS have to increment
                    if(feature["type"] == "CDS"):
                        if(feature["ID"] not in original_CDS_count):
                            original_CDS_count[feature["ID"]]=1
                        else:
                            original_CDS_count[feature["ID"]]+=1

                        feature["ID"]+="."+str(original_CDS_count[feature["ID"]])

            #Collect
            if(feature["type"] == "gene"):
                original_features_identifiers_dict[feature["ID"]]=dict()
            if(feature["type"] == "mRNA"):
                original_features_identifiers_dict[feature["Parent"]][feature["ID"]]=dict()
                original_features_parents_dict[feature["ID"]]=feature["Parent"]
            if(feature["type"] == "CDS"):
                original_features_identifiers_dict[original_features_parents_dict[feature["Parent"]]][feature["Parent"]][feature["ID"]]=1

            original_features_identifiers_list.append(feature)
            original_features_identifiers_count[feature["ID"]]=len(original_features_identifiers_list)-1

    updated_features_identifiers_dict = dict()
    updated_features_identifiers_list = list()
    updated_features_identifiers_count = dict()
    updated_features_parents_dict = dict()
    updated_CDS_count = dict()
    for gene in sorted(original_features_identifiers_dict):
#        eprint(scientific_name,release,gene,original_features_identifiers_dict[gene])

        #retrieve original object
        gene_ftr = original_features_identifiers_list[original_features_identifiers_count[gene]]
        gene_ftr = sanitize_identifiers(gene_ftr,id_cleanup)

        #store gene
        updated_features_identifiers_dict[gene_ftr["ID"]]=dict()
        updated_features_identifiers_list.append(gene_ftr)
        updated_features_identifiers_count[gene_ftr["ID"]]=len(updated_features_identifiers_list)-1

        for mRNA in sorted(original_features_identifiers_dict[gene], key=lambda x: original_features_identifiers_count[x]):
            #retrieve feature
            mRNA_ftr = original_features_identifiers_list[original_features_identifiers_count[mRNA]]

            if("PAC" in mRNA[0:3]):
                if("Name" in mRNA_ftr):
                    mRNA_ftr["ID"]=mRNA_ftr["Name"]
            mRNA_ftr = sanitize_identifiers(mRNA_ftr,id_cleanup)

            updated_features_identifiers_dict[gene_ftr["ID"]][mRNA_ftr["ID"]]=dict()
            updated_features_parents_dict[mRNA_ftr["ID"]]=mRNA_ftr["Parent"]

            updated_features_identifiers_list.append(mRNA_ftr)
            updated_features_identifiers_count[mRNA_ftr["ID"]]=len(updated_features_identifiers_list)-1

            for CDS in sorted(original_features_identifiers_dict[gene][mRNA],key=lambda x: original_features_identifiers_count[x]):
                #retrieve feature
                CDS_ftr = original_features_identifiers_list[original_features_identifiers_count[CDS]]

                if("PAC" in CDS[0:3]):
                    CDS_ftr["ID"]=mRNA_ftr["ID"]+".CDS"

                    if(CDS_ftr["ID"] not in updated_CDS_count):
                        updated_CDS_count[CDS_ftr["ID"]]=1
                    else:
                        updated_CDS_count[CDS_ftr["ID"]]+=1

                    CDS_ftr["ID"]+="."+str(updated_CDS_count[CDS_ftr["ID"]])
                    CDS_ftr["Parent"]=mRNA_ftr["ID"]

                CDS_ftr = sanitize_identifiers(CDS_ftr,id_cleanup)

                updated_features_identifiers_dict[gene_ftr["ID"]][mRNA_ftr["ID"]][CDS_ftr["ID"]]=1
                updated_features_parents_dict[CDS_ftr["ID"]]=CDS_ftr["Parent"]

                updated_features_identifiers_list.append(CDS_ftr)
                updated_features_identifiers_count[CDS_ftr["ID"]]=len(updated_features_identifiers_list)-1

#        eprint(scientific_name,release,gene_ftr["ID"],updated_features_identifiers_dict[gene_ftr["ID"]])

    features_list = list()
    for gene in sorted(updated_features_identifiers_dict):
        #retrieve updated object
        gene_ftr = updated_features_identifiers_list[updated_features_identifiers_count[gene]]

        feature_object = convert_ftr_object(gene_ftr,contigs_sequences[gene_ftr["contig"]]["sequence"])
        feature_object["type"]="gene"

        #Add ontology
        ontology_terms = dict()
        if(feature_object["id"] not in annotations):
            if(len(annotations.keys())>0):
                eprint(feature_object["id"]+" not in annotations")
        else:
            if("GO" in annotations[feature_object["id"]]):
                for go_id in annotations[feature_object["id"]]["GO"].keys():
                    if(go_id not in go_ontology):
                        eprint(go_id,"missing")
                    else:
                        if("GO" not in ontology_terms):
                            ontology_terms["GO"]=dict()
                        if(go_id not in ontology_terms["GO"]):
                            OntologyEvidence=[{"method":"GFF_Fasta_Genome_to_KBaseGenomes_Genome","timestamp":time_string,"method_version":"1.0"},
                                              {"method":"Phytozome annotation_info.txt","timestamp":time_string,"method_version":"11"}]
                            OntologyData={"id":go_id,"ontology_ref":"KBaseOntology/gene_ontology","term_name":go_ontology[go_id]["name"],"term_lineage":[],"evidence":OntologyEvidence}
                            ontology_terms["GO"][go_id]=OntologyData

#                if("KEGG/ec" in annotations[mRNA_ftr["pacid"]]):
#                    for ec in annotations[mRNA_ftr["pacid"]]["KEGG/ec"].split(","):
#                        if(ec == ""):
#                            continue
#                        if("EC" not in ontology_terms):
#                            ontology_terms["EC"]=dict()
#                        if(ec not in ontology_terms["EC"]):
#                            OntologyEvidence=[{"method":"GFF_Fasta_Genome_to_KBaseGenomes_Genome","timestamp":time_string,"method_version":"1.0"},
#                                              {"method":"Phytozome annotation_info.txt","timestamp":time_string,"method_version":"11"}]
#                            OntologyData={"id":ec,"term_lineage":[],"evidence":OntologyEvidence}
#                            ontology_terms["EC"][ec]=OntologyData

        feature_object["ontology_terms"]=ontology_terms
        features_list.append(feature_object)

        for mRNA in sorted(updated_features_identifiers_dict[gene], key=lambda x: updated_features_identifiers_count[x]):
            #retrieve updated object
            mRNA_ftr = updated_features_identifiers_list[updated_features_identifiers_count[mRNA]]
            
            feature_object = convert_ftr_object(mRNA_ftr,contigs_sequences[mRNA_ftr["contig"]]["sequence"])
            feature_object["type"]="mRNA"
            #Skipping for now
            #features_list.append(feature_object)

            #CDS aggregation needs to be done
            #Retrieve CDS list, then append locations to first one
            CDS_list = sorted(updated_features_identifiers_dict[gene][mRNA], key=lambda x: updated_features_identifiers_count[x])
            CDS_ftr = updated_features_identifiers_list[updated_features_identifiers_count[CDS_list[0]]]
            feature_object = convert_ftr_object(CDS_ftr,contigs_sequences[CDS_ftr["contig"]]["sequence"])
            #remove digit at end of CDS identifier
            feature_object["id"] = feature_object["id"].replace(".CDS.1","") 
            feature_object["type"]="CDS"

            #collect phases, and lengths of exons
            phases = [CDS_ftr["phase"]]
            exons = [len(feature_object["dna_sequence"])]

            #Remove base(s) according to phase
            feature_object["dna_sequence"] = feature_object["dna_sequence"][int(CDS_ftr["phase"]):]

            for CDS in (CDS_list[1:]):
                #retrieve updated object
                add_ftr = updated_features_identifiers_list[updated_features_identifiers_count[CDS]]
                phases.append(add_ftr["phase"])
                add_ftr_obj = convert_ftr_object(add_ftr,contigs_sequences[add_ftr["contig"]]["sequence"])
                exons.append(len(add_ftr_obj["dna_sequence"]))

                #Remove base(s) according to phase
                add_ftr_obj["dna_sequence"] = add_ftr_obj["dna_sequence"][int(add_ftr["phase"]):]

                feature_object["dna_sequence"]+=add_ftr_obj["dna_sequence"]
                feature_object["location"].append(add_ftr_obj["location"][0])

            #translate sequence
            dna_sequence = Seq(feature_object["dna_sequence"], IUPAC.ambiguous_dna)
            rna_sequence = dna_sequence.transcribe()

            #Force start codon for now
            if str(rna_sequence.upper())[:3] not in codon_table.start_codons:
                print("Forcing AUG to "+feature_object["id"],phases,exons)
                temp_seq = 'AUG'+str(rna_sequence.upper())[3:]
                rna_sequence = Seq(temp_seq, IUPAC.ambiguous_dna)

            #Append N(N) if necessary
            codon_count = len(str(rna_sequence)) % 3
            if codon_count != 0:
                temp_seq = str(rna_sequence.upper())+"N"
                if codon_count == 1:
                    temp_seq+="N"
                new_codon_count=len(temp_seq) % 3
                print("Appending N(s) to "+feature_object["id"],phases,exons)
                rna_sequence = Seq(temp_seq, IUPAC.ambiguous_dna)

            protein_sequence = Seq("")
            try:
                protein_sequence = rna_sequence.translate() #cds=True)
            except CodonTable.TranslationError as te:
                print("TranslationError for "+feature_object["id"],phases,exons,": "+str(te))

            feature_object["protein_translation"] = str(protein_sequence).upper()

            if mRNA_ftr["ID"] in proteins:
                print("Protein found for "+mRNA_ftr["ID"])
                feature_object["protein_translation"]=proteins[mRNA_ftr["ID"]]['sequence']

            feature_object["protein_translation_length"]=len(feature_object["protein_translation"])
            feature_object["md5"] = hashlib.md5(feature_object["protein_translation"]).hexdigest()

            #remove non-cds keys
            del feature_object["dna_sequence"]
            del feature_object["dna_sequence_length"]
            
            features_list.append(feature_object)

    genome = dict()
    genome["id"]=core_genome_name
    genome["scientific_name"]=scientific_name
    genome["taxon_ref"]=taxon_reference
    genome["assembly_ref"]=assembly_ref
    genome["features"]=features_list
    genome["source"]=source
    genome["domain"]="Plant"
    genome["genetic_code"]=1
    genome["gc_content"]=float("{0:.2f}".format(gc_percent))
    genome["dna_size"]=dna_size
    genome["features"]=features_list

    if(gff_handle_ref == None):
        #Upload GFF to shock
        token = os.environ.get('KB_AUTH_TOKEN') 

        shock_info = script_utils.upload_file_to_shock(logger, shock_service_url, input_gff_file, token=token)
        shock_id = shock_info["id"]
        handles = script_utils.getHandles(logger, shock_service_url, handle_service_url, [shock_id], [gff_handle_ref], token)   
        gff_handle_ref=handles[0]
    genome['gff_handle_ref'] = gff_handle_ref

    genome_string = simplejson.dumps(genome, sort_keys=True, indent=4)
    genome_file = open("Bulk_Phytozome_Upload/"+core_genome_name+'.json', 'w+')
    genome_file.write(genome_string)
    genome_file.close()

    #Provencance has a 1 MB limit.  We may want to add more like the accessions, but to be safe for now not doing that.
    genome_provenance = [{"script": __file__, "script_ver": "0.1", "description": "Genome from GFF from %s" % (source)}]

    logger.info("Attempting Genome save for %s" % (core_genome_name))
    genome_not_saved = True
    while genome_not_saved:
        try:
            genome_info =  ws_client.save_objects({"workspace":workspace_name,
                                                   "objects":[ { "type": "KBaseGenomes.Genome",
                                                                 "data": genome,
                                                                 "name": core_genome_name,
                                                                 "provenance": genome_provenance}]}) 
            genome_not_saved = False 
            logger.info("Genome saved for %s" % (core_genome_name))
        except biokbase.workspace.client.ServerError as err: 
            raise

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(prog=__file__)

    parser.add_argument('--input_gff_file', nargs='?', help='GFF file', required=True)
    parser.add_argument('--input_fasta_file', nargs='?', help='FASTA file', required=True)
    parser.add_argument('--input_annotation_file', nargs='?', help='Annotations file', required=False)
    parser.add_argument('--input_protein_file', nargs='?', help='Proteins file', required=False)
    parser.add_argument('--workspace_name', nargs='?', help='Workspace to populate', required=True)

    parser.add_argument('--shock_service_url', type=str, nargs='?', required=False, default='https://ci.kbase.us/services/shock-api/')
    parser.add_argument('--handle_service_url', type=str, nargs='?', required=False, default='https://ci.kbase.us/services/handle_service/')
    parser.add_argument('--workspace_service_url', type=str, nargs='?', required=False, default='https://ci.kbase.us/services/ws/')

    parser.add_argument('--organism', nargs='?', help='Taxon', required=True)
    parser.add_argument('--taxon_wsname', nargs='?', help='Taxon Workspace', required=False, default='ReferenceTaxons')
    parser.add_argument('--taxon_names_file', nargs='?', help='Taxon Mappings', required=False)

    parser.add_argument('--source', help="data source : examples Refseq, Genbank, Pythozyme, Gramene, etc", nargs='?', required=False, default="Phytozome") 
    parser.add_argument('--release', help="Release or version of the data.  Example Ensembl release 30", nargs='?', required=False, default = "11") 

    parser.add_argument('--id_cleanup', help="String to remove from identifiers.", nargs='?', required=False, default = "")

    args, unknown = parser.parse_known_args()

    logger = script_utils.stderrlogger(__file__)
    logger.debug(args)

    if not os.path.isfile(args.input_gff_file):
        logger.warning("{0} is not a recognizable file".format(args.input_gff_file))

    if(args.input_gff_file[-3:len(args.input_gff_file)] != '.gz'):
        logger.warning("{0} is not a gzipped file".format(args.input_gff_file))

    if not os.path.isfile(args.input_fasta_file):
        logger.warning("{0} is not a recognizable file".format(args.input_fasta_file))

    if(args.input_fasta_file[-3:len(args.input_fasta_file)] != '.gz'):
        logger.warning("{0} is not a gzipped file".format(args.input_fasta_file))

    if args.input_annotation_file != None:
        if not os.path.isfile(args.input_annotation_file):
            logger.warning("{0} is not a recognizable file".format(args.input_annotation_file))
        if args.input_annotation_file[-3:len(args.input_annotation_file)] != '.gz':
            logger.warning("{0} is not a gzipped file".format(args.input_annotation_file))

    if args.input_protein_file != None:
        if not os.path.isfile(args.input_protein_file):
            logger.warning("{0} is not a recognizable file".format(args.input_protein_file))
        if args.input_protein_file[-3:len(args.input_protein_file)] != '.gz':
            logger.warning("{0} is not a gzipped file".format(args.input_protein_file))

    if(args.id_cleanup == ""):
        args.id_cleanup=args.release
    else:
        args.id_cleanup=args.id_cleanup+","+args.release

    #Read and compile protein sequences
    print("Compiling proteins")
    proteins=dict()
    if(args.input_protein_file is not None):
        input_file_handle = gzip.open(args.input_protein_file,'rb')
        faiter = (x[1] for x in itertools.groupby(input_file_handle, lambda line: line[0] == ">"))
        for header in faiter:
            # drop the ">"
            header = header.next()[1:].strip()
            # join all sequence lines to one.
            seq = "".join(s.strip() for s in faiter.next())

            try:
                fasta_header,fasta_description = header.split(' ',1)
            except:
                fasta_header = header
                fasta_description = None

            seq=seq.upper()
            proteins[fasta_header]= {'description':fasta_description,'sequence':seq}

    print("Fetching taxonomy")
    ws_client = biokbase.workspace.client.Workspace(args.workspace_service_url)

    #Get the taxon_lookup_object
    taxon_lookup = ws_client.get_object( {'workspace':args.taxon_wsname,
                                          'id':"taxon_lookup"})['data']['taxon_lookup']
    
    #Taxon lookup dependent on full genus
    #Can use a file that has two columns, and will translate the organism text into what's in the second column
    #In order to match what's in the taxon lookup
    #Example: Athaliana    Arabidopsis thaliana
    if(args.taxon_names_file is not None):
        if(os.path.isfile(args.taxon_names_file)):
            for line in open(args.taxon_names_file):
                line=line.strip()
                array=line.split('\t')
                if(args.organism in array[0]):
                    args.organism = array[1]
                    break
        else:
            print("Warning taxon_names_file argument (%s) doesn't exist" % args.taxon_names_file)

    tax_id=0
    taxon_object_name = "unknown_taxon"
    display_sc_name = None

    if(args.organism[0:3] in taxon_lookup and args.organism in taxon_lookup[args.organism[0:3]]):
        tax_id=taxon_lookup[args.organism[0:3]][args.organism]
        taxon_object_name = "%s_taxon" % (str(tax_id))

    taxon_info = ws_client.get_objects([{"workspace": args.taxon_wsname, 
                                         "name": taxon_object_name}])[0]

    taxon_ref = "%s/%s/%s" % (taxon_info['info'][6], taxon_info['info'][0], taxon_info['info'][4])
    display_sc_name = taxon_info['data']['scientific_name']
    core_genome_name = "%s_%s_%s" % (tax_id,args.source,args.release)

    #Read and compile ontology
    print("Compiling ontology")
    annotations=dict()
    if(args.input_annotation_file is not None):
        annotation_header=list()
        input_file_handle = gzip.open(args.input_annotation_file,'rb')
        current_line = input_file_handle.readline()
        while ( current_line != '' ):
            current_line=current_line.strip()

            if(len(annotation_header)==0):
                annotation_header=current_line.split('\t')
            else:
                annotation_items=current_line.split('\t')
                identifier=annotation_items[1]
                annotation_dict=dict()
                for i in range(0,len(annotation_items)):
                    if(annotation_header[i] not in annotation_dict):
                        annotation_dict[annotation_header[i]]=dict()
                    for entry in annotation_items[i].split(","):
                        if(entry == ''):
                            continue
                        entry=entry.replace("GO:GO:","GO:")
                        annotation_dict[annotation_header[i]][entry]=1
                annotations[identifier]=annotation_dict
        
            current_line = input_file_handle.readline()

    print("Uploading genome")
    genome_type="Reference"
    upload_genome(input_gff_file=args.input_gff_file,input_fasta_file=args.input_fasta_file,workspace_name=args.workspace_name,
                  shock_service_url=args.shock_service_url,handle_service_url=args.handle_service_url,workspace_service_url=args.workspace_service_url,
                  taxon_reference=taxon_ref,source=args.source,release=args.release,core_genome_name=core_genome_name,genome_type=genome_type,scientific_name=display_sc_name,
                  logger=logger,id_cleanup=args.id_cleanup,annotations=annotations,proteins=proteins)
    sys.exit(0)
