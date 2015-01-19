#!/usr/bin/env python

# standard library imports
import os
import sys
import logging
import re
import hashlib

# 3rd party imports
import simplejson

# KBase imports
import biokbase.Transform.script_utils as script_utils


# transformation method that can be called if this module is imported
# Note the logger has different levels it could be run.  
# See: https://docs.python.org/2/library/logging.html#logging-levels
#
# The default level is set to INFO which includes everything except DEBUG
def transform(shock_service_url=None, handle_service_url=None, 
              output_file_name=None, input_directory=None, 
              working_directory=None, shock_id=None, handle_id=None, 
              input_mapping=None, fasta_reference_only=False, 
              level=logging.INFO, logger=None):
    """
    Description:
        Converts FASTA file to KBaseGenomes.ContigSet json string.  
        Note the MD5 for the contig is generated by uppercasing the sequence.
        The ContigSet MD5 is generated by taking the MD5 of joining the sorted 
        list of individual contig's MD5s with a comma separator.

    Args:
        shock_service_url: A url for the KBase SHOCK service.
        handle_service_url: A url for the KBase Handle Service.
        output_file_name: A file name where the output JSON string should be stored.  
                          If the output file name is not specified the name will default 
                          to the name of the input file appended with '_contig_set'
        input_directory: The directory the resulting json file will be written to.
        working_directory: The directory the resulting json file will be written to.
        shock_id: Shock id for the fasta file if it already exists in shock
        handle_id: Handle id for the fasta file if it already exists as a handle
        input_mapping: JSON string mapping of input files to expected types.  
                       If you don't get this you need to scan the input 
                       directory and look for your files.
        fasta_reference_only: Creates a reference to the fasta file in Shock, but does not store the sequences in the workspace object.  Not recommended unless the fasta file is larger than 1GB. This is the default behavior for files that large.
        level: Logging level, defaults to logging.INFO.
        
    Returns:
        JSON file on disk that can be saved as a KBase workspace object.

    Authors:
        Jason Baumohl, Matt Henderson
    """

    if logger is None:
        logger = script_utils.stderrlogger(__file__)
    
    logger.info("Starting conversion of FASTA to KBaseGenomes.ContigSet")
    token = os.environ.get('KB_AUTH_TOKEN')
        
    if input_mapping is None:
        logger.info("Scanning for FASTA files.")
    
        valid_extensions = [".fa",".fasta",".fna"]
    
        files = os.listdir(working_directory)
        fasta_files = [x for x in files if os.path.splitext(x)[-1] in valid_extensions]
            
        assert len(fasta_files) != 0
    
        logger.info("Found {0}".format(str(fasta_files)))

        input_file_name = files[0]
    
        if len(fasta_files) > 1:
            logger.warning("Not sure how to handle multiple FASTA files in this context. Using {0}".format(input_file_name))
    else:
        input_file_name = os.path.join(os.path.join(input_directory, "fasta_assembly"), simplejson.loads(input_mapping)["fasta_assembly"])
        
                
    logger.info("Building Object.")
 
    if not os.path.isfile(input_file_name):
        raise Exception("The input file name {0} is not a file!".format(input_file_name))        

    if not os.path.isdir(args.working_directory):
        raise Exception("The working directory {0} is not a valid directory!".format(working_directory))        

    # default if not too large
    contig_set_has_sequences = True 
    if fasta_reference_only:
        contig_set_has_sequences = False 

    fasta_filesize = os.stat(input_file_name).st_size
    if fasta_filesize > 1000000000:
        # Fasta file too large to save sequences into the ContigSet object.
        contigset_warn = """The FASTA input file seems to be too large. A ContigSet
                            object will be created without sequences, but will
                            contain a reference to the file."""
        logger.warning(contigset_warn) 
        contig_set_has_sequences = False 

    input_file_handle = open(input_file_name, 'r')
    
    fasta_header = None
    sequence_list = []
    fasta_dict = dict()
    first_header_found = False
    contig_set_md5_list = []
    # Pattern for replacing white space
    pattern = re.compile(r'\s+')
    sequence_exists = False
    
    for current_line in input_file_handle:
        if (current_line[0] == ">"):
            # found a header line
            # Wrap up previous fasta sequence
            if (not sequence_exists) and first_header_found:
                logger.error("There is no sequence related to FASTA record : {0}".format(fasta_header))        
                raise Exception("There is no sequence related to FASTA record : {0}".format(fasta_header))
                    
            if not first_header_found:
                first_header_found = True
            else:
                # build up sequence and remove all white space
                total_sequence = ''.join(sequence_list)
                total_sequence = re.sub(pattern, '', total_sequence)
                fasta_key = fasta_header.strip()
                contig_dict = dict() 
                contig_dict["id"] = fasta_key 
                contig_dict["length"] = len(total_sequence) 
                contig_dict["name"] = fasta_key 
                contig_dict["description"] = "Note MD5 is generated from uppercasing the sequence" 
                contig_md5 = hashlib.md5(total_sequence.upper()).hexdigest() 
                contig_dict["md5"] = contig_md5 
                contig_set_md5_list.append(contig_md5)
                 
                if contig_set_has_sequences: 
                    contig_dict["sequence"]= total_sequence
                else: 
                    contig_dict["sequence"]= ""
                
                fasta_dict[fasta_key] = contig_dict
               
                # get set up for next fasta sequence
                sequence_list = []
                sequence_exists = False
            
            fasta_header = current_line.replace('>','')
        else:
            sequence_list.append(current_line)
            sequence_exists = True

    input_file_handle.close()

    # wrap up last fasta sequence
    if (not sequence_exists) and first_header_found: 
        logger.error("There is no sequence related to FASTA record : {0}".format(fasta_header))        
        raise Exception("There is no sequence related to FASTA record : {0}".format(fasta_header)) 
    else: 
        # build up sequence and remove all white space      
        total_sequence = ''.join(sequence_list)
        total_sequence = re.sub(pattern, '', total_sequence)
        fasta_key = fasta_header.strip()
        contig_dict = dict()
        contig_dict["id"] = fasta_key 
        contig_dict["length"] = len(total_sequence)
        contig_dict["name"] = fasta_key
        contig_dict["description"] = "Note MD5 is generated from uppercasing the sequence" 
        contig_md5 = hashlib.md5(total_sequence.upper()).hexdigest()
        contig_dict["md5"]= contig_md5
        contig_set_md5_list.append(contig_md5)
        
        if contig_set_has_sequences: 
            contig_dict["sequence"] = total_sequence 
        else:
            contig_dict["sequence"]= ""
         
        fasta_dict[fasta_key] = contig_dict 


    if output_file_name is None:
        # default to input file name minus file extenstion adding "_contig_set" to the end
        base = os.path.basename(input_file_name)
        output_file_name = "{0}_contig_set.json".format(os.path.splitext(base)[0])
    
    contig_set_dict = dict()
    contig_set_dict["md5"] = hashlib.md5(",".join(sorted(contig_set_md5_list))).hexdigest()
    contig_set_dict["id"] = output_file_name
    contig_set_dict["name"] = output_file_name
    contig_set_dict["source"] = "KBase"
    contig_set_dict["source_id"] = os.path.basename(input_file_name) 
    contig_set_dict["contigs"] = [fasta_dict[x] for x in sorted(fasta_dict.keys())]

    if shock_id is None:
        shock_info = script_utils.upload_file_to_shock(logger, shock_service_url, input_file_name, token=token)
        shock_id = shock_info["id"]
    
    contig_set_dict["fasta_ref"] = shock_id

    # For future development if the type is updated to the handle_reference instead of a shock_reference

    # This generates the json for the object
    objectString = simplejson.dumps(contig_set_dict, sort_keys=True, indent=4)

    logger.info("ContigSet data structure creation completed.  Writing out JSON.")

    output_file_path = os.path.join(working_directory,output_file_name) 
    with open(output_file_path, "w") as outFile:
        outFile.write(objectString)
    
    logger.info("Conversion completed.")


# called only if script is run from command line
if __name__ == "__main__":
    script_details = script_utils.parse_docs(transform.__doc__)    

    import argparse

    parser = argparse.ArgumentParser(prog=__file__, 
                                     description=script_details["Description"],
                                     epilog=script_details["Authors"])
                                     
    parser.add_argument('--shock_service_url', 
                        help=script_details["Args"]["shock_service_url"],
                        action='store', type=str, nargs='?', required=True)
    parser.add_argument('--handle_service_url', 
                        help=script_details["Args"]["handle_service_url"], 
                        action='store', type=str, nargs='?', default=None, required=False)
    parser.add_argument('--input_directory', 
                        help=script_details["Args"]["input_directory"], 
                        action='store', type=str, nargs='?', required=True)
    parser.add_argument('--working_directory', 
                        help=script_details["Args"]["working_directory"], 
                        action='store', type=str, nargs='?', required=True)
    parser.add_argument('--output_file_name', 
                        help=script_details["Args"]["output_file_name"],
                        action='store', type=str, nargs='?', default=None, required=False)
    parser.add_argument('--shock_id', 
                        help=script_details["Args"]["shock_id"],
                        action='store', type=str, nargs='?', default=None, required=False)
    parser.add_argument('--handle_id', 
                        help=script_details["Args"]["handle_id"], 
                        action='store', type=str, nargs='?', default=None, required=False)

    parser.add_argument('--input_mapping', 
                        help=script_details["Args"]["input_mapping"], 
                        action='store', type=unicode, nargs='?', default=None, required=False)

    # Example of a custom argument specific to this uploader
    parser.add_argument('--fasta_reference_only', 
                        help=script_details["Args"]["fasta_reference_only"], 
                        action='store_true', required=False)

    args, unknown = parser.parse_known_args()

    logger = script_utils.stderrlogger(__file__)
    try:
        transform(shock_service_url = args.shock_service_url, 
                  handle_service_url = args.handle_service_url, 
                  output_file_name = args.output_file_name, 
                  input_directory = args.input_directory, 
                  working_directory = args.working_directory, 
                  shock_id = args.shock_id, 
                  handle_id = args.handle_id,
                  input_mapping = args.input_mapping,
                  fasta_reference_only = args.fasta_reference_only,
                  logger = logger)
    except Exception, e:
        logger.exception(e)
        sys.exit(1)
    
    sys.exit(0)

