import os
import argparse
import sys
import subprocess
import gzip
import multiprocessing
import itertools
import re
import threading
import contextlib

class Repeat:
        """Class storing repeats including their attributes"""
        def __init__(self, repeat_properties, masked=False):
                self.unit = repeat_properties["unit"]
                self.read_name = repeat_properties["read_name"]
                self.read_length = repeat_properties["read_length"]
                self.start_in_read = repeat_properties["start_in_read"]
                self.end_in_read = repeat_properties["end_in_read"]
                self.imperfections = repeat_properties["imperfections"]
                self.length = repeat_properties["length"]
                self.normalized_repeat_length = repeat_properties["normalized_repeat_length"]
                self.alignment_score = repeat_properties["alignment_score"]
                self.perfection = repeat_properties["perfection"]
                self.mismatches = repeat_properties["mismatches"]
                self.insertions = repeat_properties["insertions"]
                self.deletions = repeat_properties["deletions"]
                self.n = repeat_properties["N"]
                self.unit_offset = repeat_properties["unit_offset"]
                self.unit_length = len(repeat_properties["unit"])
                self.copy_number = self.normalized_repeat_length/self.unit_length
                self.masked = masked
                self.groups = None #groups attribute of a repeat can e.g. look like this: (('perfection', ('[', 0.0, 100.0, ')')), ('length', ('[', 0.0, 40.0, ')')))
                self.grouping_strings = None #can e.g. look like this: ['AACCCT <perfection:[100.0,100.0] length:[40.0,80.0)>'], warning: list can have multiple elements for some repeats

        def group_repeat(self, grouping, grouping_motif_setting):
                """Adds 'groups' and 'grouping_string' to repeat instance"""

                def determine_grouping_repeat_motif(grouping_motif_setting):                        
                        if grouping_motif_setting == "detected":
                                return(self.unit)
                        else:
                                reverse_complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
                                reverse_complement = "".join(reverse_complement_dict.get(base, base) for base in reversed(self.unit))
                                if grouping_motif_setting == "rc":
                                        return(reverse_complement)
                                elif grouping_motif_setting == "combine":
                                        min_fw_circ_shift = min([self.unit[-i:] + self.unit[:-i] for i in range(len(self.unit))])
                                        min_rc_circ_shift = min([reverse_complement[-i:] + reverse_complement[:-i] for i in range(len(reverse_complement))])
                                        return(min(min_fw_circ_shift, min_rc_circ_shift))
                
                def create_grouping_string(group):
                        first_attribute = True
                        for attribute in group:
                                if first_attribute:
                                        first_attribute = False
                                        grouping_string = "".join([grouping_motif, " <", attribute[0], ":", attribute[1][0], str(attribute[1][1]), ",", str(attribute[1][2]), attribute[1][3]])
                                else:
                                        grouping_string = "".join([grouping_string, " ", attribute[0], ":", attribute[1][0], str(attribute[1][1]), ",", str(attribute[1][2]), attribute[1][3]])
                        grouping_string = "".join([grouping_string, ">"])
                        return(grouping_string)        
                self.groups = []
                self.grouping_strings = []
                grouping_motif = determine_grouping_repeat_motif(grouping_motif_setting)
                if not grouping == None:
                        for group in grouping:
                                repeat_fits_bin = True
                                for attribute in group:
                                        if attribute[1][0] == "(":
                                                if attribute[1][3] == ")" and not (attribute[1][1] < getattr(self, attribute[0]) < attribute[1][2]):
                                                        repeat_fits_bin = False
                                                elif attribute[1][3] == "]" and not (attribute[1][1] < getattr(self, attribute[0]) <= attribute[1][2]):
                                                        repeat_fits_bin = False
                                        elif attribute[1][0] == "[":
                                                if attribute[1][3] == ")" and not (attribute[1][1] <= getattr(self, attribute[0]) < attribute[1][2]):
                                                        repeat_fits_bin = False
                                                elif attribute[1][3] == "]" and not (attribute[1][1] <= getattr(self, attribute[0]) <= attribute[1][2]):
                                                        repeat_fits_bin = False
                                if repeat_fits_bin:
                                        self.groups.append(group)
                                        grouping_string = create_grouping_string(group)
                                        self.grouping_strings.append(grouping_string)
                        if self.grouping_strings == []:
                                self.masked = True
                else:
                        self.grouping_strings = [grouping_motif]
               
def remove_prefix(text, prefix):
        """Replacement for removeprefix function for Python < 3.9"""
        if text.startswith(prefix):
                return text[len(prefix):]
        return text

def remove_suffix(text, suffix):
        """Replacement for removesuffix function for Python < 3.9"""
        if text.endswith(suffix):
                return text[:0-len(suffix)]
        return text

def compute_repeat_unit_offset(repeat_unit, alignment_optimal_repeat):
        """Finds repeat unit offset in repsect to repeat unit in the provided format"""
        optimal_repeat_no_gaps = alignment_optimal_repeat.replace("-","")
        repeat_unit_offset = 0
        offset_shifted_repeat_unit = repeat_unit
        for i in range(len(repeat_unit)):
                if offset_shifted_repeat_unit == optimal_repeat_no_gaps[0:len(repeat_unit)]:
                        break
                else:
                        offset_shifted_repeat_unit = offset_shifted_repeat_unit[1:] + offset_shifted_repeat_unit[0]
                        repeat_unit_offset += 1
        return(repeat_unit_offset)

def interpret_alignment(alignment_repeat_in_read, alignment_optimal_repeat, repeat_unit_offset, repeat_unit, start_in_read, end_in_read):
        """Extracts positional information of repeat imperfections and repeat unit offset from the alignment to the ideal repeat provided by Phobos output"""
        repeat_imperfections = [] #sequence variants in observed repeats in relation to optimal repeat
        insertion_offset_counter = 0
        deletion_offset_counter = 0
        current_insertion = ""
        current_deletion = ""
        for i in range(len(alignment_optimal_repeat)):
                if alignment_optimal_repeat[i] == "-":
                        current_insertion += alignment_repeat_in_read[i]
                        insertion_offset_counter += 1
                elif alignment_repeat_in_read[i] == "-":
                        current_deletion +=  alignment_optimal_repeat[i]
                        deletion_offset_counter += 1
                else: #if match or mismatch, note: indel conditions do not need to be run again after loop is done because last nucleotide in alignment is never an indel
                        if (current_insertion != ""):
                                ins_flank_end_repeat = (i + 1) - insertion_offset_counter
                                ins_flank_start_repeat =  ins_flank_end_repeat - 1 
                                ins_flank_end_repeat_unit = (ins_flank_end_repeat + repeat_unit_offset - 1)%len(repeat_unit) + 1 #not sure about +1
                                ins_flank_start_repeat_unit = (ins_flank_start_repeat + repeat_unit_offset - 1)%len(repeat_unit) + 1 #not sure about +1
                                ins_end_in_read_repeat = i - deletion_offset_counter + (start_in_read -1)
                                ins_start_in_read_repeat = ins_end_in_read_repeat - len(current_insertion) + 1
                                repeat_imperfections.append("["  + str(ins_flank_start_repeat_unit) + "_" + str(ins_flank_end_repeat_unit) + "ins" + current_insertion + "," +
                                                            str(ins_flank_start_repeat) + "_" + str(ins_flank_end_repeat) + "ins" + current_insertion + "," +
                                                            str(ins_start_in_read_repeat) + "_" + str(ins_end_in_read_repeat) +  "del" + current_insertion + "]")
                                current_insertion = ""
                        if (current_deletion != ""):
                                del_flank_end_in_read_repeat = (i + 1) - deletion_offset_counter + (start_in_read -1)
                                del_flank_start_in_read_repeat = del_flank_end_in_read_repeat - 1
                                del_end_repeat = i - insertion_offset_counter
                                del_start_repeat = del_end_repeat - (len(current_deletion)) + 1
                                del_end_repeat_unit = (del_end_repeat + repeat_unit_offset - 1)%len(repeat_unit) + 1 
                                del_start_repeat_unit = (del_start_repeat + repeat_unit_offset - 1)%len(repeat_unit) + 1
                                
                                repeat_imperfections.append("["  + str(del_start_repeat_unit) + "_" + str(del_end_repeat_unit) + "del" + current_deletion + "," +
                                                            str(del_start_repeat) + "_" + str(del_end_repeat) + "del" + current_deletion + "," +
                                                            str(del_flank_start_in_read_repeat) + "_" + str(del_flank_end_in_read_repeat) + "ins" + current_deletion + "]")
                                current_deletion = ""
                        if alignment_repeat_in_read[i] != alignment_optimal_repeat[i]:
                                mismatch_pos_repeat_unit = (i - insertion_offset_counter + repeat_unit_offset )%len(repeat_unit) + 1
                                mismatch_pos_repeat = (i + 1) - insertion_offset_counter
                                mismatch_pos_read_repeat = i + start_in_read - deletion_offset_counter
                                repeat_imperfections.append("["  + str(mismatch_pos_repeat_unit) + alignment_optimal_repeat[i] + ">" + alignment_repeat_in_read[i] + "," +
                                                                str(mismatch_pos_repeat) + alignment_optimal_repeat[i] + ">" + alignment_repeat_in_read[i] +"," +
                                                                str(mismatch_pos_read_repeat) + alignment_repeat_in_read[i] + ">" + alignment_optimal_repeat[i] + "]")
        return(repeat_imperfections)







def create_dir_if_not_present(path):
        """Creates the directory at the given path if not already present"""
        if os.path.dirname(path) != "": 
                os.makedirs(os.path.dirname(path), exist_ok=True) #requires Python >= 3.2

def determine_output_location(input_path, outputdirectory, output_prefix):
        """Determines the output path consisting of the directory and output file name (minus file type extension)"""
        if output_prefix == "":
                input_filename = os.path.basename(input_path)
                input_filename = remove_suffix(input_filename, ".gz")
                input_prefix = os.path.splitext(input_filename)[0]
                return(os.path.join(outputdirectory, input_prefix))
        else:
                return(os.path.join(outputdirectory, output_prefix))

def read_list_file(list_file_path):
        """Reads text files that are simple lists, i.e. '\n' separated strings, returns a set"""
        list_file_elements = set()
        if list_file_path != None:
                if os.path.isfile(list_file_path):
                        with open (list_file_path, "r") as list_file_open:
                                for line in list_file_open:
                                        entry = remove_suffix(line, "\n").split()[0]
                                        if entry != "":
                                                list_file_elements.add(entry)
        return(list_file_elements)

def determine_grouping_setting(grouping_setting):
        """Reads user-specified repeat grouping settings and outputs the grouping in a list"""
        if grouping_setting == None:
                return(None)
        else:
                grouping_tasks = grouping_setting.split(" ")
                grouping_dict = {}
                for grouping_task in grouping_tasks:
                        grouping_attribute = grouping_task.split(":")[0]
                        grouping_range_string = re.findall("(\(|\[)(.*?)(\)|\])", grouping_task.split(":")[1])
                        grouping_ranges = [(x[0], float(x[1].split(",")[0]), float(x[1].split(",")[1]), x[2]) for x in grouping_range_string]
                        grouping_dict[grouping_attribute] = grouping_ranges
                bin_list = [] #list, e.g.: [[('perfection', 0, 100), ('perfection', 100, 100)], [('length', 0, 30), ('length', 30, 100)]]
                for grouping_attribute in grouping_dict.keys():
                        bin_list.append([(grouping_attribute, x) for x in grouping_dict[grouping_attribute]])
                grouping_list = list(itertools.product(*bin_list)) #create all combinations from 1 item per sub-list in bin_list
                return(grouping_list)

def filter_repeats(repeat, min_perfection, max_perfection, min_rep_region_length, max_rep_region_length, min_unit_size, max_unit_size, min_copy_number, max_copy_number, read_whitelist, read_blacklist):
        """Filters repeats according to user supplied settings"""
        if (repeat.read_name in read_blacklist) or (repeat.read_name in read_whitelist) or not (
                (min_perfection <= repeat.perfection <= max_perfection) and (min_rep_region_length <= repeat.length <= max_rep_region_length) and
                (min_unit_size <= repeat.unit_length <= max_unit_size) and (min_copy_number <= repeat.copy_number <= max_copy_number)):
                repeat.masked = True

class CustomOpen:
        """Custom file open class for FASTQ(.gz) / FASTA(.gz) input files and repeat info output files"""
        def __init__(self, input_file, mode):
                self.gzipped = self.detect_gzip(input_file)
                self.file = self.open_file(input_file, mode)
                self.file_format = self.detect_file_format(input_file)
        def __enter__(self):
                if self.file_format in ("FASTA","FASTQ"): #needs to be able to return self in case of fasta/fastq such that file format attribute can be read from outside
                        return(self)
                elif self.file_format == "repeatinfo": #needs to return file such that it can be written
                        return(self.file)
        def __exit__(self, *args):
                self.file.close()
        def __iter__(self):
                return(self.file)
        def write(self, write_string):
                self.file.write(write_string)
        def close(self):
                self.file.close()
        def detect_gzip(self, input_file):
                file_extension = os.path.splitext(input_file)[1].lower()
                if file_extension in (".gz", ".gzip"):
                        with open(input_file, 'rb') as in_file_open:
                                try: #Python >= 3.7 will throw OSError when not gzip, older Python will jump into if/else
                                        if in_file_open.read(3) == b'\x1f\x8b\x08': #check for magic number
                                                return(True)
                                        else:
                                                return(False)
                                except OSError:
                                        return(False)
                else:
                        return(False)
        def detect_file_format(self, input_file):
                if self.gzipped:
                        file_gzip_extension = os.path.splitext(input_file)[1]
                        file_format_extension = os.path.splitext(remove_suffix(input_file, file_gzip_extension))[1]
                else:
                        file_format_extension = os.path.splitext(input_file)[1]
                if file_format_extension.lower() in (".fq", ".fastq"):
                        return("FASTQ")
                elif file_format_extension.lower() in (".fa", ".fasta"):
                        return("FASTA")
                elif file_format_extension.lower() == ".txt":
                        return("repeatinfo")
        def open_file(self, input_file, mode):
                if self.gzipped and mode in ("r", "rt"):
                        return(gzip.open(input_file, "rt"))
                elif not self.gzipped and mode in ("r", "rt"):
                        return(open(input_file, "rt"))
                elif self.gzipped and mode in ("w+", "wt+"):
                        return(gzip.open(input_file, "wt+"))
                elif not self.gzipped and mode in ("w+", "wt+"):
                        return(open(input_file, "wt+"))                

def add_fasta_read_identifiers(fasta):
        """Adds a unique ID (line no.) for each read to its name, avoids problems with (hypothetical) reads with same names, ID is removed later"""
        fasta_line_list = []
        read_lengths = {}
        for line_counter, line in enumerate(fasta):
                if line.strip(): #if line not empty/tabs/spaces
                        if line.startswith(">"):
                                read_name = line.split()[0][1:] + "[ID:" + str(line_counter) + "]"
                                fasta_line_list.append(">" + read_name + "\n")
                        else:
                                fasta_line_list.append(line)
                                read_lengths[read_name] = len(remove_suffix(line, "\n"))
        return(fasta_line_list, read_lengths)

def convert_fq_line_list_to_fa(fastq_line_list):
        """Creates an internal fasta file from a list of lines in fastq format and adds a unique ID (line no.) for each read to its name, avoids problems with (hypothetical) reads with same names, ID is removed later"""
        fasta_line_list = []
        read_lengths = {}
        current_read_block_lines = 0
        for line_counter, line in enumerate(fastq_line_list):
                if line.strip(): #if line not empty/tabs/spaces
                        current_read_block_lines += 1
                        if line.startswith("@") and current_read_block_lines == 1:
                                read_name = line[1:].split()[0] + "[ID:" + str(line_counter) + "]"
                                fasta_line_list.append(">" + read_name + "\n") #adds a unique ID (line no.) for each read, avoids problems with (hypothetical) reads with same names, is removed later
                        elif current_read_block_lines == 2:
                                fasta_line_list.append(line)
                                read_lengths[read_name] = len(remove_suffix(line, "\n"))
                        elif current_read_block_lines == 4:
                                current_read_block_lines = 0
        return(fasta_line_list, read_lengths)

def filter_to_longest_repeats(current_read_repeats, multi_rep_reads_setting):
        """Filters repeats detected in a single read depending on user settings"""
        if multi_rep_reads_setting == "all":
                pass
        elif multi_rep_reads_setting == "none" and len(current_read_repeats) > 1:
                for repeat in current_read_repeats:
                        repeat.masked = True
        elif multi_rep_reads_setting == "longest":
                longest_repeat = None
                for repeat in current_read_repeats:
                        if not repeat.masked:
                                if longest_repeat == None:
                                        longest_repeat = repeat
                                else:
                                        if repeat.length > longest_repeat.length or (repeat.length == longest_repeat.length and repeat.perfection > longest_repeat.perfection):
                                                longest_repeat.masked = True
                                                longest_repeat = repeat
        elif multi_rep_reads_setting == "unique_longest":
                longest_repeats = {}
                for repeat in current_read_repeats:
                        if repeat.unit in longest_repeats:
                                if repeat.length > longest_repeats[repeat.unit].length or (repeat.length == longest_repeat.length and repeat.perfection > longest_repeat.perfection):
                                        longest_repeats[repeat.unit].masked = True
                                        longest_repeats[repeat.unit] = repeat
                
def read_phobos_out_string(arguments):
        """Reads and interprets Phobos output in string format, constantly updates watcher processes with data"""
        def assign_attribute_datatypes(attribute, value):
                if attribute == "length":
                        value = int(value)
                elif attribute == "normalized_repeat_length":
                        value = int(value)
                elif attribute == "alignment_score":
                        value = int(value)
                elif attribute == "perfection":
                        value = float(value)
                elif attribute == "insertions":
                        value = int(value)
                elif attribute == "deletions":
                        value = int(value)
                elif attribute == "mismatches":
                        value = int(value)
                elif attribute == "N":
                        value = int(value)
                elif attribute == "unit":
                        value = str(value)
                elif attribute == "unit_size":
                        value = int(value)
                elif attribute == "copy_number":
                        value = float(value)
                elif attribute == "read_name":
                        value = str(value)
                return(value)

        reads, reads_file_type, phobos_path, grouping, grouping_motif_setting, min_perfection, max_perfection, min_rep_region_length, max_rep_region_length, min_unit_size, max_unit_size, min_copy_number, max_copy_number, multi_rep_reads_setting, whitelisted_reads, blacklisted_reads, countstable_file_name, repeatinfo_file_name, repeatinfo_gz_file_name, add_phobos_arguments = arguments
        attribute_mapping = {"unit":"unit", "bp":"length", "BP":"normalized_repeat_length", "pt":"alignment_score", "%":"perfection",  "mis":"mismatches", "ins":"insertions", "del":"deletions", "N":"N"}

        if reads_file_type == "FASTA":
                fasta_with_read_identifiers, read_lengths = add_fasta_read_identifiers(reads)
                fasta = "".join(fasta_with_read_identifiers)
        else:
                fasta_line_list, read_lengths = convert_fq_line_list_to_fa(reads)
                fasta = "".join(fasta_line_list)

        phobos_arguments = [phobos_path, "/dev/stdin", "/dev/stdout", "--outputFormat 1", "--reportUnit 1", "--printRepeatSeqMode 2"]

        if add_phobos_arguments != None:
                phobos_arguments.extend(add_phobos_arguments.split(";"))

        phobos_process = subprocess.Popen(phobos_arguments, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        phobos_out, phobos_err = phobos_process.communicate(input=fasta)
        phobos_output_string = phobos_out
        
        detected_repeats = []
        repeat_counts = {}
        phobos_output_string = phobos_output_string.splitlines()
        current_read = ""
        current_read_repeats = []
        for line in phobos_output_string:
                if line.startswith(">"):
                        repeat_properties = {}
                        read_name = remove_suffix(remove_prefix(line, ">"), "\n")
                        repeat_properties["read_length"] = read_lengths[read_name]
                        repeat_properties["read_name"] = read_name[:read_name.rindex("[")] #remove (internal) index from read name
                        if current_read != read_name and current_read != "":
                                filter_to_longest_repeats(current_read_repeats, multi_rep_reads_setting)
                                for repeat in current_read_repeats:
                                        if not (repeatinfo_file_name == None and repeatinfo_gz_file_name == None):
                                                detected_repeats.append(repeat)
                                        if not (countstable_file_name == None or repeat.masked):
                                                if grouping != None and repeat.groups != None:
                                                        for group_index, group in enumerate(repeat.groups):
                                                                if repeat.grouping_strings[group_index] in repeat_counts:
                                                                        repeat_counts[repeat.grouping_strings[group_index]] += 1
                                                                else:
                                                                        repeat_counts[repeat.grouping_strings[group_index]] = 1
                                                elif grouping == None:
                                                        if repeat.unit in repeat_counts:
                                                                repeat_counts[repeat.grouping_strings[0]] += 1
                                                        else:
                                                                repeat_counts[repeat.grouping_strings[0]] = 1
                        current_read_repeats = []
                        current_read = read_name                                
                else:
                        if re.search(r"(\d+) :\s+\d", line): #if output line contains the repeat details
                                alignment_repeat_in_read = ""
                                alignment_optimal_repeat = ""
                                split_attributes = remove_suffix(re.sub(" +", " ", line), "\n").split("|")
                                for i, item in enumerate(split_attributes):
                                        if (i == len(split_attributes)-1):
                                                attribute_split = item.split(" ")
                                                repeat_properties["unit"] = item.split(" ")[2]
                                        elif (i !=0):
                                                attribute_split = item.split(" ")
                                                repeat_properties[attribute_mapping[attribute_split[2]]] = assign_attribute_datatypes(attribute_mapping[attribute_split[2]], attribute_split[1])
                                        else:
                                                attribute_split = item.split(" ")
                                                repeat_properties["start_in_read"] = int(attribute_split[2])
                                                repeat_properties["end_in_read"] = int(attribute_split[4])
                                                
                        elif re.search(r"(^[\-ACGTN]+$)", line): #if output line is part of the alignment output
                                if alignment_repeat_in_read == "": #if line is the read repeat sequence part of the alignment
                                        alignment_repeat_in_read = remove_suffix(line, "\n")
                                else: #if line is the optimal repeat sequence part of the alignment
                                        alignment_optimal_repeat = remove_suffix(line, "\n")
                                        repeatunit_offset = compute_repeat_unit_offset(repeat_properties["unit"], alignment_optimal_repeat)
                                        imperfections = interpret_alignment(alignment_repeat_in_read, alignment_optimal_repeat, repeatunit_offset, repeat_properties["unit"], repeat_properties["start_in_read"], repeat_properties["end_in_read"])
                                        repeat_properties["unit_offset"] = repeatunit_offset
                                        repeat_properties["imperfections"] = imperfections
                                        current_repeat = Repeat(repeat_properties)
                                        current_repeat.group_repeat(grouping, grouping_motif_setting)
                                        filter_repeats(current_repeat, min_perfection, max_perfection, min_rep_region_length, max_rep_region_length, min_unit_size, max_unit_size, min_copy_number, max_copy_number, whitelisted_reads, blacklisted_reads) #function directly filters repeats by setting their 'masked' attribute to True, function does not return anything
                                        current_read_repeats.append(current_repeat)

        #add last repeat
        filter_to_longest_repeats(current_read_repeats, multi_rep_reads_setting)
        for repeat in current_read_repeats:
                if not (repeatinfo_file_name == None and repeatinfo_gz_file_name == None):
                        detected_repeats.append(repeat)
                if not (countstable_file_name == None or repeat.masked):
                        if grouping!=None and repeat.groups!=None:
                                for group_index, group in enumerate(repeat.groups):
                                        if repeat.grouping_strings[group_index] in repeat_counts:
                                                repeat_counts[repeat.grouping_strings[group_index]] += 1
                                        else:
                                                repeat_counts[repeat.grouping_strings[group_index]] = 1
                        elif grouping == None:
                                if repeat.unit in repeat_counts:
                                        repeat_counts[repeat.grouping_strings[0]] += 1
                                else:
                                        repeat_counts[repeat.grouping_strings[0]] = 1
                                        
        return(repeat_counts, detected_repeats)

def process_repeats(input_path, countstable_file_name, repeatinfo_file_name, repeatinfo_gz_file_name, phobos_path, processes_number, grouping, grouping_motif_setting, min_perfection, max_perfection, min_rep_region_length, max_rep_region_length, min_unit_size, max_unit_size, min_copy_number, max_copy_number, multi_rep_reads_setting, read_whitelist_file, read_blacklist_file, read_chunk_size, add_phobos_arguments):
        """Starts and kills listener processes, iterates through reads file and sends chunks of reads to read_phobos_out_string function that processes them"""

        class ThreadsafeIterator:
                """Takes an iterator/generator and makes it thread-safe by serializing call to the `next` method of given iterator/generator"""
                def __init__(self, iterator):
                        self.iterator = iterator
                        self.lock = threading.Lock()
                def __iter__(self):
                        return self
                def __next__(self):
                        with self.lock:
                                return next(self.iterator)

        def threadsafe_generator(function):
                """Decorator function that takes a generator and makes it thread-safe via ThreadsafeIterator object"""
                def generator(*a, **kw):
                        return ThreadsafeIterator(function(*a, **kw))
                return generator
       
        @threadsafe_generator
        def get_seqfile_chunks():
                with CustomOpen(input_path, "r") as seqfile_open:
                        next_n_lines = list(itertools.islice(seqfile_open, read_chunk_size))
                        while True:
                                next_line = list(itertools.islice(seqfile_open, 1))
                                if not next_n_lines or not next_line:
                                        yield(next_n_lines, seqfile_open.file_format, phobos_path, grouping, grouping_motif_setting, min_perfection, max_perfection, min_rep_region_length, max_rep_region_length, min_unit_size, max_unit_size, min_copy_number, max_copy_number, multi_rep_reads_setting, whitelisted_reads, blacklisted_reads, countstable_file_name, repeatinfo_file_name, repeatinfo_gz_file_name, add_phobos_arguments)
                                        break
                                elif seqfile_open.file_format == "FASTQ" and next_line == ["+\n"]:
                                        next_n_lines.append(next_line[0])
                                        next_n_lines.append(list(itertools.islice(seqfile_open, 1))[0])
                                        yield(next_n_lines, "FASTQ", phobos_path, grouping, grouping_motif_setting, min_perfection, max_perfection, min_rep_region_length, max_rep_region_length, min_unit_size, max_unit_size, min_copy_number, max_copy_number, multi_rep_reads_setting, whitelisted_reads, blacklisted_reads, countstable_file_name, repeatinfo_file_name, repeatinfo_gz_file_name, add_phobos_arguments)
                                        next_n_lines = list(itertools.islice(seqfile_open, read_chunk_size))
                                elif seqfile_open.file_format == "FASTA" and next_line[0].startswith(">"):
                                        yield(next_n_lines, "FASTA", phobos_path, grouping, grouping_motif_setting, min_perfection, max_perfection, min_rep_region_length, max_rep_region_length, min_unit_size, max_unit_size, min_copy_number, max_copy_number, multi_rep_reads_setting, whitelisted_reads, blacklisted_reads, countstable_file_name, repeatinfo_file_name, repeatinfo_gz_file_name, add_phobos_arguments)
                                        next_n_lines = next_line + list(itertools.islice(seqfile_open, read_chunk_size))
                                elif next_line != []:
                                        next_n_lines.append(next_line[0])

        def write_repeatinfo(repeatinfo_open, detected_repeats):
                for repeat in detected_repeats:
                        if repeat.masked:
                                repeatinfo_open.write("#")
                        repeatinfo_open.write(repeat.unit + "\t" + str(repeat.perfection) + "\t" + str(repeat.length) + "\t" + str(repeat.normalized_repeat_length) + "\t" +
                                             str(repeat.unit_offset) + "\t" + str(repeat.start_in_read) + "\t" + str(repeat.end_in_read) +
                                            "\t" + str(round(repeat.copy_number,3)) + "\t" + str(repeat.read_length) + "\t" + str(repeat.alignment_score) + "\t" + str(repeat.mismatches) + "\t" +
                                            str(repeat.insertions) + "\t" + str(repeat.deletions) + "\t" + str(repeat.n) + "\t" + repeat.read_name + "\t" +
                                            (str(repeat.grouping_strings)[1:-1].replace("'","").replace(", ", ";") if not repeat.grouping_strings == [] else "none") + "\t" +
                                            (str(repeat.imperfections)[1:-1].replace("'","").replace(", ",";") if not repeat.imperfections == [] else "none") + "\n")

        def update_counts(counter, repeat_counts):
                for group, count in repeat_counts.items():
                        if group in counter:
                                counter[group] += count
                        else:
                                counter[group] = count

        def write_counts(repeats_by_group, output_file):
                """Writes counts table for a sample"""
                with CustomOpen(output_file, "w+") as outfile_open:
                        for group, group_count in sorted(repeats_by_group.items()):
                                outfile_open.write(group + "\t" + str(group_count) + "\n")

        whitelisted_reads = read_list_file(read_whitelist_file)
        blacklisted_reads = read_list_file(read_blacklist_file)

        with multiprocessing.Pool(processes_number) as pool:
                with contextlib.ExitStack() as exitstack:
                        if countstable_file_name != None:
                                counter = {}
                                create_dir_if_not_present(countstable_file_name)
                                countstable_open = exitstack.enter_context(CustomOpen(countstable_file_name, "w+"))
                        if repeatinfo_file_name != None:
                                create_dir_if_not_present(repeatinfo_file_name)
                                repeatinfo_open = exitstack.enter_context(CustomOpen(repeatinfo_file_name, "w+"))
                                repeatinfo_open.write("unit\tperfection\tlength\tnormalized_length\tunit_offset\tstart_in_read\tend_in_read\tcopy_number\tread_length\talignment_score\tmismatches\tinsertions\tdeletions\tNs\tread_name\tgrouping\timperfections\n")
                        if repeatinfo_gz_file_name != None:
                                create_dir_if_not_present(repeatinfo_gz_file_name)
                                repeatinfo_gz_open = exitstack.enter_context(CustomOpen(repeatinfo_gz_file_name, "w+"))
                                repeatinfo_open.write("unit\tperfection\tlength\tnormalized_length\tunit_offset\tstart_in_read\tend_in_read\tcopy_number\tread_length\talignment_score\tmismatches\tinsertions\tdeletions\tNs\tread_name\tgrouping\timperfections\n")
                        for chunk in pool.imap_unordered(read_phobos_out_string, get_seqfile_chunks()):
                                chunk_repeat_counts, detected_repeats = chunk
                                if repeatinfo_file_name != None:
                                        write_repeatinfo(repeatinfo_open, detected_repeats)
                                if repeatinfo_gz_file_name != None:
                                        write_repeatinfo(repeatinfo_gz_open, detected_repeats)
                                if countstable_file_name != None:
                                        update_counts(counter, chunk_repeat_counts)
                                
        if countstable_file_name != None:
                create_dir_if_not_present(countstable_file_name)
                write_counts(counter, countstable_file_name)               
                                
def read_countstables(count_files, sample_names):
        """Reads a batch of .countstable.txt files"""
        ordered_samples = []
        repeat_counts = []
        all_repeats = set()
        for i, count_file in enumerate(count_files):
                with open(count_file, "r") as count_file_open:
                        if not sample_names == None:
                                ordered_samples.append(sample_names[i])
                        else:
                                ordered_samples.append(count_file)
                        count_file_dict = {}
                        for line in count_file_open:
                                if not line.startswith("#"):
                                        repeat_unit = line.split("\t")[0]
                                        repeat_count = remove_suffix(line.split("\t")[1], "\n")
                                        count_file_dict[repeat_unit] = repeat_count
                                        all_repeats.add(repeat_unit)
                        repeat_counts.append(count_file_dict)
                                        
        for sample_dict in repeat_counts:
                for repeat in all_repeats:
                        if not repeat in sample_dict:
                                sample_dict[repeat] = 0 #add repeats that are not present in all samples with zero values
        return(ordered_samples, repeat_counts)

def write_count_matrix(ordered_samples, repeat_counts, output_file):
        """Writes count matrix file"""
        create_dir_if_not_present(output_file)
        with open (output_file, "w+") as outfile_open:
                outfile_open.write("Repeat_group")
                for sample in ordered_samples:
                        outfile_open.write("\t" + sample)
                outfile_open.write("\n")
                for repeat in repeat_counts[0].keys():
                        outfile_open.write(repeat)
                        for i, sample in enumerate(repeat_counts):
                                outfile_open.write("\t" + str(repeat_counts[i][repeat]))
                        outfile_open.write("\n")

def parse_parameters():
        """Returns parsed command line parameters"""
        parser = argparse.ArgumentParser(description="counTR is a program that de novo detects, filters and groups tandem repeats and generates a tandem repeat count matrix for NGS samples. ")
        subparsers = parser.add_subparsers(dest="program_mode")

        parser_processrepeats = subparsers.add_parser("processrepeats", help="runs trcount from start to end")
        parser_summarizecounts = subparsers.add_parser("summarizecounts", help="summarizes 'countrepeats' output files into a single count matrix for downstream analysis")

        for subparser in [parser_processrepeats, parser_summarizecounts]:
                if subparser == parser_processrepeats:
                        subparser.add_argument("inputpath", type=str, help="path to sequencing data file in fasta(.gz) or fastq(.gz) format")
                        subparser.add_argument("outputdirectory", type=str, help="directory where the output will be written to")
                        subparser.add_argument("phobospath", type=str, help="path to Phobos executable")
                        subparser.add_argument("--outputprefix", type=str, default="", dest="output_prefix", help="prefix of output files, prefix will be taken from input file, if empty string (default: %(default)s)")
                        subparser.add_argument("--outputtype", type=str, default="c", dest="output_type", help="output to generate, countstable.txt (c), repeatinfo.txt (i), repeatinfo.txt.gz (g), concatenate the letters for multiple outputs (default: %(default)s)")
                        subparser.add_argument("--processes", type=str, default="auto", dest="processes_number", help="number of parallel processes to be used, to automatically set to maximum number of available logical cores, use 'auto' (default: %(default)s)")
                        subparser.add_argument("--grouping", type=str, default=None, dest="grouping_setting", help="repeat grouping settings, example: 'perfection:[0,100)[100,100] length:[0,30)[30,inf]' (note the single quotation marks), if 'None', repeats will be only grouped by their motif (default: %(default)s)")
                        subparser.add_argument("--groupingmotif", type=str, default="detected", dest="grouping_motif_setting", choices=['detected', 'rc', 'combine'], help="motif to report for grouping, report the detected motif as is (detected), its reverse complement (rc), or combine forward and reverse complement (combine), all motifs are reported as their lexicographically minimal string rotation (default: %(default)s)')")
                        subparser.add_argument("--minperfection", type=float, default=0, dest="min_perfection", help="minimum perfection of a repeat to be considered (default: %(default)s)")
                        subparser.add_argument("--maxperfection", type=float, default=100, dest="max_perfection", help="maximum perfection of a repeat to be considered (default: %(default)s)")
                        subparser.add_argument("--minrepeatlength", type=float, default=0, dest="min_rep_region_length", help="minimum repeat region length for a repeat to be considered (default: %(default)s)")
                        subparser.add_argument("--maxrepeatlength", type=float, default="inf", dest="max_rep_region_length", help="maximum repeat region length for a repeat to be considered (for infinite, set value to: inf) (default: %(default)s)")
                        subparser.add_argument("--minunitsize", type=float, default=0, dest="min_unit_size", help="minimum repeat unit size for a repeat to be considered (default: %(default)s)")
                        subparser.add_argument("--maxunitsize", type=float, default="inf", dest="max_unit_size", help="maximum repeat unit size for a repeat to be considered (for infinite, set value to: inf) (default: %(default)s)")
                        subparser.add_argument("--mincopynumber", type=float, default=0, dest="min_copy_number", help="minimum number of repeat unit copies in a repeat for a repeat to be considered (default: %(default)s)")
                        subparser.add_argument("--maxcopynumber", type=float, default="inf", dest="max_copy_number", help="maximum number of repeat unit copies in a repeat for a repeat to be considered (for infinite, set value to: inf) (default: %(default)s)")
                        subparser.add_argument("--multirepreads", type=str, default="all", dest="multi_rep_reads_setting", help="which repeat to consider in case of reads with multiple repeats (after other filters have been applied), either 'all' (consider all repeats for each read), 'none' (ignore multi repeat reads), 'longest' (only consider the longest repeat) or 'unique_longest' (for each unique repeat unit, only consider the longest) (default: %(default)s)")
                        subparser.add_argument("--readwhitelist", type=str, default=None, dest="read_whitelist_file", help="path to list of readnames that will not be filtered out, the rest is filtered (default: %(default)s)")
                        subparser.add_argument("--readblacklist", type=str, default=None, dest="read_blacklist_file", help="path to list of readnames that will be filtered out, the rest is kept (default: %(default)s)")
                        subparser.add_argument("--readchunksize", type=int, default=50000, dest="read_chunk_size", help="approximate number of lines that are analyzed at once in a (parallel) process (default: %(default)s)")
                        subparser.add_argument("--addphobosarguments", type=str, default=None, dest="add_phobos_arguments", help="add arguments to the default Phobos call (which is run with: --outputFormat 1 --reportUnit 1 --printRepeatSeqMode 2), example: '--indelScore -4;--mismatchScore -5' (note the single quotation marks). Warning: This command can lead to unexpected behavior and crashes, if used incorrectly (default: %(default)s)")
                elif subparser == parser_summarizecounts:
                        subparser.add_argument("outputfile", type=str, help="path to output count matrix.")
                        subparser.add_argument("inputpaths", type=str, nargs="+", help="countstable.txt files to be summarized into a count matrix")
                        subparser.add_argument("--samplenames", type=str, default=None, nargs="+", dest="sample_names", help="list of sample names to be used in the resulting header in the same order as input files. If not set, input file names will be used (default: %(default)s)")

        if len(sys.argv) < 2:
                parser.print_help()
                sys.exit(0)

        args = parser.parse_args()
        return(args)

def main():
        args = parse_parameters()
        if args.program_mode == "processrepeats":
                input_path = args.inputpath
                output_directory = args.outputdirectory
                phobos_path = args.phobospath
                output_prefix = args.output_prefix
                output_type = args.output_type
                processes_number = args.processes_number
                grouping_setting = None if args.grouping_setting in (None, "None", "none") else args.grouping_setting
                grouping_motif_setting = args.grouping_motif_setting
                min_perfection = args.min_perfection
                max_perfection = args.max_perfection
                min_rep_region_length = args.min_rep_region_length
                max_rep_region_length = args.max_rep_region_length
                min_unit_size = args.min_unit_size
                max_unit_size = args.max_unit_size
                min_copy_number = args.min_copy_number
                max_copy_number = args.max_copy_number
                multi_rep_reads_setting = args.multi_rep_reads_setting
                read_whitelist_file = args.read_whitelist_file
                read_blacklist_file = args.read_blacklist_file
                read_chunk_size = args.read_chunk_size
                add_phobos_arguments = args.add_phobos_arguments

                output_prefix = determine_output_location(input_path, output_directory, output_prefix)
                output_flags = list(output_type)

                repeatinfo_file_name = None
                repeatinfo_gz_file_name = None
                countstable_file_name = None
                if "i" in output_flags:
                        repeatinfo_file_name = "{}{}".format(output_prefix, ".repeatinfo.txt")
                if "g" in output_flags:
                        repeatinfo_gz_file_name = "{}{}".format(output_prefix, ".repeatinfo.txt.gz")
                if "c" in output_flags:
                        countstable_file_name = "{}{}".format(output_prefix, ".countstable.txt")

                if processes_number == "auto":
                        processes_number = len(os.sched_getaffinity(0)) + 2
                else:
                        processes_number = int(processes_number)
                        
                grouping = determine_grouping_setting(grouping_setting)
                process_repeats(input_path, countstable_file_name, repeatinfo_file_name, repeatinfo_gz_file_name, phobos_path, processes_number, grouping, grouping_motif_setting, min_perfection, max_perfection, min_rep_region_length, max_rep_region_length, min_unit_size, max_unit_size, min_copy_number, max_copy_number, multi_rep_reads_setting, read_whitelist_file, read_blacklist_file, read_chunk_size, add_phobos_arguments)

        elif args.program_mode == "summarizecounts":
                input_paths = args.inputpaths
                output_file = args.outputfile
                sample_names = args.sample_names
                ordered_samples, repeat_counts = read_countstables(input_paths, sample_names)
                write_count_matrix(ordered_samples, repeat_counts, output_file)
                
if __name__ == '__main__':
        main()
