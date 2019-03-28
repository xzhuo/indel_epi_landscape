'''
A script to process sv_calls result bed file from zev's DASVC pipeline. https://github.com/zeeev/DASVC
It will combine "insertion" and "deletion" in the same locus into a "replace" SV,
and get exact INDEL positions for query genome (the query starts and ends are starts and ends of alignment block in the sv_calls file).
'''

import argparse


def _get_args():
    parser = argparse.ArgumentParser(description='simple arguments')
    parser.add_argument(
        '--sv',
        '-s',
        action="store",
        dest="sv",
        help='The input sv calling file.',
    )
    parser.add_argument(
        '--bam',
        '-b',
        action="store",
        dest="bam",
        help='The input bam file.',
    )
    parser.add_argument(
        '--distance',
        '-d',
        action="store",
        dest="distance",
        help='The distance extended to both flanking region used for get_flanking_refinedSV.py.',
    )
    parser.add_argument(
        '--target_size',
        '-t',
        action="store",
        dest="target_size",
        help='The target chromosome size file used for get_flanking_refinedSV.py.',
    )
    parser.add_argument(
        '--query_size',
        '-q',
        action="store",
        dest="query_size",
        help='The query chromosome size file used for get_flanking_refinedSV.py.',
    )
    return parser.parse_args()


def main():
    args = _get_args()
    bam_file = args.bam
    with open(args.sv, 'r') as Fh:
        last_line = None
        for line in Fh:
            if line.startswith('#'):
                continue
            line = line.rstrip()
            linelist = line.split()
            if linelist[3] == "INS:BETWEEN":
                line = correct_reverse_between(bam_file, line)
                if last_line:
                    print(last_line)
                last_line = line
            elif linelist[3] == "DEL:BETWEEN":
                line = correct_reverse_between(bam_file, line)
                if last_line:
                    line = merge_line(last_line, line)
                    print(line)
                    last_line = None
                else:
                    print(line)
            else:
                if last_line:
                    print(last_line)
                    last_line = None
                target_name, target_start, target_end, sv_type, sv_length, per_id, matching_bases, query_name, query_start, query_end, sequence = splitline(line, "INTERNAL")
                query_start, query_end, target_start, target_end, sv_strand = get_query(bam_file, target_name, target_start, target_end, "no")  # if it is INTERNAL, get the precise insertion/deletion position in query.
                line = "\t".join([target_name, str(target_start), str(target_end), sv_type, str(sv_length), str(per_id), str(matching_bases), query_name, str(query_start), str(query_end), sv_strand, sequence])
                print(line)


def get_query(bam, target_name, target_start, target_end, stringent):
    '''
    Parse bam file using pysam, and get position alignment information.
    '''
    import pysam
    target_start -= 1
    samfile = pysam.AlignmentFile(bam, "rb")
    final_start = None
    final_end = None
    for num, read in enumerate(samfile.fetch(target_name, target_start, target_end)):
        position_pairs = read.get_aligned_pairs()
        # position_dict = {x[1]: x[0] for x in position_pairs}
        position_dict = dict()
        for x in position_pairs:
            if x[1] is not None:
                if x[0] is not None:
                    position_dict.update({x[1]: x[0]})
                else:
                    last_value = position_dict[x[1] - 1]  # last one!
                    position_dict.update({x[1]: last_value})

        # start_list = [x[0] for x in position_pairs if x[1] == target_start]
        # end_list = [x[0] for x in position_pairs if x[1] == target_end]

        # if len(start_list) == 1 and len(end_list) == 1:
        #     final_start = start_list[0] + read.get_tag("QS") + 1
        #     final_end = end_list[0] + read.get_tag("QS")
        # else:
        #     ipdb.set_trace()
        all_target_positions = position_dict.keys()
        target_min = min(all_target_positions)
        target_max = max(all_target_positions)
        strand_list = []
        if (stringent == "start" and target_start in position_dict) or (stringent == "end" and target_end in position_dict) or stringent == "no":
            if read.is_reverse:
                query_length = read.query_length
                strand_list.append("-")
                if target_start >= target_min:
                    final_start = query_length - position_dict.get(target_start) + read.get_tag("QS")
                else:
                    bailout_start = query_length - position_dict.get(target_min) + read.get_tag("QS")
                    alt_target_start = target_min

                if target_end <= target_max:
                    final_end = query_length - position_dict.get(target_end) + read.get_tag("QS") + 1
                else:
                    bailout_end = query_length - position_dict.get(target_max) + read.get_tag("QS") + 1
                    alt_target_end = target_max
            else:
                strand_list.append("+")
                if target_start >= target_min:
                    final_start = position_dict.get(target_start) + read.get_tag("QS") + 1
                else:
                    bailout_start = position_dict.get(target_min) + read.get_tag("QS") + 1
                    alt_target_start = target_min

                if target_end <= target_max:
                    final_end = position_dict.get(target_end) + read.get_tag("QS")
                else:
                    bailout_end = position_dict.get(target_max) + read.get_tag("QS")
                    alt_target_end = target_max

    # final_start = bailout_start if final_start is None else final_start
    # final_end = bailout_end if final_end is None else final_end
    if final_start is None:
        try:
            final_start = bailout_start
            target_start = alt_target_start
        except:
            pass

    if final_end is None:
        try:
            final_end = bailout_end
            target_end = alt_target_end
        except:
            pass
    strand_set = set(strand_list)
    return final_start, final_end, target_start + 1, target_end, ",".join(strand_set)


def splitline(line, type, SV=False):
    '''
    Split the line to a list.

    if type == "INTERNAL", per_id and matching_bases are number.
    if type == "BETWEEN", per_id and matching_bases are further splitted using comma.
    '''
    linelist = line.split()
    if SV:
        if len(linelist) == 12:
            target_name, target_start, target_end, sv_type, sv_length, per_id, matching_bases, query_name, query_start, query_end, sv_strand, sequence = linelist
        else:
            target_name, target_start, target_end, sv_type, sv_length, per_id, matching_bases, query_name, query_start, query_end, sv_strand = linelist
            sequence = ''
    elif len(linelist) == 11:  # if the last sequence column is missing go to else.
        target_name, target_start, target_end, sv_type, sv_length, per_id, matching_bases, query_name, query_start, query_end, sequence = linelist
    else:
        target_name, target_start, target_end, sv_type, sv_length, per_id, matching_bases, query_name, query_start, query_end = linelist
        sequence = ''
    target_start = int(target_start)
    target_end = int(target_end)
    query_start = int(query_start)
    query_end = int(query_end)
    sv_length = int(sv_length)
    if type == "INTERNAL":
        per_id = float(per_id)
        matching_bases = int(matching_bases)
    elif type == "BETWEEN":
        per_id = per_id.split(',')
        per_id[0] = float(per_id[0])
        per_id[1] = float(per_id[1])
        matching_bases = matching_bases.split(',')
        matching_bases[0] = int(matching_bases[0])
        matching_bases[1] = int(matching_bases[1])
    if SV:
        return target_name, target_start, target_end, sv_type, sv_length, per_id, matching_bases, query_name, query_start, query_end, sv_strand, sequence
    else:
        return target_name, target_start, target_end, sv_type, sv_length, per_id, matching_bases, query_name, query_start, query_end, sequence


def merge_line(last_line, line):
    '''
    Merge 2 lines to 1 line. the first line must be INS:BETWEEN, the 2nd line must be DEL:BETWEEN.

    The 2 lines merged and replaced with 1 line as a new sv_type: REPLACE.
    For lines cannot merge, return last_line\nline
    '''
    target_name, target_start, target_end, sv_type, sv_length, per_id, matching_bases, query_name, query_start, query_end, sv_strand, sequence = splitline(line, "BETWEEN", True)
    last_target_name, last_target_start, last_target_end, last_sv_type, last_sv_length, last_per_id, last_matching_bases, last_query_name, last_query_start, last_query_end, last_sv_strand, last_sequence = splitline(last_line, "BETWEEN", True)
    # always! last_sv_length == last_query_end - last_query_start
    if target_name == last_target_name and target_start == last_target_start:
        # print("same target position!")
        sv_strand = sv_strand if sv_strand == last_sv_strand else last_sv_strand + "," + sv_strand
        if query_name == last_query_name and (query_start == last_query_start and query_end == last_query_end):
            # print("same query position, merge!")
            sv_type = "REPLACE"
            sv_length = str(sv_length) + "," + str(last_sv_length)
            sequence = sequence + "," + last_sequence
            per_id = str(per_id[0]) + "," + str(per_id[1])
            matching_bases = str(matching_bases[0]) + "," + str(matching_bases[1])
            line = "\t".join([target_name, str(target_start), str(target_end), sv_type, str(sv_length), str(per_id), str(matching_bases), query_name, str(last_query_start), str(last_query_end), sv_strand, sequence])
        elif per_id == last_per_id and matching_bases == last_matching_bases:
            sv_type = "COMPLEX_REPLACE"
            sv_length = str(sv_length) + "," + str(last_sv_length)
            sequence = sequence + "," + last_sequence
            per_id = str(per_id[0]) + "," + str(per_id[1])
            matching_bases = str(matching_bases[0]) + "," + str(matching_bases[1])
            line = "\t".join([target_name, str(target_start), str(target_end), sv_type, str(sv_length), str(per_id), str(matching_bases), query_name, str(query_start), str(query_end), sv_strand, sequence])
        else:
            # ipdb.set_trace()
            line = "## more complicated SV!"
    else:
        line = last_line + "\n" + line
    return line

    # linelist[3] = "REPLACE"
    # linelist[8] = last_linelist[8]
    # linelist[9] = last_linelist[9]
    # linelist[4] = linelist[4] + "," + last_linelist[4]
    # linelist[10] = linelist[10] + "," + last_linelist[10]
    return line


def correct_reverse_between(bam, line):
    '''
    For del:between or ins:between lines, the query_end and query_start are reversed if it is alignment is reversed.
    This function correct query_start and query_end if the bam alignment is reverse.
    '''
    import pysam
    target_name, target_start, target_end, sv_type, sv_length, per_id, matching_bases, query_name, query_start, query_end, sequence = splitline(line, "BETWEEN")
    samfile = pysam.AlignmentFile(bam, "rb")
    read_list = []
    for read in samfile.fetch(target_name, target_start - 1, target_end + 20):  # the empty site is always < 20bp. so +- 20bp is enough here.
        read_list.append(read)

    reverse_list = [x.is_reverse for x in read_list]
    reverse_set = set(reverse_list)
    qs_list = [x.get_tag("QS") for x in read_list]
    qe_list = [x.get_tag("QE") for x in read_list]
    ts_list = [x.get_tag("TS") for x in read_list]
    te_list = [x.get_tag("TE") for x in read_list]
    joint_dict = {"frags": len(read_list),
                  "ori_qs": query_start,
                  "ori_qe": query_end,
                  "ori_ts": target_start,
                  "ori_te": target_end,
                  "reverse_set": reverse_set,
                  "qs_list": qs_list,
                  "qe_list": qe_list,
                  "ts_list": ts_list,
                  "te_list": te_list
                  }
    if len(read_list) == 1:
        sv_strand = "-" if read_list[0].is_reverse else "+"
        target_start = read_list[0].get_tag("TE")
        query_start = read_list[0].get_tag("QS") if read_list[0].is_reverse else read_list[0].get_tag("QE")
    elif len(read_list) == 2:
        if len(reverse_set) == 1:
            sv_strand = "-" if read_list[0].is_reverse else "+"
        else:
            sv_strand = ",".join(["-" if x else "+" for x in reverse_list])
        target_start = read_list[0].get_tag("TE")
        query_start = read_list[0].get_tag("QS") if read_list[0].is_reverse else read_list[0].get_tag("QE")
        target_end = read_list[1].get_tag("TS")
        query_end = read_list[1].get_tag("QE") if read_list[1].is_reverse else read_list[1].get_tag("QS")
    else:
        sv_strand = "3_frags"

    # TODO: 
    # for + strand:
    # TS = read1.TE, TE = read2.TS, QS = read1.QE, QE = read2.QS
    # for - strand:
    # TS = read1.TE, TE = read2.TS, QS = read1.QS, QE = read2.QE
    # for +,- :
    # TS = read1.TE, TE = read2.TS, QS = read1.QE, QE = read2.QE

    per_id = str(per_id[0]) + "," + str(per_id[1])
    matching_bases = str(matching_bases[0]) + "," + str(matching_bases[1])
    line = "\t".join([target_name, str(target_start), str(target_end), sv_type, str(sv_length), str(per_id), str(matching_bases), query_name, str(query_start), str(query_end), sv_strand, sequence])
    return line


if __name__ == '__main__':
    main()
