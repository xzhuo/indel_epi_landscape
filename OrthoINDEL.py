'''
A script to parse crossmap.py results for SV project
'''
import sys
import argparse
import copy


class Region:
    def __init__(self):
        self.frags = []

    def start_line(self, line):
        list = line.split()
        self.from_chr = list[0]
        self.from_start = int(list[1])
        self.from_end = int(list[2])
        self.from_summit = "." if list[3] == "." else int(list[3])
        self.from_signal = list[4]
        self.from_strand = list[5]
        self.frags = []
        frag = Frag()
        frag.digest_line(line)
        self.frags.append(frag)

    def summit_frag_list(self):  # Return a list of all frags with summit (theoratically should only be one)
        return [x for x in self.frags if x.from_chr == self.from_chr and x.from_start <= self.from_summit and x.from_end >= self.from_summit]

    def length(self, assembly):  # assembly has to be "from" or "to"
        length = 0
        for frag in self.frags:
            length += frag.length(assembly)
        return length

    def combine_frags(self, size_limit, mini):
        ''' Combine frags based on chr, strand, and within size_limit. at least one of the gaps is < mini to be considered as continuous fragments.'''
        frag_dict = {}
        region_message = ''
        for frag in self.frags:
            key = frag.to_chr
            if key in frag_dict:
                frag_dict[key].append(frag)
            else:
                frag_dict[key] = [frag]
        if len(frag_dict) > 1:
            region_message = "Split_from_chromosome"
        combined_frag_dict = {}
        for key in frag_dict:
            split_pos = []
            for i in range(1, len(frag_dict[key])):
                mergeable, message = Region.can_merge(frag_dict[key][i - 1], frag_dict[key][i], size_limit, mini, False)
                if not mergeable:
                    split_pos.append(i)
                    if len(message) > len(region_message):
                        region_message = message
            if not len(split_pos):
                combined_frag_dict[key] = frag_dict[key]
            for i in range(len(split_pos)):
                key_index = key + str(i)
                if i == 0:
                    combined_frag_dict[key_index] = frag_dict[key][:split_pos[i]]
                else:
                    combined_frag_dict[key_index] = frag_dict[key][split_pos[i - 1]:split_pos[i]]
                if i + 1 not in range(len(split_pos)):
                    key_index = key + str(i + 1)
                    combined_frag_dict[key_index] = frag_dict[key][split_pos[i]:]
        self.frags = []
        for key_index in combined_frag_dict:
            merged_frag = Region.merge_all_frags(combined_frag_dict[key_index])
            self.frags.append(merged_frag)
        return region_message

    def keep_primary(self, perc, broad):
        '''If primary fragments found (frag with summit), then delele frags that are less then given perctage of total length'''
        self.from_summit
        length = self.length("from")
        can_remove = []
        found_primary = False
        for i, frag in enumerate(self.frags):
            if frag.length("from") / length < perc:
                can_remove.append(i)
            if broad or (self.from_summit > frag.from_start and self.from_summit < frag.from_end):  # broad is True if it is working with broadPeak. then remove small fragments anyway.
                found_primary = True
        if found_primary:
            self.frags = [i for j, i in enumerate(self.frags) if j not in can_remove]

    @staticmethod
    def merge_all_frags(frag_list):
        new_frag = Frag()
        new_frag.from_chr = [x.from_chr for x in frag_list][0]
        new_frag.from_strand = [x.from_strand for x in frag_list][0]
        new_frag.from_start = min([x.from_start for x in frag_list])
        new_frag.from_end = max([x.from_end for x in frag_list])
        new_frag.to_chr = [x.to_chr for x in frag_list][0]
        new_frag.to_strand = [x.to_strand for x in frag_list][0]
        new_frag.to_start = min([x.to_start for x in frag_list])
        new_frag.to_end = max([x.to_end for x in frag_list])
        return new_frag

    @staticmethod
    def can_merge(frag1, frag2, distance, mini, noindel):
        from_gap = abs(max(frag1.from_start, frag2.from_start) - min(frag1.from_end, frag2.from_end))
        to_gap = abs(max(frag1.to_start, frag2.to_start) - min(frag1.to_end, frag2.to_end))
        if noindel:
            mergeable = (frag1.from_chr == frag2.from_chr
                         and frag1.from_strand == frag2.from_strand
                         and frag1.to_chr == frag2.to_chr
                         and frag1.to_strand == frag2.to_strand
                         and max(from_gap, to_gap) < distance
                         and max(from_gap, to_gap) < mini
                         and (frag2.to_start >= frag1.to_start if (frag1.to_strand == frag1.from_strand) else frag2.to_end <= frag1.to_end))
            return mergeable, ''
        else:
            # return (frag1.from_chr == frag2.from_chr
            #         and frag1.from_strand == frag2.from_strand
            #         and frag1.to_chr == frag2.to_chr
            #         and frag1.to_strand == frag2.to_strand
            #         and max(from_gap, to_gap) < distance
            #         and min(from_gap, to_gap) < mini
            #         and (frag2.to_start >= frag1.to_start if (frag1.to_strand == frag1.from_strand) else frag2.to_end <= frag1.to_end))

            # evaluation:
            mergeable = False
            if (frag1.from_chr != frag2.from_chr):
                message = "Split_from_chromosome"
            elif (frag1.from_strand != frag2.from_strand):
                message = "Inverse_from_chr"
            elif (frag1.to_chr != frag2.to_chr):
                message = "Split_to_chr"
            elif (frag1.to_strand != frag2.to_strand):
                message = "Inv_to_chr"
            elif max(from_gap, to_gap) >= distance:
                message = "TooBigGap"
            elif min(from_gap, to_gap) >= mini:
                message = "Replacem"
            elif ((frag2.to_start < frag1.to_start if (frag1.to_strand == frag1.from_strand) else frag2.to_end > frag1.to_end)):
                message = "Overlap"
            else:
                mergeable = True
                message = ''
            return mergeable, message


class Frag:
    def __init__(self):
        pass

    def digest_line(self, line):
        list = line.split()
        self.to_chr = list[7]
        self.to_start = int(list[8])
        self.to_end = int(list[9])
        self.to_strand = list[12]
        if list[6] == '->':
            self.from_chr = list[0]
            self.from_start = int(list[1])
            self.from_end = int(list[2])
            self.from_strand = list[5]
        else:
            from_frag_list = list[6].strip('()').split(":")
            self.from_chr = from_frag_list[1]
            self.from_start = int(from_frag_list[2])
            self.from_end = int(from_frag_list[3])
            self.from_strand = from_frag_list[4]

    def merge_frags(self, frag):
        self.from_start = min(self.from_start, frag.from_start)
        self.from_end = max(self.from_end, frag.from_end)
        self.to_start = min(self.to_start, frag.to_start)
        self.to_end = max(self.to_end, frag.to_end)

    def length(self, assembly):  # assembly has to be "from" or "to"
        if assembly == "from":
            length = self.from_end - self.from_start
        elif assembly == "to":
            length = self.to_end - self.to_start
        return length


def _get_args():
    parser = argparse.ArgumentParser(description='simple arguments')
    parser.add_argument(
        '--input',
        '-i',
        action="store",
        dest="input",
        help='The input crossmap output file.',
    )
    parser.add_argument(
        '--distance',
        '-d',
        action="store",
        dest="distance",
        default=50,
        type=int,
        help='The distance used to define continous regions. fragments within a region are considered continous if the gap is smaller than the distance. Default is 50bp',
    )
    parser.add_argument(
        '--max',
        '-m',
        action="store",
        dest="max",
        default=50000,
        type=int,
        help='The distance used as upper limit of tolerated INDELs. Fragments separated by INDELs larger than the max distance will not be merged. Default is 50,000bp',
    )
    parser.add_argument(
        '--perc',
        '-p',
        action="store",
        dest="perc",
        default=0.2,
        type=float,
        help='The threshold to remove non-primary fragments (fraction of total region length)',
    )
    parser.add_argument(
        '--broad',
        '-b',
        action="store_true",
        dest="broad",
        default=False,
        help='If the input peaks are broadpeaks (broad peaks do not have summit defined)',
    )
    parser.add_argument(
        '--stringent',
        '-s',
        action="store_true",
        dest="stringent",
        default=True,
        help='If false, return splitted regions after final combination process. Default is True',
    )
    parser.add_argument(
        '--noindel',
        '-n',
        action="store_true",
        dest="noindel",
        default=False,
        help='If True, only return regions without any indel > distance',
    )
    return parser.parse_args()


def main():
    args = _get_args()
    summary_list = []
    with open(args.input, 'r') as Fh:
        regions = []
        last_region = Region()
        for line in Fh:
            line = line.rstrip()
            split = line.split()[6]
            if split != 'Fail' and split != 'Unmap':
                if not (hasattr(last_region, 'from_chr') and line.split()[0] == last_region.from_chr and int(line.split()[1]) == last_region.from_start and int(line.split()[2]) == last_region.from_end):
                    if len(last_region.frags) > 0:
                        regions.append(copy.deepcopy(last_region))
                    last_region.start_line(line)
                else:
                    frag = Frag()
                    frag.digest_line(line)
                    mergeable, info = Region.can_merge(last_region.frags[-1], frag, args.distance, args.distance, args.noindel)
                    if mergeable:
                        last_region.frags[-1].merge_frags(frag)
                    else:
                        last_region.frags.append(frag)
            else:
                summary_list.append(line)
        regions.append(copy.deepcopy(last_region))  # append the last entry.

    for region in regions:
        if args.noindel:
            if (len(region.frags) == 1 and abs(region.frags[0].from_start - region.from_start) < args.distance
               and abs(region.frags[0].from_end - region.from_end) < args.distance):
                print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
                      (region.from_chr, region.from_start, region.from_end, region.from_strand, region.from_summit, region.from_signal,
                       region.frags[0].to_chr, region.frags[0].to_start, region.frags[0].to_end, region.frags[0].to_strand))
        else:
            message = region.combine_frags(args.max, args.distance)
            num_frags = len(region.frags)
            if args.stringent:  # only one fragment in one region allowed if stringent.
                if num_frags == 1:
                    print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
                          (region.from_chr, region.from_start, region.from_end, region.from_strand, region.from_summit, region.from_signal,
                           region.frags[0].to_chr, region.frags[0].to_start, region.frags[0].to_end, region.frags[0].to_strand))
                else:
                    outcome = message
                    line = "\t".join([region.from_chr, str(region.from_start), str(region.from_end), region.from_strand, str(region.from_summit), str(region.from_signal), outcome])
                    summary_list.append(line)
            else:
                region.keep_primary(args.perc, args.broad)
                for index, frag in enumerate(region.frags):
                    if index == 0:
                        frag.from_start = region.from_start
                    if index == num_frags - 1:
                        frag.from_end = region.from_end
                    print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
                          (frag.from_chr, frag.from_start, frag.from_end, frag.from_strand, region.from_summit, region.from_signal,
                           frag.to_chr, frag.to_start, frag.to_end, frag.to_strand))
    with open("OrthoINDELUnmapped.txt", 'w') as Fh:
        Fh.writelines("\n".join(summary_list))


if __name__ == "__main__":
    main()
