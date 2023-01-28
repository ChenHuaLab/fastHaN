#!/usr/bin/python3
# -*- coding: utf-8 -*-

#*************************************************************************
#    > File Name: Fasta2Phylip.py
#    > Author: xlzh
#    > Mail: xiaolongzhang2015@163.com 
#    > Created Time: 2023年01月26日 星期四 11时20分11秒
#*************************************************************************


import sys


def main():
    args = sys.argv

    if len(args) != 2:
        sys.stderr.write("usage: python Fasta2Phylip.py alignment.fasta > out.phylip\n")
        sys.exit(-1)

    align_fp = open(args[1], 'r')

    ''' align_dict = {'seq_name': [ATC...TTC, 'TCG...CCG'], ...}
    '''
    align_dict = {}
    cur_seq_name = None

    for line in align_fp:
        if line.startswith('>'):
            cur_seq_name = line.rstrip().strip('>')
            align_dict[cur_seq_name] = []
        else:
            align_dict[cur_seq_name].append(line.rstrip())

    # generate the phylip file
    n_align = len(align_dict)
    align_len = len(''.join(align_dict[cur_seq_name]))

    sys.stdout.write("%d %d\n" % (n_align, align_len))
    for seq_name in align_dict:
        sys.stdout.write("%s %s\n" % (seq_name, ''.join(align_dict[seq_name])))

    sys.stdout.flush()


if __name__ == '__main__':
    main()

